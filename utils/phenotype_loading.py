import csv
import subprocess
import os
import tempfile
import hail as hl


def pre_process_data_dictionary(pheno_description_raw_path, pheno_description_path):
    """
    Convert Data_Dictionary_Showcase.csv to tsv to enable load into hail

    :param str pheno_description_raw_path: Input file
    :param str pheno_description_path: Parsed tsv file
    """
    local_pheno_description_path = '/tmp/Data_Dictionary_Showcase.csv'
    local_pheno_description_out_path = '/tmp/Data_Dictionary_Showcase.tsv'
    hl.hadoop_copy(pheno_description_raw_path, f'file://{local_pheno_description_path}')
    with open(local_pheno_description_path) as f, open(local_pheno_description_out_path, 'w') as g:
        reader = csv.reader(f)
        for line in reader:
            g.write('\t'.join(line) + '\n')
    hl.hadoop_copy(f'file://{local_pheno_description_out_path}', pheno_description_path)


def get_codings():
    """
    Read codings data from Duncan's repo and load into hail Table

    :return: Hail table with codings
    :rtype: Table
    """
    root = f'{tempfile.gettempdir()}/PHESANT'
    if subprocess.check_call(['git', 'clone', 'https://github.com/astheeggeggs/PHESANT.git', root]):
        raise Exception('Could not clone repo')
    hts = []
    coding_dir = f'{root}/WAS/codings'
    for coding_file in os.listdir(f'{coding_dir}'):
        hl.hadoop_copy(f'file://{coding_dir}/{coding_file}', f'{coding_dir}/{coding_file}')
        ht = hl.import_table(f'{coding_dir}/{coding_file}')
        if 'node_id' not in ht.row:
            ht = ht.annotate(node_id=hl.null(hl.tstr), parent_id=hl.null(hl.tstr), selectable=hl.null(hl.tstr))
        ht = ht.annotate(coding_id=hl.int(coding_file.split('.')[0].replace('coding', '')))
        hts.append(ht)
    full_ht = hts[0].union(*hts[1:]).key_by('coding_id', 'coding')
    return full_ht.repartition(10)


def pheno_ht_to_mt(pheno_ht: hl.Table, data_type: str, special_fields: str = ('age', 'sex'),
                   pheno_function_type = hl.int, rekey: bool = True):
    """
    Input Hail Table with lots of phenotype row fields, distill into
    MatrixTable with either categorical or continuous data types
    as entries

    :param Table pheno_ht: Input hail Table with phenotypes as row fields
    :param str data_type: one of "categorical" or "continuous"
    :return: Hail MatrixTable with phenotypes as entries
    :rtype: MatrixTable
    """
    if data_type == 'categorical':
        filter_type = {hl.tbool}
        value_type = hl.bool
    else:
        filter_type = {hl.tint, hl.tfloat}
        value_type = hl.float

    special_fields_to_include = []
    fields = set(pheno_ht.row_value)
    for field in special_fields:
        if field in fields:
            fields.remove(field)
            special_fields_to_include.append(field)
    select_fields = {x: value_type(pheno_ht[x]) for x in fields if pheno_ht[x].dtype in filter_type}
    pheno_ht = pheno_ht.select(*special_fields_to_include, **select_fields)

    mt = pheno_ht.to_matrix_table_row_major(
        columns=list(select_fields), entry_field_name='value', col_field_name='phesant_pheno'
    )
    if rekey:
        mt = mt.key_cols_by(
            pheno=pheno_function_type(mt.phesant_pheno.split('_')[0]),
            coding=hl.case()
                .when(hl.len(mt.phesant_pheno.split('_')) == 1, mt.phesant_pheno)
                .when(hl.len(mt.phesant_pheno.split('_')) > 1, mt.phesant_pheno.split('_', 2)[1])  # TODO: fix to 1 when https://github.com/hail-is/hail/issues/7893 is fixed
                .or_error('A categorical was found not in the format of int_int')
        )
    return mt


def get_missing_codings(ht):
    """
    Download missing coding data files from UKB website

    :param Table ht: Input table with `meaning` and `coding_id`
    """
    missing_codings = ht.aggregate(hl.agg.filter(hl.is_missing(ht.meaning), hl.agg.collect_as_set(ht.coding_id)))
    import requests
    print(f'Missing: {missing_codings}')
    for coding in missing_codings:
        r = requests.post(url='http://biobank.ndph.ox.ac.uk/showcase/codown.cgi', data={'id': coding})
        with open(f'/tmp/coding{coding}.tsv', 'w') as f:
            f.write(r.text)


def get_all_codings():
    """
    Download all coding data files from UKB website
    """
    import requests
    coding_prefix = '/tmp/coding'
    all_codings = requests.post(url='http://biobank.ndph.ox.ac.uk/showcase/scdown.cgi', data={'fmt': 'txt', 'id': 2})
    all_codings = all_codings.text.strip().split('\n')[1:]
    hts = []
    for coding_list in all_codings:
        coding = coding_list.split('\t')[0]
        r = requests.post(url='http://biobank.ndph.ox.ac.uk/showcase/codown.cgi', data={'id': coding})
        req_data = r.text
        if r.status_code != 200 or not req_data or req_data.startswith('<!DOCTYPE HTML>'):
            print(f'Issue with {coding}: {r.text}')
            continue
        with open(f'{coding_prefix}{coding}.tsv', 'w') as f:
            f.write(req_data)
        hl.hadoop_copy(f'file://{coding_prefix}{coding}.tsv', f'{coding_prefix}{coding}.tsv')
        ht = hl.import_table(f'{coding_prefix}{coding}.tsv')
        if 'node_id' not in ht.row:
            ht = ht.annotate(node_id=hl.null(hl.tstr), parent_id=hl.null(hl.tstr), selectable=hl.null(hl.tstr))
        ht = ht.annotate(coding_id=hl.int(coding))
        hts.append(ht)
    full_ht = hts[0].union(*hts[1:]).key_by('coding_id', 'coding')
    return full_ht.repartition(10)


def get_phesant_reassignments(phesant_summary):
    """
    Helper function for add_coding_information.
    Parse PHESANT phenotype description data to get any coding reassignments.

    :param Table phesant_summary: Summary hail Table with PHESANT metadata
    :return: Table with reassignments
    :rtype: Table
    """
    phesant_summary = phesant_summary.annotate(
        reassign=phesant_summary['PHESANT.reassignments'].split(' ')[1].split('\|').map(lambda x: x.split('=')))
    ht = phesant_summary.explode('reassign')
    ht = ht.filter(ht.reassign[1] != 'NA')
    ht = ht.transmute(reassign_from=ht.reassign[0], reassign_to=ht.reassign[1])
    ht = ht.key_by(
        pheno=hl.int(ht.FieldID.split('_')[0]),
        coding=hl.or_missing(hl.len(ht.FieldID.split('_')) > 1, ht.FieldID.split('_')[1])
    )
    return ht.filter(ht.reassign_to == ht.coding)


def add_coding_information(mt: hl.MatrixTable, coding_ht: hl.Table, phesant_phenotype_info_path: str,
                           download_missing_codings: bool = False) -> hl.MatrixTable:
    """
    Add coding information from coding_ht as column annotations into mt

    :param MatrixTable mt: Input MT
    :param Table coding_ht: HT with coding information
    :param str phesant_phenotype_info_path: PHESANT phenotype metadata path
    :param bool download_missing_codings: Whether to download missing coding data
    :return: MT with coding information in column data
    :rtype: MatrixTable
    """
    mt = mt.annotate_cols(**coding_ht[(mt.coding_id, hl.str(mt.coding))])
    if download_missing_codings: get_missing_codings(mt.cols())
    phesant_summary = hl.import_table(phesant_phenotype_info_path, impute=True, missing='', key='FieldID')
    phesant_reassign = get_phesant_reassignments(phesant_summary)
    phesant_reassign = phesant_reassign.key_by('pheno', coding=hl.str(phesant_reassign.coding))
    mt = mt.annotate_cols(recoding=hl.or_missing(
        hl.is_missing(mt.meaning), phesant_reassign[mt.col_key].reassign_from
    ))
    return mt.annotate_cols(**hl.cond(hl.is_defined(mt.meaning),
                                      hl.struct(**{x: mt[x] for x in list(coding_ht.row_value)}),
                                      coding_ht[(mt.coding_id, hl.str(mt.recoding))]),
                            )


def combine_datasets(mt_path_dict: dict, summary_tsv_path_dict: dict = None,
                     pheno_description_path: str = None, coding_ht_path: str = None, data_type: str = 'categorical'):
    """
    Combine "both sexes", female, and male MTs into one with multiple entry fields,
    adding phenotype descriptions and coding.

    :param dict mt_path_dict: Dict of MTs (includes `both_sexes_no_sex_specific`, `females`, `males`)
    :param dict summary_tsv_path_dict: Dict of summary TSVs (includes `both_sexes_no_sex_specific`, `females`, `males`)
    :param str pheno_description_path: Phenotype description TSV path
    :param str coding_ht_path: Coding hail Table path
    :param str data_type: One of "categorical" or "continuous"
    :return: MatrixTable with all 3 entries combined
    :rtype: MatrixTable
    """
    both_mt = hl.read_matrix_table(mt_path_dict['both_sexes_no_sex_specific'])
    female_mt = hl.read_matrix_table(mt_path_dict['females'])
    male_mt = hl.read_matrix_table(mt_path_dict['males'])

    if pheno_description_path is not None:
        description_ht = hl.import_table(pheno_description_path, impute=True, missing='', key='FieldID')
        description_ht = description_ht.transmute(coding_id=description_ht.Coding)

        both_mt = both_mt.annotate_cols(**description_ht[both_mt.pheno])
        female_mt = female_mt.annotate_cols(**description_ht[female_mt.pheno])
        male_mt = male_mt.annotate_cols(**description_ht[male_mt.pheno])

    if coding_ht_path is not None and summary_tsv_path_dict is not None and data_type == 'categorical':
        coding_ht = hl.read_table(coding_ht_path)
        both_mt = add_coding_information(both_mt, coding_ht, summary_tsv_path_dict['both_sexes_no_sex_specific'])
        female_mt = add_coding_information(female_mt, coding_ht, summary_tsv_path_dict['females'])
        male_mt = add_coding_information(male_mt, coding_ht, summary_tsv_path_dict['males'])

    mt = hl.experimental.full_outer_join_mt(both_mt, female_mt)
    mt = mt.select_entries(
        both_sexes=mt.left_entry.value,
        females=mt.right_entry.value,
    ).drop('left_row', 'right_row')
    # ht = mt.cols().persist()
    # assert ht.all(ht.left_col == ht.right_col)
    mt = mt.select_cols(both_sexes_pheno=mt.left_col.drop(*mt.col_key), females_pheno=mt.right_col.drop(*mt.col_key))

    mt = hl.experimental.full_outer_join_mt(mt, male_mt)
    mt = mt.select_entries(**mt.left_entry, males=mt.right_entry.value).drop('left_row', 'right_row')
    mt = mt.select_cols(**mt.left_col.drop(*mt.col_key), males_pheno=mt.right_col.drop(*mt.col_key))
    return mt


def load_icd_data(pre_phesant_data_path, icd_codings_path, temp_directory,
                  force_overwrite_intermediate: bool = False,
                  include_dates: bool = False, icd9: bool = False):
    """
    Load raw (pre-PHESANT) phenotype data and extract ICD codes into hail MatrixTable with booleans as entries

    :param str pre_phesant_data_path: Input phenotype file
    :param str icd_codings_path: Input coding metadata
    :param str temp_directory: Temp bucket/directory to write intermediate file
    :param bool force_overwrite_intermediate: Whether to overwrite intermediate loaded file
    :param bool include_dates: Whether to also load date data (not implemented yet)
    :param bool icd9: Whether to load ICD9 data
    :return: MatrixTable with ICD codes
    :rtype: MatrixTable
    """
    if icd9:
        code_locations = {
            'primary_codes': '41203',
            'secondary_codes': '41205'
        }
    else:
        code_locations = {
            'primary_codes': '41202',
            'secondary_codes': '41204',
            'external_codes': '41201',
            'cause_of_death_codes': '40001'
        }
    date_locations = {
        'primary_codes': '41262'
    }
    ht = hl.import_table(pre_phesant_data_path, impute=not icd9, min_partitions=100, missing='', key='userId', types={'userId': hl.tint32})
    ht = ht.checkpoint(f'{temp_directory}/pre_phesant.ht', _read_if_exists=not force_overwrite_intermediate)
    all_phenos = list(ht.row_value)
    fields_to_select = {code: [ht[x] for x in all_phenos if x.startswith(f'x{loc}')] for code, loc in code_locations.items()}
    if include_dates:
        fields_to_select.update({f'date_{code}': [ht[x] for x in all_phenos if x.startswith(f'x{loc}')] for code, loc in date_locations.items()})
    ht = ht.select(**fields_to_select)
    ht = ht.annotate(
        **{code: ht[code].filter(lambda x: hl.is_defined(x)) for code in code_locations},
        # **{f'date_{code}': ht[code].filter(lambda x: hl.is_defined(x)) for code in date_locations}
    )
    # ht = ht.annotate(primary_codes_with_date=hl.dict(hl.zip(ht.primary_codes, ht.date_primary_codes)))
    all_codes = hl.sorted(hl.array(hl.set(hl.flatmap(lambda x: hl.array(x), ht.aggregate(
        [hl.agg.explode(lambda c: hl.agg.collect_as_set(c), ht[code]) for code in code_locations],
        _localize=True)))))
    ht = ht.select(
        bool_codes=all_codes.map(lambda x: hl.struct(**{code: ht[code].contains(x) for code in code_locations})))
    ht = ht.annotate_globals(all_codes=all_codes.map(lambda x: hl.struct(icd_code=x)))
    mt = ht._unlocalize_entries('bool_codes', 'all_codes', ['icd_code'])
    mt = mt.annotate_entries(any_codes=hl.any(lambda x: x, list(mt.entry.values())))
    # mt = mt.annotate_entries(date=hl.cond(mt.primary_codes, mt.primary_codes_with_date[mt.icd_code], hl.null(hl.tstr)))
    mt = mt.annotate_cols(truncated=False).annotate_globals(code_locations=code_locations)
    mt = mt.checkpoint(f'{temp_directory}/raw_icd.mt', _read_if_exists=not force_overwrite_intermediate)
    trunc_mt = mt.filter_cols((hl.len(mt.icd_code) == 3) | (hl.len(mt.icd_code) == 4))
    trunc_mt = trunc_mt.key_cols_by(icd_code=trunc_mt.icd_code[:3])
    trunc_mt = trunc_mt.group_cols_by('icd_code').aggregate_entries(
        **{code: hl.agg.any(trunc_mt[code]) for code in list(code_locations.keys()) + ['any_codes']}
    ).aggregate_cols(n_phenos_truncated=hl.agg.count()).result()
    trunc_mt = trunc_mt.filter_cols(trunc_mt.n_phenos_truncated > 1)
    trunc_mt = trunc_mt.annotate_cols(**mt.cols().drop('truncated', 'code_locations')[trunc_mt.icd_code],
                                      truncated=True).drop('n_phenos_truncated')
    mt = mt.union_cols(trunc_mt)
    coding_ht = hl.read_table(icd_codings_path)
    return mt.annotate_cols(**coding_ht[mt.col_key])


def get_full_icd_data_description(icd_codings_path, temp_path='hdfs://codings'):
    """
    Recursively get full ICD description

    :param str icd_codings_path: Input coding19.tsv
    :param str temp_path: Temporary path to write intermediate files
    :return: Table with full descriptions
    :rtype: Table
    """
    icd_ht = hl.import_table(icd_codings_path, impute=True, key='coding')
    icd_ht = icd_ht.annotate(original_node_id=icd_ht.node_id, original_parent_id=icd_ht.parent_id, short_meaning=icd_ht.meaning)
    icd_ht = icd_ht.checkpoint(f'{temp_path}/icd_codings.ht', overwrite=True)
    orig_icd_ht = hl.read_table(f'{temp_path}/icd_codings.ht').key_by('node_id')
    i = 0
    while True:
        icd_ht = icd_ht.checkpoint(f'{temp_path}/icd_codings_{i}.ht', overwrite=True)
        remaining = icd_ht.aggregate(hl.agg.count_where(icd_ht.parent_id != 0))
        print(f'Phenos remaining: {remaining}')
        if not remaining:
            break
        icd_join = orig_icd_ht[icd_ht.parent_id]
        icd_ht = icd_ht.annotate(meaning=hl.cond(icd_ht.parent_id == 0, icd_ht.meaning, icd_join.meaning + ' | ' + icd_ht.meaning),
                                 parent_id=hl.cond(icd_ht.parent_id == 0, icd_ht.parent_id, icd_join.parent_id))
        i += 1
    return icd_ht.transmute(node_id=icd_ht.original_node_id, parent_id=icd_ht.original_parent_id)


def load_prescription_data(prescription_data_tsv_path: str, prescription_mapping_tsv_path):
    ht = hl.import_table(prescription_data_tsv_path, types={'eid': hl.tint, 'data_provider': hl.tint}, key='eid')
    mapping_ht = hl.import_table(prescription_mapping_tsv_path, impute=True, key='Original_Prescription')
    ht = ht.annotate(issue_date=hl.cond(hl.len(ht.issue_date) == 0, hl.null(hl.tint64),
                                        hl.experimental.strptime(ht.issue_date + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT')),
                     **mapping_ht[ht.drug_name])
    ht = ht.filter(ht.Generic_Name != '').key_by('eid', 'Generic_Name', 'Drug_Category_and_Indication').collect_by_key()
    ht = ht.annotate(values=hl.sorted(ht.values, key=lambda x: x.issue_date))
    return ht.to_matrix_table(row_key=['eid'], col_key=['Generic_Name'], col_fields=['Drug_Category_and_Indication'])


def make_cooccurrence_ht(mt: hl.MatrixTable):
    mt = mt.annotate_cols(n_cases=hl.agg.sum(mt.any_codes))
    mt = mt.filter_cols(mt.n_cases >= 500).add_col_index()
    ht = mt.key_cols_by('col_idx').cols()
    bm = hl.linalg.BlockMatrix.from_entry_expr(mt.any_codes)
    bm = bm.T @ bm
    pheno_ht = bm.entries()
    pheno_ht = pheno_ht.annotate(i_data=ht[pheno_ht.i], j_data=ht[pheno_ht.j])
    pheno_ht = pheno_ht.annotate(prop_overlap=pheno_ht.entry / pheno_ht.i_data.n_cases)
    # pheno_ht = pheno_ht.checkpoint(f'{pheno_folder}/pheno_combo_explore/pheno_overlap.ht')


def make_correlation_ht(mt: hl.MatrixTable):
    mt = mt.filter_cols(mt.n_cases_both_sexes > 0).add_col_index()
    ht = mt.cols().key_by('col_idx')
    pheno = mt.value if 'value' in list(mt.entry) else mt.both_sexes
    bm = hl.linalg.BlockMatrix.from_entry_expr(pheno, mean_impute=True, center=True, normalize=True, axis='cols', block_size=512)
    bm = bm.T @ bm
    pheno_ht = bm.entries()
    pheno_ht = pheno_ht.annotate(i_data=ht[pheno_ht.i], j_data=ht[pheno_ht.j])
    return pheno_ht


def combine_pheno_files(pheno_file_dict: dict):
    full_mt: hl.MatrixTable = None
    for data_type, mt in pheno_file_dict.items():
        if 'pheno' in list(mt.col_key):
            mt = mt.key_cols_by(pheno=hl.str(mt.pheno), coding=mt.coding)
            criteria = mt.value if data_type == 'categorical' else hl.is_defined(mt.value)
            mt = mt.annotate_cols(n_cases=hl.agg.count_where(criteria))
            mt = mt.select_entries(value=hl.float64(mt.value))
        elif 'icd_code' in list(mt.col_key):
            mt = mt.key_cols_by(pheno=mt.icd_code, coding=mt.icd_version)
            mt = mt.filter_cols(mt.truncated)
            mt = mt.annotate_cols(n_cases=hl.agg.count_where(mt.any_codes))
            mt = mt.select_entries(value=hl.float64(mt.any_codes))
        elif 'phecode' in list(mt.col_key):
            mt = mt.key_cols_by(pheno=mt.phecode, coding=mt.phecode_sex)
            mt = mt.annotate_cols(n_cases=hl.agg.count_where(mt.case_control))
            mt = mt.select_entries(value=hl.float64(mt.case_control))
        elif 'Generic_Name' in list(mt.col_key):
            mt = mt.select_entries(value=hl.float64(hl.or_else(hl.len(mt.values) > 0, False)))
            mt2 = mt.group_cols_by(
                pheno=mt.Drug_Category_and_Indication,
                coding=mt.Drug_Category_and_Indication
            ).aggregate(value=hl.float64(hl.agg.any(mt.value > 0)))
            mt = mt.key_cols_by(pheno=mt.Generic_Name, coding=mt.Drug_Category_and_Indication).select_cols()
            mt = mt.union_cols(mt2)
            mt = mt.annotate_cols(n_cases=hl.int64(hl.agg.sum(mt.value)))
        else:
            raise ValueError('pheno or icd_code not in column key. New data type?')
        mt = mt.select_cols('n_cases', data_type=data_type,
                            n_defined=hl.agg.count_where(hl.is_defined(mt.value)))
        if full_mt is None:
            full_mt = mt
        else:
            full_mt = full_mt.union_cols(mt, row_join_type='outer' if data_type == 'prescriptions' else 'inner')
    full_mt = full_mt.unfilter_entries()
    return full_mt.select_entries(value=hl.cond(
        full_mt.data_type == 'prescriptions',
        hl.or_else(full_mt.value, hl.float64(0.0)),
        full_mt.value))


def combine_pheno_files_multi_sex(pheno_file_dict: dict, cov_ht: hl.Table):
    full_mt: hl.MatrixTable = None
    sexes = ('both_sexes', 'females', 'males')

    def compute_cases_binary(field, sex_field):
        return dict(
            n_cases_both_sexes=hl.agg.count_where(field),
            n_cases_females=hl.agg.count_where(field & (sex_field == 0)),
            n_cases_males=hl.agg.count_where(field & (sex_field == 1))
        )
    def format_entries(field, sex_field):
        return dict(
            both_sexes=hl.float64(field),
            females=hl.float64(hl.or_missing(sex_field == 0, field)),
            males=hl.float64(hl.or_missing(sex_field == 1, field))
        )

    for data_type, mt in pheno_file_dict.items():
        mt = mt.select_rows(**cov_ht[mt.row_key])
        print(data_type)
        if data_type == 'phecode':
            mt = mt.key_cols_by(pheno=mt.phecode, coding=mt.phecode_sex)
            mt = mt.annotate_cols(n_cases=hl.agg.count_where(mt.case_control))
            mt = mt.select_cols(**compute_cases_binary(mt.case_control, mt.sex),
                                data_type=data_type, meaning=mt.phecode_description,
                                path=mt.phecode_group)
            mt = mt.select_entries(**format_entries(mt.case_control, mt.sex))
        elif data_type == 'prescriptions':
            def format_prescription_name(pheno):
                return pheno.replace(',', '|').replace('/', '+').replace(' ', '')
            mt = mt.select_entries(value=hl.or_else(hl.len(mt.values) > 0, False))
            mt2 = mt.group_cols_by(
                pheno=format_prescription_name(mt.Drug_Category_and_Indication), coding='', Drug_Category_and_Indication=mt.Drug_Category_and_Indication,
            ).aggregate(value=hl.agg.any(mt.value))
            mt = mt.key_cols_by(pheno=format_prescription_name(mt.Generic_Name), coding='', Drug_Category_and_Indication=mt.Drug_Category_and_Indication).select_cols()
            mt = mt.union_cols(mt2).key_cols_by('pheno', 'coding')
            mt = mt.select_cols(**compute_cases_binary(mt.value, mt.sex),
                                data_type=data_type, meaning=mt.pheno,
                                path=mt.Drug_Category_and_Indication)
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
        elif data_type == 'biomarkers':
            mt = mt.key_cols_by(pheno=hl.str(mt.pheno), coding=mt.coding)
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(hl.is_defined(mt[sex])) for sex in sexes},
                                data_type=data_type, meaning=hl.null(hl.tstr), path=hl.null(hl.tstr))
        elif data_type == 'custom':
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
                hl.cond(mt.data_type == 'categorical', mt[sex] == 1.0, hl.is_defined(mt[sex]))
            ) for sex in sexes}, data_type=mt.data_type, meaning=hl.null(hl.tstr), path=hl.null(hl.tstr))
        elif data_type == 'additional':
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(hl.is_defined(mt[sex])) for sex in sexes},
                                data_type='continuous', meaning=hl.null(hl.tstr), path=hl.null(hl.tstr))
        elif 'pheno' in list(mt.col_key):
            mt = mt.key_cols_by(pheno=hl.str(mt.pheno), coding=mt.coding)

            def check_func(x):
                return x if data_type == 'categorical' else hl.is_defined(x)
            meaning_field = 'meaning' if data_type == 'categorical' else 'Field'
            path_field = 'Field' if data_type == 'categorical' else 'Path'
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(check_func(mt[sex])) for sex in sexes},
                                data_type=data_type,
                                meaning=hl.coalesce(*[mt[f'{sex}_pheno'][meaning_field] for sex in sexes]),
                                path=hl.coalesce(*[mt[f'{sex}_pheno'][path_field] for sex in sexes]))
            mt = mt.select_entries(**{sex: hl.float64(mt[sex]) for sex in sexes})

        elif 'icd_code' in list(mt.col_key):
            icd_version = mt.icd_version if 'icd_version' in list(mt.col) else ''
            mt = mt.key_cols_by(pheno=mt.icd_code, coding=icd_version)
            mt = mt.filter_cols(hl.len(mt.icd_code) == 3)
            mt = mt.select_cols(**compute_cases_binary(mt.any_codes, mt.sex),
                                data_type=data_type,
                                meaning=mt.short_meaning,
                                path=mt.meaning)
            mt = mt.select_entries(**format_entries(mt.any_codes, mt.sex))
        else:
            raise ValueError('pheno or icd_code not in column key. New data type?')

        if full_mt is None:
            full_mt = mt
        else:
            full_mt = full_mt.union_cols(mt, row_join_type='outer' if data_type == 'prescriptions' else 'inner')
    full_mt = full_mt.unfilter_entries()
    return full_mt.select_entries(**{sex: hl.cond(
        full_mt.data_type == 'prescriptions',
        hl.or_else(full_mt[sex], hl.float64(0.0)),
        full_mt[sex]) for sex in sexes})