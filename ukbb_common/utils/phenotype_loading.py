import csv
import subprocess
import os
import tempfile
import hail as hl
from ukbb_common.resources.generic import *

NULL_STR_KEY = ''
NULL_STR = hl.null(hl.tstr)


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


def pheno_ht_to_mt(pheno_ht: hl.Table, data_type: str, special_fields: str = ('age', 'sex'), rekey: bool = True):
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
            trait_type=data_type,
            phenocode=mt.phesant_pheno.split('_')[0],
            pheno_sex='both_sexes',
            coding=hl.case()
                .when((data_type == 'categorical') & (hl.len(mt.phesant_pheno.split('_')) > 1), mt.phesant_pheno.split('_', 2)[1])  # TODO: fix to 1 when https://github.com/hail-is/hail/issues/7893 is fixed
                .default(NULL_STR_KEY),
            modifier=hl.case()
                .when((data_type == 'continuous') & (hl.len(mt.phesant_pheno.split('_')) > 1),
                      mt.phesant_pheno.split('_', 2)[1])  # TODO: fix to 1 when https://github.com/hail-is/hail/issues/7893 is fixed
                .default(NULL_STR_KEY)
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
        pheno=ht.FieldID.split('_')[0],
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
    mt = mt.annotate_cols(recoding=hl.or_missing(
        hl.is_missing(mt.meaning), phesant_reassign[mt.col_key.select('phenocode', 'coding')].reassign_from
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
        description_ht = hl.import_table(pheno_description_path, impute=True, missing='', key='FieldID', types={'FieldID': hl.tstr})
        description_ht = description_ht.transmute(coding_id=description_ht.Coding)

        both_mt = both_mt.annotate_cols(**description_ht[both_mt.phenocode])
        female_mt = female_mt.annotate_cols(**description_ht[female_mt.phenocode])
        male_mt = male_mt.annotate_cols(**description_ht[male_mt.phenocode])

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


def make_pairwise_ht(mt: hl.MatrixTable, pheno_field, min_cases: int = 500, correlation: bool = False):
    mt = mt.annotate_entries(_pheno=pheno_field)
    if not correlation:
        mt = mt.annotate_cols(n_cases=hl.agg.sum(mt._pheno))
        mt = mt.filter_cols(mt.n_cases >= min_cases)
    else:
        mt = mt.filter_cols(mt.n_cases_both_sexes > 0)  # TODO: update this
    mt = mt.add_col_index()
    index_ht = mt.cols().key_by('col_idx')
    if correlation:
        bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, mean_impute=True, center=True, normalize=True, axis='cols', block_size=1024)
    else:
        bm = hl.linalg.BlockMatrix.from_entry_expr(mt._pheno, block_size=1024)
    bm = bm.T @ bm
    pheno_ht = bm.entries()
    pheno_ht = pheno_ht.annotate(i_data=index_ht[pheno_ht.i], j_data=index_ht[pheno_ht.j])
    if not correlation:
        pheno_ht = pheno_ht.annotate(prop_overlap=pheno_ht.entry / pheno_ht.i_data.n_cases)
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


def combine_pheno_files_multi_sex_legacy(pheno_file_dict: dict, cov_ht: hl.Table, truncated_codes_only: bool = True):
    full_mt: hl.MatrixTable = None
    sexes = ('both_sexes', 'females', 'males')

    for data_type, mt in pheno_file_dict.items():
        mt = mt.select_rows(**cov_ht[mt.row_key])
        print(data_type)
        if data_type == 'phecode':
            mt = mt.key_cols_by(trait_type=data_type, phenocode=mt.phecode, pheno_sex=mt.phecode_sex, coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            mt = mt.select_cols(**compute_cases_binary(mt.case_control, mt.sex),
                                description=mt.phecode_description, description_more=NULL_STR, coding_description=NULL_STR, category=mt.phecode_group)
            mt = mt.select_entries(**format_entries(mt.case_control, mt.sex))
        elif data_type == 'prescriptions':
            def format_prescription_name(pheno):
                return pheno.replace(',', '|').replace('/', '_')
            mt = mt.select_entries(value=hl.or_else(hl.len(mt.values) > 0, False))
            mt2 = mt.group_cols_by(
                trait_type=data_type,
                phenocode=format_prescription_name(mt.Drug_Category_and_Indication),
                pheno_sex='both_sexes',
                coding=NULL_STR_KEY, modifier=NULL_STR_KEY,
            ).aggregate(value=hl.agg.any(mt.value)).select_cols(category=NULL_STR)
            mt = mt.key_cols_by(
                trait_type=data_type,
                phenocode=format_prescription_name(mt.Generic_Name),
                pheno_sex='both_sexes',
                coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            mt = mt.select_cols(category=mt.Drug_Category_and_Indication)
            mt = mt.union_cols(mt2)
            mt = mt.select_cols(**compute_cases_binary(mt.value, mt.sex),
                                description=NULL_STR, description_more=NULL_STR, coding_description=NULL_STR, category=mt.category)
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
        elif data_type == 'custom':
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
                hl.cond(mt.trait_type == 'categorical', mt[sex] == 1.0, hl.is_defined(mt[sex]))
            ) for sex in sexes}, **{extra_col: mt[extra_col] if extra_col in list(mt.col) else NULL_STR
                                    for extra_col in PHENO_DESCRIPTION_FIELDS})
        elif data_type == 'additional':
            mt = mt.key_cols_by(trait_type='continuous', phenocode=mt.pheno, pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=mt.coding)
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(hl.is_defined(mt[sex])) for sex in sexes},
                                description=hl.coalesce(*[mt[f'{sex}_pheno'].meaning for sex in sexes]),
                                description_more=hl.coalesce(*[mt[f'{sex}_pheno'].description for sex in sexes]),
                                coding_description=NULL_STR, category=NULL_STR)
        elif data_type in ('categorical', 'continuous'):
            mt = mt.key_cols_by(trait_type=data_type, phenocode=hl.str(mt.pheno), pheno_sex='both_sexes',
                                coding=mt.coding if data_type == 'categorical' else NULL_STR_KEY,
                                modifier=NULL_STR_KEY if data_type == 'categorical' else mt.coding)

            def check_func(x):
                return x if data_type == 'categorical' else hl.is_defined(x)
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(check_func(mt[sex])) for sex in sexes},
                                description=hl.coalesce(*[mt[f'{sex}_pheno'].Field for sex in sexes]),
                                description_more=hl.coalesce(*[mt[f'{sex}_pheno'].Notes for sex in sexes]),
                                coding_description=hl.coalesce(*[mt[f'{sex}_pheno'].meaning for sex in sexes]) if
                                data_type == 'categorical' else NULL_STR,
                                category=hl.coalesce(*[mt[f'{sex}_pheno'].Path for sex in sexes]))
            mt = mt.select_entries(**{sex: hl.float64(mt[sex]) for sex in sexes})

        elif 'icd_code' in list(mt.col_key):
            icd_version = mt.icd_version if 'icd_version' in list(mt.col) else ''
            mt = mt.key_cols_by(trait_type=icd_version, phenocode=mt.icd_code, pheno_sex='both_sexes',
                                coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            if truncated_codes_only:
                mt = mt.filter_cols(hl.len(mt.icd_code) == 3)
                mt = mt.collect_cols_by_key()
                mt = mt.annotate_cols(keep=hl.if_else(
                    hl.len(mt.truncated) == 1, 0,
                    hl.zip_with_index(mt.truncated).filter(lambda x: x[1]).map(lambda x: x[0])[0]))
                mt = mt.select_entries(**{x: mt[x][mt.keep] for x in mt.entry})
                mt = mt.select_cols(**{x: mt[x][mt.keep] for x in mt.col_value if x != 'keep'})
            mt = mt.select_cols(**compute_cases_binary(mt.any_codes, mt.sex),
                                description=mt.short_meaning,
                                description_more="truncated: " + hl.str(mt.truncated) if 'truncated' in list(mt.col_value) else NULL_STR,
                                coding_description=NULL_STR,
                                category=mt.meaning)
            mt = mt.select_entries(**format_entries(mt.any_codes, mt.sex))
        elif data_type == 'icd_first_occurrence':
            mt = mt.select_entries(**format_entries(hl.is_defined(mt.value), mt.sex))
            mt = mt.select_cols(**compute_cases_binary(hl.is_defined(mt.both_sexes), mt.sex),
                                description=mt.Field, description_more=mt.Notes,
                                coding_description=NULL_STR, category=mt.Path)
        else: # 'biomarkers', 'activity_monitor'
            mt = mt.key_cols_by(trait_type=mt.trait_type if 'trait_type' in list(mt.col) else data_type,
                                phenocode=hl.str(mt.pheno), pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(hl.is_defined(mt[sex])) for sex in sexes},
                                description=mt.Field, description_more=NULL_STR, coding_description=NULL_STR, category=mt.Path)
        # else:
        #     raise ValueError('pheno or icd_code not in column key. New data type?')
        mt = mt.checkpoint(tempfile.mktemp(prefix=f'/tmp/{data_type}_', suffix='.mt'))
        if full_mt is None:
            full_mt = mt
        else:
            full_mt = full_mt.union_cols(mt, row_join_type='left')
            # full_mt = full_mt.union_cols(mt, row_join_type='outer' if data_type == 'prescriptions' else 'inner')
    full_mt = full_mt.unfilter_entries()

    # Here because prescription data was smaller than the others (so need to set the missing samples to 0)
    return full_mt.select_entries(**{sex: hl.cond(
        full_mt.trait_type == 'prescriptions',
        hl.or_else(full_mt[sex], hl.float64(0.0)),
        full_mt[sex]) for sex in sexes})



# TODO: move most of this into load functions
def combine_pheno_files_multi_sex(pheno_file_dict: dict, cov_ht: hl.Table, truncated_codes_only: bool = True):
    full_mt: hl.MatrixTable = None
    sexes = ('both_sexes', 'females', 'males')

    def counting_func(value, trait_type):
        return value if trait_type == 'categorical' else hl.is_defined(value)

    for data_type, mt in pheno_file_dict.items():
        mt = mt.select_rows(**cov_ht[mt.row_key])
        print(data_type)
        if data_type == 'custom':
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(
                hl.cond(mt.trait_type == 'categorical', mt[sex] == 1.0, hl.is_defined(mt[sex]))
            ) for sex in sexes}, **{extra_col: mt[extra_col] if extra_col in list(mt.col) else NULL_STR
                                    for extra_col in PHENO_DESCRIPTION_FIELDS})
        elif data_type in ('categorical', 'continuous'):

            def get_non_missing_field(mt, field_name):
                return hl.coalesce(*[mt[f'{sex}_pheno'][field_name] for sex in sexes])

            mt = mt.select_cols(**compute_cases_binary(counting_func(mt.both_sexes, data_type), mt.sex),
                                # **{f'n_cases_{sex}': hl.agg.count_where(counting_func(mt[sex], mt.trait_type)) for sex in sexes},
                                description=get_non_missing_field(mt, 'Field'),
                                description_more=get_non_missing_field(mt, 'Notes'),
                                coding_description=get_non_missing_field(mt, 'meaning') if
                                data_type == 'categorical' else NULL_STR,
                                category=get_non_missing_field(mt, 'Path'))
            mt = mt.select_entries(**{sex: hl.float64(mt[sex]) for sex in sexes})

        # TODO: got here - move some of this to ICD load (get icd_version as icd10 and move truncation steps there)
        elif 'icd_code' in list(mt.col_key):
            icd_version = mt.icd_version if 'icd_version' in list(mt.col) else 'icd10'
            mt = mt.key_cols_by(trait_type=icd_version, phenocode=mt.icd_code, pheno_sex='both_sexes',
                                coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            if truncated_codes_only:
                mt = mt.filter_cols(hl.len(mt.icd_code) == 3)
                mt = mt.collect_cols_by_key()
                mt = mt.annotate_cols(keep=hl.if_else(
                    hl.len(mt.truncated) == 1, 0,
                    hl.zip_with_index(mt.truncated).filter(lambda x: x[1]).map(lambda x: x[0])[0]))
                mt = mt.select_entries(**{x: mt[x][mt.keep] for x in mt.entry})
                mt = mt.select_cols(**{x: mt[x][mt.keep] for x in mt.col_value if x != 'keep'})
            mt = mt.select_cols(**compute_cases_binary(mt.any_codes, mt.sex),
                                description=mt.short_meaning,
                                description_more="truncated: " + hl.str(mt.truncated) if 'truncated' in list(mt.col_value) else NULL_STR,
                                coding_description=NULL_STR,
                                category=mt.meaning)
            mt = mt.select_entries(**format_entries(mt.any_codes, mt.sex))
        elif data_type == 'icd_first_occurrence':
            mt = mt.select_entries(**format_entries(hl.is_defined(mt.value), mt.sex))
            mt = mt.select_cols(**compute_cases_binary(mt.both_sexes == 1.0, mt.sex),
                                description=mt.Field, description_more=mt.Notes,
                                coding_description=NULL_STR, category=mt.Path)
        else: # 'biomarkers', 'activity_monitor'
            mt = mt.key_cols_by(trait_type=mt.trait_type if 'trait_type' in list(mt.col) else data_type,
                                phenocode=hl.str(mt.pheno), pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
            mt = mt.select_entries(**format_entries(mt.value, mt.sex))
            mt = mt.select_cols(**{f'n_cases_{sex}': hl.agg.count_where(hl.is_defined(mt[sex])) for sex in sexes},
                                description=mt.Field, description_more=NULL_STR, coding_description=NULL_STR, category=mt.Path)
        # else:
        #     raise ValueError('pheno or icd_code not in column key. New data type?')
        mt = mt.checkpoint(tempfile.mktemp(prefix=f'/tmp/{data_type}_', suffix='.mt'))
        if full_mt is None:
            full_mt = mt
        else:
            full_mt = full_mt.union_cols(mt, row_join_type='outer')
    full_mt = full_mt.unfilter_entries()

    # Here because prescription data was smaller than the others (so need to set the missing samples to 0)
    return full_mt.select_entries(**{sex: hl.cond(
        full_mt.trait_type == 'prescriptions',
        hl.or_else(full_mt[sex], hl.float64(0.0)),
        full_mt[sex]) for sex in sexes})

def load_dob_ht(pre_phesant_tsv_path, key_name = 'userId',
                year_field = 'x34_0_0', month_field = 'x52_0_0', recruitment_center_field = 'x54_0_0', quote=None):
    dob_ht = hl.import_table(pre_phesant_tsv_path, impute=False, min_partitions=100, missing='', key=key_name, quote=quote,
                             types={key_name: hl.tint32, recruitment_center_field: hl.tint32})
    year_field, month_field = dob_ht[year_field], dob_ht[month_field]
    month_field = hl.cond(hl.len(month_field) == 1, '0' + month_field, month_field)
    dob_ht = dob_ht.select(
        date_of_birth=hl.experimental.strptime(year_field + month_field + '15 00:00:00', '%Y%m%d %H:%M:%S', 'GMT'),
        month_of_birth=month_field,
        year_of_birth=year_field,
        recruitment_center=dob_ht[recruitment_center_field]
    )
    return dob_ht


def load_first_occurrence_data(first_exposure_and_activity_monitor_data_path, pre_phesant_tsv_path, delimiter: str = ',',
                               quote = '"'):
    ht = hl.import_table(first_exposure_and_activity_monitor_data_path, delimiter=delimiter, quote=quote,
                         missing='', impute=True, key='eid')  #, min_partitions=500)
    pseudo_dates = {'1901-01-01', '2037-07-07'}  # Pseudo date information at http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=819
    mt = filter_and_annotate_ukb_data(ht, lambda k, v: k.startswith('13') and len(k) == 10 and k.endswith('-0.0'), hl.str)

    dob_ht = load_dob_ht(pre_phesant_tsv_path, 'eid', year_field='34-0.0', month_field='52-0.0',
                         recruitment_center_field='54-0.0', quote=quote)[mt.row_key]
    dob = dob_ht.date_of_birth
    month = dob_ht.month_of_birth

    def parse_first_occurrence(x):
        return (hl.case(missing_false=True)
            .when(hl.is_defined(hl.parse_float(x)), hl.float64(x))  # Source of the first code ...
            .when(hl.literal(pseudo_dates).contains(hl.str(x)), hl.null(hl.tfloat64))  # Setting past and future dates to missing
            .when(hl.str(x) == '1902-02-02', 0.0)  # Matches DOB
            .when(hl.str(x) == '1903-03-03',  # Within year of birth (taking midpoint between month of birth and EOY)
                  (hl.experimental.strptime('1970-12-31 00:00:00', '%Y-%m-%d %H:%M:%S', 'GMT') -
                   hl.experimental.strptime('1970-' + month + '-15 00:00:00', '%Y-%m-%d %H:%M:%S',
                                            'GMT')) / 2)
            .default(hl.experimental.strptime(hl.str(x) + ' 00:00:00', '%Y-%m-%d %H:%M:%S', 'GMT') - dob
        ))
    mt = mt.annotate_entries(value=parse_first_occurrence(mt.value))

    mt = mt.key_cols_by(trait_type='icd_first_occurrence',
                        phenocode=mt.phenocode, pheno_sex='both_sexes',
                        coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
    return mt


def load_covid_data(all_samples_ht: hl.Table, covid_data_path: str, hes_main_path: str, hes_diag_path: str, death_path: str, wave: str = '01'):
    print(f'Loading COVID wave {wave}...')

    #### PATH OF DATA USED IN TESTING ####
    # def get_covid_data_path(wave: str = '20201103'):
    #     return f'gs://ukb31063/ukb31063.covid19_test_result.{wave}.txt'

    # def get_hes_main_data_path(wave: str = '20200909'):
    #     return f'gs://ukb31063/ukb31063.hesin.{wave}.txt'

    # def get_hes_diag_data_path(wave: str = '20200909'):
    #     return f'gs://ukb31063/ukb31063.hesin_diag.{wave}.txt'

    # def get_death_data_path(wave: str = '20201012'):
    #     return f'gs://ukb31063/ukb31063.death.{wave}.txt'

    covid_ht = hl.import_table(covid_data_path, delimiter='\t', missing='', impute=True, key='eid', min_partitions=100)
    hes_main_ht = hl.import_table(hes_main_path, delimiter='\t', missing='', impute=True, key=('eid', 'ins_index'), min_partitions=100)
    hes_diag_ht = hl.import_table(hes_diag_path, delimiter='\t', missing='', impute=True, key=('eid', 'ins_index'), min_partitions=100)
    death_ht = hl.import_table(death_path, delimiter='\t', missing='', impute=True, key='eid', min_partitions=100)

    death_ht = death_ht.annotate(death_date=hl.experimental.strptime(death_ht.date_of_death + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT'))  # Add Time Info to Death Register Data
    hes_main_ht = hes_main_ht.annotate(diagdate=hl.or_else(hes_main_ht.epistart, hes_main_ht.admidate))
    hes_main_ht = hes_main_ht.annotate(diag_date=hl.experimental.strptime(hes_main_ht.diagdate + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT'),
                                       admi_date=hl.experimental.strptime(hes_main_ht.admidate + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT'))  # Add Time Info to HES Inpatient Data

    hes_ht = hes_main_ht.join(hes_diag_ht)  # Join HES Inpatient Tables
    hes_pcr_pos = hes_ht.filter(hes_ht.diag_icd10 == 'U071')  # Subset Patients with Positive PCR Test Results
    hes_pcr_pos = hes_pcr_pos.select(covid_diag_date=hes_pcr_pos.diag_date, pcr_result=True)  # Select PCR-Positive Date info
    hes_ht = hes_ht.select('admi_date', 'diag_icd10')  # Select Admission Date and Diagnoses Info
    hes_ht = hes_ht.key_by('eid').join(hes_pcr_pos.key_by('eid'), how='outer')  # Create PCR-Positive information for all diagnoses

    hes_death_ht = hes_ht.join(death_ht, how='outer')  # Join Death Register Data to HES Info
    hes_death_ht = hes_death_ht.annotate(inpatient2=hes_death_ht.admi_date >= hes_death_ht.covid_diag_date,
                                         death=hes_death_ht.death_date >= hes_death_ht.covid_diag_date)  # Compare PCR-Positive Date to Death Date and Admission Date
    hes_death_ht = hes_death_ht.group_by('eid').aggregate(inpatient2=hl.agg.any(hes_death_ht.inpatient2),
                                                          death=hl.agg.any(hes_death_ht.death),
                                                          pcr_result=hl.agg.any(hes_death_ht.pcr_result),
    )

    covid_ht = covid_ht.group_by('eid').aggregate(
        origin=hl.agg.any(covid_ht.origin == 1),
        result1=hl.agg.any(covid_ht.result == 1),
        inpatient1=hl.agg.any(covid_ht.reqorg == 1),
    )
    covid_ht = covid_ht.join(hes_death_ht, how='outer')
    covid_ht = covid_ht.annotate(inpatient=covid_ht.inpatient1 | covid_ht.inpatient2,
                                 result=covid_ht.result1 | covid_ht.pcr_result)
    
    # TODO: add aoo parse to separate trait_type (covid_quantitative?)
    # dob = load_dob_ht(pre_phesant_tsv_path)[ht.key].date_of_birth
    # ht = ht.annotate(aoo=hl.or_missing(ht.result == 1, hl.experimental.strptime(ht.specdate + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT') - dob),
    #                  specdate=hl.experimental.strptime(ht.specdate + ' 00:00:00', '%d/%m/%Y %H:%M:%S', 'GMT')).drop('specdate')

    ht = all_samples_ht.annotate(**covid_ht[all_samples_ht.key])
    centers = hl.literal(ENGLAND_RECRUITMENT_CENTERS)

    analyses = {
	'A2_v2': hl.or_else(ht.death, False),
        'B1_v2': hl.or_missing(ht.result, ht.inpatient),  # fka ANA2
        'B1_v2_origin': hl.or_missing(ht.result, ht.origin),  # fka ANA2
        'C2_v2': hl.or_else(ht.result, False),  # fka ANA5
        'C2_v2_england_controls': hl.or_missing(centers.contains(ht.recruitment_center),  # fka ANA5_england_controls
                                               hl.or_else(ht.result, False)),
        'C1_v2': ht.result,  # fka ANA5_strict
        'B2_v2': hl.or_else(ht.result & ht.inpatient, False),  # fka ANA6
        'B2_v2_origin': hl.or_else(ht.result & ht.origin, False)  # fka ANA6
    }
    analysis_names = {
	'A2_v2': 'Death vs Survival (among COVID-19 positive)',
        'B1_v2': 'Hospitalized vs non-hospitalized (among COVID-19 positive)',  # fka ANA2
        'B1_v2_origin': 'Hospitalized vs non-hospitalized (among COVID-19 positive; old definition using "origin" field)',  # fka ANA2
        'C2_v2': 'COVID-19 positive (controls include untested)',  # fka ANA5
        'C2_v2_england_controls': 'COVID-19 positive (controls include untested), only patients from centers in England',  # fka ANA5_england_controls
        'C1_v2': 'COVID-19 positive (controls only COVID-19 negative)',  # fka ANA5_strict
        'B2_v2': 'Hospitalized vs non-hospitalized (controls include untested)',  # ANA6
        'B2_v2_origin': 'Hospitalized vs non-hospitalized (controls include untested; old definition using "origin" field)'  # ANA6
    }
    assert set(analyses.keys()) == set(analysis_names.keys())

    ht = ht.select(**analyses)
    mt = filter_and_annotate_ukb_data(ht, lambda k, v: True, annotate_with_showcase=False,
                                      format_col_name=lambda x: x)
    mt = mt.key_cols_by(trait_type='categorical', phenocode='COVID19', pheno_sex='both_sexes',
                        coding=mt.phenocode, modifier=wave)
    mt = mt.annotate_cols(description=hl.literal(analysis_names)[mt.coding])

    mt.annotate_cols(
        n_cases=hl.agg.count_where(mt.value == 1.0),
        n_controls=hl.agg.count_where(mt.value == 0.0)
    ).cols().show()

    return mt


def load_showcase(pheno_description_path):
    return hl.import_table(pheno_description_path, impute=True, missing='', key='FieldID', types={'FieldID': hl.tstr})


def load_activity_monitor_data(first_exposure_and_activity_monitor_data_path):
    ht = hl.import_table(first_exposure_and_activity_monitor_data_path, delimiter=',', quote='"', missing='', impute=True, key='eid')  #, min_partitions=500)
    quality_fields = ['90015-0.0', '90016-0.0', '90017-0.0']
    qual_ht = ht.select(hq=hl.is_missing(ht['90002-0.0']) & hl.all(lambda x: x == 1, [ht[x] for x in quality_fields]))
    mt = filter_and_annotate_ukb_data(ht, lambda x, v: x.startswith('90') and x.endswith('-0.0') and
                                                       v.dtype in {hl.tint32, hl.tfloat64})
    mt = mt.filter_cols(mt.ValueType == 'Continuous')
    mt = mt.annotate_rows(**qual_ht[mt.row_key])
    mt = mt.annotate_entries(value=hl.or_missing(hl.is_defined(mt.hq), mt.value))
    mt = mt.key_cols_by(trait_type='continuous', phenocode=mt.phenocode, pheno_sex='both_sexes', coding=NULL_STR_KEY, modifier=NULL_STR_KEY)
    return mt


def filter_and_annotate_ukb_data(ht, criteria, type_cast_function = hl.float64, annotate_with_showcase: bool = True,
                                 format_col_name = lambda x: x.split('-')[0]):
    fields_to_keep = {format_col_name(x): type_cast_function(v) for x, v in ht.row_value.items() if criteria(x, v)}
    ht = ht.select(**fields_to_keep)
    mt = ht.to_matrix_table_row_major(columns=list(fields_to_keep), entry_field_name='value', col_field_name='phenocode')
    if annotate_with_showcase:
        description_ht = load_showcase(pheno_description_path)
        mt = mt.annotate_cols(**description_ht[mt.phenocode])
    return mt
