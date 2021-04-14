
PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier')
PHENO_DESCRIPTION_FIELDS = ('description', 'description_more', 'coding_description', 'category')
PHENO_COLUMN_FIELDS = ('n_cases_both_sexes', 'n_cases_females', 'n_cases_males', *PHENO_DESCRIPTION_FIELDS)
PHENO_GWAS_FIELDS = ('n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized')

pheno_description_raw_path = 'gs://ukbb-exome-public/pheno_meta/Data_Dictionary_Showcase.csv'
pheno_description_path = 'gs://ukbb-exome-public/pheno_meta/Data_Dictionary_Showcase.tsv'
coding_ht_path = 'gs://ukbb-exome-public/pheno_meta/all_codings.ht'


def get_coding_path(coding_id: int, extension = 'tsv'):
    return f'gs://ukbb-exome-public/pheno_meta/codings/coding{coding_id}.{extension}'


icd_codings_tsv_path = get_coding_path(19)
icd9_codings_tsv_path = get_coding_path(87)
icd_codings_ht_path = get_coding_path(19, 'ht')
icd9_codings_ht_path = get_coding_path(87, 'ht')


PILOT_PHENOTYPES = set(map(lambda x: ('continuous', x, 'both_sexes', '', 'irnt'), {'50', '699', '23104'})).union(
    set(map(lambda x: ('categorical', x[0], 'both_sexes', x[1], ''), {('20004', '1095'), ('20004', '1479')}))).union(
    # set(map(lambda x: (x, '', 'icd'), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
    #                                    'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: ('icd10', x, 'both_sexes', '', ''), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                                       'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: ('phecode', x, 'both_sexes', '', ''), {'401', '411'}))).union(
    {('continuous','1717', 'both_sexes', '', ''),
     ('continuous', 'random', 'both_sexes', '', 'random'),
     ('continuous', 'random', 'both_sexes', '', 'random_strat'),
     ('continuous', 'whr', 'both_sexes', '', 'irnt'),
     ('categorical', '1747', 'both_sexes', '4', ''),
     ('continuous', '30040', 'both_sexes', '', 'irnt'),
     ('continuous', '30890', 'both_sexes', '', ''),
     ('prescriptions', 'HMG CoA reductase inhibitor|statin', 'both_sexes', '', ''),
     ('icd_first_occurrence', '131494', 'both_sexes', '', ''),
     ('icd_first_occurrence', '131306', 'both_sexes', '', '')})


ENGLAND_RECRUITMENT_CENTERS = {
11012, #	Barts
11021, #	Birmingham
11011, #	Bristol
11008, #	Bury
#11003	Cardiff
11024, #	Cheadle (revisit)
11020, #	Croydon
#11005	Edinburgh
#11004	Glasgow
11018, #	Hounslow
11010, #	Leeds
11016, #	Liverpool
11001, #	Manchester
11017, #	Middlesborough
11009, #	Newcastle
11013, #	Nottingham
11002, #	Oxford
11007, #	Reading
11014, #	Sheffield
10003, #	Stockport (pilot)
11006, #	Stoke
#11022	Swansea
#11023	Wrexham
11025, #	Cheadle (imaging)
11026, #	Reading (imaging)
11027, #	Newcastle (imaging)
11028 #	Bristol (imaging)
}