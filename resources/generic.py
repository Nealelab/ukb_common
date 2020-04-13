
pheno_description_raw_path = 'gs://ukbb-exome-public/pheno_meta/Data_Dictionary_Showcase.csv'
pheno_description_path = 'gs://ukbb-exome-public/pheno_meta/Data_Dictionary_Showcase.tsv'
coding_ht_path = 'gs://ukbb-exome-public/pheno_meta/all_codings.ht'


def get_coding_path(coding_id: int, extension = 'tsv'):
    return f'gs://ukbb-exome-public/pheno_meta/codings/coding{coding_id}.{extension}'


icd_codings_tsv_path = get_coding_path(19)
icd9_codings_tsv_path = get_coding_path(87)
icd_codings_ht_path = get_coding_path(19, 'ht')
icd9_codings_ht_path = get_coding_path(87, 'ht')


def get_gene_intervals_path(reference: str = 'GRCh37'):
    return f'{public_bucket}/misc/gene_intervals_{reference}.ht'


PILOT_PHENOTYPES = set(map(lambda x: (x, 'irnt', 'continuous'), {'50', '699', '23104'})).union(
    set(map(lambda x: (*x, 'categorical'), {('20004', '1095'), ('20004', '1479')}))).union(
    set(map(lambda x: (x, '', 'icd'), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                                       'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: (x, 'icd10', 'icd_all'), {'K519', 'K509', 'E109', 'E119', 'J459', 'I251',
                                       'K51', 'K50', 'E10', 'E11', 'J45', 'I25'}))).union(
    set(map(lambda x: (x, 'both_sexes', 'phecode'), {'401', '411'}))).union(
    {('1717', '1717', 'continuous'),
     ('random', 'random', 'continuous'),
     ('random', 'random_strat', 'continuous'),
     ('whr', 'whr', 'continuous'),
     ('1747', '4', 'categorical'),
     ('30040', 'irnt', 'continuous'),
     ('30890', '30890', 'biomarkers'),
     ('HMG CoA reductase inhibitor|statin', '', 'prescriptions')})
