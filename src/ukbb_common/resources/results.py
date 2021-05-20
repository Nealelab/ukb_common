

def get_ukb_sumstats_mt_path(reference: str = 'GRCh37', sex: str = 'both_sexes'):
    """
    Get UKB sumstats MatrixTable path

    :param str reference: Which reference to use (one of "GRCh37" or "GRCh38")
    :param str sex: Which sex to return results for (one of "both_sexes" (default), "female", "male")
    :return: Path to results MatrixTable
    :rtype: str
    """
    assert reference in ('GRCh37', 'GRCh38')
    assert sex in ('both_sexes', 'female', 'male')
    if reference == 'GRCh38': reference += '.liftover'
    return f'gs://hail-datasets/ukbb_imputed_v3_gwas_results_{sex}.{reference}.mt'


def get_ukb_sumstats_biomarkers_mt_path(sex: str = 'both_sexes', dilution_factor: bool = False):
    """
    Get UKB biomarker sumstats MatrixTable path

    :param str sex: Which sex to return results for (one of "both_sexes" (default), "female", "male")
    :param bool dilution_factor: Whether to return GWAS data after dilution factor applied
    :return: Path to results MatrixTable
    :rtype: str
    """
    dilution_factor = '.dilution_factor' if dilution_factor else ''
    return f'gs://ukb31063-mega-gwas/biomarkers/results-matrix-tables/ukb31063.biomarker_gwas_results.{sex}{dilution_factor}.mt'


def get_gwas_sig_path(data_format: str = 'ht', reference: str = 'GRCh37'):
    return f'gs://ukbb-mega-gwas-results-public/round2/ukbb_imputed_v3_gwas_significant.{reference}.{data_format}'


def get_top_p_path(data_format: str = 'ht', reference: str = 'GRCh37'):
    return f'gs://ukbb-mega-gwas-results-public/round2/ukbb_imputed_v3_top_p_per_variant.{reference}.{data_format}'


def check_trait_types(trait_type):
    _TRAIT_TYPES = ('icd', 'continuous', 'categorical', 'all')
    if trait_type not in _TRAIT_TYPES:
        raise ValueError(f'{trait_type} not in {_TRAIT_TYPES}')


def check_timing_type(timing_type):
    _TIMING_TYPES = ('saige', 'null_model', 'full')
    if timing_type not in _TIMING_TYPES:
        raise ValueError(f'{timing_type} not in {_TIMING_TYPES}')