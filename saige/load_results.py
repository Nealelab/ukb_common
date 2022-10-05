#!/usr/bin/env python3

__author__ = 'konradk'

from ukbb_common import *
import argparse
import tempfile

PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier')


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'load_results'), suffix='.log'),
            default_reference=args.reference)

    cases, controls = get_cases_and_controls_from_log(args.saige_run_log_format)

    quantitative_trait = args.trait_type in ('continuous', 'biomarkers')
    heritability = get_heritability_from_log(args.null_glmm_log, quantitative_trait) if args.null_glmm_log else -1.0
    inv_normalized = get_inverse_normalize_status(args.null_glmm_log) if args.null_glmm_log else 'NA'
    saige_version = get_saige_version_from_log(args.null_glmm_log) if args.null_glmm_log else 'NA'

    extension = 'single.txt' if args.analysis_type == 'gene' else 'single_variant.txt'
    pheno_key_dict = {k: getattr(args, k) for k in PHENO_KEY_FIELDS}
    if args.analysis_type == 'gene':
        load_gene_data(args.input_dir, pheno_key_dict, args.gene_map_ht_raw_path, cases, controls, heritability, saige_version, inv_normalized, args.overwrite)
    load_variant_data(args.input_dir, pheno_key_dict, args.ukb_vep_ht_path, extension, cases, controls, heritability, saige_version, inv_normalized,
                      args.log_pvalue, args.overwrite, args.legacy_annotations)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_dir', help='Input directory', required=True)
    parser.add_argument('--trait_type', help='Trait type', required=True)
    parser.add_argument('--phenocode', help='Phenotype ID', required=True)
    parser.add_argument('--pheno_sex', help='Phenotype sex', default='both_sexes')
    parser.add_argument('--coding', help='Phenotype coding', default='')
    parser.add_argument('--modifier', help='Phenotype modifier', default='')
    parser.add_argument('--null_glmm_log', help='Path to log file from null model')
    parser.add_argument('--saige_run_log_format', help='Path to log file from SAIGE test with {chr} placeholder', required=True)
    parser.add_argument('--analysis_type', help='Analysis type', choices=('gene', 'variant'), default='gene')
    parser.add_argument('--reference', help='Reference genome', default='GRCh38')
    parser.add_argument('--gene_map_ht_raw_path', help='Path to raw gene map')
    parser.add_argument('--ukb_vep_ht_path', help='Path to UKB VEP data', required=True)
    parser.add_argument('--n_threads', help='Number of threads to run', type=int, default=8)
    parser.add_argument('--legacy_annotations', help='Use old annotation picking (preferred for genotype data)', action='store_true')
    parser.add_argument('--log_pvalue', help='P-value was logged', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    args = parser.parse_args()

    main(args)