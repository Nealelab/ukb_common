import hail as hl
import argparse
import os
import tempfile
import sys
import importlib


PHENO_KEY_FIELDS = ('trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier')


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'export_pheno'), suffix='.log'),
            default_reference='GRCh38')

    sys.path.append('/')
    load_module = importlib.import_module(args.load_module)
    add_args = []
    if args.additional_args is not None:
        add_args = args.additional_args.split(',')
    mt = getattr(load_module, args.load_mt_function)(*add_args)

    mt = mt.filter_cols(hl.all(lambda x: x, [mt[k] == getattr(args, k, False) for k in PHENO_KEY_FIELDS if k != 'pheno_sex']))
    pheno_sex_mt = mt.filter_cols(mt.pheno_sex == args.pheno_sex)
    if pheno_sex_mt.count_cols() == 1:
        mt = pheno_sex_mt
    else:
        mt = mt.filter_cols(mt.pheno_sex == 'both_sexes')
    mt = mt.select_entries(value=mt[args.pheno_sex])
    if args.binary_trait:
        mt = mt.select_entries(value=hl.int(mt.value))
    if args.proportion_single_sex > 0:
        prop_female = mt.n_cases_females / (mt.n_cases_males + mt.n_cases_females)
        prop_female = prop_female.collect()[0]
        print(f'Female proportion: {prop_female}')
        if prop_female <= args.proportion_single_sex:
            print(f'{prop_female} less than {args.proportion_single_sex}. Filtering to males...')
            mt = mt.filter_rows(mt.sex == 1)
        elif prop_female >= 1 - args.proportion_single_sex:
            print(f'{prop_female} greater than {1 - args.proportion_single_sex}. Filtering to females...')
            mt = mt.filter_rows(mt.sex == 0)
    ht = mt.key_cols_by().select_cols().entries()
    ht.export(args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--load_module', help='Module to load with helper functions', required=True)
    parser.add_argument('--load_mt_function', help='Function in module to load analysis-ready MatrixTable', default='get_ukb_pheno_mt')
    parser.add_argument('--additional_args', help='Comma separated list of arguments to pass to mt_function')

    parser.add_argument('--trait_type', help='Which trait type to load', required=True)
    parser.add_argument('--phenocode', help='Pheno to output', required=True)
    parser.add_argument('--pheno_sex', help='Sex to use for pheno value', choices=('both_sexes', 'females', 'males'), default='both_sexes')
    parser.add_argument('--coding', help='Coding for pheno', default='')
    parser.add_argument('--modifier', help='Modifier for pheno', default='')
    parser.add_argument('--binary_trait', help='Whether to set trait to integer', action='store_true')
    parser.add_argument('--output_file', help='Output file', required=True)
    parser.add_argument('--proportion_single_sex', help='If set and proportion of male or female cases is less than '
                                                        'this number, then filter to females and males respectively',
                        type=float, default=0.0)
    parser.add_argument('--n_threads', help='Number of threads', type=int, default=8)
    args = parser.parse_args()

    main(args)
