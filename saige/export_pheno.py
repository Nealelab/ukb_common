import hail as hl
import argparse
import os
import tempfile


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'export_pheno'), suffix='.log'),
            default_reference='GRCh38')
    mt = hl.read_matrix_table(args.input_file)
    mt = mt.filter_cols((mt.pheno == args.pheno) & (mt.coding == args.coding))
    mt = mt.select_entries(value=mt[args.sex])
    if args.binary_pheno:
        mt = mt.select_entries(value=hl.int(mt.value))
    ht = mt.key_cols_by().entries()
    ht.export(args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_file', help='Input phenotype file', required=True)
    parser.add_argument('--binary_pheno', help='Whether to convert phenotype value to integer', action='store_true')
    parser.add_argument('--sex', help='Sex to use for pheno value', choices=('both_sexes', 'females', 'males'))
    parser.add_argument('--pheno', help='Pheno to output', required=True)
    parser.add_argument('--coding', help='Coding for pheno', default='')
    parser.add_argument('--output_file', help='Output file', required=True)
    parser.add_argument('--n_threads', help='Number of threads', type=int, default=8)
    args = parser.parse_args()

    main(args)
