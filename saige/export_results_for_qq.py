#!/usr/bin/env python3

__author__ = 'konradk'

from ukbb_common import *
import argparse
import tempfile


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'export_results_for_qq'), suffix='.log'))

    ht = hl.read_table(args.input_dir + '/variant_results.ht').key_by()
    ht = ht.select(CHR=ht.locus.contig, BP=ht.locus.position, Pvalue=ht.Pvalue)
    ht.export(args.output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_dir', help='Input directory', required=True)
    parser.add_argument('--output_file', help='Output file', required=True)
    parser.add_argument('--n_threads', help='Number of threads to run', type=int, default=8)
    args = parser.parse_args()

    main(args)