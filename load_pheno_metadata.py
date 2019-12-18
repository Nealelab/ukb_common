#!/usr/bin/env python3

__author__ = 'konradk'

from ukb_phenotypes import *


def main(args):
    hl.init(log='/load_pheno_meta.log')

    pre_process_data_dictionary(pheno_description_raw_path, pheno_description_path)
    get_all_codings().write(coding_ht_path, overwrite=args.overwrite)
    get_full_icd_data_description(icd_codings_tsv_path).write(icd_codings_ht_path, args.overwrite)
    get_full_icd_data_description(icd9_codings_tsv_path).write(icd9_codings_ht_path, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
