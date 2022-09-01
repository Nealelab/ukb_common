#!/usr/bin/env python3

__author__ = 'konradk'

import argparse
from pprint import pprint
from gnomad.utils import *
from ukbb_common import *

threshold: float = 1e-6
gwas_sig_threshold = 5e-8


def main(args):
    hl.init(log='/create_gwas_sig.log')

    if args.create_gwas_sig_file:
        hts = []
        for sex in ('both_sexes', 'female', 'male'):
            mt = hl.read_matrix_table(get_ukb_sumstats_mt_path(sex=sex))
            entry_fields = list(mt.entry)
            col_fields = list(mt.col)
            mt = mt.filter_entries((mt.pval <= threshold) & ~mt.low_confidence_variant)
            hts.append(mt.entries().annotate(sex=sex).naive_coalesce(1000).persist())
        ht = hts[0].union(*hts[1:]).key_by()
        # TODO: fix key, add biomarker data
        ht = ht.select(chrom=ht.locus.contig, pos=ht.locus.position, ref=ht.alleles[0], alt=ht.alleles[1], rsid=ht.rsid,
                       sex=ht.sex, AF=ht.AF, consequence=ht.consequence, info=ht.info,
                       **{f: ht[f] for f in entry_fields}, **{f: ht[f] for f in col_fields})
        ht = ht.checkpoint(get_gwas_sig_path(), overwrite=args.overwrite)
        ht.export(get_gwas_sig_path('tsv.bgz'))
        ht.drop('PHESANT_transformation', 'notes', 'ytx', 'tstat', 'rsid').export('/tmp/data.tsv.bgz')

    if args.create_top_p_file:
        hts = []
        for sex in ('both_sexes', 'female', 'male'):
            mt = hl.read_matrix_table(get_ukb_sumstats_mt_path(sex=sex))
            mt = mt.filter_entries(~mt.low_confidence_variant)
            ht = get_top_p_from_mt(mt, -hl.abs(mt.tstat)).annotate(sex=sex)
            hts.append(ht.filter(hl.is_defined(ht.pval)).naive_coalesce(1000).persist())
        ht = hts[0].union(*hts[1:])
        ht = ht.checkpoint(get_top_p_path(), overwrite=args.overwrite)
        row_fields = list(ht.row)
        ht = ht.key_by()
        ht = ht.select(chrom=ht.locus.contig, pos=ht.locus.position, ref=ht.alleles[0], alt=ht.alleles[1], *row_fields)
        ht.export(get_top_p_path('tsv.bgz'))

        ht = hl.read_table(get_top_p_path())
        print(f'Found {ht.aggregate(hl.agg.count_where((ht.pval < 5e-8) & (ht.sex == "both_sexes")))} significant hits '
              f'out of {int(ht.count() / 3)} variants.')

    ht = hl.read_table(get_gwas_sig_path())
    ht = ht.filter(ht.pval < gwas_sig_threshold)
    pprint(ht.aggregate(hl.struct(
        total_associations=hl.agg.count(),
        unique_phenos=hl.len(hl.agg.collect_as_set(ht.phenotype))
    )))
    ht = ht.group_by(ht.phenotype).aggregate(n_sig=hl.agg.count())
    ht.order_by(-ht.n_sig).show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--create_gwas_sig_file', action='store_true')
    parser.add_argument('--create_top_p_file', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
