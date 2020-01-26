#!/usr/bin/env python3

__author__ = 'konradk'

from gnomad_hail import *
from ukb_common import *

threshold: float = 1e-6

def main(args):
    hl.init(log='/create_gwas_sig.log')
    hts = []
    for sex in ('both_sexes', 'female', 'male'):
        mt = hl.read_matrix_table(get_ukb_sumstats_mt_path(sex=sex))
        entry_fields = list(mt.entry)
        mt = mt.filter_entries(mt.pval <= threshold)
        hts.append(mt.entries().annotate(sex=sex).persist())
    ht = hts[0].union(*hts[1:]).key_by()
    ht = ht.select(chrom=ht.locus.contig, pos=ht.locus.position, ref=ht.alleles[0], alt=ht.alleles[1], rsid=ht.rsid,
                   sex=ht.sex, AF=ht.AF, consequence=ht.consequence, info=ht.info, **{f: ht[f] for f in entry_fields})
    ht = ht.checkpoint(get_gwas_sig_path(), overwrite=args.overwrite)
    ht.export(get_gwas_sig_path('tsv.bgz'))

    hts = []
    for sex in ('both_sexes', 'female', 'male'):
        mt = hl.read_matrix_table(get_ukb_sumstats_mt_path(sex=sex))
        ht = mt.annotate_rows(top_p=hl.agg.take(mt.entry.annotate(**mt.col), 1, ordering=mt.pval)[0]).rows()
        hts.append(ht.transmute(**ht.top_p, sex=sex).persist())
    ht = hts[0].union(*hts[1:])
    ht.checkpoint(get_top_p_path(), overwrite=args.overwrite)
    row_fields = list(ht.row)
    ht = ht.key_by()
    ht = ht.select(chrom=ht.locus.contig, pos=ht.locus.position, ref=ht.alleles[0], alt=ht.alleles[1], *row_fields)
    ht.export(get_top_p_path('tsv.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
