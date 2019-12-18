import hail as hl
import argparse
import os
import sys
import tempfile


def gt_to_gp(mt, location: str = 'gp'):
    return mt.annotate_entries(**{location: hl.or_missing(
        hl.is_defined(mt.GT),
        hl.map(lambda i: hl.cond(mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0),
               hl.range(0, hl.triangle(hl.len(mt.alleles)))))})


def impute_missing_gp(mt, location: str = 'gp', mean_impute: bool = True):
    mt = mt.annotate_entries(_gp = mt.location)
    if mean_impute:
        mt = mt.annotate_rows(_mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp))
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop('_gp')


def main(args):
    hl.init(master=f'local[{args.n_threads}]',
            log=hl.utils.timestamp_path(os.path.join(tempfile.gettempdir(), 'extract_vcf'), suffix='.log'),
            default_reference='GRCh38')
    gene_ht = hl.read_table(args.gene_map_ht_path)
    mt = hl.read_matrix_table(args.mt_path)
    meta_ht = hl.read_table(args.meta_ht_path)
    meta_ht = meta_ht.filter(~meta_ht.is_filtered & (meta_ht.hybrid_pop == '12'))
    sample_mapping_ht = hl.import_table(args.sample_mapping_file, key='eid_sample', delimiter=',')
    qual_ht = hl.read_table(args.qual_ht_path)
    qual_ht = qual_ht.filter(hl.len(qual_ht.filters) == 0)

    if args.gene is not None:
        gene_ht = gene_ht.filter(gene_ht.gene_symbol == args.gene)
        interval = gene_ht.aggregate(hl.agg.take(gene_ht.interval, 1), _localize=False)
    else:
        interval = [hl.parse_locus_interval(args.interval)]
        gene_ht = hl.filter_intervals(gene_ht, interval)

    gene_ht = gene_ht.filter(hl.set(args.groups.split(',')).contains(gene_ht.annotation))
    gene_ht.select(group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation, variant=hl.delimit(gene_ht.variants, '\t')
                   ).key_by().drop('start').export(args.group_output_file, header=False)

    # TODO: possible minor optimization: filter output VCF to only variants in `gene_ht.variants`
    if not args.no_adj:
        mt = mt.filter_entries(mt.adj)

    mt = hl.filter_intervals(mt, interval).select_entries('GT')
    mt = mt.annotate_cols(exome_id=sample_mapping_ht[mt.s.split('_')[1]].eid_26041, meta=meta_ht[mt.col_key]).key_cols_by('exome_id')
    mt = mt.filter_cols(hl.is_defined(mt.meta) & hl.is_defined(mt.exome_id))
    mt = mt.filter_rows(hl.is_defined(qual_ht[mt.row_key]) & (hl.agg.count_where(mt.GT.is_non_ref()) > 0))
    mt = mt.annotate_rows(rsid=mt.locus.contig + ':' + hl.str(mt.locus.position) + '_' + mt.alleles[0] + '/' + mt.alleles[1])

    if args.callrate_filter:
        mt = mt.filter_rows(hl.agg.fraction(hl.is_defined(mt.GT)) > 0.95)

    if args.export_bgen:
        mt = gt_to_gp(mt)
        mt = impute_missing_gp(mt, mean_impute=args.mean_impute_missing)
        hl.export_bgen(mt, args.output_file, gp=mt.gp, varid=mt.rsid)
    else:
        mt = mt.annotate_entries(GT=hl.or_else(mt.GT, hl.call(0, 0)))
        # Note: no mean-imputation for VCF
        hl.export_vcf(mt, args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--gene_map_ht_path', help='Path to gene map HT', required=True)
    parser.add_argument('--mt_path', help='Path to full MT', required=True)
    parser.add_argument('--sample_mapping_file', help='Path to sample mapping file', required=True)
    parser.add_argument('--meta_ht_path', help='Path to sample meta HT', required=True)
    parser.add_argument('--qual_ht_path', help='Path to variant filters HT', required=True)

    parser.add_argument('--groups', help='Which variant groupings to use', required=True)
    parser.add_argument('--gene', help='Gene to export')
    parser.add_argument('--interval', help='Interval to export')

    parser.add_argument('--no_adj', help='Use all genotypes instead of only high-quality ones', action='store_true')
    parser.add_argument('--mean_impute_missing', help='Whether to mean impute missing genotypes (BGEN only) '
                                                      '(default: set to hom ref)', action='store_true')
    parser.add_argument('--export_bgen', help='Export BGEN instead of VCF', action='store_true')
    parser.add_argument('--callrate_filter', help='Impose 95% callrate filter', action='store_true')

    parser.add_argument('--group_output_file', help='Output file for variant groupings', required=True)
    parser.add_argument('--output_file', help='Output file', required=True)
    parser.add_argument('--n_threads', help='Number of threads', type=int, default=8)
    args = parser.parse_args()

    if int(args.gene is not None) + int(args.interval is not None) != 1:
        sys.exit('Error: One and only one of --gene or --interval must be specified')

    main(args)

