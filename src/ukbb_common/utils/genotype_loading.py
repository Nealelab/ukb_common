import hail as hl


def load_all_mfi_data():
    # Unverified: this file appears to have actual minor allele freqs, and then the next column seems to be minor allele
    # See https://gnomad.broadinstitute.org/variant/1-757640-G-A?dataset=gnomad_r2_1
    # 1:757640_G_A    rs3115853       757640  G       A       0.145479        G       0.995208
    # 1:757658_T_C    rs543105093     757658  T       C       0.000261269     C       0.23176
    header_fields = ('varid', 'rsid', 'pos', 'ref', 'alt', 'maf', 'minor_allele', 'info')
    types = {'f2': hl.tint32, 'f5': hl.tfloat64, 'f7': hl.tfloat64}
    header_rename = {f'f{i}': v for i, v in enumerate(header_fields)}
    hts = []
    for chromosome in CHROMOSOMES:
        ht = hl.import_table(ukb_imputed_info_path.format(chromosome), no_header=True, types=types).rename(header_rename)
        locus_chrom = 'X' if chromosome == 'XY' else chromosome
        ht = ht.key_by(locus=hl.locus(locus_chrom, ht.pos, reference_genome=REFERENCE_GENOME),
                       alleles=[ht.ref, ht.alt])
        hts.append(ht.annotate(original_chromosome=chromosome))
    return hts[0].union(*hts[1:])


def mac_category_case_builder(call_stats_expr):
    return (hl.case()
            .when(call_stats_expr.AC <= 5, call_stats_expr.AC)
            .when(call_stats_expr.AC <= 10, 10)
            .when(call_stats_expr.AC <= 20, 20)
            .when(call_stats_expr.AF <= 0.001, 0.001)
            .when(call_stats_expr.AF <= 0.01, 0.01)
            .when(call_stats_expr.AF <= 0.1, 0.1)
            .default(0.99))


def filter_ht_for_plink(ht: hl.Table, n_samples: int, min_call_rate: float = 0.95,
                         variants_per_mac_category: int = 2000,
                         variants_per_maf_category: int = 10000):
    from gnomad.utils.filtering import filter_to_autosomes
    ht = filter_to_autosomes(ht)
    ht = ht.filter((ht.call_stats.AN >= n_samples * 2 * min_call_rate) &
                   (ht.call_stats.AC > 0))
    ht = ht.annotate(mac_category=mac_category_case_builder(ht.call_stats))
    category_counter = ht.aggregate(hl.agg.counter(ht.mac_category))
    print(category_counter)
    ht = ht.annotate_globals(category_counter=category_counter)
    return ht.filter(hl.rand_unif(0, 1) <
                     hl.cond(ht.mac_category >= 1, variants_per_mac_category, variants_per_maf_category) / ht.category_counter[ht.mac_category]
                     )
