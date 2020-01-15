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
