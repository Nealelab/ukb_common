from gnomad_hail import *


def annotation_case_builder(worst_csq_by_gene_canonical_expr):
    return (hl.case(missing_false=True)
            .when(worst_csq_by_gene_canonical_expr.lof == 'HC', 'pLoF')
            .when(worst_csq_by_gene_canonical_expr.lof == 'LC', 'LC')
            .when((worst_csq_by_gene_canonical_expr.most_severe_consequence == 'missense_variant') |
                  (worst_csq_by_gene_canonical_expr.most_severe_consequence == 'inframe_insertion') |
                  (worst_csq_by_gene_canonical_expr.most_severe_consequence == 'inframe_deletion'), 'missense')
            .when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
            .or_missing())


def get_vep_formatted_data(ukb_vep_path: str):
    ht = hl.read_table(ukb_vep_path)
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    return ht.select(
        gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))


def load_variant_data(directory: str, pheno: str, coding: str, trait_type: str, ukb_vep_path: str, extension: str = 'single.txt',
                      n_cases: int = -1, n_controls: int = -1, overwrite: bool = False):
    output_ht_path = f'{directory}/variant_results.ht'
    ht = hl.import_table(f'{directory}/*.{extension}', delimiter=' ', impute=True)
    print(f'Loading: {directory}/*.{extension} ...')
    marker_id_col = 'markerID' if extension == 'single.txt' else 'SNPID'
    locus_alleles = ht[marker_id_col].split('_')
    ht = ht.key_by(locus=hl.parse_locus(locus_alleles[0]), alleles=locus_alleles[1].split('/'),
                   pheno=pheno, coding=coding, trait_type=trait_type).distinct().naive_coalesce(50)
    if marker_id_col == 'SNPID':
        ht = ht.drop('CHR', 'POS', 'rsid', 'Allele1', 'Allele2')
    ht = ht.transmute(Pvalue=ht['p.value']).annotate_globals(n_cases=n_cases, n_controls=n_controls)
    ht = ht.annotate(**get_vep_formatted_data(ukb_vep_path)[
        hl.struct(locus=ht.locus, alleles=ht.alleles)])  # TODO: fix this for variants that overlap multiple genes
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite).drop('n_cases', 'n_controls')
    mt = ht.to_matrix_table(['locus', 'alleles'], ['pheno', 'coding', 'trait_type'],
                            [marker_id_col, 'gene', 'annotation'], []).annotate_cols(n_cases=n_cases, n_controls=n_controls)
    mt.checkpoint(output_ht_path.replace('.ht', '.mt'), overwrite=overwrite, _read_if_exists=not overwrite)


def load_gene_data(directory: str, pheno: str, coding: str, trait_type: str, gene_ht_map_path: str,
                   n_cases: int = -1, n_controls: int = -1, overwrite: bool = False):
    output_ht_path = f'{directory}/gene_results.ht'
    print(f'Loading: {directory}/*.gene.txt ...')
    types = {f'Nmarker_MACCate_{i}': hl.tint32 for i in range(1, 9)}
    types.update({x: hl.tfloat64 for x in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_skato_NA', 'Pvalue_burden_NA', 'Pvalue_skat_NA')})
    ht = hl.import_table(f'{directory}/*.gene.txt', delimiter=' ', impute=True, types=types)

    fields = ht.Gene.split('_')
    gene_ht = hl.read_table(gene_ht_map_path).select('interval').distinct()
    ht = ht.key_by(gene_id=fields[0], gene_symbol=fields[1], annotation=fields[2],
                   pheno=pheno, coding=coding, trait_type=trait_type).drop('Gene').naive_coalesce(10).annotate_globals(n_cases=n_cases, n_controls=n_controls)
    ht = ht.annotate(total_variants=hl.sum([v for k, v in list(ht.row_value.items()) if 'Nmarker' in k]),
                     interval=gene_ht.key_by('gene_id')[ht.gene_id].interval)
    ht = ht.checkpoint(output_ht_path, overwrite=overwrite, _read_if_exists=not overwrite).drop('n_cases', 'n_controls')
    mt = ht.to_matrix_table(['gene_symbol', 'gene_id', 'annotation', 'interval'],
                            ['pheno', 'coding', 'trait_type'], [], []).annotate_cols(n_cases=n_cases, n_controls=n_controls)
    mt.checkpoint(output_ht_path.replace('.ht', '.mt'), overwrite=overwrite, _read_if_exists=not overwrite)


def get_cases_and_controls_from_log(log_prefix):
    cases = controls = -1
    for chrom in range(10, 23):
        try:
            with hl.hadoop_open(f'{log_prefix}_chr{chrom}_000000001.gene.log') as f:
                for line in f:
                    if line.startswith('Analyzing'):
                        fields = line.split()
                        if len(fields) == 6:
                            try:
                                cases = int(fields[1])
                                controls = int(fields[4])
                                break
                            except ValueError:
                                logger.warn(f'Could not load number of cases or controls from {line}.')
            return cases, controls
        except:
            pass
    return cases, controls


def union_mts_by_tree(all_mts, temp_dir):
    chunk_size = int(len(all_mts) ** 0.5) + 1
    outer_mts = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_mts): break
        mt = all_mts[i * chunk_size]
        for j in range(1, chunk_size):
            if i * chunk_size + j >= len(all_mts): break
            # try:
            mt = mt.union_cols(all_mts[i * chunk_size + j], row_join_type='outer')
            # except:
            #     print(f'problem with {i * chunk_size} and {i * chunk_size + j}')
            #     mt.describe()
            #     all_mts[i * chunk_size + j].describe()
            #     sys.exit(1)
        outer_mts.append(mt.checkpoint(f'{temp_dir}/temp_output_{i}.mt', overwrite=True))
    mt = outer_mts[0]
    for next_mt in outer_mts[1:]:
        mt = mt.union_cols(next_mt, row_join_type='outer')
    return mt