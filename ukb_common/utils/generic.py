import hail as hl
import uuid

UKB_GNOMAD_POP_MAPPING = {'AFR': 'afr',
                          'AMR': 'amr',
                          'CSA': 'sas',
                          'EAS': 'eas',
                          'EUR': 'nfe_nwe',
                          'MID': 'mid'}


def pull_out_fields_from_entries(mt, shared_fields, index='rows'):
    func = mt.annotate_rows if index == 'rows' else mt.annotate_cols
    mt = func(**{f'_{field}': hl.agg.take(mt[field], 1)[0] for field in shared_fields})
    return mt.drop(*shared_fields).rename({f'_{field}': field for field in shared_fields})


def create_broadcast_dict(key, value = None):
    """
    Create broadcast join (local dictionary from key -> value)
    from a Hail Table.

    :param Expression key: Key Expression
    :param Expression value: Value Expression
    :return: Hail DictExpression (without an index)
    :rtype: DictExpression
    """
    if isinstance(key, hl.Table):
        key = key.key
    ht = key._indices.source
    if value is None:
        value = ht.row_value
    return hl.dict(ht.aggregate(hl.agg.collect((key, value)), _localize=False))


def all_axis_join(left_mt, right_mt, row_join: str = '', col_join: str = '',
                  entry_join: str = '', global_join: str = ''):
    """
    Annotate `left_mt` with possible fields from `right_mt`. Note that the
    default parameters will overwrite any fields in the `left_mt`. To put
    these fields into a subfield, set the *_join parameter to a specified
    string.

    :param MatrixTable left_mt: Left MatrixTable
    :param MatrixTable right_mt: Right MatrixTable to annotate onto `left_mt`
    :param str row_join: Row field name for right_mt's row data. If empty string, splat into all subfields. If None, do no row field joining.
    :param str col_join: Col field name for right_mt's col data. If empty string, splat into all subfields. If None, do no col field joining.
    :param str entry_join: Entry field name for right_mt's entry data. If empty string, splat into all subfields. If None, do no entry field joining.
    :param str global_join: Global field name for right_mt's entry data. If empty string, splat into all subfields. If None, do no global field joining.
    :return: Fully joined MatrixTable
    :rtype: MatrixTable
    """
    row_data = {}
    if row_join is not None:
        row_data = right_mt.rows()[left_mt.row_key]
        row_data = hl.struct(**row_data) if row_join == '' else hl.struct(**{row_join: row_data})

    col_data = {}
    if col_join is not None:
        col_data = right_mt.cols()[left_mt.col_key]
        col_data = hl.struct(**col_data) if col_join == '' else hl.struct(**{col_join: col_data})

    entry_data = {}
    if entry_join is not None:
        entry_data = right_mt[left_mt.row_key, left_mt.col_key]
        entry_data = hl.struct(**entry_data) if entry_join == '' else hl.struct(**{entry_join: entry_data})

    global_data = {}
    if global_join is not None:
        global_data = right_mt.index_globals()
        global_data = hl.struct(**global_data) if global_join == '' else hl.struct(**{global_join: global_data})

    return left_mt._annotate_all(row_exprs=row_data, col_exprs=col_data,
                                 entry_exprs=entry_data, global_exprs=global_data)


def _load_gencode_gtf(gtf_file=None, reference_genome=None) -> hl.Table:
    """
    Get Gencode GTF (from file or reference genome)

    Parameters
    ----------
    reference_genome : :obj:`str` or :class:`.ReferenceGenome`, optional
       Reference genome to use (passed along to import_gtf).
    gtf_file : :obj:`str`
       GTF file to load. If none is provided, but `reference_genome` is one of
       `GRCh37` or `GRCh38`, a default will be used (on Google Cloud Platform).

    Returns
    -------
    :class:`.Table`
    """
    GTFS = {
        'GRCh37': 'gs://hail-common/references/gencode/gencode.v19.annotation.gtf.bgz',
        'GRCh38': 'gs://hail-common/references/gencode/gencode.v29.annotation.gtf.bgz',
    }
    if reference_genome is None:
        reference_genome = hl.default_reference().name
    else:
        reference_genome = reference_genome.name
    if gtf_file is None:
        gtf_file = GTFS.get(reference_genome)
        if gtf_file is None:
            raise ValueError(
                'get_gene_intervals requires a GTF file, or the reference genome be one of GRCh37 or GRCh38 (when on Google Cloud Platform)')
    ht = hl.experimental.import_gtf(gtf_file, reference_genome=reference_genome,
                                    skip_invalid_contigs=True, min_partitions=12)
    ht = ht.annotate(gene_id=ht.gene_id.split(f'\\.')[0],
                     transcript_id=ht.transcript_id.split('\\.')[0])
    return ht


def create_genome_intervals_file() -> hl.Table:
    # Load GTF file
    tmp_path = f'/tmp_{uuid.uuid4()}.ht'
    ht = _load_gencode_gtf()
    ht.filter((ht.feature == 'gene') & (ht.gene_type == 'protein_coding')).write(tmp_path, True)

    # Scan to get bounds, create intervals
    tmp_path2 = f'/tmp/tmp_{uuid.uuid4()}.ht'
    ht = hl.read_table(tmp_path)
    ht = ht.filter((ht.feature == 'gene') & (ht.gene_type == 'protein_coding'))
    ht = ht.select('gene_id', 'gene_name')
    last_locus = hl.scan.take(ht.row, 1, ordering=-ht.interval.start.global_position())
    intergenic_region = hl.or_missing((hl.len(last_locus) > 0) &
                                      (last_locus[0].interval.end.contig == ht.interval.start.contig),
                                      hl.interval(last_locus[0].interval.end, ht.interval.start))
    ht = ht.annotate(
        last_locus=last_locus,
        intergenic_region=intergenic_region
    )
    intergenic_length = ht.intergenic_region.end.position - ht.intergenic_region.start.position
    intergenic_region = hl.or_missing(intergenic_length > 0, ht.intergenic_region)
    intergenic_dist = hl.int((intergenic_region.end.position - intergenic_region.start.position) / 2)
    chrom = ht.interval.start.contig
    def interval(pos1, pos2):
        return hl.interval(hl.locus(chrom, pos1), hl.locus(chrom, pos2))
    ht = ht.transmute(
        intergenic_region1=hl.or_missing(hl.is_defined(intergenic_region),
                                         interval(intergenic_region.start.position + 1,  # gene interval is closed
                                                  intergenic_region.start.position + intergenic_dist)),
        intergenic_region2=hl.or_missing(hl.is_defined(intergenic_region),
                                         interval(intergenic_region.start.position + intergenic_dist,
                                                  intergenic_region.end.position))
    ).key_by()
    regions = hl.array([hl.struct(interval=ht.interval, gene_id=ht.gene_id, gene_name=ht.gene_name, within_gene=True)])
    regions = hl.if_else(hl.is_defined(ht.intergenic_region1), regions.extend([
        hl.struct(interval=ht.intergenic_region1, gene_id=ht.last_locus[0].gene_id, gene_name=ht.last_locus[0].gene_name, within_gene=False),
        hl.struct(interval=ht.intergenic_region2, gene_id=ht.gene_id, gene_name=ht.gene_name, within_gene=False)
    ]), regions)
    ht = ht.annotate(regions=regions).explode('regions')
    ht = ht.select(**ht.regions)
    return ht.key_by('interval')


def downsample_table_by_x_y(ht, x, y, label: dict = None, x_field_name: str = 'x', y_field_name: str = 'y',
                            n_divisions: int = 500):
    res = ht.aggregate(hl.agg.downsample(x, y, label=hl.array(list(label.values())), n_divisions=n_divisions), _localize=False)
    ht = hl.utils.range_table(1).annotate(data=res).explode('data')
    # ht = ht.drop('idx')  # TODO: add once https://github.com/hail-is/hail/issues/8751 is fixed
    ht = ht.select(**{x_field_name: ht.data[0], y_field_name: ht.data[1]},
                   **{l: ht.data[2][i] for i, l in enumerate(label.keys())})
    return ht


def locus_alleles_to_chr_pos_ref_alt(t, unkey_drop_and_add_as_prefix: bool = False):
    annotation = {
        'chrom': t.locus.contig, 'pos': t.locus.position,
        'ref': t.alleles[0], 'alt': t.alleles[1],
    }
    if unkey_drop_and_add_as_prefix:
        if isinstance(t, hl.MatrixTable):
            t = t.annotate_rows(**annotation).key_rows_by()
            return t.select_rows(*annotation, *t.row.drop('locus', 'alleles', *annotation))
        else:
            t = t.annotate(**annotation).key_by()
            return t.select(*annotation, *t.row.drop('locus', 'alleles', *annotation))
    else:
        return t.annotate_rows(**annotation) if isinstance(t, hl.MatrixTable) else t.annotate(**annotation)

