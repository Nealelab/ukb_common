import hail as hl
import uuid


def pull_out_fields_from_entries(mt, shared_fields, index='rows'):
    func = mt.annotate_rows if index == 'rows' else mt.annotate_cols
    mt = func(**{f'_{field}': hl.agg.take(mt[field], 1)[0] for field in shared_fields})
    return mt.drop(*shared_fields).rename({f'_{field}': field for field in shared_fields})


def create_broadcast_dict(key, value):
    """
    Create broadcast join (local dictionary from key -> value)
    from a Hail Table.

    :param Expression key: Key Expression
    :param Expression value: Value Expression
    :return: Hail DictExpression (without an index)
    :rtype: DictExpression
    """
    ht = key._indices.source
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