import hail as hl

def combine_phenotypes(mt: hl.MatrixTable, column_field, entry_field, lists_of_columns,
                       new_col_name='grouping', new_entry_name='new_entry', grouping_function=hl.agg.any):
    """
    Group by non-unique fields and apply grouping_function in order to combine entries in MatrixTable.

    Example:

    mt = hl.balding_nichols_model(1, 4, 10)
    mt = mt.annotate_entries(pheno=hl.rand_bool(0.5))
    lists_of_columns = [[0, 1], [0, 3]]
    entry_field = mt.pheno
    column_field = mt.sample_idx

    :param MatrixTable mt: Input MatrixTable
    :param Expression column_field: Column-indexed Expression to group by
    :param Expression entry_field: Entry-indexed Expression to which to apply `grouping_function`
    :param list of list lists_of_columns: Entry in this list should be the same type as `column_field`
    :param str new_col_name: Name for new column key (default 'grouping')
    :param str new_entry_name: Name for new entry expression (default 'new_entry')
    :param function grouping_function: Aggregator function to apply to `entry_field` (default hl.agg.any)
    :return: Re-grouped MatrixTable
    :rtype: MatrixTable
    """
    lists_of_columns = hl.literal(lists_of_columns)
    mt = mt._annotate_all(col_exprs={'_col_expr': column_field}, entry_exprs={'_entry_expr': entry_field})
    mt = mt.annotate_cols(**{new_col_name: lists_of_columns.filter(lambda x: x.contains(mt._col_expr))})
    mt = mt.explode_cols(new_col_name)
    return mt.group_cols_by(new_col_name).aggregate(**{new_entry_name: grouping_function(mt._entry_expr)})


def conditional_phenotypes(mt: hl.MatrixTable, column_field, entry_field, lists_of_columns,
                           new_col_name='grouping', new_entry_name='new_entry'):
    """
    Create a conditional phenotype by setting phenotype1 to missing for any individual without phenotype2.

    Pheno1 Pheno2 new_pheno
    T      T      T
    T      F      NA
    F      F      NA
    F      T      F

    `lists_of_columns` should be a list of lists (of length 2 for the inner list).
    The first element corresponds to the phenotype to maintain, except for setting to missing when the
    phenotype coded by the second element is False.

    new_entry = Pheno1 conditioned on having Pheno2

    Example:

    mt = hl.balding_nichols_model(1, 3, 10).drop('GT')
    mt = mt.annotate_entries(pheno=hl.rand_bool(0.5))
    lists_of_columns = [[0, 1], [2, 1]]
    entry_field = mt.pheno
    column_field = mt.sample_idx

    :param MatrixTable mt: Input MatrixTable
    :param Expression column_field: Column-indexed Expression to group by
    :param Expression entry_field: Entry-indexed Expression to which to apply `grouping_function`
    :param list of list lists_of_columns: Entry in this list should be the same type as `column_field`
    :param str new_col_name: Name for new column key (default 'grouping')
    :param str new_entry_name: Name for new entry expression (default 'new_entry')
    :return: Re-grouped MatrixTable
    :rtype: MatrixTable
    """
    assert all([len(x) == 2 for x in lists_of_columns])
    lists_of_columns = hl.literal(lists_of_columns)
    mt = mt._annotate_all(col_exprs={'_col_expr': column_field}, entry_exprs={'_entry_expr': entry_field})
    mt = mt.annotate_cols(_col_expr=lists_of_columns.filter(lambda x: x.contains(mt._col_expr)).map(lambda y: (y, y[0] == mt._col_expr)))
    mt = mt.explode_cols('_col_expr')
    # if second element (~mt._col_expr[1]) is false (~mt._entry_expr), then return missing
    # otherwise, get actual element (either true if second element, or actual first element)
    bool_array = hl.agg.collect(
        hl.if_else(~mt._col_expr[1] & ~mt._entry_expr, hl.null(hl.tbool), mt._entry_expr)
    )
    # if any element is missing, return missing. otherwise return first element
    return mt.group_cols_by(**{new_col_name: mt._col_expr[0]}).aggregate(
        **{new_entry_name: hl.if_else(hl.any(lambda x: hl.is_missing(x), bool_array),
                                      hl.null(hl.tbool),
                                      bool_array[0] & bool_array[1])})
