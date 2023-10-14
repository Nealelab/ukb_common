""" Generate IHME GBD MatrixTable

This script is meant to be used with the R script 
210106_construct_combined_IHME_table.R, which outputs the zipped flat file to be used
as input here.

This script outputs a MatrixTable with rows representing diseases, columns representing
age groups, and entries representing various measures of morbidity. With default parameters,
the entry schema is a set of structs, each representing a measure (e.g., DALYs). Each
struct contains sexes, confidence intervals, and the values of the measure.

Optional arguments can be used to output a Table for just one of the age groups. Note
from the GBD guide (https://www.ncbi.nlm.nih.gov/books/NBK11818/):
    An age-standardized rate is a weighted average of the age-specific rates, where 
    the weights are the proportions of a standard population in the corresponding age 
    groups (q.v.). The potential confounding effect of age is removed when comparing
    age-standardized rates computed using the same standard population.

This can be imported to use the following functions:
convert_ihme_mt_to_table(mt: hl.MatrixTable, age_groups: list, explode_icd: str) -> hl.Table

@author: rgupta
"""
import hail as hl
import numpy as np
from copy import deepcopy
import argparse, re
__author__ = 'rgupta'


def get_entry_cols(entry_substr: str, field_names: list):
    entry_substr_l = entry_substr.split(',')
    matching_list_entries = []
    for ele in entry_substr_l:
        matching_items = [field for field in field_names if re.search(ele,field)]
        matching_list_entries = matching_list_entries + matching_items
    return list(set(matching_list_entries))


def get_column_cols(col_key: str, col_other: str, field_names: list):
    col_other_l = col_other.split(',')
    if not all([col in field_names for col in col_other_l]):
        raise ValueError('All column fields must be in the columns of the read table.')
    if not col_key in field_names:
        raise ValueError('Column key must be in set of fields.')
    return col_key, list(np.setdiff1d(list(set(col_other_l)),[col_key]))


def get_row_cols(row_key: str, row_other: list, col_fields: list, 
                 entry_fields: list, field_names: list):
    if not row_key in field_names:
        raise ValueError('Row key must be in set of fields.')
    if(len(row_other) == 0):
        row_other_l = np.setdiff1d(np.setdiff1d(field_names, col_fields), entry_fields)
    else:
        row_other_l = row_other.split(',')
        if not all([row in field_names for row in row_other_l]):
            raise ValueError('All row fields must be in the columns of the read table.')
    return row_key, list(np.setdiff1d(list(set(row_other_l)),[row_key]))


def recurse_struct(string_array: list, depth: int, ordering: list, 
                   elements: list, ht: hl.Table):
    """A recursive function to produce a nested struct

    Parameters
    ----------
    string_array : list
        An array that has the same length as the depth of the desired nested structure. Should contain
        None for any positions that haven't been visited yet. Once this array is filled (i.e.,
        at the deepest position), will be concatenated into a period-delimited string which should
        be found in the Hail Table.
    
    depth : int
        The depth of the strucure to currently iterate in.
    
    ordering : list
        This is used to map the elements of the column names (period delimited) to
        the struct. This function reverses this order, taking elements (which is 
        ordered in the order of the nested structure) and converting it back to the
        order of the original column name.
    
    elements : list
        A list of lists. Each element forms the names of that depth of the nested struct.
    
    ht : Hail Table
        Should contain columns that will be transformed into the new struct.

    Returns
    -------
    Hail Table
    """
    current_struct = hl.struct()
    struct_holder = dict()
    for this_ele in elements[depth]:
        string_array_m = deepcopy(string_array)
        pos_this = [idx for idx, ele in enumerate(ordering) if ele == depth+1][0]
        string_array_m[pos_this] = this_ele
        if depth == len(ordering)-1:
            # we are at the terminal level
            column_from = '.'.join(string_array_m)
            struct_holder.update({this_ele:ht[column_from]})
        else:
            depth_m = depth + 1
            new_struct = recurse_struct(string_array_m, depth_m, ordering, elements, ht)
            struct_holder.update({this_ele:new_struct})
    current_struct_m = current_struct.annotate(**struct_holder)
    return current_struct_m


def make_entries_into_struct(ht: hl.Table, entry_fields: list, 
                             entry_segment_order: list, drop_old_entries=True):
    entry_order_l = [int(ele) for ele in entry_segment_order.split(',')]
    entry_fields_spl = [field.split('.') for field in entry_fields]

    if not all([len(spl) == len(entry_order_l) for spl in entry_fields_spl]):
        raise ValueError('All entry columns must have the same number of segments ' +
                         '(period delimited), and this number must match the number ' +
                         'of (comma delimited) segments in the --entry-segment-order.')
    order_list_cp = deepcopy(entry_order_l)
    order_list_cp.sort()
    if not all([ordered == idx for ordered,idx in zip(order_list_cp, range(1, len(order_list_cp)+1))]):
        raise ValueError('--entry-segment-order must be a comma delimited set of ' +
                         'integers that are consecutive, starting from one.')

    elements = [] # indexed by the order
    for idx in range(0, max(entry_order_l)):
        this_vec = [ele for vec in entry_fields_spl for ele, order in zip(vec, entry_order_l) if order == idx+1]
        this_unique = list(set(this_vec))
        this_unique.sort()
        elements.append(this_unique)

    this_array = [None] * len(elements)
    resultant_struct = recurse_struct(string_array=this_array, depth=0, 
                                      ordering=entry_order_l, elements=elements, ht=ht)
    entry_holder = {k:v for k,v in resultant_struct.items()}

    ht_m = ht.annotate(**entry_holder)
    if drop_old_entries:
        ht_m = ht_m.drop(*entry_fields)
    
    return ht_m


def obtain_ICD_fields(row_fields: list):
    return [field for field in row_fields if re.search('ICD',field)]


def convert_ihme_mt_to_table(mt: hl.MatrixTable, age_groups: list, explode_icd: str = ''):
    """Convert the IHME MT to a Table.

    Parameters
    ----------
    mt : MatrixTable
        IHME MatrixTable.
    
    age_groups : list
        List of age_groups to filter the MT by. All elements must be present in the set of
        age_groups included in the MatrixTable.
    
    explode_icd : str
        A row field to explode. Must be an ICD code field. If provided, will remove all
        other ICD code fields.

    Returns
    -------
    Hail Table
    """
    age_vals = list(set(mt.age_name.collect()))
    age_groups = list(set(age_groups))
    age_groups_not_in_vals = np.setdiff1d(age_groups,age_vals)
    if len(age_groups_not_in_vals) == 0:
        mt_fc = mt_f.filter_cols(hl.literal(age_groups).contains(mt_f.age_name))
        if len(explode_icd) > 0:
            ICD_fields = obtain_ICD_fields(list(mt.row.keys()))
            if explode_icd in ICD_fields:
                ICD_non_explode = [ICD for ICD in ICD_fields if ICD != explode_icd]
                mt_explICD = mt_fc.drop(*ICD_non_explode)
                mt_fin = mt_explICD.explode_rows(mt_explICD[explode_icd])
            else:
                raise ValueError('If exploding the final table by one of the ICD columns,' +
                                 ' the provided column must actually contain ICD codes.')
        else:
            mt_fin = mt_fc
        return mt_fin.entries()
    else:
        raise ValueError('If trimming MatrixTable by an age group, provided age groups ' +
                         'must be present in the MatrixTable column fields.')


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Location of the table containing the data for merging. Expects this to contain data for a single year.')
parser.add_argument('--output', type=str, help='Output MatrixTable.')
parser.add_argument('--year', type=int, help='Year from which data was obtained.', default=2019)
parser.add_argument('--col-key', type=str, help='The column that will be the column field key.', default='age_name')
parser.add_argument('--col-other', type=str, help='Other columns to include as column fields, comma delimited.', default='age')
parser.add_argument('--row-key', type=str, help='The column that will be the column field key.', default='cause_name')
parser.add_argument('--row-other', type=str, help='Other columns to include as column fields, comma delimited. If blank, will include all other columns. Will convert any row fields matching ICD to an array by splitting on ", ".', default='')
parser.add_argument('--entry-substr', type=str, help='Substrings, column delimited, to define an entry column. Uses OR between provided substraings.', default='val,upper,lower')
parser.add_argument('--entry-segment-order', type=str, help='Entry columns are those that contain the substrings in --entry-subsr. They are period-delimited to represent hail structs. This comma-delimited string defines nesting order.', default='3,2,1')
parser.add_argument('--output-table', type=str, help='Output a derivative Hail Table to a directory. If not provided, will not output. If provided, expects --filter-specific-column to select a column for the final table.', default='')
parser.add_argument('--filter-specific-column', type=str, help='Filter the MatrixTable to a specific column for conversion to a Table. Only considered if --output-table is provided. Note that the column "Age-standardized" only contains rates.', default='Age-standardized')
parser.add_argument('--explode-icd', type=str, help='Explode an ICD column. If not provided, will not explode. If provided, will remove all other ICD columns. Only considered if --output-table is provided.', default='')


if __name__ == '__main__':
    args = parser.parse_args()
    
    # Import
    ht = hl.import_table(args.input, force=True, impute=True)
    field_names = list(ht.row.keys())

    entry_cols = get_entry_cols(args.entry_substr, field_names)
    column_key, column_other_fields = get_column_cols(args.col_key, args.col_other, field_names)
    row_key, row_other_fields = get_row_cols(args.row_key, args.row_other, 
                                             col_fields=[column_key] + column_other_fields,
                                             entry_fields=entry_cols, field_names=field_names)
    
    # Ensure that row, column, and entry schema contains appropriate names in the Hail Table.
    all_cols = [column_key, row_key] + row_other_fields + column_other_fields + entry_cols
    all_cols.sort()
    if not all([col in field_names for col in all_cols]):
        raise ValueError('All selected fields must be in the columns of the read table.')
    if any(all_cols.count(x) > 1 for x in all_cols):
        raise ValueError('There are duplicate items in the final set of columns. ' + \
                        'Ensure that row items are only assigned to rows, column ' + \
                        'items only to columns, etc.')

    # Make MatrixTable with proper entry structure.
    ht_m = make_entries_into_struct(ht, entry_cols, args.entry_segment_order, drop_old_entries=True)
    mt = ht_m.to_matrix_table(row_key=[row_key], row_fields=row_other_fields,
                              col_key=[column_key], col_fields=column_other_fields,
                              n_partitions=100)
    mt = mt.annotate_globals(year=args.year)

    # Splitting any row fields that have ICD in the name. Will not modify the row key.
    ICD_fields = obtain_ICD_fields(row_other_fields)
    modifications = {field: mt[field].split(', ') for field in ICD_fields}
    mt_f = mt.annotate_rows(**modifications)

    # Output MatrixTable
    mt_f.write(args.output, overwrite=True)

    # Now construct a derivative Hail table if the argument is enabled.
    if len(args.output_table) > 0:
        if len(args.filter_specific_column) > 0:
            print(args.filter_specific_column)
            ht_out = convert_ihme_mt_to_table(mt_f, age_groups=[args.filter_specific_column], 
                                              explode_icd=args.explode_icd)
            ht_out.write(args.output_table, overwrite=True)
        else:
            raise ValueError('If trimming MatrixTable by an age group, an age-group ' + 
                             'must be provided via --filter-specific-column.')

