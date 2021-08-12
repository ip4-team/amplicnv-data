import argparse
import configparser
import glob
import re
import sys
from typing import List

import pandas as pd
from bedhandler.handler import BedFileLoader
from cnvfinder.nrrhandler import NRRList


def get_files(param):
    param_list = param.strip().split('\n')
    return sum([glob.glob(path) for path in param_list], [])


def build_type_map(text_fields: List = None, number_fields: List = None, self_fields: List = None) -> dict:
    text_fields = [] if text_fields is None else text_fields
    number_fields = [] if number_fields is None else number_fields
    self_fields = [] if self_fields is None else self_fields

    text_dict = {field: str for field in text_fields}
    number_dict = {field: int for field in number_fields}
    self_dict = {field: lambda obj: obj for field in self_fields}
    return {**text_dict, **number_dict, **self_dict}


def build_column_map(columns: List):
    return dict(zip(columns, range(0, len(columns))))


def build_dataframe(lines: list, columns: list, column_map: dict, type_map: dict):
    targets = [[type_map[column](line[column_map[column]]) for column in columns] for line in lines]
    return pd.DataFrame(targets, columns=columns)


def build_df_for_attr(attr: str, data: List) -> pd.DataFrame:
    ids = ['S1']
    columns = ['chrom', 'chrom_start', 'chrom_end']
    new_df = data[0][columns + [attr]].copy()
    new_df.column_names = columns + ids
    for i in range(1, len(data)):
        new_id = f'S{i + 1}'
        new_df[new_id] = data[i][attr]
        ids.append(new_id)
        new_df.column_names = columns + ids
    return new_df.set_index(columns)


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Build target stats table')
    parser.add_argument('-m', '--overall-mappability', type=str, required=True, help='overall mappability BED file')
    parser.add_argument('-r', '--cov-files-regex', type=str, required=True,
                        help='amplicon coverage files regex (glob)')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    cov_files = {re.search(r'[a-zA-Z]{2}-.*-\d+', filename).group(0): filename
                 for filename in get_files(parsed_args.cov_files_regex)}
    bed_columns = ['chrom', 'chrom_start', 'chrom_end']
    mappability_table = pd.read_table(parsed_args.overall_mappability,
                                      header=None, names=bed_columns + ['mappability']).set_index(bed_columns)
    tables = []
    for sample_id, filename in cov_files.items():
        print(f'Loading sample {sample_id}: {filename}')
        cov_data = BedFileLoader(filename)
        df = build_dataframe(cov_data.expand_columns(), cov_data.columns,
                             build_column_map(cov_data.columns),
                             build_type_map(text_fields=['chrom', 'region_id', 'gene'],
                                            number_fields=['chrom_start', 'chrom_end', 'gc_count', 'overlaps',
                                                           'fwd_e2e', 'rev_e2e', 'total_reads', 'fwd_reads',
                                                           'rev_reads', 'cov20x', 'cov100x', 'cov500x'],
                                            self_fields=['pools']))
        tables.append(df)

    sample_ids = [f'S{i + 1}' for i in range(0, len(tables))]

    # read count
    read_count_table = build_df_for_attr('total_reads', tables)
    min_read_list = [20, 10, 5]
    results = []
    column_names = []

    stats = read_count_table.quantile([0.25, 0.75]).T.rename(columns={0.25: 'Q1', 0.75: 'Q3'})
    stats['IQR'] = stats['Q3'] - stats['Q1']
    stats['Q3 + 1.5 * IQR'] = stats['Q3'] + 1.5 * stats['IQR']
    stats['Q1 - 1.5 * IQR'] = stats['Q1'] - 1.5 * stats['IQR']

    above_outliers = read_count_table.apply(lambda c: c > stats.loc[c.name, 'Q3 + 1.5 * IQR'])
    stats['# outliers (IQR)'] = above_outliers.sum()
    stats['% outliers (IQR)'] = (stats['# outliers (IQR)'] / len(above_outliers) * 100).round(2)
    n_above_outliers = above_outliers.sum(axis='columns')
    percent_above_outliers = ((n_above_outliers / len(sample_ids)) * 100).round(2)
    results.append(list(percent_above_outliers))
    column_names.append('total reads > Q3 + 1.5 * IQR (%)')

    for min_read in min_read_list:
        below_outliers = read_count_table.apply(lambda c: c < min_read)
        stats[f'# outliers (< {min_read})'] = below_outliers.sum()
        stats[f'% outliers (< {min_read})'] = (stats[f'# outliers (< {min_read})'] / len(below_outliers) * 100).round(2)
        n_below_outliers = below_outliers.sum(axis='columns')
        percent_below_outliers = ((n_below_outliers / len(sample_ids)) * 100).round(2)
        results.append(list(percent_below_outliers))
        column_names.append(f'total reads < {min_read} (%)')

    nrr_list = NRRList(covfiles=list(cov_files.values()))
    nrr_list.compute_metrics()
    mad_div_median = [round((mad / median) * 100, 2) if median != 0 else 'NA'
                      for mad, median in zip(nrr_list.mad, nrr_list.normalized_median)]
    results.append(mad_div_median)
    column_names.append('mad / median (%)')

    for column_name, column_data in zip(column_names, results):
        mappability_table[column_name] = column_data

    mappability_table.to_csv(parsed_args.output, sep='\t')
