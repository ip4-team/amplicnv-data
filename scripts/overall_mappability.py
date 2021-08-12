import argparse
import sys
from typing import List

import pandas


class GenomicRegion:
    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __str__(self):
        return f'{self.chrom}:{self.start}-{self.end}'


class Mappability:
    def __init__(self, value: float, region: GenomicRegion):
        self.value = value
        self.region = region

    def __str__(self):
        return f'{self.region}={self.value}'


def build_query(region: GenomicRegion) -> str:
    full_overlap = f'start >= {region.start} & end <= {region.end}'
    left_overlap = f'start < {region.start} & end > {region.start}'
    right_overlap = f'start < {region.end} & end > {region.end}'
    return f'chrom == "{region.chrom}" & (({full_overlap}) | ({left_overlap}) | ({right_overlap}))'


def create_mappability_obj(data) -> Mappability:
    return Mappability(
        data['mappability'],
        GenomicRegion(data['chrom'], data['start'], data['end'])
    )


def find_all_overlapping(region: GenomicRegion, data: pandas.DataFrame) -> List[Mappability]:
    return list(map(create_mappability_obj, data.query(build_query(region)).to_dict('records')))


def compute_overall_mappability(region: GenomicRegion, overlapping_mappability: List[Mappability]) -> float:
    base_counter = 0
    value_accumulator = 0
    for mappability in overlapping_mappability:
        region_start = region.start if region.start >= mappability.region.start else mappability.region.start
        region_end = region.end if region.end <= mappability.region.end else mappability.region.end
        region_len = region_end - region_start
        value_accumulator += value_accumulator + region_len * mappability.value
        base_counter += base_counter + region_len
    return value_accumulator / base_counter


def overall_mappability(region: GenomicRegion, data: pandas.DataFrame) -> float:
    print(f'Computing overall mappability for {region.chrom}:{region.start}-{region.end}')
    overlapping_mappability = find_all_overlapping(region, data)
    return compute_overall_mappability(region, overlapping_mappability)


def intervals_len(df):
    return df.apply(lambda row: abs(row['end'] - row['start']), axis='columns')


def get_mappability_table(targets_bed: str, mappability_bedgraph: str, k: int):
    bed_columns = ['chrom', 'start', 'end']
    targets_table = pandas.read_csv(targets_bed, sep='\t', skiprows=1, usecols=[0, 1, 2], names=bed_columns)
    map_table = pandas.read_table(mappability_bedgraph, names=bed_columns + ['mappability'])
    targets_table['mappability'] = targets_table.apply(
        lambda row: overall_mappability(
            GenomicRegion(row['chrom'], row['start'], row['end'] - (k - 1)), map_table), axis='columns')
    return targets_table


def parse_args(args) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Compute targets overall mappability')
    parser.add_argument('-t', '--targets', type=str, required=True, help='target BED file')
    parser.add_argument('-m', '--mappability', type=str, required=True, help='mappability bedgraph file')
    parser.add_argument('-k', '--k-mer-length', type=int, required=True,
                        help='k-mer len used when computing mappability')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file')
    return parser.parse_args(args)


if __name__ == '__main__':
    parsed_args = parse_args(sys.argv[1:])
    get_mappability_table(parsed_args.targets, parsed_args.mappability, parsed_args.k_mer_length)\
        .to_csv(parsed_args.output, sep='\t', index=False, header=False)
