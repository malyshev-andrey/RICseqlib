#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import pandas as pd


BED_COLUMNS = ['chr', 'start', 'end', 'name', 'score', 'strand']
columns = BED_COLUMNS + [
    'junction',
    'cigar'
] + [
    f'gene_{c}' for c in BED_COLUMNS
] + [
    'source',
    'gene_id',
    'count',
    'coverage',
    'intersection'
]


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-p', '--pairwise_output', required=True)
parser.add_argument('-e', '--extended', action='store_true', default=False)
parser.add_argument('-s', '--preferred_source')


args = parser.parse_args()

if not args.extended:
    columns.remove('source')
    columns.remove('gene_id')
else:
    assert args.preferred_source

data = pd.read_csv(
    sys.stdin if args.input == 'stdin' else args.input,
    sep='\t',
    header=None,
    names=columns
)

data[['read_id', 'n']] = data['name'].str.split('/', expand=True)
data['n'] = data['n'].astype('int')
reads = data['read_id'].nunique()
names = data['name'].nunique()

# choose more frequent genes for particular read
data['read_gene_count'] = data.groupby(['read_id', 'gene_name'], sort=False).transform('size')
max_count = data.groupby('name', sort=False)['read_gene_count'].transform('max')
data = data[data['read_gene_count'] == max_count][:]

# choose preferred source
if args.extended:
    data['is_preferred'] = data['source'] == args.preferred_source
    preferred_available = data.groupby('name', sort=False)['is_preferred'].transform('any')
    data = data[data['is_preferred'] == preferred_available][:]

# choose gene with highest intersection fraction
data['intersection_frac'] = data['intersection'] / (data['end'] - data['start'])
max_frac = data.groupby('name', sort=False)['intersection_frac'].transform('max')
data = data[data['intersection_frac'] == max_frac][:]

# coverage voting
max_coverage = data.groupby('name', sort=False)['coverage'].transform('max')
data = data[data['coverage'] == max_coverage][:]

print('Inevitable duplicates:', data['name'].duplicated().sum())

# just first annotation for each part
data = data.groupby('name', sort=False, as_index=False).first()

assert data['name'].is_unique
assert data.groupby('read_id', sort=False)['n'].transform(lambda s: s == np.arange(len(s)) + 1).all()
assert (reads, names) == (data['read_id'].nunique(), data['name'].nunique())

data.to_csv(args.output, sep='\t', index=False, header=False, columns=columns)

data['gene'] = data['gene_name'].str.cat(
    others=[
        data['gene_chr'],
        data['gene_start'].astype('str'),
        data['gene_end'].astype('str'),
        data['gene_strand']
    ] + ([data['source']] if args.extended else []),
    sep=':'
)

if not args.extended:
    data['gene_id'] = data['gene_name']

columns = [
    'chr', 'start', 'end', 'name', 'score', 'strand', 'junction',
    'gene_id', 'gene_score', 'gene'
]

data = data[columns + ['read_id', 'n']]

pairwise = data.join(
    data.shift(-1),
    rsuffix='2',
    how='left',
    validate='one_to_one'
)

del data

pairwise = pairwise[pairwise['read_id'] == pairwise['read_id2']][:]
assert (pairwise['n'] + 1 == pairwise['n2']).all()

pairwise = pairwise[pairwise['gene'] != pairwise['gene2']]
pairwise = pairwise[pairwise['gene_id'] != pairwise['gene_id2']]

columns += [f'{c}2' for c in columns]

columns.remove('junction2')

pairwise.to_csv(args.pairwise_output, sep='\t', index=False, columns = columns)
