import argparse
import sys
import numpy as np
import pandas as pd
import scipy.stats as ss

INPUT_COLUMNS = [
    'chr',
    'start',
    'end',
    'name',
    'score',
    'strand',
    'junction',
    'gene_id',
    'gene_type',
    'gene'
]

GFF3_COLUMNS = [
    'seqid',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attributes'
]


def fisher(count12, count1, count2, count):
    table = [[count12, count1], [count2, count]]
    return ss.fisher_exact(table, alternative='greater').pvalue


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True)
parser.add_argument('-o', '--output', required=True)
parser.add_argument('-e', '--extended', action='store_true')
parser.add_argument('-a', '--annotation')


args = parser.parse_args()

columns = [f'{c}1' for c in INPUT_COLUMNS] + [f'{c}2' for c in INPUT_COLUMNS]
columns.remove('junction2')
data = pd.read_csv(
    sys.stdin if args.input == 'stdin' else args.input,
    sep='\t',
    header=None,
    names=columns
)

assert (data['gene_id1'] != data['gene_id2']).all()

gene_columns = ('gene_id', 'gene_type', 'gene')

new_columns = []
swap = {}
for n in 1, 2:
    for c in gene_columns:
        swap[f'{c}{n}'] = f'{c}{3-n}'
        new_columns.append(f'{c}{n}')

data = data[new_columns][:]
data = data[(data['gene_id1'] != '.') & (data['gene_id2'] != '.')][:]

data = pd.concat([data, data.rename(columns=swap)])

data['count1'] = data.groupby(['gene_id1', 'gene_type1', 'gene1'], sort=False).transform('size')
data['count2'] = data.groupby(['gene_id2', 'gene_type2', 'gene2'], sort=False).transform('size')

data = data[data['gene_id1'] > data['gene_id2']]

data = data.groupby(
    list(data.columns),
    sort=False,
    as_index=False
).agg('size').rename(columns={'size': 'observed'})

data['p_value'] = np.vectorize(fisher)(
    data['observed'],
    data['count1'],
    data['count2'],
    data['observed'].sum() - data['count2'] - data['count1'] - data['observed']
)

data['p_adj'] = ss.false_discovery_control(data['p_value'], method='bh')

data.sort_values(['p_adj', 'observed'], ascending=(True, False), inplace=True)

columns = ['gene_name', 'seqid', 'start', 'end', 'strand']
columns += ['source'] if args.extended else []
for n in 1, 2:
    data[[f'{c}{n}' for c in columns]] = data[f'gene{n}'].str.split(':', expand=True)

data['name'] = data['gene_id1'] + '__' + data['gene_id2']


if args.annotation:
    annotation = pd.read_csv(
        args.annotation, sep='\t', comment='#',
        header=None, names=GFF3_COLUMNS
    )
    annotation = annotation[annotation['type'] == 'gene'][:]
    pairs = annotation['attributes'].str.findall(r'([^=]+)=([^;]+);')
    annotation.reset_index().join(
        pd.json_normalize(pairs.apply(dict)),
        validate='one_to_one'
    )
    attributes = pd.json_normalize(pairs.apply(dict)).set_index(pairs.index)
    annotation = annotation.join(attributes).set_index('gene_id')

    for n in 1, 2:
        data[f'gene_type{n}'] = data.join(
            annotation['gene_type'],
            on=f'gene_id{n}',
            how='left',
            validate='many_to_one'
        )['gene_type']


data.to_csv(
    args.output,
    sep='\t',
    index=False,
    columns=[
        'seqid1', 'start1', 'end1', 'seqid2', 'start2', 'end2',
        'name', 'observed', 'strand1', 'strand2',
        'gene_name1', 'gene_type1', 'gene_name2', 'gene_type2',
        'p_value', 'p_adj'
    ]
)
