#!/usr/bin/env python3
import argparse
import sys
import warnings
from typing import Generator, Callable
import pandas as pd

warnings.filterwarnings("ignore")

CHUNKSIZE = 10**7
OUTPUT_BUFFER = 536870912


def adjacent_groupby2(*args, chunk_processing: Callable = None, 
                      name_column, **kwargs) -> Generator[pd.DataFrame, None, None]:
    if 'chunksize' not in kwargs:
        kwargs['chunksize'] = CHUNKSIZE
        
    reader = pd.read_csv(*args, **kwargs)
    remaining = next(reader)
    remaining['read_id'] = remaining[name_column].str.split('/', expand=True)[0]
    
    for chunk in reader:
        chunk['read_id'] = chunk[name_column].str.split('/', expand=True)[0]
        first = (chunk['read_id'] == chunk['read_id'].iloc[0])
        yield pd.concat(
            [remaining, chunk[first]],
            verify_integrity=True
        )
        remaining = chunk[~first]
    yield remaining


def classify_reads(df: pd.DataFrame, name_column, gene_column):
    
    df['without_gene'] = (df[gene_column] == '.')
    df['without_gene'] = df.groupby('read_id', sort=False)['without_gene'].transform('all')

    df['gene_name_f'] = 1 / df.groupby([name_column, gene_column], sort=False).transform('size')
    gene_count = df.groupby(['read_id', gene_column], sort=False)['gene_name_f'].transform('sum')
    names_count = df.groupby('read_id', sort=False)[name_column].transform('nunique')

    df['intra'] = (gene_count == names_count)
    df['intra'] = df.groupby('read_id', sort=False)['intra'].transform('any')


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True)
parser.add_argument('-p', '--prefix', required=True)
parser.add_argument('-n', '--name_column', required=True, type=int)
parser.add_argument('-g', '--gene_column', required=True, type=int)

args = parser.parse_args()

reader = adjacent_groupby2(
    sys.stdin if args.input == 'stdin' else args.input,
    sep='\t',
    header=None,
    name_column=args.name_column,
    chunksize=CHUNKSIZE
)

category = {
    'intermolecular': lambda df: df[~df['intra']],
    'intramolecular': lambda df: df[(df['intra'] & ~df['without_gene'])],
    'without_gene': lambda df: df[df['without_gene']]
}

output = {}
for suffix in category:
    output[suffix] = (
        open(f'{args.prefix}.{suffix}.bed', 'wt', buffering=OUTPUT_BUFFER),
        category[suffix]
    )

for data in reader:
    columns = list(data.columns)
    columns.remove('read_id')
    classify_reads(data, args.name_column, args.gene_column)
    for file, func in output.values():
        func(data).to_csv(file, sep='\t', index=False, header=False, columns=columns)

for file, _ in output.values():
    file.close()
