#!/usr/bin/env python3
from __future__ import annotations
import argparse
from dataclasses import dataclass
from multiprocessing import Pool
import re
import sys
import numpy as np
import pandas as pd


OUTPUT_BUFFER = 536870912


SAM_COLUMNS = [
    'QNAME',
    'FLAG',
    'RNAME',
    'POS',
    'MAPQ',
    'CIGAR',
    'RNEXT',
    'PNEXT',
    'TLEN',
    'SEQ',
    'QUAL',
    'jM'
]

JUNCTION_CODES = {
    'NONE': -1,
    'CHIMERIC': 30,
    'MATES': 40,
    'AMBIGUOUS': 50
}

OFFSET_PATTERNS = {
    '+': re.compile(r'^(\d+)S'),
    '-': re.compile(r'(\d+)S$')
}

STRANDS = ('+', '-')
MATES = (1, 2)

# cigar operations
OPERATIONS = {}
for cigar_operation in {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}:
    OPERATIONS[cigar_operation] = {
        'consumes_query': cigar_operation in {'M', 'I', 'S', '=', 'X'},
        'consumes_reference': cigar_operation in {'M', 'D', 'N', '=', 'X'},
    }
OPERATION_PATTERN = re.compile(r'(\d+)(\D)')


@dataclass
class Tag:
    chr: str = '.'
    start: int = 0
    end: int = 0
    strand: str = '.'
    junction: int = JUNCTION_CODES['NONE']
    cigar: str = '.'
    score: int = 0

    def __post_init__(self):
        assert self.start < self.end or (self.start == self.end == 0)
        assert self.strand in STRANDS or self.strand == '.'

    def intersection(self, other: Tag) -> bool:
        if self.chr != other.chr or self.strand != other.strand:
            return False
        elif (self.end <= other.start) or (other.end <= self.start):
            return False
        return True

    def merge(self, other: Tag) -> Tag:
        if not self.intersection(other):
            return Tag()
        if self.junction == other.junction:
            junction = self.junction
        elif self.junction == JUNCTION_CODES['NONE']:
            junction = other.junction
        elif other.junction == JUNCTION_CODES['NONE']:
            junction = self.junction
        else:
            junction = JUNCTION_CODES['AMBIGUOUS']
        equality = (self.start == other.start) and (self.end == other.end)
        result = Tag(
            chr=self.chr,
            strand=self.strand,
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            cigar=f'{self.cigar}/{other.cigar}',
            score=1 if equality else 0.5,
            junction=junction
        )
        return result

    def __bool__(self) -> bool:
        return self.chr != '.'


def collapse_series(series: pd.Series) -> any:
    assert series.nunique() == 1, f'Value is not unique: {series}'
    return series.iloc[0]


def mate(sam_flag: int) -> int:
    forward, reverse = bool(sam_flag % 128 // 64), bool(sam_flag % 256 // 128)
    assert forward ^ reverse, 'Either 0x40 or 0x80 flag must be set for read!'
    return 1 if forward else 2


def strand(sam_flag: int) -> str:
    return '-' if sam_flag % 32 // 16 else '+'


def cigar2offset(sam_cigar: str, sam_strand: str) -> int:
    offset = re.findall(OFFSET_PATTERNS[sam_strand], sam_cigar)
    assert len(offset) <= 1, f'Incorrect fragment cigar: {sam_cigar}'
    return int(offset[0]) if offset else 0


def tags_split(fragment_data: pd.DataFrame) -> list[Tag]:
    start = pos = fragment_data['POS'] - 1  # to 0-based
    junctions = map(int, fragment_data['jM'].removeprefix('jM:B:c,').split(','))
    result = []
    current_cigar = []
    for n, operation in re.findall(OPERATION_PATTERN, fragment_data['CIGAR']):
        n = int(n)
        if operation == 'N':
            result.append(Tag(
                chr=fragment_data['RNAME'],
                start=start,
                end=pos,
                strand=fragment_data['strand'],
                cigar=''.join(current_cigar),
                junction=next(junctions)
            ))
            start = pos + n
            current_cigar = []
        else:
            current_cigar.append(f'{n}{operation}')
        if OPERATIONS[operation]['consumes_reference']:
            pos += n
    if current_cigar:
        result.append(Tag(
            chr=fragment_data['RNAME'],
            start=start,
            end=pos,
            strand=fragment_data['strand'],
            cigar=''.join(current_cigar),
            junction=JUNCTION_CODES['NONE']
        ))
    return result


def reverse_tags(tags: list[Tag]) -> None:
    tags.reverse()
    for i in range(len(tags)-1):
        tags[i].junction = tags[i+1].junction
    tags[-1].junction = JUNCTION_CODES['NONE']


def assemble_mate(mate_data: pd.DataFrame) -> list[Tag]:
    mate_n = collapse_series(mate_data['mate'])
    mate_tags = []
    for _, fragment_data in mate_data.sort_values('offset').iterrows():
        assert fragment_data['strand'] in STRANDS
        fragment_tags = fragment_data['tags'][:]
        if fragment_data['strand'] == '-':
            reverse_tags(fragment_tags)
        if mate_tags:
            last_tag = mate_tags[-1]
            assert last_tag.junction == JUNCTION_CODES['NONE']
            last_tag.junction = JUNCTION_CODES['CHIMERIC']
        mate_tags.extend(fragment_tags)
    assert mate_n in MATES
    if mate_n == 2:
        reverse_tags(mate_tags)
    return mate_tags


def align_mates(mate1: list[Tag], mate2: list[Tag]) -> list[Tag]:
    for tag in mate2:
        tag.strand = '+' if tag.strand == '-' else '-'
    matrix = []
    for tag1 in mate1:
        matrix.append([tag1.intersection(tag2) for tag2 in mate2])
    matrix = np.array(matrix, dtype=np.ushort)
    i_max, j_max = matrix.shape
    for i in range(1, i_max):
        for j in range(1, j_max):
            matrix[i][j] += matrix[i - 1][j - 1]
    i, j = np.unravel_index(matrix.argmax(), matrix.shape)
    not_aligned = (matrix[i][j] == 0)
    if not_aligned:
        assert mate1[-1].junction == JUNCTION_CODES['NONE']
        mate1[-1].junction = JUNCTION_CODES['MATES']
        return mate1 + mate2

    result = []
    i, j = i - max(i, j), j - max(i, j)
    while i < len(mate1) or j < len(mate2):
        tag1 = mate1[i] if 0 <= i < len(mate1) else {}
        tag2 = mate2[j] if 0 <= j < len(mate2) else {}
        if tag1 and tag2:
            result.append(tag1.merge(tag2))
        else:
            result.append(tag1 or tag2)
        i, j = i+1, j+1
    return result


def assemble_read(read_data: pd.DataFrame) -> list[Tag]:
    result = {}
    for mate_n, mate_data in read_data.groupby('mate'):
        result[mate_n] = assemble_mate(mate_data)
    if len(result) < 2:
        _, result = result.popitem()
    else:
        result = align_mates(result[1], result[2])
    return result


def chunk_processing(chunk: pd.DataFrame) -> None:
    chunk['mate'] = chunk['FLAG'].apply(mate)
    chunk['strand'] = chunk['FLAG'].apply(strand)
    chunk['offset'] = chunk.apply(
        lambda df: cigar2offset(df['CIGAR'], df['strand']),
        axis=1
    )
    chunk['tags'] = chunk.apply(tags_split, axis=1)


def adjacent_groupby(*args, groupby_column: str, processing: callable, **kwargs):
    if 'chunksize' not in kwargs:
        kwargs['chunksize'] = 700
    with pd.read_csv(*args, **kwargs) as chunk_iterator:
        chunk = next(chunk_iterator)
        processing(chunk)
        *mid, remaining = chunk.groupby(groupby_column, sort=False)
        yield from mid
        for chunk in chunk_iterator:
            processing(chunk)
            try:
                first, *mid, last = chunk.groupby(groupby_column, sort=False)
            except ValueError:
                first, *mid = chunk.groupby(groupby_column, sort=False)
                last = None
            if first[0] == remaining[0]:
                yield first[0], pd.concat([remaining[1], first[1]])
            else:
                yield remaining
                yield first
            yield from mid
            remaining = last
        if remaining is not None:
            yield remaining


def pmap(func, iterable, threads: int):
    with Pool(threads) as pool:
        yield from pool.imap(func, iterable)


def get_tags_names(input_filename: str) -> list[str]:
    with open(input_filename, 'rt') as input_file:
        for line in input_file:
            if line[0] == '@':
                continue
            tags = line.split('\t')[len(SAM_COLUMNS):]
            return [s.split(':')[0] for s in tags]


def read_data2bed_lines(read_data: tuple[str, pd.DataFrame]):
    assembly = assemble_read(read_data[1])
    ambiguous = False
    result = []
    for n, tag in enumerate(assembly):
        bed_line = [tag.chr, tag.start, tag.end, f'{read_data[0]}/{n+1}',
                    tag.score, tag.strand, tag.junction, tag.cigar]
        if not tag:
            ambiguous = True
        result.append(bed_line)
    return ambiguous, result


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-a', '--ambiguous', required=True)
    parser.add_argument('-T', '--threads', default=8, type=int)

    args = parser.parse_args()

    reader = adjacent_groupby(
        sys.stdin if args.input == 'stdin' else args.input,
        names=SAM_COLUMNS,
        sep='\t',
        groupby_column='QNAME',
        processing=chunk_processing,
        comment='@'
    )
    with (open(args.output, mode='wt', buffering=OUTPUT_BUFFER) as output,
          open(args.ambiguous, mode='wt', buffering=OUTPUT_BUFFER) as ambiguous):
        for flag, lines in pmap(read_data2bed_lines, reader, args.threads):
            output_file = ambiguous if flag else output
            for line in lines:
                print(*line, sep='\t', file=output_file)


if __name__ == '__main__':
    main()
