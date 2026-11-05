"""Equivalence cases for hicexplorer.hicMergeLoops (apieq.py).

mergeLoops gets what the 3.7.6 main() hands it: the loop files of the
hicMergeLoops tests (gm12878 at 5, 10 and 25 kb) read with pandas,
concatenated, sorted by chromosome and start (a stable pandas sort in place
of pybedtools, so that the cases need no bedtools), duplicates dropped, and
the x and y anchors as tuples. Constructed frames cover a chromosome that
reappears in the tuples, a first loop on an unknown chromosome, the y tree
looked up under the x chromosome, and stale overlap sets. readFile is
compared as the DataFrame's columns, dtypes and values (driver 'frame').
"""
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data', 'hicMergeLoops')
M = 'hicexplorer.hicMergeLoops.'


def read(name):
    import pandas as pd
    return pd.read_csv(os.path.join(TEST_DATA, name), sep='\t', header=None)


def prepared(names, resolution, limit=None):
    import pandas as pd

    def make():
        frame = pd.concat([read(n) for n in names])
        frame = frame.sort_values([0, 1], kind='mergesort').reset_index(drop=True)
        frame = frame.drop_duplicates(keep=False)
        if limit is not None:
            frame = frame.iloc[:limit]
        tuples_x = [tuple(x) for x in frame[[0, 1, 2]].values]
        tuples_y = [tuple(x) for x in frame[[3, 4, 5]].values]
        return (frame, resolution, tuples_x, tuples_y)
    return make


def constructed(rows, tuples_x=None, tuples_y=None, resolution=10000):
    def make():
        import pandas as pd
        frame = pd.DataFrame(rows)
        tx = tuples_x if tuples_x is not None else [tuple(x) for x in frame[[0, 1, 2]].values]
        ty = tuples_y if tuples_y is not None else [tuple(x) for x in frame[[3, 4, 5]].values]
        return (frame, resolution, tx, ty)
    return make


ROWS = [
    ('1', 10000, 20000, '1', 250000, 260000, 0.001),
    ('1', 10000, 15000, '1', 250000, 255000, 0.002),
    ('1', 12000, 17000, '1', 252000, 257000, 0.003),
    ('2', 50000, 60000, '2', 90000, 100000, 0.01),
    ('2', 52000, 57000, '2', 92000, 97000, 0.01),
]

REAPPEARING = ROWS + [('1', 11000, 16000, '1', 251000, 256000, 0.004)]

INTER = [
    ('1', 10000, 20000, '1', 250000, 260000, 0.001),
    ('1', 10000, 15000, '1', 250000, 255000, 0.002),
    ('3', 10000, 20000, '1', 250000, 260000, 0.001),
]

STALE = [
    ('1', 10000, 20000, '1', 250000, 260000, 0.001),
    ('1', 10000, 15000, '1', 250000, 255000, 0.002),
    ('1', 10000, 15000, 'Q', 250000, 255000, 0.002),
]


def frame(target, *args):
    df = target(*args)
    return (type(df).__name__, [type(c).__name__ for c in df.columns], list(df.columns),
            [str(df[c].dtype) for c in df.columns], [df[c].values for c in df.columns],
            type(df.index).__name__, len(df.index))


TSV = (b'1\t10000\t20000\t1\t250000\t260000\t0.001\n'
       b'X\t5\t6\tX\t7\t8\t1e-5\n'
       b'\n'
       b'2\t50000\t60000\t2\t90000\t100000\tNA\n')
STRINGS = b'chr1\t10\t20\tfoo\nchr2\t30\t40\t\nchr3\t50\t60\tbar baz\n'
FLOATS = b'1\t0.5\t2\n2\t1.25\t3\n3\t0.0019430592210407135\t4.0\n'

CASES = [
    # mergeLoops -------------------------------------------------------------
    Case('merge/all_resolutions', M + 'mergeLoops',
         args=prepared(['gm12878_5kb.bedgraph', 'gm12878_10kb.bedgraph', 'gm12878_25kb.bedgraph'], 25000)),
    Case('merge/10kb', M + 'mergeLoops', args=prepared(['gm12878_10kb.bedgraph'], 10000)),
    Case('merge/5kb_lowest_25kb', M + 'mergeLoops', args=prepared(['gm12878_5kb.bedgraph'], 25000)),
    Case('merge/merged_file', M + 'mergeLoops', args=prepared(['gm12878_all.bedgraph'], 25000, limit=3000)),
    Case('merge/resolution_str', M + 'mergeLoops', args=prepared(['gm12878_25kb.bedgraph'], '25000', limit=500)),
    Case('merge/constructed', M + 'mergeLoops', args=constructed(ROWS)),
    Case('merge/reappearing_chromosome', M + 'mergeLoops', args=constructed(REAPPEARING)),
    Case('merge/inter_chromosomal_key_error', M + 'mergeLoops', args=constructed(INTER)),
    Case('merge/stale_y_set', M + 'mergeLoops', args=constructed(STALE)),
    Case('merge/first_loop_unknown_chromosome', M + 'mergeLoops',
         args=constructed(ROWS, tuples_x=[('9', 0, 1)], tuples_y=[('9', 0, 1)])),
    Case('merge/y_unknown', M + 'mergeLoops',
         args=constructed(ROWS, tuples_y=[('9', 0, 1)])),
    Case('merge/empty_tuples', M + 'mergeLoops', args=constructed(ROWS, tuples_x=[], tuples_y=[])),
    Case('merge/too_few_columns', M + 'mergeLoops',
         args=lambda: (constructed(ROWS)()[0][[0, 1, 2, 3]], 10000, [('1', 0, 1)], [('1', 0, 1)])),
    # readFile ---------------------------------------------------------------
    Case('read/real_10kb', M + 'readFile', args=(os.path.join(TEST_DATA, 'gm12878_10kb.bedgraph'),), driver='frame'),
    Case('read/real_all', M + 'readFile', args=(os.path.join(TEST_DATA, 'gm12878_all.bedgraph'),), driver='frame'),
    Case('read/mixed_and_na', M + 'readFile', args=(File('a.tsv'),), files={'a.tsv': TSV}, driver='frame'),
    Case('read/strings', M + 'readFile', args=(File('a.tsv'),), files={'a.tsv': STRINGS}, driver='frame'),
    Case('read/floats', M + 'readFile', args=(File('a.tsv'),), files={'a.tsv': FLOATS}, driver='frame'),
    Case('read/empty', M + 'readFile', args=(File('a.tsv'),), files={'a.tsv': b''}, driver='frame'),
    Case('read/missing', M + 'readFile', args=(File('missing.tsv'),), driver='frame'),
]
