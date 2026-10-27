"""Equivalence cases for hicexplorer.hicValidateLocations (apieq.py).

The inputs are the tool's real test data (loops_1.bedgraph, the SMC3
narrowPeak and CTCF broadPeak files, the TAD boundaries, the RAD21 cool),
used the way the tool's main uses them, plus small tables for the edge cases.
DataFrames and Series are compared by their labels, dtypes, values and index.
bedtools' input sort, which main does through pybedtools, is done here with
pandas (chromosome, then start), which is the order bedtools merge requires.

Deviations from 3.7.6 (cpp/python/src/api/hicValidateLocations.cpp):
  * DataFrame columns that are neither int, float, bool nor object of str,
    non-numeric coordinates, a bin size that is not an int, and frames whose
    labels 1 and 2 are not the columns at positions 1 and 2 raise TypeError
    naming HiCExplorer 4 (3.7.6 computed with a float bin size, for one);
  * the errors pybedtools raises when bedtools rejects its input (unsorted
    input to merge, non-integer coordinates) and the errors of a cool file
    the C++ cool reader cannot open are RuntimeError with the C++ message;
  * read_csv's ParserError carries the C++ tokenizer's line count.
"""
import io
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'test_data', 'hicValidateLocations')
LOOPS = os.path.join(DATA, 'loops_1.bedgraph')
SMC3 = os.path.join(DATA, 'GSM935376_hg19_Gm12878_Smc3.narrowPeak')
CTCF = os.path.join(DATA, 'GSM733752_hg19_ctcf_GM12878.broadPeak')
TADS = os.path.join(DATA, 'untreated_R1_boundaries.bed')
RAD21 = os.path.join(DATA, 'GSM1436265_RAD21_ENCFF002EMQ_10kb.cool')
V = 'hicexplorer.hicValidateLocations.'


def frame(df):
    return (list(df.columns), [str(d) for d in df.dtypes], [df[c].tolist() for c in df.columns],
            type(df.index).__name__, df.index.tolist())


def series(s):
    return (s.name, str(s.dtype), s.tolist(), type(s.index).__name__, s.index.tolist())


def module_of(target):
    import sys
    return sys.modules[target.__module__]


def sort_bed(df):
    return df.sort_values([df.columns[0], df.columns[1]], kind='stable').reset_index(drop=True)


def read(target, *args):
    return frame(target(*args))


def read_handle(target, path, prefix):
    with open(path) as handle:
        return frame(target(handle, prefix))


def read_text(target, text, prefix):
    return frame(target(io.StringIO(text), prefix))


def overlap_loops(target, protein_path, prefix_loops, binned):
    m = module_of(target)
    loops = m.readLoopFile(LOOPS, prefix_loops)
    protein = sort_bed(m.readThreeColumn(protein_path, None))
    if binned:
        protein = m.applyBinning(protein, 10000)
    return series(target(loops, protein))


def overlap_tads(target, columns=3):
    m = module_of(target)
    tads = m.readThreeColumn(TADS, None)
    if columns > 3:
        tads[3] = 'x'
    protein = m.applyBinning(sort_bed(m.readThreeColumn(CTCF, None)), 20000)
    return series(target(m.applyBinning(tads, 20000, False), protein))


def correlate(target, prefix, rows=None, mutate=None):
    m = module_of(target)
    loops = m.readLoopFile(LOOPS, prefix)
    if rows:
        loops = loops.head(rows)
    if mutate == 'float':
        loops[1] = loops[1].astype(float)
    elif mutate == 'narrow':
        loops = loops[[0, 1, 2, 3]]
    elif mutate == 'out_of_range':
        loops.loc[0, 2] = 10 ** 10
        loops.loc[1, 0] = 'chrNope'
    peaks, selection = target(RAD21, loops)
    return (type(peaks).__name__, int(peaks)), series(selection)


def binning(target, path, size, merge, sort=True, head=None):
    m = module_of(target)
    df = m.readThreeColumn(path, None)
    if sort:
        df = sort_bed(df)
    if head:
        df = df.head(head)
    return frame(target(df, size, merge))


def binning_frame(target, data, size, merge=True):
    import pandas as pd
    return frame(target(pd.DataFrame(data), size, merge))


def write_loops(target, directory, as_handle=False):
    m = module_of(target)
    loops = m.readLoopFile(LOOPS, 'add').head(300)
    path = os.path.join(directory, 'out.txt')
    if as_handle:
        with open(path, 'w') as handle:
            result = target(handle, loops)
    else:
        result = target(path, loops)
    with open(path) as handle:
        return result, handle.read()


SMALL = b"chr1\t100\t200\nchr1\t150\t300\nchr2\t5\t9\n"
RAGGED = b"chr1\t100\t200\nchr1\t150\t300\t4\n"
INT_CHROMS = b"1\t100\t200\n2\t150\t300\n"

CASES = [
    Case('readThreeColumn/narrowPeak', V + 'readThreeColumn', args=(SMC3, None), driver='read'),
    Case('readThreeColumn/broadPeak_add', V + 'readThreeColumn', args=(CTCF, 'add'), driver='read'),
    Case('readThreeColumn/boundaries_remove', V + 'readThreeColumn', args=(TADS, 'remove'), driver='read'),
    Case('readThreeColumn/handle', V + 'readThreeColumn', args=(TADS, None), driver='read_handle'),
    Case('readThreeColumn/text_add', V + 'readThreeColumn', args=("chr1\t1\t2\tx\n3\t4\t5\ty\n", 'add'),
         driver='read_text'),
    Case('readThreeColumn/int_chroms_add', V + 'readThreeColumn', args=(File('a.bed'), 'add'),
         files={'a.bed': INT_CHROMS}, driver='read'),
    Case('readThreeColumn/int_chroms_remove', V + 'readThreeColumn', args=(File('a.bed'), 'remove'),
         files={'a.bed': INT_CHROMS}, driver='read'),
    Case('readThreeColumn/two_columns', V + 'readThreeColumn', args=(File('a.bed'), None),
         files={'a.bed': b"chr1\t1\n"}, driver='read'),
    Case('readThreeColumn/empty', V + 'readThreeColumn', args=(File('a.bed'), None), files={'a.bed': b""},
         driver='read'),
    Case('readThreeColumn/missing', V + 'readThreeColumn', args=(File('missing.bed'), None), driver='read'),

    Case('readLoopFile/plain', V + 'readLoopFile', args=(LOOPS, None), driver='read'),
    Case('readLoopFile/add', V + 'readLoopFile', args=(LOOPS, 'add'), driver='read'),
    Case('readLoopFile/remove', V + 'readLoopFile', args=(File('l.txt'), 'remove'),
         files={'l.txt': b"chr1\t1\t2\tchrX\t3\t4\nhr2\t1\t2\tcrr\t3\t4\n"}, driver='read'),
    Case('readLoopFile/three_columns_add', V + 'readLoopFile', args=(File('a.bed'), 'add'), files={'a.bed': SMALL},
         driver='read'),
    Case('readLoopFile/unknown_prefix_mode', V + 'readLoopFile', args=(File('a.bed'), 'None'),
         files={'a.bed': SMALL}, driver='read'),

    Case('overlapLoop/smc3_binned', V + 'overlapLoop', args=(SMC3, 'add', True), driver='overlap_loops'),
    Case('overlapLoop/ctcf_unbinned', V + 'overlapLoop', args=(CTCF, 'add', False), driver='overlap_loops'),
    Case('overlapLoop/no_prefix', V + 'overlapLoop', args=(SMC3, None, True), driver='overlap_loops'),

    Case('overlapTAD/ctcf', V + 'overlapTAD', driver='overlap_tads'),
    Case('overlapTAD/text_name_column', V + 'overlapTAD', args=(4,), driver='overlap_tads'),

    Case('correlateCool/rad21', V + 'correlateCool', args=('add',), driver='correlate'),
    Case('correlateCool/no_prefix', V + 'correlateCool', args=(None, 400), driver='correlate'),
    Case('correlateCool/float_coordinates', V + 'correlateCool', args=('add', 400, 'float'), driver='correlate'),
    Case('correlateCool/four_columns', V + 'correlateCool', args=('add', 50, 'narrow'), driver='correlate'),
    Case('correlateCool/out_of_range', V + 'correlateCool', args=('add', 50, 'out_of_range'), driver='correlate'),

    Case('applyBinning/merge', V + 'applyBinning', args=(SMC3, 10000, True), driver='binning'),
    Case('applyBinning/no_merge', V + 'applyBinning', args=(CTCF, 20000, False, False), driver='binning'),
    Case('applyBinning/no_merge_no_duplicates', V + 'applyBinning', args=(TADS, 1, False, False, 30),
         driver='binning'),
    Case('applyBinning/merge_default', V + 'applyBinning', args=(TADS, 20000, 1), driver='binning'),
    Case('applyBinning/extra_columns', V + 'applyBinning',
         args=({0: ['a', 'a', 'b', 'a'], 1: [5, 6, 1, 15000], 2: [7, 8, 2, 18000], 3: [1.5, 2.5, 3.5, 4.5]},
               10000, False), driver='binning_frame'),
    Case('applyBinning/floats', V + 'applyBinning',
         args=({0: ['a', 'a', 'b'], 1: [5.5, 19999.9, 1.0], 2: [7.0, 25000.5, 2.0]}, 10000), driver='binning_frame'),
    Case('applyBinning/text_column', V + 'applyBinning', args=({0: ['a'], 1: ['x'], 2: [3]}, 10),
         driver='binning_frame'),

    Case('writeLoopFile/path', V + 'writeLoopFile', args=(File(),), driver='write_loops'),
    Case('writeLoopFile/handle', V + 'writeLoopFile', args=(File(), True), driver='write_loops'),
]
