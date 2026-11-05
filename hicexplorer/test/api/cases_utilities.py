"""Equivalence cases for hicexplorer.utilities (apieq.py).

The matrices are real chromosome blocks of the repository's test data
(small_test_matrix_50kb_res.h5, read with PyTables so that both
interpreters build the same input), plus small constructed ones for the edge cases the
real data does not reach: explicit zeros, non-square shapes, other dtypes.
"""
import gzip
import os
import pathlib

import numpy as np
import scipy.sparse as sp

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
U = 'hicexplorer.utilities.'


# ---------------------------------------------------------------------------
# inputs

def real_block(chromosome_index=0, dtype=np.float64):
    """The upper triangle block of one chromosome of small_test_matrix_50kb_res.h5,
    symmetrised as hiCMatrix does, as a csr_matrix."""
    import tables
    with tables.open_file(os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'), 'r') as f:
        m = sp.csr_matrix((f.root.matrix.data.read(), f.root.matrix.indices.read(), f.root.matrix.indptr.read()),
                          shape=tuple(f.root.matrix.shape.read()))
        chroms = [c.decode() if isinstance(c, bytes) else c for c in f.root.intervals.chr_list.read()]
    names = []
    for c in chroms:
        if c not in names:
            names.append(c)
    first = chroms.index(names[chromosome_index])
    last = len(chroms) - chroms[::-1].index(names[chromosome_index])
    block = m[first:last, first:last]
    block = block + sp.triu(block, 1).T
    block = sp.csr_matrix(block)
    if dtype in (np.int32, np.int64):
        block.data = np.rint(block.data).astype(dtype)
    else:
        block.data = block.data.astype(dtype)
    return block


def real_cut_intervals(limit=60):
    import tables
    with tables.open_file(os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'), 'r') as f:
        chroms = [c.decode() if isinstance(c, bytes) else c for c in f.root.intervals.chr_list.read()]
        starts = f.root.intervals.start_list.read()
        ends = f.root.intervals.end_list.read()
    rows = [(chroms[i], starts[i], ends[i], 1.0) for i in range(len(chroms))]
    # drop every fifth bin, as masking does, so that the gaps must be closed
    kept = [row for i, row in enumerate(rows) if i % 5 != 2]
    return kept[:limit] + kept[-limit:]


def explicit_zero_matrix(dtype=np.float64):
    data = np.array([1.0, 0.0, 2.0, 3.0, 0.0, 4.0, 5.0], dtype=dtype)
    indices = np.array([0, 1, 2, 1, 2, 2, 3], dtype=np.int32)
    indptr = np.array([0, 3, 5, 6, 7], dtype=np.int32)
    return sp.csr_matrix((data, indices, indptr), shape=(4, 4))


def small_matrix(dtype=np.float64):
    dense = np.array([[4, 2, 0, 1],
                      [2, 6, 3, 0],
                      [0, 3, 5, 2],
                      [1, 0, 2, 7]], dtype=dtype)
    return sp.csr_matrix(dense)


def nan_inf_matrix(dtype=np.float64):
    data = np.array([1.0, np.nan, np.inf, -np.inf, 2.5, np.nan], dtype=dtype)
    return sp.csr_matrix((data, np.array([0, 1, 2, 0, 1, 2]), np.array([0, 3, 6])), shape=(2, 3))


def rng_values(n, dtype, seed=7, outliers=True):
    values = np.random.default_rng(seed).normal(10.0, 2.0, n)
    if outliers:
        values[::17] *= 9.0
    if dtype in (np.int32, np.int64, np.uint8):
        return np.rint(values).astype(dtype)
    return values.astype(dtype)


# ---------------------------------------------------------------------------
# drivers

def queue_result(target, *args):
    """expected_interactions_thread puts its result into a queue."""
    import queue
    q = queue.Queue()
    target(*args, q)
    result = q.get()
    if isinstance(result, str):
        # the traceback text of a failed worker is not part of the contract
        return result.split('Traceback')[0]
    return result


def thread_message(target, *args, **kwargs):
    result = target(*args, **kwargs)
    if isinstance(result, str):
        return result.split('Traceback')[0]
    return result


def read_file(target, *args):
    """opener returns a file object: compare its type and its content."""
    handle = target(*args)
    try:
        return type(handle).__name__, handle.read()
    finally:
        handle.close()


def call_and_list_state(target, *args):
    """For functions that change a list in place and return it."""
    result = target(*args)
    return result, args[0]


BED = (b"# a comment\n"
       b"chr1\t100\t200\tgeneA\n"
       b"chr1\t300\t400\n"
       b"chr2 500 600 geneB\n"
       b"chrX\t1\t2\t3\t4\t5\n"
       b"#chr3\t1\t2\n")


CASES = [
    # readBed -----------------------------------------------------------------
    Case('readBed/real', U + 'readBed', args=(os.path.join(TEST_DATA, 'cHi-C', 'referencePoints.bed'),)),
    Case('readBed/real_domains', U + 'readBed', args=(os.path.join(TEST_DATA, 'domains.bed'),)),
    Case('readBed/mixed', U + 'readBed', args=(File('a.bed'),), files={'a.bed': BED}),
    Case('readBed/crlf_and_cr', U + 'readBed', args=(File('a.bed'),),
         files={'a.bed': b"chr1\t1\t2\r\nchr2\t3\t4\rchr3\t5\t6"}),
    Case('readBed/unicode', U + 'readBed', args=(File('a.bed'),),
         files={'a.bed': "﻿chré1\t1\t2\nchr北　 7 8\n".encode('utf-8')}),
    Case('readBed/spaces_only', U + 'readBed', args=(File('a.bed'),), files={'a.bed': b"chr1 1 2 x\nchr1  3\t 4\n"}),
    Case('readBed/blank_line', U + 'readBed', args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\n\nchr1\t3\t4\n"}),
    Case('readBed/two_fields', U + 'readBed', args=(File('a.bed'),), files={'a.bed': b"chr1\t1\n"}),
    Case('readBed/empty', U + 'readBed', args=(File('a.bed'),), files={'a.bed': b""}),
    Case('readBed/invalid_utf8', U + 'readBed', args=(File('a.bed'),), files={'a.bed': b"chr1\t1\t2\nchr\xff\t3\t4\n"}),
    Case('readBed/invalid_utf8_second_chunk', U + 'readBed', args=(File('a.bed'),),
         files={'a.bed': b"chr1\t1\t2\n" * 1000 + b"chr\xe2\x82\t3\t4\n"}),
    Case('readBed/truncated_utf8_at_end', U + 'readBed', args=(File('a.bed'),),
         files={'a.bed': b"chr1\t1\t2\nchr\xe2\x82"}),
    Case('readBed/utf8_across_chunk', U + 'readBed', args=(File('a.bed'),),
         files={'a.bed': b"c" * 8191 + "€\t1\t2\n".encode('utf-8')}),
    Case('readBed/missing', U + 'readBed', args=(File('missing.bed'),)),
    Case('readBed/directory', U + 'readBed', args=(File(),)),
    Case('readBed/pathlib', U + 'readBed', args=lambda: (pathlib.Path(TEST_DATA) / 'regions.bed',)),
    Case('readBed/bytes_path', U + 'readBed', args=(os.path.join(TEST_DATA, 'regions.bed').encode(),)),
    Case('readBed/embedded_null', U + 'readBed', args=('a\0b',)),
    Case('readBed/not_a_path', U + 'readBed', args=(3.5,)),

    # readTargetBed -----------------------------------------------------------
    Case('readTargetBed/real_4col', U + 'readTargetBed',
         args=(os.path.join(TEST_DATA, 'cHi-C', 'chicAggregateStatistic', 'target_list_4col.bed'),)),
    Case('readTargetBed/real_3col', U + 'readTargetBed',
         args=(os.path.join(TEST_DATA, 'cHi-C', 'chicAggregateStatistic', 'target_list_3col.bed'),)),
    Case('readTargetBed/mixed', U + 'readTargetBed', args=(File('t.bed'),),
         files={'t.bed': b"#c\nchr1\t1\t2\tg1\nchr1\t3\t4\nchr2 5 6 g2 extra\nchr3\t7\t8\tg1\n"}),
    Case('readTargetBed/three_first', U + 'readTargetBed', args=(File('t.bed'),), files={'t.bed': b"chr1\t1\t2\n"}),
    Case('readTargetBed/too_few', U + 'readTargetBed', args=(File('t.bed'),), files={'t.bed': b"chr1\t1\t2\tg\nchr1 2\n"}),
    Case('readTargetBed/missing', U + 'readTargetBed', args=(File('missing.bed'),)),

    # writableFile ------------------------------------------------------------
    Case('writableFile/new', U + 'writableFile', args=(File('out.txt'),)),
    Case('writableFile/existing', U + 'writableFile', args=(File('out.txt'),), files={'out.txt': b'content'}),
    Case('writableFile/missing_directory', U + 'writableFile', args=(File('no/such/out.txt'),)),
    Case('writableFile/directory', U + 'writableFile', args=(File(),)),
    Case('writableFile/embedded_null', U + 'writableFile', args=('a\0b',)),

    # remove_outliers ---------------------------------------------------------
    Case('remove_outliers/float64', U + 'remove_outliers', args=lambda: (rng_values(200, np.float64),)),
    Case('remove_outliers/float32', U + 'remove_outliers', args=lambda: (rng_values(201, np.float32),)),
    Case('remove_outliers/int64', U + 'remove_outliers', args=lambda: (rng_values(150, np.int64),)),
    Case('remove_outliers/int32', U + 'remove_outliers', args=lambda: (rng_values(151, np.int32),)),
    Case('remove_outliers/uint8', U + 'remove_outliers', args=lambda: (rng_values(60, np.uint8),)),
    Case('remove_outliers/bool', U + 'remove_outliers', args=lambda: (np.array([True, False, True, True]),)),
    Case('remove_outliers/max_deviation', U + 'remove_outliers', args=lambda: (rng_values(200, np.float64), 1.3)),
    Case('remove_outliers/max_deviation_float32', U + 'remove_outliers',
         args=lambda: (rng_values(200, np.float32), 1.3)),
    Case('remove_outliers/2d', U + 'remove_outliers', args=lambda: (rng_values(60, np.float64).reshape(6, 10),)),
    Case('remove_outliers/strided', U + 'remove_outliers', args=lambda: (rng_values(200, np.float64)[::3],)),
    Case('remove_outliers/nan', U + 'remove_outliers', args=lambda: (np.array([1.0, np.nan, 3.0, 50.0]),)),
    Case('remove_outliers/constant', U + 'remove_outliers', args=lambda: (np.ones(9),)),
    Case('remove_outliers/empty', U + 'remove_outliers', args=lambda: (np.array([]),)),
    Case('remove_outliers/list_mad_zero', U + 'remove_outliers', args=lambda: ([1.0, 1.0, 1.0],)),
    Case('remove_outliers/list', U + 'remove_outliers', args=lambda: ([1.0, 2.0, 3.0, 400.0],)),

    # convertNansToZeros / convertInfsToZeros ---------------------------------
    Case('convertNansToZeros/float64', U + 'convertNansToZeros', args=lambda: (nan_inf_matrix(),)),
    Case('convertNansToZeros/float32', U + 'convertNansToZeros', args=lambda: (nan_inf_matrix(np.float32),)),
    Case('convertNansToZeros/int', U + 'convertNansToZeros', args=lambda: (small_matrix(np.int32),)),
    Case('convertNansToZeros/empty', U + 'convertNansToZeros', args=lambda: (sp.csr_matrix((3, 3)),)),
    Case('convertNansToZeros/real', U + 'convertNansToZeros', args=lambda: (real_block(),)),
    Case('convertInfsToZeros/float64', U + 'convertInfsToZeros', args=lambda: (nan_inf_matrix(),)),
    Case('convertInfsToZeros/float32', U + 'convertInfsToZeros', args=lambda: (nan_inf_matrix(np.float32),)),
    Case('convertInfsToZeros/empty', U + 'convertInfsToZeros', args=lambda: (sp.csr_matrix((2, 2)),)),

    # convertInfsToZeros_ArrayFloat / convertNansToOnes -----------------------
    Case('convertInfsToZeros_ArrayFloat/zero', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([1.0, np.nan, np.inf, -np.inf, 2.0]),)),
    Case('convertInfsToZeros_ArrayFloat/epsilon', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([1.0, np.nan, np.inf, -np.inf, 2.0]), True)),
    Case('convertInfsToZeros_ArrayFloat/float32_epsilon', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([np.nan, 1.0, np.inf], dtype=np.float32), 1)),
    Case('convertInfsToZeros_ArrayFloat/float16', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([np.nan, 1.0, np.inf], dtype=np.float16), True)),
    Case('convertInfsToZeros_ArrayFloat/int', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.arange(5),)),
    Case('convertInfsToZeros_ArrayFloat/empty', U + 'convertInfsToZeros_ArrayFloat', args=lambda: (np.array([]),)),
    Case('convertInfsToZeros_ArrayFloat/2d_rows', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([[np.nan, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]),),
         note='numpy indexes the first axis with the flat positions'),
    Case('convertInfsToZeros_ArrayFloat/2d_out_of_bounds', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array([[1.0, 2.0, 3.0], [4.0, np.nan, 5.0]]),)),
    Case('convertInfsToZeros_ArrayFloat/list_clean', U + 'convertInfsToZeros_ArrayFloat', args=lambda: ([1.0, 2.0],)),
    Case('convertInfsToZeros_ArrayFloat/list_nan', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: ([1.0, float('nan')],)),
    Case('convertInfsToZeros_ArrayFloat/object', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.array(['a', 'b'], dtype=object),)),
    Case('convertInfsToZeros_ArrayFloat/read_only', U + 'convertInfsToZeros_ArrayFloat',
         args=lambda: (np.frombuffer(np.array([1.0, np.nan]).tobytes()),)),
    Case('convertNansToOnes/float64', U + 'convertNansToOnes', args=lambda: (np.array([np.nan, 2.0, np.inf, np.nan]),)),
    Case('convertNansToOnes/complex', U + 'convertNansToOnes', args=lambda: (np.array([complex(np.nan, 1), 2 + 0j]),)),
    Case('convertNansToOnes/empty', U + 'convertNansToOnes', args=lambda: (np.array([], dtype=np.float32),)),

    # myAverage ---------------------------------------------------------------
    Case('myAverage/mean', U + 'myAverage', args=lambda: (np.array([1.0, np.nan, 3.0, 4.5]),)),
    Case('myAverage/median', U + 'myAverage', args=lambda: (np.array([1.0, np.nan, 3.0, 4.5]), 'median')),
    Case('myAverage/median_odd', U + 'myAverage', args=lambda: (rng_values(101, np.float64), 'anything')),
    Case('myAverage/mean_float32', U + 'myAverage', args=lambda: (rng_values(10001, np.float32),)),
    Case('myAverage/median_float32', U + 'myAverage', args=lambda: (rng_values(1000, np.float32), 'median')),
    Case('myAverage/mean_int', U + 'myAverage', args=lambda: (rng_values(20000, np.int64),)),
    Case('myAverage/median_int', U + 'myAverage', args=lambda: (rng_values(20, np.int32), 'median')),
    Case('myAverage/mean_large', U + 'myAverage', args=lambda: (rng_values(30000, np.float64),)),
    Case('myAverage/all_nan', U + 'myAverage', args=lambda: (np.array([np.nan, np.nan]),)),
    Case('myAverage/2d', U + 'myAverage', args=lambda: (rng_values(30, np.float64).reshape(5, 6), 'median')),
    Case('myAverage/list', U + 'myAverage', args=lambda: ([1.0, 2.0],)),

    # enlarge_bins ------------------------------------------------------------
    Case('enlarge_bins/doctest', U + 'enlarge_bins',
         args=lambda: ([('chr1', 10, 50, 1), ('chr1', 50, 80, 2), ('chr2', 10, 60, 3), ('chr2', 70, 90, 4)],)),
    Case('enlarge_bins/real', U + 'enlarge_bins', args=lambda: (real_cut_intervals(),)),
    Case('enlarge_bins/odd_gaps_and_overlaps', U + 'enlarge_bins',
         args=lambda: ([('c', 5, 10, 'a'), ('c', 13, 20, 'b'), ('c', 18, 30, 'c'), ('d', 3, 4, 'd'),
                        ('e', 7, 9, None)],)),
    Case('enlarge_bins/single', U + 'enlarge_bins', args=lambda: ([('chr1', 10, 50, 1)],)),
    Case('enlarge_bins/numpy_types', U + 'enlarge_bins',
         args=lambda: ([('chr1', np.int64(10), np.int64(50), 1), ('chr1', np.int64(55), np.int64(80), 2),
                        ('chr1', 81, 90, 3), ('chr1', np.int32(95), 100, 4)],)),
    Case('enlarge_bins/lists', U + 'enlarge_bins', args=lambda: ([['chr1', 10, 50, 1], ['chr1', 60, 80, 2]],)),
    Case('enlarge_bins/empty', U + 'enlarge_bins', args=lambda: ([],)),
    Case('enlarge_bins/short_tuple', U + 'enlarge_bins', args=lambda: ([('chr1', 10, 50)],)),
    Case('enlarge_bins/tuple_input', U + 'enlarge_bins', args=lambda: ((('chr1', 10, 50, 1), ('chr1', 60, 80, 2)),)),

    # genomicRegion -----------------------------------------------------------
] + [Case('genomicRegion/{}'.format(i), U + 'genomicRegion', args=(text,)) for i, text in enumerate([
    'chr1:1,000-2,000', ' chr1 : 10 - 20 ', '', '   \t', ',,;|!{}()', 'chrX', '　chr1:5 7',
    'chr1-1-2', '{chr1:(10)}', 'chré:1-2', '\ud800x'])] + [
    Case('genomicRegion/bytes_empty', U + 'genomicRegion', args=(b'',)),
    Case('genomicRegion/bytes', U + 'genomicRegion', args=(b'chr1:1-2',)),
    Case('genomicRegion/int', U + 'genomicRegion', args=(5,)),

    # getUserRegion -----------------------------------------------------------
    Case('getUserRegion/doctest1', U + 'getUserRegion', args=({'chr2': 1000}, 'chr2:10:1001')),
    Case('getUserRegion/doctest2', U + 'getUserRegion', args=({'chr2': 200000}, 'chr2:10:123344:3')),
    Case('getUserRegion/chrom_only', U + 'getUserRegion', args=({'chr2': 5000000}, 'chr2')),
    Case('getUserRegion/start_only', U + 'getUserRegion', args=({'chr2': 5000}, 'chr2:300')),
    Case('getUserRegion/unknown', U + 'getUserRegion', args=({'chr2': 1000, 'chrX': 3}, 'chr1:1:2')),
    Case('getUserRegion/start_after_end', U + 'getUserRegion', args=({'chr2': 1000}, 'chr2:500:100')),
    Case('getUserRegion/negative_start', U + 'getUserRegion', args=({'chr2': 1000}, 'chr2:-5:100')),
    Case('getUserRegion/bad_int', U + 'getUserRegion', args=({'chr2': 1000}, 'chr2:1x:100')),
    Case('getUserRegion/python_int_syntax', U + 'getUserRegion', args=({'chr2': 100000}, 'chr2: 1_0:+2000')),
    Case('getUserRegion/pairs', U + 'getUserRegion', args=([('chr1', 10), ('chr2', 7000)], 'chr2:0:6000:7')),
    Case('getUserRegion/numpy_size', U + 'getUserRegion', args=lambda: ({'chr2': np.int64(7000)}, 'chr2:10:9000')),
    Case('getUserRegion/numpy_size_tiled', U + 'getUserRegion',
         args=lambda: ({'chr2': np.int64(7000)}, 'chr2:10:9000:300')),
    Case('getUserRegion/chunk_limit_float', U + 'getUserRegion', args=({'chr2': 5000000}, 'chr2:0:4000000:7', 1e5)),
    Case('getUserRegion/chunk_limit_int', U + 'getUserRegion', args=({'chr2': 5000000}, 'chr2:0:4000000:7', 100000)),
    Case('getUserRegion/chunk_limit_no_tile', U + 'getUserRegion', args=({'chr2': 5000000}, 'chr2:0:4000000', 12345.7)),
    Case('getUserRegion/tile_zero', U + 'getUserRegion', args=({'chr2': 5000}, 'chr2:10:20:0')),
    Case('getUserRegion/negative_tile', U + 'getUserRegion', args=({'chr2': 5000}, 'chr2:10:20:-7')),
    Case('getUserRegion/bytes', U + 'getUserRegion', args=({'chr2': 5000}, b'chr2:10:20')),

    # expected interactions ---------------------------------------------------
    Case('expected_interactions_in_distance/real', U + 'expected_interactions_in_distance',
         args=lambda: (1000, 2, real_block())),
    Case('expected_interactions_in_distance/real_float32', U + 'expected_interactions_in_distance',
         args=lambda: (1000, 2, real_block(1, np.float32))),
    Case('expected_interactions_in_distance/int32', U + 'expected_interactions_in_distance',
         args=lambda: (10, 1, small_matrix(np.int32))),
    Case('expected_interactions_in_distance/float_count', U + 'expected_interactions_in_distance',
         args=lambda: (9.7, 2.5, small_matrix())),
    Case('expected_interactions_in_distance/explicit_zeros', U + 'expected_interactions_in_distance',
         args=lambda: (5, 1, explicit_zero_matrix())),
    Case('expected_interactions_in_distance/non_square', U + 'expected_interactions_in_distance',
         args=lambda: (5, 1, sp.csr_matrix(np.array([[0, 0, 1.0], [0, 0, 0]])))),
    Case('expected_interactions_in_distance/csc', U + 'expected_interactions_in_distance',
         args=lambda: (1000, 2, sp.csc_matrix(real_block()))),
    Case('expected_interactions_in_distance/coo', U + 'expected_interactions_in_distance',
         args=lambda: (1000, 2, sp.coo_matrix(real_block()))),
    Case('expected_interactions_non_zero/real', U + 'expected_interactions_non_zero', args=lambda: (real_block(),)),
    Case('expected_interactions_non_zero/float32', U + 'expected_interactions_non_zero',
         args=lambda: (real_block(1, np.float32),)),
    Case('expected_interactions_non_zero/explicit_zeros', U + 'expected_interactions_non_zero',
         args=lambda: (explicit_zero_matrix(),)),
    Case('expected_interactions_non_zero/nan', U + 'expected_interactions_non_zero', args=lambda: (nan_inf_matrix(),)),
    Case('expected_interactions/real', U + 'expected_interactions', args=lambda: (real_block(),)),
    Case('expected_interactions/real_float32', U + 'expected_interactions', args=lambda: (real_block(1, np.float32),)),
    Case('expected_interactions/int64', U + 'expected_interactions', args=lambda: (real_block(0, np.int64),)),
    Case('expected_interactions/int32', U + 'expected_interactions', args=lambda: (real_block(1, np.int32),)),
    Case('expected_interactions/threads', U + 'expected_interactions', args=lambda: (real_block(), 3)),
    Case('expected_interactions/negative_threads', U + 'expected_interactions', args=lambda: (small_matrix(), -1)),
    Case('expected_interactions/float_threads', U + 'expected_interactions', args=lambda: (small_matrix(), 2.0)),
    Case('expected_interactions/empty', U + 'expected_interactions', args=lambda: (sp.csr_matrix((4, 4)),)),
    Case('expected_interactions/explicit_zeros', U + 'expected_interactions', args=lambda: (explicit_zero_matrix(),)),
    Case('expected_interactions/explicit_zeros_threads', U + 'expected_interactions',
         args=lambda: (explicit_zero_matrix(), 2), driver='thread_message'),
    Case('expected_interactions/non_square', U + 'expected_interactions',
         args=lambda: (sp.csr_matrix(np.array([[0, 0, 0, 1.0], [0, 0, 0, 0]])),)),
    Case('expected_interactions_thread/sums', U + 'expected_interactions_thread',
         args=lambda: (np.array([1.0, 2.0, 3.0, 4.0, 5.5]), np.array([0, 1, 1, 3, 0]), 0, 4, 5),
         driver='queue_result'),
    Case('expected_interactions_thread/float32', U + 'expected_interactions_thread',
         args=lambda: (np.arange(20, dtype=np.float32) / 3, np.arange(20) % 4, 1, 3, 6), driver='queue_result'),
    Case('expected_interactions_thread/mismatch', U + 'expected_interactions_thread',
         args=lambda: (np.array([1.0, 2.0]), np.array([0, 1, 1]), 0, 2, 3), driver='queue_result'),
    Case('expected_interactions_thread/out_of_bounds', U + 'expected_interactions_thread',
         args=lambda: (np.array([1.0, 2.0]), np.array([0, 5]), 0, 6, 3), driver='queue_result'),

    # compute_zscore ----------------------------------------------------------
    Case('compute_zscore/real', U + 'compute_zscore', args=lambda: (real_block(), None, None)),
    Case('compute_zscore/real_depth', U + 'compute_zscore', args=lambda: (real_block(), 10, None)),
    Case('compute_zscore/real_depth_large', U + 'compute_zscore', args=lambda: (real_block(1), 10000, None)),
    Case('compute_zscore/float_depth', U + 'compute_zscore', args=lambda: (small_matrix(), 2.5, None)),
    Case('compute_zscore/float32', U + 'compute_zscore', args=lambda: (real_block(1, np.float32), 20, None)),
    Case('compute_zscore/int', U + 'compute_zscore', args=lambda: (real_block(0, np.int64), None, None)),
    Case('compute_zscore/empty', U + 'compute_zscore', args=lambda: (sp.csr_matrix((3, 3)), None, None)),
    Case('compute_zscore/threads', U + 'compute_zscore', args=lambda: (small_matrix(), None, 2)),

    # obs/exp -----------------------------------------------------------------
    Case('obs_exp_matrix_lieberman/real', U + 'obs_exp_matrix_lieberman', args=lambda: (real_block(), 2000, 3)),
    Case('obs_exp_matrix_lieberman/float32', U + 'obs_exp_matrix_lieberman',
         args=lambda: (real_block(1, np.float32), 500, 1)),
    Case('obs_exp_matrix_lieberman/int32', U + 'obs_exp_matrix_lieberman',
         args=lambda: (real_block(0, np.int32), 300, 2)),
    Case('obs_exp_matrix_lieberman/explicit_zeros', U + 'obs_exp_matrix_lieberman',
         args=lambda: (explicit_zero_matrix(), 4, 1)),
    Case('obs_exp_matrix_lieberman/empty', U + 'obs_exp_matrix_lieberman', args=lambda: (sp.csr_matrix((3, 3)), 4, 1)),
    Case('obs_exp_matrix_non_zero/real', U + 'obs_exp_matrix_non_zero', args=lambda: (real_block(),)),
    Case('obs_exp_matrix_non_zero/ligation', U + 'obs_exp_matrix_non_zero', args=lambda: (real_block(), True)),
    Case('obs_exp_matrix_non_zero/ligation_float32', U + 'obs_exp_matrix_non_zero',
         args=lambda: (real_block(1, np.float32), True)),
    Case('obs_exp_matrix_non_zero/ligation_int32', U + 'obs_exp_matrix_non_zero',
         args=lambda: (real_block(0, np.int32), True)),
    Case('obs_exp_matrix_non_zero/copy', U + 'obs_exp_matrix_non_zero',
         args=lambda: (real_block(), False, False, True)),
    Case('obs_exp_matrix_non_zero/csc_ligation', U + 'obs_exp_matrix_non_zero',
         args=lambda: (sp.csc_matrix(real_block()), True)),
    Case('obs_exp_matrix_non_zero/explicit_zeros', U + 'obs_exp_matrix_non_zero',
         args=lambda: (explicit_zero_matrix(), True)),
    Case('obs_exp_matrix_non_zero/nan', U + 'obs_exp_matrix_non_zero', args=lambda: (nan_inf_matrix(), False, True, True)),
    Case('obs_exp_matrix/real', U + 'obs_exp_matrix', args=lambda: (real_block(),)),
    Case('obs_exp_matrix/copy_epsilon', U + 'obs_exp_matrix', args=lambda: (real_block(1), False, True)),
    Case('obs_exp_matrix/float32', U + 'obs_exp_matrix', args=lambda: (real_block(1, np.float32),)),
    Case('obs_exp_matrix/int32', U + 'obs_exp_matrix', args=lambda: (real_block(0, np.int32),)),
    Case('obs_exp_matrix/distance', U + 'obs_exp_matrix', args=lambda: (real_block(), True, False, None, 15)),
    Case('obs_exp_matrix/distance_copy', U + 'obs_exp_matrix', args=lambda: (real_block(), False, False, None, 7)),
    Case('obs_exp_matrix/threads', U + 'obs_exp_matrix', args=lambda: (small_matrix(), True, False, 2)),
    Case('obs_exp_matrix/empty', U + 'obs_exp_matrix', args=lambda: (sp.csr_matrix((3, 3)),)),
    Case('obs_exp_matrix/explicit_zeros', U + 'obs_exp_matrix', args=lambda: (explicit_zero_matrix(),)),
    Case('obs_exp_matrix/coo', U + 'obs_exp_matrix', args=lambda: (sp.coo_matrix(small_matrix()),)),

    # strings and bytes -------------------------------------------------------
] + [Case('toString/{}'.format(name), U + 'toString', args=value) for name, value in [
    ('str', ('abc',)), ('bytes', (b'abc',)), ('bytes_non_ascii', (b'\xc3\xa9',)),
    ('npbytes', (np.bytes_(b'chr1'),)), ('npstr', (np.str_('chr1'),)),
    ('list', ([b'a', 'b', [b'c', 1]],)), ('tuple', ((b'a',),)), ('int', (5,)), ('none', (None,)),
    ('array_S', (np.array([b'chr1', b'chrX10', b''], dtype='S6'),)),
    ('array_U', (np.array(['a', 'bc']),)),
    ('array_S_non_ascii', (np.array([b'\xffab']),)),
    ('array_int32', (np.array([1, -20, 300], dtype=np.int32),)),
    ('array_int64', (np.arange(4).reshape(2, 2),)),
    ('array_uint8', (np.array([0, 255], dtype=np.uint8),)),
    ('array_bool', (np.array([True, False]),)),
    ('array_float64', (np.array([0.1, 1e16, 1e-5, -0.0, np.nan, np.inf, 123456789.0, 1.5e300]),)),
    ('array_float32', (np.array([0.1, 1e16, 16777216.0, 3.4e38, 1e-40], dtype=np.float32),)),
    ('array_object', (np.array(['x', 12, None], dtype=object),)),
]] + [Case('toBytes/{}'.format(name), U + 'toBytes', args=value) for name, value in [
    ('str', ('abc',)), ('bytes', (b'abc',)), ('non_ascii', ('é',)), ('list', (['a', b'b', ['c']],)),
    ('array', (np.array(['a']),)), ('int', (3,)),
]] + [
    Case('check_chrom_str_bytes/str_keys', U + 'check_chrom_str_bytes', args=({'chr1': 1}, b'chr1')),
    Case('check_chrom_str_bytes/bytes_keys', U + 'check_chrom_str_bytes', args=({b'chr1': 1}, 'chr1')),
    Case('check_chrom_str_bytes/npbytes_keys', U + 'check_chrom_str_bytes',
         args=lambda: ([np.bytes_(b'chr1')], ['chr1', 'chr2'])),
    Case('check_chrom_str_bytes/list_obj', U + 'check_chrom_str_bytes', args=(['chr1'], [b'chr1', b'x'])),
    Case('check_chrom_str_bytes/same_type', U + 'check_chrom_str_bytes', args=(['chr1'], 'chrX')),
    Case('check_chrom_str_bytes/ints', U + 'check_chrom_str_bytes', args=([1, 2], 'chrX')),
    Case('check_chrom_str_bytes/iterator', U + 'check_chrom_str_bytes',
         args=lambda: (iter([b'a', 'b', b'c']), 'x'), compare_args=False),
    Case('check_chrom_str_bytes/empty', U + 'check_chrom_str_bytes', args=([], 'x')),
] + [Case('change_chrom_names/{}'.format(name), U + 'change_chrom_names', args=value) for name, value in [
    ('ucsc', ('chr1',)), ('ensembl', ('X',)), ('bytes', (b'chr2',)), ('npstr', (np.str_('chrM'),)),
    ('empty', ('',)), ('unicode', ('chré',)), ('list', (['chr1'],)), ('int', (1,)),
]] + [
    Case('opener/plain', U + 'opener', args=(File('f.txt'),), files={'f.txt': b'line1\nline2\n'}, driver='read_file'),
    Case('opener/gzip', U + 'opener', args=(File('f.gz'),), files={'f.gz': gzip.compress(b'zipped\n', mtime=0)},
         driver='read_file'),
    Case('opener/one_byte', U + 'opener', args=(File('f'),), files={'f': b'\x1f'}, driver='read_file'),
    Case('opener/missing', U + 'opener', args=(File('missing'),), driver='read_file'),
] + [Case('remove_non_ascii/{}'.format(name), U + 'remove_non_ascii', args=value) for name, value in [
    ('ascii', ('plain text',)), ('latin', ('Knösós Æ',)), ('cjk', ('北亰',)),
    ('greek', ('Κνωσός',)), ('astral', ('a\U0001f600b\U000e0001',)),
    ('surrogate', ('x\udc80y',)), ('bytes', (b'abc',)), ('empty', ('',)),
]] + [
    Case('check_cooler/cool_name', U + 'check_cooler', args=('missing.cool',)),
    Case('check_cooler/real_cool', U + 'check_cooler', args=(os.path.join(TEST_DATA, 'small_test_matrix.cool'),)),
    Case('check_cooler/real_h5', U + 'check_cooler', args=(os.path.join(TEST_DATA, 'small_test_matrix.h5'),)),
    Case('check_cooler/real_mcool', U + 'check_cooler', args=(os.path.join(TEST_DATA, 'matrix.mcool'),)),
    Case('check_cooler/mcool_uri', U + 'check_cooler',
         args=(os.path.join(TEST_DATA, 'hicBuildMatrix', 'multi_small_test_matrix.mcool') + '::/resolutions/10000',)),
    Case('check_cooler/missing_mcool', U + 'check_cooler', args=('missing.mcool',)),
    Case('check_cooler/missing_h5', U + 'check_cooler', args=('missing.h5',)),
    Case('check_cooler/bad_uri', U + 'check_cooler', args=('a::b::c',)),
    Case('check_cooler/text_file', U + 'check_cooler', args=(File('a.txt'),), files={'a.txt': b'not hdf5'}),
] + [Case('in_units/{}'.format(i), U + 'in_units', args=(value,)) for i, value in enumerate([
    0, 12.345, 1500, 1500.0001, 1501, 23456.5, 1.5e6, 1500000.0001, 2.5e7, -3e6, '2000', float('nan'),
    float('inf'), -float('inf'), np.float32(1234567.0), 2 ** 70, 999.995])
]
