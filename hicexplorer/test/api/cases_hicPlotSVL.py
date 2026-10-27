"""Equivalence cases for hicexplorer.hicPlotSVL (apieq.py).

compute_relation_short_long_range puts its result on a queue and returns
None; the driver passes a queue.Queue and returns what was put, together
with the bin size the function cached on the hiCMatrix object. The matrices
are the repository's: small_test_matrix_50kb_res (int32 cool, int64 h5),
Li_et_al_2015 (float64 with restriction fragment bins), and float32 copies
made in the driver.

Deviations from 3.7.6 (cpp/python/src/api/hicPlotSVL.cpp): a cool file the
C++ reader cannot open raises OSError with the reader's message; the matrix
of a hiCMatrix object must be csr.
"""
import os
import queue

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
T = 'hicexplorer.hicPlotSVL.compute_relation_short_long_range'

COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
LI_COOL = os.path.join(TEST_DATA, 'Li_et_al_2015.cool')
LI_H5 = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')


def load(path):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


def on_object(target, path, chromosomes, distance, dtype=None, explicit_zero=False):
    hic = load(path)
    if dtype is not None:
        hic.matrix.data = hic.matrix.data.astype(dtype)
    if explicit_zero:
        first, last = hic.getChrBinRange(chromosomes[0])
        row = first + 3
        for k in range(hic.matrix.indptr[row], hic.matrix.indptr[row + 1]):
            if first <= hic.matrix.indices[k] < last:
                hic.matrix.data[k] = 0
                break
    if chromosomes is None:
        chromosomes = list(hic.chrBinBoundaries)
    q = queue.Queue()
    returned = target(hic, chromosomes, distance, False, q)
    return returned, q.get_nowait(), q.empty(), hic.bin_size


def on_file(target, path, chromosomes, distance):
    q = queue.Queue()
    returned = target(path, chromosomes, distance, True, q)
    return returned, q.get_nowait(), q.empty()


ALL_COOL = ['chr2RHet', 'chr3RHet', 'chr2LHet', 'chr4', 'chrYHet', 'chr3L', 'chr2L', 'chrU', 'chrX', 'chrXHet',
            'chr2R', 'chr3R', 'chrUextra', 'chrM', 'chr3LHet']

CASES = [
    Case('object/h5_all', T, driver='on_object', args=(H5, None, 2000000)),
    Case('object/h5_short', T, driver='on_object', args=(H5, ['chr2L', 'chrX', 'chr4'], 100000)),
    Case('object/cool_int32', T, driver='on_object', args=(COOL, None, 2000000)),
    Case('object/float_distance', T, driver='on_object', args=(COOL, ['chr3R', 'chr2R'], 333333.3)),
    Case('object/distance_zero', T, driver='on_object', args=(H5, ['chr3L'], 0)),
    Case('object/float32', T, driver='on_object', args=(H5, ['chr2L', 'chr3R', 'chrX'], 1000000, np.float32)),
    Case('object/float64', T, driver='on_object', args=(H5, ['chr2L', 'chr3R'], 1000000, np.float64)),
    Case('object/li_fragments', T, driver='on_object', args=(LI_H5, ['X'], 50000)),
    Case('object/li_cool', T, driver='on_object', args=(LI_COOL, ['X'], 50000)),
    Case('object/tuple_of_chromosomes', T, driver='on_object', args=(H5, ('chr2L', 'chr2L'), 2000000)),
    Case('object/empty', T, driver='on_object', args=(H5, [], 2000000)),
    Case('object/unknown_chromosome', T, driver='on_object', args=(H5, ['chr2L', 'chr9'], 2000000)),
    Case('object/explicit_zero', T, driver='on_object', args=(H5, ['chr2L'], 2000000, None, True)),
    Case('file/cool_all', T, driver='on_file', args=(COOL, ALL_COOL, 2000000)),
    Case('file/cool_short', T, driver='on_file', args=(COOL, ['chr2L', 'chrX'], 150000)),
    Case('file/li', T, driver='on_file', args=(LI_COOL, ['X'], 1000000)),
    Case('file/single_bin_chromosome', T, driver='on_file', args=(COOL, ['chrM'], 2000000)),
    Case('file/unknown_chromosome', T, driver='on_file', args=(COOL, ['chr9'], 2000000)),
]
