"""Equivalence cases for hicexplorer.hicAverageRegions (apieq.py).

The functions query a hiCMatrix, so every case is a driver that loads one of
the repository's matrices with hicmatrix in its own interpreter and calls the
function on it: small_test_matrix_50kb_res as cool (Python int bins) and as h5
(numpy.int64 bins), and Li_et_al_2015 (restriction fragment bins).

Deviations from 3.7.6 (cpp/python/src/api/hicAverageRegions.cpp,
hicmatrix_queries.hpp):
  * ranges that are neither Python ints nor signed numpy integers raise
    TypeError naming HiCExplorer 4 (3.7.6 computed floats from float ranges);
  * a position that int() cannot convert raises int()'s ValueError; hicmatrix
    logs it with a '%d' format that fails inside logging.
"""
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicAverageRegions.'

COOL = 'small_test_matrix_50kb_res.cool'
H5 = 'small_test_matrix_50kb_res.h5'
LI = 'Li_et_al_2015.cool'


def load(name):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(os.path.join(TEST_DATA, name))


def on_matrix(target, matrix, *args):
    """target(hiCMatrix of `matrix`, *args); the matrix itself is not returned
    (its repr holds an address)."""
    return target(load(matrix), *args)


def regions_file(target, matrix, path, rng, mapping, bins):
    """The loop of hicAverageRegions.main over a BED file of the test data."""
    hic = load(matrix)
    out = []
    with open(os.path.join(TEST_DATA, path)) as handle:
        for line in handle:
            fields = line.strip().split('\t')
            viewpoint = (fields[0], fields[1], fields[2])
            try:
                out.append(target(hic, viewpoint, rng, mapping))
            except Exception as error:  # compared as part of the result
                out.append((type(error).__name__, str(error)))
    return out


def case(id, name, matrix, *args, **kw):
    return Case(id, M + name, driver='on_matrix', args=(matrix,) + args, **kw)


CASES = [
    # calculateViewpointRange
    case('range/start', 'calculateViewpointRange', COOL, ('chr2L', '2500000', '3500000'), [100000, 200000], 'start'),
    case('range/end', 'calculateViewpointRange', COOL, ('chr2L', '2500000', '3500000'), [100000, 200000], 'end'),
    case('range/center', 'calculateViewpointRange', COOL, ('chr2L', '2500001', '3500000'), [100000, 200000], 'center'),
    case('range/center_float_text', 'calculateViewpointRange', COOL, ('chr2L', '2500001.5', '3500000'), [10, 20],
         'center'),
    case('range/clamp_start', 'calculateViewpointRange', COOL, ('chr2L', '50000', '60000'), [100000, 200000], 'start'),
    case('range/clamp_end', 'calculateViewpointRange', COOL, ('chr2RHet', '3200000', '3250000'), [100000, 2000000],
         'end'),
    case('range/clamp_both', 'calculateViewpointRange', H5, ('chr4', '10', '20'), [100000000, 100000000], 'start'),
    case('range/h5', 'calculateViewpointRange', H5, ('chr3L', '1000000', '1500000'), [100000, 200000], 'center'),
    case('range/numpy_range', 'calculateViewpointRange', COOL, ('chr2L', '2500000', '3500000'),
         np.array([100000, 200000]), 'start'),
    case('range/numpy_int32_range', 'calculateViewpointRange', COOL, ('chr2L', '2500000', '3500000'),
         [np.int32(1), np.int32(2)], 'end'),
    case('range/int_viewpoint', 'calculateViewpointRange', COOL, ('chr2L', 2500000, 3500000), (5, 7), 'start'),
    case('range/li', 'calculateViewpointRange', LI, ('X', '10000', '12000'), [5000, 5000], 'center'),
    case('range/bad_mapping', 'calculateViewpointRange', COOL, ('chr2L', '2500000', '3500000'), [1, 2], 'middle'),
    case('range/unknown_chrom', 'calculateViewpointRange', COOL, ('chr9', '2500000', '3500000'), [1, 2], 'start'),
    case('range/not_a_number', 'calculateViewpointRange', COOL, ('chr2L', '25x', '3500000'), [1, 2], 'start'),
    Case('range/regions_file', M + 'calculateViewpointRange', driver='regions_file',
         args=(COOL, 'hicAverageRegions/regions_multi.bed', [500000, 500000], 'center', False)),

    # getBinIndices
    case('bins/inside', 'getBinIndices', COOL, ('chr2L', 2500000, 3500000)),
    case('bins/text_positions', 'getBinIndices', H5, ('chr3R', '2500000', '3500000')),
    case('bins/same_bin', 'getBinIndices', COOL, ('chrX', 10, 20)),
    case('bins/outside', 'getBinIndices', COOL, ('chr2L', 2500000, 350000000)),
    case('bins/negative', 'getBinIndices', COOL, ('chr2L', -5, 100)),
    case('bins/bytes_name', 'getBinIndices', COOL, (b'chr2L', 0, 100)),
    case('bins/unknown_chrom', 'getBinIndices', H5, ('chrM', 0, 100)),
    case('bins/li_fragments', 'getBinIndices', LI, ('X', 1000, 100000)),
    case('bins/float_positions', 'getBinIndices', COOL, ('chr2L', 2500000.7, 3500000.2)),

    # calculateViewpointRangeBins
    case('range_bins/start', 'calculateViewpointRangeBins', COOL, ('chr2L', '2500000', '3500000'), [10, 20], 'start'),
    case('range_bins/end', 'calculateViewpointRangeBins', COOL, ('chr2L', '2500000', '3500000'), [10, 20], 'end'),
    case('range_bins/center', 'calculateViewpointRangeBins', H5, ('chr2L', '2500000', '3500001'), [10, 20],
         'center'),
    case('range_bins/other_mapping_is_center', 'calculateViewpointRangeBins', COOL,
         ('chr2L', '2500000', '3500001'), [10, 20], 'start_end'),
    case('range_bins/clamp_start', 'calculateViewpointRangeBins', COOL, ('chr3L', '0', '100000'), [10, 20], 'start'),
    case('range_bins/clamp_end', 'calculateViewpointRangeBins', H5, ('chr3LHet', '2500000', '2550000'), [10, 200],
         'end'),
    case('range_bins/numpy_range', 'calculateViewpointRangeBins', COOL, ('chr2L', '2500000', '3500000'),
         np.array([3, 4], dtype=np.int64), 'start'),
    case('range_bins/outside', 'calculateViewpointRangeBins', COOL, ('chr2L', '2500000', '350000000'), [3, 4],
         'start'),
    case('range_bins/unknown_chrom', 'calculateViewpointRangeBins', COOL, ('chr9', '1', '2'), [3, 4], 'start'),
    case('range_bins/li', 'calculateViewpointRangeBins', LI, ('X', '20000', '30000'), [30, 30], 'end'),
    Case('range_bins/regions_file', M + 'calculateViewpointRangeBins', driver='regions_file',
         args=(H5, 'hicAverageRegions/regions_multi.bed', [50, 50], 'start', True)),
]
