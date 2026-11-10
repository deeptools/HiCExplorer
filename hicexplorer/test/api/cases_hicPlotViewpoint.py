"""Equivalence cases for hicexplorer.hicPlotViewpoint (apieq.py).

getViewpointValues loads the matrix itself (hm.hiCMatrix in 3.7.6, the C++
matrix layer here) and returns the bins in hiCMatrix's types; the cases use a
cool file (str and int bins, float64 values), an h5 file (numpy.str_ and
numpy.int64 bins, int64 values) and the cHi-C matrix, with and without the
interaction list and the chromosome filter (a list, and a str, which
keepOnlyTheseChr iterates by character).
"""
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicPlotViewpoint.'

LI = os.path.join(TEST_DATA, 'Li_et_al_2015.cool')
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
FL = os.path.join(TEST_DATA, 'cHi-C', 'FL-E13-5_chr1.cool')


def exits(target, *args, **kwargs):
    try:
        return target(*args, **kwargs)
    except SystemExit as error:
        return 'SystemExit', error.code


CASES = [
    Case('relabelTicks/kb', M + 'relabelTicks', args=(500,)),
    Case('relabelTicks/below_mb', M + 'relabelTicks', args=(999999,)),
    Case('relabelTicks/mb', M + 'relabelTicks', args=(1e6,)),
    Case('relabelTicks/float', M + 'relabelTicks', args=(2500000.7,)),
    Case('relabelTicks/numpy', M + 'relabelTicks', args=lambda: (np.int64(1234567),)),
    Case('relabelTicks/negative', M + 'relabelTicks', args=(-1500,)),
    Case('relabelTicks/str', M + 'relabelTicks', args=('5',)),

    Case('getViewpointValues/cool_point', M + 'getViewpointValues',
         args=(LI, ['X', '3000000'], 'X', 2500000, 3500000)),
    Case('getViewpointValues/cool_region_interactions', M + 'getViewpointValues',
         args=(LI, ['X', '3000000', '3010000'], 'X', 2990000, 3020000, True)),
    Case('getViewpointValues/cool_chromosome_str', M + 'getViewpointValues',
         args=(LI, ['X', '3000000'], 'X', 2900000, 3100000, None, 'X')),
    Case('getViewpointValues/cool_chromosome_str_fails', M + 'getViewpointValues',
         args=(LI, ['2L', '3000000'], '2L', 2900000, 3100000, None, '2L')),
    Case('getViewpointValues/cool_chromosome_list', M + 'getViewpointValues',
         args=(LI, ['X', '3000000'], 'X', 2950000, 3050000, [], ['X'])),
    Case('getViewpointValues/h5_interactions', M + 'getViewpointValues',
         args=(H5, ['chr2L', '1000000', '1100000'], 'chr2L', 800000, 1300000, [])),
    Case('getViewpointValues/h5_chromosome_list', M + 'getViewpointValues',
         args=(H5, ['chrX', '2000000'], 'chrX', 1000000, 3000000, None, ['chr3L', 'chrX'])),
    Case('getViewpointValues/chic_matrix', M + 'getViewpointValues',
         args=(FL, ['chr1', '14300280'], 'chr1', 14200000, 14400000)),
    Case('getViewpointValues/reference_outside', M + 'getViewpointValues',
         args=(LI, ['X', '900000000'], 'X', 2500000, 3500000)),
    Case('getViewpointValues/region_outside', M + 'getViewpointValues',
         args=(LI, ['X', '3000000'], 'X', 2500000, 900000000)),
    Case('getViewpointValues/unknown_chromosome', M + 'getViewpointValues',
         args=(LI, ['chrZ', '3000000'], 'chrZ', 2500000, 3500000)),
    Case('getViewpointValues/reversed_region', M + 'getViewpointValues',
         args=(LI, ['X', '3000000'], 'X', 3500000, 2500000)),
    Case('getViewpointValues/bad_reference', M + 'getViewpointValues', driver='exits',
         args=(LI, ['X'], 'X', 2500000, 3500000)),
]
