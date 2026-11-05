"""Equivalence cases for hicexplorer.chicViewpoint (apieq.py).

compute_viewpoint runs on the cHi-C test matrices with the repository's
background model and reference point files, as chicViewpoint's main() calls
it (a Viewpoint with the matrix, the parsed arguments, a queue); the driver
returns what it put into the queue. adjustViewpointData is called where
chicViewpoint calls it, on a reference point near the chromosome end, and on
constructed inputs; compute_x_fold on the dtypes and shapes it meets.
"""
import argparse
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
CHIC = os.path.join(TEST_DATA, 'cHi-C')
M = 'hicexplorer.chicViewpoint.'

FL = os.path.join(CHIC, 'FL-E13-5_chr1.cool')
MB = os.path.join(CHIC, 'MB-E10-5_chr1.cool')
BACKGROUND = os.path.join(CHIC, 'background.txt')

_MATRICES = {}


def _viewpoint(matrix=None):
    """A Viewpoint on hiCMatrix(matrix), the matrix loaded once per interpreter
    (no function here changes it)."""
    import hicmatrix.HiCMatrix as hm
    from hicexplorer.lib.viewpoint import Viewpoint
    if matrix is not None and matrix not in _MATRICES:
        _MATRICES[matrix] = hm.hiCMatrix(matrix)
    return Viewpoint(_MATRICES[matrix] if matrix is not None else None)


def near_end_adjust(target, matrix, reference_point, prange):
    """adjustViewpointData where compute_viewpoint calls it."""
    viewpoint = _viewpoint(matrix)
    model = viewpoint.readBackgroundDataFile(BACKGROUND, prange, 500000, pMean=True)
    region_start, region_end, _range = viewpoint.calculateViewpointRange(reference_point, prange)
    data, _ = viewpoint.computeViewpoint(reference_point, reference_point[0], region_start, region_end)
    background = viewpoint.interactionBackgroundData(model, _range).flatten()
    return len(data), len(background), target(viewpoint, data, background, reference_point, region_start, region_end)


def on_matrix(target, matrix, *args):
    return target(_viewpoint(matrix), *args)


def compute(target, matrix, reference_file, prange, fixate, average, threads_put=None):
    import queue
    viewpoint = _viewpoint(matrix)
    reference_points, genes = viewpoint.readReferencePointFile(reference_file)
    background = viewpoint.readBackgroundDataFile(BACKGROUND, prange, fixate)
    mean = viewpoint.readBackgroundDataFile(BACKGROUND, prange, fixate, pMean=True)
    args = argparse.Namespace(fixateRange=fixate, range=prange, averageContactBin=average)
    q = queue.Queue()
    target(viewpoint, args, q, reference_points, genes, matrix, background, mean)
    result = q.get()
    assert q.empty()
    return result


def compute_points(target, matrix, reference_points, genes, prange, fixate, average):
    import queue
    viewpoint = _viewpoint(matrix)
    background = viewpoint.readBackgroundDataFile(BACKGROUND, prange, fixate)
    mean = viewpoint.readBackgroundDataFile(BACKGROUND, prange, fixate, pMean=True)
    args = argparse.Namespace(fixateRange=fixate, range=prange, averageContactBin=average)
    q = queue.Queue()
    target(viewpoint, args, q, reference_points, genes, matrix, background, mean)
    return q.get()


CASES = [
    # adjustViewpointData ------------------------------------------------------
    Case('adjustViewpointData/near_end', M + 'adjustViewpointData', driver='near_end_adjust',
         args=(FL, ('chr1', '197145200', '197145200'), (200000, 200000))),
    Case('adjustViewpointData/near_end_wide_range', M + 'adjustViewpointData', driver='near_end_adjust',
         args=(MB, ('chr1', '197100000', '197100000'), (300000, 300000))),
    Case('adjustViewpointData/data_longer', M + 'adjustViewpointData', driver='on_matrix',
         args=(FL, [0.5, 1.5, 2.5, 3.5, 4.5, 5.5], [0.1, 0.2, 0.3], ('chr1', '14300280', '14300280'),
               14298000, 14302000)),
    Case('adjustViewpointData/background_longer', M + 'adjustViewpointData', driver='on_matrix',
         args=lambda: (FL, np.array([0.5, 1.5]), np.array([0.1, 0.2, 0.3, 0.4, 0.5]),
                       ('chr1', '14300280', '14300280'), 14298000, 14302000)),
    Case('adjustViewpointData/empty_range', M + 'adjustViewpointData', driver='on_matrix',
         args=(FL, [1.0], [2.0], ('chr1', '14300280', '14300280'), 14302000, 14298000)),

    # compute_x_fold -------------------------------------------------------------
    Case('compute_x_fold/arrays', M + 'compute_x_fold',
         args=lambda: (np.array([0.5, 0.25, 0.0, 1.0]), np.array([0.25, 0.5, 0.0, 0.0]))),
    Case('compute_x_fold/array_and_list', M + 'compute_x_fold',
         args=lambda: (np.array([0.5, 0.25, 0.125]), [0.25, np.float64(0.5), 2.0])),
    Case('compute_x_fold/float32_and_list', M + 'compute_x_fold',
         args=lambda: (np.array([0.5, 0.3], dtype=np.float32), [0.7, 0.1])),
    Case('compute_x_fold/float32_and_float32', M + 'compute_x_fold',
         args=lambda: (np.array([0.5, 0.3], dtype=np.float32), np.array([0.7, 0.1], dtype=np.float32))),
    Case('compute_x_fold/ints', M + 'compute_x_fold', args=lambda: (np.array([1, 2, 3]), np.array([2, 2, 0]))),
    Case('compute_x_fold/scalar', M + 'compute_x_fold', args=lambda: (np.array([1.0, 2.0]), 4.0)),
    Case('compute_x_fold/lists', M + 'compute_x_fold', args=([1.0], [2.0])),
    Case('compute_x_fold/shapes', M + 'compute_x_fold', args=lambda: (np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0]))),

    # compute_viewpoint ----------------------------------------------------------
    Case('compute_viewpoint/three_genes', M + 'compute_viewpoint', driver='compute',
         args=(FL, os.path.join(CHIC, 'referencePoints_chicViewpoint.bed'), [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/no_smoothing', M + 'compute_viewpoint', driver='compute',
         args=(MB, os.path.join(CHIC, 'referencePoints_chicViewpoint.bed'), [100000, 150000], 500000, 0)),
    Case('compute_viewpoint/one_gene', M + 'compute_viewpoint', driver='compute',
         args=(MB, os.path.join(CHIC, 'referencePoints_chicViewpoint_one_gene.bed'), [200000, 200000], 500000, 4)),
    Case('compute_viewpoint/near_end', M + 'compute_viewpoint', driver='compute_points',
         args=(FL, [('chr1', '197145200', '197145200')], ['NearEnd'], [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/near_end_no_smoothing', M + 'compute_viewpoint', driver='compute_points',
         args=(MB, [('chr1', '197145200', '197145200')], ['NearEnd'], [200000, 200000], 500000, 0)),
    Case('compute_viewpoint/near_start', M + 'compute_viewpoint', driver='compute_points',
         args=(FL, [('chr1', '150500', '150500')], ['NearStart'], [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/wide', M + 'compute_viewpoint', driver='compute_points',
         args=(MB, [('chr1', '24202042', '24205042')], ['Wide'], [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/unknown_chromosome', M + 'compute_viewpoint', driver='compute_points',
         args=(FL, [('chr1', '4487435', '4487435'), ('chr2', '1', '1')], ['a', 'b'], [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/missing_gene', M + 'compute_viewpoint', driver='compute_points',
         args=(FL, [('chr1', '4487435', '4487435')], [], [200000, 200000], 500000, 5)),
    Case('compute_viewpoint/empty', M + 'compute_viewpoint', driver='compute_points',
         args=(FL, [], [], [200000, 200000], 500000, 5)),
]
