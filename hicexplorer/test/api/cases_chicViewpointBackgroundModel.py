"""Equivalence cases for hicexplorer.chicViewpointBackgroundModel (apieq.py).

compute_background runs as chicViewpointBackgroundModel's main() calls it,
on the cHi-C test matrices with the repository's reference point files, with
and without smoothing; the driver returns what it put into the queue.
"""
import argparse
import os

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicViewpointBackgroundModel.'

FL = os.path.join(CHIC, 'FL-E13-5_chr1.cool')
MB = os.path.join(CHIC, 'MB-E10-5_chr1.cool')

_MATRICES = {}


def _viewpoint(matrix):
    """A Viewpoint on hiCMatrix(matrix), loaded once per interpreter."""
    import hicmatrix.HiCMatrix as hm
    from hicexplorer.lib.viewpoint import Viewpoint
    if matrix not in _MATRICES:
        _MATRICES[matrix] = hm.hiCMatrix(matrix)
    return Viewpoint(_MATRICES[matrix])


def background(target, matrix, reference_points, fixate, average):
    import queue
    viewpoint = _viewpoint(matrix)
    if isinstance(reference_points, str):
        reference_points, _ = viewpoint.readReferencePointFile(reference_points)
    q = queue.Queue()
    target(reference_points, viewpoint, argparse.Namespace(fixateRange=fixate, averageContactBin=average), q)
    return q.get()


CASES = [
    Case('compute_background/reference_points', M + 'compute_background', driver='background',
         args=(FL, os.path.join(CHIC, 'referencePoints.bed'), 50000, 5)),
    Case('compute_background/no_smoothing', M + 'compute_background', driver='background',
         args=(MB, os.path.join(CHIC, 'referencePoints_chicViewpoint.bed'), 20000, 0)),
    Case('compute_background/even_window', M + 'compute_background', driver='background',
         args=(MB, os.path.join(CHIC, 'referencePoints_chicViewpoint.bed'), 10000, 4)),
    Case('compute_background/near_borders', M + 'compute_background', driver='background',
         args=(FL, [('chr1', '197145200', '197145200'), ('chr1', '150500', '150500'),
                    ('chr1', '24202042', '24205042')], 200000, 5)),
    Case('compute_background/unknown_chromosome', M + 'compute_background', driver='background',
         args=(FL, [('chr1', '4487435', '4487435'), ('chrX', '1', '1')], 20000, 5)),
    Case('compute_background/position_outside', M + 'compute_background', driver='background',
         args=(FL, [('chr1', '300000000', '300000000')], 20000, 5)),
    Case('compute_background/empty', M + 'compute_background', driver='background', args=(FL, [], 20000, 5)),
]
