"""Equivalence cases for hicexplorer.chicQualityControl (apieq.py).

compute_sparsity runs as chicQualityControl's main() calls it, on the cHi-C
test matrices with the quality control reference points of the repository
(which include reference points the matrix cannot place); the driver returns
what it put into the queue.
"""
import argparse
import os

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicQualityControl.'

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


def sparsity(target, matrix, reference_points, fixate):
    import queue
    viewpoint = _viewpoint(matrix)
    if isinstance(reference_points, str):
        reference_points, _ = viewpoint.readReferencePointFile(reference_points)
    q = queue.Queue()
    target(reference_points, viewpoint, argparse.Namespace(fixateRange=fixate), q)
    return q.get()


CASES = [
    Case('compute_sparsity/qc_reference_points', M + 'compute_sparsity', driver='sparsity',
         args=(FL, os.path.join(CHIC, 'referencePoints_qc.bed'), 500000)),
    Case('compute_sparsity/qc_reference_points_mb', M + 'compute_sparsity', driver='sparsity',
         args=(MB, os.path.join(CHIC, 'referencePoints_qc.bed'), 200000)),
    Case('compute_sparsity/edge_cases', M + 'compute_sparsity', driver='sparsity',
         args=(FL, os.path.join(CHIC, 'referencePoints_edge_cases.bed'), 500000)),
    Case('compute_sparsity/none_and_other_chromosome', M + 'compute_sparsity', driver='sparsity',
         args=(FL, [None, ('chr2', '1', '1'), ('chr1', '4487435', '4487435')], 1000)),
    Case('compute_sparsity/outside', M + 'compute_sparsity', driver='sparsity',
         args=(FL, [('chr1', '300000000', '300000000'), ('chr1', '70931031', '70093803')], 1000)),
    Case('compute_sparsity/bad_number', M + 'compute_sparsity', driver='sparsity',
         args=(FL, [('chr1', 'x', '1')], 1000)),
]
