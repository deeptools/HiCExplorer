"""Equivalence cases for hicexplorer.hicMergeMatrixBins (apieq.py).

The hiCMatrix objects are loaded with hicmatrix from the repository's test
matrices in both interpreters: small_test_matrix_50kb_res.h5 (int64, with
NaN bins, so remove_nans_if_needed masks them), its cool version (int32),
the ICE corrected matrix of the hicCorrectMatrix tests (float64), and a
float64 copy of the first with random bias factors, whose sums show the
summation order of scipy. The functions change the hiCMatrix in place and
return it; the drivers return whether the result is the argument and the
state of the object (matrix, bin table, nan_bins, correction factors,
masked bin bookkeeping, chrBinBoundaries and the interval trees), because
a hiCMatrix has no value representation of its own. reduce_matrix's timing
debug records are compared with the time masked.
"""
import logging
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicMergeMatrixBins.'


def load(name, biased=False, dtype=None, nan_bins=None):
    from hicmatrix import HiCMatrix as hm
    hic = hm.hiCMatrix(os.path.join(TEST_DATA, name))
    if biased:
        rng = np.random.default_rng(5)
        bias = rng.uniform(0.2, 4.0, hic.matrix.shape[0])
        coo = hic.matrix.tocoo()
        coo.data = coo.data * bias[coo.row] * bias[coo.col] / 3.0
        hic.matrix = coo.tocsr()
    if dtype is not None:
        hic.matrix.data = hic.matrix.data.astype(dtype)
    if nan_bins is not None:
        hic.nan_bins = nan_bins
    return hic


def small_block(dtype=np.int64):
    """setMatrix with the first two chromosomes of the 50 kb matrix."""
    from hicmatrix import HiCMatrix as hm
    full = hm.hiCMatrix(os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5'))
    n = full.chrBinBoundaries['chr3RHet'][1]
    hic = hm.hiCMatrix()
    hic.nan_bins = []
    hic.setMatrix(full.matrix[:n, :n].astype(dtype), list(full.cut_intervals[:n]))
    return hic


class _MaskTimes(logging.Filter):
    def filter(self, record):
        message = record.getMessage()
        if message.startswith('time '):
            record.msg = message.split(':')[0] + ': <masked>'
            record.args = ()
        return True


def state(hic):
    trees = getattr(hic, 'interval_trees', None)
    boundaries = getattr(hic, 'chrBinBoundaries', None)
    return {
        'matrix': hic.matrix,
        'cut_intervals': hic.cut_intervals,
        'nan_bins': hic.nan_bins,
        'correction_factors': hic.correction_factors,
        'distance_counts': hic.distance_counts,
        'orig_bin_ids': hic.orig_bin_ids,
        'orig_cut_intervals': hic.orig_cut_intervals,
        'chrBinBoundaries': None if boundaries is None else list(boundaries.items()),
        'interval_trees': None if trees is None else
        [(name, sorted((iv.begin, iv.end, iv.data) for iv in tree)) for name, tree in trees.items()],
    }


def run(target, hic, *args):
    logger = logging.getLogger('hicexplorer.reduceMatrix')
    mask = _MaskTimes()
    logger.addFilter(mask)
    try:
        result = target(hic, *args)
    finally:
        logger.removeFilter(mask)
    return result is hic, state(result)


def run_state_on_error(target, hic, *args):
    try:
        return run(target, hic, *args)
    except Exception as error:
        return type(error).__name__, str(error), state(hic)


def h5(**kw):
    return lambda: (load('small_test_matrix_50kb_res.h5', **kw),)


def h5_with(*extra, **kw):
    return lambda: (load('small_test_matrix_50kb_res.h5', **kw),) + extra


R = dict(driver='run', compare_args=False)
CASES = [
    # remove_nans_if_needed --------------------------------------------------
    Case('remove_nans/h5_nan_bins', M + 'remove_nans_if_needed', args=h5(), **R),
    Case('remove_nans/cool_nan_bins', M + 'remove_nans_if_needed',
         args=lambda: (load('small_test_matrix_50kb_res.cool'),), **R),
    Case('remove_nans/no_nan_bins', M + 'remove_nans_if_needed', args=h5(nan_bins=[]), **R),
    Case('remove_nans/none', M + 'remove_nans_if_needed', args=h5(nan_bins=None), **R),
    Case('remove_nans/ice_corrected', M + 'remove_nans_if_needed',
         args=lambda: (load('hicCorrectMatrix/small_test_matrix_ICEcorrected_chrUextra_chr3LHet.h5'),), **R),
    # merge_bins ------------------------------------------------------------
    Case('merge_bins/h5_2', M + 'merge_bins', args=h5_with(2), **R),
    Case('merge_bins/h5_3', M + 'merge_bins', args=h5_with(3), **R),
    Case('merge_bins/h5_10', M + 'merge_bins', args=h5_with(10), **R),
    Case('merge_bins/h5_1', M + 'merge_bins', args=h5_with(1), **R),
    Case('merge_bins/h5_numpy_int', M + 'merge_bins', args=h5_with(np.int64(4)), **R),
    Case('merge_bins/cool_5', M + 'merge_bins', args=lambda: (load('small_test_matrix_50kb_res.cool'), 5), **R),
    Case('merge_bins/biased_float64', M + 'merge_bins', args=h5_with(3, biased=True), **R),
    Case('merge_bins/biased_float32', M + 'merge_bins', args=h5_with(3, biased=True, dtype=np.float32), **R),
    Case('merge_bins/ice_corrected', M + 'merge_bins',
         args=lambda: (load('hicCorrectMatrix/small_test_matrix_ICEcorrected_chrUextra_chr3LHet.h5'), 4), **R),
    Case('merge_bins/set_matrix_block', M + 'merge_bins', args=lambda: (small_block(), 4), **R),
    Case('merge_bins/set_matrix_block_int32', M + 'merge_bins', args=lambda: (small_block(np.int32), 7), **R),
    Case('merge_bins/larger_than_chromosomes', M + 'merge_bins', args=lambda: (small_block(), 100), **R),
    Case('merge_bins/zero', M + 'merge_bins', args=lambda: (small_block(), 0), driver='run_state_on_error',
         compare_args=False),
    # running_window_merge ---------------------------------------------------
    Case('running_window/h5_3', M + 'running_window_merge', args=h5_with(3), **R),
    Case('running_window/h5_5', M + 'running_window_merge', args=h5_with(5), **R),
    Case('running_window/h5_1', M + 'running_window_merge', args=h5_with(1), **R),
    Case('running_window/cool_3', M + 'running_window_merge',
         args=lambda: (load('small_test_matrix_50kb_res.cool'), 3), **R),
    Case('running_window/biased_float64_3', M + 'running_window_merge', args=h5_with(3, biased=True), **R),
    Case('running_window/biased_float64_5', M + 'running_window_merge', args=h5_with(5, biased=True), **R),
    Case('running_window/biased_float32_3', M + 'running_window_merge',
         args=h5_with(3, biased=True, dtype=np.float32), **R),
    Case('running_window/set_matrix_block', M + 'running_window_merge', args=lambda: (small_block(), 3), **R),
    Case('running_window/even', M + 'running_window_merge', args=lambda: (small_block(), 4),
         driver='run_state_on_error', compare_args=False),
]
