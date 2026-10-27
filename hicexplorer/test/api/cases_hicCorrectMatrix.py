"""Equivalence cases for hicexplorer.hicCorrectMatrix (apieq.py).

The matrices are the repository's test matrices loaded with hicmatrix in
both interpreters: small_test_matrix_50kb_res.h5 (int64 counts with NaN
bins and runs of them), the same with random bias factors (float64), its
cool version (int32) and the ICE corrected matrix of the hicCorrectMatrix
tests. The MAD cases use the coverage vectors filter_by_zscore computes
from them. A hiCMatrix has no value representation, so the cases that take
one compare what the function returns and, where it changes the object,
the state the drivers return.

correct_subparser is compared as its help text (with a fixed prog) and the
attributes of every action; plot_total_contact_dist as the PNG it writes,
the matrix it leaves in the hiCMatrix and its log records.
"""
import hashlib
import os

import numpy as np
import scipy.sparse as sp

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
C = 'hicexplorer.hicCorrectMatrix.'


def load(name='small_test_matrix_50kb_res.h5', biased=False, dtype=None):
    from hicmatrix import HiCMatrix as hm
    hic = hm.hiCMatrix(os.path.join(TEST_DATA, name))
    if biased:
        rng = np.random.default_rng(3)
        bias = rng.uniform(0.2, 4.0, hic.matrix.shape[0])
        coo = hic.matrix.tocoo()
        coo.data = coo.data * bias[coo.row] * bias[coo.col] / 3.0
        hic.matrix = coo.tocsr()
    if dtype is not None:
        hic.matrix.data = hic.matrix.data.astype(dtype)
    return hic


class Holder(object):
    """An object with a matrix attribute, what fill_gaps reads."""

    def __init__(self, matrix):
        self.matrix = matrix

    def __repr__(self):
        return 'Holder({!r}, nnz={})'.format(self.matrix.shape, self.matrix.nnz)


def coverage(hic):
    row_sum = np.asarray(hic.matrix.sum(axis=1)).flatten()
    return row_sum - hic.matrix.diagonal()


def namespace(**kw):
    import argparse
    return argparse.Namespace(**kw)


# ---------------------------------------------------------------------------
# drivers

def parser_summary(target):
    parser = target()
    parser.prog = 'hicCorrectMatrix'
    actions = [(a.option_strings, a.dest, a.nargs, a.const, a.default,
                None if a.type is None else a.type.__name__, a.choices, a.required, a.help, a.metavar,
                type(a).__name__, getattr(a, 'version', None)) for a in parser._actions]
    groups = [(g.title, [a.dest for a in g._group_actions]) for g in parser._action_groups]
    return parser.format_help(), actions, groups, parser.add_help


def mad_summary(target, points):
    mad = target(points)
    out = [mad.mad_b_value, mad.median, mad.med_abs_deviation, mad.modified_z_score,
           mad.get_motified_zscores(), mad.is_outlier(-1.5, 5.0), mad.is_outlier(-0.3, 0.3)]
    for value in (0.0, 3, 12.5, np.float64(7.0), np.array([0.0, 100.0]), np.array([[1.0], [2.0]]),
                  np.arange(4)):
        out.append(mad.value_to_mad(value))
        out.append(mad.mad_to_value(value))
    out.append(mad.value_to_mad([1.0, 2.0]))
    return out


def mad_list_to_value(target, points):
    return target(points).mad_to_value([1.0, 2.0])


def mad_zero_deviation(target):
    mad = target(np.array([5.0, 5.0, 5.0, 1.0]))
    return mad.med_abs_deviation, mad.value_to_mad(7.0), mad.value_to_mad(np.array([1.0, 9.0]))


def mad_changed_attributes(target, points):
    mad = target(points)
    mad.median = np.float64(2.5)
    mad.med_abs_deviation = np.float64(4.0)
    mad.modified_z_score = np.array([-3.0, 0.0, 3.0])
    return mad.value_to_mad(1.0), mad.mad_to_value(1.0), mad.is_outlier(-1.0, 1.0), mad.get_motified_zscores()


def filter_driver(target, hic, *args, **kwargs):
    return target(hic, *args, **kwargs)


def plot_driver(target, hic, args):
    target(hic, args)
    with open(args.plotName, 'rb') as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    return digest, hic.matrix


def plot_exit(target, hic, args):
    try:
        return plot_driver(target, hic, args)
    except SystemExit as error:
        return 'SystemExit', error.code, os.path.exists(args.plotName)


# ---------------------------------------------------------------------------
# inputs

def fill_args(name='small_test_matrix_50kb_res.h5', biased=False, dtype=None, bins=None, as_list=False):
    def make():
        hic = load(name, biased, dtype)
        failed = np.asarray(hic.nan_bins if bins is None else bins, dtype=np.int64)
        return (Holder(hic.matrix), list(failed) if as_list else failed)
    return make


def symmetric_block(dtype=np.float64):
    hic = load(biased=True)
    m = hic.matrix[:300, :300]
    m = sp.csr_matrix((m + m.T) / 2.0)
    m.data = m.data.astype(dtype)
    return m


def zeroed():
    hic = load()
    hic.matrix.data[:] = 0
    return hic


def zero_bins():
    hic = load()
    return [i for i in range(hic.matrix.shape[0]) if hic.matrix[i].nnz == 0]


NAN_BINS_SOME = [3, 7, 8, 9, 40, 41, 100, 250, 600, 601, 602, 1500]

CASES = [
    # correct_subparser ---------------------------------------------------------
    Case('correct_subparser', C + 'correct_subparser', driver='parser_summary'),
    # iterative_correction ------------------------------------------------------
    Case('iterative_correction/real', C + 'iterative_correction',
         args=lambda: (symmetric_block(), namespace(iterNum=50, verbose=False))),
    Case('iterative_correction/iter_num', C + 'iterative_correction',
         args=lambda: (symmetric_block(), namespace(iterNum=3, verbose=False))),
    Case('iterative_correction/int', C + 'iterative_correction',
         args=lambda: (sp.csr_matrix(np.rint(symmetric_block()).astype(np.int64)), namespace(iterNum=20, verbose=0))),
    Case('iterative_correction/missing_attribute', C + 'iterative_correction',
         args=lambda: (symmetric_block(), namespace(verbose=False))),
    # fill_gaps ------------------------------------------------------------------
    Case('fill_gaps/nan_bins_int64', C + 'fill_gaps', args=fill_args(bins=NAN_BINS_SOME)),
    Case('fill_gaps/nan_bins_biased', C + 'fill_gaps', args=fill_args(biased=True, bins=NAN_BINS_SOME)),
    Case('fill_gaps/int32_cool', C + 'fill_gaps', args=fill_args('small_test_matrix_50kb_res.cool', bins=NAN_BINS_SOME)),
    Case('fill_gaps/float32', C + 'fill_gaps', args=fill_args(biased=True, dtype=np.float32, bins=NAN_BINS_SOME)),
    Case('fill_gaps/real_nan_bins', C + 'fill_gaps', args=fill_args(biased=True)),
    Case('fill_gaps/first_and_single', C + 'fill_gaps', args=fill_args(biased=True, bins=[0, 5, 20])),
    Case('fill_gaps/runs_only', C + 'fill_gaps', args=fill_args(bins=[10, 11, 12, 30, 31])),
    Case('fill_gaps/empty', C + 'fill_gaps', args=fill_args(bins=[])),
    Case('fill_gaps/last_bin_index_error', C + 'fill_gaps', args=fill_args(bins=[5, 2793])),
    Case('fill_gaps/contiguous_index_error', C + 'fill_gaps', args=fill_args(bins=NAN_BINS_SOME),
         kwargs={'fill_contiguous': True}),
    Case('fill_gaps/contiguous_truthy_is_not_true', C + 'fill_gaps', args=fill_args(bins=NAN_BINS_SOME),
         kwargs={'fill_contiguous': 1}),
    Case('fill_gaps/list_type_error', C + 'fill_gaps', args=fill_args(bins=[3, 7, 8], as_list=True)),
    Case('fill_gaps/ice_corrected', C + 'fill_gaps',
         args=fill_args('hicCorrectMatrix/small_test_matrix_ICEcorrected_chrUextra_chr3LHet.h5', bins=[4, 9, 10, 300])),
    # MAD ------------------------------------------------------------------------
    Case('MAD/coverage_int64', C + 'MAD', args=lambda: (coverage(load()),), driver='mad_summary'),
    Case('MAD/coverage_float64', C + 'MAD', args=lambda: (coverage(load(biased=True)),), driver='mad_summary'),
    Case('MAD/coverage_int32', C + 'MAD', args=lambda: (coverage(load('small_test_matrix_50kb_res.cool')),),
         driver='mad_summary'),
    Case('MAD/column', C + 'MAD', args=lambda: (coverage(load())[:, None],), driver='mad_summary'),
    Case('MAD/three_columns', C + 'MAD', args=lambda: (np.random.default_rng(1).gamma(2, 3, (200, 3)),),
         driver='mad_summary'),
    Case('MAD/with_nan', C + 'MAD', args=lambda: (np.array([1.0, np.nan, 3.0, 4.0, 10.0]),), driver='mad_summary'),
    Case('MAD/odd_even', C + 'MAD', args=lambda: (np.array([4.0, 1.0, 2.0, 8.0, 0.0, -1.0]),), driver='mad_summary'),
    Case('MAD/empty', C + 'MAD', args=lambda: (np.array([], dtype=np.float64),), driver='mad_summary'),
    Case('MAD/no_positive', C + 'MAD', args=lambda: (np.array([0.0, -1.0, -2.0]),), driver='mad_summary'),
    Case('MAD/zero_deviation', C + 'MAD', driver='mad_zero_deviation'),
    Case('MAD/changed_attributes', C + 'MAD', args=lambda: (np.array([1.0, 2.0, 3.0]),),
         driver='mad_changed_attributes'),
    Case('MAD/list_points', C + 'MAD', args=([1.0, 2.0, 3.0],)),
    Case('MAD/mad_to_value_list', C + 'MAD', args=lambda: (np.array([1.0, 2.0, 3.0]),), driver='mad_list_to_value'),
    # filter_by_zscore -------------------------------------------------------------
    Case('filter_by_zscore/genome', C + 'filter_by_zscore', args=lambda: (load(), -1.5, 5.0),
         driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/genome_biased', C + 'filter_by_zscore', args=lambda: (load(biased=True), -1.0, 3.0),
         driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/perchr', C + 'filter_by_zscore', args=lambda: (load(), -1.5, 5.0),
         kwargs={'perchr': True}, driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/perchr_biased', C + 'filter_by_zscore', args=lambda: (load(biased=True), -2.0, 2.0),
         kwargs={'perchr': True}, driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/perchr_cool_int32', C + 'filter_by_zscore',
         args=lambda: (load('small_test_matrix_50kb_res.cool'), -1.5, 5.0), kwargs={'perchr': True},
         driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/perchr_ice', C + 'filter_by_zscore',
         args=lambda: (load('hicCorrectMatrix/small_test_matrix_ICEcorrected_chrUextra_chr3LHet.h5'), -1.5, 5.0),
         kwargs={'perchr': True}, driver='filter_driver', compare_args=False),
    Case('filter_by_zscore/nothing_removed', C + 'filter_by_zscore', args=lambda: (load(), -1000.0, 1000.0),
         kwargs={'perchr': True}, driver='filter_driver', compare_args=False),
    # plot_total_contact_dist ------------------------------------------------------
    Case('plot/genome', C + 'plot_total_contact_dist',
         args=lambda: (load(), namespace(perchr=False, xMax=None, plotName=File('plot.png'))),
         driver='plot_driver', compare_args=False),
    Case('plot/genome_xmax', C + 'plot_total_contact_dist',
         args=lambda: (load(), namespace(perchr=False, xMax=5000, plotName=File('plot.png'))),
         driver='plot_driver', compare_args=False),
    Case('plot/perchr', C + 'plot_total_contact_dist',
         args=lambda: (load(), namespace(perchr=True, xMax=None, plotName=File('plot.png'))),
         driver='plot_driver', compare_args=False),
    Case('plot/perchr_biased', C + 'plot_total_contact_dist',
         args=lambda: (load(biased=True), namespace(perchr=True, xMax=None, plotName=File('plot.png'))),
         driver='plot_driver', compare_args=False),
    Case('plot/empty_exit', C + 'plot_total_contact_dist',
         args=lambda: (zeroed(), namespace(perchr=False, xMax=None, plotName=File('plot.png'))),
         driver='plot_exit', compare_args=False),
]
