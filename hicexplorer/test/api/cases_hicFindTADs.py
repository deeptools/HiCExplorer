"""Equivalence cases for hicexplorer.hicFindTADs (apieq.py).

The pipeline cases run HicFindTads the way main() of 3.7.6 did (construct,
compute_spectra_matrix, find_boundaries, save the bedgraph matrix, the
domains and the boundaries) on real matrices of the test data, and compare
the instance's state, the state of its hiCMatrix and every file written. The
other cases call the module functions and the static methods directly, on
real data and on the edge cases the real data does not reach.

A hiCMatrix given to a function is loaded in each interpreter with
hicmatrix (the input), and its state after the call is compared through
hic_state(), which lists the attributes the 3.7.6 code reads and writes.

Deviations from 3.7.6, kept on purpose (the C++ side documents them in
cpp/python/src/api/hicFindTADs.cpp):

* HicFindTads(<file name>) loads the matrix with the C++ matrix layer into a
  hiCMatrix whose matrixFileHandler is None, so a later hic_ma.save() picks
  its writer from the output name instead of from the input.
* get_incremental_step_size raises ValueError where 3.7.6 never returns
  (start_step_len <= 0 with min_win_size <= max_win_size, a NaN argument or
  an infinite maximum); no case can pin an endless loop.
* Window lengths that are not whole numbers, y_axis arrays of float32 and
  matrices that are neither integer nor float raise TypeError naming
  HiCExplorer 4.
* compute_spectra_matrix uses num_processors threads instead of processes.
  The result does not depend on the count, and a count that leaves a worker
  without a single usable bin still raises the ValueError of 3.7.6
  ('not enough values to unpack', pinned by too_many_processors).
* hiCMatrix matrices are compared with sorted column indices (sorted_matrix);
  the order of the stored entries within a row is not compared.
* The log records of hicmatrix's own logger (hicmatrix.HiCMatrix) are not
  part of the comparison, which covers the hicexplorer loggers.
* The CTCF filter of the hicFindTADs tool (--ctcfFile) is a HiCExplorer 4
  extension with the divergent orientation rule ('-' before, '+' after a
  boundary); no API name filters by CTCF.
"""
import argparse
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicFindTADs.'
T = M + 'HicFindTads'

SMALL_50KB = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
SMALL_50KB_COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
SMALL = os.path.join(TEST_DATA, 'small_test_matrix.h5')


# ---------------------------------------------------------------------------
# state

def sorted_matrix(matrix):
    """A copy with sorted column indices: scipy's fancy indexing leaves the
    columns of a row in the order of the source row, which is not part of the
    comparison (the entries, explicit zeros included, are)."""
    if matrix is None:
        return None
    copy = matrix.copy()
    copy.sort_indices()
    return copy


def hic_state(h):
    trees = getattr(h, 'interval_trees', None)
    return {
        'matrix': sorted_matrix(h.matrix),
        'cut_intervals': list(h.cut_intervals),
        'nan_bins': h.nan_bins,
        'correction_factors': h.correction_factors,
        'distance_counts': h.distance_counts,
        'orig_bin_ids': h.orig_bin_ids,
        'orig_cut_intervals': list(h.orig_cut_intervals),
        'chrBinBoundaries': getattr(h, 'chrBinBoundaries', None),
        'interval_trees': None if trees is None else
        [(k, sorted((iv.begin, iv.end, iv.data) for iv in v)) for k, v in trees.items()],
        'bin_size': h.bin_size,
    }


ATTRIBUTES = ('num_processors', 'max_depth', 'min_depth', 'step', 'delta', 'min_boundary_distance',
              'use_zscore', 'binsize', 'correct_for_multiple_testing', 'threshold_comparisons', 'pvalueFDR',
              'bedgraph_matrix', 'boundaries')


def tads_state(t):
    out = {name: getattr(t, name, '<missing>') for name in ATTRIBUTES}
    out['hic_ma'] = hic_state(t.hic_ma)
    return out


def files_of(directory):
    out = {}
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if os.path.isfile(path):
            with open(path, 'rb') as handle:
                out[name] = handle.read()
    return out


def load(path):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


# ---------------------------------------------------------------------------
# drivers

def exits(function):
    """exit() of the 3.7.6 code raises SystemExit, which apieq does not catch:
    the drivers return its code instead."""
    def driver(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except SystemExit as error:
            return ('SystemExit', error.code)
    driver.__name__ = function.__name__
    return driver


@exits
def pipeline(cls, directory, matrix, kwargs, as_object=False, perchr=True):
    """main() of 3.7.6: compute the TAD-separation score, find and save the boundaries."""
    import hicexplorer.hicFindTADs as module
    if as_object:
        matrix = load(matrix)
    t = cls(matrix, **kwargs)
    after_init = tads_state(t)
    t.compute_spectra_matrix(perchr=perchr)
    t.save_bedgraph_matrix(os.path.join(directory, 'out_tad_score.bm'))
    t.find_boundaries()
    t.save_domains_and_boundaries(os.path.join(directory, 'out'))
    return {'init': after_init, 'final': tads_state(t), 'files': files_of(directory),
            'global_is_hic_ma': module.hic_ma is t.hic_ma}


@exits
def construct(cls, matrix, kwargs, as_object=False):
    if as_object:
        matrix = load(matrix)
    return tads_state(cls(matrix, **kwargs))


@exits
def spectra(cls, matrix, kwargs, as_object=False):
    if as_object:
        matrix = load(matrix)
    t = cls(matrix, **kwargs)
    t.compute_spectra_matrix()
    return tads_state(t)


@exits
def from_bedgraph(cls, directory, matrix, bm_name, kwargs, chromosomes=None):
    """The --TAD_sep_score_prefix route: a saved score is loaded and the boundaries called."""
    t = cls(matrix, **kwargs)
    t.load_bedgraph_matrix(os.path.join(directory, bm_name), pChromosomes=chromosomes)
    state = tads_state(t)
    return state


@exits
def bedgraph_roundtrip(cls, directory, matrix, kwargs):
    t = cls(matrix, **kwargs)
    t.compute_spectra_matrix()
    path = os.path.join(directory, 'score.bm')
    t.save_bedgraph_matrix(path)
    t.bedgraph_matrix = None
    t.load_bedgraph_matrix(path)
    t.find_boundaries()
    t.save_domains_and_boundaries(os.path.join(directory, 'again'))
    return {'loaded': tads_state(t), 'files': files_of(directory)}


@exits
def with_bedgraph(cls, directory, matrix, kwargs, bedgraph, boundaries, attributes, method, *args):
    """Calls one method on an instance whose bedgraph matrix and boundaries are
    those the named input functions return."""
    t = cls(matrix, **kwargs)
    if bedgraph is not None:
        t.bedgraph_matrix = globals()[bedgraph]()
    if boundaries is not None:
        t.boundaries = globals()[boundaries]()
    for name, value in attributes.items():
        setattr(t, name, value)
    result = getattr(t, method)(*[os.path.join(directory, a[5:]) if isinstance(a, str) and a.startswith('FILE:')
                                  else a for a in args])
    return {'result': result, 'state': tads_state(t), 'files': files_of(directory)}


def module_function(target, matrix_path, *args):
    """A module function called with a hiCMatrix loaded from `matrix_path`."""
    h = load(matrix_path)
    result = target(h, *args)
    return {'result': result, 'hic': hic_state(h)}


def module_function_matrix(target, matrix_path, *args):
    """A module function called with the scipy matrix of `matrix_path`."""
    h = load(matrix_path)
    return target(h.matrix, *args)


def compute_with_global(target, matrix_path, bins, *args):
    import hicexplorer.hicFindTADs as module
    module.hic_ma = load(matrix_path)
    try:
        return target(list(bins), *args)
    finally:
        module.hic_ma = None


def compute_wrapper(target, matrix_path, args):
    import hicexplorer.hicFindTADs as module
    module.hic_ma = load(matrix_path)
    try:
        return target(args)
    finally:
        module.hic_ma = None


def call_static(cls, name, *args, **kwargs):
    return getattr(cls, name)(*args, **kwargs)


def clustering(cls, *args):
    z, clusters = cls.hierarchical_clustering(None, *args)
    return z, clusters


def save_linkage(cls, directory, Z, name):
    result = cls.save_linkage(Z, os.path.join(directory, name))
    return result, files_of(directory)


def save_clusters(cls, directory, clusters):
    result = cls.save_clusters(clusters, os.path.join(directory, 'clusters'))
    return result, files_of(directory)


def print_args(target, values):
    return target(argparse.Namespace(**values))


# ---------------------------------------------------------------------------
# inputs

def score_matrix(rows=60, columns=4, seed=3):
    rng = np.random.RandomState(seed)
    base = np.sin(np.arange(rows) * 0.45)[:, None] + 0.1 * rng.randn(rows, columns)
    return base


def chromosomes(rows=60, split=35):
    return np.array(['chr2'] * split + ['chr1'] * (rows - split))


def score_mean(rows=60, seed=3):
    return score_matrix(rows, seed=seed).mean(axis=1)


BM = (b'#{"step":100000,"minDepth":150000,"maxDepth":500000,"binsize":50000}\n'
      + b''.join('chr{}\t{}\t{}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(
          'X' if i < 30 else '2L', (i % 30) * 50000, (i % 30 + 1) * 50000,
          np.sin(i * 0.5), np.cos(i * 0.3), np.sin(i * 0.2) + 0.1).encode() for i in range(60)))

# chromosome names out of sorted order: np.unique sorts them, np.split then
# cuts the starts at unsorted points and the assert of 3.7.6 fails (pinned)
BOUNDARIES_UNSORTED = [('chr1', 0, 1.0), ('chr1', 100, 0.5), ('chr1', 250, 0.9), ('chr1', 400, 0.2),
                       ('chr1', 600, 1.0), ('chrX', 50, 0.3), ('chrX', 300, 0.8), ('chrX', 520, 0.4),
                       ('chrX', 700, 0.7), ('chr10', 0, 0.1), ('chr10', 90, 0.6), ('chr10', 180, 0.6)]
BOUNDARIES = [('chr1', 0, 1.0), ('chr1', 100, 0.5), ('chr1', 250, 0.9), ('chr1', 400, 0.2),
              ('chr1', 600, 1.0), ('chr2', 50, 0.3), ('chr2', 300, 0.8), ('chr2', 520, 0.4),
              ('chr2', 700, 0.7), ('chr2', 810, 0.75), ('chr3', 0, 0.1), ('chr3', 90, 0.6),
              ('chr3', 180, 0.6), ('chr3', 260, 0.35)]

KW = {'min_depth': 150000, 'max_depth': 500000, 'step': 100000, 'delta': 0.01}


def small_bedgraph():
    return {'chrom': np.array(['chrA'] * 12 + ['chrB'] * 12),
            'chr_start': np.array([i * 100 for i in range(12)] * 2),
            'chr_end': np.array([i * 100 + 100 for i in range(12)] * 2),
            'matrix': score_matrix(24, 3, seed=5)}


def small_boundaries():
    return {'min_idx': (3, 6, 15, 19), 'delta': {3: 0.2, 6: np.float64(0.005), 15: 0.4},
            'pvalues': {3: np.float64(0.001), 6: np.float64(0.001), 15: np.float64(0.2), 19: 0.0}}


CASES = [
    # the pipeline, real matrices
    Case('pipeline_50kb_fdr', T, args=(File(), SMALL_50KB, dict(KW, p_correct_for_multiple_testing='fdr',
                                                                p_threshold_comparisons=0.05)),
         driver='pipeline'),
    Case('pipeline_50kb_default_depths', T, args=(File(), SMALL_50KB, {}), driver='pipeline'),
    Case('pipeline_50kb_bonferroni', T,
         args=(File(), SMALL_50KB, dict(KW, p_correct_for_multiple_testing='bonferroni',
                                        p_threshold_comparisons=0.5)),
         driver='pipeline'),
    Case('pipeline_50kb_none', T,
         args=(File(), SMALL_50KB, dict(KW, p_correct_for_multiple_testing='None', p_threshold_comparisons=0.3,
                                        delta=0.0, min_boundary_distance=100000)),
         driver='pipeline'),
    Case('pipeline_50kb_cool', T, args=(File(), SMALL_50KB_COOL, dict(KW, num_processors=3)), driver='pipeline'),
    Case('pipeline_50kb_object_chromosomes', T,
         args=(File(), SMALL_50KB, dict(KW, pChromosomes=['chr3R', 'chrX', 'chrNone'])),
         kwargs={'as_object': True}, driver='pipeline'),
    Case('pipeline_50kb_no_zscore', T, args=(File(), SMALL_50KB, dict(KW, use_zscore=False)),
         driver='pipeline'),
    Case('pipeline_50kb_genome_wide_zscore', T, args=(File(), SMALL_50KB, KW), kwargs={'perchr': False},
         driver='pipeline'),
    Case('pipeline_small_10kb_cool', T,
         args=(File(), os.path.join(TEST_DATA, 'small_test_matrix.cool'),
               {'min_depth': 30000, 'max_depth': 100000, 'step': 10000, 'num_processors': 4,
                'pChromosomes': ['chrX', 'chr3L']}),
         driver='pipeline'),
    Case('too_many_processors', T, args=(SMALL_50KB, dict(KW, num_processors=400)), driver='spectra'),
    Case('bedgraph_roundtrip', T, args=(File(), SMALL_50KB, KW), driver='bedgraph_roundtrip'),
    # construction and set_variables
    Case('init_object_defaults', T, args=(SMALL_50KB, {}), kwargs={'as_object': True}, driver='construct'),
    Case('init_max_not_above_min', T, args=(SMALL_50KB, {'min_depth': 500000, 'max_depth': 500000}),
         driver='construct'),
    Case('init_max_depth_too_small', T, args=(SMALL_50KB, {'max_depth': 100000}), driver='construct'),
    Case('init_min_depth_too_small', T, args=(SMALL_50KB, {'min_depth': 100000}), driver='construct'),
    Case('init_step_too_small', T, args=(SMALL_50KB, {'step': 10000}), driver='construct'),
    Case('init_unknown_chromosome_only', T, args=(SMALL_50KB, {'pChromosomes': ['chrX', 'nope']}),
         driver='construct'),
    # load_bedgraph_matrix
    Case('load_bedgraph', T, args=(File(), SMALL_50KB, 'a.bm', KW), files={'a.bm': BM}, driver='from_bedgraph'),
    Case('load_bedgraph_chromosomes', T, args=(File(), SMALL_50KB, 'a.bm', KW, ['chr2L']), files={'a.bm': BM},
         driver='from_bedgraph'),
    Case('load_bedgraph_without_header', T, args=(File(), SMALL_50KB, 'a.bm', KW),
         files={'a.bm': BM.split(b'\n', 1)[1]}, driver='from_bedgraph'),
    Case('load_bedgraph_ragged', T, args=(File(), SMALL_50KB, 'a.bm', KW),
         files={'a.bm': BM + b'chrX\t1\t2\t0.5\n'}, driver='from_bedgraph'),
    Case('load_bedgraph_bad_number', T, args=(File(), SMALL_50KB, 'a.bm', KW),
         files={'a.bm': BM + b'chrX\t1\t2\tx\t1\t1\n'}, driver='from_bedgraph'),
    Case('load_bedgraph_missing_file', T, args=(File(), SMALL_50KB, 'missing.bm', KW), driver='from_bedgraph'),
    # methods on a given state
    Case('save_domains_given_state_fdr', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', 'small_boundaries', {'pvalueFDR': 0.01, 'delta': 0.01},
               'save_domains_and_boundaries', 'FILE:x'),
         driver='with_bedgraph'),
    Case('save_domains_given_state_bonferroni', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', 'small_boundaries',
               {'correct_for_multiple_testing': 'bonferroni', 'threshold_comparisons': 0.5, 'delta': 0.0},
               'save_domains_and_boundaries', 'FILE:x'),
         driver='with_bedgraph'),
    Case('save_domains_given_state_unknown_mode', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', 'small_boundaries',
               {'correct_for_multiple_testing': 'holm'}, 'save_domains_and_boundaries', 'FILE:x'),
         driver='with_bedgraph'),
    Case('save_domains_without_pvalue_fdr', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', 'small_boundaries', {}, 'save_domains_and_boundaries',
               'FILE:x'),
         driver='with_bedgraph'),
    Case('save_bedgraph_numpy_step', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', None, {'step': np.int64(5)}, 'save_bedgraph_matrix',
               'FILE:x.bm'),
         driver='with_bedgraph'),
    Case('save_bedgraph_float_depth', T,
         args=(File(), SMALL_50KB, KW, 'small_bedgraph', None, {'min_depth': 1.5e5, 'binsize': None},
               'save_bedgraph_matrix', 'FILE:x.bm'),
         driver='with_bedgraph'),
    Case('find_boundaries_given_scores', T,
         args=(File(), SMALL_50KB, KW, None, None, {}, 'load_bedgraph_matrix', 'FILE:a.bm'),
         files={'a.bm': BM}, driver='with_bedgraph'),
    # static methods
    Case('peakdetect_default', T, args=('peakdetect', score_mean()), kwargs={'chrom': chromosomes()},
         driver='call_static'),
    Case('peakdetect_x_axis_delta', T, args=('peakdetect', score_mean(), np.arange(60) * 10 + 5),
         kwargs={'lookahead': 2, 'delta': 0.05, 'chrom': chromosomes()}, driver='call_static'),
    Case('peakdetect_lookahead_zero', T, args=('peakdetect', score_mean()),
         kwargs={'lookahead': 0, 'chrom': chromosomes()}, driver='call_static'),
    Case('peakdetect_lookahead_long', T, args=('peakdetect', score_mean()),
         kwargs={'lookahead': 9, 'chrom': chromosomes()}, driver='call_static'),
    Case('peakdetect_negative_delta', T, args=('peakdetect', score_mean()),
         kwargs={'delta': -1, 'chrom': chromosomes()}, driver='call_static'),
    Case('peakdetect_length_mismatch', T, args=('peakdetect', score_mean(), np.arange(5)),
         kwargs={'chrom': chromosomes()}, driver='call_static'),
    Case('peakdetect_no_chrom', T, args=('peakdetect', score_mean()), driver='call_static'),
    Case('peakdetect_list', T, args=('peakdetect', [0.5, 0.2, 0.1, 0.3, 0.6, 0.4, 0.1, 0.0, 0.2]),
         kwargs={'chrom': ['c'] * 9, 'lookahead': 2}, driver='call_static'),
    Case('peakdetect_infinity', T, args=('peakdetect', np.array([0.1, 0.5, np.inf, 0.2, 0.3, 0.1])),
         kwargs={'chrom': ['c'] * 6, 'lookahead': 1}, driver='call_static'),
    Case('delta_wrt_window', T, args=('delta_wrt_window', [3, 12, 20, 33, 36, 47, 58], score_mean(),
                                      chromosomes()), driver='call_static'),
    Case('delta_wrt_window_numpy_keys', T,
         args=('delta_wrt_window', tuple(np.array([12, 20, 47])), score_mean(), chromosomes()),
         kwargs={'window_len': 4}, driver='call_static'),
    Case('find_consensus_minima', T, args=('find_consensus_minima', score_matrix()),
         kwargs={'chrom': chromosomes()}, driver='call_static'),
    Case('find_consensus_minima_none', T, args=('find_consensus_minima', np.ones((10, 3))),
         kwargs={'chrom': ['c'] * 10, 'lookahead': 2}, driver='call_static'),
    Case('get_domains', T, args=('get_domains', [('a', 0, 1), ('a', 10, 2), ('a', None, 3), ('a', 30, 4),
                                                 ('b', 0, 5), ('b', 40, 6)]), driver='call_static'),
    Case('get_domains_empty', T, args=('get_domains', []), driver='call_static'),
    Case('clustering', T, args=(BOUNDARIES, [0.3, 0.6, 0.6]), driver='clustering'),
    Case('clustering_one_cutoff', T, args=(BOUNDARIES, np.array([0.5])), driver='clustering'),
    Case('clustering_unsorted_chromosomes', T, args=(BOUNDARIES_UNSORTED, [0.3]), driver='clustering'),
    Case('clustering_no_cutoff', T, args=(BOUNDARIES,), driver='clustering'),
    Case('clustering_integer_values', T,
         args=([('c', 0, 3), ('c', 10, 1), ('c', 30, 2), ('c', 45, 5), ('c', 70, 1)], [2]),
         driver='clustering'),
    Case('clustering_empty', T, args=([],), driver='clustering'),
    Case('save_linkage', T, args=(File(), {'chr1': [(1, 2, 0.9, 2, 175.5, 325), (0, 4, np.float64(0.5), 3, 50, 250)],
                                           'chrX': [(0, 1, 1, 2, 10, 20)]}, 'z.txt'), driver='save_linkage'),
    Case('save_linkage_unwritable', T, args=(File(), {'chr1': []}, 'missing/z.txt'), driver='save_linkage'),
    Case('save_linkage_short_tuple', T, args=(File(), {'chr1': [(1, 2, 0.9, 2, 175, 325), (1, 2)]}, 'z.txt'),
         driver='save_linkage'),
    Case('save_clusters', T, args=(File(), {'0.3': [('chr1', 0, 100), ('chr1', 100, 400)], '0.6': []}),
         driver='save_clusters'),
    # module functions
    Case('step_sizes_ints', M + 'get_incremental_step_size', args=(8, 50, 2)),
    Case('step_sizes_bp', M + 'get_incremental_step_size', args=(150000, 500000, 100000)),
    Case('step_sizes_floats', M + 'get_incremental_step_size', args=(2.5, 30.0, 1.5)),
    Case('step_sizes_numpy', M + 'get_incremental_step_size', args=(np.int64(10), 40, np.int64(3))),
    Case('step_sizes_empty', M + 'get_incremental_step_size', args=(10, 5, 2)),
    Case('step_sizes_nan_step', M + 'get_incremental_step_size', args=(10, 50, float('nan'))),
    Case('cut_weight_mean', M + 'get_cut_weight', args=(SMALL_50KB, 30, 200000, True), driver='module_function'),
    Case('cut_weight_values', M + 'get_cut_weight', args=(SMALL_50KB, 30, 150000), driver='module_function'),
    Case('cut_weight_first_bin', M + 'get_cut_weight', args=(SMALL_50KB, 0, 150000, True),
         driver='module_function'),
    Case('cut_weight_past_end', M + 'get_cut_weight', args=(SMALL_50KB, 10000, 150000, True),
         driver='module_function'),
    Case('cut_weight_negative', M + 'get_cut_weight', args=(SMALL_50KB, -1, 150000, True),
         driver='module_function'),
    Case('cut_weight_by_bin_mean', M + 'get_cut_weight_by_bin_id', args=(SMALL_50KB, 40, 6, True),
         driver='module_function_matrix'),
    Case('cut_weight_by_bin_values', M + 'get_cut_weight_by_bin_id', args=(SMALL_50KB, 40, 4),
         driver='module_function_matrix'),
    Case('cut_weight_by_bin_edge', M + 'get_cut_weight_by_bin_id', args=(SMALL_50KB, 0, 4, True),
         driver='module_function_matrix'),
    Case('cut_weight_by_bin_none', M + 'get_cut_weight_by_bin_id', args=(SMALL_50KB, -3, 4, True),
         driver='module_function_matrix'),
    Case('bins_at_distance', M + 'get_idx_of_bins_at_given_distance', args=(SMALL_50KB, 50, 300000),
         driver='module_function'),
    Case('bins_at_distance_chromosome_start', M + 'get_idx_of_bins_at_given_distance',
         args=(SMALL_50KB, 66, 300000), driver='module_function'),
    Case('bins_at_distance_past_end', M + 'get_idx_of_bins_at_given_distance', args=(SMALL_50KB, 99999, 50000),
         driver='module_function'),
    Case('triangle', M + 'get_triangle', args=(SMALL_50KB, 50, 200000), driver='module_function'),
    Case('compute_matrix', M + 'compute_matrix', args=(SMALL_50KB, range(0, 140, 3), 150000, 500000, 100000),
         driver='compute_with_global'),
    Case('compute_matrix_nothing_usable', M + 'compute_matrix', args=(SMALL_50KB, [-1, -2], 150000, 500000, 100000),
         driver='compute_with_global'),
    Case('compute_matrix_wrapper', M + 'compute_matrix_wrapper', args=(SMALL_50KB, ([10, 20, 30], 100000, 300000, 50000)),
         driver='compute_wrapper'),
    Case('print_args', M + 'print_args', args=({'matrix': 'm.h5', 'delta': 0.01, 'chromosomes': None},),
         driver='print_args'),
]


# Cases whose floats differ from 3.7.6 in the last bits since the special
# functions (normal tail, incomplete beta and gamma functions, log-gamma) are
# Boost.Math's and no longer the routines scipy uses: the p-values agree
# to a relative 2e-14 or better. The cases compare every float at rtol 1e-12;
# integers, strings, coordinates and the structure of the result are still
# compared exactly, so a changed call, a moved TAD or a different row count
# fails.
_BOOST_LAST_BITS = {
    'bedgraph_roundtrip',
    'pipeline_50kb_bonferroni',
    'pipeline_50kb_cool',
    'pipeline_50kb_default_depths',
    'pipeline_50kb_fdr',
    'pipeline_50kb_genome_wide_zscore',
    'pipeline_50kb_no_zscore',
    'pipeline_50kb_none',
    'pipeline_50kb_object_chromosomes',
    'pipeline_small_10kb_cool',
}
for _case in CASES:
    if _case.id in _BOOST_LAST_BITS:
        _case.rtol = 1e-12
        _case.note = (_case.note + '; ' if _case.note else '') + (
            'Boost.Math special functions round the last bits of the p-values differently from scipy 1.14')
