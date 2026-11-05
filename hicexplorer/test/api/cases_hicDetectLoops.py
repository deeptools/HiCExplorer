"""Equivalence cases for hicexplorer.hicDetectLoops (apieq.py).

compute_loops runs on the matrices of the hicDetectLoops tests with the
options those tests use: small_test_matrix.h5 (chr2L, every --expected,
--pValuePreselection as a number and as a threshold file, the window
defaults), Li_et_al_2015.h5 (restriction fragment bins) and, with a queue,
GSE63525_GM12878_insitu_primary_2_5mb.cool, which compute_loops then loads
itself. The functions below it get the inputs compute_loops hands them,
built here with numpy from the chr2L block: the upper triangle without its
diagonal and an obs/exp matrix truncated to integers, as obs_exp_matrix
leaves an integer count matrix.

Functions that put their result into a queue are called with a queue.Queue
and compared as what they put (driver 'queued'); a 'Fail: ' string is
compared up to the traceback, whose frames are not part of the contract
(the C++ names its own). A hiCMatrix has no value representation, so the
cases that take one compare the results, the namespace after the call and
whether the hiCMatrix still has its matrix (compute_long_range_contacts
deletes it). The timing of the 3.7.6 queue polling (one second sleeps) is
not compared.
"""
import os

import numpy as np
import scipy.sparse as sp

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
D = 'hicexplorer.hicDetectLoops.'
MATRIX_H5 = os.path.join(TEST_DATA, 'small_test_matrix.h5')
MATRIX_COOL = os.path.join(TEST_DATA, 'hicDetectLoops', 'GSE63525_GM12878_insitu_primary_2_5mb.cool')
LI_H5 = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')
THRESHOLD_FILE = os.path.join(TEST_DATA, 'hicCreateThresholdFile', 'thresholdFile_loose_pValue.txt')


# ---------------------------------------------------------------------------
# inputs

def load(path=MATRIX_H5):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


def namespace(**kw):
    import argparse
    args = dict(matrix=MATRIX_H5, maxLoopDistance=2000000, windowSize=5, peakWidth=2, expected='mean',
                threadsPerChromosome=1, pValuePreselection=0.5, pValue=0.5, peakInteractionsThreshold=1,
                obsExpThreshold=1.5)
    args.update(kw)
    return argparse.Namespace(**args)


_BLOCK = {}


def chr2l_block(first=6486, size=600):
    """triu of a part of chr2L of small_test_matrix.h5 without its diagonal,
    the matrix compute_long_range_contacts sees."""
    key = (first, size)
    if key not in _BLOCK:
        import tables
        with tables.open_file(MATRIX_H5, 'r') as f:
            m = sp.csr_matrix((f.root.matrix.data.read(), f.root.matrix.indices.read(), f.root.matrix.indptr.read()),
                              shape=tuple(f.root.matrix.shape.read()))
        block = sp.triu(m[first:first + size, first:first + size], k=1, format='csr')
        block.eliminate_zeros()
        _BLOCK[key] = block
    return _BLOCK[key].copy()


def obs_exp_int(block):
    """count / (sum of the distance / (n + 1 - d)) at distance ceil(d / 2),
    truncated to the integer dtype, as obs_exp_matrix leaves an int matrix
    (a truncated zero is set to one, so that the structure stays)."""
    n = block.shape[0]
    coo = block.tocoo()
    d = np.abs(coo.row - coo.col)
    sums = np.bincount(d, weights=coo.data, minlength=n)
    expected = sums / np.arange(n + 1, 1, -1)
    lookup = np.ceil(d / 2).astype(np.int32)
    values = coo.data.astype(np.float32) / expected[lookup]
    values[~np.isfinite(values)] = 1e-6
    values = values.astype(block.dtype)
    # keep the structure of the count matrix, as compute_loops requires
    values[values == 0] = 1
    return sp.csr_matrix((values, (coo.row, coo.col)), shape=block.shape)


def obs_exp_float(block):
    out = obs_exp_int(block.astype(np.float64))
    out.data = out.data / 37.0
    return out


class Hic(object):
    """What compute_long_range_contacts and the mapping read of a hiCMatrix."""

    def __init__(self, matrix, first=6486):
        self.matrix = matrix
        self.cut_intervals = [('chr2L', i * 5000, (i + 1) * 5000, 1.0) for i in range(matrix.shape[0])]
        self.first = first

    def getBinSize(self):
        return 5000

    def getBinPos(self, binIndex):
        if binIndex < len(self.cut_intervals):
            return self.cut_intervals[binIndex]
        raise ValueError(f"binIndex: {binIndex} not found")

    def __repr__(self):
        return 'Hic({})'.format(self.matrix.shape if hasattr(self, 'matrix') else 'deleted')


def distribution_inputs(obs=None):
    obs = obs_exp_int(chr2l_block()) if obs is None else obs
    rows, cols = obs.nonzero()
    distances = np.absolute(rows - cols)
    return obs.data, distances


def distributions(obs=None):
    data, distances = distribution_inputs(obs)
    positions, obs_exp = {}, {}
    for distance in range(distances.min(), distances.max() + 1):
        mask = distances == distance
        obs_exp[distance] = data[mask]
        positions[distance] = np.argwhere(mask).flatten()
    return obs_exp, positions


def candidates_of(obs, top=400):
    order = np.argsort(-obs.data, kind='stable')[:top]
    rows, cols = obs.nonzero()
    return np.array([*zip(rows[order], cols[order])])


# ---------------------------------------------------------------------------
# drivers

def _strip(value):
    if isinstance(value, str) and 'Fail: ' in value:
        return value.split('Traceback')[0]
    if isinstance(value, (list, tuple)):
        return type(value)(_strip(v) for v in value)
    return value


def queued(target, *args):
    import queue
    q = queue.Queue()
    returned = target(*args, q)
    items = []
    while not q.empty():
        items.append(_strip(q.get()))
    return _strip(returned), items


def stripped(target, *args, **kwargs):
    return _strip(target(*args, **kwargs))


def loops(target, hic, region, args, is_cooler):
    result = target(hic, region, args, is_cooler)
    return _strip(result), vars(args), hasattr(hic, 'matrix')


def loops_queued(target, region, args, is_cooler):
    import queue
    q = queue.Queue()
    returned = target(None, region, args, is_cooler, q)
    items = []
    while not q.empty():
        items.append(_strip(q.get()))
    return returned, items, vars(args)


def long_range(target, hic, *args):
    result = target(hic, *args)
    return _strip(result), hasattr(hic, 'matrix')


def line_number(target):
    return target()


def bedgraph(target, loops_, path, *args):
    target(loops_, path, *args)
    with open(path) as handle:
        return handle.read()


def bedgraph_partial(target, loops_, path):
    try:
        target(loops_, path)
    except Exception as error:
        with open(path) as handle:
            return type(error).__name__, str(error), handle.read()


def mapping_hic():
    return load()


LOOPS = [('chr2L', 485000, 490000, 'chr2L', 490000, 495000, 0.4385623894062053),
         ('chr2L', 960000, 965000, 'chr2L', 1045000, 1050000, np.float64(0.4016704795074202)),
         ('chr3L', np.int64(5), np.int64(10), 'chr3L', np.int64(15), np.int64(20), 1e-5)]

THRESHOLDS = (b'# a comment\n# Standard threshold 0.5\n-5000\t0.5\n0\t0.25\n5000\t1e-3\n'
              b'10000\t0.125\n\n15000\t0.9\n')


def pvalue_inputs(keys=None, threshold=0.5, obs=None):
    def make():
        obs_exp, positions = distributions(obs)
        key_list = sorted(obs_exp) if keys is None else keys
        return (obs_exp, key_list, threshold, positions, 5000, 1, 1.5)
    return make


def long_range_inputs(threads=1, preselection=0.5, obs_kind='int', explicit_zero=False):
    def make():
        block = chr2l_block()
        obs = obs_exp_int(block) if obs_kind == 'int' else obs_exp_float(block)
        if explicit_zero:
            obs.data[5] = 0
        return (Hic(block), obs, 5, 0.5, 2, preselection, 1, 1.5, threads)
    return make


def merge_inputs(threads=1, top=400):
    def make():
        obs = obs_exp_int(chr2l_block())
        return (candidates_of(obs, top), 5, obs, threads)
    return make


def region_inputs(threads=1, peak=2, window=5):
    def make():
        obs = obs_exp_int(chr2l_block())
        return (obs, candidates_of(obs, 300), window, 0.5, peak, threads)
    return make


L = dict(driver='loops', compare_args=False)
CASES = [
    # get_linenumber --------------------------------------------------------------
    Case('get_linenumber', D + 'get_linenumber', driver='line_number'),
    # read_threshold_file ---------------------------------------------------------
    Case('read_threshold_file/real', D + 'read_threshold_file', args=(THRESHOLD_FILE,)),
    Case('read_threshold_file/constructed', D + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': THRESHOLDS}),
    Case('read_threshold_file/one_field', D + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b'5000\t0.5\n10000\n'}),
    Case('read_threshold_file/three_fields', D + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b'5000\t0.5\t1\n'}),
    Case('read_threshold_file/bad_float', D + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b'5000\tabc\n'}),
    Case('read_threshold_file/bad_int', D + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b'5e3\t0.5\n'}),
    Case('read_threshold_file/missing', D + 'read_threshold_file', args=(File('missing.txt'),)),
    # write_bedgraph ----------------------------------------------------------------
    Case('write_bedgraph/all', D + 'write_bedgraph', args=(LOOPS, File('o.bedgraph')), driver='bedgraph'),
    Case('write_bedgraph/region', D + 'write_bedgraph', args=(LOOPS, File('o.bedgraph'), 400000, 1100000),
         driver='bedgraph'),
    Case('write_bedgraph/region_zero_start', D + 'write_bedgraph', args=(LOOPS, File('o.bedgraph'), 0, 1100000),
         driver='bedgraph'),
    Case('write_bedgraph/empty', D + 'write_bedgraph', args=([], File('o.bedgraph')), driver='bedgraph'),
    Case('write_bedgraph/short_tuple', D + 'write_bedgraph', args=(LOOPS[:1] + [('a', 1, 2)], File('o.bedgraph')),
         driver='bedgraph_partial'),
    Case('write_bedgraph/long_tuple', D + 'write_bedgraph', args=([LOOPS[0] + (1,)], File('o.bedgraph')),
         driver='bedgraph_partial'),
    Case('write_bedgraph/list_item', D + 'write_bedgraph', args=([list(LOOPS[0])], File('o.bedgraph')),
         driver='bedgraph_partial'),
    Case('write_bedgraph/missing_directory', D + 'write_bedgraph', args=(LOOPS, File('no/such/o.bedgraph'))),
    # cluster_to_genome_position_mapping ---------------------------------------------
    Case('mapping/real', D + 'cluster_to_genome_position_mapping',
         args=lambda: (mapping_hic(), np.array([[6500, 6510], [6600, 6900], [100, 104]], dtype=np.int32),
                       np.array([0.1, 0.2, 0.3]), None), compare_args=False),
    Case('mapping/max_distance', D + 'cluster_to_genome_position_mapping',
         args=lambda: (mapping_hic(), np.array([[6500, 6510], [6600, 6900], [100, 104]], dtype=np.int32),
                       np.array([0.1, 0.2, 0.3]), 100000), compare_args=False),
    Case('mapping/out_of_range', D + 'cluster_to_genome_position_mapping',
         args=lambda: (mapping_hic(), [(6500, 99999)], [0.1], None), compare_args=False),
    # create_distance_distribution ----------------------------------------------------
    Case('create_distance_distribution/real', D + 'create_distance_distribution',
         args=lambda: distribution_inputs() + (5, 1, 40), driver='queued'),
    Case('create_distance_distribution/empty_range', D + 'create_distance_distribution',
         args=lambda: distribution_inputs() + (5, 40, 40), driver='queued'),
    Case('create_distance_distribution/mismatch', D + 'create_distance_distribution',
         args=lambda: (distribution_inputs()[0][:-3], distribution_inputs()[1], 5, 1, 4), driver='queued'),
    Case('create_distance_distribution/float_bounds', D + 'create_distance_distribution',
         args=lambda: distribution_inputs() + (5, 1.0, 4), driver='queued'),
    # compute_p_values_mask -----------------------------------------------------------
    Case('compute_p_values_mask/float_threshold', D + 'compute_p_values_mask', args=pvalue_inputs(),
         driver='queued'),
    Case('compute_p_values_mask/tight', D + 'compute_p_values_mask', args=pvalue_inputs(threshold=0.01),
         driver='queued'),
    Case('compute_p_values_mask/some_keys', D + 'compute_p_values_mask', args=pvalue_inputs(keys=[3, 1, 40, 2]),
         driver='queued'),
    Case('compute_p_values_mask/no_keys', D + 'compute_p_values_mask', args=pvalue_inputs(keys=[]),
         driver='queued'),
    Case('compute_p_values_mask/dict_threshold', D + 'compute_p_values_mask',
         args=pvalue_inputs(keys=[1, 2, 3, 4], threshold={5000: 0.5, 10000: 0.01, 15000: 0.9, 20000: 0.2}),
         driver='queued'),
    Case('compute_p_values_mask/dict_missing_key', D + 'compute_p_values_mask',
         args=pvalue_inputs(keys=[1, 2], threshold={5000: 0.5}), driver='queued'),
    Case('compute_p_values_mask/int_threshold', D + 'compute_p_values_mask',
         args=pvalue_inputs(keys=[1, 2], threshold=1), driver='queued'),
    Case('compute_p_values_mask/missing_distance', D + 'compute_p_values_mask',
         args=pvalue_inputs(keys=[1, 100000]), driver='queued'),
    # neighborhood_merge(_thread) ----------------------------------------------------
    Case('neighborhood_merge_thread/real', D + 'neighborhood_merge_thread',
         args=lambda: merge_inputs()()[:3], driver='queued'),
    Case('neighborhood_merge_thread/empty', D + 'neighborhood_merge_thread',
         args=lambda: ([],) + merge_inputs()()[1:3], driver='queued'),
    Case('neighborhood_merge_thread/out_of_range', D + 'neighborhood_merge_thread',
         args=lambda: ([(5, 10000)],) + merge_inputs()()[1:3], driver='queued'),
    Case('neighborhood_merge/one_thread', D + 'neighborhood_merge', args=merge_inputs(1)),
    Case('neighborhood_merge/three_threads', D + 'neighborhood_merge', args=merge_inputs(3)),
    Case('neighborhood_merge/more_threads_than_candidates', D + 'neighborhood_merge', args=merge_inputs(8, top=5)),
    Case('neighborhood_merge/tuples', D + 'neighborhood_merge',
         args=lambda: ([tuple(int(v) for v in c) for c in merge_inputs()()[0]],) + merge_inputs()()[1:]),
    Case('neighborhood_merge/fail', D + 'neighborhood_merge',
         args=lambda: ([(5, 10000), (6, 7)],) + merge_inputs(2)()[1:], driver='stripped'),
    Case('neighborhood_merge/zero_threads', D + 'neighborhood_merge', args=merge_inputs(0)),
    # candidate_region_test(_thread) ---------------------------------------------------
    Case('candidate_region_test_thread/real', D + 'candidate_region_test_thread',
         args=lambda: region_inputs()()[:5], driver='queued'),
    Case('candidate_region_test_thread/empty', D + 'candidate_region_test_thread',
         args=lambda: (region_inputs()()[0], []) + region_inputs()()[2:5], driver='queued'),
    Case('candidate_region_test_thread/peak_wider', D + 'candidate_region_test_thread',
         args=lambda: region_inputs(peak=6)()[:5], driver='queued'),
    Case('candidate_region_test_thread/window_3', D + 'candidate_region_test_thread',
         args=lambda: region_inputs(window=3, peak=1)()[:5], driver='queued'),
    Case('candidate_region_test/one_thread', D + 'candidate_region_test', args=region_inputs(1)),
    Case('candidate_region_test/four_threads', D + 'candidate_region_test', args=region_inputs(4)),
    Case('candidate_region_test/empty', D + 'candidate_region_test',
         args=lambda: (region_inputs()()[0], []) + region_inputs()()[2:]),
    Case('candidate_region_test/out_of_range', D + 'candidate_region_test',
         args=lambda: (region_inputs()()[0], [(5, 10000)]) + region_inputs()()[2:], driver='stripped'),
    # compute_long_range_contacts --------------------------------------------------------
    Case('compute_long_range_contacts/one_thread', D + 'compute_long_range_contacts', args=long_range_inputs(),
         driver='long_range', compare_args=False),
    Case('compute_long_range_contacts/three_threads', D + 'compute_long_range_contacts',
         args=long_range_inputs(3), driver='long_range', compare_args=False),
    Case('compute_long_range_contacts/threshold_dict', D + 'compute_long_range_contacts',
         args=long_range_inputs(2, preselection={k * 5000: 0.4 for k in range(0, 700)}), driver='long_range',
         compare_args=False),
    Case('compute_long_range_contacts/threshold_dict_missing', D + 'compute_long_range_contacts',
         args=long_range_inputs(2, preselection={5000: 0.4}), driver='long_range', compare_args=False),
    Case('compute_long_range_contacts/explicit_zero', D + 'compute_long_range_contacts',
         args=long_range_inputs(explicit_zero=True), driver='long_range', compare_args=False),
    Case('compute_long_range_contacts/float_obs_exp', D + 'compute_long_range_contacts',
         args=long_range_inputs(obs_kind='float'), driver='long_range', compare_args=False),
    # compute_loops ------------------------------------------------------------------------
    Case('compute_loops/h5_chr2L', D + 'compute_loops', args=lambda: (load(), 'chr2L', namespace(), False), **L),
    Case('compute_loops/h5_chr2L_mean_nonzero', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(expected='mean_nonzero'), False), **L),
    Case('compute_loops/h5_chr2L_window_defaults', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(windowSize=None, peakWidth=None), False), **L),
    Case('compute_loops/h5_chr2L_tight', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(pValuePreselection='0.01', maxLoopDistance=100000), False), **L),
    Case('compute_loops/h5_chr2L_threshold_file', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(pValuePreselection=THRESHOLD_FILE, maxLoopDistance=200000), False),
         **L),
    Case('compute_loops/h5_chr2L_threshold_file_fail', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(pValuePreselection=THRESHOLD_FILE), False), **L),
    Case('compute_loops/h5_chrM_too_small', D + 'compute_loops', args=lambda: (load(), 'chrM', namespace(), False), **L),
    Case('compute_loops/h5_unknown_expected', D + 'compute_loops',
         args=lambda: (load(), 'chr2L', namespace(expected='median'), False), **L),
    Case('compute_loops/h5_unknown_chromosome', D + 'compute_loops',
         args=lambda: (load(), 'chrZ', namespace(), False), **L),
    Case('compute_loops/li_x', D + 'compute_loops',
         args=lambda: (load(LI_H5), 'X', namespace(matrix=LI_H5, maxLoopDistance=100000), False), **L),
    Case('compute_loops/cool_queue', D + 'compute_loops',
         args=lambda: ('1', namespace(matrix=MATRIX_COOL, maxLoopDistance=50000000, pValuePreselection=0.55), True),
         driver='loops_queued', compare_args=False),
    Case('compute_loops/cool_queue_diagonal_only', D + 'compute_loops',
         args=lambda: ('2', namespace(matrix=MATRIX_COOL, maxLoopDistance=3000000, pValuePreselection=0.55), True),
         driver='loops_queued', compare_args=False),
    Case('compute_loops/h5_queue', D + 'compute_loops',
         args=lambda: ('chr2L', namespace(maxLoopDistance=1000000), False), driver='loops_queued',
         compare_args=False),
]


# Cases whose floats differ from 3.7.6 in the last bits since the special
# functions (normal tail, incomplete beta and gamma functions, log-gamma) are
# Boost.Math's and no longer the routines scipy uses: the p-values agree
# to a relative 2e-14 or better. The cases compare every float at rtol 1e-12;
# integers, strings, coordinates and the structure of the result are still
# compared exactly, so a changed call, a moved TAD or a different row count
# fails.
_BOOST_LAST_BITS = {
    'candidate_region_test/four_threads',
    'candidate_region_test/one_thread',
    'candidate_region_test_thread/real',
    'candidate_region_test_thread/window_3',
    'compute_long_range_contacts/one_thread',
    'compute_long_range_contacts/three_threads',
    'compute_long_range_contacts/threshold_dict',
    'compute_loops/cool_queue',
    'compute_loops/h5_chr2L',
    'compute_loops/h5_chr2L_mean_nonzero',
    'compute_loops/h5_chr2L_threshold_file',
    'compute_loops/h5_chr2L_tight',
    'compute_loops/h5_chr2L_window_defaults',
    'compute_loops/h5_queue',
    'compute_loops/li_x',
}
for _case in CASES:
    if _case.id in _BOOST_LAST_BITS:
        _case.rtol = 1e-12
        _case.note = (_case.note + '; ' if _case.note else '') + (
            'Boost.Math special functions round the last bits of the p-values differently from scipy 1.14')
