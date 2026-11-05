"""Equivalence cases for hicexplorer.chicSignificantInteractions (apieq.py).

The functions run as chicSignificantInteractions' main() runs them: on the
interaction file of the cHi-C test data (chicViewpoint/two_matrices.hdf5)
with the repository's background model and threshold files, in dual and
single combination mode, with an x-fold or a loose p-value preselection or
none. The output files are compared by an h5py dump; errorLog.txt, which
compute_interaction_file appends to in the working directory, is read back
from the case directory.

call_multi_core runs compute_interaction_file in this process
(cpp/python/src/api/chicSignificantInteractions.cpp); its results are those
of 3.7.6's processes. The records the workers log reach the caller's
handlers now, which 3.7.6's worker processes did not deliver, so the
call_multi_core cases use inputs on which the workers log nothing. A failure's traceback text is not compared (the frames
are not 3.7.6's): the drivers keep the text before 'Traceback'.
"""
import argparse
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicSignificantInteractions.'

INTERACTIONS = os.path.join(CHIC, 'chicViewpoint', 'two_matrices.hdf5')
BACKGROUND = os.path.join(CHIC, 'background.txt')
XFOLD = os.path.join(CHIC, 'thresholdFile_xfold.txt')
PVALUE = os.path.join(CHIC, 'thresholdFile_pValue.txt')
LOOSE = os.path.join(CHIC, 'thresholdFile_loose_pValue.txt')
FL, MB = 'FL-E13-5_chr1', 'MB-E10-5_chr1'
GENES = ('Eya1', 'Sox17', 'Tfap2d')
DUAL = [[[FL, 'chr1', gene], [MB, 'chr1', gene]] for gene in GENES]
SINGLE = [[[sample, 'chr1', gene]] for sample in (FL, MB) for gene in GENES]
RESOLUTION = np.int64(1000)


def dump_h5(path):
    """Every group, attribute and dataset of an HDF5 file with h5py."""
    import h5py
    out = []

    def attributes(obj):
        return sorted((k, repr(v), type(v).__name__) for k, v in obj.attrs.items())

    with h5py.File(path, 'r') as f:
        out.append(('/', attributes(f)))

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                value = obj[()]
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                out.append((name, 'dataset', obj.dtype.str, obj.shape, obj.compression, obj.compression_opts,
                            obj.chunks, value, attributes(obj)))
            else:
                out.append((name, 'group', attributes(obj)))
        f.visititems(visit)
    return out


def _threshold(value):
    if isinstance(value, str):
        from hicexplorer import chicSignificantInteractions
        return chicSignificantInteractions.read_threshold_file(value)
    return value


def _args(pvalue=0.2, xfold=1.5, loose=None, truncate=False, peak=5, threads=1, mode='dual'):
    return argparse.Namespace(pValue=_threshold(pvalue), xFoldBackground=_threshold(xfold),
                              loosePValue=_threshold(loose), truncateZeroPvalues=truncate,
                              peakInteractionsThreshold=peak, threads=threads, range=[200000, 200000],
                              combinationMode=mode, fixateRange=500000)


def _background():
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint().readBackgroundDataFile(BACKGROUND, [200000, 200000], 500000)


def _data(sample):
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint().readInteractionFile(INTERACTIONS, sample)


def _short(result):
    if isinstance(result, str):
        return result.split('Traceback')[0]
    return result


class _InDirectory(object):
    def __init__(self, directory):
        self.directory = directory

    def __enter__(self):
        self.previous = os.getcwd()
        os.chdir(self.directory)

    def __exit__(self, *exc):
        os.chdir(self.previous)


def _error_log(directory):
    path = os.path.join(directory, 'errorLog.txt')
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return handle.read()


def threshold(target, path):
    return target(path)


def filter_data(target, sample, pvalue, peak):
    data = _data(sample)
    return target(data[0], _threshold(pvalue), data[1], peak), data[0]


def merge_x(target, sample, xfold):
    from hicexplorer.lib.viewpoint import Viewpoint
    return target(_threshold(xfold), _data(sample), Viewpoint(), RESOLUTION)


def merge_loose(target, sample, loose, truncate):
    from hicexplorer.lib.viewpoint import Viewpoint
    return target(_threshold(loose), _data(sample), Viewpoint(), RESOLUTION, truncate)


def new_p_values(target, sample, xfold, pvalue, peak, drop_background_key=None):
    from hicexplorer import chicSignificantInteractions
    from hicexplorer.lib.viewpoint import Viewpoint
    scores, merged = chicSignificantInteractions.merge_neighbors_x_fold(_threshold(xfold), _data(sample),
                                                                        Viewpoint(), RESOLUTION)
    background = _background()
    if drop_background_key is not None:
        background.pop(drop_background_key, None)
    result = target(scores, background, _threshold(pvalue), merged, peak)
    return result, scores


def interaction_file(target, directory, files, use_queue=True, path=INTERACTIONS, **settings):
    import queue
    from hicexplorer.lib.viewpoint import Viewpoint
    q = queue.Queue() if use_queue else None
    with _InDirectory(directory):
        result = target(files, _args(**settings), Viewpoint(), _background(), path, RESOLUTION, q)
    out = _short(q.get()) if use_queue else None
    return _short(result), out, _error_log(directory)


class _WithoutTraceback(object):
    """Cuts a logged failure message before its traceback text."""

    def filter(self, record):
        if isinstance(record.msg, str) and 'Traceback' in record.msg:
            record.msg = record.msg.split('Traceback')[0]
        return True


def multi_core(target, directory, files, **settings):
    import logging
    from hicexplorer.lib.viewpoint import Viewpoint
    logger = logging.getLogger('hicexplorer.chicSignificantInteractions')
    cut = _WithoutTraceback()
    logger.addFilter(cut)
    try:
        with _InDirectory(directory):
            try:
                result = target(files, _args(**settings), Viewpoint(), _background(), INTERACTIONS, RESOLUTION)
            except SystemExit as error:
                result = ('SystemExit', error.code)
    finally:
        logger.removeFilter(cut)
    return result, _error_log(directory)


def write_files(target, directory, out, files, **settings):
    from hicexplorer import chicSignificantInteractions as s
    from hicexplorer.lib.viewpoint import Viewpoint
    args = _args(**settings)
    with _InDirectory(directory):
        significant, keys, target_data, target_keys, points_target, points_significant = s.call_multi_core(
            files, args, Viewpoint(), _background(), INTERACTIONS, RESOLUTION)
    if target.__name__ == 'writeSignificantHDF':
        try:
            target(out, significant, keys, Viewpoint(), points_significant, args)
        except Exception as error:
            return 'raised', type(error).__name__, str(error), dump_h5(out) if os.path.exists(out) else None
    else:
        try:
            target(out, target_data, target_keys, Viewpoint(), RESOLUTION, points_target, args)
        except Exception as error:
            return 'raised', type(error).__name__, str(error), dump_h5(out) if os.path.exists(out) else None
    return dump_h5(out)


SCORES = {np.int64(-1000): ['chr1', 100, 200, 'g', 5.0, -1000, 0.1, 0.5, 1.5, 7.0],
          np.int64(0): ['chr1', 200, 300, 'g', 5.0, 0, 0.2, 0.01, 2.5, 3.0],
          np.int64(1000): ['chr1', 300, 400, 'g', 5.0, 1000, 0.3, 0.3, 3.5, 9.0]}

CASES = [
    # read_threshold_file ----------------------------------------------------------
    Case('read_threshold_file/xfold', M + 'read_threshold_file', driver='threshold', args=(XFOLD,)),
    Case('read_threshold_file/pvalue', M + 'read_threshold_file', driver='threshold', args=(PVALUE,)),
    Case('read_threshold_file/blank_line_stops', M + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b"# c\n-1000\t0.5\n\n1000\t0.25\n"}),
    Case('read_threshold_file/bad_line', M + 'read_threshold_file', args=(File('t.txt'),),
         files={'t.txt': b"-1000\t0.5\t3\n"}),
    Case('read_threshold_file/missing', M + 'read_threshold_file', args=(File('missing.txt'),)),

    # filter_by_pvalue -----------------------------------------------------------
    Case('filter_by_pvalue/float', M + 'filter_by_pvalue', driver='filter_data', args=([FL, 'chr1', 'Eya1'], 0.2, 1.0)),
    Case('filter_by_pvalue/threshold_file', M + 'filter_by_pvalue', driver='filter_data',
         args=([MB, 'chr1', 'Sox17'], PVALUE, 1.0)),
    Case('filter_by_pvalue/int_threshold', M + 'filter_by_pvalue', driver='filter_data',
         args=([MB, 'chr1', 'Sox17'], 1, 1.0)),
    Case('filter_by_pvalue/constructed', M + 'filter_by_pvalue',
         args=lambda: ({k: [v[6], v[7], v[9], v[8]] for k, v in SCORES.items()}, 0.4, SCORES, 2.0)),

    # merge_neighbors_x_fold / merge_neighbors_loose_p_value ---------------------
    Case('merge_neighbors_x_fold/float', M + 'merge_neighbors_x_fold', driver='merge_x', args=([FL, 'chr1', 'Eya1'], 1.5)),
    Case('merge_neighbors_x_fold/threshold_file', M + 'merge_neighbors_x_fold', driver='merge_x',
         args=([MB, 'chr1', 'Tfap2d'], XFOLD)),
    Case('merge_neighbors_x_fold/nothing_accepted', M + 'merge_neighbors_x_fold', driver='merge_x',
         args=([MB, 'chr1', 'Tfap2d'], 1e9)),
    Case('merge_neighbors_x_fold/int', M + 'merge_neighbors_x_fold', driver='merge_x', args=([MB, 'chr1', 'Tfap2d'], 2)),
    Case('merge_neighbors_loose_p_value/float', M + 'merge_neighbors_loose_p_value', driver='merge_loose',
         args=([FL, 'chr1', 'Sox17'], 0.3, False)),
    Case('merge_neighbors_loose_p_value/truncate', M + 'merge_neighbors_loose_p_value', driver='merge_loose',
         args=([FL, 'chr1', 'Sox17'], 0.3, True)),
    Case('merge_neighbors_loose_p_value/threshold_file', M + 'merge_neighbors_loose_p_value', driver='merge_loose',
         args=([MB, 'chr1', 'Eya1'], LOOSE, True)),

    # compute_new_p_values -----------------------------------------------------------
    Case('compute_new_p_values/float', M + 'compute_new_p_values', driver='new_p_values',
         args=([FL, 'chr1', 'Eya1'], 1.5, 0.2, 5)),
    Case('compute_new_p_values/threshold_file', M + 'compute_new_p_values', driver='new_p_values',
         args=([MB, 'chr1', 'Tfap2d'], XFOLD, PVALUE, 1)),
    Case('compute_new_p_values/missing_background_key', M + 'compute_new_p_values', driver='new_p_values',
         args=([FL, 'chr1', 'Eya1'], 1.5, 0.2, 5, -1000)),
    Case('compute_new_p_values/int_threshold', M + 'compute_new_p_values', driver='new_p_values',
         args=([FL, 'chr1', 'Eya1'], 1.5, 1, 5)),

    # compute_interaction_file -----------------------------------------------------
    Case('compute_interaction_file/dual_xfold', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL)),
    Case('compute_interaction_file/single_loose', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), SINGLE), kwargs={'xfold': None, 'loose': 0.3, 'truncate': True, 'mode': 'single'}),
    Case('compute_interaction_file/no_preselection', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[:1]), kwargs={'xfold': None, 'pvalue': 0.05, 'peak': 1.2}),
    Case('compute_interaction_file/threshold_files', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[1:]), kwargs={'xfold': XFOLD, 'pvalue': PVALUE}),
    Case('compute_interaction_file/nothing_accepted', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[:2]), kwargs={'xfold': 1e9}),
    Case('compute_interaction_file/no_queue', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[:1], False)),
    Case('compute_interaction_file/failure', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[:1], True, os.path.join(CHIC, 'missing.hdf5'))),
    Case('compute_interaction_file/failure_no_queue', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), DUAL[:1], False, os.path.join(CHIC, 'missing.hdf5'))),
    Case('compute_interaction_file/empty_combination', M + 'compute_interaction_file', driver='interaction_file',
         args=(File(), [[]])),

    # call_multi_core ----------------------------------------------------------------
    Case('call_multi_core/dual_one_thread', M + 'call_multi_core', driver='multi_core', args=(File(), DUAL)),
    Case('call_multi_core/dual_two_threads', M + 'call_multi_core', driver='multi_core', args=(File(), DUAL),
         kwargs={'threads': 2}),
    Case('call_multi_core/single_three_threads', M + 'call_multi_core', driver='multi_core', args=(File(), SINGLE),
         kwargs={'threads': 3, 'xfold': None, 'loose': 0.2, 'mode': 'single'}),
    Case('call_multi_core/repeated_samples', M + 'call_multi_core', driver='multi_core',
         args=(File(), DUAL + DUAL[:1]), kwargs={'threads': 2}),
    Case('call_multi_core/zero_threads', M + 'call_multi_core', driver='multi_core', args=(File(), DUAL),
         kwargs={'threads': 0}),
    Case('call_multi_core/negative_threads', M + 'call_multi_core', driver='multi_core', args=(File(), DUAL),
         kwargs={'threads': -1}),
    Case('call_multi_core/failure', M + 'call_multi_core', driver='multi_core',
         args=(File(), DUAL[:1] + [[]])),

    # writeSignificantHDF / writeTargetHDF -------------------------------------------
    Case('writeSignificantHDF/dual_xfold', M + 'writeSignificantHDF', driver='write_files',
         args=(File(), File('s.hdf5'), DUAL)),
    Case('writeSignificantHDF/single_loose', M + 'writeSignificantHDF', driver='write_files',
         args=(File(), File('s.hdf5'), SINGLE), kwargs={'xfold': None, 'loose': 0.3, 'mode': 'single'}),
    Case('writeSignificantHDF/no_preselection', M + 'writeSignificantHDF', driver='write_files',
         args=(File(), File('s.hdf5'), DUAL), kwargs={'xfold': None, 'pvalue': 0.05, 'peak': 1.2, 'truncate': True}),
    Case('writeSignificantHDF/threshold_dict', M + 'writeSignificantHDF', driver='write_files',
         args=(File(), File('s.hdf5'), DUAL[:1]), kwargs={'xfold': XFOLD}),
    Case('writeSignificantHDF/repeated_samples', M + 'writeSignificantHDF', driver='write_files',
         args=(File(), File('s.hdf5'), DUAL + DUAL[:1]), kwargs={'threads': 2}),
    Case('writeTargetHDF/dual_xfold', M + 'writeTargetHDF', driver='write_files',
         args=(File(), File('t.hdf5'), DUAL)),
    Case('writeTargetHDF/single_loose', M + 'writeTargetHDF', driver='write_files',
         args=(File(), File('t.hdf5'), SINGLE), kwargs={'xfold': None, 'loose': 0.3, 'mode': 'single'}),
    Case('writeTargetHDF/no_preselection', M + 'writeTargetHDF', driver='write_files',
         args=(File(), File('t.hdf5'), DUAL), kwargs={'xfold': None, 'pvalue': 0.05, 'peak': 1.2}),
    Case('writeTargetHDF/repeated_genes', M + 'writeTargetHDF', driver='write_files',
         args=(File(), File('t.hdf5'), DUAL + DUAL[:2]), kwargs={'threads': 2}),
    Case('writeTargetHDF/missing_directory', M + 'writeTargetHDF', driver='write_files',
         args=(File(), File('no/t.hdf5'), DUAL)),
]
