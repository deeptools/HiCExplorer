"""Equivalence cases for hicexplorer.chicAggregateStatistic (apieq.py).

The scores are the interaction file of the cHi-C test data
(chicViewpoint/two_matrices.hdf5), the targets the three kinds main() knows:
an HDF5 target file of chicSignificantInteractions, a bed4 target dict
(utilities.readTargetBed) and a bed3 file (as a path for
run_target_list_compilation, as an interval tree for
filter_scores_target_list). The aggregate files are compared by an h5py dump.

call_multi_core runs the worker in this process (cpp/python/src/api/
chicAggregateStatistic.cpp): its results are those of 3.7.6's processes, and
the cases cover the splitting into threads, one and several targets and a
failing worker (exit(1) after logging the message).
"""
import argparse
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicAggregateStatistic.'

INTERACTIONS = os.path.join(CHIC, 'chicViewpoint', 'two_matrices.hdf5')
TARGET_DUAL = os.path.join(CHIC, 'chicSignificantInteractions', 'targetFile_dual.hdf5')
BED3 = os.path.join(CHIC, 'chicAggregateStatistic', 'target_list_3col.bed')
BED4 = os.path.join(CHIC, 'chicAggregateStatistic', 'target_list_4col.bed')

FL, MB = 'FL-E13-5_chr1', 'MB-E10-5_chr1'
PAIRS = [[[FL, 'chr1', gene], [MB, 'chr1', gene]] for gene in ('Eya1', 'Sox17', 'Tfap2d')]
TARGETS = [[FL, MB, 'genes', gene] for gene in ('Eya1', 'Sox17', 'Tfap2d')]


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


def scores(sample):
    from hicexplorer.lib.viewpoint import Viewpoint
    return Viewpoint().readInteractionFile(INTERACTIONS, sample)[1]


def filter_hdf5(target, sample, target_list, target_file=TARGET_DUAL):
    data = scores(sample)
    result = target(data, 'hdf5', None, pTargetList=target_list, pTargetFile=target_file)
    return result, data


def filter_bed4(target, sample, gene):
    from hicexplorer import utilities
    _, _, positions = utilities.readTargetBed(BED4)
    data = scores(sample)
    return target(data, 'bed4', positions, pTargetList=['c_adj_norm', 't_adj_norm', 'genes', gene]), positions


def filter_bed3(target, sample, path=BED3):
    import hicmatrix.HiCMatrix as hm
    from hicexplorer import utilities
    tree = hm.hiCMatrix().intervalListToIntervalTree(utilities.readBed(path))[0]
    return target(scores(sample), 'bed3', None, pTargetIntervalTree=tree)


def filter_constructed(target, data, *args, **kwargs):
    result = target(data, *args, **kwargs)
    return result, data


def run_compilation(target, files, targets, ftype, one_target, use_queue=True, positions=False):
    import queue
    from hicexplorer import utilities
    from hicexplorer.lib.viewpoint import Viewpoint
    pos = utilities.readTargetBed(BED4)[2] if positions else None
    args = argparse.Namespace(interactionFile=INTERACTIONS, targetFile=TARGET_DUAL)
    q = queue.Queue() if use_queue else None
    result = target(files, targets, ftype, pos, args, Viewpoint(), q, one_target)
    if q is None:
        return result
    out = q.get()
    if isinstance(out, str):
        return out.split('Traceback')[0]
    return result, out


def multi_core(target, files, targets, ftype, threads, positions=False):
    from hicexplorer import chicAggregateStatistic, utilities
    from hicexplorer.lib.viewpoint import Viewpoint
    pos = utilities.readTargetBed(BED4)[2] if positions else None
    args = argparse.Namespace(interactionFile=INTERACTIONS, targetFile=TARGET_DUAL, threads=threads)
    try:
        result = target(files, targets, ftype, pos, chicAggregateStatistic.run_target_list_compilation, args,
                        Viewpoint())
    except SystemExit as error:
        return 'SystemExit', error.code, args.threads
    return result, args.threads


def failing_worker(pInteractionFilesList, pTargetList, pTargetFType, pTargetPosDict, pArgs, pViewpointObj, pQueue,
                   pOneTarget):
    """A worker that fails on its second slice (the traceback text 3.7.6's
    workers append is not part of the comparison)."""
    if pInteractionFilesList and pInteractionFilesList[0][0][2] != 'Eya1':
        pQueue.put('Fail: the worker of {} failed'.format(len(pInteractionFilesList)))
    else:
        pQueue.put([[['name', len(pInteractionFilesList), pOneTarget, pTargetList]], [{'x': pTargetFType}]])


def multi_core_worker(target, files, targets, threads):
    args = argparse.Namespace(threads=threads)
    try:
        result = target(files, targets, 'bed3', None, failing_worker, args, None)
    except SystemExit as error:
        return 'SystemExit', error.code, args.threads
    return result, args.threads


def write_aggregate(target, out, files, targets, ftype, positions=False):
    from hicexplorer import chicAggregateStatistic, utilities
    from hicexplorer.lib.viewpoint import Viewpoint
    pos = utilities.readTargetBed(BED4)[2] if positions else None
    args = argparse.Namespace(interactionFile=INTERACTIONS, targetFile=TARGET_DUAL, threads=2)
    names, accepted = chicAggregateStatistic.call_multi_core(
        files, targets, ftype, pos, chicAggregateStatistic.run_target_list_compilation, args, Viewpoint())
    target(out, names, accepted, args)
    return dump_h5(out)


def write_constructed(target, out, names, accepted):
    target(out, names, accepted, None)
    return dump_h5(out)


def lines():
    """Score lines with targets that two, one and no position reach."""
    return ({np.int64(-2000): ['chr1', 100, 200, 'g', 5.0, -2000, 0.1, 0.5, 1.5, 3.0],
             np.int64(-1000): ['chr1', 200, 300, 'g', 5.0, -1000, 0.2, 0.25, 2.5, 4.0],
             np.int64(0): ['chr1', 300, 400, 'g', 5.0, 0, 0.3, 0.125, 3.5, 5.0],
             np.int64(1000): ['chr2', 300, 400, 'g', 5.0, 1000, 0.4, 0.0625, 4.5, 6.0],
             np.int64(3000): ['chr1', 900, 950, 'g', 5.0, 3000, 0.5, 0.03125, 5.5, 7.0]},)


POSITIONS = {'g': {'chromosome': 'chr1', 'start_list': ['150', '350', '980'], 'end_list': ['250', '360', '990']},
             'h': {'chromosome': 'chr1', 'start_list': [], 'end_list': []},
             'n': {'chromosome': 'chr1', 'start_list': ['5'], 'end_list': ['5']}}

AGGREGATED_NAMES = [[[FL, 'chr1', 'g'], [MB, 'chr1', 'g']]]

CASES = [
    # filter_scores_target_list -----------------------------------------------
    Case('filter_scores_target_list/hdf5_eya1', M + 'filter_scores_target_list', driver='filter_hdf5',
         args=([FL, 'chr1', 'Eya1'], TARGETS[0])),
    Case('filter_scores_target_list/hdf5_tfap2d_mb', M + 'filter_scores_target_list', driver='filter_hdf5',
         args=([MB, 'chr1', 'Tfap2d'], TARGETS[2])),
    Case('filter_scores_target_list/hdf5_missing_target', M + 'filter_scores_target_list', driver='filter_hdf5',
         args=([MB, 'chr1', 'Tfap2d'], [FL, MB, 'genes', 'Nope'])),
    Case('filter_scores_target_list/bed4', M + 'filter_scores_target_list', driver='filter_bed4',
         args=([FL, 'chr1', 'Sox17'], 'Sox17')),
    Case('filter_scores_target_list/bed4_other_gene', M + 'filter_scores_target_list', driver='filter_bed4',
         args=([FL, 'chr1', 'Sox17'], 'Eya1')),
    Case('filter_scores_target_list/bed3_tree', M + 'filter_scores_target_list', driver='filter_bed3',
         args=([FL, 'chr1', 'Eya1'],)),
    Case('filter_scores_target_list/bed3_tree_mb', M + 'filter_scores_target_list', driver='filter_bed3',
         args=([MB, 'chr1', 'Tfap2d'],)),
    Case('filter_scores_target_list/constructed_bed4', M + 'filter_scores_target_list', driver='filter_constructed',
         args=lambda: lines() + ('bed4', POSITIONS), kwargs={'pTargetList': ['x', 'g']}),
    Case('filter_scores_target_list/constructed_empty_targets', M + 'filter_scores_target_list',
         driver='filter_constructed', args=lambda: lines() + ('bed4', POSITIONS), kwargs={'pTargetList': ['h']}),
    Case('filter_scores_target_list/constructed_null_interval', M + 'filter_scores_target_list',
         driver='filter_constructed', args=lambda: lines() + ('bed4', POSITIONS), kwargs={'pTargetList': ['n']}),
    Case('filter_scores_target_list/bed3_without_tree', M + 'filter_scores_target_list', driver='filter_constructed',
         args=lambda: lines() + ('bed3', None)),
    Case('filter_scores_target_list/bed3_with_list', M + 'filter_scores_target_list', driver='filter_constructed',
         args=lambda: lines() + ('bed3', None), kwargs={'pTargetList': ['g']}),
    Case('filter_scores_target_list/unknown_type', M + 'filter_scores_target_list', driver='filter_constructed',
         args=lambda: lines() + ('bed5', None)),

    # run_target_list_compilation -----------------------------------------------
    Case('run_target_list_compilation/hdf5', M + 'run_target_list_compilation', driver='run_compilation',
         args=(PAIRS, TARGETS, 'hdf5', False)),
    Case('run_target_list_compilation/bed4', M + 'run_target_list_compilation', driver='run_compilation',
         args=(PAIRS[:2], [['c_adj_norm', 't_adj_norm', 'genes', 'Eya1'], ['c_adj_norm', 't_adj_norm', 'genes', 'Sox17']],
               'bed4', False, True, True)),
    Case('run_target_list_compilation/bed3', M + 'run_target_list_compilation', driver='run_compilation',
         args=(PAIRS, BED3, 'bed3', True)),
    Case('run_target_list_compilation/bed3_missing', M + 'run_target_list_compilation', driver='run_compilation',
         args=(PAIRS, os.path.join(CHIC, 'missing.bed'), 'bed3', True)),
    Case('run_target_list_compilation/bed3_not_one_target', M + 'run_target_list_compilation',
         driver='run_compilation', args=(PAIRS[:1], BED3, 'bed3', False)),
    Case('run_target_list_compilation/no_queue', M + 'run_target_list_compilation', driver='run_compilation',
         args=(PAIRS[:1], TARGETS[:1], 'hdf5', False, False)),
    Case('run_target_list_compilation/no_files', M + 'run_target_list_compilation', driver='run_compilation',
         args=([], TARGETS, 'hdf5', False)),

    # call_multi_core -------------------------------------------------------------
    Case('call_multi_core/hdf5_two_threads', M + 'call_multi_core', driver='multi_core',
         args=(PAIRS, TARGETS, 'hdf5', 2)),
    Case('call_multi_core/hdf5_more_threads_than_files', M + 'call_multi_core', driver='multi_core',
         args=(PAIRS, TARGETS, 'hdf5', 8)),
    Case('call_multi_core/bed3_one_target', M + 'call_multi_core', driver='multi_core',
         args=(PAIRS, [BED3], 'bed3', 2)),
    Case('call_multi_core/failing_worker', M + 'call_multi_core', driver='multi_core_worker',
         args=(PAIRS, ['t1', 't2', 't3'], 2)),
    Case('call_multi_core/worker_one_thread', M + 'call_multi_core', driver='multi_core_worker',
         args=(PAIRS, ['t1'], 1)),
    Case('call_multi_core/no_files', M + 'call_multi_core', driver='multi_core',
         args=([], TARGETS, 'hdf5', 2)),

    # writeAggregateHDF -----------------------------------------------------------
    Case('writeAggregateHDF/hdf5', M + 'writeAggregateHDF', driver='write_aggregate',
         args=(File('a.hdf5'), PAIRS, TARGETS, 'hdf5')),
    Case('writeAggregateHDF/bed3', M + 'writeAggregateHDF', driver='write_aggregate',
         args=(File('a.hdf5'), PAIRS[1:], [BED3], 'bed3')),
    Case('writeAggregateHDF/constructed_genes_key', M + 'writeAggregateHDF', driver='write_constructed',
         args=lambda: (File('a.hdf5'), [[[FL, 'genes', 'g'], [MB, 'chr1', 'g']]],
                       [[{0: ['chr1', 1, 2, 'g', 3.0, 0, 0.5, 0.1, 1.5, 4.0]},
                         {0: ['chr1', 1, 2, 'g', 3.0, 0, 0.5, 0.1, 1.5, np.float64(4.0)]}]])),
    Case('writeAggregateHDF/constructed_empty', M + 'writeAggregateHDF', driver='write_constructed',
         args=lambda: (File('a.hdf5'), AGGREGATED_NAMES, [[{}, {}]])),
    Case('writeAggregateHDF/constructed_twice', M + 'writeAggregateHDF', driver='write_constructed',
         args=lambda: (File('a.hdf5'), AGGREGATED_NAMES * 2,
                       [[{0: ['chr1', 1, 2, 'g', 3.0, 0, 4.0]}, {}]] * 2)),
]
