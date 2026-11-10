"""Equivalence cases for hicexplorer.chicExportData (apieq.py).

exportData runs as chicExportData's main() calls it, on the interaction,
significant, target, aggregate and differential files of the cHi-C test
data, for text and bigwig output; the driver returns what it put into the
queue. A failure puts 'Fail: ' + str(exp) + traceback.format_exc() into the
queue: the traceback names the frames the error passed, which are not
3.7.6's, so the driver compares the text before 'Traceback'.
"""
import argparse
import os

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicExportData.'

INTERACTIONS = os.path.join(CHIC, 'chicViewpoint', 'two_matrices.hdf5')
SIGNIFICANT = os.path.join(CHIC, 'chicSignificantInteractions', 'significantInteractions_dual.hdf5')
TARGET = os.path.join(CHIC, 'chicSignificantInteractions', 'targetFile_dual.hdf5')
AGGREGATE = os.path.join(CHIC, 'chicAggregateStatistic', 'aggregate.hdf5')
DIFFERENTIAL = os.path.join(CHIC, 'chicDifferentialTest', 'differential.hdf5')
BACKGROUND = os.path.join(CHIC, 'background.txt')


def export(target, file_list, path, file_type, output='txt', decimals=12, value='p-value', background=False,
           sizes=None):
    import queue
    from hicexplorer.lib.viewpoint import Viewpoint
    viewpoint = Viewpoint()
    background_data = None
    if background:
        background_data = viewpoint.readBackgroundDataFile(BACKGROUND, [200000, 200000], 200000, pMean=True)
    args = argparse.Namespace(file=path, outputFileType=output, outputValueBigwig=value,
                              backgroundModelFile=BACKGROUND if background else None)
    q = queue.Queue()
    target(file_list, args, viewpoint, decimals, sizes, background_data, file_type, q)
    result = q.get()
    if isinstance(result, str):
        return result.split('Traceback')[0]
    return result


SAMPLES = [[['FL-E13-5_chr1', 'chr1', 'Eya1']], [['MB-E10-5_chr1', 'genes', 'Sox17']]]
SIZES = {'chr1': 197195432}

CASES = [
    Case('exportData/interactions_txt', M + 'exportData', driver='export',
         args=(SAMPLES, INTERACTIONS, 'interactions')),
    Case('exportData/interactions_txt_decimals', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1', 'chr1', 'Tfap2d'], ['MB-E10-5_chr1', 'chr1', 'Tfap2d']]], INTERACTIONS,
               'interactions', 'txt', 3)),
    Case('exportData/significant_txt', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1', 'chr1', 'Sox17']], [['MB-E10-5_chr1', 'chr1', 'Tfap2d']]], SIGNIFICANT,
               'significant', 'txt', 5)),
    Case('exportData/interactions_missing_sample', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1', 'chr1', 'Nope']]], INTERACTIONS, 'interactions')),
    Case('exportData/interactions_bigwig', M + 'exportData', driver='export',
         args=(SAMPLES, INTERACTIONS, 'interactions', 'bigwig', 12, 'raw', False, SIZES)),
    Case('exportData/interactions_bigwig_background', M + 'exportData', driver='export',
         args=(SAMPLES, INTERACTIONS, 'interactions', 'bigwig', 12, 'relative-interactions', True, SIZES)),
    Case('exportData/interactions_bigwig_xfold', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1', 'chr1', 'Eya1'], ['MB-E10-5_chr1', 'chr1', 'Eya1']]], INTERACTIONS,
               'interactions', 'bigwig', 12, 'x-fold', True, SIZES)),
    Case('exportData/interactions_bigwig_no_size', M + 'exportData', driver='export',
         args=(SAMPLES, INTERACTIONS, 'interactions', 'bigwig', 12, 'p-value', False, {})),
    Case('exportData/interactions_bigwig_empty_sample', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1', 'chr1', 'Nope']]], INTERACTIONS, 'interactions', 'bigwig', 12, 'raw', False,
               SIZES)),
    Case('exportData/interactions_other_format', M + 'exportData', driver='export',
         args=(SAMPLES, INTERACTIONS, 'interactions', 'png')),
    Case('exportData/no_samples', M + 'exportData', driver='export', args=([[]], INTERACTIONS, 'interactions')),
    Case('exportData/target', M + 'exportData', driver='export',
         args=([['FL-E13-5_chr1', 'MB-E10-5_chr1', 'genes', 'Eya1'], ['FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1',
                                                                      'Tfap2d']], TARGET, 'target')),
    Case('exportData/target_not_strings', M + 'exportData', driver='export',
         args=([['FL-E13-5_chr1', 'chr1', 'Eya1']], INTERACTIONS, 'target')),
    Case('exportData/target_missing', M + 'exportData', driver='export',
         args=([['FL-E13-5_chr1', 'x']], TARGET, 'target')),
    Case('exportData/aggregate', M + 'exportData', driver='export',
         args=([[['FL-E13-5_chr1_MB-E10-5_chr1', 'FL-E13-5_chr1', 'chr1', 'Eya1'],
                 ['FL-E13-5_chr1_MB-E10-5_chr1', 'MB-E10-5_chr1', 'genes', 'Eya1']]], AGGREGATE, 'aggregate', 'txt', 4)),
    Case('exportData/differential', M + 'exportData', driver='export',
         args=([['FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1', 'Eya1'], ['FL-E13-5_chr1', 'MB-E10-5_chr1', 'genes', 'Sox17'],
                ['FL-E13-5_chr1', 'MB-E10-5_chr1', 'chr1', 'Nope']], DIFFERENTIAL, 'differential', 'txt', 5)),
    Case('exportData/unknown_type', M + 'exportData', driver='export', args=(SAMPLES, INTERACTIONS, 'other')),
    Case('exportData/missing_file', M + 'exportData', driver='export',
         args=([['a']], os.path.join(CHIC, 'missing.hdf5'), 'target')),
]
