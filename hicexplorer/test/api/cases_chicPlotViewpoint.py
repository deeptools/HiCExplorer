"""Equivalence cases for hicexplorer.chicPlotViewpoint (apieq.py).

plot_images draws with matplotlib (it is on the drawing allowlist); the data
it draws come from the Viewpoint object and from hicexplorer._hicx. The
cases run it as chicPlotViewpoint's main() does on the interaction,
significant and differential files of the cHi-C test data and compare the
PNG images it puts into the queue byte for byte (both interpreters draw with
matplotlib 3.8.4), together with the file names.
"""
import argparse
import os

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
CHIC = os.path.join(HERE, '..', 'test_data', 'cHi-C')
M = 'hicexplorer.chicPlotViewpoint.'

INTERACTIONS = os.path.join(CHIC, 'chicViewpoint', 'two_matrices.hdf5')
SIGNIFICANT = os.path.join(CHIC, 'chicSignificantInteractions', 'significantInteractions_dual.hdf5')
DIFFERENTIAL = os.path.join(CHIC, 'chicDifferentialTest', 'differential.hdf5')
BACKGROUND = os.path.join(CHIC, 'background.txt')
FL, MB = 'FL-E13-5_chr1', 'MB-E10-5_chr1'


def plot(target, files, prange, background=True, pvalue=True, differential=False, significant=False,
         truncate=False, min_p=None, max_p=None, levels=None, xfold=None, plot_significant=True, use_queue=True):
    import queue
    import numpy as np
    from hicexplorer.lib.viewpoint import Viewpoint
    viewpoint = Viewpoint()
    background_data = None
    if background:
        background_data = viewpoint.readBackgroundDataFile(BACKGROUND, prange, prange[1], pMean=True)
    args = argparse.Namespace(interactionFile=INTERACTIONS, range=prange, dpi=60, colorList=['g', 'b'],
                              pValue=pvalue, differentialTestResult=DIFFERENTIAL if differential else None,
                              significantInteractions=SIGNIFICANT if significant else None,
                              plotSignificantInteractions=plot_significant, truncateZeroPvalues=truncate,
                              minPValue=min_p, maxPValue=max_p, colorMapPvalue='RdYlBu',
                              pValueSignificanceLevels=levels, xFold=xfold, outputFormat='png')
    differential_list = [[f[0][0], f[1][0], f[1][1], f[1][2], 'rejected'] if len(f) > 1 else [] for f in files]
    significant_list = [[list(sample) for sample in f] for f in files]
    q = queue.Queue() if use_queue else None
    result = target(files, differential_list, background_data, args, viewpoint, significant_list,
                    np.int64(1000), q)
    if q is None:
        return result
    out = q.get()
    if isinstance(out, str):
        return out.split('Traceback')[0]
    images, names = out
    return [image.getvalue() for image in images], names


DUAL = [[[FL, 'chr1', gene], [MB, 'chr1', gene]] for gene in ('Eya1', 'Sox17')]
SINGLE = [[[FL, 'chr1', 'Tfap2d']]]

CASES = [
    Case('plot_images/dual_background', M + 'plot_images', driver='plot', args=(DUAL[:1], [200000, 200000])),
    Case('plot_images/single_no_pvalue', M + 'plot_images', driver='plot',
         args=(SINGLE, [100000, 100000], False, False)),
    Case('plot_images/differential_significant', M + 'plot_images', driver='plot',
         args=(DUAL[:1], [200000, 200000]),
         kwargs={'differential': True, 'significant': True, 'truncate': True, 'min_p': 0.01, 'max_p': 0.5,
                 'xfold': 1.5}),
    Case('plot_images/small_range_kb_ticks', M + 'plot_images', driver='plot', args=(DUAL[1:], [5000, 4000]),
         kwargs={'max_p': 0.2}),
    Case('plot_images/wide_range', M + 'plot_images', driver='plot', args=(SINGLE, [1500000, 1500000], False)),
    Case('plot_images/no_data', M + 'plot_images', driver='plot', args=([[[FL, 'chr1', 'Nope']]], [200000, 200000])),
    Case('plot_images/no_queue', M + 'plot_images', driver='plot', args=(SINGLE, [20000, 20000]),
         kwargs={'use_queue': False}),
    Case('plot_images/failure', M + 'plot_images', driver='plot', args=([[[FL, 'chr1']]], [20000, 20000])),
]
