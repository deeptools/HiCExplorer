"""Equivalence cases for hicexplorer.hicPlotMatrix (apieq.py).

The data functions run on hiCMatrix objects of the repository's matrices
(small_test_matrix_50kb_res as cool and h5, Li_et_al_2015, a cool region
whose bins do not start at 0). The drawing functions plotPerChr and
plotBigwig stay matplotlib code but take their data from C++, so their cases
compare what ends up on the axes (the heatmap values, the bigwig polygons,
titles and limits) and the log records; both interpreters draw with the same
matplotlib. exit() raises SystemExit, which the drivers return as a value.

Deviations from 3.7.6 (cpp/python/src/api/hicPlotMatrix.cpp):
translate_region and the bigwig file names take str (other types raise the
TypeError or AttributeError 3.7.6 raised); a negative number of bigwig bins
raises ValueError (pyBigWig does not return).
"""
import argparse
import os

import numpy as np

from apieq import Case

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicPlotMatrix.'
COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
LI = os.path.join(TEST_DATA, 'Li_et_al_2015.cool')
OBS_EXP = os.path.join(TEST_DATA, 'hicPCA', 'obsexp_norm.h5')
TRACK = os.path.join(TEST_DATA, 'bigwig_chrx_2e6_5e6.bw')
DM3_TRACK = os.path.join(TEST_DATA, 'hicPCA', 'pca1.bw')


def load(path, chromosomes=None, keep=None):
    from hicmatrix import HiCMatrix as hm
    hic = hm.hiCMatrix(path, pChrnameList=chromosomes)
    if keep is not None:
        hic.keepOnlyTheseChr(keep)
    return hic


def guarded(function):
    try:
        return function()
    except SystemExit as error:
        return ('SystemExit', error.code)


def on_matrix(target, path, *args, **kwargs):
    chromosomes = kwargs.pop('chromosomes', None)
    keep = kwargs.pop('keep', None)
    hic = load(path, chromosomes, keep)
    return guarded(lambda: target(hic, *args))


def translate(target, path, region, chromosomes=None):
    hic = load(path, chromosomes)
    return guarded(lambda: target(region, hic))


def region_of(target, path, region, region2=None, matrix_name=None):
    hic = load(path)
    args = argparse.Namespace(region=region, matrix=matrix_name or path)
    if region2 is not None:
        args.region2 = region2
    result = guarded(lambda: target(args, hic))
    return result, args.region


def plot_args(**kwargs):
    defaults = dict(bigwig=None, increaseFigureHeight=0.5, increaseFigureWidth=0.5, dpi=72,
                    bigwigAdditionalVerticalAxis=False, log=False, log1p=False, vMin=None, vMax=None,
                    region=None, region2=None, matrix=H5, title=None, rotationX=0, rotationY=0, fontsize=10,
                    scoreName=None, chromosomeOrder=None, flipBigwigSign=False, scaleFactorBigwig=1.0,
                    vMinBigwig=None, vMaxBigwig=None)
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def figure_content(fig):
    """What was drawn: per axes its title, labels, limits, the arrays of its
    heatmaps and the vertices of its filled areas."""
    out = []
    for ax in fig.axes:
        meshes = [np.asarray(c.get_array()) for c in ax.collections if type(c).__name__ == 'QuadMesh']
        polygons = [p.vertices for c in ax.collections if type(c).__name__ in ('PolyCollection', 'FillBetweenPolyCollection')
                    for p in c.get_paths()]
        out.append({'title': ax.get_title(), 'xlabel': ax.get_xlabel(), 'ylabel': ax.get_ylabel(),
                    'xlim': tuple(float(v) for v in ax.get_xlim()), 'ylim': tuple(float(v) for v in ax.get_ylim()),
                    'meshes': meshes, 'polygons': polygons})
    return out


def per_chromosome(target, path, keep, bigwig=None, **kwargs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    hic = load(path, keep=keep)
    args = plot_args(matrix=path, bigwig=None if bigwig is None else [bigwig], **kwargs)
    try:
        fig = target(hic, 'viridis', args, bigwig is not None, 50000)
    except SystemExit as error:
        return ('SystemExit', error.code)
    content = figure_content(fig)
    plt.close(fig)
    return content, args.region


def bigwig(target, files, sizes_matrix=OBS_EXP, region=None, vertical=False, resolution=5000, keep=None, **kwargs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import importlib
    module = importlib.import_module('hicexplorer.hicPlotMatrix')
    sizes = module.chromosome_start_end(load(sizes_matrix, keep=keep)) if sizes_matrix else None
    fig, axes = plt.subplots(max(len(files), 1), 1, squeeze=False)
    axes = list(axes[:, 0])
    try:
        returned = target(axes, files, pChromosomeSizes=sizes, pRegion=region, pVertical=vertical,
                          pResolution=resolution, **kwargs)
    except SystemExit as error:
        returned = ('SystemExit', error.code)
    content = [(a.get_frame_on(), a.xaxis.get_visible(), a.yaxis.get_visible()) for a in axes]
    content = (returned, content, figure_content(fig))
    plt.close(fig)
    return content


CASES = [
    # relabel_ticks
    Case('ticks/mbp', M + 'relabel_ticks', args=lambda: (np.array([0., 1e6, 2e6, 3.5e6]),)),
    Case('ticks/kbp', M + 'relabel_ticks', args=lambda: (np.array([0., 500., 1000., 2000.]),)),
    Case('ticks/bp', M + 'relabel_ticks', args=lambda: (np.array([-10.5, 0., 100.125, 1000.]),)),
    Case('ticks/rounding', M + 'relabel_ticks', args=lambda: (np.array([1500.5, 2500.5, 1600500., 1700500.]),)),
    Case('ticks/list_of_ints', M + 'relabel_ticks', args=([0, 2000000, 4000000],)),
    Case('ticks/nan_inf', M + 'relabel_ticks', args=lambda: (np.array([np.nan, -np.inf, np.inf, 2e6]),)),
    Case('ticks/one', M + 'relabel_ticks', args=lambda: (np.array([2e6]),)),
    Case('ticks/empty_array', M + 'relabel_ticks', args=lambda: (np.array([]),)),
    Case('ticks/empty_list', M + 'relabel_ticks', args=([],)),

    # translate_region
    Case('region/cool', M + 'translate_region', driver='translate', args=(COOL, 'chr2L:1,000,000-2,000,000')),
    Case('region/h5', M + 'translate_region', driver='translate', args=(H5, 'chr2L:1000000-2000000')),
    Case('region/chromosome_only_h5', M + 'translate_region', driver='translate', args=(H5, 'chr3R')),
    Case('region/chromosome_only_cool', M + 'translate_region', driver='translate', args=(COOL, 'chrX')),
    Case('region/start_only', M + 'translate_region', driver='translate', args=(H5, 'chrX:500')),
    Case('region/other_name', M + 'translate_region', driver='translate', args=(COOL, '2L:10;0-2!00')),
    Case('region/li', M + 'translate_region', driver='translate', args=(LI, 'chrX')),
    Case('region/cool_region_load', M + 'translate_region', driver='translate',
         args=(COOL, 'chr2L', ['chr2L:1000000-2000000'])),
    Case('region/unknown', M + 'translate_region', driver='translate', args=(COOL, 'chrZ:1-2')),
    Case('region/not_a_number', M + 'translate_region', driver='translate', args=(COOL, 'chr2L:1x-2')),
    Case('region/bytes', M + 'translate_region', driver='translate', args=(COOL, b'chr2L')),
    Case('region/int', M + 'translate_region', driver='translate', args=(COOL, 5)),

    # getRegion
    Case('get_region/cool', M + 'getRegion', driver='region_of', args=(COOL, 'chr2L:1000001-2000001')),
    Case('get_region/h5', M + 'getRegion', driver='region_of', args=(H5, 'chr2L:1000001-2000001')),
    Case('get_region/chromosome', M + 'getRegion', driver='region_of', args=(H5, 'chr4')),
    Case('get_region/region2_cool', M + 'getRegion', driver='region_of',
         args=(COOL, 'chr2L:1000000-2000000', 'chr3R:999999-1500001')),
    Case('get_region/region2_h5', M + 'getRegion', driver='region_of',
         args=(H5, 'chr2L:1000000-2000000', 'chr3R:999999-1500001')),
    Case('get_region/region2_other_name', M + 'getRegion', driver='region_of', args=(H5, 'chrX', '3R')),
    Case('get_region/region2_unknown', M + 'getRegion', driver='region_of', args=(H5, 'chrX', 'chrZ')),
    Case('get_region/region2_empty', M + 'getRegion', driver='region_of', args=(H5, 'chrX', '')),
    Case('get_region/li_cool', M + 'getRegion', driver='region_of', args=(LI, 'X:100000-200000')),
    Case('get_region/no_bins', M + 'getRegion', driver='region_of', args=(H5, 'chrX:10-20')),
    Case('get_region/cool_named_h5', M + 'getRegion', driver='region_of',
         args=(COOL, 'chrX:10-20', None, 'matrix.h5')),

    # chromosome_start_end
    Case('start_end/cool', M + 'chromosome_start_end', driver='on_matrix', args=(COOL,)),
    Case('start_end/h5', M + 'chromosome_start_end', driver='on_matrix', args=(H5,)),
    Case('start_end/li', M + 'chromosome_start_end', driver='on_matrix', args=(LI,)),
    Case('start_end/region', M + 'chromosome_start_end', driver='on_matrix', args=(COOL,),
         kwargs={'chromosomes': ['chr2L:1000000-2000000']}),

    # make_start_pos_array
    Case('start_pos/cool', M + 'make_start_pos_array', driver='on_matrix', args=(COOL,)),
    Case('start_pos/h5', M + 'make_start_pos_array', driver='on_matrix', args=(H5,)),
    Case('start_pos/li', M + 'make_start_pos_array', driver='on_matrix', args=(LI,)),
    Case('start_pos/region', M + 'make_start_pos_array', driver='on_matrix', args=(COOL,),
         kwargs={'chromosomes': ['chr2L:1000000-2000000']}),
    Case('start_pos/subset', M + 'make_start_pos_array', driver='on_matrix', args=(H5,),
         kwargs={'keep': ['chrX', 'chr2L', 'chr4']}),

    # plotPerChr (drawing on C++ data)
    Case('per_chr/plain', M + 'plotPerChr', driver='per_chromosome', args=(H5, ['chr4', 'chrXHet', 'chr2LHet'])),
    Case('per_chr/log1p', M + 'plotPerChr', driver='per_chromosome', args=(H5, ['chr4', 'chr2LHet']),
         kwargs={'log1p': True}),
    Case('per_chr/log', M + 'plotPerChr', driver='per_chromosome', args=(COOL, ['chr4', 'chrM', 'chrYHet']),
         kwargs={'log': True}),
    Case('per_chr/log_with_data', M + 'plotPerChr', driver='per_chromosome', args=(COOL, ['chr4', 'chr2LHet']),
         kwargs={'log': True}),
    Case('per_chr/many', M + 'plotPerChr', driver='per_chromosome',
         args=(COOL, ['chr4', 'chrXHet', 'chr2LHet', 'chrYHet', 'chrM', 'chr3LHet'])),
    Case('per_chr/bigwig', M + 'plotPerChr', driver='per_chromosome', args=(H5, ['chr4', 'chr2LHet'], DM3_TRACK)),
    Case('per_chr/bigwig_vertical', M + 'plotPerChr', driver='per_chromosome',
         args=(H5, ['chr4'], DM3_TRACK), kwargs={'bigwigAdditionalVerticalAxis': True, 'log1p': True}),

    # plotBigwig (drawing on C++ data)
    Case('bigwig/region', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 3000000]}),
    Case('bigwig/region_beyond_chromosome', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 22000000, 10 ** 15]}),
    Case('bigwig/region_vertical', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 3000000], 'vertical': True, 'pFlipBigwigSign': True}),
    Case('bigwig/region_other_name', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['X', 2500000, 2600000], 'sizes_matrix': LI, 'resolution': 1000}),
    Case('bigwig/region_scaled_clipped', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 4000000], 'pScaleFactorBigwig': 2.5, 'pValueMin': 1.0,
                 'pValueMax': 6.0, 'resolution': 20000}),
    Case('bigwig/region_min_only', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 4000000], 'pValueMin': 2.0, 'resolution': 20000}),
    Case('bigwig/region_too_small', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 2001000]}),
    Case('bigwig/genome', M + 'plotBigwig', driver='bigwig', args=([TRACK],), kwargs={'resolution': 50000}),
    Case('bigwig/genome_one_chromosome', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'resolution': 50000, 'keep': ['chrX'], 'pValueMax': 3.0}),
    Case('bigwig/genome_dm3', M + 'plotBigwig', driver='bigwig', args=([DM3_TRACK],),
         kwargs={'sizes_matrix': H5, 'resolution': 100000, 'pFlipBigwigSign': True}),
    Case('bigwig/genome_vertical', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'resolution': 50000, 'vertical': True}),
    Case('bigwig/genome_missing_chromosome', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'sizes_matrix': H5, 'resolution': 100000}),
    Case('bigwig/two_files', M + 'plotBigwig', driver='bigwig', args=([DM3_TRACK, TRACK],),
         kwargs={'sizes_matrix': H5, 'resolution': 100000}),
    Case('bigwig/two_files_too_coarse', M + 'plotBigwig', driver='bigwig', args=([DM3_TRACK, TRACK],),
         kwargs={'sizes_matrix': H5, 'resolution': 500000}),
    Case('bigwig/no_sizes', M + 'plotBigwig', driver='bigwig', args=([TRACK],), kwargs={'sizes_matrix': None}),
    Case('bigwig/not_bigwig', M + 'plotBigwig', driver='bigwig', args=([TRACK, 'track.bedgraph'],)),
    Case('bigwig/no_files', M + 'plotBigwig', driver='bigwig', args=([],)),
    Case('bigwig/missing_file', M + 'plotBigwig', driver='bigwig', args=([os.path.join(TEST_DATA, 'none.bw')],)),
    Case('bigwig/no_resolution', M + 'plotBigwig', driver='bigwig', args=([TRACK],),
         kwargs={'region': ['chrX', 2000000, 3000000], 'resolution': None}),
]
