from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)
from scipy.sparse import load_npz
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate
from mpl_toolkits.axes_grid1 import make_axes_locatable


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
        hicPlotAverage regions plots the data computed by hicAverageRegions. It shows the summed up and averaged regions around
        all given reference points. This is useful to determine the difference between samples if the TAD configuration is equal or changed.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The averaged regions file computed by hicAverageRegions, ends npz.',
                                required=True)
    parserRequired.add_argument('--outputFile', '-o',
                                help='The averaged regions plot.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--log1p',
                           help='Plot the log1p of the matrix values.',
                           action='store_true')

    parserOpt.add_argument('--log',
                           help='Plot the log of the matrix values.',
                           action='store_true')
    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='hot_r')
    parserOpt.add_argument('--vMin',
                           help='Minimum score value.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--vMax',
                           help='Maximum score value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--dpi',
                           help='Resolution for the image in case the'
                           'ouput is a raster graphics image (e.g png, jpg).',
                           type=int,
                           default=300)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    matrix = load_npz(args.matrix)

    matrix = matrix.toarray()
    matrix = np.triu(matrix)
    matrix = rotate(matrix, 45)
    matrix_shapes = matrix.shape
    matrix = matrix[:matrix_shapes[0] // 2, :]
    norm = None
    if args.vMax is not None or args.vMin is not None:
        matrix = matrix.clip(min=args.vMin, max=args.vMax)
    if args.log1p or args.log:

        mask = matrix == 0
        mask_nan = np.isnan(matrix)
        mask_inf = np.isinf(matrix)

        try:
            matrix[mask] = np.nanmin(matrix[mask == False])
            matrix[mask_nan] = np.nanmin(matrix[mask_nan == False])
            matrix[mask_inf] = np.nanmin(matrix[mask_inf == False])
        except Exception:
            log.warning('Clearing of matrix failed. Plotting can fail.')
        if args.log1p:
            matrix += 1
            norm = LogNorm()

        elif args.log:
            norm = LogNorm()

    fig = plt.figure()
    axis = plt.gca()
    matrix_axis = axis.matshow(matrix, cmap=args.colorMap, norm=norm)
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)

    fig.colorbar(matrix_axis, cax=cax)
    plt.tight_layout()
    plt.savefig(args.outputFile, dpi=args.dpi)
