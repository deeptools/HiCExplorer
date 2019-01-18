from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=UserWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import argparse
import os
import numpy as np
from builtins import range
from past.builtins import map
from scipy.sparse import triu
from scipy.stats import pearsonr, spearmanr

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
# for plotting
from matplotlib import use as mplt_use
mplt_use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Computes pairwise correlations between Hi-C matrices data. '
        'The correlation is computed taking the values from each pair '
        'of matrices and discarding values that are zero in both matrices.'
        'Parameters that strongly affect correlations are bin size of the Hi-C '
        'matrices and the considered range. The smaller the bin size of the '
        'matrices, the finer differences you score. The --range parameter should '
        'be selected at a meaningful genomic scale according to, for example, the '
        'mean size of the TADs in the organism you work with.')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrices', '-m',
                                help='Matrices to correlate (usually .h5 but other formats are allowed). '
                                'hicCorrelate is better used on un-corrected matrices in order to '
                                'exclude any changes introduced by the correction.',
                                nargs='+',
                                required=True)

    parserHeatmap = parser.add_argument_group('Heatmap arguments',
                                              description='Options for generating the correlation heatmap \r')

    # define the arguments
    parserHeatmap.add_argument('--zMin', '-min',
                               default=None,
                               help='Minimum value for the heatmap intensities. '
                               'If not specified the value is set automatically.',
                               type=float)

    parserHeatmap.add_argument('--zMax', '-max',
                               default=None,
                               help='Maximum value for the heatmap intensities.'
                               'If not specified the value is set automatically.',
                               type=float)

    parserHeatmap.add_argument(
        '--colorMap', default='jet',
        metavar='',
        help='Color map to use for the heatmap. Available values can be '
             'seen here: '
             'http://matplotlib.org/examples/color/colormaps_reference.html')

    parserHeatmap.add_argument('--plotFileFormat',
                               metavar='FILETYPE',
                               help='Image format type. If given, this option '
                               'overrides the image format based on the plotFile '
                               'ending. The available options are: png, emf, '
                               'eps, pdf and svg.',
                               choices=['png', 'pdf', 'svg', 'eps', 'emf'])

    parserHeatmap.add_argument('--plotNumbers',
                               help='If set, then the correlation number is plotted '
                               'on top of the heatmap.',
                               action='store_true',
                               required=False)

    parserOpt = parser.add_argument_group('Optional arguments')

    # define the arguments
    parserOpt.add_argument('--method',
                           help='Correlation method to use.',
                           choices=['pearson', 'spearman'],
                           default='pearson')

    parserOpt.add_argument('--log1p',
                           help='If set, then the log1p of the matrix values is used. This parameter has no '
                           'effect for Spearman correlations but changes the output of Pearson correlation '
                           'and, for the scatter plot, if set, the visualization of the values is easier.',
                           action='store_true')

    parserOpt.add_argument('--labels', '-l',
                           metavar='sample1 sample2',
                           help='User defined labels instead of default labels '
                           'from file names. '
                           'Multiple labels have to be separated by space, e.g. '
                           '--labels sample1 sample2 sample3',
                           nargs='+')

    parserOpt.add_argument('--range',
                           help='In bp with the format low_range:high_range, '
                           'for example 1000000:2000000. If --range is given '
                           'only counts within this range are considered. '
                           'The range should be adjusted to the size of interacting  '
                           'domains in the genome you are working with.')

    parserOpt.add_argument('--outFileNameHeatmap', '-oh',
                           help='File name to save the resulting heatmap plot.',
                           required=True)

    parserOpt.add_argument('--outFileNameScatter', '-os',
                           help='File name to save the resulting scatter plot.',
                           required=True)

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the '
                           'correlation.',
                           default=None,
                           nargs='+')

    parserOpt.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module. Is only used with \'cool\' matrix format.'
                           ' One master process which is used to read the input file into the buffer and one process which is merging '
                           'the output bam files of the processes into one output bam file. All other threads do the actual computation.',
                           required=False,
                           default=4,
                           type=int
                           )

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def plot_correlation(corr_matrix, labels, plot_filename, vmax=None,
                     vmin=None, colormap='Reds', image_format=None):
    import scipy.cluster.hierarchy as sch
    num_rows = corr_matrix.shape[0]

    # set the minimum and maximum values
    if vmax is None:
        vmax = 1
    if vmin is None:
        vmin = 0 if corr_matrix.min() >= 0 else -1

    # Compute and plot dendrogram.
    fig = plt.figure(figsize=(10.5, 9.5))
    axdendro = fig.add_axes([0.02, 0.1, 0.1, 0.7])
    axdendro.set_axis_off()
    y_var = sch.linkage(corr_matrix, method='complete')
    z_var = sch.dendrogram(y_var, orientation='left',
                           link_color_func=lambda k: 'black')
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    cmap = plt.get_cmap(colormap)

    # this line simply makes a new cmap, based on the original
    # colormap that goes from 0.0 to 0.9
    # This is done to avoid colors that
    # are too dark at the end of the range that do not offer
    # a good contrast between the correlation numbers that are
    # plotted on black.
    cmap = cmap.from_list(colormap + "clipped", cmap([0.0, 0.8]))
    # Plot distance matrix.
    axmatrix = fig.add_axes([0.13, 0.1, 0.6, 0.7])
    index = z_var['leaves']
    corr_matrix = corr_matrix[index, :]
    corr_matrix = corr_matrix[:, index]
    img_mat = axmatrix.matshow(corr_matrix, aspect='equal', origin='lower',
                               cmap=cmap, extent=(0, num_rows, 0, num_rows),
                               vmax=vmax, vmin=vmin)
    axmatrix.yaxis.tick_right()
    axmatrix.set_yticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_yticklabels(np.array(labels).astype('str')[index],
                             fontsize=14)

    axmatrix.set_xticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_xticklabels(np.array(labels).astype('str')[index],
                             fontsize=14,
                             rotation=45,
                             ha='left')

#    axmatrix.set_xticks([])
    # Plot colorbar.
    axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
    plt.colorbar(img_mat, cax=axcolor, orientation='horizontal')
    for row in range(num_rows):
        for col in range(num_rows):
            axmatrix.text(row + 0.5, col + 0.5,
                          "{:.2f}".format(corr_matrix[row, col]),
                          ha='center', va='center')

    fig.savefig(plot_filename, format=image_format)


def get_vectors(mat1, mat2):
    """
    Uses sparse matrix tricks to convert
    into a vector the matrix values such
    that zero values that appear in only
    one of the matrices is kept. But
    zeros in two matrices are removed

    Requires two sparse matrices as input
    """
    assert mat1.shape == mat2.shape, "Matrices have different shapes. "\
        "Computation of correlation is not possible."

    # create a new matrix that is the sum of the two
    # matrices to compare. The goal is to have
    # a matrix that contains all the positions
    # that are non-zero in both matrices
    _mat = mat1 + mat2

    # add one to each element in the new matrix
    _mat.data += 1

    # get a vector of the values in mat1 from
    # _mat
    values1 = (_mat - mat1).data - 1

    # get a vector of the values in mat2 from
    # _mat
    values2 = (_mat - mat2).data - 1

    return values1, values2


def main(args=None):

    args = parse_arguments().parse_args(args)

    if args.labels and len(args.matrices) != len(args.labels):
        log.error("The number of labels does not match the number of matrices.")
        exit(0)
    if not args.labels:
        args.labels = map(lambda x: os.path.basename(x), args.matrices)

    num_files = len(args.matrices)
    map(lambda x: os.path.basename(x), args.matrices)
    # initialize results matrix
    results = np.zeros((num_files, num_files), dtype='float')

    rows, cols = np.triu_indices(num_files)
    correlation_opts = {'spearman': spearmanr,
                        'pearson': pearsonr}
    hic_mat_list = []
    max_value = None
    min_value = None
    all_mat = None
    all_nan = []

    for i, matrix in enumerate(args.matrices):
        log.info("loading hic matrix {}\n".format(matrix))

        if (check_cooler(args.matrices[i])) and args.chromosomes is not None and len(args.chromosomes) == 1:
            _mat = hm.hiCMatrix(matrix, pChrnameList=args.chromosomes)
        else:
            _mat = hm.hiCMatrix(matrix)
            if args.chromosomes:
                _mat.keepOnlyTheseChr(args.chromosomes)
            _mat.filterOutInterChrCounts()

        _mat.diagflat(0)
        log.info("restore masked bins {}\n".format(matrix))
        bin_size = _mat.getBinSize()
        all_nan = np.unique(np.concatenate([all_nan, _mat.nan_bins]))

        _mat = triu(_mat.matrix, k=0, format='csr')
        if args.range:
            min_dist, max_dist = args.range.split(":")
            min_dist = int(min_dist)
            max_dist = int(max_dist)
            if max_dist < bin_size:
                log.error("Please specify a max range that is larger than bin size ({})".format(bin_size))
                exit()
            max_depth_in_bins = int(max_dist / bin_size)
            max_dist = int(max_dist) // bin_size
            min_dist = int(min_dist) // bin_size
            # work only with the upper matrix
            # and remove all pixels that are beyond
            # max_depth_in_bis
            # (this is done by subtracting a second sparse matrix
            # that contains only the upper matrix that wants to be removed.
            _mat = triu(_mat, k=0, format='csr') - triu(_mat, k=max_depth_in_bins, format='csr')

            _mat.eliminate_zeros()

            _mat_coo = _mat.tocoo()
            dist = _mat_coo.col - _mat_coo.row
            keep = np.flatnonzero((dist <= max_dist) & (dist >= min_dist))
            _mat_coo.data = _mat_coo.data[keep]
            _mat_coo.row = _mat_coo.row[keep]
            _mat_coo.col = _mat_coo.col[keep]
            _mat = _mat_coo.tocsr()
        else:
            _mat = triu(_mat, k=0, format='csr')

        if args.log1p:
            _mat.data = np.log1p(_mat.data)
        if all_mat is None:
            all_mat = _mat
        else:
            all_mat = all_mat + _mat

        if max_value is None or max_value < _mat.data.max():
            max_value = _mat.data.max()
        if min_value is None or min_value > _mat.data.min():
            min_value = _mat.data.min()

        hic_mat_list.append(_mat)

    # remove nan bins
    rows_keep = cols_keep = np.delete(list(range(all_mat.shape[1])), all_nan)
    all_mat = all_mat[rows_keep, :][:, cols_keep]

    # make large matrix to correlate by
    # using sparse matrix tricks

    big_mat = None
    for mat in hic_mat_list:
        mat = mat[rows_keep, :][:, cols_keep]
        sample_vector = (mat + all_mat).data - all_mat.data
        if big_mat is None:
            big_mat = sample_vector
        else:
            big_mat = np.vstack([big_mat, sample_vector])

    # take the transpose such that columns represent each of the samples
    big_mat = np.ma.masked_invalid(big_mat).T

    grids = gridspec.GridSpec(num_files, num_files)
    grids.update(wspace=0, hspace=0)
    fig = plt.figure(figsize=(2 * num_files, 2 * num_files))
    plt.rcParams['font.size'] = 8.0

    min_value = int(big_mat.min())
    max_value = int(big_mat.max())
    if (min_value % 2 == 0 and max_value % 2 == 0) or \
            (min_value % 1 == 0 and max_value % 2 == 1):
        # make one value odd and the other even
        max_value += 1

    if args.log1p:
        major_locator = FixedLocator(list(range(min_value, max_value, 2)))
        minor_locator = FixedLocator(list(range(min_value, max_value, 1)))

    for index in range(len(rows)):
        row = rows[index]
        col = cols[index]
        if row == col:
            results[row, col] = 1

            # add titles as
            # empty plot in the diagonal
            ax = fig.add_subplot(grids[row, col])
            ax.text(0.6, 0.6, args.labels[row],
                    verticalalignment='center',
                    horizontalalignment='center',
                    fontsize=10, fontweight='bold',
                    transform=ax.transAxes)
            ax.set_axis_off()
            continue

        log.info("comparing {} and {}\n".format(args.matrices[row],
                                                args.matrices[col]))

        # remove cases in which both are zero or one is zero and
        # the other is one
        _mat = big_mat[:, [row, col]]
        _mat = _mat[_mat.sum(axis=1) > 1, :]
        vector1 = _mat[:, 0]
        vector2 = _mat[:, 1]

        results[row, col] = correlation_opts[args.method](vector1, vector2)[0]

        # scatter plots
        ax = fig.add_subplot(grids[row, col])
        if args.log1p:
            ax.xaxis.set_major_locator(major_locator)
            ax.xaxis.set_minor_locator(minor_locator)
            ax.yaxis.set_major_locator(major_locator)
            ax.yaxis.set_minor_locator(minor_locator)

        ax.text(0.2, 0.8, "{}={:.2f}".format(args.method,
                                             results[row, col]),
                horizontalalignment='left',
                transform=ax.transAxes)
        ax.get_yaxis().set_tick_params(
            which='both',
            left='off',
            right='off',
            direction='out')

        ax.get_xaxis().set_tick_params(
            which='both',
            top='off',
            bottom='off',
            direction='out')

        if col != num_files - 1:
            ax.set_yticklabels([])
        else:
            ax.yaxis.tick_right()
            ax.get_yaxis().set_tick_params(
                which='both',
                left='off',
                right='on',
                direction='out')
        if col - row == 1:
            ax.xaxis.tick_bottom()
            ax.get_xaxis().set_tick_params(
                which='both',
                top='off',
                bottom='on',
                direction='out')
        else:
            ax.set_xticklabels([])

        ax.hist2d(vector1, vector2, bins=150, cmin=0.1)
    fig.tight_layout()
    log.info("saving {}".format(args.outFileNameScatter))
    fig.savefig(args.outFileNameScatter, bbox_inches='tight')

    results = results + np.triu(results, 1).T
    plot_correlation(results, args.labels,
                     args.outFileNameHeatmap,
                     args.zMax,
                     args.zMin,
                     args.colorMap,
                     image_format=args.plotFileFormat)
#                    plot_numbers=args.plotNumbers)
