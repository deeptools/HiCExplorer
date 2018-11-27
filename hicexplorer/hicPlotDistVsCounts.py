from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
import numpy as np
import pandas as pd
import argparse

from hicmatrix import HiCMatrix
from hicexplorer._version import __version__

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from collections import OrderedDict
from builtins import range
from past.builtins import zip
from six import iteritems

from scipy.sparse import triu


import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        add_help=False,
        description='This program creates distance vs. Hi-C counts plots. It can use several matrix files to compare '
                    'them at once. If the `--perchr` option is given, each chromosome is plotted independently. '
                    'When plotting multiple matrices, denser matrices are scaled down to match the sum of the smallest matrix.')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrices', '-m',
                                help='Hi-C normalized (corrected) matrices. Each path should be separated by a space.',
                                nargs="+",
                                required=True)

    parserRequired.add_argument('--plotFile', '-o',
                                help='File name to save the file. The given file '
                                'ending will be used '
                                'to determine the image format. '
                                'The available options are: .png, .emf, '
                                '.eps, .pdf and .svg.',
                                type=argparse.FileType('w'),
                                metavar='file name',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--labels',
                           help='Label to assign to each matrix file. Each label should be separated by a space. Quote '
                           'labels that contain spaces: E.g. --labels label1 "labels 2". If no labels are given '
                           'then the file name is used.',
                           nargs="+")

    parserOpt.add_argument('--skipDiagonal', '-s',
                           help='If set, diagonal counts are not included.',
                           action='store_true')

    parserOpt.add_argument('--maxdepth',
                           help='Maximum distance from diagonal to use. In other words, distances up to maxDeph are '
                           'computed. Default is 3 million bp.',
                           metavar='INT bp',
                           type=int,
                           default=int(3e6))

    parserOpt.add_argument('--perchr',
                           help='If given, computes and display distance versus Hi-C counts plots for each chromosome stored '
                           'in the matrices passed to --matrices.',
                           action='store_true')

    parserOpt.add_argument('--chromosomeExclude',
                           help='Exclude the given list of chromosomes. This is useful for example to exclude '
                           'the Y chromosome. The names of the chromosomes should be separated by space.',
                           nargs='+')

    parserOpt.add_argument('--outFileData',
                           help='If given, the data underlying the plots is saved on this file.',
                           type=argparse.FileType('w'),
                           )

    parserOpt.add_argument('--plotsize',
                           help='Width and height of the plot (in inches). Default is 6*number of cols, 4 * number of '
                           'rows. The maximum number of rows is 4. Example: --plotsize 6 5',
                           nargs=2,
                           type=float
                           )

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_distance_mean(hicmat, maxdepth=None, perchr=False):
    """
    Converts a corrected counts matrix into a
    obs / expected matrix or z-scores fast.

    The caveat is that the obs/exp or z-score are only
    computed for non-zero values, although zero values that
    are not part of the sparse matrix are considered.

    For each diagonal the mean (and std when computing z-scores) are
    calculated and then each non-zero value of the sparse matrix is
    replaced by the obs/exp or z-score.

    Parameters
    ----------
    hicmat: HiCMatrix object
    maxdepth: maximum distance from the diagonal to consider. All contacts beyond this distance will not
                     be considered.
    perchr: bool to indicate if computations should be perform per chromosome


    Returns
    -------
    observed / expected sparse matrix

    >>> from scipy.sparse import csr_matrix, dia_matrix
    >>> row, col = np.triu_indices(5)
    >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    >>> hic = HiCMatrix.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 1,  8,  5, 3, 0],
    ... [ 0,  4, 15, 5, 1],
    ... [ 0,  0,  0, 7, 2],
    ... [ 0,  0,  0, 0, 1],
    ... [ 0,  0,  0, 0, 0]])

    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals)
    >>> hic.convert_to_obs_exp_matrix().todense()
    matrix([[1. , 0.8, 1. , 1. , 0. ],
            [0. , 4. , 1.5, 1. , 1. ],
            [0. , 0. , 0. , 0.7, 2. ],
            [0. , 0. , 0. , 0. , 1. ],
            [0. , 0. , 0. , 0. , 0. ]])

    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.convert_to_obs_exp_matrix(maxdepth=20).todense()
    matrix([[1. , 0.8, 1. , 0. , 0. ],
            [0. , 4. , 1.5, 1. , 0. ],
            [0. , 0. , 0. , 0.7, nan],
            [0. , 0. , 0. , 0. , nan],
            [0. , 0. , 0. , 0. , 0. ]])

    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.convert_to_obs_exp_matrix(zscore=True).todense()
    matrix([[ 0.        , -0.56195149,         nan,         nan, -1.41421356],
            [ 0.        ,  1.93649167,  1.40487872,         nan,  0.        ],
            [ 0.        ,  0.        , -0.64549722, -0.84292723,  1.41421356],
            [ 0.        ,  0.        ,  0.        , -0.64549722,  0.        ],
            [ 0.        ,  0.        ,  0.        ,  0.        , -0.64549722]])

    nans occur where the standard deviation is zero
    """

    binsize = hicmat.getBinSize()

    if maxdepth:
        if maxdepth < binsize:
            exit("Please specify a maxDepth larger than bin size ({})".format(binsize))

        max_depth_in_bins = int(float(maxdepth * 1.5) / binsize)
        # work only with the upper matrix
        # and remove all pixels that are beyond
        # max_depth_in_bis
        # (this is done by subtracting a second sparse matrix
        # that contains only the upper matrix that wants to be removed.
        hicmat.matrix = triu(hicmat.matrix, k=0, format='csr') - \
            triu(hicmat.matrix, k=max_depth_in_bins, format='csr')
    else:
        hicmat.matrix = triu(hicmat.matrix, k=0, format='csr')

    hicmat.matrix.eliminate_zeros()

    chr_submatrix = OrderedDict()
    cut_intervals = OrderedDict()
    chrom_sizes = OrderedDict()
    chrom_range = OrderedDict()
    if perchr:
        for chrname in hicmat.getChrNames():
            chr_range = hicmat.getChrBinRange(chrname)
            chr_submatrix[chrname] = hicmat.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]].tocoo()
            cut_intervals[chrname] = [hicmat.cut_intervals[x] for x in range(chr_range[0], chr_range[1])]
            chrom_sizes[chrname] = [chr_submatrix[chrname].shape[0]]
            chrom_range[chrname] = (chr_range[0], chr_range[1])

    else:
        chr_submatrix['all'] = hicmat.matrix.tocoo()
        cut_intervals['all'] = hicmat.cut_intervals
        chrom_sizes['all'] = np.array([v[1] - v[0] for k, v in iteritems(hicmat.chrBinBoundaries)])
        chrom_range['all'] = (0, hicmat.matrix.shape[0])

    mean_dict = {}

    for chrname, submatrix in iteritems(chr_submatrix):
        log.info("processing chromosome {}\n".format(chrname))

        dist_list, chrom_list = hicmat.getDistList(submatrix.row, submatrix.col,
                                                   HiCMatrix.hiCMatrix.fit_cut_intervals(cut_intervals[chrname]))

        # to get the sum of all values at a given distance I use np.bincount which
        # is quite fast. However, the input of bincount is positive integers. Moreover
        # it returns the sum for every consecutive integer, even if this is not on the list.
        # Thus, dist_list, which contains the distance in bp between any two bins is
        # converted to bin distance.

        # Because positive integers are needed we add +1 to all bin distances
        # such that the value of -1 (which means different chromosomes) can now be used

        dist_list[dist_list == -1] = -binsize
        # divide by binsize to get a list of bin distances and add +1 to remove negative values
        dist_list = (np.array(dist_list).astype(float) / binsize).astype(int) + 1

        # for each distance, return the sum of all values
        sum_counts = np.bincount(dist_list, weights=submatrix.data)
        distance_len = np.bincount(dist_list)
        # compute the average for each distance
        mat_size = submatrix.shape[0]
        # compute mean value for each distance
        mu = {}
        zero_value_bins = []
        for bin_dist_plus_one, sum_value in enumerate(sum_counts):
            if maxdepth and bin_dist_plus_one == 0:  # this is for intra chromosomal counts
                # when max depth is set, the computation
                # of the total_intra is not accurate and is safer to
                # output np.nan
                mu[bin_dist_plus_one] = np.nan
                continue

            if bin_dist_plus_one == 0:
                total_intra = mat_size ** 2 - sum([size ** 2 for size in chrom_sizes[chrname]])
                diagonal_length = total_intra / 2
            else:
                # to compute the average counts per distance we take the sum_counts and divide
                # by the number of values on the respective diagonal
                # which is equal to the size of each chromosome - the diagonal offset (for those
                # chromosome larger than the offset)
                # In the following example with two chromosomes
                # the first (main) diagonal has a size equal to the matrix (6),
                # while the next has 1 value less for each chromosome (4) and the last one has only 2 values

                # 0 1 2 . . .
                # - 0 1 . . .
                # - - 0 . . .
                # . . . 0 1 2
                # . . . - 0 1
                # . . . - - 0

                # idx - 1 because earlier the values where
                # shifted.
                diagonal_length = sum([size - (bin_dist_plus_one - 1) for size in chrom_sizes[chrname]
                                       if size > (bin_dist_plus_one - 1)])

            # the diagonal length should contain the number of values at a certain distance.
            # If the matrix is dense, the distance_len[bin_dist_plus_one] correctly contains the number of values
            # If the matrix is equally spaced, then, the diagonal_length as computed before is accurate.
            # But, if the matrix is both sparse and with unequal bins, then none of the above methods is
            # accurate but the the diagonal_length as computed before will be closer.
            diagonal_length = max(diagonal_length, distance_len[bin_dist_plus_one])

            if diagonal_length == 0:
                mu[bin_dist_plus_one] = np.nan
            else:
                mu[bin_dist_plus_one] = np.float64(sum_value) / diagonal_length
                if sum_value == 0:
                    zero_value_bins.append(bin_dist_plus_one)
                    log.info("zero value for {}, diagonal len: {}\n".format(bin_dist_plus_one, diagonal_length))
                if len(zero_value_bins) > 10:
                    diff = np.diff(zero_value_bins)
                    if len(diff[diff == 1]) > 10:
                        # if too many consecutive bins with zero are found that means that probably no
                        # further counts will be found
                        log.info("skipping rest of chromosome {}. Too many emtpy diagonals\n".format(chrname))
                        break
            if np.isnan(sum_value):
                log.info("nan value found for distance {}\n".format((bin_dist_plus_one - 1) * binsize))

        if maxdepth is None:
            maxdepth = np.inf
        mean_dict[chrname] = OrderedDict([((k - 1) * binsize, v) for k, v in iteritems(mu) if k > 0 and
                                          (k - 1) * binsize <= maxdepth])
        # mean_dict[chrname]['intra_chr'] = mu[0]

    return mean_dict


def main(args=None):
    """
    for each distance, compare the
    distribution of two samples,
    report number of cases were they differ
    """

    args = parse_arguments().parse_args(args)
    mean_dict = OrderedDict()
    matrix_sum = {}
    if args.labels is None:
        labels = OrderedDict([(x, os.path.basename(x)) for x in args.matrices])
    else:
        labels = OrderedDict(zip(args.matrices, args.labels))

    chroms = set()
    for matrix_file in args.matrices:
        hic_ma = HiCMatrix.hiCMatrix(matrix_file)
        matrix_sum[matrix_file] = hic_ma.matrix.sum()
        if args.chromosomeExclude is None:
            args.chromosomeExclude = []

        chrtokeep = [x for x in list(hic_ma.interval_trees) if x not in args.chromosomeExclude]
        hic_ma.keepOnlyTheseChr(chrtokeep)

        mean_dict[matrix_file] = compute_distance_mean(hic_ma, maxdepth=args.maxdepth, perchr=args.perchr)
        chroms = chroms.union([k for k in list(mean_dict[matrix_file]) if len(mean_dict[matrix_file][k]) > 1])

    # compute scale factors such that values are comparable
    min_sum = min(matrix_sum.values())
    scale_factor = dict([(matrix_file, float(min_sum) / mat_sum) for matrix_file, mat_sum in iteritems(matrix_sum)])
    log.info("The scale factors used are: {}".format(scale_factor))
    if len(args.matrices) > 1 and args.perchr:
        # in this case, for each chromosome a plot is made that combines the data from the
        # hic matrices
        max_cols = 4
        num_rows = int(np.ceil(float(len(chroms)) / max_cols))
        num_cols = min(len(chroms), max_cols)

    else:
        num_cols = num_rows = 1

    if args.plotsize is None:
        width = 6
        height = 4
    else:
        width, height = args.plotsize
    fig = plt.figure(figsize=(width * num_cols, height * num_rows))

    axs = np.empty((num_rows, num_cols), dtype='object')
    for matrix_file in args.matrices:
        idx = 0
        for chrom, mean_values in iteritems(mean_dict[matrix_file]):
            if len(mean_values) <= 1:
                log.debug("No values found for: {}, chromosome: {}\n".format(matrix_file, chrom))
                continue
            x, y = zip(*[(k, v) for k, v in iteritems(mean_values) if v > 0])
            if len(x) <= 1:
                log.debug("No values found for: {}, chromosome: {}\n".format(matrix_file, chrom))
                continue
            if args.perchr and len(args.matrices) == 1:
                col = 0
                row = 0
            else:
                col = idx % num_cols
                row = idx // num_cols
            if axs[row, col] is None:
                ax = plt.subplot2grid((num_rows, num_cols), (row, col))
                ax.set_xlabel('genomic distance')
                ax.set_ylabel('corrected Hi-C counts')
                try:
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                except ValueError:
                    continue
            else:
                ax = axs[row, col]
            y = np.array(y) * scale_factor[matrix_file]
            if args.perchr and len(args.matrices) > 1:
                label = labels[matrix_file]
                ax.set_title(chrom)
            elif args.perchr:
                label = chrom
            else:
                label = labels[matrix_file]

            ax.plot(x, y, label=label)
            axs[row, col] = ax
            idx += 1
            if args.outFileData is not None:
                x_vals = np.stack(x).T
                y_vals = np.stack(y).T
                table_to_export = pd.DataFrame({'Matrix': labels[matrix_file],
                                                'Chromosome': chrom,
                                                'Distance': x_vals,
                                                'Contacts': y_vals})
                table_to_export.to_csv(args.outFileData, sep='\t')

    for ax in axs.reshape(-1):
        if ax is None:
            continue
        ax.legend(prop={'size': 'small'})
        ax.set_xlim(0, args.maxdepth)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(args.plotFile.name, bbox_inches='tight', bbox_extra_artists=(lgd,))
    plt.close(fig)
