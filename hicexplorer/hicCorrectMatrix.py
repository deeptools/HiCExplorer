"""hicCorrectMatrix: a thin layer over the C++ executable hicCorrectMatrix.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
import sys
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from past.builtins import zip
from scipy.sparse import lil_matrix

from hicexplorer.iterativeCorrection import iterativeCorrection
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.utilities import check_cooler

# Knight-Ruiz algorithm:
from krbalancing import *

# packages needed for plotting:
from matplotlib import use
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib
import numpy as np
debug = 0

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicCorrectMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCorrectMatrix')


def correct_subparser():
    # define the arguments
    parser = argparse.ArgumentParser(add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Name of the Hi-C matrix to correct in .h5 format.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. '
                                'The output is a .h5 file.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--correctionMethod',
                           help='Method to be used for matrix correction. It'
                           ' can be set to KR or ICE'
                           ' (Default: %(default)s).',
                           type=str,
                           metavar='STR',
                           default='KR',
                           choices=['KR', 'ICE'])

    parserOpt.add_argument('--filterThreshold', '-t',
                           help='Removes bins of low or large coverage. '
                                'Usually these bins do not contain valid '
                                'Hi-C data or represent regions that '
                                'accumulate reads and thus must be '
                                'discarded. Use hicCorrectMatrix '
                                'diagnostic_plot to identify the modified '
                                'z-value thresholds. A lower and upper '
                                'threshold are required separated by '
                                'space, e.g. --filterThreshold -1.5 5. Applied'
                                ' only for ICE!',
                           type=float,
                           nargs=2,
                           default=None)

    parserOpt.add_argument('--iterNum', '-n',
                           help='Number of iterations to compute.'
                           'only for ICE!'
                           ' (Default: %(default)s).',
                           type=int,
                           metavar='INT',
                           default=500)

    parserOpt.add_argument('--inflationCutoff',
                           help='Value corresponding to the maximum number of '
                           'times a bin can be scaled up during the iterative '
                           'correction. For example, an inflation cutoff of 3 '
                           'will filter out all bins that were expanded 3 '
                           'times or more during the iterative correctionself.'
                           'Only for ICE!',
                           type=float)

    parserOpt.add_argument('--transCutoff', '-transcut',
                           help='Clip high counts in the top -transcut trans '
                           'regions (i.e. between chromosomes). A usual value '
                           'is 0.05. Only for ICE! ',
                           type=float)

    parserOpt.add_argument('--sequencedCountCutoff',
                           help='Each bin receives a value indicating the '
                           'fraction that is covered by reads. A cutoff of '
                           '0.5 will discard all those bins that have less '
                           'than half of the bin covered. Only for ICE!',
                           default=None,
                           type=float)

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the '
                           'iterative correction. The order of the given '
                           'chromosomes will be then kept for the resulting '
                           'corrected matrix',
                           default=None,
                           nargs='+')

    parserOpt.add_argument('--skipDiagonal', '-s',
                           help='If set, diagonal counts are not included.'
                           ' Only for ICE!',
                           action='store_true')

    parserOpt.add_argument('--perchr',
                           help='Normalize each chromosome separately. This is'
                           ' useful for samples from cells with uneven number '
                           'of chromosomes and/or translocations.',
                           action='store_true')

    parserOpt.add_argument('--filteredBed',
                           help='Print bins filtered our by  --filterThreshold to this file')

    parserOpt.add_argument('--verbose',
                           help='Print processing status.',
                           action='store_true')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def iterative_correction(matrix, args):
    corrected_matrix, correction_factors = iterativeCorrection(matrix,
                                                               M=args.iterNum,
                                                               verbose=args.verbose)

    return corrected_matrix, correction_factors


def fill_gaps(hic_ma, failed_bins, fill_contiguous=False):
    """ try to fill-in the failed_bins the matrix by adding the
    average values of the neighboring rows and cols. The idea
    for the iterative correction is that is best to put
    something in contrast to not put anything

    hic_ma:  Hi-C matrix object
    failed_bins: list of bin ids
    fill_contiguous: If True, stretches of masked rows/cols are filled.
                     Otherwise, these cases are skipped

    """
    log.debug("starting fill gaps")
    mat_size = hic_ma.matrix.shape[0]
    fill_ma = hic_ma.matrix.copy().tolil()
    if fill_contiguous is True:
        discontinuous_failed = failed_bins
        consecutive_failed_idx = np.array([])
    else:
        # find stretches of consecutive failed regions
        consecutive_failed_idx = np.flatnonzero(np.diff(failed_bins) == 1)
        # the banned list of indices is equal to the actual list
        # and the list plus one, to identify consecutive failed regions.
        # for [1,2,5,10] the np.diff is [1,3,5]. The consecutive id list
        # is [0], for '1', in the original list, but we are missing the '2'
        # thats where the consecutive_failed_idx+1 comes.
        consecutive_failed_idx = np.unique(np.sort(
            np.concatenate([consecutive_failed_idx,
                            consecutive_failed_idx + 1])))
        # find the failed regions that are not consecutive
        discontinuous_failed = [x for idx, x in enumerate(failed_bins)
                                if idx not in consecutive_failed_idx]

    log.debug("Filling {} failed bins\n".format(
        len(discontinuous_failed)))

    """
    for missing_bin in discontinuous_failed:
        if 0 < missing_bin < mat_size - 1:
            for idx in range(1, mat_size - 2):
                if idx % 100 == 0:
                    log.info(".")
                # the new row value is the mean between the upper
                # and lower bins corresponding to the same diagonal
                fill_ma[missing_bin, idx :] = \
                    (hic_ma.matrix[missing_bin-1, idx-1] +
                     hic_ma.matrix[missing_bin+1, idx+1]) / 2

                # same for cols
                fill_ma[idx, missing_bin] = \
                    (hic_ma.matrix[idx-1, missing_bin-1] +
                     hic_ma.matrix[idx+1, missing_bin+1]) / 2

    """
    for missing_bin in discontinuous_failed:
        if 0 < missing_bin < mat_size - 1:
            # the new row value is the mean between the upper
            # and lower rows
            fill_ma[missing_bin, 1:mat_size - 1] = \
                (hic_ma.matrix[missing_bin - 1, :mat_size - 2] +
                 hic_ma.matrix[missing_bin + 1, 2:]) / 2

            # same for cols
            fill_ma[1:mat_size - 1, missing_bin] = \
                (hic_ma.matrix[:mat_size - 2, missing_bin - 1] +
                 hic_ma.matrix[2:, missing_bin + 1]) / 2

    # identify the intersection points of the failed regions because they
    # neighbors get wrong values
    for bin_a in discontinuous_failed:
        for bin_b in discontinuous_failed:
            if 0 < bin_a < mat_size and \
                    0 < bin_b < mat_size:
                # the fill value is the average over the
                # neighbors that do have a value

                fill_value = np.mean([
                    hic_ma.matrix[bin_a - 1, bin_b - 1],
                    hic_ma.matrix[bin_a - 1, bin_b + 1],
                    hic_ma.matrix[bin_a + 1, bin_b - 1],
                    hic_ma.matrix[bin_a + 1, bin_b + 1],
                ])

                fill_ma[bin_a - 1, bin_b] = fill_value
                fill_ma[bin_a + 1, bin_b] = fill_value
                fill_ma[bin_a, bin_b - 1] = fill_value
                fill_ma[bin_a, bin_b + 1] = fill_value

    # return the matrix and the bins that continue to be failed regions
    return fill_ma.tocsr(), np.sort(failed_bins[consecutive_failed_idx])


class MAD(object):

    def __init__(self, points):
        """
        Returns a boolean array with True if points are outliers and False
        otherwise.

        :param points: An array

        :returns: A numobservations-length boolean array.

        :references: Boris Iglewicz and David Hoaglin (1993), "Volume 16:
         How to Detect and Handle Outliers", The ASQC Basic References in
         Quality Control:
         Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
        """

        self.mad_b_value = 0.6745
        if len(points.shape) == 1:
            points = points[:, None]
        self.median = np.median(points[points > 0], axis=0)

        diff = np.sum((points - self.median), axis=-1)

        self.med_abs_deviation = np.median(np.abs(diff))
        self.modified_z_score = np.multiply(self.mad_b_value, np.divide(diff, self.med_abs_deviation))

    def get_motified_zscores(self):

        return self.modified_z_score

    def is_outlier(self, lower_threshold, upper_threshold):
        """
        Returns a boolean list of outliers

        :param lower_threshold: Lower median absolute deviation
        :param upper_threshold: upper median absolute deviation

        :return: boolean array
        """

        return (self.modified_z_score < lower_threshold) | \
               (self.modified_z_score > upper_threshold)

    def value_to_mad(self, value):
        """
        return the mad value for a given value
        based on the data
        """
        log.debug("self.median: {}".format(self.median))
        diff = value - self.median
        log.debug("diff: {}".format(diff))
        log.debug("self.med_abs_deviation: {}".format(self.med_abs_deviation))
        log.debug("self.mad_b_value: {}".format(self.mad_b_value))
        log.debug("all together: {}".format(
            self.mad_b_value * diff / self.med_abs_deviation))
        # workaround for 'Axis limits cannot be NaN or Inf' bug in version 2.1.1
        # prevent dividing by 0!!!
        if self.med_abs_deviation == 0.0:
            return self.mad_b_value * diff

        return self.mad_b_value * diff / self.med_abs_deviation

    def mad_to_value(self, mad):
        """
        return the numeric value for a given mad score
        based on the data

        z = b_v * (x - median) / mad
        z * mad / b_v = x - median
        (z * mad / b_v) + median = x
        """

        return (mad * self.med_abs_deviation / self.mad_b_value) + self.median


def plot_total_contact_dist(hic_ma, args):
    """
    Plots the distribution of number of contacts (excluding self contacts)
    Outliers with a high number are removed for the plot

    :param hic_ma: sparse matrix
    :return:
    """
    use('Agg')

    majorlocator = MultipleLocator(1)
    majorformatter = FormatStrFormatter('%d')
    minorlocator = MultipleLocator(0.2)

    def plot_histogram(row_sum_values, mad_values, ax1, title=None):

        if args.xMax:
            ax1.set_xlim(ax1.get_xlim()[0], args.xMax)
            row_sum_values = row_sum_values[row_sum_values < args.xMax]

        ax1.set_xlabel("total counts per bin")
        ax1.set_ylabel("frequency")
    #    ax1.xaxis.grid(True)
        ax1.patch.set_visible(False)
        dist, bin_s, __ = ax1.hist(row_sum_values, 100, color='green')

        # add second axis on top
        ax2 = ax1.twiny()
        ax2.set_xlabel("modified z-score")
        ax2.xaxis.set_major_locator(majorlocator)
        ax2.xaxis.set_major_formatter(majorformatter)
        ax2.xaxis.grid(True, which='minor')
        # for the minor ticks, use no labels; default NullFormatter
        ax2.xaxis.set_minor_locator(minorlocator)

        # update second axis values by mapping the min max
        # of the main axis to the translated values
        # into modified z score.

        # workaround for 'Axis limits cannot be NaN or Inf' bug in version 2.1.1
        log.debug("ax1.get_xlim(): {}".format(ax1.get_xlim()))
        log.debug("np.array(ax1.get_xlim()): {}".format(np.array(
                                                        ax1.get_xlim())))
        log.debug("mad_values.value_to_mad(np.array(ax1.get_xlim())): {}".format(
            mad_values.value_to_mad(np.array(ax1.get_xlim()))))

        ax2.set_xlim(mad_values.value_to_mad(np.array(ax1.get_xlim())))

        # get first local mininum value
        local_min = [x for x, y in enumerate(dist) if 1 <= x < len(
            dist) - 1 and dist[x - 1] > y < dist[x + 1]]

        if len(local_min) > 0:
            threshold = bin_s[local_min[0]]
        else:
            threshold = None

        if threshold:
            mad_threshold = mad_values.value_to_mad(threshold)
            ymin, ymax = ax2.get_ylim()
            ax2.vlines(mad_threshold, ymin, ymax)
            if title:
                log.info("{}: mad threshold {}".format(title, mad_threshold))
            else:
                log.info("mad threshold {}".format(mad_threshold))

    # replace nan by 0
    # hic_ma.matrix.data[np.isnan(hic_ma.matrix.data)] = 0
    hic_ma.matrix = convertNansToZeros(hic_ma.matrix)
    hic_ma.matrix = convertInfsToZeros(hic_ma.matrix)

    if args.perchr:
        chroms = hic_ma.getChrNames()
        if len(chroms) > 30:
            log.warning("The matrix contains {} chromosomes. It is not "
                        "practical to plot each. Try using --chromosomes to "
                        "select some chromosomes or plot a single histogram.")
        num_rows = int(np.ceil(float(len(chroms)) / 5))
        num_cols = min(len(chroms), 5)
        grids = gridspec.GridSpec(num_rows, num_cols)
        fig = plt.figure(figsize=(6 * num_cols, 5 * num_rows))
        ax = {}
        for plot_num, chrname in enumerate(chroms):
            log.info("Plotting chromosome {}".format(chrname))

            chr_range = hic_ma.getChrBinRange(chrname)
            chr_submatrix = hic_ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]

            row_sum = np.asarray(chr_submatrix.sum(axis=1)).flatten()
            row_sum = row_sum - chr_submatrix.diagonal()
            mad = MAD(row_sum)
            modified_z_score = mad.get_motified_zscores()

            # high remove outliers
            row_sum = row_sum[modified_z_score < 5]
            if len(row_sum) == 0:
                log.error('No data for chromosome {}. Continue with next chromosome.'.format(chrname))
            col = plot_num % num_cols
            row = plot_num // num_cols
            ax[chrname] = fig.add_subplot(grids[row, col])

            plot_histogram(row_sum, mad, ax[chrname], title=chrname)
            ax[chrname].set_title(chrname)
    else:
        fig = plt.figure()
        row_sum = np.asarray(hic_ma.matrix.sum(axis=1)).flatten()
        row_sum = row_sum - hic_ma.matrix.diagonal()
        mad = MAD(row_sum)
        modified_z_score = mad.get_motified_zscores()

        # high remove outliers
        row_sum = row_sum[modified_z_score < 5]
        if len(row_sum) == 0:
            log.error('No data available, exit.')
            exit(1)
        ax = fig.add_subplot(111)
        plot_histogram(row_sum, mad, ax)

    plt.tight_layout()
    plt.savefig(args.plotName)
    plt.close()


def filter_by_zscore(hic_ma, lower_threshold, upper_threshold, perchr=False):
    """
    The method defines thresholds per chromosome
    to avoid introducing bias due to different chromosome numbers

    """
    log.info("filtering by z-score")
    to_remove = []
    if perchr:
        for chrname in list(hic_ma.interval_trees):
            chr_range = hic_ma.getChrBinRange(chrname)
            chr_submatrix = hic_ma.matrix[chr_range[0]:chr_range[1],
                                          chr_range[0]:chr_range[1]]

            # replace nan values by zero
            chr_submatrix.data[np.isnan(chr_submatrix.data)] = 0
            row_sum = np.asarray(chr_submatrix.sum(axis=1)).flatten()
            # subtract from row sum, the diagonal
            # to account for interactions with other bins
            # and not only self interactions that are the dominant count
            row_sum = row_sum - chr_submatrix.diagonal()
            mad = MAD(row_sum)
            problematic = np.flatnonzero(
                mad.is_outlier(lower_threshold, upper_threshold))

            # because the problematic indices are specific for the given chromosome
            # they need to be updated to match the large matrix indices
            problematic += chr_range[0]

            if len(problematic) == 0:
                log.warn("Warning. No bins removed for chromosome {} using thresholds {} {}"
                         "\n".format(chrname, lower_threshold, upper_threshold))

            to_remove.extend(problematic)

    else:
        row_sum = np.asarray(hic_ma.matrix.sum(axis=1)).flatten()
        # subtract from row sum, the diagonal
        # to account for interactions with other bins
        # and not only self interactions that are the dominant count
        row_sum = row_sum - hic_ma.matrix.diagonal()
        mad = MAD(row_sum)
        to_remove = np.flatnonzero(mad.is_outlier(
            lower_threshold, upper_threshold))
    return sorted(to_remove)


def main(args=None):
    """Run the C++ implementation of hicCorrectMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicCorrectMatrix')(args)
