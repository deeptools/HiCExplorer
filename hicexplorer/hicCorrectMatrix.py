from __future__ import division
import warnings
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

import numpy as np
debug = 0

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        conflict_handler='resolve',
        description="""
Iterative correction for a Hi-C matrix (see Imakaev et al. 2012
Nature Methods for details). For the method to work correctly, bins with
zero reads assigned to them should be removed as they can not be corrected.
Also, bins with low number of reads should be removed,
otherwise, during the correction step, the counts associated with
those bins will be amplified (usually, zero and low coverage bins
tend contain repetitive regions).  Bins with extremely high number
of reads can also be removed from the correction as they may represent
copy number variations.

To aid in the identification of bins with low and high read coverage, the
histogram of the number of reads can be plotted together with the
Median Absolute Deviation (MAD).

It is recommended to run hicCorrectMatrix as follows:

    $ hicCorrectMatrix diagnostic_plot --matrix hic_matrix.h5 -o plot_file.png

Then, after revising the plot and deciding the threshold values:

    $ hicCorrectMatrix correct --matrix hic_matrix.h5 \r
    --filterThreshold <lower threshold> <upper threshold> -o corrected_matrix

For a more in-depth review of how to determine the threshold values, please visit:
http://hicexplorer.readthedocs.io/en/latest/content/example_usage.html#correction-of-hi-c-matrix
"""
    )

    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        title="Options",
        dest='command',
        metavar='',
        help="""To get detailed help on each of the options: \r

    $ hicCorrectMatrix diagnostic_plot -h \r

    $ hicCorrectMatrix correct -h
    """)
    plot_mode = subparsers.add_parser(
        'diagnostic_plot',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="""Plots a histogram of the coverage per bin together with the modified
z-score based on the median absolute deviation method
(see Boris Iglewicz and David Hoaglin 1993, Volume 16: How to Detect
and Handle Outliers The ASQC Basic References in Quality Control:
Statistical Techniques, Edward F. Mykytka, Ph.D., Editor).
        """,
        usage='%(prog)s '
              '--matrix hic_matrix.h5 '
              '-o file.png')
    plot_modeRequired = plot_mode.add_argument_group('Required arguments')
    plot_modeRequired.add_argument('--matrix', '-m',
                                   help='Name of the Hi-C matrix to correct in .h5 format.',
                                   required=True)

    plot_modeRequired.add_argument('--plotName', '-o',
                                   help='File name to save the diagnostic plot.',
                                   required=True)

    plot_modeOpt = plot_mode.add_argument_group('Optional arguments')
    plot_modeOpt.add_argument('--chromosomes',
                              help='List of chromosomes to be included in the iterative '
                              'correction. The order of the given chromosomes will be then '
                              'kept for the resulting corrected matrix.',
                              default=None,
                              nargs='+')

    plot_modeOpt.add_argument('--xMax',
                              help='Max value for the x-axis in counts per bin.',
                              default=None,
                              type=float)

    plot_modeOpt.add_argument(
        '--perchr',
        help='Compute histogram per chromosome. For samples from cells with uneven number '
        'of chromosomes and/or translocations it is advisable to check the histograms '
        'per chromosome to find the most conservative `filterThreshold`.',
        action='store_true')

    plot_modeOpt.add_argument('--verbose',
                              help='Print processing status.',
                              action='store_true')

    subparsers.add_parser(
        'correct',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[correct_subparser()],
        help="""Run the iterative correction.""",
        usage='%(prog)s '
              '--matrix hic_matrix.h5 '
              '--filterThreshold -1.2 5 '
              '-out corrected_matrix.h5 \n')

    return parser


def correct_subparser():
    # define the arguments
    parser = argparse.ArgumentParser(add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Name of the Hi-C matrix to correct in .h5 format.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. The '
                                'output is a .h5 file.',
                                required=True)

    parserRequired.add_argument('--filterThreshold', '-t',
                                help='Removes bins of low or large coverage. '
                                'Usually these bins do not contain valid Hi-C data or represent '
                                'regions that accumulate reads and thus must be discarded. '
                                'Use hicCorrectMatrix diagnostic_plot '
                                'to identify the modified z-value thresholds. A lower and upper '
                                'threshold are required separated by space, e.g. --filterThreshold '
                                '-1.5 5',
                                type=float,
                                nargs=2,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--iterNum', '-n',
                           help='Number of iterations to compute.',
                           type=int,
                           metavar='INT',
                           default=500)

    parserOpt.add_argument('--inflationCutoff',
                           help='Value corresponding to the maximum number of times a bin '
                           'can be scaled up during the iterative correction. For example, '
                           'an inflation cutoff of 3 will filter out all bins that were '
                           'expanded 3 times or more during the iterative correction.',
                           type=float)

    parserOpt.add_argument('--transCutoff', '-transcut',
                           help='Clip high counts in the top -transcut trans '
                           'regions (i.e. between chromosomes). A usual value '
                           'is 0.05 ',
                           type=float)

    parserOpt.add_argument('--sequencedCountCutoff',
                           help='Each bin receives a value indicating the '
                           'fraction that is covered by reads. A cutoff of '
                           '0.5 will discard all those bins that have less '
                           'than half of the bin covered.',
                           default=None,
                           type=float)

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the iterative '
                           'correction. The order of the given chromosomes will be then '
                           'kept for the resulting corrected matrix',
                           default=None,
                           nargs='+')

    parserOpt.add_argument('--skipDiagonal', '-s',
                           help='If set, diagonal counts are not included',
                           action='store_true')

    parserOpt.add_argument('--perchr',
                           help='Normalize each chromosome separately. This is useful for '
                           'samples from cells with uneven number of chromosomes and/or translocations.',
                           action='store_true')

    parserOpt.add_argument('--verbose',
                           help='Print processing status',
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

        :references: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
            Handle Outliers", The ASQC Basic References in Quality Control:
            Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
        """

        self.mad_b_value = 0.6745
        if len(points.shape) == 1:
            points = points[:, None]
        self.median = np.median(points[points > 0], axis=0)

        diff = np.sum((points - self.median), axis=-1)

        self.med_abs_deviation = np.median(np.abs(diff))
        self.modified_z_score = self.mad_b_value * diff / self.med_abs_deviation

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
        log.debug("all together: {}".format(self.mad_b_value * diff / self.med_abs_deviation))
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
    from matplotlib import use
    use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
        log.debug("np.array(ax1.get_xlim()): {}".format(np.array(ax1.get_xlim())))
        log.debug("mad_values.value_to_mad(np.array(ax1.get_xlim())): {}".format(mad_values.value_to_mad(np.array(ax1.get_xlim()))))

        ax2.set_xlim(mad_values.value_to_mad(np.array(ax1.get_xlim())))

        # get first local mininum value
        local_min = [x for x, y in enumerate(dist) if 1 <= x < len(dist) - 1 and
                     dist[x - 1] > y < dist[x + 1]]

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
            log.warning("The matrix contains {} chromosomes. It is not practical to plot "
                        "each. Try using --chromosomes to select some chromosomes or "
                        "plot a single histogram.")
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
            problematic = np.flatnonzero(mad.is_outlier(lower_threshold, upper_threshold))

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
        to_remove = np.flatnonzero(mad.is_outlier(lower_threshold, upper_threshold))

    return sorted(to_remove)


def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.verbose:
        log.setLevel(logging.INFO)

    # args.chromosomes
    if check_cooler(args.matrix) and args.chromosomes is not None and len(args.chromosomes) == 1:
        ma = hm.hiCMatrix(args.matrix, pChrnameList=toString(args.chromosomes))
    else:
        ma = hm.hiCMatrix(args.matrix)

        if args.chromosomes:
            ma.reorderChromosomes(toString(args.chromosomes))

    # mask all zero value bins
    row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()
    log.info("Removing {} zero value bins".format(sum(row_sum == 0)))
    ma.maskBins(np.flatnonzero(row_sum == 0))
    matrix_shape = ma.matrix.shape
    ma.matrix = convertNansToZeros(ma.matrix)
    ma.matrix = convertInfsToZeros(ma.matrix)

    if 'plotName' in args:
        plot_total_contact_dist(ma, args)
        log.info("Saving diagnostic plot {}\n".format(args.plotName))
        return

    log.info("matrix contains {} data points. Sparsity {:.3f}.".format(
        len(ma.matrix.data),
        float(len(ma.matrix.data)) / (ma.matrix.shape[0] ** 2)))

    if args.skipDiagonal:
        ma.diagflat(value=0)

    outlier_regions = filter_by_zscore(ma, args.filterThreshold[0], args.filterThreshold[1], perchr=args.perchr)
    # compute and print some statistics
    pct_outlier = 100 * float(len(outlier_regions)) / ma.matrix.shape[0]
    ma.printchrtoremove(outlier_regions, label="Bins that are MAD outliers ({:.2f}%) "
                                               "out of".format(pct_outlier, ma.matrix.shape[0]),
                        restore_masked_bins=False)

    assert matrix_shape == ma.matrix.shape
    # mask filtered regions
    ma.maskBins(outlier_regions)
    total_filtered_out = set(outlier_regions)

    if args.sequencedCountCutoff and 0 < args.sequencedCountCutoff < 1:
        chrom, _, _, coverage = zip(*ma.cut_intervals)

        assert type(coverage[0]) == np.float64

        failed_bins = np.flatnonzero(
            np.array(coverage) < args.sequencedCountCutoff)

        ma.printchrtoremove(failed_bins, label="Bins with low coverage", restore_masked_bins=False)
        ma.maskBins(failed_bins)
        total_filtered_out = set(failed_bins)
        """
        ma.matrix, to_remove = fill_gaps(ma, failed_bins)
        log.warning("From {} failed bins, {} could "
                         "not be filled\n".format(len(failed_bins),
                                                  len(to_remove)))
        ma.maskBins(to_remove)
        """

    if args.transCutoff and 0 < args.transCutoff < 100:
        cutoff = float(args.transCutoff) / 100
        # a usual cutoff is 0.05
        ma.truncTrans(high=cutoff)

    pre_row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()
    correction_factors = []
    if args.perchr:
        corrected_matrix = lil_matrix(ma.matrix.shape)
        # normalize each chromosome independently
        for chrname in list(ma.interval_trees):
            chr_range = ma.getChrBinRange(chrname)
            chr_submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            _matrix, _corr_factors = iterative_correction(chr_submatrix, args)
            corrected_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = _matrix
            correction_factors.append(_corr_factors)
        correction_factors = np.concatenate(correction_factors)

    else:
        corrected_matrix, correction_factors = iterative_correction(ma.matrix, args)

    ma.setMatrixValues(corrected_matrix)
    ma.setCorrectionFactors(correction_factors)
    log.info("Correction factors {}".format(correction_factors[:10]))
    if args.inflationCutoff and args.inflationCutoff > 0:
        after_row_sum = np.asarray(corrected_matrix.sum(axis=1)).flatten()
        # identify rows that were expanded more than args.inflationCutoff times
        to_remove = np.flatnonzero(after_row_sum / pre_row_sum >= args.inflationCutoff)
        ma.printchrtoremove(to_remove,
                            label="inflated >={} "
                            "regions".format(args.inflationCutoff), restore_masked_bins=False)
        total_filtered_out = total_filtered_out.union(to_remove)

        ma.maskBins(to_remove)

    ma.printchrtoremove(sorted(list(total_filtered_out)),
                        label="Total regions to be removed", restore_masked_bins=False)

    ma.save(args.outFileName, pApplyCorrection=False)
