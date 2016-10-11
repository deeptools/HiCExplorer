import sys, argparse
from scipy.sparse import lil_matrix
import logging

from hicexplorer.iterativeCorrection import iterativeCorrection
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__

import numpy as np
debug = 0


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""
Iterative correction for a hic matrix (see Imakaev et al. 2012
Nature Methods for details). For the method to work correctly, bins
with low or too high coverage need to be filtered. For this, it is
recommended to first run some diagnostic plots to determine the
modified z-score cut off.

It is recommended to run hicCorrectMatrix as follows:

    $ hicCorrectMatrix diagnostic_plot --matrix hic_matrix -o plot_file.png

Then, after revising the plot and deciding the threshold values:

    $ hicCorrectMatrix correct --matrix hic_matrix \\
         --filterThreshold <lower threshold> <upper threshold> -o corrected_matrix


To get detailed help on each of the options:

    $ hicCorrectMatrix diagnostic_plot -h
    $ hicCorrectMatrix correct -h


"""
    )

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        metavar='')

    correct_mode = subparsers.add_parser(
        'correct',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[correct_subparser()],
        help="Run the iterative correction.",
        usage='%(prog)s correct '
              '--matrix hic_matrix.npz '
              '--filterThreshold -1.2 5'
              '-out corrected_matrix.npz \n')

    plot_mode = subparsers.add_parser(
        'diagnostic_plot',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Plots a histogram of the coverage per bin together with the modified z-score "
             "based on the median absolute deviation method (see Boris Iglewicz and David Hoaglin 1993, "
             "Volume 16: How to Detect and Handle Outliers The ASQC Basic References in Quality Control:  "
             "Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. ",
        usage='%(prog)s '
              '--matrix hic_matrix.npz '
              '-o file.png')
    plot_mode.add_argument('--matrix', '-m',
                           help='Hi-C matrix.',
                           required=True)

    plot_mode.add_argument('--plotName', '-o',
                           help='File name to save the diagnostic plot.',
                           required=True)

    plot_mode.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the iterative '
                           'correction. The order of the given chromosomes will be then '
                           'kept for the resulting corrected matrix',
                           default=None,
                           nargs='+')

    plot_mode.add_argument('--xMax',
                        help='Max value for the x-axis in counts per bin',
                        default=None,
                        type=float)


    merge_mode = subparsers.add_parser(
        'merge_failed',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Merges together failed bins to rescue some of the information instead of discarding it. This option "
             "is mostly useful with processing small restriction fragment size bins.",
        usage='%(prog)s '
              '--matrix hic_matrix.npz '
              '--outMatrixFile hic_matrix_merged_failed.npz '
              '-o file.png')
    merge_mode.add_argument('--matrix', '-m',
                            help='Hi-C matrix.',
                            required=True)

    merge_mode.add_argument('--outMatrixFile',
                            help='Name to save the resulting matrix.',
                            required=True)

    merge_mode.add_argument('--plotName', '-o',
                            help='File name to save the diagnostic plot.',
                            required=True)

    merge_mode.add_argument('--filterThreshold', '-t',
                        help='Bins of low coverage or large coverage need to be removed. '
                             'Usually they do not contain valid Hi-C data of represent '
                             'regions that accumulate reads. Use hicCorrectMatrix diagnostic_plot '
                             'to identify the modified z-value thresholds. A lower and upper '
                             'threshold are required separated by space. Eg. --filterThreshold '
                             '-1.5 5',
                        type=float,
                        nargs=2,
                        required=True)

    merge_mode.add_argument('--chromosomes',
                            help='List of chromosomes to be included in the iterative '
                            'correction. The order of the given chromosomes will be then '
                            'kept for the resulting corrected matrix',
                            default=None,
                            nargs='+')

    merge_mode.add_argument('--xMax',
                            help='Max value for the X field in counts per bin',
                            default=None,
                            type=float)
    return parser


def correct_subparser():
    # define the arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--matrix', '-m',
                        help='Hi-C matrix.',
                        required=True)

    parser.add_argument('--iterNum', '-n',
                        help='number of iterations',
                        type=int,
                        metavar='INT',
                        default=500)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. The '
                             'output is a .npz file.',
                        required=True)

    parser.add_argument('--filterThreshold', '-t',
                        help='Bins of low coverage or large coverage need to be removed. '
                             'Usually they do not contain valid Hi-C data of represent '
                             'regions that accumulate reads. Use hicCorrectMatrix diagnostic_plot '
                             'to identify the modified z-value thresholds. A lower and upper '
                             'threshold are required separated by space. Eg. --filterThreshold '
                             '-1.5 5',
                        type=float,
                        nargs=2,
                        required=True)

    parser.add_argument('--inflationCutoff',
                        help='Value corresponding to the maximum number of times a bin '
                        'can be scaled up during the iterative correction. For example, '
                        'a inflation Cutoff of 3 will filter out all bins that were '
                        'expanded 3 times or more during the iterative correction.',
                        type=float)

    parser.add_argument('--transCutoff', '-transcut',
                        help='Clip high counts in the top -transcut trans '
                        'regions (i.e. between chromosomes). A usual value '
                        'is 0.05 ',
                        type=float)

    parser.add_argument('--sequencedCountCutoff',
                        help='Each bin receives a value indicating the '
                        'fraction that is covered by reads. A cutoff of '
                        '0.5 will discard all those bins that have less '
                        'than half of the bin covered.',
                        default=None,
                        type=float)

    parser.add_argument('--chromosomes',
                        help='List of chromosomes to be included in the iterative '
                        'correction. The order of the given chromosomes will be then '
                        'kept for the resulting corrected matrix',
                        default=None,
                        nargs='+')

    parser.add_argument('--skipDiagonal', '-s',
                        help='If set, diagonal counts are not included',
                        action='store_true')

    parser.add_argument('--perchr',
                        help='Normalize each chromosome separately',
                        action='store_true')

    parser.add_argument('--verbose',
                        help='Print processing status',
                        action='store_true')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def iterative_correction(matrix, args):
    corrected_matrix, correction_factors = iterativeCorrection(matrix,
                                                               M=args.iterNum,
                                                               verbose=args.verbose)

    return corrected_matrix, correction_factors


def merge_failed_bins(hic_matrix, failed_bins):
    """
    Merges the failed bins instead of removing them
    :param hic_matrix: hicMatrix object
    :param failed_bins: list of failed bins
    :return: hicMatrix object
    """

    import hicexplorer.reduceMatrix
    hic_matrix.restoreMaskedBins()
    # get the bins to merge
    ref_name_list, start_list, end_list, coverage_list = zip(*hic_matrix.cut_intervals)
    new_bins = []
    bins_to_merge = []
    coverage_list = np.array(coverage_list)
    def get_merged_bin(bin_list):
        if len(bin_list) < 2:
            print idx, bin_list, consecutive, ref_name_list[idx - 1], ref_name_list[idx]
            import ipdb;ipdb.set_trace()
        assert len(bin_list) > 1, "Error, bin_list length has less than 2 elements."
        coverage = np.mean(coverage_list[bin_list])
        return ref_name_list[bin_list[0]], start_list[bin_list[0]], end_list[bin_list[-1]], coverage

    # if consecutive failed bins are found
    # merge them together
    # otherwise, merge them to the smallest neighboring bin
    consecutive = []
    for idx in range(len(ref_name_list)):
        if idx in failed_bins:
            # chromosome name change
            if idx > 1 and ref_name_list[idx - 1] != ref_name_list[idx]:
                if consecutive:
                    if len(consecutive) == 1:
                        new_bins.append(hic_matrix.cut_intervals[consecutive[0]])
                        bins_to_merge.append([consecutive[0]])
                    else:
                        new_bins.append(get_merged_bin(consecutive))
                        bins_to_merge.append(consecutive)
                    consecutive = []
                if idx + 1 in failed_bins:
                    consecutive.append(idx)
                    continue

            elif idx + 1 in failed_bins:
                consecutive.append(idx)
                continue
            elif len(consecutive):
                consecutive.append(idx)
                new_bins.append(get_merged_bin(consecutive))
                bins_to_merge.append(consecutive)
                consecutive = []
                continue

            if idx == 0 or ref_name_list[idx - 1] != ref_name_list[idx]:
                # can only merge to the right bin
                new_bins.append(get_merged_bin([idx, idx + 1]))
                bins_to_merge.append([idx, idx + 1])
                continue
            elif idx + 1 == len(ref_name_list):
                # can only merge to the left bin, but since this should have been
                # already added, then it is updated
                new_bins[-1] = get_merged_bin([idx -1, idx])
                bins_to_merge[-1] = [idx -1, idx]
                continue

            # merge to the shorter neighboring bin
            prev_bin_len = end_list[idx -1] - start_list[idx -1]
            next_bin_len = end_list[idx +1] - start_list[idx -+2]
            if prev_bin_len < next_bin_len:
                new_bins[-1] = get_merged_bin([idx -1, idx])
                bins_to_merge[-1] = [idx -1, idx]
            else:
                new_bins.append(get_merged_bin([idx, idx + 1]))
                bins_to_merge.append([idx, idx + 1])

        else:
            # skip if the bin was already added in the previous loop
            if idx == 0 or idx not in bins_to_merge[-1]:
                bins_to_merge.append([idx,])
                new_bins.append(hic_matrix.cut_intervals[idx])

    diff = np.diff(np.concatenate(bins_to_merge))
    if len(np.flatnonzero(diff > 1)) != 0:
        import ipdb;ipdb.set_trace()
    assert len(np.flatnonzero(diff > 1)) == 0, "Some indexes are missing"
    hic_matrix.update_matrix(hicexplorer.reduceMatrix.reduce_matrix(hic_matrix.matrix, bins_to_merge, diagonal=True),
                             new_bins)
    return hic_matrix


def fill_gaps(hic_ma, failed_bins, fill_contiguous=False):
    """ try to fill-in the failed_bins the matrix by adding the
    average values of the neighboring rows and cols. The idea
    for the iterative correction is that is best to put
    something in contrast to not put anything

    hic_ma: hic matrix object
    failed_bins: list of bin ids
    fill_contiguous: If True, stretches of masked rows/cols are filled.
                     Otherwise, these cases are skipped

    """
    logging.info("starting fill gaps")
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
                                consecutive_failed_idx+1])))
        # find the failed regions that are not consecutive
        discontinuous_failed = [x for idx, x in enumerate(failed_bins)
                                if idx not in consecutive_failed_idx]

    sys.stderr.write("Filling {} failed bins\n".format(
            len(discontinuous_failed)))

    """
    for missing_bin in discontinuous_failed:
        if 0 < missing_bin < mat_size - 1:
            for idx in range(1, mat_size - 2):
                if idx % 100 == 0:
                    sys.stderr.write(".")
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
            fill_ma[missing_bin, 1:mat_size-1] = \
                (hic_ma.matrix[missing_bin - 1, :mat_size-2] +
                 hic_ma.matrix[missing_bin + 1, 2:]) / 2

            # same for cols
            fill_ma[1:mat_size-1, missing_bin] = \
                (hic_ma.matrix[:mat_size-2, missing_bin - 1] +
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
                        hic_ma.matrix[bin_a-1, bin_b-1],
                        hic_ma.matrix[bin_a-1, bin_b+1],
                        hic_ma.matrix[bin_a+1, bin_b-1],
                        hic_ma.matrix[bin_a+1, bin_b+1],
                        ])

                fill_ma[bin_a-1, bin_b] = fill_value
                fill_ma[bin_a+1, bin_b] = fill_value
                fill_ma[bin_a, bin_b-1] = fill_value
                fill_ma[bin_a, bin_b+1] = fill_value

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

        diff = value - self.median
        return self.mad_b_value * diff / self.med_abs_deviation


def plot_total_contact_dist(hic_ma, args):
    """
    Plots the distribution of number of contacts (excluding self contacts)
    Outliers with a high number are removed for the plot

    :param hic_ma: sparse matrix
    :return:
    """

    # replace nan by 0
    hic_ma.data[np.isnan(hic_ma.data)] = 0
    row_sum = np.asarray(hic_ma.sum(axis=1)).flatten()
    row_sum = row_sum - hic_ma.diagonal()
    mad = MAD(row_sum)
    modified_z_score = mad.get_motified_zscores()

    # high remove outliers
    row_sum = row_sum[modified_z_score < 5]

    from matplotlib import use
    use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    majorLocator = MultipleLocator(1)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(0.2)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    if args.xMax:
        ax1.set_xlim(ax1.get_xlim()[0], args.xMax)
        row_sum = row_sum[row_sum <  args.xMax]


    ax1.set_xlabel("total counts per bin")
    ax1.set_ylabel("frequency")
#    ax1.xaxis.grid(True)
    ax1.patch.set_visible(False)
    dist, bin_s, __ = ax1.hist(row_sum, 100, color='green')

    # add second axis on top
    ax2 = ax1.twiny()
    ax2.set_xlabel("modified z-score")
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.grid(True, which='minor')
    # for the minor ticks, use no labels; default NullFormatter
    ax2.xaxis.set_minor_locator(minorLocator)

    # update second axis values by mapping the min max
    # of the main axis to the translated values
    # into modified z score.
    ax2.set_xlim(mad.value_to_mad(np.array(ax1.get_xlim())))

    # get first local miminum value
    local_min = [x for x, y in enumerate(dist) if 1 <= x < len(dist) - 1 and
                 dist[x-1] > y < dist[x+1]]

    if len(local_min) > 0:
        threshold = bin_s[local_min[0]]
    else:
        threshold = None

    if threshold:
        mad_threshold = mad.value_to_mad(threshold)
        ymin, ymax = ax2.get_ylim()
        ax2.vlines(mad_threshold, ymin, ymax)
        print "mad threshold {}".format(mad_threshold)

    plt.savefig(args.plotName)
    plt.close()


def filter_by_zscore(hic_ma, lower_threshold, upper_threshold):
    """
    The method defines thresholds per chromosome
    to avoid introducing bias due to different chromosome numbers

    """
    to_remove = []
    for chrname in hic_ma.interval_trees.keys():
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
            sys.stderr.write("Warning. No bins removed for chromosome {} "
                             "using thresholds {} {}"
                             "\n".format(chrname, lower_threshold, upper_threshold))

        to_remove.extend(problematic)

    return sorted(to_remove)


def main():
    args = parse_arguments().parse_args()
    ma = hm.hiCMatrix(args.matrix)

    if args.chromosomes:
        ma.reorderChromosomes(args.chromosomes)

    if 'outMatrixFile' in args:
        # get below threshold outliers by using an extremely high upper threshold
        outlier_regions = filter_by_zscore(ma, args.filterThreshold[0], 1e6)
        print len(outlier_regions), ma.matrix.shape
        # compute and print some statistics
        pct_outlier = 100 * float(len(outlier_regions)) / ma.matrix.shape[0]
        ma.printchrtoremove(outlier_regions, label="Bins that are MAD outliers ({:.2f}%)".format(pct_outlier))
        # try to recover some of the outliers by merging them
        ma = merge_failed_bins(ma, outlier_regions)
        ma.save(args.outMatrixFile)
        plot_total_contact_dist(ma.matrix, args)
        sys.stderr.write("Saving diagnostic plot {}\n".format(args.plotName))
        exit()

    elif 'plotName' in args:
        plot_total_contact_dist(ma.matrix, args)
        sys.stderr.write("Saving diagnostic plot {}\n".format(args.plotName))
        exit()

    if args.verbose:
        print "matrix contains {} data points. Sparsity {:.3f}.".format(
            len(ma.matrix.data),
            float(len(ma.matrix.data))/(ma.matrix.shape[0]**2))

    if args.skipDiagonal:
        ma.diagflat(value=0)

    outlier_regions = filter_by_zscore(ma, args.filterThreshold[0], args.filterThreshold[1])
    # compute and print some statistics
    pct_outlier = 100 * float(len(outlier_regions)) / ma.matrix.shape[0]
    ma.printchrtoremove(outlier_regions, label="Bins that are MAD outliers after merge ({:.2f}%) "
                                               "out of".format(pct_outlier, ma.matrix.shape[0]))

    # mask filtered regions
    ma.maskBins(outlier_regions)
    total_filtered_out = set(outlier_regions)

    if args.sequencedCountCutoff and 0 < args.sequencedCountCutoff < 1:
        chrom, _, _, coverage = zip(*ma.cut_intervals)

        assert type(coverage[0]) == np.float64

        failed_bins = np.flatnonzero(
            np.array(coverage) < args.sequencedCountCutoff)

        ma.printchrtoremove(failed_bins, label="Bins with low coverage")
        ma.maskBins(failed_bins)
        total_filtered_out = set(failed_bins)
        """
        ma.matrix, to_remove = fill_gaps(ma, failed_bins)
        sys.stderr.write("From {} failed bins, {} could "
                         "not be filled\n".format(len(failed_bins),
                                                  len(to_remove)))
        ma.maskBins(to_remove)
        """



    if args.transCutoff and 0 < args.transCutoff < 100:
        cutoff = float(args.transCutoff)/100
        # a usual cutoff is 0.05
        ma.truncTrans(high=cutoff)

    pre_row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()
    correction_factors = []
    if args.perchr:
        corrected_matrix = lil_matrix(ma.matrix.shape)
        # normalize each chromosome independently
        for chrname in ma.interval_trees.keys():
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
    if args.inflationCutoff and args.inflationCutoff > 0:
        after_row_sum = np.asarray(corrected_matrix.sum(axis=1)).flatten()
        # identify rows that were expanded more than args.inflationCutoff times
        to_remove = np.flatnonzero(after_row_sum / pre_row_sum >= args.inflationCutoff)
        ma.printchrtoremove(to_remove,
                            label="inflated >={} "
                            "regions".format(args.inflationCutoff))
        total_filtered_out = total_filtered_out.union(to_remove)

        ma.maskBins(to_remove)

    ma.printchrtoremove(sorted(list(total_filtered_out)),
                        label="Total regions to be removed")

    ma.save(args.outFileName)
