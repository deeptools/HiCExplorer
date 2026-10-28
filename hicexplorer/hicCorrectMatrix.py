"""hicCorrectMatrix: the tool is the C++ executable hicCorrectMatrix.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names keep their 3.7.6 signatures and behaviour
and call hicexplorer._hicx.hicCorrectMatrix
(cpp/python/src/api/hicCorrectMatrix.cpp); plot_total_contact_dist draws with
matplotlib from data computed there.
"""
import warnings
import sys  # noqa: F401
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse  # noqa: F401
from scipy.sparse import lil_matrix  # noqa: F401

from hicexplorer.iterativeCorrection import iterativeCorrection  # noqa: F401 (iterative_correction calls it)
from hicmatrix import HiCMatrix as hm  # noqa: F401
from hicexplorer._version import __version__  # noqa: F401
from hicexplorer.utilities import toString  # noqa: F401
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros  # noqa: F401 (the plot data use them)
from hicexplorer.utilities import check_cooler  # noqa: F401


# packages needed for plotting:
from matplotlib import use
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib  # noqa: F401
import numpy as np
debug = 0

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicCorrectMatrix as _impl


def parse_arguments(args=None):
    """The argparse parser of hicCorrectMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCorrectMatrix')


def correct_subparser():
    return _impl.correct_subparser()


def iterative_correction(matrix, args):
    return _impl.iterative_correction(matrix, args)


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
    return _impl.fill_gaps(hic_ma, failed_bins, fill_contiguous)


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
        _impl.MAD___init__(self, points)

    def get_motified_zscores(self):
        return _impl.MAD_get_motified_zscores(self)

    def is_outlier(self, lower_threshold, upper_threshold):
        """
        Returns a boolean list of outliers

        :param lower_threshold: Lower median absolute deviation
        :param upper_threshold: upper median absolute deviation

        :return: boolean array
        """
        return _impl.MAD_is_outlier(self, lower_threshold, upper_threshold)

    def value_to_mad(self, value):
        """
        return the mad value for a given value
        based on the data
        """
        return _impl.MAD_value_to_mad(self, value)

    def mad_to_value(self, mad):
        """
        return the numeric value for a given mad score
        based on the data

        z = b_v * (x - median) / mad
        z * mad / b_v = x - median
        (z * mad / b_v) + median = x
        """
        return _impl.MAD_mad_to_value(self, mad)


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
            row_sum_values = _impl.values_below(row_sum_values, args.xMax)

        ax1.set_xlabel("total counts per bin")
        ax1.set_ylabel("frequency")
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
        threshold = _impl.histogram_threshold(row_sum_values)

        if threshold:
            mad_threshold = mad_values.value_to_mad(threshold)
            ymin, ymax = ax2.get_ylim()
            ax2.vlines(mad_threshold, ymin, ymax)
            if title:
                log.info("{}: mad threshold {}".format(title, mad_threshold))
            else:
                log.info("mad threshold {}".format(mad_threshold))

    # replace nan and inf by 0
    _impl.prepare_total_contact_dist(hic_ma)

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
            # the row sums without the self contacts, outliers removed
            row_sum, mad = _impl.chromosome_contact_dist(hic_ma, chrname)
            col = plot_num % num_cols
            row = plot_num // num_cols
            ax[chrname] = fig.add_subplot(grids[row, col])

            plot_histogram(row_sum, mad, ax[chrname], title=chrname)
            ax[chrname].set_title(chrname)
    else:
        fig = plt.figure()
        row_sum, mad = _impl.genome_contact_dist(hic_ma)
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
    return _impl.filter_by_zscore(hic_ma, lower_threshold, upper_threshold, perchr)


def main(args=None):
    """Run the C++ implementation of hicCorrectMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicCorrectMatrix')(args)
