#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
import sys
import logging
import argparse
import json
from collections import OrderedDict
from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import enlarge_bins
from scipy import sparse
import numpy as np
import multiprocessing
from hicexplorer._version import __version__
from hicexplorer.utilities import toString, toBytes, check_chrom_str_bytes

# python 2 / 3 compatibility
from past.builtins import zip
from six import iteritems
from builtins import range
from past.builtins import map

log = logging.getLogger(__name__)

# this is a holder vor the
hic_ma = None


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""
Uses a measure called TAD-separation score to identify the degree of separation between
the left and right regions at each Hi-C matrix bin. This is done for a
running window of different sizes. Then, TADs are called as those
positions having a local TAD-separation score minimum. The TAD-separation score is
measured using the z-score of the Hi-C matrix and is defined as the mean zscore of all
the matrix contacts between the left and right regions (diamond).

To find the TADs, the program  needs to compute first the
TAD scores at different window sizes. Then, the results of that computation
are used to call the TADs. This is convenient to test different filtering criteria quickly
as the demanding step is the computation of TAD-separation scores.

 A simple example usage is:

$ hicFindTads -m hic_matrix.h5 --outPrefix TADs --correctForMultipleTesting fdr

The bedgraph file produced by this tool can be used to plot the so-called insulation score
along the genome or at specific regions. This score is much more reliable across samples
than the number of TADs or the TADs width that can vary depending on the sequencing depth because of the lack
of information at certain bins, and depending on the parameters used with this tool.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Corrected Hi-C matrix to use for the computations',
                                required=True)

    parserRequired.add_argument('--outPrefix',
                                help='File prefix to save the resulting files: 1. <prefix>_tad_separation.bm '
                                'The format of the output file is chrom start end TAD-sep1 TAD-sep2 TAD-sep3 .. etc. '
                                'We call this format a bedgraph matrix and can be plotted using '
                                '`hicPlotTADs`. Each of the TAD-separation scores in the file corresponds to '
                                'a different window length starting from --minDepth to --maxDepth. '
                                '2. <prefix>_zscore_matrix.h5, the zscore matrix used for the computation of '
                                'the TAD-separation score.  3. < prefix > _boundaries.bed, which'
                                'contains the positions of boundaries. The genomic coordinates in this file '
                                'correspond to the resolution used. Thus, for Hi-C bins of '
                                '10.000bp the boundary position is 10.000bp long. For restriction fragment '
                                'matrices the boundary position varies depending on the fragment length '
                                'at the boundary. 4. <prefix>_domains.bed '
                                'contains the TADs positions. This is a non-overlapping set of genomic '
                                'positions. 5. <prefix>_boundaries.gff Similar to the boundaries bed file '
                                'but with extra information (pvalue, delta). 6. <prefix>_score.bedgraph file '
                                'contains the TAD-separation score '
                                'measured at each Hi-C bin coordinate. Is useful to visualize in a genome '
                                'browser. The delta and pvalue settings are saved as part of the name.',
                                required=True)

    parserRequired.add_argument('--correctForMultipleTesting',
                                help='Select the bonferroni or false discovery rate for a multiple comparison. Bonferroni '
                                'controlls the familywise error rate (FWER) and needs a p-value. The false discovery rate '
                                '(FDR) controls the likelyhood of type I errors and needs a q-value. As a third option '
                                'it is possible to not use a multiple comparison method at all.',
                                type=str,
                                default="fdr",
                                choices=['fdr', 'bonferroni', 'None'],
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--minDepth',
                           help='Minimum window length (in bp) to be considered to the left and to the right '
                           'of each Hi-C bin. This number should be at least 3 times '
                           'as large as the bin size of the Hi-C matrix.',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--maxDepth',
                           help='Maximum window length to be considered to the left and to the right '
                           'of the cut point in bp. This number should around 6-10 times '
                           'as large as the bin size of the Hi-C matrix.',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--step',
                           help='Step size when moving from --minDepth to --maxDepth. Note, the step size'
                           'grows exponentially as '
                           '`maxDeph + (step * int(x)**1.5) for x in [0, 1, ...]` until  it '
                           'reaches `maxDepth`. For example, selecting  step=10,000, minDepth=20,000 '
                           'and maxDepth=150,000 will compute TAD-scores for window sizes: '
                           '20,000, 30,000, 40,000, 70,000 and 100,000',
                           metavar='INT bp',
                           type=int)

    parserOpt.add_argument('--TAD_sep_score_prefix',
                           help='Sometimes it is useful to change some of the parameters without recomputing the '
                           'z-score matrix and the TAD-separation score. For this case, the prefix containing the '
                           'TAD separation score and the z-score matrix can be given. If this option is given, '
                           'new boundaries will be computed but the values of --minDepth, --maxDepth and --step will '
                           'not be used.',
                           required=False)

    parserOpt.add_argument('--thresholdComparisons',
                           help='P-value threshold for the bonferroni correction / q-value for FDR. '
                           'The probability of a local minima to be a boundary '
                           'is estimated by comparing the distribution (Wilcoxon ranksum) of '
                           'the  zscores between the left and right '
                           'regions (diamond) at the local minimum with the matrix zscores for a '
                           'diamond at --minDepth to the left and a diamond --minDepth to the right. '
                           'If --correctForMultipleTesting is \'None\' the threshold is applied on the '
                           'raw p-values without any multiple testing correction. Set it to \'1\' if no threshold should be used.',
                           type=float,
                           default=0.01)

    parserOpt.add_argument('--delta',
                           help='Minimum threshold of the difference between the TAD-separation score of a '
                           'putative boundary and the mean of the TAD-sep. score of surrounding bins. '
                           'The delta value reduces spurious boundaries that are shallow, which usually '
                           'occur at the center of large TADs when the TAD-sep. score is flat. Higher '
                           'delta threshold values produce more conservative boundary estimations. By '
                           'default a value of 0.01 is used.',
                           type=float,
                           default=0.01)

    parserOpt.add_argument('--minBoundaryDistance',
                           help='Minimum distance between boundaries (in bp). This parameter can be '
                           'used to reduce spurious boundaries caused by noise.',
                           type=int)

    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes and order in which the '
                           'chromosomes should be plotted. This option '
                           'overrides --region and --region2.',
                           nargs='+')

    parserOpt.add_argument('--numberOfProcessors', '-p',
                           help='Number of processors to use ',
                           type=int,
                           default=1)

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_matrix_wrapper(args):
    return compute_matrix(*args)


def get_cut_weight_by_bin_id(matrix, cut, depth, return_mean=False):
    """
    like get_cut_weight which is the 'diamond' representing the counts
    from a region -dept to +depth from the given bin position (cut):

    """
    if cut < 0 or cut > matrix.shape[0]:
        return None
    # the range [start:i] should have running window
    # length elements (i is excluded from the range)
    start = max(0, cut - depth)
    # same for range [i+1:end] (i is excluded from the range)
    end = min(matrix.shape[0], cut + depth)

    # the idea is to evaluate the interactions
    # between the upstream neighbors with the
    # down stream neighbors. In other words
    # the inter-domain interactions
    if return_mean is True:
        return matrix[start:cut, cut:end].mean()
    else:
        return matrix[start:cut, cut:end].todense().A1


def get_idx_of_bins_at_given_distance(hic_matrix, idx, window_len):
    """
    Returns the right and left indices that are at a given distance in bp
    from `idx`
    Parameters
    ----------
    idx reference index
    window_len distance in bp

    Returns
    -------
    tuple, with left and right bin indices
    """

    # find the bin that is closer to window_len
    # the range [start:i] should have running window
    # length elements (i is excluded from the range)
    chrom, cut_start, cut_end, _ = hic_matrix.getBinPos(idx)

    left_start = max(0, cut_start - window_len)
    left_idx = hic_matrix.getRegionBinRange(chrom, left_start, left_start + 1)[0]
    # chr_end_pos = hic_matrix.get_chromosome_sizes()[chrom]
    # if ?ring(chrom)
    chromosome_size = hic_matrix.get_chromosome_sizes()
    chrom = check_chrom_str_bytes(chromosome_size, chrom)
    # if type(next(iter(chromosome_size))) != type(chrom):
    #     if type(next(iter(chromosome_size))) is str:
    #         chrom = toString(chrom)
    #     elif type(next(iter(chromosome_size))) is bytes:
    #         chrom = toBytes(chrom)
    chr_end_pos = chromosome_size[chrom]

    right_end = min(chr_end_pos, cut_end + window_len) - 1
    right_idx = hic_matrix.getRegionBinRange(chrom, right_end, right_end)[0]
    return left_idx, right_idx


def get_cut_weight(hic_matrix, cut, window_len, return_mean=False):
    """
    'diamond' (when matrix is rotated 45 degrees) representing the counts
    from a region -window_len to +window_len from the given bin position (cut):

    """
    if cut < 0 or cut > hic_matrix.matrix.shape[0]:
        return None

    try:
        left_idx, right_idx = get_idx_of_bins_at_given_distance(hic_matrix, cut, window_len)
    except TypeError:
        # log.warn("Problem with cut: {}, window length: {}".format(cut, window_len))
        return None

    if return_mean is True:
        if hic_matrix.matrix[left_idx:cut, cut:right_idx].nnz == 0:
            # i.e. if the submatrix is emtpy. Calling the mean on an empty matrix
            # triggers a RuntimeWarning error
            return 0
        else:
            return hic_matrix.matrix[left_idx:cut, cut:right_idx].todense().mean()
    else:
        return hic_matrix.matrix[left_idx:cut, cut:right_idx].todense().A1


def get_triangle(hic_matrix, cut, window_len, return_mean=False):
    """
    like get_cut_weight which is the 'diamond' representing the counts
    from a region -window_len to +window_len from the given bin position (cut):

    """
    if cut < 0 or cut > hic_matrix.matrix.shape[0]:
        return None

    left_idx, right_idx = get_idx_of_bins_at_given_distance(hic_matrix, cut, window_len)

    def remove_lower_triangle(matrix):
        """
        remove all values in the lower triangle of a matrix
        """
        return matrix[np.triu_indices_from(matrix)].A1

    edges_left = remove_lower_triangle(hic_matrix.matrix[left_idx:cut, :][:, left_idx:cut].todense())
    edges_right = remove_lower_triangle(hic_matrix.matrix[cut:right_idx, :][:, cut:right_idx].todense())
#    if cut > 1000:
#        import ipdb;ipdb.set_trace()
    return np.concatenate([edges_left, edges_right])


def get_incremental_step_size(min_win_size, max_win_size, start_step_len):
    """
    generates a list of incremental windows sizes (measured in bins)

    :param min_win_size: starting window size
    :param max_win_size: end window size
    :param start_step_len: start length
    :return: incremental_step list of bin lengths
    """
    incremental_step = []
    step = -1
    while 1:
        step += 1
        inc_step = min_win_size + int(start_step_len * (step ** 1.5))
        if step > 1 and inc_step == incremental_step[-1]:
            continue
        if inc_step > max_win_size:
            break
        incremental_step.append(inc_step)
    return incremental_step


def compute_matrix(bins_list, min_win_size=8, max_win_size=50, step_len=2):
    """
    Receives a number of bins for which the tad-score should be computed.

    Parameters
    ----------
    hic_ma Hi-C matrix object
    bins_list list of bins to process
    min_win_size
    max_win_size
    step_len

    Returns
    -------

    """

    """
    :param hic_ma: Hi-C matrix object from HiCMatrix
    :param outfile: String, path of a file to save the conductance
                matrix in *bedgraph matrix* format
    :return: (chrom, start, end, matrix)
    """
    global hic_ma
    positions_array = []
    cond_matrix = []
    incremental_step = get_incremental_step_size(min_win_size, max_win_size, step_len)
    # if type(next(iter(self.interval_trees))) is np.bytes_:
    #     chrname = toBytes(chrname)
    # else:
    #     chrname = toString(chrname)
    # print("cut_intervals", hic_ma.cut_intervals)
    for cut in bins_list:
        chrom, chr_start, chr_end, _ = hic_ma.cut_intervals[cut]
        # get conductance
        # for multiple window lengths at a time
        mult_matrix = [get_cut_weight(hic_ma, cut, depth, return_mean=True) for depth in incremental_step]
        if any(x is None or np.isnan(x) for x in mult_matrix):
            # skip problematic cases
            continue

        cond_matrix.append(mult_matrix)

        positions_array.append((chrom, chr_start, chr_end))
    chrom, chr_start, chr_end = zip(*positions_array)
    cond_matrix = np.vstack(cond_matrix)

    return chrom, chr_start, chr_end, cond_matrix


class HicFindTads(object):

    def __init__(self, matrix, num_processors=1, max_depth=None, min_depth=None, step=None, delta=0.01,
                 min_boundary_distance=None, use_zscore=True, p_correct_for_multiple_testing="fdr", p_threshold_comparisons=0.01,
                 pChromosomes=None):
        """

        Parameters
        ----------
        matrix  Either a filename or a Hi-C Matrix object
        num_processors
        max_depth max window distance to consider (total window length is 2* max depth)
        min_depth min window to consider (total window length is 2* max min)
        step progression step from min_depth to max depth. The value give is for the first step, wich iteratively grows
                exponentially (exponent = 1.5) until it reaches the max_depth value.
        delta
        min_boundary_distance
        use_zscore boolean. By default is true. Set to other option
        pCorrectForMultipleTesting Multiple comparisons method: FDR, Bonferroni or None
        pThresholdComparisons The threshold for the Multiple comparisons. It is used as p-value for Bonferroni or as q-value for FDR.
        pChromosomes The chromomes that should be included for the analysis.
        """

        # if matrix is string, loaded, else, assume is a HiCMatrix object

        self.set_matrix(matrix, pChromosomes)
        if max_depth is not None and min_depth is not None and max_depth <= min_depth:
            log.error("Please check that maxDepth is larger than minDepth.")
            exit()

        self.num_processors = num_processors
        self.max_depth = max_depth
        self.min_depth = min_depth
        self.step = step
        self.delta = delta
        self.min_boundary_distance = min_boundary_distance
        self.use_zscore = use_zscore
        self.binsize = self.hic_ma.getBinSize()
        self.bedgraph_matrix = None
        self.boundaries = None
        self.set_variables()
        self.correct_for_multiple_testing = p_correct_for_multiple_testing
        self.threshold_comparisons = p_threshold_comparisons

    def set_matrix(self, pMatrix, pChromosomes):
        if isinstance(pMatrix, str):
            self.hic_ma = hm.hiCMatrix(pMatrix)
        else:
            self.hic_ma = pMatrix

        if pChromosomes is not None:
            valid_chromosomes = []
            invalid_chromosomes = []
            log.debug('args.chromosomeOrder: {}'.format(pChromosomes))
            log.debug("ma.chrBinBoundaries {}".format(self.hic_ma.chrBinBoundaries))
            if sys.version_info[0] == 3:
                pChromosomes = toBytes(pChromosomes)
            for chrom in toString(pChromosomes):
                if chrom in self.hic_ma.chrBinBoundaries:
                    valid_chromosomes.append(chrom)
                else:
                    invalid_chromosomes.append(chrom)

            if len(invalid_chromosomes) > 0:
                log.warning("WARNING: The following chromosome/scaffold names were not found. Please check"
                            "the correct spelling of the chromosome names. \n")
                log.warning("\n".join(invalid_chromosomes))
            self.hic_ma.reorderChromosomes(valid_chromosomes)

    def set_variables(self):
        """
        Checks the value of the max_depth, min_depth and step variables, setting default parameters
        if neccesary.

        Returns
        -------
        None
        """

        if self.max_depth is None:
            if self.binsize < 1000:
                self.max_depth = self.binsize * 60
            elif 1000 <= self.binsize < 20000:
                self.max_depth = self.binsize * 40
            else:
                self.max_depth = self.binsize * 10
        elif self.max_depth < self.binsize * 5:
            log.error("Please specify a --maxDepth that is at least 5 times larger than the matrix bin size")
            exit(1)

        if self.min_depth is None:
            if self.binsize < 1000:
                self.min_depth = self.binsize * 30
            elif 1000 <= self.binsize < 20000:
                self.min_depth = self.binsize * 10
            else:
                self.min_depth = self.binsize * 5
        elif self.min_depth < self.binsize * 3:
            log.error("Please specify a --minDepth that is at least 3 times larger than the matrix bin size")
            exit(1)

        if self.step is None:
            if self.binsize < 1000:
                self.step = self.binsize * 4
            else:
                self.step = self.binsize * 2

        elif self.step < self.binsize:
            log.error("Please specify a --step that is at least the size of the matrix bin size")
            exit(1)

        # print parameters used
        log.debug("max depth:\t{}\n".format(self.max_depth))
        log.debug("min depth:\t{}\n".format(self.min_depth))
        log.debug("step:\t{}\n".format(self.step))
        log.debug("bin size:\t{}\n".format(self.binsize))

    @staticmethod
    def peakdetect(y_axis, x_axis=None, lookahead=3, delta=0, chrom=None):
        """
        Based on the MATLAB script at:
        http://billauer.co.il/peakdet.html

        function for detecting local maximum and minimum in a signal.
        Discovers peaks by searching for values which are surrounded by lower
        or larger values for maximum and minimum respectively

        keyword arguments:
        :param: y_axis -- A list containing the signal over which to find peaks
        :param: x_axis -- (optional) A x-axis whose values correspond to the y_axis list
            and is used in the return to specify the position of the peaks. If
            omitted an index of the y_axis is used. (default: None)
        :param: lookahead -- (optional) distance to look ahead from a peak candidate to
            determine if it is the actual peak
        :param: delta -- (optional) this specifies a minimum difference between a peak and
            the following points, before a peak may be considered a peak. Useful
            to hinder the function from picking up false peaks towards to end of
            the signal. To work well delta should be set to delta >= RMSnoise * 5.


        :return: -- two lists [max_peaks, min_peaks] containing the positive and
            negative peaks respectively. Each cell of the lists contains a tuple
            of: (position, peak_value)
            to get the average peak value do: np.mean(max_peaks, 0)[1] on the
            results to unpack one of the lists into x, y coordinates do:
            x, y = zip(*tab)
        """
        max_peaks = []
        min_peaks = []
        dump = []   # Used to pop the first hit which almost always is false

        # check input data
        if x_axis is None:
            x_axis = np.arange(len(y_axis))

        if len(y_axis) != len(x_axis):
            raise ValueError('Input vectors y_axis and x_axis must have same length')

        # store data length for later use

        if not (np.isscalar(delta) and delta >= 0):
            raise ValueError("delta must be a positive number")

        # maximum and minimum candidates are temporarily stored in
        # min_x and min_y respectively
        min_y, max_y = np.Inf, -np.Inf
        max_pos, min_pos = None, None
        search_for = None
        # Only detect peak if there is 'lookahead' amount of points after it
        prev_chrom = None
        for index, (x, y) in enumerate(zip(x_axis[:-lookahead], y_axis[:-lookahead])):
            assert -np.inf < y < np.inf, "Error, infinity value detected for value at position {}".format(index)
            if prev_chrom is None:
                prev_chrom = chrom[index]
            if prev_chrom != chrom[index]:
                # reset variables at the start of a chromosome
                min_y, max_y = np.Inf, -np.Inf
                max_pos, min_pos = None, None
                search_for = None

            prev_chrom = chrom[index]

            if y > max_y:
                max_y = y
                max_pos = x
            if y < min_y:
                min_y = y
                min_pos = x

            # look for max
            if y < max_y - delta and max_y != np.Inf and search_for != 'min':
                # Maximum peak candidate found
                # look ahead in signal to ensure that this is a peak and not jitter
                if y_axis[index:index + lookahead].max() < max_y:
                    max_peaks.append([max_pos, max_y])
                    dump.append(True)
                    # set algorithm to only find minimum now
                    max_y = y
                    min_y = y
                    min_pos = x
                    search_for = 'min'
                    continue

            # look for min
            if y > min_y + delta and min_y != -np.Inf and search_for != 'max':
                # Minimum peak candidate found
                # look ahead in signal to ensure that this is a peak and not jitter
                if y_axis[index:index + lookahead].min() > min_y:
                    min_peaks.append([min_pos, min_y])
                    dump.append(False)
                    # set algorithm to only find maximum now
                    min_y = y
                    max_y = y
                    max_pos = x
                    search_for = 'max'

        # Remove the false hit on the first value of the y_axis
        try:
            if dump[0]:
                max_peaks.pop(0)
            else:
                min_peaks.pop(0)
            del dump
        except IndexError:
            # no peaks were found, should the function return empty lists?
            pass

        return [max_peaks, min_peaks]

    @staticmethod
    def delta_wrt_window(min_idx_list, matrix_avg, chrom, window_len=10):
        """
        Computes the local 'delta' for each minima identified. The local
        delta is defined as the difference between the TAD-separation score at the minimum with respect
        to each of the TAD-separation scores of the `window_len` to the left, plus the `window_len` to
        the right of each minimum. The minimum itself is excluded.

        Parameters
        ----------
        matrix TAD score matrix for different window sizes
        chrom list of chrom names for each matrix row
        window_len the length of the window to look for the local TAD-score average

        Returns
        -------
        A dict of each minima delta to the mean of the TAD-separation score surrounding the minima. The keys
        are the min_idx

        """

        # compute the start and end points of the chromosomes
        # and store it as tuples in chrom_ranges list (e.g. chrom_ranges = [(0, 29254), (29254, 60006), ...]
        unique_chroms, chr_start_idx = np.unique(chrom, return_index=True)
        chr_start_idx = np.concatenate([chr_start_idx, [len(chrom) - 1]])
        chr_start_idx = np.sort(chr_start_idx)
        chrom_ranges = [(chr_start_idx[x], chr_start_idx[x + 1]) for x in range(len(chr_start_idx) - 1)]

        delta_to_mean = {}
        for min_idx in min_idx_list:
            # check that the min_idx is not to close to any of the chromosome boundaries
            close_to_chrom_border = True
            for start_range, end_range in chrom_ranges:
                if start_range < min_idx < end_range:
                    if min_idx - window_len >= start_range and min_idx + window_len < end_range:
                        close_to_chrom_border = False
                        continue
            if close_to_chrom_border is True:
                delta_to_mean[min_idx] = np.nan

            else:
                local_tad_score = np.concatenate([matrix_avg[min_idx - window_len:min_idx + 3],
                                                  matrix_avg[min_idx + 4: min_idx + window_len]])
                delta_to_mean[min_idx] = local_tad_score.mean() - matrix_avg[min_idx]

        return delta_to_mean

    @staticmethod
    def find_consensus_minima(tad_score_matrix, lookahead=3, chrom=None):
        """
        Finds the minimum over the average values per column
        :param tad_score_matrix: TAD-separation score matrix
        :return:
        """

        tad_score_matrix_avg = tad_score_matrix.mean(axis=1)

        # compute local minima for the matrix average
        _max, _min = HicFindTads.peakdetect(tad_score_matrix_avg, lookahead=lookahead, chrom=chrom)
        if _min:
            min_idx, _ = zip(*_min)
        else:
            min_idx = []

        # get the delta for each boundary
        delta_to_mean = HicFindTads.delta_wrt_window(min_idx, tad_score_matrix_avg, chrom)

        return min_idx, delta_to_mean

    def hierarchical_clustering(self, boundary_list, clusters_cutoff=[]):
        """
        :param boundary_list: is a list of tuples each containing
        the location of a boundary. The order should be sorted
        and contain the following values:
            (chrom, start, value)
        :param clusters_cutoff: List of values to separate clusters. The
            clusters found at those value thresholds are returned.


        :return: z_value, clusters

        For z_value, the format used is similar as the scipy.cluster.hierarchy.linkage() function
        which is described as follows:

        A 4 by :math:`(n-1)` matrix ``z_value`` is returned. At the
        :math:`i`-th iteration, clusters with indices ``z_value[i, 0]`` and
        ``z_value[i, 1]`` are combined to form cluster :math:`n + i`. A
        cluster with an index less than :math:`n` corresponds to one of
        the :math:`n` original observations. The distance between
        clusters ``z_value[i, 0]`` and ``z_value[i, 1]`` is given by ``z_value[i, 2]``. The
        fourth value ``z_value[i, 3]`` represents the number of original
        observations in the newly formed cluster.

        The difference is that instead of a 4 times n-1 array, a
        6 times n-1 array is returned. Where positions 4, and 5
        correspond to the genomic coordinates of ``z_value[i, 0]`` and ``z_value[i, 1]``

        """
        # run the hierarchical clustering per chromosome
        if clusters_cutoff:
            # sort in reverse order
            clusters_cutoff = np.sort(np.unique(clusters_cutoff))[::-1]

        chrom, start, value = zip(*boundary_list)

        unique_chr, indices = np.unique(chrom, return_index=True)
        indices = indices[1:]  # the first element is not needed
        start_per_chr = np.split(start, indices)
        value_per_chr = np.split(value, indices)
        z_value = {}

        def get_domain_positions(boundary_position):
            """
            returns for each boundary a start,end position
            corresponding to each TAD
            :param boundary_position: list of boundary chromosomal positions
            :return: list of (start, end) tuples.
            """
            start_ = None
            domain_list = []
            for position in boundary_position:
                if start_ is None:
                    start_ = position
                    continue
                domain_list.append((start_, position))
                start_ = position

            return domain_list

        def find_in_clusters(clusters_, search_id):
            """
            Given a list of clusters (each cluster defined as as set,
            the function returns the position in which an id is found
            :param clusters_:
            :param search_id:
            :return:
            """
            for set_idx, set_of_ids in enumerate(clusters_):
                if search_id in set_of_ids:
                    return set_idx

        def cluster_to_regions(clusters_, chrom_name):
            """
            Transforms a list of sets of ids from the hierarchical
            clustering to genomic positions
            :param clusters_: cluster ids
            :param chrom_name: chromosome name
            :return: list of tuples with (chrom_name, start, end)

            Example:

            clusters = [set(1,2,3), set(4,5,10)]

            """
            start_list = []
            end_list = []
            for set_ in clusters_:
                if len(set_) == 0:
                    continue

                # the ids in the sets are created in such a
                # that the min id is the one with the smaller start position
                start_list.append(domains[min(set_)][0])
                end_list.append(domains[max(set_)][1])

            start_list = np.array(start_list)
            end_list = np.array(end_list)
            order = np.argsort(start_list)

            return zip([chrom_name] * len(order), start_list[order], end_list[order])

        return_clusters = {}  # collects the genomic positions of the clusters per chromosome
        # The values are a list, one for each cutoff.
        for chrom_idx, chrom_name in enumerate(unique_chr):
            z_value[chrom_name] = []
            return_clusters[chrom_name] = []
            clust_cutoff = clusters_cutoff[:]
            domains = get_domain_positions(start_per_chr[chrom_idx])
            clusters = [{x} for x in range(len(domains))]

            # initialize the cluster_x with the genomic position of domain centers
            cluster_x = [int(d_start + float(d_end - d_start) / 2) for d_start, d_end in domains]
            # number of domains should be equal to the number of values minus 1
            assert len(domains) == len(value_per_chr[chrom_idx]) - 1, "error"

            """
            domain:id
                 0            1               2            3
             |---------|---------------|----------------|----|
            values:id
             0         1               3                3    4
            values id after removing flanks
                       0               1                2
             """
            values = value_per_chr[chrom_idx][1:-1]  # remove flanking values that do not join TADs

            # from highest to lowest merge neighboring domains
            order = np.argsort(values)[::-1]
            for idx, order_idx in enumerate(order):
                if len(clust_cutoff) and idx + 1 < len(order) and \
                        values[order_idx] >= clust_cutoff[0] > values[order[idx + 1]]:
                    clust_cutoff = clust_cutoff[1:]  # remove first element
                    return_clusters[chrom_name].append(cluster_to_regions(clusters, chrom_name))
                # merge domains order_idx - 1 and order_idx
                left = find_in_clusters(clusters, order_idx)
                right = find_in_clusters(clusters, order_idx + 1)
                z_value[chrom_name].append((left, right, values[order_idx],
                                            len(clusters[left]) + len(clusters[right]),
                                            cluster_x[left], cluster_x[right]))

                # set as new cluster position the center between the two merged
                # clusters
                gen_dist = int(float(abs(cluster_x[left] - cluster_x[right])) / 2)
                cluster_x.append(min(cluster_x[left], cluster_x[right]) + gen_dist)

                clusters.append(clusters[left].union(clusters[right]))
                clusters[left] = set()
                clusters[right] = set()

        # convert return_clusters from a per chromosome dictionary to
        # a per cut_off dictionary merging all chromosomes in to one list.
        ret_ = {}  # dictionary to hold the clusters per cutoff. The key of
        # each item is the str(cutoff)

        for idx, cutoff in enumerate(clusters_cutoff):
            cutoff = str(cutoff)
            ret_[cutoff] = []
            for chr_name in return_clusters:
                try:
                    ret_[cutoff].extend(return_clusters[chr_name][idx])
                except IndexError:
                    pass

        return z_value, ret_

    def save_linkage(Z, file_name):
        """

        :param Z: Z has a format similar to the scipy.cluster.linkage matrix (see function
                    hierarchical_clustering).
        :param file_name: File name to save the results
        :return: None
        """

        try:
            file_h = open(file_name, 'w')
        except IOError:
            log.error("Can't save linkage file:\n{}".format(file_name))
            return

        count = 0
        for chrom, values in iteritems(Z):
            for id_a, id_b, distance, num_clusters, pos_a, pos_b in values:
                count += 1
                file_h.write('{}\t{}\t{}\tclust_{}'
                             '\t{}\t.\t{}\t{}\t{}\n'.format(chrom,
                                                            int(pos_a),
                                                            int(pos_b),
                                                            count,
                                                            distance,
                                                            id_a, id_b,
                                                            num_clusters))

    def get_domains(boundary_list):
        """
        returns for each boundary a chrom, start,end position
        corresponding to each TAD
        :param boundary_position: list of boundary chromosomal positions
        :return: list of (chrom, start, end, value) tuples.
        """
        prev_start = None
        prev_chrom = boundary_list[0][0]
        domain_list = []
        for chrom, start, value in boundary_list:
            if start is None:
                prev_start = start
                prev_chrom = chrom
                continue
            if prev_chrom != chrom:
                prev_chrom = chrom
                prev_start = None
                continue
            domain_list.append((chrom, prev_start, start, value))
            prev_start = start
            prev_chrom = chrom

        return domain_list

    def save_bedgraph_matrix(self, outfile):
        """
        Save matrix as chrom, start, end ,row, values separated by tab
        I call this a bedgraph matrix and the ending is .bm

        Returns
        -------
        None
        """
        # get params to save as part of the bedgraph file
        params = OrderedDict()
        params['step'] = self.step
        params['minDepth'] = self.min_depth
        params['maxDepth'] = self.max_depth
        params['binsize'] = self.binsize
        params_str = json.dumps(params, separators=(',', ':'))

        with open(outfile, 'w') as f:
            f.write("#" + params_str + "\n")
            for idx in range(len(self.bedgraph_matrix['chrom'])):
                matrix_values = "\t".join(np.char.mod('%f', self.bedgraph_matrix['matrix'][idx, :]))

                f.write("{}\t{}\t{}\t{}\n".format(toString(self.bedgraph_matrix['chrom'][idx]),
                                                  toString(self.bedgraph_matrix['chr_start'][idx]),
                                                  toString(self.bedgraph_matrix['chr_end'][idx]),
                                                  toString(matrix_values)))

    def save_clusters(clusters, file_prefix):
        """

        :param clusters: is a dictionary whose key is the cut of used to create it.
                         the value is a list of tuples, each representing
                          a genomec interval as ('chr', start, end).
        :param file_prefix: file prefix to save the resulting bed files
        :return: list of file names created
        """
        for cutoff, intervals in iteritems(clusters):
            fileh = open("{}_{}.bed".format(file_prefix, cutoff), 'w')
            for chrom, start, end in intervals:
                fileh.write("{}\t{}\t{}\t.\t0\t.\n".format(chrom, start, end))

    def save_domains_and_boundaries(self, prefix):

        # a boundary is added to the start and end of each chromosome
        # np.unique return index is used to quickly get
        # the indices at which the name of the chromosome changes (chrom, start, end should be  sorted)
        chrom = self.bedgraph_matrix['chrom']
        chr_start = self.bedgraph_matrix['chr_start']
        chr_end = self.bedgraph_matrix['chr_end']
        matrix = self.bedgraph_matrix['matrix']

        min_idx = self.boundaries['min_idx']
        delta_of_min = self.boundaries['delta']
        pvalue_of_min = self.boundaries['pvalues']

        unique_chroms, chr_start_idx = np.unique(chrom, return_index=True)

        # the trick to get the start positions using np.unique only works when
        # more than one chromosome is present
        if len(unique_chroms) == 1:
            chr_start_idx = [0]
            chr_end_idx = [len(chrom) - 1]
        else:
            chr_end_idx = chr_start_idx
            # the indices for the end of the chromosomes
            # are the the start indices - 1, with the exception
            # of the idx == 0 that is transformed into the length
            # of the chromosome to get the idx for the end of the
            # last chromosome
            chr_end_idx[chr_end_idx == 0] = len(chrom)
            chr_end_idx -= 1

        # put all indices together and sort
        min_idx = np.sort(np.concatenate([chr_start_idx, chr_end_idx, min_idx]))

        mean_mat_all = matrix.mean(axis=1)

        filtered_min_idx = []

        for idx in min_idx:
            # filter by delta and pvalue_thresholds
            if idx not in delta_of_min:
                delta_of_min[idx] = np.nan
            if self.correct_for_multiple_testing == 'fdr':
                if delta_of_min[idx] >= self.delta and idx in pvalue_of_min and pvalue_of_min[idx] <= self.pvalueFDR:
                    filtered_min_idx += [idx]
            elif self.correct_for_multiple_testing == 'bonferroni':
                if delta_of_min[idx] >= self.delta and idx in pvalue_of_min and pvalue_of_min[idx] <= self.threshold_comparisons:
                    filtered_min_idx += [idx]
            elif self.correct_for_multiple_testing == 'None':
                if delta_of_min[idx] >= self.delta and idx in pvalue_of_min and pvalue_of_min[idx] <= self.threshold_comparisons:
                    filtered_min_idx += [idx]

        if self.correct_for_multiple_testing == 'fdr':
            log.info("FDR correction. Number of boundaries for delta {}, qval {}: {}".format(self.delta, self.threshold_comparisons,
                                                                                             len(filtered_min_idx)))
        elif self.correct_for_multiple_testing == 'bonferroni':
            log.info("Bonferroni correction. Number of boundaries for delta {} and pval {}: {}".format(self.delta, self.threshold_comparisons,
                                                                                                       len(filtered_min_idx)))
        else:
            log.info("No multiple testing correction. Number of boundaries for delta {}: {}, used threshold: {}".format(self.delta, len(filtered_min_idx), self.threshold_comparisons))

        count = 1
        with open(prefix + '_boundaries.bed', 'w') as file_boundary_bin, open(prefix + '_domains.bed', 'w') as file_domains, open(prefix + '_boundaries.gff', 'w') as gff:
            for idx, min_bin_id in enumerate(filtered_min_idx):
                # skip if the start of the boundary
                # is the end of the chromosome
                if min_bin_id in chr_end_idx:
                    continue

                if min_bin_id not in delta_of_min:
                    delta_of_min[min_bin_id] = np.nan

                # 1. save boundaries at 1bp position
                right_bin_center = chr_start[min_bin_id] + int((chr_end[min_bin_id] - chr_start[min_bin_id]) / 2)

                left_bin_center = chr_start[min_bin_id - 1] + int((chr_end[min_bin_id - 1] - chr_start[min_bin_id - 1]) / 2)

                if chrom[min_bin_id] != chrom[min_bin_id - 1]:
                    continue

                # 2. save the position of the boundary range
                file_boundary_bin.write("{}\t{}\t{}\tB{:05d}\t{:.12f}\t.\n".format(toString(chrom[min_bin_id]),
                                                                                   left_bin_center,
                                                                                   right_bin_center,
                                                                                   min_bin_id,
                                                                                   mean_mat_all[min_bin_id]))

                # safe gff file that can contain more information
                gff.write("{chrom}\tHiCExplorer\tboundary\t{start}\t{end}\t{score:.12f}"
                          "\t.\t.\tID=B{id:05d};delta={delta:.12f};pvalue={pvalue:.12f};"
                          "tad_sep={score:.12f}\n".format(chrom=toString(chrom[min_bin_id]),
                                                          start=left_bin_center,
                                                          end=right_bin_center,
                                                          delta=delta_of_min[min_bin_id],
                                                          pvalue=pvalue_of_min[min_bin_id],
                                                          score=mean_mat_all[min_bin_id],
                                                          id=min_bin_id))

                start = chr_start[min_bin_id]
                # check that the next boundary exists and is in the same chromosome
                if idx + 1 == len(filtered_min_idx) or chrom[min_bin_id] != chrom[filtered_min_idx[idx + 1]]:
                    continue

                end = chr_start[filtered_min_idx[idx + 1]]

                # 2. save domain intervals
                if count % 2 == 0:
                    rgb = '51,160,44'
                else:
                    rgb = '31,120,180'

                file_domains.write("{0}\t{1}\t{2}\tID_{6}_{3}\t{4:.12f}\t.\t{1}\t{2}\t{5}\n".format(toString(chrom[min_bin_id]),
                                                                                                    start, end, count,
                                                                                                    mean_mat_all[min_bin_id],
                                                                                                    rgb, self.delta))

                count += 1

        # save track with mean values in bedgraph format
        with open(prefix + '_score.bedgraph', 'w') as tad_score:
            for idx in range(1, len(chrom)):
                right_bin_center = chr_start[idx] + int((chr_end[idx] - chr_start[idx]) / 2)
                left_bin_center = chr_start[idx - 1] + int((chr_end[idx - 1] - chr_start[idx - 1]) / 2)
                if right_bin_center <= left_bin_center:
                    # this condition happens at chromosome borders
                    continue
                tad_score.write("{}\t{}\t{}\t{:.12f}\n".format(toString(chrom[idx]), left_bin_center, right_bin_center,
                                                               mean_mat_all[idx]))

    def compute_spectra_matrix(self, perchr=True):
        """
        Uses multiple processors to compute the TAD-score

        Returns
        -------
        bed matrix (chrom start end value1 value2 ... value_n)

        """

        # remove self counts
        log.info('removing diagonal values\n')
        self.hic_ma.diagflat(value=0)

        # mask bins without any information
        self.hic_ma.maskBins(self.hic_ma.nan_bins)
        orig_intervals = self.hic_ma.cut_intervals[:]

        if self.use_zscore:
            # use zscore matrix
            log.info("Computing z-score matrix...\n")
            self.hic_ma.convert_to_zscore_matrix(maxdepth=self.max_depth * 2.5, perchr=perchr)

        # extend remaining bins to remove gaps in
        # the matrix
        new_intervals = enlarge_bins(self.hic_ma.cut_intervals)

        # rebuilt bin positions if necessary
        if new_intervals != orig_intervals:
            self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = self.hic_ma.intervalListToIntervalTree(new_intervals)
            self.hic_ma.cut_intervals = new_intervals
            self.hic_ma.orig_cut_intervals = new_intervals
            self.hic_ma.orig_bin_ids = list(range(len(new_intervals)))
            self.hic_ma.nan_bins = []

        if self.min_depth % self.hic_ma.getBinSize() != 0:
            log.warn('Warning. specified *depth* is not multiple of the '
                     'Hi-C matrix bin size ({})\n'.format(self.hic_ma.getBinSize()))
        if self.step % self.hic_ma.getBinSize() != 0:
            log.warn('Warning. specified *step* is not multiple of the '
                     'Hi-C matrix bin size ({})\n'.format(self.hic_ma.getBinSize()))

        self.binsize = self.hic_ma.getBinSize()

        log.info("Computing TAD-separation scores...\n")
        min_depth_in_bins = int(self.min_depth / self.binsize)
        max_depth_in_bins = int(self.max_depth / self.binsize)
        step_in_bins = int(self.step / self.binsize)
        if step_in_bins == 0:
            log.error("Please select a step size larger than {}".format(self.binsize))
            exit(1)

        incremental_step = get_incremental_step_size(self.min_depth, self.max_depth, self.step)

        log.info("computing spectrum for window sizes between {} ({} bp)"
                 "and {} ({} bp) at the following window sizes {} {}\n".format(min_depth_in_bins,
                                                                               self.binsize * min_depth_in_bins,
                                                                               max_depth_in_bins,
                                                                               self.binsize * max_depth_in_bins,
                                                                               step_in_bins, incremental_step))
        if min_depth_in_bins <= 1:
            log.error('ERROR\nminDepth length too small. Use a value that is at least '
                      'twice as large as the bin size which is: {}\n'.format(self.binsize))
            exit(0)

        if max_depth_in_bins <= 1:
            log.error('ERROR\nmaxDepth length too small. Use a value that is larger '
                      'than the bin size which is: {}\n'.format(self.binsize))
            exit(0)
        # work only with the upper matrix
        # and remove all pixels that are beyond
        # 2 * max_depth_in_bis which are not required
        # (this is done by subtracting a second sparse matrix
        # that contains only the upper matrix that wants to be removed.
        limit = 2 * max_depth_in_bins
        self.hic_ma.matrix = sparse.triu(self.hic_ma.matrix, k=0, format='csr') - sparse.triu(self.hic_ma.matrix,
                                                                                              k=limit, format='csr')
        self.hic_ma.matrix.eliminate_zeros()

        func = compute_matrix_wrapper
        TASKS = []
        bins_to_consider = []
        # to speed up parallel computation the self.hic_ma (HiCMatrix object) is converted into a global object.
        global hic_ma
        hic_ma = self.hic_ma
        for chrom in list(self.hic_ma.chrBinBoundaries):
            bins_to_consider.extend(list(range(*self.hic_ma.chrBinBoundaries[chrom])))

        for idx_array in np.array_split(bins_to_consider, self.num_processors):
            TASKS.append((idx_array, self.min_depth, self.max_depth, self.step))

        if self.num_processors > 1:
            pool = multiprocessing.Pool(self.num_processors)
            log.info("Using {} processors\n".format(self.num_processors))
            res = pool.map_async(func, TASKS).get(9999999)
            pool.close()
        else:
            res = map(func, TASKS)

        chrom = []
        chr_start = []
        chr_end = []
        matrix = []
        for _chrom, _chr_start, _chr_end, _matrix in res:
            chrom.extend(toString(_chrom))
            chr_start.extend(_chr_start)
            chr_end.extend(_chr_end)
            matrix.append(_matrix)

        matrix = np.vstack(matrix)

        self.bedgraph_matrix = {'chrom': np.array(chrom),
                                'chr_start': np.array(chr_start).astype(int),
                                'chr_end': np.array(chr_end).astype(int),
                                'matrix': matrix}

    def load_bedgraph_matrix(self, filename, pChromosomes=None):
        # load spectrum matrix:
        matrix = []
        chrom_list = []
        start_list = []
        end_list = []
        with open(filename, 'r') as fh:
            for line in fh:
                # if type(line)
                if line.startswith("#"):
                    # recover the parameters used to generate the spectrum_matrix
                    parameters = json.loads(line[1:].strip())
                    continue
                fields = line.strip().split('\t')
                chrom, start, end = fields[0:3]
                if pChromosomes is not None:
                    if chrom in pChromosomes:
                        chrom_list.append(chrom)
                        start_list.append(int(float(start)))
                        end_list.append(int(float(end)))
                        matrix.append(map(float, fields[3:]))
                else:
                    chrom_list.append(chrom)
                    start_list.append(int(float(start)))
                    end_list.append(int(float(end)))
                    matrix.append(map(float, fields[3:]))
        self.min_depth = parameters['minDepth']
        self.max_depth = parameters['maxDepth']
        self.step = parameters['step']
        self.binsize = parameters['binsize']

        matrix = np.vstack(matrix)
        self.bedgraph_matrix = {'chrom': np.array(chrom_list),
                                'chr_start': np.array(start_list).astype(int),
                                'chr_end': np.array(end_list).astype(int),
                                'matrix': matrix}

    def min_pvalue(self, min_idx):
        """
        For each putative local minima, find the -window_len diammond and the +window_len diamond
        and compare with the local minima using wilcoxon rank sum.

        Parameters
        ----------
        min_idx list of local minima

        Returns
        -------
        list of p-values per each local minima

        """

        log.info("Computing p-values for window length: {}\n".format(self.min_depth))
        pvalues = []
        chrom = self.bedgraph_matrix['chrom']
        chr_start = self.bedgraph_matrix['chr_start']
        chr_end = self.bedgraph_matrix['chr_end']
        window_len = self.min_depth

        from scipy.stats import ranksums

        new_min_idx = []
        for idx in min_idx:

            matrix_idx = self.hic_ma.getRegionBinRange(chrom[idx], chr_start[idx], chr_end[idx])
            if matrix_idx is None:
                continue
            else:
                matrix_idx = matrix_idx[0]

            new_min_idx += [idx]
            min_chr, min_start, min_end, _ = self.hic_ma.getBinPos(matrix_idx)
            assert toString(chrom[idx]) == toString(min_chr) and chr_start[idx] == min_start and chr_end[idx] == min_end
            left_idx, right_idx = get_idx_of_bins_at_given_distance(self.hic_ma, matrix_idx, window_len)

            left = get_cut_weight(self.hic_ma, left_idx, window_len)
            right = get_cut_weight(self.hic_ma, right_idx, window_len)
            boundary = get_cut_weight(self.hic_ma, matrix_idx, window_len)

            if left is None:
                left = []
            if right is None:
                right = []

            if len(left) == 0 and len(right) == 0:
                pval = np.nan

            elif boundary is None or len(boundary) == 0 or len(left) == 0 or len(right) == 0:
                pval = np.nan

            else:

                try:
                    pval1 = ranksums(boundary, left)[1]
                    pval2 = ranksums(boundary, right)[1]
                    pval = min(pval1, pval2)
                except ValueError:
                    # this condition happens when boundary, left and right are only composed of 'zeros'
                    pval = np.nan
            pvalues += [pval]

        assert len(pvalues) == len(new_min_idx)

        # fdr
        if self.correct_for_multiple_testing == 'fdr':

            pvalues = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
            pvalues_ = sorted(pvalues)
            largest_p_i = 0
            for i, p in enumerate(pvalues_):
                if p <= (self.threshold_comparisons * (i + 1) / len(pvalues_)):
                    if p >= largest_p_i:
                        largest_p_i = p
            self.pvalueFDR = largest_p_i
        elif self.correct_for_multiple_testing == 'bonferroni':
            # bonferroni correction
            pvalues = np.array(pvalues) * len(pvalues)
            to_one_index_values = np.array([e > 1 if ~np.isnan(e) else False for e in pvalues])
            if len(to_one_index_values) > 0:
                pvalues[to_one_index_values] = 1

        return OrderedDict(zip(new_min_idx, pvalues))

    def find_boundaries(self):

        # perform some checks
        avg_bin_size = np.median(self.bedgraph_matrix['chr_end'] - self.bedgraph_matrix['chr_start'])

        # compute lookahead (in number of bins)
        if self.min_boundary_distance is None:
            self.min_boundary_distance = avg_bin_size * 4

        lookahead = int(self.min_boundary_distance / avg_bin_size)
        if lookahead < 1:
            raise ValueError("minBoundaryDistance must be '1' or above in value")

        min_idx, delta = HicFindTads.find_consensus_minima(self.bedgraph_matrix['matrix'], lookahead=lookahead,
                                                           chrom=self.bedgraph_matrix['chrom'])

        pvalues = self.min_pvalue(min_idx)

        if len(min_idx) <= 10:
            mat_mean = self.bedgraph_matrix['matrix'].mean(axis=1)
            m_mean = mat_mean.mean()
            m_median = np.median(mat_mean)
            m_75 = np.percentile(mat_mean, 75)
            m_25 = np.percentile(mat_mean, 25)

            msg = ("Please check the parameters:\n"
                   " delta: {}\n"
                   " minBoundaryDistance: {}\n\n"
                   "TAD Score values:\n"
                   " mean: {:.3f}\n"
                   " median: {:.3f}\n"
                   " 1st quartile: {:.3f}\n"
                   " 3rd quartile: {:.3f}\n".format(self.delta, self.min_boundary_distance,
                                                    m_mean, m_median, m_25, m_75))

            if len(min_idx) == 0:
                log.error("\n*ERROR*\nNo boundaries were found. {}".format(msg))
                exit(1)
            else:
                log.info("Only {} boundaries found. {}".format(len(min_idx), msg))

        self.boundaries = {'min_idx': min_idx,
                           'delta': delta,
                           'pvalues': pvalues}


def print_args(args):
    """
    Print to stderr the parameters used
    Parameters
    ----------
    args

    Returns
    -------

    """
    for key, value in args._get_kwargs():
        log.info("{}:\t{}\n".format(key, value))


def main(args=None):

    args = parse_arguments().parse_args(args)
    ft = HicFindTads(args.matrix, num_processors=args.numberOfProcessors, max_depth=args.maxDepth,
                     min_depth=args.minDepth, step=args.step, delta=args.delta,
                     min_boundary_distance=args.minBoundaryDistance, use_zscore=True,
                     p_correct_for_multiple_testing=args.correctForMultipleTesting, p_threshold_comparisons=args.thresholdComparisons,
                     pChromosomes=args.chromosomes)

    tad_score_file = args.outPrefix + "_tad_score.bm"
    zscore_matrix_file = args.outPrefix + "_zscore_matrix.h5"

    if args.TAD_sep_score_prefix is not None:
        tad_score_file = args.TAD_sep_score_prefix + "_tad_score.bm"
        zscore_matrix_file = args.TAD_sep_score_prefix + "_zscore_matrix.h5"
        # check that the given file exists
        if not os.path.isfile(tad_score_file):
            log.error("The given TAD_sep_score_prefix does not contain a valid TAD-separation score. Please check.\n"
                      "Could not find file {}".format(tad_score_file))
            exit(1)
        if not os.path.isfile(zscore_matrix_file):
            log.error("The given TAD_sep_score_prefix does not contain a valid z-score matrix. Please check.\n"
                      "Could not find file {}".format(zscore_matrix_file))
            exit(1)
        log.info("\nUsing existing TAD-separation score file: {}\n".format(tad_score_file))
        ft.set_matrix(zscore_matrix_file, args.chromosomes)
        ft.load_bedgraph_matrix(tad_score_file, args.chromosomes)

    elif not os.path.isfile(tad_score_file):
        ft.compute_spectra_matrix()
        # save z-score matrix that is needed for find TADs algorithm
        ft.hic_ma.save(args.outPrefix + "_zscore_matrix.h5")
        ft.save_bedgraph_matrix(tad_score_file)
    else:
        log.info("\nFound existing TAD-separation score file: {}\n".format(tad_score_file))
        log.info("This file will be used\n")
        ft.set_matrix(zscore_matrix_file, args.chromosomes)
        # ft.hic_ma = hm.hiCMatrix(zscore_matrix_file)
        ft.load_bedgraph_matrix(tad_score_file, args.chromosomes)

    ft.find_boundaries()
    ft.save_domains_and_boundaries(args.outPrefix)

    # turn of hierarchical clustering which is apparently not working.

    # if 2 == 1:
    #     boundary_list = [(hic_ma.cut_intervals[min_][0], hic_ma.cut_intervals[min_][2], mean_mat[min_]) for min_ in min_idx]

    #     Z, clusters = hierarchical_clustering(boundary_list, clusters_cutoff=[0.4, 0.3, 0.2])

    #     save_linkage(Z, args.outPrefix + '_linkage.bed')
    #     save_clusters(clusters, args.outPrefix)
