"""hicPlotDistVsCounts: a thin layer over the C++ executable hicPlotDistVsCounts.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
import numpy as np
import pandas as pd
import argparse
import sys

from hicmatrix import HiCMatrix
from hicexplorer._version import __version__

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

from collections import OrderedDict
from past.builtins import zip

from scipy.sparse import triu

from .utilities import change_chrom_names

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicPlotDistVsCounts, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotDistVsCounts')


def compute_distance_mean(hicmat, maxdepth=None, perchr=False, custom_cut_intervals=None):
    """
    Compute average of values for each distance (only intra-chr or intra-unit defined in custom_cut_intervals)

    The caveat is that the average are only
    computed for non-zero values, although zero values that
    are not part of the sparse matrix are considered.

    For each diagonal the mean are calculated

    All unit whose names begin with _ignore_ are ignored

    Parameters
    ----------
    hicmat: HiCMatrix object
    maxdepth: maximum distance from the diagonal to consider. All contacts beyond this distance will not
                     be considered.
    perchr: bool to indicate if computations should be perform per chromosome
    custom_cut_intervals: by default the distance mean is computed for intra-chr
                          but another unit can be defined for example TADs.
                          then custom_cut_intervals should contains a list of tuple
                          with the same length as hic.cut_intervals for example:
                          if cut_intervals are [('a', 0, 10, 1), ('a', 10, 20, 1), ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
                          custom could be [('tad1', 0, 10, 1), ('tad1', 10, 20, 1), ('tad2', 0, 10, 1), ('tad2', 10, 20, 1), ('tad3', 0, 10, 1)]


    Returns
    -------
    dictionary where keys are either 'all' or the chromosome names when perchr is set to True
               and values are ordered dictionary where keys are distances and values are a
               2-tuple with the mean and the number of bins considered.

    """

    binsize = hicmat.getBinSize()

    if maxdepth:
        if maxdepth < binsize:
            exit("Please specify a maxDepth larger than bin size ({})".format(binsize))

        max_depth_in_bins = int(float(maxdepth * 1.5) / binsize)
        # work only with the upper matrix
        # and remove all pixels that are beyond
        # max_depth_in_bins
        # (this is done by subtracting a second sparse matrix
        # that contains only the upper matrix that wants to be removed.
        hicmat.matrix = triu(hicmat.matrix, k=0, format='csr') - \
            triu(hicmat.matrix, k=max_depth_in_bins, format='csr')
    else:
        hicmat.matrix = triu(hicmat.matrix, k=0, format='csr')

    hicmat.matrix.eliminate_zeros()

    chr_submatrix = OrderedDict()
    cut_intervals = OrderedDict()
    unit_sizes = OrderedDict()
    chrom_range = OrderedDict()

    if custom_cut_intervals is None:
        cut_intervals_genome_wide = hicmat.cut_intervals
    else:
        assert len(custom_cut_intervals) == len(hicmat.cut_intervals), \
            "The custom_cut_intervals should have the same size as hicmat.cut_intervals"
        cut_intervals_genome_wide = custom_cut_intervals

    if perchr:
        for chrname in hicmat.getChrNames():
            chr_range = hicmat.getChrBinRange(chrname)
            chr_submatrix[chrname] = hicmat.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]].tocoo()
            cut_intervals[chrname] = [cut_intervals_genome_wide[x] for x in range(chr_range[0], chr_range[1])]
            unit_sizes[chrname] = np.array([v[1] - v[0]
                                            for k, v in hicmat.intervalListToIntervalTree(cut_intervals[chrname])[1].items()
                                            if not k.startswith('_ignore_')])
            chrom_range[chrname] = (chr_range[0], chr_range[1])

    else:
        chr_submatrix['all'] = hicmat.matrix.tocoo()
        cut_intervals['all'] = cut_intervals_genome_wide
        unit_sizes['all'] = np.array([v[1] - v[0]
                                      for k, v in hicmat.intervalListToIntervalTree(cut_intervals_genome_wide)[1].items()
                                      if not k.startswith('_ignore_')])
        chrom_range['all'] = (0, hicmat.matrix.shape[0])

    mean_dict = {}

    for chrname, submatrix in chr_submatrix.items():
        log.info("processing chromosome {}\n".format(chrname))

        dist_list, chrom_list = hicmat.getDistList(submatrix.row, submatrix.col, HiCMatrix.hiCMatrix.fit_cut_intervals(cut_intervals[chrname]))
        # We filter out the interactions where the chrom starts with _ignore_
        dist_list = dist_list[[not chrom.startswith('_ignore_') for chrom in chrom_list]]

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
        sum_counts = np.bincount(dist_list, weights=submatrix.data[[not chrom.startswith('_ignore_') for chrom in chrom_list]])
        distance_len = np.bincount(dist_list)
        # compute the average for each distance
        mat_size = submatrix.shape[0]
        # compute mean value for each distance
        mu = {}
        # also store the number of bins used to do the average
        nb = {}
        zero_value_bins = []
        for bin_dist_plus_one, sum_value in enumerate(sum_counts):
            if maxdepth and bin_dist_plus_one == 0:  # this is for inter chromosomal counts
                # when max depth is set, the computation
                # of the total_inter is not accurate and is safer to
                # output np.nan
                mu[bin_dist_plus_one] = np.nan
                nb[bin_dist_plus_one] = np.nan
                continue

            if bin_dist_plus_one == 0:
                total_intra = mat_size ** 2 - sum([size ** 2 for size in unit_sizes[chrname]])
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
                diagonal_length = sum([size - (bin_dist_plus_one - 1) for size in unit_sizes[chrname]
                                       if size > (bin_dist_plus_one - 1)])

            # the diagonal length should contain the number of values at a certain distance.
            # If the matrix is dense, the distance_len[bin_dist_plus_one] correctly contains the number of values
            # If the matrix is equally spaced, then, the diagonal_length as computed before is accurate.
            # But, if the matrix is both sparse and with unequal bins, then none of the above methods is
            # accurate but the the diagonal_length as computed before will be closer.
            diagonal_length = max(diagonal_length, distance_len[bin_dist_plus_one])

            if diagonal_length == 0:
                mu[bin_dist_plus_one] = np.nan
                nb[bin_dist_plus_one] = 0
            else:
                mu[bin_dist_plus_one] = np.float64(sum_value) / diagonal_length
                nb[bin_dist_plus_one] = diagonal_length
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
        mean_dict[chrname] = OrderedDict([((k - 1) * binsize, (v, nb[k])) for k, v in mu.items() if k > 0 and
                                          (k - 1) * binsize <= maxdepth])
        # mean_dict[chrname]['intra_chr'] = mu[0]

    return mean_dict


def from_bed_to_cut_interval(hicmat, fh):
    """
    Generate cut_interval (list of 4-uple) compatible with a HiCMatrix
    where 'chrom' are units from the bed file

    Parameters
    ----------
    hicmat: HiCMatrix object
    bedfile: bed file with non-overlapping intervals and not more than 2 features should overlap a bin (for example domains.bed)

    Returns
    -------
    list of tuple like hicmat.cut_intervals

    >>> from scipy.sparse import csr_matrix
    >>> row, col = np.triu_indices(5)
    >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('a', 30, 40, 1),
    ... ('b', 20, 30, 1), ('b', 40, 50, 1)]
    >>> hic = HiCMatrix.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 1,  8,  5, 3,  0,  0],
    ... [ 0,  4, 15, 5,  1,  0],
    ... [ 0,  0,  0, 7,  2,  1],
    ... [ 0,  0,  0, 0,  1,  2],
    ... [ 0,  0,  0, 0, 10, 20],
    ... [ 0,  0,  0, 0,  0,  5]])

    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals)
    >>> tad_line='a\\t0\\t30\\nb\\t20\\t50'
    >>> with open('/tmp/test.bed', 'w') as fh:
    ...     fh.write(tad_line)
    14
    >>> fh = open('/tmp/test.bed', 'r')
    >>> from_bed_to_cut_interval(hic, fh)
    [('0', 0, 10, 1), ('0', 10, 20, 1), ('0', 20, 30, 1), ('_ignore_3', 0, 10, 1), ('1', 0, 10, 1), ('1', 20, 30, 1)]
    """
    original_cut_intervals = hicmat.cut_intervals
    new_cut_intervals = [()] * len(original_cut_intervals)
    chrom_list = hicmat.getChrNames()
    id = 0
    for line in fh:
        if line[0] == "#":
            continue
        fields = line.strip().split()
        if fields[0] not in chrom_list:
            if change_chrom_names(fields[0]) in chrom_list:
                fields[0] = change_chrom_names(fields[0])
            else:
                continue
        chrom, start, end = fields[:3]
        overlapping_bins = [i for s, e, i in hicmat.interval_trees[chrom].overlap(int(start), int(end))]
        original_start_bin_pos = original_cut_intervals[min(overlapping_bins)][1]
        for i in overlapping_bins:
            # Check it does not overlap another region
            if len(new_cut_intervals[i]) != 0:
                raise Exception("2 features must not overlap the same bin."
                                "{} is overlapped twice".format(original_cut_intervals[i]))
            # Fill the new_cut_interval
            new_cut_intervals[i] = (str(id), original_cut_intervals[i][1] - original_start_bin_pos,
                                    original_cut_intervals[i][2] - original_start_bin_pos,
                                    original_cut_intervals[i][3])
        id += 1
    assert id > 0, \
        "No region overlapped with bins."
    # Fill all intervals which are not in the fh with a 'chromosome' which starts with '_ignore_'
    for i in [i for i, interval in enumerate(new_cut_intervals) if len(interval) == 0]:
        new_cut_intervals[i] = ('_ignore_{}'.format(i), 0, original_cut_intervals[i][2] - original_cut_intervals[i][1], original_cut_intervals[i][3])
    return new_cut_intervals


if sys.version_info >= (3, 12):
    compute_distance_mean.__doc__ += """
    >>> from scipy.sparse import csr_matrix
    >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('a', 30, 40, 1),
    ... ('b', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]
    >>> hic = HiCMatrix.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 1,  8,  5,  3,  0,  0,  0],
    ... [ 0,  4, 15,  5,  1,  0,  0],
    ... [ 0,  0,  0,  7,  2,  0,  1],
    ... [ 0,  0,  0,  0,  1,  0,  2],
    ... [ 0,  0,  0,  0, 10,  0, 20],
    ... [ 0,  0,  0,  0,  0,  0,  0],
    ... [ 0,  0,  0,  0,  0,  0,  6]])

    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals)
    >>> compute_distance_mean(hic)
    {'all': OrderedDict({0: (3.0, 7), 10: (6.0, 5), 20: (10.0, 3), 30: (3.0, 1)})}
    >>> compute_distance_mean(hic, perchr=True)
    {'a': OrderedDict({0: (1.25, 4), 10: (10.0, 3), 20: (5.0, 2), 30: (3.0, 1)}), 'b': OrderedDict({0: (5.333333333333333, 3), 10: (0.0, 2), 20: (20.0, 1)})}
    >>> custom_cut = [('tad1', 0, 10, 1), ('tad1', 10, 20, 1), ('tad2', 0, 10, 1),
    ... ('tad2', 10, 20, 1), ('tad3', 0, 10, 1), ('tad3', 10, 20, 1), ('tad3', 20, 30, 1)]
    >>> compute_distance_mean(hic, custom_cut_intervals=custom_cut)
    {'all': OrderedDict({0: (3.0, 7), 10: (3.75, 4), 20: (20.0, 1)})}
    >>> compute_distance_mean(hic, perchr=True, custom_cut_intervals=custom_cut)
    {'a': OrderedDict({0: (1.25, 4), 10: (7.5, 2)}), 'b': OrderedDict({0: (5.333333333333333, 3), 10: (0.0, 2), 20: (20.0, 1)})}
    >>> custom_cut = [('_ignore_0', 0, 10, 1), ('0', 0, 10, 1),
    ... ('0', 10, 20, 1), ('_ignore_3', 0, 10, 1),
    ... ('1', 0, 10, 1), ('1', 10, 20, 1), ('1', 20, 30, 1)]
    >>> compute_distance_mean(hic, custom_cut_intervals=custom_cut)
    {'all': OrderedDict({0: (4.0, 5), 10: (5.0, 3), 20: (20.0, 1)})}
    """
else:
    compute_distance_mean.__doc__ += """
    >>> from scipy.sparse import csr_matrix
    >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('a', 30, 40, 1),
    ... ('b', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]
    >>> hic = HiCMatrix.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 1,  8,  5,  3,  0,  0,  0],
    ... [ 0,  4, 15,  5,  1,  0,  0],
    ... [ 0,  0,  0,  7,  2,  0,  1],
    ... [ 0,  0,  0,  0,  1,  0,  2],
    ... [ 0,  0,  0,  0, 10,  0, 20],
    ... [ 0,  0,  0,  0,  0,  0,  0],
    ... [ 0,  0,  0,  0,  0,  0,  6]])
    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals)
    >>> compute_distance_mean(hic)
    {'all': OrderedDict([(0, (3.0, 7)), (10, (6.0, 5)), (20, (10.0, 3)), (30, (3.0, 1))])}
    >>> compute_distance_mean(hic, perchr=True)
    {'a': OrderedDict([(0, (1.25, 4)), (10, (10.0, 3)), (20, (5.0, 2)), (30, (3.0, 1))]), 'b': OrderedDict([(0, (5.333333333333333, 3)), (10, (0.0, 2)), (20, (20.0, 1))])}
    >>> custom_cut = [('tad1', 0, 10, 1), ('tad1', 10, 20, 1), ('tad2', 0, 10, 1),
    ... ('tad2', 10, 20, 1), ('tad3', 0, 10, 1), ('tad3', 10, 20, 1), ('tad3', 20, 30, 1)]
    >>> compute_distance_mean(hic, custom_cut_intervals=custom_cut)
    {'all': OrderedDict([(0, (3.0, 7)), (10, (3.75, 4)), (20, (20.0, 1))])}
    >>> compute_distance_mean(hic, perchr=True, custom_cut_intervals=custom_cut)
    {'a': OrderedDict([(0, (1.25, 4)), (10, (7.5, 2))]), 'b': OrderedDict([(0, (5.333333333333333, 3)), (10, (0.0, 2)), (20, (20.0, 1))])}
    >>> custom_cut = [('_ignore_0', 0, 10, 1), ('0', 0, 10, 1),
    ... ('0', 10, 20, 1), ('_ignore_3', 0, 10, 1),
    ... ('1', 0, 10, 1), ('1', 10, 20, 1), ('1', 20, 30, 1)]
    >>> compute_distance_mean(hic, custom_cut_intervals=custom_cut)
    {'all': OrderedDict([(0, (4.0, 5)), (10, (5.0, 3)), (20, (20.0, 1))])}
    """


def main(args=None):
    """Run the C++ implementation of hicPlotDistVsCounts with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotDistVsCounts')(args)
