"""hicPlotDistVsCounts, implemented in C++.

main() runs the C++ executable hicPlotDistVsCounts and parse_arguments() is
built from its --help-json (hicexplorer._cpp). compute_distance_mean and
from_bed_to_cut_interval keep the name, the signature and the behaviour they
had in HiCExplorer 3.7.6; their bodies call hicexplorer._hicx.hicPlotDistVsCounts
(cpp/python/src/api/hicPlotDistVsCounts.cpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# numpy and HiCMatrix are the names the doctest of from_bed_to_cut_interval uses
import numpy as np  # noqa: F401

from hicmatrix import HiCMatrix  # noqa: F401

import matplotlib
matplotlib.use('Agg')

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicPlotDistVsCounts as _impl


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
    return _impl.compute_distance_mean(hicmat, maxdepth, perchr, custom_cut_intervals)


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
    return _impl.from_bed_to_cut_interval(hicmat, fh)


def main(args=None):
    """Run the C++ implementation of hicPlotDistVsCounts with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotDistVsCounts')(args)
