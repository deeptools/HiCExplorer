"""hicMergeMatrixBins: the tool is the C++ executable hicMergeMatrixBins.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). remove_nans_if_needed, running_window_merge and merge_bins
keep their 3.7.6 signatures and behaviour and call
hicexplorer._hicx.hicMergeMatrixBins (cpp/python/src/api/hicMergeMatrixBins.cpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import numpy as np  # noqa: F401 (the doctests use it)
import logging
log = logging.getLogger(__name__)

from hicmatrix import HiCMatrix as hm  # noqa: F401 (the doctests use it)
from hicexplorer.reduceMatrix import reduce_matrix  # noqa: F401 (the C++ hicMergeMatrixBins looks it up in this module)

from hicexplorer._hicx import hicMergeMatrixBins as _impl


def parse_arguments(args=None):
    """The argparse parser of hicMergeMatrixBins, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicMergeMatrixBins')


def remove_nans_if_needed(hic):
    return _impl.remove_nans_if_needed(hic)


def running_window_merge(hic_matrix, num_bins):
    """Creates a 'running window' merge without changing the
    original resolution of the matrix. The window size is
    defined by the num_bins that are merged. Num bins
    had to be an odd number such that equal amounts of left and
    right bins can be merged.
       a | b | c
       ---------
       d | e | f
       ---------
       g | h | i
    In this matrix, using a merge of num_bins 3,
    the merge is done as follows, a = a + b + d + e,
    e = a + b + c + d + e + f etc,

    >>> from scipy.sparse import csr_matrix, dia_matrix
    >>> row, col = np.triu_indices(5)
    >>> cut_intervals = [('a', 0, 10, 0.5), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('b', 40, 50, 1)]
    >>> hic = hm.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 1, 1 ],
    ... [ 1, 1 ]])

    make the matrix symmetric:
    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals[:2])
    >>> merge_matrix = running_window_merge(hic, 3)
    >>> print(merge_matrix.matrix.todense())
    [[3 3]
     [3 3]]

    >>> matrix = np.array([
    ... [ 1, 1, 1, 1 ],
    ... [ 1, 1, 1, 1 ],
    ... [ 1, 1, 1, 1 ],
    ... [ 1, 1, 1, 1 ]])

    make the matrix symmetric:
    >>> hic.matrix = csr_matrix(matrix)
    >>> hic.setMatrix(hic.matrix, cut_intervals)
    >>> merge_matrix = running_window_merge(hic, 3)
    >>> print(merge_matrix.matrix.todense())
    [[3 5 6 4]
     [5 6 8 6]
     [6 8 6 5]
     [4 6 5 3]]
    """
    return _impl.running_window_merge(hic_matrix, num_bins)


def merge_bins(hic, num_bins):
    """
    Merge the bins using the specified number of bins. This
    functions takes care to make new intervals

    Parameters
    ----------

    hic : HiCMatrix object

    num_bins : number of consecutive bins to merge.

    Returns
    -------

    A sparse matrix.

    Set up a Hi-C test matrix
    >>> from scipy.sparse import csr_matrix
    >>> row, col = np.triu_indices(5)
    >>> cut_intervals = [('a', 0, 10, 0.5), ('a', 10, 20, 1),
    ... ('a', 20, 30, 1), ('a', 30, 40, 0.1), ('b', 40, 50, 1)]
    >>> hic = hm.hiCMatrix()
    >>> hic.nan_bins = []
    >>> matrix = np.array([
    ... [ 50, 10,  5,  3,   0],
    ... [  0, 60, 15,  5,   1],
    ... [  0,  0, 80,  7,   3],
    ... [  0,  0,  0, 90,   1],
    ... [  0,  0,  0,  0, 100]], dtype=np.int32)

    make the matrix symmetric:
    >>> from scipy.sparse import dia_matrix

    >>> dia = dia_matrix(([matrix.diagonal()], [0]), shape=matrix.shape)
    >>> hic.matrix = csr_matrix(matrix + matrix.T - dia)
    >>> hic.setMatrix(hic.matrix, cut_intervals)

    run merge_matrix
    >>> merge_matrix = merge_bins(hic, 2)
    >>> merge_matrix.cut_intervals
    [('a', 0, 20, 0.75), ('a', 20, 40, 0.55), ('b', 40, 50, 1.0)]
    >>> merge_matrix.matrix.todense()
    matrix([[120,  28,   1],
            [ 28, 177,   4],
            [  1,   4, 100]], dtype=int32)
    """
    return _impl.merge_bins(hic, num_bins)


def main(args=None):
    """Run the C++ implementation of hicMergeMatrixBins with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicMergeMatrixBins')(args)
