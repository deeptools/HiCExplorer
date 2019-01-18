from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from scipy.sparse import coo_matrix, dia_matrix, triu
import numpy as np
import time

import logging
log = logging.getLogger(__name__)


def reduce_matrix(matrix, bins_to_merge, use_triu=True, diagonal=False):
    """
    This function sums the rows and columns corresponding
    to the bins_to_merge, returning a new sparse
    matrix of size len(bins_to_merge).

    The function uses sparse matrices tricks to work very fast.
    The idea is that all the cells in the new reduced matrix
    whose sum needs to be computed are expressed as a vector
    of bins.  Using the numpy bincount function the sum of each bin
    is obtained. Then this sum vector is converted to a
    sparse matrix.

    Parameters
    ----------

    bins_to_merge : A list of lists. The values of the lists should
        correspond to the indices of the matrix.

    use_triu : use only the upper triangle of the matrix. This is to avoid double
               counting values close to the diagonal.
               E.g. the matrix [[1, 2],
                                 2, 1]]

               when merging the two bins, should be reduced to [4] and not [6], because
               the number '2' should be counted only once.
               When this option is used, the result is converted to a symmetric matrix

    diagonal : If set to true, then the main diagonal is preserved (not deleted). Only works
               when use_triu is set

    Returns
    -------

    A sparse matrix.

    >>> from scipy.sparse import csr_matrix
    >>> A = csr_matrix(np.array([[1,0],[0,1]]), dtype=np.int32)
    >>> reduce_matrix(A, [(0,1)], diagonal=True).todense()
    matrix([[2]], dtype=int32)
    >>> A = csr_matrix(np.array([[5,5,2,2,0],[0,5,2,2,1],
    ... [0,0,1,1,0], [0,0,0,1,0], [0,0,0,0,0]]), dtype=np.int32)
    >>> A.todense()
    matrix([[5, 5, 2, 2, 0],
            [0, 5, 2, 2, 1],
            [0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0]], dtype=int32)
    >>> ll = [(0,1), (2,3), (4,)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[15,  8,  1],
            [ 0,  3,  0],
            [ 0,  0,  0]], dtype=int32)

    >>> ll = [(0,1,2), (3,4)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[20,  6],
            [ 0,  1]], dtype=int32)

    >>> ll = [(0,1), (2,3), (4,)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=True).todense()
    matrix([[15,  8,  1],
            [ 8,  3,  0],
            [ 1,  0,  0]], dtype=int32)

    Test symmetric matrix
    >>> A = csr_matrix(np.array([[2,2,1,1,1,1],[2,2,1,1,1,1],
    ... [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1]]), dtype=np.int32)
    >>> A.todense()
    matrix([[2, 2, 1, 1, 1, 1],
            [2, 2, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1]], dtype=int32)
    >>> ll = [(0,1), (2,3), (4,5)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[8, 4, 4],
            [4, 4, 4],
            [4, 4, 4]], dtype=int32)

    >>> reduce_matrix(A, ll, diagonal=True, use_triu=True).todense()
    matrix([[6, 4, 4],
            [4, 3, 4],
            [4, 4, 3]], dtype=int32)

    >>> ll = [(0,1,2), (3,4,5)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[13,  9],
            [ 9,  9]], dtype=int32)

    Test symmetric matrix non consecutive
    >>> A = csr_matrix(np.array([[5,1,5,1,0],[0,1,2,1,1],
    ... [0,0,5,1,0], [0,0,0,1,0], [0,0,0,0,0]]), dtype=np.int32)
    >>> dia = dia_matrix(([A.diagonal()], [0]), shape=A.shape)
    >>> A= csr_matrix(A + A.T - dia, dtype=np.int32)
    >>> A.todense()
    matrix([[5, 1, 5, 1, 0],
            [1, 1, 2, 1, 1],
            [5, 2, 5, 1, 0],
            [1, 1, 1, 1, 0],
            [0, 1, 0, 0, 0]], dtype=int32)
    >>> ll = [(0,2), (1,3), (4,)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[20,  5,  0],
            [ 5,  4,  1],
            [ 0,  1,  0]], dtype=int32)

    Test removal of row/columns when the scaffold list does
    not contains all indices
    >>> ll = [(1,2), (3,)]
    >>> reduce_matrix(A, ll, diagonal=True, use_triu=False).todense()
    matrix([[10,  2],
            [ 2,  1]], dtype=int32)

    Test with float and nan
    >>> A = csr_matrix(np.array([[0.1,0.1,0.2,0.2,np.nan],
    ... [0,0.1,0.2,0.2,1.1],
    ... [0,0,0.2,0.2,0], [0,0,0,0.1,0], [0,0,0,0,0]]))
    >>> dia = dia_matrix(([A.diagonal()], [0]), shape=A.shape)
    >>> A= csr_matrix(A + A.T - dia)
    >>> print(A.todense())
    [[0.1 0.1 0.2 0.2 nan]
     [0.1 0.1 0.2 0.2 1.1]
     [0.2 0.2 0.2 0.2 0. ]
     [0.2 0.2 0.2 0.1 0. ]
     [nan 1.1 0.  0.  0. ]]

    >>> ll = [(0,1), (2,3), (4,)]
    >>> print(reduce_matrix(A, ll, diagonal=True, use_triu=False).todense())
    [[0.4 0.8 nan]
     [0.8 0.7 0. ]
     [nan 0.  0. ]]
    """

    if use_triu:
        ma = triu(matrix, k=0, format='coo')
    else:
        ma = matrix.tocoo()
    M = len(bins_to_merge)
    start_time = time.time()

    if M == ma.shape[0]:
        return matrix

    # check for nans and make a warning in case they are
    # present
    num_nan = len(np.flatnonzero(np.isnan(np.array(ma.data))))
    if num_nan > 0:
        log.warning("*Warning*\nmatrix contains {} NaN values.".format(num_nan))

    # each original col and row index is converted
    # to a new index based on the bins_to_merge.
    # For example, if rows 1 and 10 are to be merged
    # then all rows whose value is 1 or 10 are given
    # as new value the index in the bins_to_merge list.

    map_ = np.zeros(ma.shape[0], dtype=int) - 1     # -1 such that all cases not replaced by the next loop
    # can be identified later. Those cases that remain as -1
    # are for the rows/cols not appearing in the bins_to_merge
    for k, v in enumerate(bins_to_merge):
        for x in v:
            map_[x] = k

    new_row = np.take(map_, ma.row)
    new_col = np.take(map_, ma.col)

    # remove rows and cols with -1. Those
    # correspond to regions that are not part of
    # the bins_to_merge to be merged
    keep = (new_row > -1) & (new_col > -1)
    ma.data = ma.data[keep]
    new_row = new_row[keep]
    new_col = new_col[keep]

    # The following line converts each combination
    # of row and col into a list of bins. For example,
    # for each case in which row = 1 and col =3 or
    # row = 3 and col = 1 a unique bin id is given
    # The trick I use to get the unique bin ids is
    # to order and convert the pair (row,col) into
    # a complex number  ( for the pair 3,1 the
    # complex es 1,3j.
    uniq, bin_array = np.unique(new_row + 1j * new_col, return_inverse=True)

    elapsed_time = time.time() - start_time
    start_time = time.time()
    log.debug("time complex unique: {:.5f}".format(elapsed_time))
    # sum the bins. For each bin, all values associated
    # to the bin members are added. Those values are in
    # the ma.data array.
    sum_array = np.bincount(bin_array, weights=ma.data)

    # To reconstruct the new matrix the row and col
    # positions need to be obtained. They correspond
    # to the positions of the unique values in the bin_array.
    uniq, ind = np.unique(bin_array, return_index=True)
    new_row_ = new_row[ind]
    new_col_ = new_col[ind]

    result = coo_matrix((sum_array, (new_row_, new_col_)),
                        shape=(M, M), dtype=ma.dtype)
    elapsed_time = time.time() - start_time

    log.debug("time create sparse matrix: {:.5f}".format(elapsed_time))

    if use_triu:
        diagmatrix = dia_matrix(([result.diagonal()], [0]),
                                shape=(M, M), dtype=ma.dtype)

        # to make set the main diagonal to zero I
        # multiply the diagonal by 2
        if diagonal is False:
            diagmatrix *= 2
        # make result symmetric
        matrix = result + result.T - diagmatrix
    else:
        matrix = result
    matrix.eliminate_zeros()

    return matrix
