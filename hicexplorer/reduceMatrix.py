"""hicexplorer.reduceMatrix, implemented in C++.

reduce_matrix keeps its 3.7.6 signature and results; the body calls
hicexplorer._hicx.reduceMatrix (cpp/python/src/api/reduceMatrix.cpp,
cpp/core/include/hicx/reduce_matrix_entries.hpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
# the doctests below use these names
from scipy.sparse import dia_matrix  # noqa: F401 (the doctests use it)
import numpy as np  # noqa: F401 (the doctests use it)

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import reduceMatrix as _impl


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
    return _impl.reduce_matrix(matrix, bins_to_merge, use_triu, diagonal)
