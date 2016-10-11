import sys
from scipy.sparse import coo_matrix, dia_matrix, triu
import numpy as np
import time
import logging


def reduce_matrix(matrix, bins_to_merge, diagonal=False):
    """
    This function sums the rows and columns corresponding
    to the bins_to_merge, returning a new sparse
    matrix of size len(bins_to_merge).

    The function uses sparse matrices tricks to work very fast.
    The idea is that all the cells in the new reduced matrix
    whose sum needs to be computed are expressed as a vector
    of bins.  Using the bincount function the sum of each bin
    is obtained. Then this sum vector is converted to a
    sparse matrix.


    Parameters
    ----------

    bins_to_merge : A list of lists. The values of the lists should
        correspond to the indices of the matrix.

    diagonal : If set to true, then the main diagonal is preserved (not deleted)

    Returns
    -------

    A sparse matrix.

    >>> from scipy.sparse import csr_matrix
    >>> A = csr_matrix(np.array([[1,0],[0,1]]))
    >>> reduce_matrix(A, [(0,1)], diagonal=True).todense()
    matrix([[2]])
    >>> A = csr_matrix(np.array([[1,0,0],[0,1,0],[0,0,0]]))
    >>> reduce_matrix(A, [(0,1), (2,)], diagonal=True).todense()
    matrix([[2, 0],
            [0, 0]])
    >>> A = csr_matrix(np.array([[12,5,3,2,0],[0,11,4,1,1],
    ... [0,0,9,6,0], [0,0,0,10,0], [0,0,0,0,0]]))
    >>> A.todense()
    matrix([[12,  5,  3,  2,  0],
            [ 0, 11,  4,  1,  1],
            [ 0,  0,  9,  6,  0],
            [ 0,  0,  0, 10,  0],
            [ 0,  0,  0,  0,  0]])
    >>> ll = [(0,2), (1,3), (4,)]
    >>> reduce_matrix(A, ll, diagonal=True).todense()
    matrix([[24, 17,  0],
            [17, 22,  1],
            [ 0,  1,  0]])
    >>> ll = [(0,2,4), (1,3)]
    >>> reduce_matrix(A, ll, diagonal=True).todense()
    matrix([[24, 18],
            [18, 22]])

    Test removal of row/columns when the scaffold list does
    not contains all indices
    >>> ll = [(1,2), (3,)]
    >>> reduce_matrix(A, ll, diagonal=True).todense()
    matrix([[24,  7],
            [ 7, 10]])

    Test with float and nan
    >>> A = csr_matrix(np.array([[12.4,5.3,3.1,2.2,np.nan],
    ... [0,11.4,4.4,1.1,1.1],
    ... [0,0,9,6,0], [0,0,0,10,0], [0,0,0,0,0]]))
    >>> A.todense()
    matrix([[ 12.4,   5.3,   3.1,   2.2,   nan],
            [  0. ,  11.4,   4.4,   1.1,   1.1],
            [  0. ,   0. ,   9. ,   6. ,   0. ],
            [  0. ,   0. ,   0. ,  10. ,   0. ],
            [  0. ,   0. ,   0. ,   0. ,   0. ]])

    >>> ll = [(0,2), (1,3), (4,)]
    >>> reduce_matrix(A, ll, diagonal=True).todense()
    matrix([[ 24.5,  17.9,   nan],
            [ 17.9,  22.5,   1.1],
            [  nan,   1.1,   0. ]])
    """

    #use only the upper triangle of the matrix
    # and convert to coo
    #logging.info("reducing matrix")
    """
    try:
        if sum([len(x) for x in bins_to_merge]) != matrix.shape[0]:
            raise Exception("bins_to_merge length different than "
                        "matrix length")
    except:
        import ipdb;ipdb.set_trace()
    """

    start_time = time.time()
    ma = triu(matrix, k=0, format='coo')
    M = len(bins_to_merge)
    elapsed_time = time.time() - start_time
    start_time = time.time()
    #logging.debug("time triu: {:.5f}".format(elapsed_time))
    # place holder for the new row and col
    # vectors

    if M == ma.shape[0]:
        return matrix

    # check for nans and make a warning in case they are
    # present
    num_nan = len(np.flatnonzero(np.isnan(np.array(ma.data))))
    if num_nan > 0:
        sys.stderr.write(
            "*Warning*\nmatrix contains {} NaN values.\n".format(num_nan))

    # each original col and row index is converted
    # to a new index based on the bins_to_merge.
    # For example, if rows 1 and 10 are to be merged
    # then all rows whose value is 1 or 10 are given
    # as new value the index in the bins_to_merge list.

    map_ = np.zeros(ma.shape[0], dtype=int) -1 # -1 such that all 
                                               # cases not replaced
                                               # by the next loop
                                               # can be identified later.
                                               # Those cases that remain as -1
                                               # are for the rows/cols not
                                               # appearing in the bins_to_merge
    for k, v in enumerate(bins_to_merge):
        for x in v:
            map_[x] = k

    new_row = np.take(map_, ma.row)
    new_col = np.take(map_, ma.col)

    elapsed_time = time.time() - start_time
    start_time = time.time()
    # logging.debug("time mapping: {:.5f}".format(elapsed_time))

    # remove rows and cols with -1. Those
    # correspond to regions that are not part of
    # the bins_to_merge to be merged
    keep = (new_row > -1) & (new_col >-1)
    ma.data = ma.data[keep]
    new_row = new_row[keep]
    new_col = new_col[keep]

    elapsed_time = time.time() - start_time
    start_time = time.time()
    # logging.debug( "time removing unused: {:.5f}".format(elapsed_time))

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
    logging.debug("time complex unique: {:.5f}".format(elapsed_time))
    # sum the bins. For each bin, all values associated
    # to the bin members are added. Those values are in
    # the ma.data array.
    sum_array = np.bincount(bin_array, weights=ma.data)
    elapsed_time = time.time() - start_time
    start_time = time.time()
    # logging.debug( "time bincount: {:.5f}".format(elapsed_time))

    # To reconstruct the new matrix the row and col
    # positions need to be obtained. They correspond
    # to the positions of the unique values in the bin_array.
    uniq, ind = np.unique(bin_array, return_index=True)
    new_row_ = new_row[ind]
    new_col_ = new_col[ind]

    elapsed_time = time.time() - start_time
    start_time = time.time()

    # logging.debug( "time rebuild matrix guts: {:.5f}".format(elapsed_time))
    result = coo_matrix((sum_array, (new_row_, new_col_)),
                        shape=(M,M), dtype=ma.dtype)
    elapsed_time = time.time() - start_time
    start_time = time.time()
    logging.debug( "time create sparse matrix: {:.5f}".format(elapsed_time))

    diagmatrix = dia_matrix(([result.diagonal()], [0]),
                            shape=(M, M), dtype=ma.dtype)

    # to make set the main diagonal to zero I
    # multiply the diagonal by 2
    if diagonal is False:
        diagmatrix *= 2

    elapsed_time = time.time() - start_time
    start_time = time.time()
    #logging.debug( "time make diagonal: {:.5f}".format(elapsed_time))
    # make result symmetric
    matrix = result + result.T - diagmatrix 
    matrix.eliminate_zeros()

    elapsed_time = time.time() - start_time
    start_time = time.time()
    #logging.debug( "time final matrix: {:.5f}".format(elapsed_time))

    return matrix

