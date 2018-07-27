from __future__ import division
# from __future__ import unicode_literals
from builtins import range
from past.builtins import zip
from six import iteritems

import numpy as np
# import os
import sys
from collections import OrderedDict
from scipy.sparse import csr_matrix, dia_matrix, triu, tril
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
import tables
from intervaltree import IntervalTree, Interval
from .utilities import toBytes
from .utilities import toString
from .utilities import check_chrom_str_bytes
from .lib import MatrixFileHandler
# import gzip

import cooler

import logging
log = logging.getLogger(__name__)

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings(action="ignore", message="numpy.dtype size changed")
warnings.filterwarnings(action="ignore", message="numpy.ndarray size changed")
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=ImportWarning)
warnings.simplefilter(action='ignore', category=tables.exceptions.FlavorWarning)


# try to import pandas if exists
# try:
#     import pandas as pd
#     pandas = True
# except ImportError:
#     pandas = False


class hiCMatrix:
    """
    Class to handle Hi-C matrices
    contains routines to get intrachromosomal distances
    get sub matrices by chrname.
    """

    def __init__(self, pMatrixFile=None, pChrnameList=None, pCooler_only_init=None,
                 pApplyCorrectionCooler=None):
        # self.correction_factors = None  # this value is set in case a matrix was iteratively corrected
        self.non_homogeneous_warning_already_printed = False
        # self.distance_counts = None  # only defined when getCountsByDistance is called
        self.bin_size = None
        self.bin_size_homogeneous = None  # track if the bins are equally spaced or not
        self.uncorrected_matrix = None

        self.matrix = None
        self.cut_intervals = None
        self.nan_bins = None
        self.correction_factors = None
        self.distance_counts = None
        # # when NaN bins are masked, this variable becomes contains the bin index
        # # needed to put the masked bins back into the matrix.
        self.orig_bin_ids = []
        self.orig_cut_intervals = []  # similar to orig_bin_ids. Used to identify the position of masked nan bins
        self.matrixFileHandler = None

        if pMatrixFile is not None:
            log.debug('Load self.matrixFileHandler')
            fileType = 'cool'
            if pMatrixFile.endswith('.h5'):
                fileType = 'h5'
            self.matrixFileHandler = MatrixFileHandler(pFileType=fileType, pMatrixFile=pMatrixFile, pChrnameList=pChrnameList, pCooler_only_init=pCooler_only_init,
                                                       pApplyCorrectionCooler=pApplyCorrectionCooler)
            self.matrix, self.cut_intervals, self.nan_bins, \
                self.correction_factors, self.distance_counts = self.matrixFileHandler.load()

            # log.info('self.matrix {}'.format(self.matrix))
            # log.info('self.cut_intervals {}'.format(self.cut_intervals))
            # log.info('self.nan_bins {}'.format(self.nan_bins))
            # log.info('self.correction_factors {}'.format(self.correction_factors))
            # log.info('self.distance_counts {}'.format(self.distance_counts))
            # log.info('self.matrix {}'.format(self.matrix))
            if self.nan_bins is None:
                self.nan_bins = np.array([])
            # log.debug('Apply fillLowerTriangle')
            # log.debug('self.matix before: {}'.format(self.matrix))
            # log.debug('self.matix.shape before: {}'.format(self.matrix.shape))

            self.fillLowerTriangle()
            # log.debug('self.matix after: {}'.format(self.matrix))

            self.restoreMaskedBins()
            self.interval_trees, self.chrBinBoundaries = \
                self.intervalListToIntervalTree(self.cut_intervals)
        elif pCooler_only_init:
            self.matrixFileHandler = MatrixFileHandler()
        elif pMatrixFile is None:
            # self.matrixFileHandler = MatrixFileHandler()

            log.info('Only init object, no matrix given.')
        else:
            log.error('matrix file not given')
            sys.exit(1)

    def save(self, pMatrixName, pSymmetric=True, pApplyCorrection=False):
        """ As an output format cooler and mcooler are supported.
        """
        if self.matrixFileHandler is None:
            fileType = 'cool'
            if pMatrixName.endswith('h5'):
                fileType = 'h5'
            self.matrixFileHandler = MatrixFileHandler(pFileType=fileType)

        self.restoreMaskedBins()
        self.matrixFileHandler.set_matrix_variables(self.matrix, self.cut_intervals, self.nan_bins,
                                                    self.correction_factors, self.distance_counts)
        if pMatrixName.endswith('cool') or pMatrixName.endswith('h5'):
            self.matrixFileHandler.save(pMatrixName, pSymmetric=pSymmetric, pApplyCorrection=pApplyCorrection)

    def fillLowerTriangle(self):
        """
        checks if the matrix is complete or if only half of the matrix was saved.
        Returns a whole matrix.

        >>> from scipy.sparse import csr_matrix
        >>> import numpy as np
        >>> A = csr_matrix(np.array([[12,5,3,2,0],[0,11,4,1,1],
        ... [0,0,9,6,0], [0,0,0,10,0], [0,0,0,0,0]]), dtype=np.int32)
        >>> hm = hiCMatrix()
        >>> hm.matrix = A
        >>> hm.fillLowerTriangle()
        >>> B = hm.matrix
        >>> B.todense()
        matrix([[12,  5,  3,  2,  0],
                [ 5, 11,  4,  1,  1],
                [ 3,  4,  9,  6,  0],
                [ 2,  1,  6, 10,  0],
                [ 0,  1,  0,  0,  0]], dtype=int32)

        """
        # log.debug('sum of tril: {}'.format(tril(self.matrix, k=-1).sum()))
        if tril(self.matrix, k=-1).sum() == 0:
            # this case means that the lower triangle of the
            # symmetric matrix (below the main diagonal)
            # is zero. In this case, replace the lower
            # triangle using the upper triangle
            self.matrix = self.matrix + triu(self.matrix, 1).T

        # return matrix

    def setCutIntervals(self, cut_intervals):
        """
        Replace the cut_intervals of a matrix
        """

        # check that the matrix is squared
        if len(cut_intervals) != self.matrix.shape[0]:
            raise Exception("Length of cut_intervals {} does not match the "
                            "matrix size {}".format(len(cut_intervals), self.matrix.shape))

        self.cut_intervals = cut_intervals
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def setMatrix(self, matrix, cut_intervals):
        """
        Initialize a matrix with a given matrix
        and cut_intervals. Mostly useful for
        testing.
        """

        # check that the matrix is squared
        if matrix.shape[0] != matrix.shape[1]:
            raise Exception("Matrix is not squared. Shape is {}".format(matrix.shape))
        if len(cut_intervals) != matrix.shape[0]:
            raise Exception("Length of cut_intervals {} does not match the matrix size {}".format(len(cut_intervals),
                                                                                                  matrix.shape))

        self.matrix = matrix
        self.cut_intervals = cut_intervals
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def getBinSize(self):
        """
        estimates the bin size. In case the bin size
        is not equal for all bins (maybe except for the
        bin at the en of the chromosomes) a warning is issued.
        In case of uneven bins, the median is returned.
        """
        if self.bin_size is None:
            chrom, start, end, extra = zip(*self.cut_intervals)
            median = int(np.median(np.diff(start)))
            diff = np.array(end) - np.array(start)

            # check if the bin size is
            # homogeneous
            if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
                self.bin_size_homogeneous = False
                if self.non_homogeneous_warning_already_printed is False:
                    log.warning('Bin size is not homogeneous. \
                                      Median {}\n'.format(median))
                    self.non_homogeneous_warning_already_printed = True
            self.bin_size = median
        return self.bin_size

    def getMatrix(self):
        matrix = self.matrix.todense()
        if len(self.nan_bins):
            # to set NaN values the matrix type has to be
            # float. Corrected matrices are of float
            # type while uncorrected matrices are of
            # of int type
            if np.issubdtype(self.matrix, 'float') is False:
                matrix = matrix.astype(float)
            matrix[self.nan_bins, :] = np.nan
            matrix[:, self.nan_bins] = np.nan

        return matrix

    def getChrBinRange(self, chrName):
        """
        Given a chromosome name,
        This functions return the start and end bin indices in the matrix
        """

        return self.chrBinBoundaries[chrName]

    def getChrNames(self):
        """
        returns the names of the chromosomes
        present in the matrix
        """
        return list(self.chrBinBoundaries)

    def getBinPos(self, binIndex):
        """
        given a bin, it returns the chromosome name,
        start position and end position
        """
        return self.cut_intervals[binIndex]

    def getRegionBinRange(self, chrname, startpos, endpos):
        """
        Given a chromosome region, this function returns
        the bin indices that overlap with such region.
        """

        try:
            # chromosome_size = hic_matrix.get_chromosome_sizes()
            # chrname = check_chrom_str_bytes(self.interval_trees, chrname)
            if type(next(iter(self.interval_trees))) != type(chrname):
                if type(next(iter(self.interval_trees))) is str:
                    chrname = toString(chrname)
                elif type(next(iter(self.interval_trees))) is bytes:
                    chrname = toBytes(chrname)
                elif type(next(iter(self.interval_trees))) is np.bytes_:
                    chrname = toBytes(chrname)
            # chr_end_pos = chromosome_size[chrname]
            self.interval_trees[chrname]
        except KeyError:

            log.exception("chromosome: {} name not found in matrix".format(chrname))
            log.exception("valid names are:")
            log.exception(self.interval_trees.keys())
            exit(1)
        try:
            startpos = int(startpos)
            endpos = int(endpos)
        except ValueError:
            log.exeption("{} or {}  are not valid "
                         "position values.".format(startpos, endpos))
            exit(1)

        try:
            startbin = sorted(self.interval_trees[chrname][startpos:startpos + 1])[0].data
            endbin = sorted(self.interval_trees[chrname][endpos:endpos + 1])[0].data
        except IndexError:
            log.exception("Index error")
            return None

        return startbin, endbin

    @staticmethod
    def getDistList(rows, cols, cut_intervals):
        """
            Given a list of rows and cols
            an array is returned containing
            the genomic distance between
            each element of the row array
            with each element of the col array.
            -1 is returned for inter-chromosomal
            interactions.

            A matching list containing the chromosome name
            is also returned

        >>> from scipy.sparse import coo_matrix
        >>> import numpy as np
        >>> row, col = np.triu_indices(5)
        >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
        ... ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
        >>> dist_list, chrom_list = hiCMatrix.getDistList(row, col,
        ... cut_intervals)
        >>> coo_matrix((dist_list, (row, col)), shape=(5,5), dtype=np.int32).todense()
        matrix([[ 0, 10, 20, 30, -1],
                [ 0,  0, 10, 20, -1],
                [ 0,  0,  0, 10, -1],
                [ 0,  0,  0,  0, -1],
                [ 0,  0,  0,  0,  0]], dtype=int32)
        >>> chrom_list.tolist()
        ['a', 'a', 'a', 'a', '', 'a', 'a', 'a', '', 'a', 'a', '', 'a', '', 'b']
        """
        chrnamelist, startlist, endlist, extralist = zip(*cut_intervals)
        # now the distance between any two points
        # is computed and arranged such that for each
        # element of the data array, a corespondent distance is stored
        start_row = np.take(startlist, rows)
        start_col = np.take(startlist, cols)
        dist_list = start_col - start_row

        # now  all distances that are between chromosomes are removed
        # to do this I convert the array of chromosomes to
        # a array of indices. Then, when subtracting the
        # values that correspond to matrix.row and matrix.col
        # using the array of indices, any value other
        # than 0 means inter-chromosomal row,col combination.

        # chr_id_list is based on a trick using np.unique
        # to get from a list of strings
        # a list of integers
        chr_id_list = np.unique(chrnamelist, return_inverse=True)[1]

        chr_row = np.take(chr_id_list, rows)
        chr_col = np.take(chr_id_list, cols)
        chr_diff = chr_row - chr_col
        # set in dist_list array '-1' for all interchromosomal values
        dist_list[chr_diff != 0] = -1

        # make a corresponding chromosome name list
        # if filtering per chromosome is required
        chrom_list = np.take(chrnamelist, rows)
        chrom_list[chr_diff != 0] = ''

        return dist_list, chrom_list

    @staticmethod
    def fit_cut_intervals(cut_intervals):
        # check that the matrix has bins of same size
        # otherwise try to adjust the bins to
        # to match a regular binning
        if len(cut_intervals) <= 1:
            # do nothing if there is only one interval
            return cut_intervals
        chrom, start, end, extra = zip(*cut_intervals)

        median = int(np.median(np.diff(start)))
        diff = np.array(end) - np.array(start)
        # check if the bin size is homogeneous
        if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
            # set the start position of a bin to the closest multiple
            # of the median
            def snap_nearest_multiple(start_x, m):
                resi = [-1 * (start_x % m), -start_x % m]
                return start_x + resi[np.argmin(np.abs(resi))]
            start = [snap_nearest_multiple(x, median) for x in start]
            end = [snap_nearest_multiple(x, median) for x in end]
            cut_intervals = zip(chrom, start, end, extra)
            log.info('[getCountsByDistance] Bin size is not '
                     'homogeneous, setting \n'
                     'the bin distance to the median: {}\n'.format(median))
        return cut_intervals

    def convert_to_zscore_matrix(self, maxdepth=None, perchr=False):
        return self.convert_to_obs_exp_matrix(maxdepth=maxdepth, zscore=True, perchr=perchr)

    def convert_to_obs_exp_matrix(self, maxdepth=None, zscore=False, perchr=False):
        """
        Converts a corrected counts matrix into a
        obs / expected matrix or z-scores fast.

        The caveat is that the obs/exp or z-score are only
        computed for non-zero values, although zero values that
        are not part of the sparse matrix are considered.

        For each diagonal the mean (and std when computing z-scores) are
        calculated and then each non-zero value of the sparse matrix is
        replaced by the obs/exp or z-score.

        Parameters
        ----------
        maxdepth: maximum distance from the diagonal to consider. All contacts beyond this distance will not
                         be considered.
        zscore: if a zscore wants to be returned instead of obs/exp


        Returns
        -------
        observed / expected sparse matrix

        >>> from scipy.sparse import csr_matrix, dia_matrix
        >>> import numpy as np
        >>> row, col = np.triu_indices(5)
        >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
        ... ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
        >>> hic = hiCMatrix()
        >>> hic.nan_bins = []
        >>> matrix = np.array([
        ... [ 1,  8,  5, 3, 0],
        ... [ 0,  4, 15, 5, 1],
        ... [ 0,  0,  0, 7, 2],
        ... [ 0,  0,  0, 0, 1],
        ... [ 0,  0,  0, 0, 0]])

        >>> hic.matrix = csr_matrix(matrix)
        >>> hic.setMatrix(hic.matrix, cut_intervals)
        >>> hic.convert_to_obs_exp_matrix().todense()
        matrix([[1. , 0.8, 1. , 1. , 0. ],
                [0. , 4. , 1.5, 1. , 1. ],
                [0. , 0. , 0. , 0.7, 2. ],
                [0. , 0. , 0. , 0. , 1. ],
                [0. , 0. , 0. , 0. , 0. ]])

        >>> hic.matrix = csr_matrix(matrix)
        >>> hic.convert_to_obs_exp_matrix(maxdepth=20).todense()
        matrix([[1. , 0.8, 1. , 0. , 0. ],
                [0. , 4. , 1.5, 1. , 0. ],
                [0. , 0. , 0. , 0.7, nan],
                [0. , 0. , 0. , 0. , nan],
                [0. , 0. , 0. , 0. , 0. ]])

        >>> hic.matrix = csr_matrix(matrix)
        >>> hic.convert_to_obs_exp_matrix(zscore=True).todense()
        matrix([[ 0.        , -0.56195149,         nan,         nan, -1.41421356],
                [ 0.        ,  1.93649167,  1.40487872,         nan,  0.        ],
                [ 0.        ,  0.        , -0.64549722, -0.84292723,  1.41421356],
                [ 0.        ,  0.        ,  0.        , -0.64549722,  0.        ],
                [ 0.        ,  0.        ,  0.        ,  0.        , -0.64549722]])

        nans occur where the standard deviation is zero
        """

        binsize = self.getBinSize()
        max_depth_in_bins = None

        if maxdepth:
            if maxdepth < binsize:
                exit("Please specify a maxDepth larger than bin size ({})".format(binsize))

            max_depth_in_bins = int(float(maxdepth * 1.5) / binsize)
            # work only with the upper matrix
            # and remove all pixels that are beyond
            # max_depth_in_bis
            # (this is done by subtracting a second sparse matrix
            # that contains only the upper matrix that wants to be removed.
            self.matrix = triu(self.matrix, k=0, format='csr') - \
                triu(self.matrix, k=max_depth_in_bins, format='csr')
        else:
            self.matrix = triu(self.matrix, k=0, format='csr')

        self.matrix.eliminate_zeros()
        depth = None
        if zscore is True:
            from scipy.sparse import diags
            m_size = self.matrix.shape[0]
            if max_depth_in_bins is not None:
                depth = max_depth_in_bins
            else:
                depth = m_size
                estimated_size_dense_matrix = m_size ** 2 * 8
                if estimated_size_dense_matrix > 100e6:
                    log.info("To compute z-scores a dense matrix is required. This will use \n"
                             "{} Mb of memory.\n To reduce memory use the maxdeph option."
                             "".format(estimated_size_dense_matrix / 1e6))

            # to compute zscore the zero values need to be accounted and the matrix
            # need to become dense. This is only practical if only up to certain distance
            # wants to be evaluated, otherwise the dense matrix is too large.
            # To make the matrix dense and keep the same computations as when
            # the matrix is sparse the following is done:
            # A sparse diagonal matrix of shape = matrix.shape is created with ones
            # (only upper triangle contains diagonals up to maxdeph)
            # This  sparse matrix is then added to self.matrix
            # then, -1 is subtracted to the self.matrix.data, thus effectively
            # adding zeros.
            diag_mat_ones = diags(np.repeat([1], m_size * depth).reshape(depth, m_size), list(range(depth)))

            self.matrix += diag_mat_ones

        from scipy.sparse import lil_matrix
        trasf_matrix = lil_matrix(self.matrix.shape)

        chr_submatrix = OrderedDict()
        cut_intervals = OrderedDict()
        chrom_sizes = OrderedDict()
        chrom_range = OrderedDict()
        if perchr:
            for chrname in self.getChrNames():
                chr_range = self.getChrBinRange(chrname)
                chr_submatrix[chrname] = self.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]].tocoo()
                cut_intervals[chrname] = [self.cut_intervals[x] for x in range(chr_range[0], chr_range[1])]
                chrom_sizes[chrname] = [chr_submatrix[chrname].shape[0]]
                chrom_range[chrname] = (chr_range[0], chr_range[1])

        else:
            chr_submatrix['all'] = self.matrix.tocoo()
            cut_intervals['all'] = self.cut_intervals
            chrom_sizes['all'] = np.array([v[1] - v[0] for k, v in iteritems(self.chrBinBoundaries)])
            chrom_range['all'] = (0, self.matrix.shape[0])

        for chrname, submatrix in iteritems(chr_submatrix):
            log.info("processing chromosome {}\n".format(chrname))
            if zscore is True:
                # this step has to be done after tocoo()
                submatrix.data -= 1

            dist_list, chrom_list = self.getDistList(submatrix.row, submatrix.col,
                                                     hiCMatrix.fit_cut_intervals(cut_intervals[chrname]))

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
            sum_counts = np.bincount(dist_list, weights=submatrix.data)
            distance_len = np.bincount(dist_list)
            # compute the average for each distance
            mat_size = submatrix.shape[0]
            mu = {}
            std = {}
            # compute mean value for each distance

            for bin_dist_plus_one, sum_value in enumerate(sum_counts):
                if maxdepth and bin_dist_plus_one == 0:  # this is for intra chromosomal counts
                    # when max depth is set, the computation
                    # of the total_intra is not accurate and is safer to
                    # output np.nan
                    mu[bin_dist_plus_one] = np.nan
                    std[bin_dist_plus_one] = np.nan
                    continue

                if bin_dist_plus_one == 0:
                    total_intra = mat_size ** 2 - sum([size ** 2 for size in chrom_sizes[chrname]])
                    diagonal_length = int(total_intra / 2)
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
                    diagonal_length = sum([size - (bin_dist_plus_one - 1) for size in chrom_sizes[chrname] if size > (bin_dist_plus_one - 1)])
                    log.debug("Type of diagonal_length {}".format(type(diagonal_length)))

                # the diagonal length should contain the number of values at a certain distance.
                # If the matrix is dense, the distance_len[bin_dist_plus_one] correctly contains the number of values
                # If the matrix is equally spaced, then, the diagonal_length as computed before is accurate.
                # But, if the matrix is both sparse and with unequal bins, then none of the above methods is
                # accurate but the the diagonal_length as computed before will be closer.
                diagonal_length = max(diagonal_length, distance_len[bin_dist_plus_one])
                log.debug("Type of diagonal_length {}".format(type(diagonal_length)))

                if diagonal_length == 0:
                    mu[bin_dist_plus_one] = np.nan
                else:
                    mu[bin_dist_plus_one] = np.float64(sum_value) / diagonal_length

                if np.isnan(sum_value):
                    log.info("nan value found for distance {}\n".format((bin_dist_plus_one - 1) * binsize))

                # if zscore is needed, compute standard deviation: std = sqrt(mean(abs(x - x.mean())**2))
                if zscore:
                    values_sqrt_diff = \
                        np.abs((submatrix.data[dist_list == bin_dist_plus_one] - mu[bin_dist_plus_one]) ** 2)
                    # the standard deviation is the sum of the differences with mu squared (value variable)
                    # plus all zeros that are not included in the sparse matrix
                    # for which the standard deviation is
                    # (0 - mu)**2 = (mu)**2
                    # The number of zeros is the diagonal length - the length of the non zero values
                    zero_values_sqrt_diff_sum = (diagonal_length - len(values_sqrt_diff)) * mu[bin_dist_plus_one] ** 2

                    _std = np.sqrt((values_sqrt_diff.sum() + zero_values_sqrt_diff_sum) / diagonal_length)
                    std[bin_dist_plus_one] = _std

            # use the expected values to compute obs/exp
            transf_ma = np.zeros(len(submatrix.data))
            for idx, value in enumerate(submatrix.data):
                if depth is not None and dist_list[idx] > depth + 1:
                    continue
                if zscore:
                    if std[dist_list[idx]] == 0:
                        transf_ma[idx] = np.nan
                    else:
                        transf_ma[idx] = (value - mu[dist_list[idx]]) / std[dist_list[idx]]
                else:
                    transf_ma[idx] = value / mu[dist_list[idx]]

            submatrix.data = transf_ma
            trasf_matrix[chrom_range[chrname][0]:chrom_range[chrname][1], chrom_range[chrname][0]:chrom_range[chrname][1]] = submatrix.tolil()

        self.matrix = trasf_matrix.tocsr()

        return self.matrix

    @staticmethod
    def dist_list_to_dict(data, dist_list):
        """
        splits data, into numeric groups defined by dist_list
        Return a dictionary containing, for
        each unique distance a dictionary
        """

        order = np.argsort(dist_list)
        dist_list = dist_list[order]
        data = data[order]

        # having the dist_list sorted, np.split
        # is used to divide the data into
        # groups that lie at the same distance, for this
        # np.diff together with np.flatnonzero is used to
        # find the indices where the distance changes.
        # the '+1' is needed because the np.diff array is
        # one element smaller than the original array, thus
        # the indices based no the np.diff array are off by 1
        # with respect to the original array
        groups = np.split(data, np.flatnonzero(np.diff(dist_list)) + 1)

        # because the dist_list is sorted
        # the order of the unique values
        # corresponds to that of the groups.
        # In other words, group[0]
        # has distance_unique[0]
        # np.sort after np.unique  in theory
        # is not needed, but just in case...
        distance_unique = np.sort(np.unique(dist_list))

        # convert to dictionary having as key
        # the distance
        distance = {}
        for index in range(len(distance_unique)):
            distance[distance_unique[index]] = groups[index]

        return distance

    def keepOnlyTheseChr(self, chromosome_list):
        """
        given a list of chromosome names,
        these are kept, while any other is removed
        from the matrix
        """
        chromosome_list = check_chrom_str_bytes(self.interval_trees, chromosome_list)

        try:
            [self.chrBinBoundaries[x] for x in chromosome_list]
        except KeyError as e:
            raise Exception("Chromosome name {} not in matrix.".format(e))

        self.restoreMaskedBins()
        size = self.matrix.shape
        # initialize a 1D array containing the columns (and rows) to
        # select. By default none are selected
        sel = np.empty(size[0], dtype=np.bool)
        sel[:] = False

        for chrName in list(self.interval_trees):
            if chrName not in chromosome_list:
                continue

            # identify start and end rows
            # of chromosomes that wants to be
            # kept
            index_start, index_end = self.getChrBinRange(chrName)
            sel[index_start:index_end] = True

        sel_id = np.flatnonzero(sel)
        mat = self.matrix[sel_id, :][:, sel_id]

        # update bin ids
        self.cut_intervals = [self.cut_intervals[x] for x in sel_id]

        # update correction factors
        if self.correction_factors is not None:
            self.correction_factors = [self.correction_factors[x] for x in sel_id]

        # keep track of nan bins
        if len(self.nan_bins):
            _temp = np.zeros(size[0])
            _temp[self.nan_bins] = 1
            _temp = _temp[sel_id]
            self.nan_bins = np.flatnonzero(_temp == 1)
        else:
            self.nan_bins = []

        self.numCols = len(sel_id)

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)
        # remove distanceCounts
        try:
            self.distance_counts = None
        except AttributeError:
            pass
        self.matrix = mat
        return self.matrix

# <<<<<<< HEAD
# =======
#     def save_bing_ren(self, fileName):
#         """
#         Saves the matrix using bing ren's
#         method which is chrom_name\tstart_bin\tend_bin\tvalues...
#         """
#         if os.path.isdir(fileName):
#             exit('Please specify a file name to save the data. The given file name is a folder: \n{}\n'.format(fileName))

#         if fileName[-3:] != '.gz':
#             fileName += '.gz'

#         try:
#             fileh = gzip.open(fileName, 'wt')
#         except IOError:
#             msg = "{} file can't be opened for writing".format(fileName)
#             raise msg

#         colNames = []
#         for x in range(self.matrix.shape[0]):
#             chrom, start, end = self.cut_intervals[x][0:3]
#             colNames.append("{}\t{}\t{}".format(chrom, start, end))

#         for row in range(self.matrix.shape[0]):
#             values = [str(x) for x in self.matrix[row, :].toarray().flatten()]
#             fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values)))

#         fileh.close()

#     def save_dekker(self, fileName):
#         """
#         Saves the matrix using dekker format
#         """
#         if os.path.isdir(fileName):
#             exit('Please specify a file name to save the data. The given file name is a folder: \n{}\n'.format(fileName))

#         if fileName[-3:] != '.gz':
#             fileName += '.gz'

#         try:
#             fileh = gzip.open(fileName, 'wt')
#         except IOError:
#             msg = "{} file can't be opened for writing".format(fileName)
#             raise msg

#         colNames = []
#         for x in range(self.matrix.shape[0]):
#             chrom, start, end = self.cut_intervals[x][0:3]
#             colNames.append("{}|--|{}:{}-{}".format(x, toString(chrom), start, end))  # adds dm3 to the end (?problem..)

#         fileh.write("#converted from hicexplorer\n")
#         fileh.write("\t" + "\t".join(colNames) + "\n")
#         for row in range(self.matrix.shape[0]):
#             values = [str(x) for x in self.matrix[row, :].toarray().flatten()]
#             fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values)))

#         fileh.close()

#     def save_lieberman(self, fileName):
#         """
#         Saves the matrix using lieberman format. Given an output directory name and resolution of the matrix.

#         In contrast to other methods, when saving using liebermans format a folder is required
#         where the data is saved per chromosome
#         """

#         if os.path.isfile(fileName):
#             exit("a file with the same name as the given path exists ({}). Please choose a folder.".format(fileName))

#         if not os.path.isdir(fileName):
#             os.mkdir(fileName)

#         resolution = self.getBinSize()
#         if not self.bin_size_homogeneous:
#             log.warning("*WARNING* The hic matrix contains bins of difference size but\n"
#                         "the 'lieberman' format requires equally spaced bins. The program\n"
#                         "will proceed but the results may be unreliable.\n")

#         for chrom in list(self.interval_trees):
#             chrstart, chrend = self.getChrBinRange(chrom)
#             chrwise_mat = self.matrix[chrstart:chrend, chrstart:chrend]
#             if len(chrwise_mat.data) > 0:
#                 log.info("Saving chromosome {}...\n".format(chrom))
#                 chrwise_mat_coo = triu(chrwise_mat, k=0, format='csr').tocoo()
#                 start = chrwise_mat_coo.row * resolution
#                 end = chrwise_mat_coo.col * resolution
#                 fileh = gzip.open("{}/chr{}.gz".format(fileName, chrom), 'wt')
#                 fileh.write("#converted from HiCExplorer format\n")
#                 for idx in range(len(start)):
#                     fileh.write("{}\t{}\t{}\n".format(start[idx], end[idx], chrwise_mat_coo.data[idx]))
#                 fileh.close()

#     def save_GInteractions(self, fileName):
#         """
#         Saves the matrix using bioconductor's GInteraction format. `bin_pos1 , bin_pos2, number of interactions`
#         """
#         self.restoreMaskedBins()
#         log.debug(self.matrix.shape)
#         mat_coo = triu(self.matrix, k=0, format='csr').tocoo()
#         fileh = open("{}.tsv".format(fileName), 'w')
#         for idx, counts in enumerate(mat_coo.data):
#             chr_row, start_row, end_row, _ = self.cut_intervals[mat_coo.row[idx]]
#             chr_col, start_col, end_col, _ = self.cut_intervals[mat_coo.col[idx]]
#             fileh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_row, int(start_row), int(end_row),
#                                                               chr_col, int(start_col), int(end_col), counts))
#         fileh.close()

#     def create_empty_cool_file(self, pFileName):
#         bins_data_frame = pd.DataFrame(columns=['chrom', 'start', 'end', 'weight'])
#         matrix_data_frame = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'count'])
#         cooler.io.create(cool_uri=pFileName,
#                          bins=bins_data_frame,
#                          pixels=matrix_data_frame)

#     def save_cooler(self, pFileName, pSymmetric=True, pApplyCorrection=True):
#         # for value in self.nan_bins:
#         self.restoreMaskedBins()
#         self.matrix = self.matrix.tolil()
#         if self.nan_bins is not None:
#             self.matrix[self.nan_bins, :] = 0
#             self.matrix[:, self.nan_bins] = 0
#         self.matrix = self.matrix.tocsr()

#         for i in range(len(self.matrix.data)):
#             if np.isnan(self.matrix.data[i]):
#                 self.matrix.data[i] = 0

#         self.matrix.eliminate_zeros()
#         # save only the upper triangle of the
#         if pSymmetric:
#             # symmetric matrix
#             matrix = triu(self.matrix, k=0, format='csr')
#         else:
#             matrix = self.matrix

#         cut_intervals_ = []
#         for value in self.cut_intervals:
#             cut_intervals_.append(tuple((value[0], value[1], value[2])))

#         bins_data_frame = pd.DataFrame(cut_intervals_, columns=['chrom', 'start', 'end'])

#         # append correction factors if they exist
#         if self.correction_factors is not None and pApplyCorrection:
#             weight = convertNansToOnes(np.array(self.correction_factors).flatten())
#             bins_data_frame = bins_data_frame.assign(weight=weight)

#         # get only the upper triangle of the matrix to save to disk
#         # upper_triangle = triu(self.matrix, k=0, format='csr')
#         # create a tuple list and use it to create a data frame

#         # save correction factors and original matrix

#         # revert correction to store orginal matrix
#         if self.correction_factors is not None and pApplyCorrection:
#             log.info("Reverting correction factors on matrix...")
#             instances, features = matrix.nonzero()
#             self.correction_factors = np.array(self.correction_factors)

#             # do not apply if correction factors are just 1's
#             if np.sum(self.correction_factors) != len(self.correction_factors):
#                 instances_factors = self.correction_factors[instances]
#                 features_factors = self.correction_factors[features]

#                 instances_factors *= features_factors
#                 matrix.data *= instances_factors
#                 instances_factors = None
#                 features_factors = None

#                 matrix.data = np.rint(matrix.data)
#                 matrix.data = matrix.data.astype(int)

#             data = matrix.data.tolist()

#         else:

#             instances, features = matrix.nonzero()
#             data = matrix.data.tolist()

#             if matrix.dtype not in [np.int32, int]:
#                 log.warning("Writing non-standard cooler matrix. Datatype of matrix['count'] is: {}".format(matrix.dtype))
#                 cooler._writer.COUNT_DTYPE = matrix.dtype

#         if len(instances) == 0 and len(features) == 0:
#             exit('No data present. Exit.')
#         else:
#             matrix_tuple_list = zip(instances.tolist(), features.tolist(), data)
#             matrix_data_frame = pd.DataFrame(matrix_tuple_list, columns=['bin1_id', 'bin2_id', 'count'])

#             cooler.io.create(cool_uri=pFileName,
#                              bins=bins_data_frame,
#                              pixels=matrix_data_frame,
#                              append=False)

#     def save_hdf5(self, filename, pSymmetric=True):
#         """
#         Saves a matrix using hdf5 format
#         :param filename:
#         :return: None
#         """
#         self.restoreMaskedBins()
#         if not filename.endswith(".h5"):
#             filename += ".h5"

#         # if the file name already exists
#         # try to find a new suitable name
#         if os.path.isfile(filename):
#             log.warning("*WARNING* File already exists {}\n "
#                         "Overwriting ...\n".format(filename))

#             from os import unlink
#             unlink(filename)
#         try:
#             nan_bins = np.array(self.nan_bins)
#         except Exception:
#             nan_bins = np.array([])
#         # save only the upper triangle of the
#         if pSymmetric:
#             # symmetric matrix
#             matrix = triu(self.matrix, k=0, format='csr')
#         else:
#             matrix = self.matrix
#         filters = tables.Filters(complevel=5, complib='blosc')
#         with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
#             matrix_group = h5file.create_group("/", "matrix", )
#             # save the parts of the csr matrix
#             for matrix_part in ('data', 'indices', 'indptr', 'shape'):
#                 arr = np.array(getattr(matrix, matrix_part))
#                 atom = tables.Atom.from_dtype(arr.dtype)
#                 ds = h5file.create_carray(matrix_group, matrix_part, atom,
#                                           shape=arr.shape,
#                                           filters=filters)
#                 ds[:] = arr

#             # save the matrix intervals
#             intervals_group = h5file.create_group("/", "intervals", )
#             chr_list, start_list, end_list, extra_list = zip(*self.cut_intervals)
#             for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
#                 arr = np.array(eval(interval_part))
#                 atom = tables.Atom.from_dtype(arr.dtype)
#                 ds = h5file.create_carray(intervals_group, interval_part, atom,
#                                           shape=arr.shape,
#                                           filters=filters)
#                 ds[:] = arr

#             # save nan bins
#             if len(nan_bins):
#                 atom = tables.Atom.from_dtype(nan_bins.dtype)
#                 ds = h5file.create_carray(h5file.root, 'nan_bins', atom,
#                                           shape=nan_bins.shape,
#                                           filters=filters)
#                 ds[:] = nan_bins

#             # save corrections factors
#             if self.correction_factors is not None and len(self.correction_factors):
#                 self.correction_factors = np.array(self.correction_factors)
#                 atom = tables.Atom.from_dtype(self.correction_factors.dtype)
#                 ds = h5file.create_carray(h5file.root, 'correction_factors', atom,
#                                           shape=self.correction_factors.shape,
#                                           filters=filters)
#                 ds[:] = np.array(self.correction_factors)

#             # save distance counts
#             if self.distance_counts is not None and len(self.distance_counts):
#                 atom = tables.Atom.from_dtype(self.distance_counts.dtype)
#                 ds = h5file.create_carray(h5file.root, 'distance_counts', atom,
#                                           shape=self.distance_counts.shape,
#                                           filters=filters)
#                 ds[:] = np.array(self.distance_counts)

#     def save_npz(self, filename):
#         """
#         saves using the numpy npz format. Adequate for small samples.
#         """
#         self.restoreMaskedBins()
#         chrNameList, startList, endList, extraList = zip(*self.cut_intervals)
#         try:
#             nan_bins = self.nan_bins
#         except Exception:
#             nan_bins = np.array([])
#         # save only the upper triangle of the
#         # symmetric matrix
#         matrix = triu(self.matrix, k=0, format='csr')
#         try:
#             np.savez(
#                 filename, matrix=matrix, chrNameList=chrNameList,
#                 startList=startList, endList=endList, extraList=extraList,
#                 nan_bins=nan_bins, correction_factors=self.correction_factors)
#         except Exception as e:
#             log.debug("error saving matrix: {}".format(e))
#             try:
#                 log.info("Matrix can not be saved because is too big!")
#                 log.version_info("Eliminating entries with only one count.")

#                 # try to remove noise by deleting 1
#                 matrix.data = matrix.data - 1
#                 matrix.eliminate_zeros()
#                 np.savez(
#                     filename, matrix=matrix, chrNameList=chrNameList,
#                     startList=startList, endList=endList, extraList=extraList,
#                     nan_bins=nan_bins)
#             except Exception:
#                 log.exception("Try failed. Matrix can not be saved because it is too big!")
#             exit()

#     def save(self, pMatrixName, pSymmetric=True, pApplyCorrection=False):
#         """To save please specifiy the ending of your format i.e. 'output_matrix.format' Supported are: 'dekker', 'ren',
#             'lieberman', 'npz', 'GInteractions', 'h5' and 'cool'.
#         """
#         if pMatrixName.endswith('dekker') == 'dekker':
#             self.save_dekker(pMatrixName)
#         elif pMatrixName.endswith('ren'):
#             self.save_bing_ren(pMatrixName)
#         elif pMatrixName.endswith('lieberman'):
#             self.save_lieberman(pMatrixName)
#         elif pMatrixName.endswith('npz'):
#             self.save_npz(pMatrixName)
#         elif pMatrixName.endswith('GInteractions'):
#             self.save_GInteractions(pMatrixName)
#         elif pMatrixName.endswith('cool'):
#             self.save_cooler(pMatrixName, pSymmetric=pSymmetric, pApplyCorrection=pApplyCorrection)
#         else:
#             self.save_hdf5(pMatrixName, pSymmetric=pSymmetric)

# >>>>>>> develop
    def diagflat(self, value=np.nan):
        """
        sets
        the matrix diagonal to np.nan
        """
        M = self.matrix.shape[0]
        diagmatrix = dia_matrix((np.repeat(value, M), 0), shape=(M, M))
        self_diag = dia_matrix(([self.matrix.diagonal()], [0]), shape=(M, M))
        # take matrix, subtract the values of the diagonal such that
        # it becomes all zeros, replace with new values by adding them
        self.matrix = self.matrix - self_diag + diagmatrix
        return self.matrix

    def filterOutInterChrCounts(self):
        """
        set all inter chromosomal counts to np.nan
        >>> from scipy.sparse import coo_matrix
        >>> import numpy as np
        >>> row, col = np.triu_indices(5)
        >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
        ... ('a', 20, 30, 1), ('b', 30, 40, 1), ('b', 40, 50, 1)]
        >>> hic = hiCMatrix()
        >>> hic.nan_bins = []
        >>> matrix = np.array([
        ... [ 0, 10,  5, 3, 0],
        ... [ 0,  0, 15, 5, 1],
        ... [ 0,  0,  0, 7, 3],
        ... [ 0,  0,  0, 0, 1],
        ... [ 0,  0,  0, 0, 0]])

        make the matrix symmetric:
        >>> hic.matrix = csr_matrix(matrix + matrix.T)
        >>> hic.setMatrix(csr_matrix(matrix + matrix.T, dtype=np.int32), cut_intervals)
        >>> hic.filterOutInterChrCounts().todense()
        matrix([[ 0, 10,  5,  0,  0],
                [10,  0, 15,  0,  0],
                [ 5, 15,  0,  0,  0],
                [ 0,  0,  0,  0,  1],
                [ 0,  0,  0,  1,  0]], dtype=int32)
        """

        ma_coo = self.matrix.tocoo()
        dist_list, _ = hiCMatrix.getDistList(ma_coo.row, ma_coo.col,
                                             self.cut_intervals)

        # set to zero all cases in which dist_list is zero
        ma_coo.data[dist_list == -1] = 0

        self.matrix = ma_coo.tocsr()
        self.matrix.eliminate_zeros()
        return self.matrix

    def setMatrixValues(self, newMatrix):
        """
        replace the current matrix values
        by the given matrix values. The
        shapes have to coincide
        """
        assert self.matrix.shape == newMatrix.shape,\
            "Given matrix has different shape. New " \
            "values need to have the same shape as previous matrix."

        self.matrix = csr_matrix(newMatrix)

    def setCorrectionFactors(self, correction_factors):
        assert len(correction_factors) == self.matrix.shape[0], \
            "length of correction factors and length of matrix are different."
        self.correction_factors = correction_factors

    def reorderChromosomes(self, new_chr_order):
        new_order = []
        new_chr_order = check_chrom_str_bytes(self.chrBinBoundaries, new_chr_order)

        for chrName in new_chr_order:
            # check that the chromosome names are valid
            if chrName not in self.chrBinBoundaries:
                exit("Chromosome name '{}' not found. Please check the correct spelling "
                     "of the chromosomes and try again".format(chrName))
            orig = self.chrBinBoundaries[chrName]
            new_order.extend(list(range(orig[0], orig[1])))
        self.reorderBins(new_order)

    def reorderBins(self, new_order):
        """
        reorders the rows and colums of the
        matrix according to the new order.
        The new order can be smaller
        than the original matrix. In that
        case, the ids not in the
        new order are removed.
        """
        orig_num_rows = self.matrix.shape[0]
        self.matrix = self.matrix[new_order, :][:, new_order]
        self.cut_intervals = [self.cut_intervals[x] for x in new_order]
        # reorder the masked bins
        # keep track of nan bins
        if len(self.nan_bins):
            _temp = np.zeros(orig_num_rows)
            _temp[self.nan_bins] = 1
            _temp = _temp[new_order]
            self.nan_bins = np.flatnonzero(_temp == 1)
        else:
            self.nan_bins = []

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def maskBins(self, bin_ids=None):
        """
        Mask the list of bins given. Mask means
        to remove the bins from the matrix,
        and keep the information about the intervals
        as masked
        """
        # print("self.cut_intervalsMASKBINS___START", self.cut_intervals)

        if bin_ids is None or len(bin_ids) == 0:
            return
        self.printchrtoremove(bin_ids, restore_masked_bins=False)
        try:
            # check if a masked bin already exists
            if len(self.orig_bin_ids) > 0:
                M = self.matrix.shape[0]
                previous_bin_ids = self.orig_bin_ids[M:]
                # merge new and old masked bins
                bin_ids = np.unique(np.concatenate([previous_bin_ids, self.orig_bin_ids[bin_ids]]))
                np.sort(bin_ids)
                self.restoreMaskedBins()
        except Exception:
            pass

        # join with existing nan_bins
        if self.nan_bins is not None and len(self.nan_bins) > 0:
            log.info("found existing {} nan bins that will be "
                     "included for masking ".format(len(self.nan_bins)))
            bin_ids = np.unique(np.concatenate([self.nan_bins, bin_ids]))
            self.nan_bins = []
        rows = cols = np.delete(list(range(self.matrix.shape[1])), bin_ids)

        self.matrix = self.matrix[rows, :][:, cols]

        # to keep track of removed bins
        # I add their ids to the end of the rows vector
        # to reverse the changes, I just need to do an argsort
        # to put the removed bins in place
        # log.debug("bins_ids {}".format(bin_ids))
        self.orig_bin_ids = np.concatenate([rows, bin_ids])

        new_cut_intervals = [self.cut_intervals[x] for x in rows]

        self.orig_cut_intervals = new_cut_intervals + [self.cut_intervals[x] for x in bin_ids]

        self.cut_intervals = new_cut_intervals

        self.interval_trees, self.chrBinBoundaries = self.intervalListToIntervalTree(self.cut_intervals)

        if self.correction_factors is not None:
            self.correction_factors = self.correction_factors[rows]

    def update_matrix(self, new_matrix, new_cut_intervals):
        """
        give a new matrix and list of cut intervals, the matrix,  cut intervals and
        the respective tree are updated
        :param new_matrix: now values for the sparse matrix
        :param new_cut_intervals: list of cut intervals, each entry being a tuple of the form
        (chrom, start, end, coverage)
        :return:
        """
        if len(self.orig_bin_ids) > 0:
            raise Exception("matrix contains masked bins. Restore masked bins first")

        assert len(new_cut_intervals) == new_matrix.shape[0], "matrix shape and len of cut intervals do not match"

        self.matrix = new_matrix
        self.cut_intervals = new_cut_intervals

        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

        self.nan_bins = np.flatnonzero(self.matrix.sum(0).A == 0)

    def restoreMaskedBins(self):
        """
        Puts backs into the matrix the bins
        removed


        Examples
        --------
        >>> from scipy.sparse import coo_matrix
        >>> row, col = np.triu_indices(5)
        >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
        ... ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
        >>> hic = hiCMatrix()
        >>> hic.nan_bins = []
        >>> matrix = np.array([
        ... [ 0, 10,  5, 3, 0],
        ... [ 0,  0, 15, 5, 1],
        ... [ 0,  0,  0, 7, 3],
        ... [ 0,  0,  0, 0, 1],
        ... [ 0,  0,  0, 0, 0]], dtype=np.int32)

        make the matrix symmetric:
        >>> hic.matrix = csr_matrix(matrix + matrix.T)
        >>> hic.setMatrix(csr_matrix(matrix + matrix.T), cut_intervals)

        Add masked bins masked bins
        >>> hic.maskBins([3])
        >>> hic.matrix.todense()
        matrix([[ 0, 10,  5,  0],
                [10,  0, 15,  1],
                [ 5, 15,  0,  3],
                [ 0,  1,  3,  0]], dtype=int32)
        >>> hic.cut_intervals
        [('a', 0, 10, 1), ('a', 10, 20, 1), ('a', 20, 30, 1), ('b', 40, 50, 1)]

        >>> hic.restoreMaskedBins()
        >>> hic.matrix.todense()
        matrix([[ 0., 10.,  5.,  0.,  0.],
                [10.,  0., 15.,  0.,  1.],
                [ 5., 15.,  0.,  0.,  3.],
                [ 0.,  0.,  0.,  0.,  0.],
                [ 0.,  1.,  3.,  0.,  0.]])

        >>> hic.cut_intervals
        [('a', 0, 10, 1), ('a', 10, 20, 1), ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
        """
        if len(self.orig_bin_ids) == 0:
            return
        # the rows to add are
        # as an empty sparse matrix
        M = self.matrix.shape[0]
        N = len(self.orig_bin_ids) - M
        rows_mat = csr_matrix((N, M))
        # cols to add
        cols_mat = csr_matrix((M + N, N))

        # add the rows and cols at the end of the
        # current matrix
        self.matrix = sparse_vstack([self.matrix, rows_mat])
        self.matrix = sparse_hstack([self.matrix, cols_mat], format='csr')

        # the new matrix has the right number of cols and rows, now
        # they need to be reordered to be back in their original places
        rows = cols = np.argsort(self.orig_bin_ids)
        self.matrix = self.matrix[rows, :][:, cols]
        self.cut_intervals = [self.orig_cut_intervals[x] for x in rows]
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)
        # set as nan_bins the masked bins that were restored
        self.nan_bins = self.orig_bin_ids[M:]

        if self.correction_factors is not None:
            # add missing values as nans at end of array
            self.correction_factors = np.concatenate([self.correction_factors,
                                                      np.repeat(np.nan, N)])
            # reorder array
            self.correction_factors = self.correction_factors[rows]

        # reset orig bins ids and cut intervals
        self.orig_bin_ids = []
        self.orig_cut_intervals = []
        log.info("masked bins were restored\n")

    def reorderMatrix(self, orig, dest):
        """
        Given a matrix, a region over the diagonal is moved from
        its origin to a new destination. With this method a
        new order of the chromosomes can be produced.
        :param orig: a tuple containing the indices of the region to be moved
        :param dest: the index of the region into which to insert
                     the section moved
        """

        rows = np.delete(list(range(self.matrix.shape[1])), range(orig[0], orig[1]))

        if dest > orig[1]:
            dest = dest - (orig[1] - orig[0])

        rows = cols = np.insert(
            rows, np.repeat(dest, orig[1] - orig[0]), list(range(orig[0], orig[1])))
        self.matrix = self.matrix[rows, :][:, cols]
        self.cut_intervals = [self.cut_intervals[x] for x in rows]
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

        if self.correction_factors is not None:
            self.correction_factors = self.correction_factors[rows]
        return

    def truncTrans(self, high=0.05):
        """Truncates trans contacts to remove blowouts
        Clip high counts in trans regions (i.e. between
        chromosomes) to the max value found in the 1-high*100
        percentile

        :param:  high : float, 0<high<1, optional
            Fraction of top trans interactions to be removed
        """
        mat = self.matrix.tocoo()
        dist_list = hiCMatrix.getDistList(mat.row, mat.col, self.cut_intervals)
        if np.count_nonzero(dist_list == -1) == 0:
            return
        max_inter = np.percentile(mat.data[dist_list == -1], (100 - high))
        mat.data[(mat.data >= max_inter) & (dist_list == -1)] == max_inter

        self.setMatrixValues(mat)

    def printchrtoremove(self, to_remove, label="Number of poor regions to remove", restore_masked_bins=True):
        """
        prints out the number of bin per chromosomes
        that will be removed
        """
        cnt = {}
        try:
            self.prev_to_remove
        except Exception:
            log.debug("No self.prev_to_remove defined, defining it now.")
            self.prev_to_remove = np.array([])

        # if the same information was already printed don't
        # show it again.
        if np.array_equal(self.prev_to_remove, to_remove):
            return

        if restore_masked_bins:
            try:
                # check if a masked bin already exists
                if len(self.orig_bin_ids) > 0:
                    log.info("Masked bins already present")
                    self.restoreMaskedBins()
            except Exception:
                pass
        for idx in to_remove:
            chrom = self.cut_intervals[idx][0]
            if chrom not in cnt:
                cnt[chrom] = 0
            cnt[chrom] += 1

        log.info('{}: {} {}'.format(label, len(to_remove), cnt))
        self.prev_to_remove = to_remove

    def get_chromosome_sizes(self):
        chrom_sizes = OrderedDict()
        for chrom, (start_bin, end_bin) in iteritems(self.chrBinBoundaries):
            chrom, start, end, _ = self.cut_intervals[end_bin - 1]
            chrom_sizes[chrom] = end

        return chrom_sizes

    def intervalListToIntervalTree(self, interval_list):
        """
        given an ordered list of (chromosome name, start, end)
        this is transformed to a number of interval trees,
        one for each chromosome
        """

        assert len(interval_list) > 0, "Interval list is empty"
        cut_int_tree = {}
        chrbin_boundaries = OrderedDict()
        intval_id = 0
        chr_start_id = 0
        previous_chrom = None
        for intval in interval_list:
            chrom, start, end = intval[0:3]
            start = int(start)
            end = int(end)
            if previous_chrom != chrom:
                if previous_chrom is None:
                    previous_chrom = chrom

                chrbin_boundaries[previous_chrom] = \
                    (chr_start_id, intval_id)
                chr_start_id = intval_id
                cut_int_tree[chrom] = IntervalTree()
                previous_chrom = chrom

            cut_int_tree[chrom].add(Interval(start, end, intval_id))

            intval_id += 1
        chrbin_boundaries[chrom] = (chr_start_id, intval_id)

        return cut_int_tree, chrbin_boundaries


def check_cooler(pFileName):
    if pFileName.endswith('.cool') or cooler.io.is_cooler(pFileName) or'.mcool' in pFileName:
        return True
    return False
