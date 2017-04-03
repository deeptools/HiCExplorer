import numpy as np
import sys
import os
from collections import OrderedDict
from scipy.sparse import csr_matrix, dia_matrix, coo_matrix
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
from scipy.sparse import triu, tril
import tables
from intervaltree import IntervalTree, Interval

import gzip

# try to import pandas if exists
try:
    import pandas as pd
    pandas = True
except ImportError:
    pandas = False


class hiCMatrix:
    """
    Class to handle HiC matrices
    contains routines to get intrachromosomal distances
    get sub matrices by chrname.
    """

    def __init__(self, matrixFile=None, file_format=None, skiprows=None, chrnameList=None, bplimit=None):
        self.correction_factors = None  # this value is set in case a matrix was iteratively corrected
        self.non_homogeneous_warning_already_printed = False
        self.distance_counts = None  # only defined when getCountsByDistance is called
        self.bin_size = None
        self.bin_size_homogeneous = None  # track if the bins are equally spaced or not

        if matrixFile:
            self.nan_bins = np.array([])
            if not file_format:
                if matrixFile.endswith(".npz"):
                    file_format = 'npz'
                elif matrixFile.endswith(".h5"):
                    file_format = 'h5'
                elif matrixFile.endswith('.gz'):
                    file_format = 'dekker'
                # by default assume that the matrix file_format is hd5
                else:
                    file_format = 'h5'

            if file_format == 'h5' or file_format == 'hicexplorer':
                self.matrix, self.cut_intervals, self.nan_bins, self.distance_counts, self.correction_factors = \
                    hiCMatrix.load_h5(matrixFile)
                self.restoreMaskedBins()

            elif file_format == 'npz':
                self.matrix, self.cut_intervals, self.nan_bins, self.distance_counts, self.correction_factors = \
                    hiCMatrix.load_npz(matrixFile)
                self.restoreMaskedBins()

            elif file_format == 'dekker':
                self.cut_intervals = self.getDekkerBins(matrixFile)
                if not skiprows:
                    # +1 to also skip the column labels
                    skiprows = len(self.header) + 1
                self.matrix = csr_matrix(
                    np.loadtxt(matrixFile,
                               skiprows=skiprows,
                               usecols=range(1, len(self.cut_intervals) + 1)))
                """
                # convert nans to zeros
                self.matrix.data[np.isnan(self.matrix.data)] = 0
                self.matrix.eliminate_zeros()
                # mask zero count bins:
                row_sum = np.asarray(self.matrix.sum(axis=1)).flatten()
                self.maskBins(np.flatnonzero(row_sum==0))
                """
            elif file_format == 'lieberman':  # lieberman format needs additional arguments : chrnameList
                lieberman_data = self.getLiebermanBins(filenameList=matrixFile, chrnameList=chrnameList)
                self.cut_intervals = lieberman_data['cut_intervals']
                self.matrix = lieberman_data['matrix']

            else:
                exit("matrix format not known.")

            self.interval_trees, self.chrBinBoundaries = \
                self.intervalListToIntervalTree(self.cut_intervals)

    @staticmethod
    def load_h5(matrix_filename):
        """
        Loads a matrix stored in h5 format
        :param matrix_filename:
        :return: matrix, cut_intervals, nan_bins, distance_counts, correction_factors
        """
        with tables.open_file(matrix_filename) as f:
            parts = {}
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()

            matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                                shape=parts['shape'])
            matrix = hiCMatrix.fillLowerTriangle(matrix)
            # get intervals
            intvals = {}
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                intvals[interval_part] = getattr(f.root.intervals, interval_part).read()
            cut_intervals = zip(intvals['chr_list'], intvals['start_list'], intvals['end_list'], intvals['extra_list'])
            assert len(cut_intervals) == matrix.shape[0], \
                "Error loading matrix. Length of bin intervals ({}) is different than the " \
                "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

            # get nan_bins
            if hasattr(f.root, 'nan_bins'):
                nan_bins = f.root.nan_bins.read()
            else:
                nan_bins = np.array([])

            # get correction factors
            if hasattr(f.root, 'correction_factors'):
                correction_factors = f.root.correction_factors.read()
                assert len(correction_factors) == matrix.shape[0], \
                    "Error loading matrix. Length of correction factors does not" \
                    "match size of matrix"
            else:
                correction_factors = None

            # get correction factors
            if hasattr(f.root, 'distance_counts'):
                distance_counts = f.root.correction_factors.read()
            else:
                distance_counts = None

            return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    @staticmethod
    def load_npz(matrixFile):
        _ma = np.load(matrixFile)
        matrix = hiCMatrix.fillLowerTriangle(_ma['matrix'].tolist())
        if 'dist_counts' not in _ma:
            distance_counts = None
        else:
            distance_counts = _ma['dist_counts'].tolist()

        cut_intervals = zip(_ma['chrNameList'], _ma['startList'],
                            _ma['endList'], _ma['extraList'])

        assert len(cut_intervals) == matrix.shape[0], \
            "Corrupted matrix file. Matrix size and " \
            "matrix bin definitions do not correspond"
        if 'nan_bins' in _ma.keys():
            nan_bins = _ma['nan_bins']
        else:
            nan_bins = np.array([])

        if 'correction_factors' in _ma.keys():
            try:
                # None value
                # for correction_factors is saved by numpy
                # as: array(None, dtype=object)
                # Thus, to get the original None value
                # the first item of the array is taken.
                _ma['correction_factors'][0]
                assert len(_ma['correction_factors']) == matrix.shape[0], \
                    "length of correction factors and length of matrix are different."
                correction_factors = _ma['correction_factors']
            except IndexError:
                correction_factors = None
        else:
            correction_factors = None

        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    @staticmethod
    def fillLowerTriangle(matrix):
        """
        checks if the matrix is complete or if only half of the matrix was saved.
        Returns a whole matrix.

        >>> from scipy.sparse import csr_matrix
        >>> A = csr_matrix(np.array([[12,5,3,2,0],[0,11,4,1,1],
        ... [0,0,9,6,0], [0,0,0,10,0], [0,0,0,0,0]]))
        >>> B = hiCMatrix.fillLowerTriangle(A)
        >>> B.todense()
        matrix([[12,  5,  3,  2,  0],
                [ 5, 11,  4,  1,  1],
                [ 3,  4,  9,  6,  0],
                [ 2,  1,  6, 10,  0],
                [ 0,  1,  0,  0,  0]])

        """
        if tril(matrix, k=-1).sum() == 0:
            # this case means that the lower triangle of the
            # symmetric matrix (below the main diagonal)
            # is zero. In this case, replace the lower
            # triangle using the upper triangle
            matrix = matrix + triu(matrix, 1).T

        return matrix

    def setMatrix(self, matrix, cut_intervals):
        """
        Initialize a matrix with a given matrix
        and cut_intervals. Mostly useful for
        testing.
        """

        # check that the matrix is squared
        if matrix.shape[0] != matrix.shape[1]:
            raise Exception("Matrix is not squared. Shape is {}".format(matrix.shape))
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
            # check if the bin size is homogeneous
            if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
                self.bin_size_homogeneous = False
                if self.non_homogeneous_warning_already_printed is False:
                    sys.stderr.write('WARNING: bin size is not homogeneous. Median {}\n'.format(median))
                    self.non_homogeneous_warning_already_printed = True
            self.bin_size = median
        return self.bin_size

    def getDekkerBins(self, fileName):
        """
        Reads a gziped HiC matrix in Decker's format
        This format has row and column headers and comments
        with the prefix '#'. The following code skips the
        header while saving it in the self.header list
        and then reads and returns the column headers.
        The headers have the format bin_id|species|chrom:start-end
        for example 'bin1000000|dm3|chr4:1-50001'
        """

        i = 0
        self.header = []
        try:
            for line in gzip.open(fileName, 'r').readlines():
                if line.startswith("#"):
                    self.header.append(line)
                    continue
                i += 1
                if i == 1:  # this is the first line that is not
                            # a comment which contains the column headers
                            # that are equal to the row headers
                    colLabels = line.split("\t")
                    nameList = []
                    startList = []
                    endList = []
                    binIdList = []

                    for lab in colLabels[1:]:
                        (binId, species, position) = lab.split("|")
                        (chrName, pos) = position.split(":")
                        (start, end) = pos.split("-")
                        binIdList.append(binId)
                        nameList.append(chrName)
                        startList.append(int(start))
                        endList.append(int(end))
                    break

        except IOError:
            print "Error reading {}.\nDoes the file exists? "
            "Is it gzipped?".format(fileName)
            exit(1)

        return zip(nameList, startList, endList, binIdList)

    def getLiebermanBins(self, filenameList, chrnameList, pandas=pandas):
        """
        Reads a list of txt file in liberman's format and returns
        cut intervals and matrix. Each file is seperated by chr name
        and contains: locus1,locus2,and contact score seperated by tab.
        """

        # Create empty row, col, value for the matrix

        row = np.array([]).astype("int")
        col = np.array([]).astype("int")
        value = np.array([])
        cut_intervals = []
        dim = 0
        # for each chr, append the row, col, value to the first one. Extend the dim
        for i in range(0, len(filenameList)):
            if pandas is True:
                chrd = pd.read_csv(filenameList[i], sep="\t", header=None)
                chrdata = chrd.as_matrix()
            else:
                print "Pandas unavailable. Reading files using numpy (slower).."
                chrdata = np.loadtxt(filenameList[i])

            # define resolution as the median of the difference of the rows
            # in the data table.

            resolution = np.median(np.diff(np.unique(np.sort(chrdata[:, 1]))))

            chrcol = (chrdata[:, 1] / resolution).astype(int)
            chrrow = (chrdata[:, 0] / resolution).astype(int)

            chrdim = max(max(chrcol), max(chrrow)) + 1
            row = np.concatenate([row, chrrow + dim])
            col = np.concatenate([col, chrcol + dim])
            value = np.concatenate([value, chrdata[:, 2]])
            dim += chrdim

            for _bin in range(chrdim):
                cut_intervals.append((chrnameList[i], _bin * resolution, (_bin + 1) * resolution, 0))

        final_mat = coo_matrix((value, (row, col)), shape=(dim, dim))
        lieberman_data = dict(cut_intervals=cut_intervals, matrix=final_mat)
        return lieberman_data

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
        return self.chrBinBoundaries.keys()

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
            self.interval_trees[chrname]
        except KeyError:
            """
            print "chromosome: {} name not found in matrix".format(chrName)
            print "valid names are:"
            print interval_trees.keys()
            """
            return None
        try:
            startpos = int(startpos)
            endpos = int(endpos)
        except:
            print "{} or {}  are not valid " \
                "position values.".format(startpos, endpos)
            exit()

        try:
            startbin = sorted(self.interval_trees[chrname][startpos:startpos + 1])[0].data
            endbin = sorted(self.interval_trees[chrname][endpos:endpos + 1])[0].data
        except IndexError:
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
        >>> row, col = np.triu_indices(5)
        >>> cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
        ... ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
        >>> dist_list, chrom_list = hiCMatrix.getDistList(row, col,
        ... cut_intervals)
        >>> coo_matrix((dist_list, (row, col)), shape=(5,5)).todense()
        matrix([[ 0, 10, 20, 30, -1],
                [ 0,  0, 10, 20, -1],
                [ 0,  0,  0, 10, -1],
                [ 0,  0,  0,  0, -1],
                [ 0,  0,  0,  0,  0]])
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
            sys.stderr.write('[getCountsByDistance] Bin size is not '
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
        matrix([[ 1. ,  0.8,  1. ,  1. ,  0. ],
                [ 0. ,  4. ,  1.5,  1. ,  1. ],
                [ 0. ,  0. ,  0. ,  0.7,  2. ],
                [ 0. ,  0. ,  0. ,  0. ,  1. ],
                [ 0. ,  0. ,  0. ,  0. ,  0. ]])

        >>> hic.matrix = csr_matrix(matrix)
        >>> hic.convert_to_obs_exp_matrix(maxdepth=20).todense()
        matrix([[ 1. ,  0.8,  1. ,  0. ,  0. ],
                [ 0. ,  4. ,  1.5,  1. ,  0. ],
                [ 0. ,  0. ,  0. ,  0.7,  nan],
                [ 0. ,  0. ,  0. ,  0. ,  nan],
                [ 0. ,  0. ,  0. ,  0. ,  0. ]])

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
                estimated_size_dense_matrix = m_size**2 * 8
                if estimated_size_dense_matrix > 100e6:
                    sys.stderr.write("To compute z-scores a dense matrix is required. This will use \n"
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
            diag_mat_ones = diags(np.repeat([1], m_size * depth).reshape(depth, m_size), range(depth))

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
            chrom_sizes['all'] = np.array([v[1] - v[0] for k, v in self.chrBinBoundaries.iteritems()])
            chrom_range['all'] = (0, self.matrix.shape[0])

        for chrname, submatrix in chr_submatrix.iteritems():
            sys.stderr.write("processing chromosome {}\n".format(chrname))
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
                    diagonal_length = sum([size - (bin_dist_plus_one - 1) for size in chrom_sizes[chrname] if size > (bin_dist_plus_one - 1)])

                # the diagonal length should contain the number of values at a certain distance.
                # If the matrix is dense, the distance_len[bin_dist_plus_one] correctly contains the number of values
                # If the matrix is equally spaced, then, the diagonal_length as computed before is accurate.
                # But, if the matrix is both sparse and with unequal bins, then none of the above methods is
                # accurate but the the diagonal_length as computed before will be closer.
                diagonal_length = max(diagonal_length, distance_len[bin_dist_plus_one])

                if diagonal_length == 0:
                    mu[bin_dist_plus_one] = np.nan
                else:
                    mu[bin_dist_plus_one] = np.float64(sum_value) / diagonal_length

                if np.isnan(sum_value):
                    sys.stderr.write("nan value found for distance {}\n".format((bin_dist_plus_one - 1) * binsize))

                # if zscore is needed, compute standard deviation: std = sqrt(mean(abs(x - x.mean())**2))
                if zscore:
                    values_sqrt_diff = \
                        np.abs((submatrix.data[dist_list == bin_dist_plus_one] - mu[bin_dist_plus_one])**2)
                    # the standard deviation is the sum of the differences with mu squared (value variable)
                    # plus all zeros that are not included in the sparse matrix
                    # for which the standard deviation is
                    # (0 - mu)**2 = (mu)**2
                    # The number of zeros is the diagonal length - the length of the non zero values
                    zero_values_sqrt_diff_sum = (diagonal_length - len(values_sqrt_diff)) * mu[bin_dist_plus_one]**2

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

    def getCountsByDistance(self, mean=False, per_chr=False):
        """
        computes counts for each intrachromosomal distance.
        better used with a corrected matrix


        Parameters
        ----------
        mean : if set to true (default) the mean of the distance value is returned instead of a list with each of the
                elements
        per_chr: set to true if the computation should be done per chromosome

        Returns
        -------
        returns a dictionary having as key the distance
        and as value an array containing the matrix values
        corresponding to that distance


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
        ... [ 0,  0,  0, 0, 0]])

        make the matrix symmetric:
        >>> hic.matrix = csr_matrix(matrix + matrix.T)
        >>> hic.setMatrix(csr_matrix(matrix + matrix.T), cut_intervals)
        >>> hic.getCountsByDistance()
        {0: array([0, 0, 0, 0, 0]), 10: array([10, 15,  7]), \
20: array([5, 5]), 30: array([3]), -1: array([0, 1, 3, 1])}

        Test get distance counts per chromosome
        >>> hic.distance_counts = None
        >>> hic.getCountsByDistance(per_chr=True)
        {'a': {0: array([0, 0, 0, 0]), 10: array([10, 15,  7]), \
20: array([5, 5]), 30: array([3])}, 'b': {0: array([0])}}

        Test the removal of masked bins
        >>> hic.nan_bins = [3]
        >>> hic.distance_counts = None
        >>> hic.getCountsByDistance()
        {0: array([0, 0, 0, 0]), 10: array([10, 15]), 20: array([5]), \
-1: array([0, 1, 3])}

        Test bins that are of different size
        >>> cut_intervals = [('a', 0, 12, 1), ('a', 12, 25, 1),
        ... ('a', 25, 35, 1), ('a', 35, 41, 1), ('b', 41, 50, 1)]
        >>> hic = hiCMatrix()
        >>> hic.nan_bins = []
        >>> hic.matrix = csr_matrix(matrix + matrix.T)
        >>> hic.setMatrix(csr_matrix(matrix + matrix.T), cut_intervals)
        >>> hic.getCountsByDistance()
        {0: array([0, 0, 0, 0, 0]), 33: array([3]), 11: array([10, 15,  7]), \
22: array([5, 5]), -1: array([0, 1, 3, 1])}
        """

        if self.distance_counts:
            return self.distance_counts

        cut_intervals = hiCMatrix.fit_cut_intervals(self.cut_intervals)

        M, N = self.matrix.shape
        # get the row and col indices of the upper triangle
        # of a square matrix of size M
        rows, cols = np.triu_indices(M)

        # put a negative number to all masked bins
        #
        # the upper triangle of the sparse data matrix is
        # converted into a dense matrix and represented
        # as a vector
        data = np.asarray(self.matrix.todense()[(rows, cols)]).flatten()

        # convert nans to zeros. Otherwise the computations will fail
        if np.any(np.isnan(data)):
            num_nan = len(np.flatnonzero(np.isnan(data)))
            sys.stderr.write("converting {} ({:.2f}) nans "
                             "to zeros".format(num_nan,
                                               float(num_nan) / len(data)))
            data[np.isnan(data)] = 0

        # get a vector of all distances. The order is the same
        # as the 'data' vector
        dist_list, chrom_list = hiCMatrix.getDistList(rows, cols,
                                                      cut_intervals)

        """
        remove all zeros corresponding to masked bins, needs to be done
        after the dist_list computation to get the distances right
        """
        # compute the overlap between rows and self.nan_bins
        # The row_res vector is boolean and contains 'False'
        # for all row ids that are in self.nan_bin
        row_res = np.invert(np.in1d(rows, self.nan_bins))
        col_res = np.invert(np.in1d(cols, self.nan_bins))

        # get data values for which row_res and col_res are both true
        data = data[row_res & col_res]
        dist_list = dist_list[row_res & col_res]
        chrom_list = chrom_list[row_res & col_res]

        if per_chr:
            distance = {}
            for chr_name in self.getChrNames():
                _chr_data = data[chrom_list == chr_name]
                _chr_dist_list = dist_list[chrom_list == chr_name]
                distance[chr_name] = hiCMatrix.dist_list_to_dict(_chr_data,
                                                                 _chr_dist_list)
        else:
            distance = hiCMatrix.dist_list_to_dict(data, dist_list)
        self.distance_counts = distance
        if mean:
            return [np.mean(distance[x]) for x in range(len(distance.keys()))]
        else:
            return distance

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
        for index in xrange(len(distance_unique)):
            distance[distance_unique[index]] = groups[index]

        return distance

    @staticmethod
    def getUnwantedChrs():
        # too drosophila specific. Should remove this function
        unwanted_chr = {'chrM', 'chrYHet', 'chrXHet', 'chrUextra', 'chrU', 'chr3RHet', 'chr3LHet', 'chr2RHet',
                        'chr2LHet'}
        return unwanted_chr

    def filterUnwantedChr(self, chromosome=None):
        if chromosome:
            self.keepOnlyTheseChr([chromosome])
        else:
            wanted = set(self.getChrNames()).difference(self.getUnwantedChrs())
            self.keepOnlyTheseChr(wanted)

    def keepOnlyTheseChr(self, chromosome_list):
        """
        given a list of chromosome names,
        these are kept, while any other is removed
        from the matrix
        """
        if isinstance(chromosome_list, str):
            chromosome_list = [chromosome_list]

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

        for chrName in self.interval_trees.keys():
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

    def save_bing_ren(self, fileName):
        """
        Saves the matrix using bing ren's
        method which is chrom_name\tstart_bin\tend_bin\tvalues...
        """
        if os.path.isdir(fileName):
            exit('Please specify a file name to save the data. The given file name is a folder: \n{}\n'.format(fileName))

        if fileName[-3:] != '.gz':
            fileName += '.gz'

        try:
            fileh = gzip.open(fileName, 'w')
        except:
            msg = "{} file can't be opened for writing".format(fileName)
            raise msg

        colNames = []
        for x in range(self.matrix.shape[0]):
            chrom, start, end = self.cut_intervals[x][0:3]
            colNames.append("{}\t{}\t{}".format(chrom, start, end))

        for row in range(self.matrix.shape[0]):
            values = [str(x) for x in self.matrix[row, :].toarray().flatten()]
            fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values)))

        fileh.close()

    def save_dekker(self, fileName):
        """
        Saves the matrix using dekker format
        """
        if os.path.isdir(fileName):
            exit('Please specify a file name to save the data. The given file name is a folder: \n{}\n'.format(fileName))

        if fileName[-3:] != '.gz':
            fileName += '.gz'

        try:
            fileh = gzip.open(fileName, 'w')
        except:
            msg = "{} file can't be opened for writing".format(fileName)
            raise msg

        colNames = []
        for x in range(self.matrix.shape[0]):
            chrom, start, end = self.cut_intervals[x][0:3]
            colNames.append("{}|--|{}:{}-{}".format(x, chrom, start, end))  # adds dm3 to the end (?problem..)

        fileh.write("#converted from hicexplorer\n")
        fileh.write("\t" + "\t".join(colNames) + "\n")
        for row in range(self.matrix.shape[0]):
            values = [str(x) for x in self.matrix[row, :].toarray().flatten()]
            fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values)))

        fileh.close()

    def save_lieberman(self, fileName):
        """
        Saves the matrix using lieberman format. Given an output directory name and resolution of the matrix.

        In contrast to other methods, when saving using liebermans format a folder is required
        where the data is saved per chromosome
        """

        if os.path.isfile(fileName):
            exit("a file with the same name as the given path exists ({}). Please choose a folder.".format(fileName))

        if not os.path.isdir(fileName):
            os.mkdir(fileName)

        resolution = self.getBinSize()
        if not self.bin_size_homogeneous:
            sys.stderr.write("*WARNING* The hic matrix contains bins of difference size but\n"
                             "the 'lieberman' format requires equally spaced bins. The program\n"
                             "will proceed but the results may be unreliable.\n")

        for chrom in self.interval_trees.keys():
            chrstart, chrend = self.getChrBinRange(chrom)
            chrwise_mat = self.matrix[chrstart:chrend, chrstart:chrend]
            if len(chrwise_mat.data) > 0:
                sys.stderr.write("Saving chromosome {}...\n".format(chrom))
                chrwise_mat_coo = triu(chrwise_mat, k=0, format='csr').tocoo()
                start = chrwise_mat_coo.row * resolution
                end = chrwise_mat_coo.col * resolution
                fileh = gzip.open("{}/chr{}.gz".format(fileName, chrom), 'w')
                fileh.write("#converted from HiCExplorer format\n")
                for idx in range(len(start)):
                    fileh.write("{}\t{}\t{}\n".format(start[idx], end[idx], chrwise_mat_coo.data[idx]))
                fileh.close()

    def save_GInteractions(self, fileName):
        self.restoreMaskedBins()
        print self.matrix.shape
        mat_coo = triu(self.matrix, k=0, format='csr').tocoo()
        fileh = open("{}.tsv".format(fileName), 'w')
        for idx, counts in enumerate(mat_coo.data):
            chr_row, start_row, end_row, _ = self.cut_intervals[mat_coo.row[idx]]
            chr_col, start_col, end_col, _ = self.cut_intervals[mat_coo.col[idx]]
            fileh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_row, int(start_row), int(end_row),
                                                              chr_col, int(start_col), int(end_col), counts))
        fileh.close()

    def save(self, filename):
        """
        Saves a matrix using hdf5 format
        :param filename:
        :return: None
        """
        self.restoreMaskedBins()
        if not filename.endswith(".h5"):
            filename += ".h5"

        # if the file name already exists
        # try to find a new suitable name
        if os.path.isfile(filename):
            sys.stderr.write("*WARNING* File already exists {}\n "
                             "Overwriting ...\n".format(filename))

            from os import unlink
            unlink(filename)
        try:
            nan_bins = np.array(self.nan_bins)
        except:
            nan_bins = np.array([])
        # save only the upper triangle of the
        # symmetric matrix
        matrix = triu(self.matrix, k=0, format='csr')
        filters = tables.Filters(complevel=5, complib='blosc')
        with tables.open_file(filename, mode="w", title="HiCExplorer matrix") as h5file:
            matrix_group = h5file.create_group("/", "matrix", )
            # save the parts of the csr matrix
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                arr = np.array(getattr(matrix, matrix_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(matrix_group, matrix_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save the matrix intervals
            intervals_group = h5file.create_group("/", "intervals", )
            chr_list, start_list, end_list, extra_list = zip(*self.cut_intervals)
            for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
                arr = np.array(eval(interval_part))
                atom = tables.Atom.from_dtype(arr.dtype)
                ds = h5file.create_carray(intervals_group, interval_part, atom,
                                          shape=arr.shape,
                                          filters=filters)
                ds[:] = arr

            # save nan bins
            if len(nan_bins):
                atom = tables.Atom.from_dtype(nan_bins.dtype)
                ds = h5file.create_carray(h5file.root, 'nan_bins', atom,
                                          shape=nan_bins.shape,
                                          filters=filters)
                ds[:] = nan_bins

            # save corrections factors
            if self.correction_factors is not None and len(self.correction_factors):
                self.correction_factors = np.array(self.correction_factors)
                atom = tables.Atom.from_dtype(self.correction_factors.dtype)
                ds = h5file.create_carray(h5file.root, 'correction_factors', atom,
                                          shape=self.correction_factors.shape,
                                          filters=filters)
                ds[:] = np.array(self.correction_factors)

            # save distance counts
            if self.distance_counts is not None and len(self.distance_counts):
                atom = tables.Atom.from_dtype(self.distance_counts.dtype)
                ds = h5file.create_carray(h5file.root, 'distance_counts', atom,
                                          shape=self.distance_counts.shape,
                                          filters=filters)
                ds[:] = np.array(self.distance_counts)

    def save_npz(self, filename):
        """
        saves using the numpy npz format. Adequate for small samples.
        """
        self.restoreMaskedBins()
        chrNameList, startList, endList, extraList = zip(*self.cut_intervals)
        try:
            nan_bins = self.nan_bins
        except:
            nan_bins = np.array([])
        # save only the upper triangle of the
        # symmetric matrix
        matrix = triu(self.matrix, k=0, format='csr')
        try:
            np.savez(
                filename, matrix=matrix, chrNameList=chrNameList,
                startList=startList, endList=endList, extraList=extraList,
                nan_bins=nan_bins, correction_factors=self.correction_factors)
        except Exception as e:
            print "error saving matrix: {}".format(e)
            try:
                print "Matrix can not be saved because is too big!"
                print "Eliminating entries with only one count."

                # try to remove noise by deleting 1
                matrix.data = matrix.data - 1
                matrix.eliminate_zeros()
                np.savez(
                    filename, matrix=matrix, chrNameList=chrNameList,
                    startList=startList, endList=endList, extraList=extraList,
                    nan_bins=nan_bins)
            except:
                print "Matrix can not be saved because is too big!"
            exit()

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
        >>> hic.setMatrix(csr_matrix(matrix + matrix.T), cut_intervals)
        >>> hic.filterOutInterChrCounts().todense()
        matrix([[ 0, 10,  5,  0,  0],
                [10,  0, 15,  0,  0],
                [ 5, 15,  0,  0,  0],
                [ 0,  0,  0,  0,  1],
                [ 0,  0,  0,  1,  0]])
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

    def reorderChromosomes_old(self, new_chr_order):
        if len(new_chr_order) != len(self.chrBinBoundaries):
            return
        dest = 0
        for chrName in new_chr_order:
            orig = self.chrBinBoundaries[chrName]

            self.reorderMatrix(orig, dest)
            dest += orig[1] - orig[0]

    def reorderChromosomes(self, new_chr_order):
        new_order = []
        for chrName in new_chr_order:
            # check that the chromosome names are valid
            if chrName not in self.chrBinBoundaries:
                exit("Chromosome name '{}' not found. Please check the correct spelling "
                     "of the chromosomes and try again".format(chrName))
            orig = self.chrBinBoundaries[chrName]
            new_order.extend(range(orig[0], orig[1]))
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

    def removeBins(self, bin_ids):
        """ given an array of ids, all rows and columns
        matching those ids are removed
        """
        rows = cols = np.delete(range(self.matrix.shape[1]), bin_ids)

        self.matrix = self.matrix[rows, :][:, cols]
        self.cut_intervals = [self.cut_intervals[x] for x in rows]
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

    def maskBins(self, bin_ids):
        """
        Mask the list of bins given. Mask means
        to remove the bins from the matrix,
        and keep the information about the intervals
        as masked
        """
        if len(bin_ids) == 0:
            return

        self.printchrtoremove(bin_ids, restore_masked_bins=False)
        try:
            # check if a masked bin already exists
            if len(self.orig_bin_ids) > 0:
                print "Masked bins already present"
                M = self.matrix.shape[0]
                previous_bin_ids = self.orig_bin_ids[M:]
                # merge new and old masked bins
                bin_ids = np.unique(np.concatenate([previous_bin_ids, self.orig_bin_ids[bin_ids]]))
                np.sort(bin_ids)
                self.restoreMaskedBins()
        except:
            pass

        # join with existing nan_bins
        if len(self.nan_bins) > 0:
            print "found existing {} nan bins that will be " \
                "included for masking ".format(len(self.nan_bins))
            bin_ids = np.unique(np.concatenate([self.nan_bins, bin_ids]))
            self.nan_bins = []
        rows = cols = np.delete(range(self.matrix.shape[1]), bin_ids)
        self.matrix = self.matrix[rows, :][:, cols]

        # to keep track of removed bins
        # I add their ids to the end of the rows vector
        # to reverse the changes, I just need to do an argsort
        # to put the removed bins in place
        self.orig_bin_ids = np.concatenate([rows, bin_ids])

        new_cut_intervals = [self.cut_intervals[x] for x in rows]

        self.orig_cut_intervals = \
            new_cut_intervals + [self.cut_intervals[x] for x in bin_ids]

        self.cut_intervals = new_cut_intervals
#        self.nan_bins = np.intersect1d(self.nan_bins, rows)
        self.interval_trees, self.chrBinBoundaries = \
            self.intervalListToIntervalTree(self.cut_intervals)

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
        if hasattr(self, 'orig_bin_ids'):
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
        """
        try:
            self.orig_bin_ids
        except AttributeError:
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
        del(self.orig_cut_intervals, self.orig_bin_ids)
        sys.stderr.write("masked bins were restored\n")

    def reorderMatrix(self, orig, dest):
        """
        Given a matrix, a region over the diagonal is moved from
        its origin to a new destination. With this method a
        new order of the chromosomes can be produced.
        :param orig: a tuple containing the indices of the region to be moved
        :param dest: the index of the region into which to insert
                     the section moved
        """

        rows = np.delete(range(self.matrix.shape[1]), range(orig[0], orig[1]))

        if dest > orig[1]:
            dest = dest - (orig[1] - orig[0])

        rows = cols = np.insert(
            rows, np.repeat(dest, orig[1] - orig[0]), range(orig[0], orig[1]))
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

    def truncTrans_bk(self, high=0.0005):
        """Truncates trans contacts to remove blowouts

        :param:  high : float, 0<high<1, optional
            Fraction of top trans interactions to be removed
        """
        M, N = self.matrix.shape
        rows, cols = np.triu_indices(M, 1)
        dist_list = hiCMatrix.getDistList(rows, cols)
        data = np.asarray(self.matrix.todense()[(rows, cols)]).flatten()

        max_inter = np.percentile(data[dist_list == -1], 100. * (1 - high))

        to_clip = np.flatnonzero((data >= max_inter) & (dist_list == -1))
        rows_to_clip = rows[to_clip]
        cols_to_clip = cols[to_clip]
        mat = self.matrix.tolil()
        mat[(rows_to_clip, cols_to_clip)] = max_inter
        mat[(cols_to_clip, rows_to_clip)] = max_inter

        self.setMatrixValues(mat)

    def removePoorRegions(self, cutoff=2):
        """Removes "cutoff" percent of bins with least counts.
        Also, removes bins that have less than sequencedFraction
        unique reads  with respect
        to the total number of reads found in the bin. (which
        includes multi reads).

        :param:  cutoff : int, 0<cutoff<100
        Percent of lowest-counts bins to be removed

        """
        # remove all nan values
        self.matrix.data[np.isnan(self.matrix.data)] = 0
        # using the row sum, find and remove those
        # bins that have the botton 1%
        row_sum = np.asarray(self.matrix.sum(axis=1)).flatten()
        min_threshold = np.percentile(row_sum, cutoff)
        to_remove = np.flatnonzero(row_sum <= min_threshold)
        # code to count were the bins are going to be removed
        chrom, start, end, cov = zip(*self.cut_intervals)
        chrom = np.sort(np.array(chrom)[to_remove])
        chr_count = np.diff(np.concatenate(
            [np.unique(np.sort(chrom), return_index=True)[1],
             [len(chrom)]]))
        chr_dict = dict(zip(chrom[chr_count], chr_count))
        sys.stderr.write('num poor regions to remove {}\n{}\n'.format(
            len(to_remove),
            chr_dict))
        self.maskBins(to_remove)

    def printchrtoremove(self, to_remove, label="Number of poor regions to remove", restore_masked_bins=True):
        """
        prints out the number of bin per chromosomes
        that will be removed
        """
        cnt = {}
        try:
            self.prev_to_remove
        except:
            self.prev_to_remove = np.array([])

        # if the same information was already printed don't
        # show it again.
        if np.array_equal(self.prev_to_remove, to_remove):
            return

        if restore_masked_bins:
            try:
                # check if a masked bin already exists
                if len(self.orig_bin_ids) > 0:
                    print "Masked bins already present"
                    self.restoreMaskedBins()
            except:
                pass
        for idx in to_remove:
            chrom = self.cut_intervals[idx][0]
            if chrom not in cnt:
                cnt[chrom] = 0
            cnt[chrom] += 1

        sys.stderr.write('{}: {}\n{}\n'.format(label, len(to_remove), cnt))
        self.prev_to_remove = to_remove

    def removeBySequencedCount(self, sequencedFraction=0.5):
        """
        Removes bins that have less than sequencedFraction
        unique reads  with respect
        to the total number of reads found in the bin. (which
        includes multi reads).

        :param:   sequencedFraction: float, optional, 0<x<1
            Fraction of the bin that needs to be sequenced in order
            to keep the bin

        """
        chrNameList, startList, endList, coverage = zip(*self.cut_intervals)
        if type(coverage[0]) != np.float64:
            return
        to_remove = np.flatnonzero(np.array(coverage) < sequencedFraction)
        self.maskBins(to_remove)
        return to_remove

    def get_chromosome_sizes(self):
        chrom_sizes = OrderedDict()
        for chrom, (start_bin, end_bin) in self.chrBinBoundaries.iteritems():
            chrom, start, end, _ = self.cut_intervals[end_bin - 1]
            chrom_sizes[chrom] = end
        return chrom_sizes

    @staticmethod
    def intervalListToIntervalTree(interval_list):
        """
        given an ordered list of (chromosome name, start, end)
        this is transformed to a number of interval trees,
        one for each chromosome
        """

        assert len(interval_list) > 0, "List is empty"
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
