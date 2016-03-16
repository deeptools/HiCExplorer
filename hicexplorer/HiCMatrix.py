import numpy as np
import sys
from collections import OrderedDict
from scipy.sparse import csr_matrix, dia_matrix, coo_matrix
from scipy.sparse import vstack as sparse_vstack
from scipy.sparse import hstack as sparse_hstack
from scipy.sparse import triu, tril

## try to import pandas if exists
try:
    import pandas as pd
    pandas = True
except ImportError:
    pandas = False
    

from bx.intervals.intersection import IntervalTree, Interval
import gzip


class hiCMatrix:
    """
    Class to handle HiC matrices
    contains routines to get intrachromosomal distances
    get sub matrices by chrname.
    """

    def __init__(self, matrixFile=None, format=None, skiprows=None, chrnameList=None):
        self.correction_factors = None # this value is set in case a matrix was iteratively corrected
        self.non_homogeneous_warning_already_printed = False
        self.distanceCounts = None # only defined when getCountsByDistance is called

        if matrixFile:
            self.nan_bins = np.array([])
            if not format:
                if matrixFile[-4:] == ".npz":
                    format = 'npz'
                elif matrixFile[-3:] == '.gz':
                    format = 'dekker'
                # by default assume that the matrix format is .npz
                else:
                    format = 'npz'

            if format == 'npz':
                _ma = np.load(matrixFile)
                self.matrix = hiCMatrix.fillLowerTriangle(
                    _ma['matrix'].tolist())
                if 'dist_counts' not in _ma:
                    self.distance_counts = None
                else:
                    self.distance_counts = _ma['dist_counts'].tolist()

                self.cut_intervals = zip(_ma['chrNameList'], _ma['startList'],
                                         _ma['endList'], _ma['extraList'])

                assert len(self.cut_intervals) == self.matrix.shape[0], \
                       "Corrupted matrix file. Matrix size and " \
                       "matrix bin definitions do not correspond"
                if 'nan_bins' in _ma.keys():
                    self.nan_bins = _ma['nan_bins']
                    self.restoreMaskedBins()
                if 'correction_factors' in _ma.keys():
                    try:
                        # None value
                        # for correction_factors is saved by numpy
                        # as: array(None, dtype=object)
                        # Thus, to get the original None value
                        # the first item of the array is taken.
                        _ma['correction_factors'][0]
                        self.setCorrectionFactors(_ma['correction_factors'])
                    except IndexError:
                        pass

            elif format == 'dekker':
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
            elif format == 'lieberman': # lieberman format needs additional arguments : chrnameList
                lieberman_data = self.getLiebermanBins(filenameList = matrixFile, chrnameList = chrnameList)
                self.cut_intervals = lieberman_data['cut_intervals']
                self.matrix = lieberman_data['matrix']
            else:
                exit("matrix format not known.")

            self.interval_trees, self.chrBinBoundaries = \
                self.intervalListToIntervalTree(self.cut_intervals)


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
            # symetric matrix (below the main diagonal)
            # is cero. In this case, replace the lower
            # triangle using the upper triangle
            matrix = matrix + triu(matrix, 1).T

        return matrix

    def setMatrix(self, matrix, cut_intervals):
        """
        Initialize a matrix with a given matrix
        and cut_intervals. Mostly useful for
        testing.
        """

        self.setMatrixValues = matrix
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
        chrom, start, end, extra = zip(*self.cut_intervals)
        median = int(np.median(np.diff(start)))
        diff = np.array(end) - np.array(start)
        # check if the bin size is homogeneous
        if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
            if self.non_homogeneous_warning_already_printed is False:
                sys.stderr.write('WARNING: bin size is not homogeneous. Median {}\n'.format(median))
                self.non_homogeneous_warning_already_printed = True
#            raise Exception('bin size is not homogeneous')
        return median

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

    def getLiebermanBins(self, filenameList, chrnameList, pandas = pandas):
        """
        Reads a list of txt file in liberman's format and returns
        cut intervals and matrix. Each file is seperated by chr name
        and contains: locus1,locus2,and contact score seperated by tab.
        """

        ## Create empty row, col, value for the matrix

        row = np.array([]).astype("int")
        col = np.array([]).astype("int")
        value = np.array([])
        cut_intervals = []
        dim = 0
        ## for each chr, append the row, col, value to the first one. Extend the dim
        for i in range(0, len(filenameList)):
            if pandas == True:
                chrd = pd.read_csv(filenameList[i], sep = "\t", header=None)
                chrdata = chrd.as_matrix()
            else:
                print "Pandas unavailable. Reading files using numpy (slower).."
                chrdata = np.loadtxt(filenameList[i])

            # define resolution as the median of the difference of the rows
            # in the data table.

            resolution = np.median(np.diff(np.unique(np.sort(chrdata[:,1]))))

            chrcol = (chrdata[:, 1] / resolution).astype(int)
            chrrow = (chrdata[:, 0] / resolution).astype(int)

            chrdim = max(max(chrcol), max(chrrow)) + 1
            row = np.concatenate([row, chrrow + dim])
            col = np.concatenate([col, chrcol + dim])
            value = np.concatenate([value, chrdata[:, 2]])
            dim = dim + chrdim

            for _bin in range(chrdim):
                cut_intervals.append((chrnameList[i], _bin * resolution, (_bin + 1) * resolution,0))

        final_mat = coo_matrix((value, (row, col)), shape=(dim,dim))
        lieberman_data = dict(cut_intervals = cut_intervals, matrix = final_mat)
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

    def getBinPos(self, binIndex, includeBinId=False):
        """
        given a bin, it returns the chromosome name,
        start position and end position
        """
        if includeBinId:
            chrom, start, end, extra = self.cut_intervals[binIndex]
            # I dont understand this part. There was
            # a typo in original implementation?
            ret = (chrom, start, end, chrom)
#            ret = (self.nameList[binIndex], self.startList[binIndex], self.endList[binIndex], self.nameList[binIndex] )
        else:
            ret = self.cut_intervals[binIndex]
        return ret

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
            startbin = self.interval_trees[chrname].find(
                startpos, startpos + 1)[0].value
            endbin = self.interval_trees[chrname].find(
                endpos, endpos + 1)[0].value
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

    def getCountsByDistance(self, mean=False, per_chr=False):
        """
        computes counts for each intrachromosomal distance.
        better used with a corrected matrix

        returns a dictionary having as key the distance
        and as value an array containing the matrix values
        corresponding to that distance

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
        >>> del hic.distanceCounts
        >>> hic.getCountsByDistance(per_chr=True)
        {'a': {0: array([0, 0, 0, 0]), 10: array([10, 15,  7]), \
20: array([5, 5]), 30: array([3])}, 'b': {0: array([0])}}

        Test the removal of masked bins
        >>> hic.nan_bins = [3]
        >>> del hic.distanceCounts
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
            return self.distanceCounts

        # check that the matrix has bins of same size
        # otherwise the computations will fail
        chrom, start, end, extra = zip(*self.cut_intervals)
        median = int(np.median(np.diff(start)))
        diff = np.array(end) - np.array(start)
        # check if the bin size is homogeneous
        if len(np.flatnonzero(diff != median)) > (len(diff) * 0.01):
            # set the start position of a bin to the closest multiple
            # of the median
            def snap_nearest_multiple(start_x, m):
                resi = [-1*(start_x % m), -start_x % m]
                return start_x + resi[np.argmin(np.abs(resi))]
            start = [snap_nearest_multiple(x, median) for x in start]
            end = [snap_nearest_multiple(x, median) for x in end]
            cut_intervals = zip(chrom, start, end, extra)
            sys.stderr.write('[getCountsByDiscance] Bin size is not '
                             'homogeneous, setting \n'
                             'the bin distance to the median: {}\n'.format(median))

        else:
            cut_intervals = self.cut_intervals

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
                                               float(num_nan)/len(data)))
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
        # The row_res vector is boleean and contains 'False'
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
        self.distanceCounts = distance
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
        # is used to devide the data into
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
            del(self.distanceCounts)
        except AttributeError:
            pass
        self.matrix = mat
        return self.matrix


    def save_bing_ren(self, fileName):
        """
        Saves the matrix using bing ren's
        method which is chrom_name\tstart_bin\tend_bin\tvalues...
        """
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
            fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values) ) )

        fileh.close()

    def save_dekker(self, fileName):
        """
        Saves the matrix using dekker format
        """
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
            colNames.append("{}|dm3|{}:{}-{}".format(x, chrom, start, end)) # adds dm3 to the end (?problem..)

        fileh.write("#converted from npz\n")
        fileh.write("\t" + "\t".join(colNames) + "\n")
        for row in range(self.matrix.shape[0]):
            values = [str(x) for x in self.matrix[row, :].toarray().flatten()]
            fileh.write("{}\t{}\n".format(colNames[row], "\t".join(values) ) )

        fileh.close()

    def save_lieberman(self, fileName):
        """
        Saves the matrix using lieberman format. Given an output directory name and resolution of the matrix.
        """
        try:
            os.mkdir(fileName)
        except:
            print "Directory {} exists! Files will be overwritten.".format(fileName)

        lib_mat = self.matrix
        resolution = self.getBinSize()

        for chrom in self.interval_trees.keys():
                fileh = gzip.open("{}/chr{}_{}.gz".format(fileName,chr,resolution), 'w')
                rowNames = []
                chrstart, chrend = lib_mat.getChrBinRange(chrom)
                chrwise_mat = lib_mat.matrix[chrstart:chrend, chrstart:chrend]
                chrwise_mat_coo = triu(chrwise_mat, k=0, format='csr').tocoo()
                for x in range(chrwise_mat_coo.shape[0]):
                    start = chrwise_mat_coo.row[x,]*resolution
                    end = chrwise_mat_coo.col[x,]*resolution
                    data = chrwise_mat_coo.data[x,]
                    rowNames.append("{}\t{}\t{}".format(start, end,data))

                fileh.write("#converted from npz")
                fileh.write("\n" + "\n".join(rowNames) + "\n")
                fileh.close()

    def save(self, filename):
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
        assert len(correction_factors)==self.matrix.shape[0], \
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

        self.printchrtoremove(bin_ids)
        try:
            # check if a masked bin already exists
            if len(self.orig_bin_ids) > 0:
                print "Masked bins already present"
                M = self.matrix.shape[0]
                previous_bin_ids = self.orig_bin_ids[M:]
                # merge new and old masked bins
                bin_ids = np.unique(
                    np.concatenate(
                        [previous_bin_ids, self.orig_bin_ids[bin_ids]]))
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

    def printchrtoremove(self, to_remove,
                         label="Number of poor regions to remove"):
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

        sys.stderr.write('{}: {}\n{}\n'.format(label,
                len(to_remove),
                cnt))
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

            cut_int_tree[chrom].insert_interval(
                Interval(start, end, value=intval_id))

            intval_id += 1

        chrbin_boundaries[chrom] = (chr_start_id, intval_id)

        return (cut_int_tree, chrbin_boundaries)
