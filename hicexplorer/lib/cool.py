import cooler
import logging
import numpy as np
from scipy.sparse import triu
import pandas as pd
from past.builtins import zip
from builtins import super
log = logging.getLogger(__name__)
from .matrixFile import MatrixFile
from hicexplorer.utilities import toString
from hicexplorer.utilities import convertNansToOnes


class Cool(MatrixFile, object):

    def __init__(self, pMatrixFile=None, pCooler_only_init=None):
        super().__init__(pMatrixFile)
        self.chrnameList = None
        self.correctionFactorTable = 'weight'
        self.correctionOperator = '*'

    def only_init(self):
        self.cooler_file = cooler.Cooler(self.matrixFileName)

    def getInformationCoolerBinNames(self):
        return cooler.Cooler(self.matrixFileName).bins().columns.values

    def load(self, pApplyCorrection=None, pMatrixOnly=None):
        log.debug('Load in cool format')
        log.debug('self.chrnameList {}'.format(self.chrnameList))
        if self.matrixFileName is None:
            log.info('No matrix is initalized')
        if pApplyCorrection is None:
            pApplyCorrection = True
        try:
            cooler_file = cooler.Cooler(self.matrixFileName)
        except Exception:
            log.info("Could not open cooler file. Maybe the path is wrong or the given node is not available.")
            log.info('The following file was tried to open: {}'.format(self.matrixFileName))
            log.info("The following nodes are available: {}".format(cooler.io.ls(self.matrixFileName.split("::")[0])))
            exit()

        if self.chrnameList is None:
            matrix = cooler_file.matrix(balance=False, sparse=True)[:].tocsr()
        else:
            if len(self.chrnameList) == 1:
                try:
                    matrix = cooler_file.matrix(balance=False, sparse=True).fetch(self.chrnameList[0]).tocsr()
                except ValueError:
                    exit("Wrong chromosome format. Please check UCSC / ensembl notation.")
            else:
                exit("Operation to load more as one region is not supported.")

        cut_intervals_data_frame = None
        correction_factors_data_frame = None

        if self.chrnameList is not None:
            if len(self.chrnameList) == 1:
                cut_intervals_data_frame = cooler_file.bins().fetch(self.chrnameList[0])

                if self.correctionFactorTable in cut_intervals_data_frame:
                    correction_factors_data_frame = cut_intervals_data_frame[self.correctionFactorTable]
            else:
                exit("Operation to load more than one chr from bins is not supported.")
        else:
            if pApplyCorrection and self.correctionFactorTable in cooler_file.bins():
                correction_factors_data_frame = cooler_file.bins()[[self.correctionFactorTable]][:]

            cut_intervals_data_frame = cooler_file.bins()[['chrom', 'start', 'end']][:]

        correction_factors = None
        # log.debug("{} {}".format(correction_factors_data_frame, pApplyCorrection))

        if correction_factors_data_frame is not None and pApplyCorrection:
            log.debug("Apply correction factors")
            # apply correction factors to matrix
            # a_i,j = a_i,j * c_i *c_j
            matrix.eliminate_zeros()
            matrix.data = matrix.data.astype(float)

            correction_factors = convertNansToOnes(np.array(correction_factors_data_frame.values).flatten())
            # apply only if there are not only 1's
            if np.sum(correction_factors) != len(correction_factors):
                instances, features = matrix.nonzero()
                instances_factors = correction_factors[instances]
                features_factors = correction_factors[features]
                instances_factors *= features_factors

                if self.correctionOperator == '*':
                    matrix.data *= instances_factors
                elif self.correctionOperator == '/':
                    matrix.data /= instances_factors

        cut_intervals = []

        for values in cut_intervals_data_frame.values:
            cut_intervals.append(tuple([toString(values[0]), values[1], values[2], 1.0]))

        # try to restore nan_bins.
        try:
            shape = matrix.shape[0] if matrix.shape[0] < matrix.shape[1] else matrix.shape[1]
            nan_bins = np.array(range(shape))
            nan_bins = np.setxor1d(nan_bins, matrix.indices)

            i = 0
            while i < len(nan_bins):
                if nan_bins[i] >= shape:
                    break
                i += 1
            nan_bins = nan_bins[:i]

        except Exception:
            nan_bins = None

        distance_counts = None

        # matrix = hiCMatrix.fillLowerTriangle(matrix)

        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors

    # def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
    #     super().set_matrix_variables(pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts)

    # def __create_empty_cool_file(self, pFileName):
    #     bins_data_frame = pd.DataFrame(columns=['chrom', 'start', 'end', 'weight'])
    #     matrix_data_frame = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'count'])
    #     cooler.io.create(cool_uri=pFileName,
    #                      bins=bins_data_frame,
    #                      pixels=matrix_data_frame)

    def save(self, pFileName, pSymmetric=True, pApplyCorrection=True):
        log.debug('Save in cool format')

        # log.info('self.matrix {}'.format(self.matrix))
        # log.info('self.nan_bins {}'.format(self.nan_bins))
        # log.info('self.cut_intervals {}'.format(self.cut_intervals))
        # log.info('self.correction_factors {}'.format(self.correction_factors))
        # log.info('pApplyCorrection {}'.format(pApplyCorrection))

        # for value in self.nan_bins:
        #

        self.matrix = self.matrix.tolil()
        if self.nan_bins is not None:
            self.matrix[self.nan_bins, :] = 0
            self.matrix[:, self.nan_bins] = 0
        self.matrix = self.matrix.tocsr()
        # log.info('self.matrix after nan handling{}'.format(self.matrix))

        for i in range(len(self.matrix.data)):
            if np.isnan(self.matrix.data[i]):
                self.matrix.data[i] = 0
        # log.info('self.matrix after nan handling II {}'.format(self.matrix))

        self.matrix.eliminate_zeros()
        # log.info('self.matrix after eliminate zeros{}'.format(self.matrix))

        # save only the upper triangle of the
        if pSymmetric:
            # symmetric matrix
            matrix = triu(self.matrix, format='csr')
            # log.debug('Symmetric {}'.format(pSymmetric))
        else:
            matrix = self.matrix
            # log.debug('SymmetricELSEs {}'.format(pSymmetric))

        # log.info('matrix after symmetric{}'.format(matrix))

        cut_intervals_ = []
        for value in self.cut_intervals:
            cut_intervals_.append(tuple((value[0], value[1], value[2])))

        bins_data_frame = pd.DataFrame(cut_intervals_, columns=['chrom', 'start', 'end'])

        # append correction factors if they exist
        if self.correction_factors is not None and pApplyCorrection:
            weight = convertNansToOnes(np.array(self.correction_factors).flatten())
            bins_data_frame = bins_data_frame.assign(weight=weight)

        # get only the upper triangle of the matrix to save to disk
        # upper_triangle = triu(self.matrix, k=0, format='csr')
        # create a tuple list and use it to create a data frame

        # save correction factors and original matrix

        # revert correction to store orginal matrix
        if self.correction_factors is not None and pApplyCorrection:
            log.info("Reverting correction factors on matrix...")
            instances, features = matrix.nonzero()
            self.correction_factors = np.array(self.correction_factors)

            # do not apply if correction factors are just 1's
            if np.sum(self.correction_factors) != len(self.correction_factors):
                instances_factors = self.correction_factors[instances]
                features_factors = self.correction_factors[features]

                instances_factors *= features_factors
                matrix.data /= instances_factors
                instances_factors = None
                features_factors = None

                matrix.data = np.rint(matrix.data)
                matrix.data = matrix.data.astype(int)

            data = matrix.data.tolist()

        else:

            instances, features = matrix.nonzero()
            data = matrix.data.tolist()

            if matrix.dtype not in [np.int32, int]:
                log.warning("Writing non-standard cooler matrix. Datatype of matrix['count'] is: {}".format(matrix.dtype))
                cooler._writer.COUNT_DTYPE = matrix.dtype

        if len(instances) == 0 and len(features) == 0:
            exit('No data present. Exit.')
        else:
            # log.debug(' data: {}'.format(data))
            matrix_tuple_list = zip(instances.tolist(), features.tolist(), data)
            # log.debug('Save cool, data: {}'.format(matrix_tuple_list))
            matrix_data_frame = pd.DataFrame(matrix_tuple_list, columns=['bin1_id', 'bin2_id', 'count'])

            cooler.io.create(cool_uri=pFileName,
                             bins=bins_data_frame,
                             pixels=matrix_data_frame,
                             append=False)
