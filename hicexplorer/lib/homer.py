from .matrixFile import MatrixFile
from scipy.sparse import csr_matrix
import logging
log = logging.getLogger(__name__)


class Homer(MatrixFile):

    def __init__(self, pMatrixFile):
        super().__init__(pMatrixFile)

    def load(self):
        # instances = []
        # features = []
        # data = []
        cut_intervals = []
        # x = 0
        # y = 0
        with open(self.matrixFileName, 'r') as matrix_file:
            values = matrix_file.readline()
            values = values.split('\t')

            # get bin size
            start_first = int(values[2].strip().split('-')[1])
            start_second = int(values[3].strip().split('-')[1])
            bin_size = start_second - start_first
            for i, value in enumerate(values[2:]):
                chrom, start = value.strip().split('-')
                log.debug('chrom {} start {}'.format(chrom, start))

                cut_intervals.append((chrom, int(start), int(start) + bin_size, 1))

            matrix_dense = []
            for line in matrix_file:
                values = line.split('\t')
                data = []
                for i, value in enumerate(values[2:]):
                    data.append(float(value))
                matrix_dense.append(data)
            log.debug('data {}'.format(data))

        matrix = csr_matrix(matrix_dense)
        nan_bins = None
        distance_counts = None
        correction_factors = None
        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors
