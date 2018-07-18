from .matrixFile import MatrixFile
from scipy.sparse import csr_matrix

import logging
log = logging.getLogger(__name__)


class Hicpro(MatrixFile):

    def __init__(self, pMatrixFile, pBedFile):
        super().__init__(pMatrixFileName=pMatrixFile, pBedFile=pBedFile)
        # self.bedFile = pBedFile

    def load(self):
        instances = []
        features = []
        data = []
        with open(self.matrixFileName, 'r') as matrix_file:
            for line in matrix_file:
                x, y, value = line.strip().split('\t')
                instances.append(int(x) - 1)
                features.append(int(y) - 1)
                data.append(float(value))
        cut_intervals = []
        with open(self.bedFile, 'r') as bed_file:
            for line in bed_file:
                chrom, start, end, value = line.strip().split('\t')
                cut_intervals.append((chrom, int(start), int(end), int(value)))
        shape = len(cut_intervals)
        log.debug('instances {}'.format(len(instances)))
        log.debug('max instances {}'.format(max(instances)))
        log.debug('min instances {}'.format(min(instances)))

        log.debug('features {}'.format(len(features)))
        log.debug('max features {}'.format(max(features)))
        log.debug('min features {}'.format(min(features)))

        log.debug('data {}'.format(len(data)))
        log.debug('cut_intervals {}'.format(len(cut_intervals)))
        matrix = csr_matrix((data, (instances, features)), shape=(shape, shape))

        nan_bins = None
        distance_counts = None
        correction_factors = None
        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors
