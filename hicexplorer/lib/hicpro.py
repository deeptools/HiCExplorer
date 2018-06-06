from .matrixFile import MatrixFile
import logging
log = logging.getLogger(__name__)
class Hicpro(MatrixFile):

    def __init__(self, pMatrixFile, pBedFile):
        super().__init__(pMatrixFile)
        self.__bed_file = pBedFile
    
    
    def load(self):
        instances = []
        featues = []
        data = []
        with open(pMatrixFile, 'r') as matrix_file:
            for line in matrix_file:
                x, y, value = line.split('\t')
                instances.append(int(x))
                features.append(int(y))
                data.append(float(data))
        cut_intervals = []
        with open(pBedFile, 'r') as bed_file:
            for line in bed_file:
                chrom, start, end, value = line.split('\t')
                cut_intervals.append((chrom, start, end, value))
        shape = len(cut_intervals)
        matrix = csr_matrix((data, (instances, features)), shape=(shape, shape))

        nan_bins = None
        distance_counts = None
        correction_factors = None
        return matrix, cut_intervals, nan_bins, distance_counts, correction_factors
