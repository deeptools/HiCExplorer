from .matrixFile import MatrixFile
import logging
log = logging.getLogger(__name__)
class Homer(MatrixFile):

    def __init__(self, pMatrixFile):
        super().__init__(pMatrixFile)
    
    def load(self):
        instances = []
        features = []
        data = []
        cut_intervals = []
        x = 0
        y = 0
        with open(pMatrixFile, 'r') as matrix_file:
            values = matrix_file.readline()
            for i, value in enumerate(values[2:]):
                chrom, start = value.split('-')

                cut_intervals.append((chrom, int(start), int(value[i + 1]), 1))

            for line in matrix_file:
                values = line.split('\t')
                for i, value in enumerate(values[2:]):
                    instances.append(x)
                    features.append(y)
                    data.append(int(value))
                    y += 1
                x += 1
                y = 0