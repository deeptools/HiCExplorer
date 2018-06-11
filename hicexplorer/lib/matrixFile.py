import logging
log = logging.getLogger(__name__)

class MatrixFile():

    def __init__(self, pMatrixFileName):
        self.matrixFileName = pMatrixFileName
        self.matrix = None
        self.cut_intervals = None
        self.nan_bins = None
        self.correction_factors = None
        self.distance_counts = None

    def load(self):
        pass
    
    def save(self):
        pass
    
    def is_of_type(self):
        return False

    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrix = pMatrix
        self.cut_intervals = pCutIntervals
        self.nan_bins = pNanBins
        self.correction_factors = pCorrectionFactors
        self.distance_counts = pDistanceCounts