import logging
log = logging.getLogger(__name__)


class MatrixFile():

    def __init__(self, pMatrixFileName=None, pBedFile=None):
        self.matrixFileName = pMatrixFileName
        self.matrix = None
        self.cut_intervals = None
        self.nan_bins = None
        self.correction_factors = None
        self.distance_counts = None
        self.bedFile = pBedFile

    def load(self):
        log.error('Not implemented')

    def save(self):
        log.error('Not implemented')

    def is_of_type(self):
        log.error('Not implemented')

    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrix = pMatrix
        self.cut_intervals = pCutIntervals
        self.nan_bins = pNanBins
        self.correction_factors = pCorrectionFactors
        self.distance_counts = pDistanceCounts
