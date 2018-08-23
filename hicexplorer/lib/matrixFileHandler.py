import importlib
import logging
# from .matrixFile import MatrixFile
log = logging.getLogger(__name__)


class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pFileType='cool', pMatrixFile=None, pChrnameList=None, pCooler_only_init=None,
                 pApplyCorrectionCooler=None, pBedFileHicPro=None, pCorrectionFactorTable=None,
                 pCorrectionOperator=None):

        self.class_ = getattr(importlib.import_module('.' + pFileType.lower(), package='hicexplorer.lib'), pFileType.title())
        if pFileType == 'cool' and pCooler_only_init:
            self.matrixFile = self.class_(pMatrixFile, pCooler_only_init)
            self.matrixFile.chrnameList = pChrnameList
        elif pFileType == 'hicpro':
            self.matrixFile = self.class_(pMatrixFile=pMatrixFile, pBedFile=pBedFileHicPro)
        else:
            self.matrixFile = self.class_(pMatrixFile)
            if pFileType == 'cool':
                self.matrixFile.chrnameList = pChrnameList
                if pCorrectionFactorTable is not None:
                    self.matrixFile.correctionFactorTable = pCorrectionFactorTable
                if pCorrectionOperator is not None:
                    self.matrixFile.correctionOperator = pCorrectionOperator

    def load(self):

        return self.matrixFile.load()

    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrixFile.set_matrix_variables(pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts)

    def save(self, pName, pSymmetric, pApplyCorrection):
        self.matrixFile.save(pName, pSymmetric, pApplyCorrection)

    def load_init(self):
        pass
