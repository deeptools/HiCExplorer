
import sys
import importlib
import logging
# from .matrixFile import MatrixFile
log = logging.getLogger(__name__)

class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pFileType='Cool', pMatrixFile=None, pChrnameList=None, pCooler_only_init=None, 
                    pApplyCorrectionCooler=None):
        
        self.class_ = getattr(importlib.import_module('.' + pFileType.lower(), package='hicexplorer.lib'), pFileType)
        if pFileType == 'Cool' and pCooler_only_init:
            self.matrixFile = self.class_(pMatrixFile, pCooler_only_init)
        else:
            self.matrixFile = self.class_(pMatrixFile)

    def load(self):

        return self.matrixFile.load()
        
    
    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrixFile.set_matrix_variables(pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts)

    def save(self, pName, pSymmetric, pApplyCorrection):
        self.matrixFile.save(pName, pSymmetric, pApplyCorrection)
        
    def load_init(self):
        pass    

    # def __check_file_format(self, pMatrix, p):
        
    
    # def __check_cooler(pFileName):
    #     if pFileName.endswith('.cool') or cooler.io.is_cooler(pFileName) or'.mcool' in pFileName:
    #         return True
    #     return False


    # def __convertNansToOnes(pArray):
    #     nan_elements = np.flatnonzero(np.isnan(pArray))
    #     if len(nan_elements) > 0:
    #         pArray[nan_elements] = 1.0
    #     return pArray
