
import sys
import importlib
import logging
# from .matrixFile import MatrixFile
log = logging.getLogger(__name__)

class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pFileType='cool', pMatrixFile=None, pChrnameList=None, pCooler_only_init=None, 
                    pApplyCorrectionCooler=None, pBedFileHicPro=None):
        
        self.class_ = getattr(importlib.import_module('.' + pFileType.lower(), package='hicexplorer.lib'), pFileType.title())
        if pFileType == 'cool' and pCooler_only_init:
            self.matrixFile = self.class_(pMatrixFile, pCooler_only_init)
        elif pFileType == 'hicpro':
            self.matrixFile = self.class_(pMatrixFile=pMatrixFile, pBedFile=pBedFileHicPro)
        else:
            self.matrixFile = self.class_(pMatrixFile)

    def load(self):

        return self.matrixFile.load()
        
    
    def set_matrix_variables(self, pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts):
        self.matrixFile.set_matrix_variables(pMatrix, pCutIntervals, pNanBins, pCorrectionFactors, pDistanceCounts)
        log.info('self.matrix {}'.format(self.matrixFile.matrix))
        log.info('self.nan_bins {}'.format(self.matrixFile.nan_bins))
        log.info('self.cut_intervals {}'.format(self.matrixFile.cut_intervals))
        log.info('self.correction_factors {}'.format(self.matrixFile.correction_factors)) 

    def save(self, pName, pSymmetric, pApplyCorrection):
        log.debug('Save data, call save function of specific implementatin')
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
