
import sys
import importlib
import logging
log = logging.getLogger(__name__)

class MatrixFileHandler():
    """
    This class handles the load and save of the different Hi-C contact matrix formats.
    """

    def __init__(self, pMatrixFile, pFileType='Cool', pChrnameList=None, pCooler_only_init=None, 
                    pIntraChromosomalOnly=None, pApplyCorrectionCooler=None):
        
        self.__class = getattr(importlib.import_module('.' + pFileType.lower(), package='hicexplorer.matrixFileHandlerLib'), pFileType)
        self.__matrixFile = self.__class(pMatrixFile)

    def load(self):

        return self.__matrixFile.load()
        
    def load_init(self):
        pass    
    
    def save(self, pName):
       self.__matrixFile.save(pName)
        
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
