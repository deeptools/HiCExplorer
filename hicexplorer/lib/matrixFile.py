import logging
log = logging.getLogger(__name__)

class MatrixFile():

    def __init__(self, pMatrixFileName):
        self._matrixFileName = pMatrixFileName

    def load(self):
        pass
    
    def save(self):
        pass
    
    def is_of_type(self):
        return False