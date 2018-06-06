
from .matrixFile import MatrixFile
import logging
log = logging.getLogger(__name__)
class Ginteractions(MatrixFile):

    def __init__(self, pMatrixFile):
        super().__init__(pMatrixFile)
    

    def __save_GInteractions(self, fileName):
        self.restoreMaskedBins()
        log.debug(self.matrix.shape)
        mat_coo = triu(self.matrix, k=0, format='csr').tocoo()
        fileh = open("{}.tsv".format(fileName), 'w')
        for idx, counts in enumerate(mat_coo.data):
            chr_row, start_row, end_row, _ = self.cut_intervals[mat_coo.row[idx]]
            chr_col, start_col, end_col, _ = self.cut_intervals[mat_coo.col[idx]]
            fileh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_row, int(start_row), int(end_row),
                                                              chr_col, int(start_col), int(end_col), counts))
        fileh.close()