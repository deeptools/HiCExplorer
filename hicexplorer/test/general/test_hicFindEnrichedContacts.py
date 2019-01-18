""" Testsuite for hicFindEnrichedContacts """
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicFindEnrichedContacts as hicfec
from hicmatrix import HiCMatrix as hm
import os
import numpy as np
import numpy.testing as nt
from scipy.sparse import csr_matrix

# ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
# original_matrix = ROOT + "small_test_matrix_50kb_res"

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"
dpnii_file = ROOT + "DpnII.bed"

from psutil import virtual_memory
mem = virtual_memory()
memory = mem.total / 2**30


# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200


def test_getPearson():
    cut_intervals = [('a', 0, 10, 1), ('a', 10, 20, 1),
                     ('a', 20, 30, 1), ('a', 30, 40, 1), ('b', 40, 50, 1)]
    hic = hm.hiCMatrix()
    hic.nan_bins = []
    matrix = np.array([[1, 8, 5, 3, 0],
                       [0, 4, 15, 5, 1],
                       [0, 0, 0, 0, 2],
                       [0, 0, 0, 0, 1],
                       [0, 0, 0, 0, 0]])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    pearson = hicfec.getPearson(hic.matrix)

    nt.assert_equal(pearson.shape, (5, 5))
    for row in pearson:
        for value in row:
            assert -1 <= value.all() <= 1

    hic = hm.hiCMatrix()
    hic.nan_bins = []
    matrix = np.array([[2, 2, 3, 4, 5],
                       [1, 2, 3, 4, 5],
                       [1, 2, 2, 4, 5],
                       [1, 2, 3, 2, 5],
                       [1, 2, 3, 4, 2]
                       ])

    hic.matrix = csr_matrix(matrix)
    hic.setMatrix(hic.matrix, cut_intervals)

    pearson = hicfec.getPearson(hic.matrix)

    nt.assert_equal(pearson.shape, (5, 5))
    nt.assert_equal(pearson, np.zeros(shape=pearson.shape))
    for row in pearson:
        for value in row:
            assert -1 <= value.all() <= 1
