import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicDetectLongRangeContacts
import pytest
from psutil import virtual_memory
import numpy.testing as nt
import numpy as np

from scipy.sparse import csr_matrix
mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True
# DIFF = 60


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
tolerance = 13  # default matplotlib pixed difference tolerance


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_compute_zscore_matrix():
    # outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    # outfile.close()

    instances = np.array([0, 1, 1, 1, 2, 4, 6, 6, 6])
    features = np.array([0, 1, 2, 3, 1, 4, 1, 2, 7])
    data = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9])

    matrix = csr_matrix((data, (instances, features)))
    distances = np.absolute(instances-features)
    sum_per_distance_test = np.zeros(8)
    sum_per_distance_test = hicDetectLongRangeContacts._sum_per_distance(sum_per_distance_test, data, distances, 9)
    # sum per distance
    sum_per_distance_expected = np.array([8, 17, 4 ,0, 8, 7, 0, 0])
    
    nt.assert_equal(sum_per_distance_test, sum_per_distance_expected)    

    z_score = hicDetectLongRangeContacts.compute_zscore_matrix(matrix)
    print('foo:', z_score)