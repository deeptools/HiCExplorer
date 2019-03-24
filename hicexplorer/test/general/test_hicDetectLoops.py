import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicDetectLoops
import pytest
from psutil import virtual_memory
import numpy.testing as nt
import numpy as np

from scipy.sparse import csr_matrix
mem = virtual_memory()
memory = mem.total / 2**30

import logging
log = logging.getLogger(__name__)

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True
# DIFF = 60

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")

tolerance = 13  # default matplotlib pixed difference tolerance


@pytest.mark.skipif(MID_MEMORY > memory,
                    reason="Travis has too less memory to run it.")
def test_sum_per_distance():
    # outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    # outfile.close()

    instances = np.array([0, 1, 1, 1, 2, 4, 6, 6, 6])
    features = np.array([0, 1, 2, 3, 1, 4, 1, 2, 7])
    data = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9])

    matrix = csr_matrix((data, (instances, features)))
    distances = np.absolute(instances - features)
    sum_per_distance_test = np.zeros(8)
    sum_per_distance_test, distance_count_test = hicDetectLoops._sum_per_distance(
        sum_per_distance_test, data, distances)
    # sum per distance
    sum_per_distance_expected = np.array([8, 17, 4, 0, 8, 7, 0, 0])
    distance_count_expected = np.array([3, 3, 1, 0, 1, 1, 0, 0])
    nt.assert_equal(sum_per_distance_test, sum_per_distance_expected)
    nt.assert_equal(distance_count_test, distance_count_expected)

    # z_score = hicDetectLoops.compute_zscore_matrix(matrix)
    # print('foo:', z_score)


def test_compute_zscore_matrix():
    instances = np.array([0, 1, 1, 1, 2, 4, 6, 6, 6])
    features = np.array([0, 1, 2, 3, 1, 4, 1, 2, 7])
    data = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9])

    # z_score_test = hicDetectLoops.compute_zscore_matrix(instances, features, data, 8)

    # z_score_expected = np.array([])
    # log.debug('{}'.format(z_score))


def test_filter_duplicates():
    pass


def test_candidate_region_test():
    pass


def test_cluster_to_genome_position_mapping():
    pass


def test_write_bedgraph():
    pass


def test_smoothInteractionValues():
    pass


def test_main_h5():
    outfile_loop_h5 = NamedTemporaryFile(suffix='bedgraph', delete=True)

    args = "--matrix {} -o {} -d 1e-10 -pit 1 -p 1".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)


def test_main_cool():
    # def test_correlate():
    outfile_loop_cool = NamedTemporaryFile(suffix='bedgraph', delete=True)

    # outfile_heatmap = NamedTemporaryFile(suffix='heatmap.png', prefix='hicexplorer_test', delete=False)
    # outfile_scatter = NamedTemporaryFile(suffix='scatter.png', prefix='hicexplorer_test', delete=False)

    args = "--matrix {} -o {} -d 1e-10 -pit 1 -p 1".format(
        ROOT + "small_test_matrix.cool", outfile_loop_cool.name).split()
    hicDetectLoops.main(args)

    # res = compare_images(ROOT + "hicCorrelate" + '/heatmap.png', outfile_heatmap.name, tol=40)
    # assert res is None, res

    # res = compare_images(ROOT + "hicCorrelate" + '/scatter.png', outfile_scatter.name, tol=40)
    # assert res is None, res
    # os.remove(outfile_heatmap.name)
    # os.remove(outfile_scatter.name)
