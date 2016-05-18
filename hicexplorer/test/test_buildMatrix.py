from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer import HiCMatrix as hm

import os.path
from os import unlink
import numpy.testing as nt


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"


def test_build_matrix():
    outfile = '/tmp/matrix.h5'
    args = "-s {} {} -o {} -bs 5000 -b /tmp/test.bam".format(sam_R1, sam_R2, outfile).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT + "small_test_matrix.h5")
    new = hm.hiCMatrix(outfile)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    assert test.cut_intervals == new.cut_intervals
    unlink(outfile)