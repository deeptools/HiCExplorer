from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer import HiCMatrix as hm
from tempfile import NamedTemporaryFile, mkdtemp
import shutil
import os
import numpy.testing as nt


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"


def test_build_matrix():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} -o {} -bs 5000 -b /tmp/test.bam --QCfolder {}".format(sam_R1, sam_R2,
                                                                           outfile.name,
                                                                           qc_folder).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT + "small_test_matrix.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    assert test.cut_intervals == new.cut_intervals

    print set(os.listdir(ROOT + "QC/"))
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)