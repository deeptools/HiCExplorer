from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer import HiCMatrix as hm
from tempfile import NamedTemporaryFile, mkdtemp
import shutil
import os
import numpy.testing as nt


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"
dpnii_file = ROOT + "DpnII.bed"


def are_files_equal(file1, file2):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                equal = False
                break
    return equal


def test_build_matrix():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b /tmp/test.bam --QCfolder {} --threads 4".format(sam_R1, sam_R2,
                                                                                                  outfile.name,
                                                                                                  qc_folder).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT + "small_test_matrix.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print set(os.listdir(ROOT + "QC/"))
    assert are_files_equal(ROOT + "QC/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))
    assert abs(os.path.getsize(ROOT + "small_test_matrix_result.bam") - os.path.getsize("/tmp/test.bam")) < 1000
    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
    os.unlink("/tmp/test.bam")


def test_build_matrix_rf():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} -rs {} --outFileName {}  --QCfolder {} " \
           "--restrictionSequence GATC " \
           "--danglingSequence GATC " \
           "--minDistance 150 " \
           "--maxDistance 1500".format(sam_R1, sam_R2, dpnii_file,
                                       outfile.name,
                                       qc_folder).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT + "small_test_rf_matrix.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print set(os.listdir(ROOT + "QC_rc/"))
    assert are_files_equal(ROOT + "QC_rc/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT + "QC_rc/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
