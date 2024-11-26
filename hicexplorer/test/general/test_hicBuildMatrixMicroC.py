import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicBuildMatrixMicroC, hicInfo
from hicmatrix import HiCMatrix as hm
from tempfile import NamedTemporaryFile, mkdtemp
import shutil
import os
import numpy.testing as nt
import pytest
from hicexplorer.test.test_compute_function import compute

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"
delta = 80000


def are_files_equal(file1, file2, delta=1):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        equal = False
                        break
                else:
                    equal = False
                    break
    return equal


def test_build_matrix(capsys):
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4".format(sam_R1, sam_R2,
                                                                                       outfile.name, outfile_bam.name,
                                                                                       qc_folder).split()
    compute(hicBuildMatrixMicroC.main, args, 5)
    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel_one_rc.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)
    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT + "QC_no_restriction_site/")))
    assert are_files_equal(ROOT + "QC_no_restriction_site/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT + "QC_no_restriction_site/")) == set(os.listdir(qc_folder))

    # accept delta of 80 kb, file size is around 4.5 MB
    assert abs(os.path.getsize(ROOT + "small_test_matrix_result.bam") - os.path.getsize(outfile_bam.name)) < delta

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
