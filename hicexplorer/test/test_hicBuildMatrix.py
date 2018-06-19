from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer import HiCMatrix as hm
from tempfile import NamedTemporaryFile, mkdtemp
import shutil
import os
import numpy.testing as nt


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/hicBuildMatrix/"
sam_R1 = ROOT + "R1_1000.bam"
sam_R2 = ROOT + "R2_1000.bam"
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
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    ROOT_DEFAULT = ROOT + 'default/'

    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 100000 -b /tmp/test.bam --QCfolder {} --threads 4".format(sam_R1, sam_R2,
                                                                                                  outfile.name,
                                                                                                  qc_folder).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT_DEFAULT + "R1_R2.cool")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)
    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT_DEFAULT + "QC/")))
    assert are_files_equal(ROOT_DEFAULT + "QC/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT_DEFAULT + "QC/")) == set(os.listdir(qc_folder))
    # accept delta of 60 kb, file size is around 4.5 MB
    assert abs(os.path.getsize(ROOT_DEFAULT + "R1_R2.bam") - os.path.getsize("/tmp/test.bam")) < 64000

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
    os.unlink("/tmp/test.bam")


# def test_build_matrix_cooler():
#     outfile = NamedTemporaryFile(suffix='.cool', delete=False)
#     outfile.close()
#     qc_folder = mkdtemp(prefix="testQC_")
#     args = "-s {} {} --outFileName {} -bs 5000 -b /tmp/test.bam --QCfolder {} --threads 4".format(sam_R1, sam_R2,
#                                                                                                   outfile.name,
#                                                                                                   qc_folder).split()
#     hicBuildMatrix.main(args)

#     test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel.h5")
#     new = hm.hiCMatrix(outfile.name)

#     nt.assert_equal(test.matrix.data, new.matrix.data)
#     # nt.assert_equal(test.cut_intervals, new.cut_intervals)
#     nt.assert_equal(len(new.cut_intervals), len(test.cut_intervals))
#     cut_interval_new_ = []
#     cut_interval_test_ = []
#     for x in new.cut_intervals:
#         cut_interval_new_.append(x[:3])
#     for x in test.cut_intervals:
#         cut_interval_test_.append(x[:3])

#     nt.assert_equal(cut_interval_new_, cut_interval_test_)
#     # print(set(os.listdir(ROOT + "QC/")))
#     assert are_files_equal(ROOT + "QC/QC.log", qc_folder + "/QC.log")
#     assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

#     os.unlink(outfile.name)
#     shutil.rmtree(qc_folder)


def test_build_matrix_rf():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    ROOT_RF = ROOT + 'restrictionFile/'
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} -rs {} --outFileName {}  --QCfolder {} " \
           "--restrictionSequence GATC " \
           "--danglingSequence GATC " \
           "--minDistance 150 " \
           "--maxLibraryInsertSize 1500 --threads 4".format(sam_R1, sam_R2, dpnii_file,
                                                            outfile.name,
                                                            qc_folder).split()
    hicBuildMatrix.main(args)

    test = hm.hiCMatrix(ROOT_RF + "R1_R2_restrictionFile.cool")
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(set(os.listdir(ROOT_RF + "QC/")))
    assert are_files_equal(ROOT_RF + "QC/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT_RF + "QC/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
