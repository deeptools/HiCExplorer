from hicexplorer import hicFindTADs
from hicexplorer import HiCMatrix as hm
from tempfile import mkdtemp
import shutil
import os
import numpy.testing as nt


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def are_files_equal(file1, file2, pDifference=10):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                count = sum(1 for a, b in zip(x, y) if a != b)
                if count > pDifference:
                    equal = False
                    break
    return equal


def test_find_TADs_fdr():

    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_fdr")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 4 --step 20000 \
    --outPrefix {}/test_multiFDR --minBoundaryDistance 20000 \
    --multipleComparisons fdr --thresholdComparisons 0.1".format(matrix, tad_folder).split()

    hicFindTADs.main(args)

    new = hm.hiCMatrix(tad_folder + "/test_multiFDR_zscore_matrix.h5")
    test = hm.hiCMatrix(ROOT + 'find_TADs/multiFDR_zscore_matrix.h5')
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(tad_folder + "/test_multiFDR_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/multiFDR_boundaries.bed", tad_folder + "/test_multiFDR_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/multiFDR_domains.bed", tad_folder + "/test_multiFDR_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/multiFDR_tad_score.bm", tad_folder + "/test_multiFDR_tad_score.bm")
    assert are_files_equal(ROOT + "find_TADs/multiFDR_boundaries.gff", tad_folder + "/test_multiFDR_boundaries.gff")
    assert are_files_equal(ROOT + "find_TADs/multiFDR_score.bedgraph", tad_folder + "/test_multiFDR_score.bedgraph")

    shutil.rmtree(tad_folder)


def test_find_TADs_bonferroni():

    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_bonferroni")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 4 --step 20000 \
    --outPrefix {}/test_multiBonferroni --minBoundaryDistance 20000 \
    --multipleComparisons bonferroni --thresholdComparisons 0.1".format(matrix, tad_folder).split()

    hicFindTADs.main(args)

    new = hm.hiCMatrix(tad_folder + "/test_multiBonferroni_zscore_matrix.h5")
    test = hm.hiCMatrix(ROOT + 'find_TADs/multiBonferroni_zscore_matrix.h5')
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(tad_folder + "/test_multiBonferroni_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/multiBonferroni_boundaries.bed", tad_folder + "/test_multiBonferroni_boundaries.bed", pDifference=10)
    assert are_files_equal(ROOT + "find_TADs/multiBonferroni_domains.bed", tad_folder + "/test_multiBonferroni_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/multiBonferroni_tad_score.bm", tad_folder + "/test_multiBonferroni_tad_score.bm")
    assert are_files_equal(ROOT + "find_TADs/multiBonferroni_boundaries.gff", tad_folder + "/test_multiBonferroni_boundaries.gff")
    assert are_files_equal(ROOT + "find_TADs/multiBonferroni_score.bedgraph", tad_folder + "/test_multiBonferroni_score.bedgraph")

    shutil.rmtree(tad_folder)


def test_find_TADs_none():

    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_none")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 4 --step 20000 \
    --outPrefix {}/test_multiNone --minBoundaryDistance 20000 \
    --multipleComparisons None".format(matrix, tad_folder).split()

    hicFindTADs.main(args)

    new = hm.hiCMatrix(tad_folder + "/test_multiNone_zscore_matrix.h5")
    test = hm.hiCMatrix(ROOT + 'find_TADs/multiNone_zscore_matrix.h5')
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    assert are_files_equal(ROOT + "find_TADs/multiNone_boundaries.bed", tad_folder + "/test_multiNone_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/multiNone_domains.bed", tad_folder + "/test_multiNone_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/multiNone_tad_score.bm", tad_folder + "/test_multiNone_tad_score.bm")
    assert are_files_equal(ROOT + "find_TADs/multiNone_boundaries.gff", tad_folder + "/test_multiNone_boundaries.gff")
    assert are_files_equal(ROOT + "find_TADs/multiNone_score.bedgraph", tad_folder + "/test_multiNone_score.bedgraph")

    shutil.rmtree(tad_folder)
