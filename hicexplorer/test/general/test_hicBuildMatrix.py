import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicBuildMatrix, hicInfo
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
dpnii_file = ROOT + "DpnII.bed"
# QC.log references written by the fixed HiCExplorer 4 tools (Q1 to Q6 in
# cpp/tools/hicBuildMatrix.cpp). The QC folders in test_data/QC* hold the
# 3.7.6 numbers and serve the QC report tools as input, so they stay as they
# are; only the file names of their content are compared below.
QC_V4 = ROOT + "hicBuildMatrix/hicx4_qc/"


def assert_out_bam_matches_the_contacts(bam_path, qc_log):
    """The records of --outBam are two per Hi-C contact, in the order of the
    input, and carry the insert size of the pair (Q5, fixed in HiCExplorer 4:
    3.7.6 wrote 0 in every record)."""
    import pysam
    contacts = None
    with open(qc_log) as handle:
        for line in handle:
            if line.startswith('Hi-C contacts'):
                contacts = int(line.split('\t')[1])
    with pysam.Samfile(bam_path, 'rb') as handle:
        records = list(handle)
    assert len(records) == 2 * contacts
    for first, second in zip(records[0::2], records[1::2]):
        assert first.qname == second.qname
        if first.reference_id == second.reference_id:
            assert first.isize == second.pos - first.pos
            assert second.isize == first.pos - second.pos


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
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 \
            --restrictionSequence GATC --danglingSequence GATC -rs {}".format(sam_R1, sam_R2,
                                                                              outfile.name, outfile_bam.name,
                                                                              qc_folder, dpnii_file).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)
    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel_one_rc.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)
    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT + "QC/")))
    assert are_files_equal(QC_V4 + "QC_one.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    assert_out_bam_matches_the_contacts(outfile_bam.name, qc_folder + "/QC.log")

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
    # os.unlink("/tmp/test.bam")


def test_build_matrix_restriction_enzyme(capsys):
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 --danglingSequence GATC AGCT --restrictionSequence GATC AAGCTT -rs {} {}".format(sam_R1, sam_R2,
                                                                                                                                                                outfile.name, outfile_bam.name,
                                                                                                                                                                qc_folder, dpnii_file, ROOT + 'hicFindRestSite/hindIII.bed').split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel_two_rc.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)
    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT + "QC_multi_restriction/")))
    assert are_files_equal(QC_V4 + "QC_multi_restriction.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC_multi_restriction/")) == set(os.listdir(qc_folder))

    assert_out_bam_matches_the_contacts(outfile_bam.name, qc_folder + "/QC.log")

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
    # os.unlink("/tmp/test.bam")


def test_build_matrix_restriction_enzyme_region(capsys):
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 --danglingSequence GATC AGCT --restrictionSequence GATC AAGCTT -rs {} {} --region {}".format(sam_R1, sam_R2,
                                                                                                                                                                            outfile.name, outfile_bam.name,
                                                                                                                                                                            qc_folder, dpnii_file, ROOT + 'hicFindRestSite/hindIII.bed', 'chr3R').split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel_two_rc_chr3R.h5")
    new = hm.hiCMatrix(outfile.name)

    # print('test.cut_intervals {}'.format(test.cut_intervals))
    # print('new.cut_intervals {}'.format(new.cut_intervals))
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(len(test.cut_intervals), len(new.cut_intervals))

    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT + "QC_region/")))
    assert are_files_equal(QC_V4 + "QC_region.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC_region/")) == set(os.listdir(qc_folder))

    assert_out_bam_matches_the_contacts(outfile_bam.name, qc_folder + "/QC.log")

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


def test_build_matrix_chromosome_sizes(capsys):
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 --chromosomeSizes {}  \
            --restrictionSequence GATC --danglingSequence GATC -rs {}".format(sam_R1, sam_R2,
                                                                              outfile.name, outfile_bam.name,
                                                                              qc_folder, ROOT + 'hicBuildMatrix/dm3.chrom.sizes', dpnii_file).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicBuildMatrix/chromosome_sizes/small_test_chromosome_size.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)
    # print("MATRIX NAME:", outfile.name)
    print(set(os.listdir(ROOT + "QC/")))
    assert are_files_equal(QC_V4 + "QC_one.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    assert_out_bam_matches_the_contacts(outfile_bam.name, qc_folder + "/QC.log")

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)
    # os.unlink("/tmp/test.bam")


def test_build_matrix_cooler():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4  \
            --restrictionSequence GATC --danglingSequence GATC -rs {}".format(sam_R1, sam_R2,
                                                                              outfile.name, outfile_bam.name,
                                                                              qc_folder, dpnii_file).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel.h5")
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(test.matrix.data, new.matrix.data)
    # nt.assert_equal(test.cut_intervals, new.cut_intervals)
    nt.assert_equal(len(new.cut_intervals), len(test.cut_intervals))
    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)
    # print(set(os.listdir(ROOT + "QC/")))
    assert are_files_equal(QC_V4 + "QC_one.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


def test_build_matrix_cooler_metadata():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 --genomeAssembly dm3  \
            --restrictionSequence GATC --danglingSequence GATC -rs {}".format(sam_R1, sam_R2,
                                                                              outfile.name, outfile_bam.name,
                                                                              qc_folder, dpnii_file).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_matrix_parallel.h5")
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(test.matrix.data, new.matrix.data)
    # nt.assert_equal(test.cut_intervals, new.cut_intervals)
    nt.assert_equal(len(new.cut_intervals), len(test.cut_intervals))
    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)
    # print(set(os.listdir(ROOT + "QC/")))
    assert are_files_equal(QC_V4 + "QC_one.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    outfile_metadata = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile_metadata.close()
    args = "-m {} -o {}".format(outfile.name, outfile_metadata.name).split()
    hicInfo.main(args)
    assert are_files_equal(ROOT + "hicBuildMatrix/metadata.txt", outfile_metadata.name, delta=7)

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


def test_build_matrix_cooler_multiple():
    outfile = NamedTemporaryFile(suffix='.mcool', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 10000 20000 -b {} --QCfolder {} --threads 4  \
            --restrictionSequence GATC --danglingSequence GATC -rs {}".format(sam_R1, sam_R2,
                                                                              outfile.name, outfile_bam.name,
                                                                              qc_folder, dpnii_file).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test_5000 = hm.hiCMatrix(ROOT + "hicBuildMatrix/multi_small_test_matrix_hicx4.mcool::/resolutions/5000")
    test_10000 = hm.hiCMatrix(ROOT + "hicBuildMatrix/multi_small_test_matrix_hicx4.mcool::/resolutions/10000")
    test_20000 = hm.hiCMatrix(ROOT + "hicBuildMatrix/multi_small_test_matrix_hicx4.mcool::/resolutions/20000")

    new_5000 = hm.hiCMatrix(outfile.name + '::/resolutions/5000')
    new_10000 = hm.hiCMatrix(outfile.name + '::/resolutions/10000')
    new_20000 = hm.hiCMatrix(outfile.name + '::/resolutions/20000')

    nt.assert_equal(test_5000.matrix.data, new_5000.matrix.data)
    nt.assert_equal(test_10000.matrix.data, new_10000.matrix.data)
    nt.assert_equal(test_20000.matrix.data, new_20000.matrix.data)

    # nt.assert_equal(test.cut_intervals, new.cut_intervals)
    nt.assert_equal(len(new_5000.cut_intervals), len(test_5000.cut_intervals))
    nt.assert_equal(len(new_10000.cut_intervals), len(test_10000.cut_intervals))
    nt.assert_equal(len(new_20000.cut_intervals), len(test_20000.cut_intervals))

    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new_5000.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test_5000.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)

    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new_10000.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test_10000.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)

    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new_20000.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test_20000.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)
    # print(set(os.listdir(ROOT + "QC/")))
    assert are_files_equal(QC_V4 + "QC_one.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


def test_build_matrix_rf():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} -rs {} --outFileName {}  --QCfolder {} " \
           "--restrictionSequence GATC " \
           "--danglingSequence GATC " \
           "--minDistance 150 " \
           "--maxLibraryInsertSize 1500 --threads 4".format(sam_R1, sam_R2, dpnii_file,
                                                            outfile.name,
                                                            qc_folder).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_rf_matrix.h5")
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(set(os.listdir(ROOT + "QC_rc/")))
    assert are_files_equal(QC_V4 + "QC_rc.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC_rc/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


def test_build_matrix_rf_multi():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} -rs {} {} --outFileName {}  --QCfolder {} " \
           "--restrictionSequence GATC AAGCTT " \
           "--danglingSequence GATC AGCT " \
           "--minDistance 150 " \
           "--maxLibraryInsertSize 1500 --threads 4".format(sam_R1, sam_R2, dpnii_file, ROOT + 'hicFindRestSite/hindIII.bed',
                                                            outfile.name,
                                                            qc_folder).split()
    # --danglingSequence GATC AGCT --restrictionSequence GATC AAGCTT
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)

    test = hm.hiCMatrix(ROOT + "small_test_rf_matrix_multiple.h5")
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(set(os.listdir(ROOT + "QC_rc_multiple/")))
    assert are_files_equal(QC_V4 + "QC_rc_multiple.log", qc_folder + "/QC.log", delta=0)
    assert set(os.listdir(ROOT + "QC_rc_multiple/")) == set(os.listdir(qc_folder))

    os.unlink(outfile.name)
    shutil.rmtree(qc_folder)


@pytest.mark.xfail
def test_build_matrix_fail(capsys):
    outfile = NamedTemporaryFile(suffix='', delete=False)
    outfile.close()
    outfile_bam = NamedTemporaryFile(suffix='.bam', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --outFileName {} -bs 5000 -b {} --QCfolder {} --threads 4 ".format(sam_R1, sam_R2, outfile_bam.name,
                                                                                        outfile.name,
                                                                                        qc_folder).split()
    # hicBuildMatrix.main(args)
    compute(hicBuildMatrix.main, args, 5)


# ---------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port (cpp/AGENTS_CONTRACT.md
# rule 1). What the tests above this line actually constrain is narrower than
# it looks:
#
#   * nt.assert_equal(test.matrix.data, new.matrix.data) compares the CSR value
#     array alone. Neither `indices` nor `indptr` is looked at, so a matrix
#     whose counts are correct but sit in the wrong cells passes.
#   * are_files_equal skips every line starting with "File", tolerates `delta`
#     mismatching lines (1 by default, 2 in test_build_matrix) and zips the two
#     files, so trailing lines present in only one of them are invisible. The
#     stored references under test_data/QC* are in fact stale: they say
#     "dangling end GATC" where the current tool writes "dangling end GATC
#     (restriction sequence GATC)", and the delta is what hides it.
#   * the QC folder is compared by file *name* set only, so the content of
#     QC_table.txt and the four derived tables is unchecked, although
#     hicPrepareQCreport, which has no test file at all, produces them.
#   * the output BAM is only checked by byte size within 80,000.
#   * --maxDistance, --keepSelfLigation, --keepSelfCircles, --minMappingQuality,
#     --skipDuplicationCheck, --doTestRun, --doTestRunLines and
#     --inputBufferSize are never asserted on.
#
# The tests below pin all of that at full precision on real data. Several of
# them pin a defect rather than a feature; those say so.
#
# The fixture is R1_1000.bam / R2_1000.bam against the HindIII cut sites, which
# is 983 read pairs against 60,516 sites and runs in about fourteen seconds,
# against the eighty a DpnII run on the 99,983 pair library costs.

import hashlib
import json
import numpy as np
import pysam

RS_HINDIII = ROOT + "hicFindRestSite/hindIII.bed"
R1_1000 = ROOT + "R1_1000.bam"
R2_1000 = ROOT + "R2_1000.bam"


def _small_fixture(tmp_dir, extra="", suffix=".h5", out_bam=None):
    """One hicBuildMatrix run on the small fixture. Returns (matrix, QCfolder)."""
    outfile = NamedTemporaryFile(suffix=suffix, delete=False, dir=tmp_dir)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_", dir=tmp_dir)
    args = ("-s {} {} --outFileName {} -bs 100000 --QCfolder {} --threads 4 "
            "--restrictionSequence AAGCTT --danglingSequence AGCT -rs {} {}").format(
        R1_1000, R2_1000, outfile.name, qc_folder, RS_HINDIII, extra)
    if out_bam is not None:
        args += " -b {}".format(out_bam)
    compute(hicBuildMatrix.main, args.split(), 5)
    return outfile.name, qc_folder


def _matrix_fingerprint(path):
    """Every number the matrix holds, at full precision.

    The sha256 covers `data`, `indices` and `indptr`, so unlike the assertions
    above it fails when a count lands in the wrong cell, not only when its
    value changes.
    """
    matrix = hm.hiCMatrix(path)
    csr = matrix.matrix.tocsr()
    digest = hashlib.sha256()
    for array in (csr.data.astype('<i8'), csr.indices.astype('<i8'),
                  csr.indptr.astype('<i8')):
        digest.update(array.tobytes())
    extra = np.array([interval[3] for interval in matrix.cut_intervals], dtype=float)
    return {
        "shape": csr.shape,
        "nnz": int(csr.nnz),
        "sum": int(csr.sum()),
        "dtype": str(csr.data.dtype),
        "sha256": digest.hexdigest(),
        "bins": len(matrix.cut_intervals),
        "nan_bins": int(np.isnan(extra).sum()),
        "coverage_sum": float(np.nansum(extra)),
        "coverage_max": float(np.nanmax(extra)),
        "boundaries": list(matrix.chrBinBoundaries.items()),
    }


def _qc_body(qc_folder):
    """QC.log without the line that carries the temporary output file name."""
    with open(qc_folder + "/QC.log") as handle:
        return "".join(line for line in handle if not line.startswith("File\t"))


def _table_body(qc_folder, name):
    """One QC table with the File column, which is a temporary path, removed."""
    rows = []
    with open(qc_folder + "/" + name) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            rows.append("\t".join(fields[1:]))
    return "\n".join(rows) + "\n"


# HiCExplorer 4 output. Against 3.7.6 the fixes of Q1 and Q2 change it: the 13
# pairs with both mates at quality 0 are not unique only since Q1 (3.7.6
# charged the pairs whose first mate had quality 0, which here are 13 pairs
# whose second mate is mappable, and moved them to "Low mapping quality" now),
# and the 20 self circles are removed (Q2), so 195 instead of 215 contacts, 20
# fewer outward pairs and 20 fewer short range pairs.
EXPECTED_QC_LOG = """\
Sequenced reads\t983\t\t
Min rest. site distance\t300\t\t
Max library insert size\t1000\t\t

#\tcount\t(percentage w.r.t. total sequenced reads)
Pairs mappable, unique and high quality\t300\t(30.52)
Hi-C contacts\t195\t(19.84)
One mate unmapped\t506\t(51.48)
One mate not unique\t0\t(0.00)
Low mapping quality\t177\t(18.01)

#\tcount\t(percentage w.r.t. mappable, unique and high quality pairs)
dangling end AGCT (restriction sequence AAGCTT)\t1\t(0.33)
self ligation (removed)\t5\t(1.67)
One mate not close to rest site\t0\t(0.00)
same fragment\t79\t(26.33)
self circle\t20\t(6.67)
duplicated pairs\t0\t(0.00)

#\tcount\t(percentage w.r.t. total valid pairs used)
inter chromosomal\t24\t(12.31)
Intra short range (< 20kb)\t41\t(21.03)
Intra long range (>= 20kb)\t130\t(66.67)
Read pair type: inward pairs\t42\t(21.54)
Read pair type: outward pairs\t46\t(23.59)
Read pair type: left pairs\t48\t(24.62)
Read pair type: right pairs\t35\t(17.95)
"""


def test_build_matrix_pins_the_whole_matrix_and_bin_table(tmp_path):
    """The counts, where they sit, and the per bin coverage the tool derives."""
    matrix_file, qc_folder = _small_fixture(str(tmp_path))
    fingerprint = _matrix_fingerprint(matrix_file)

    assert fingerprint["shape"] == (1697, 1697)
    assert fingerprint["bins"] == 1697
    assert fingerprint["nnz"] == 325
    assert fingerprint["sum"] == 332
    assert fingerprint["dtype"] == "int64"
    # The whole CSR, indices and row offsets included.
    assert fingerprint["sha256"] == (
        "0b0b0f3486c085fb66ff24182a9c5e41a9508d7ffc96a141f42926dc52ae0f9d")
    # The fourth element of every cut interval is the per bin maximum of the
    # coverage vector, or NaN where the bin saw no read. 1,414 of the 1,697
    # bins are empty and the ones that are not sum to 284 with a maximum of 2.
    # (Q3 and Q4, fixed in HiCExplorer 4: each mate is counted in its own bin
    # and the last cell of a bin is counted.)
    assert fingerprint["nan_bins"] == 1414
    assert fingerprint["coverage_sum"] == 284.0
    assert fingerprint["coverage_max"] == 2.0
    assert fingerprint["boundaries"][:3] == [
        ('chr2L', (0, 231)), ('chr2LHet', (231, 235)), ('chr2R', (235, 447))]

    matrix = hm.hiCMatrix(matrix_file)
    assert matrix.cut_intervals[0][:3] == ('chr2L', 0, 100000)
    assert np.isnan(matrix.cut_intervals[0][3])
    assert matrix.cut_intervals[-1][:3] == ('chrYHet', 300000, 347038)
    shutil.rmtree(qc_folder)


def test_build_matrix_pins_the_qc_log_exactly(tmp_path):
    """QC.log is the input of hicPrepareQCreport, so its text is a contract."""
    _, qc_folder = _small_fixture(str(tmp_path))
    body = _qc_body(qc_folder)
    # The log opens with an empty line, which the reference above does not
    # repeat because splitting on it would be noise.
    assert body.startswith("\n")
    assert body[1:] == EXPECTED_QC_LOG
    shutil.rmtree(qc_folder)


def test_build_matrix_pins_the_qc_tables(tmp_path):
    """The five tables hicPrepareQCreport writes, at full float precision.

    hicPrepareQCreport has no test file of its own, and hicBuildMatrix calls
    it, so these are the only assertions in the suite on its output. The
    percentages are pandas float reprs, that is the shortest string that round
    trips, and they are compared as text.
    """
    _, qc_folder = _small_fixture(str(tmp_path))

    assert _table_body(qc_folder, "QC_table.txt") == (
        "Sequenced reads\tMin rest. site distance\tMax library insert size\t"
        "Pairs mappable, unique and high quality\tHi-C contacts\t"
        "One mate unmapped\tOne mate not unique\tLow mapping quality\t"
        "dangling end AGCT (restriction sequence AAGCTT)\t"
        "self ligation (removed)\tOne mate not close to rest site\t"
        "same fragment\tself circle\tduplicated pairs\tinter chromosomal\t"
        "Intra short range (< 20kb)\tIntra long range (>= 20kb)\t"
        "Read pair type: inward pairs\tRead pair type: outward pairs\t"
        "Read pair type: left pairs\tRead pair type: right pairs\n"
        "983\t300\t1000\t300\t195\t506\t0\t177\t1\t5\t0\t79\t20\t0\t24\t41\t"
        "130\t42\t46\t48\t35\n")

    assert _table_body(qc_folder, "unmapable_table.txt") == (
        "Hi-C contacts\tHi-C contacts_%\tLow mapping quality\t"
        "Low mapping quality_%\tOne mate not unique\tOne mate not unique_%\t"
        "One mate unmapped\tOne mate unmapped_%\n"
        "195\t0.19837232960325535\t177\t0.18006103763987794\t0\t"
        "0.0\t506\t0.5147507629704985\n")

    # 'dangling end' is not a column name, only a prefix of one, so
    # make_figure_pairs_discarded emits the percentage for it and no count.
    assert _table_body(qc_folder, "discarded_table.txt") == (
        "One mate not close to rest site\tOne mate not close to rest site %\t"
        "dangling end AGCT (restriction sequence AAGCTT) %\t"
        "duplicated pairs\tduplicated pairs %\tsame fragment\tsame fragment %\t"
        "self circle\tself circle %\tself ligation (removed)\t"
        "self ligation (removed) %\n"
        "0\t0.0\t0.0033333333333333335\t0\t0.0\t79\t0.2633333333333333\t20\t"
        "0.06666666666666667\t5\t0.016666666666666666\n")

    assert _table_body(qc_folder, "distance_table.txt") == (
        "inter chromosomal\tinter chromosomal %\tIntra short range (< 20kb)\t"
        "Intra short range (< 20kb) %\tIntra long range (>= 20kb)\t"
        "Intra long range (>= 20kb) %\n"
        "24\t0.12307692307692308\t41\t0.21025641025641026\t130\t"
        "0.6666666666666666\n")

    # These percentages are relative to the sum of the four orientations, not
    # to the contact count, so they add up to one.
    assert _table_body(qc_folder, "read_orientation_table.txt") == (
        "Read pair type: inward pairs\tRead pair type: inward pairs %\t"
        "Read pair type: outward pairs\tRead pair type: outward pairs %\t"
        "Read pair type: left pairs\tRead pair type: left pairs %\t"
        "Read pair type: right pairs\tRead pair type: right pairs %\n"
        "42\t0.24561403508771928\t46\t0.26900584795321636\t48\t"
        "0.2807017543859649\t35\t0.2046783625730994\n")
    shutil.rmtree(qc_folder)


def test_build_matrix_keep_self_circles(tmp_path):
    """HiCExplorer 4 (Q2): --keepSelfCircles decides about the self circles.

    buildMatrixMethods.py:759-762 counts a self circle and then, when
    --keepSelfCircles is off, executes `continue`. That `continue` belongs to
    the enclosing `for restrictionSequence` loop, not to the loop over read
    pairs, so 3.7.6 never dropped the pair and the flag changed nothing. Now
    the 20 self circles of the fixture are removed, and kept with the flag: the
    matrix with the flag is the matrix 3.7.6 wrote in either case.
    """
    without_file, without_qc = _small_fixture(str(tmp_path))
    with_file, with_qc = _small_fixture(str(tmp_path), extra="--keepSelfCircles")

    without = _matrix_fingerprint(without_file)
    with_flag = _matrix_fingerprint(with_file)
    assert without["sum"] == 332
    assert with_flag["sum"] == 353
    # the pixels of the 3.7.6 matrix (its coverage column differs, Q3 and Q4)
    assert with_flag["sha256"] == (
        "da5aed0d813d1d472ef041794f4065ec545f29b69ada3ccb43762f2dc20baaf1")
    assert "self circle\t20\t(6.67)\n" in _qc_body(without_qc)
    assert "self circle\t20\t(6.67)\n" in _qc_body(with_qc)
    assert "Hi-C contacts\t195\t(19.84)\n" in _qc_body(without_qc)
    assert "Hi-C contacts\t215\t(21.87)\n" in _qc_body(with_qc)
    shutil.rmtree(without_qc)
    shutil.rmtree(with_qc)


def test_build_matrix_two_restriction_sequences_count_a_self_circle_once(tmp_path):
    """HiCExplorer 4 (Q2): with two restriction sequences 3.7.6 counted the
    same self circle once per sequence and kept the pair."""
    def run(extra):
        outfile = NamedTemporaryFile(suffix='.h5', delete=False, dir=str(tmp_path))
        outfile.close()
        qc_folder = mkdtemp(prefix="testQC_", dir=str(tmp_path))
        args = ("-s {} {} --outFileName {} -bs 100000 --QCfolder {} --threads 4 "
                "--restrictionSequence AAGCTT GATC --danglingSequence AGCT GATC "
                "-rs {} {} {}").format(R1_1000, R2_1000, outfile.name, qc_folder,
                                       RS_HINDIII, dpnii_file, extra)
        compute(hicBuildMatrix.main, args.split(), 5)
        return _qc_body(qc_folder)

    removed = run("")
    kept = run("--keepSelfCircles")
    assert "self circle\t2\t(0.67)\n" in removed
    assert "self circle\t2\t(0.67)\n" in kept
    assert "Hi-C contacts\t213\t(21.67)\n" in removed
    assert "Hi-C contacts\t215\t(21.87)\n" in kept


def test_build_matrix_keep_self_ligation_adds_the_self_ligated_pairs(tmp_path):
    """--keepSelfLigation, which no test exercised, does work.

    Five pairs are self ligations here, and keeping them raises the contacts
    from 195 to 200 and moves all five into the short range and inward
    buckets. The QC label changes with it.
    """
    matrix_file, qc_folder = _small_fixture(str(tmp_path), extra="--keepSelfLigation")
    fingerprint = _matrix_fingerprint(matrix_file)
    assert fingerprint["nnz"] == 329
    assert fingerprint["sum"] == 337
    assert fingerprint["sha256"] == (
        "4320b9791a7de2fb9eea0190a3a474826c3caef34f594c6c32fb1c77d7cd4505")
    assert fingerprint["nan_bins"] == 1412
    assert fingerprint["coverage_sum"] == 286.0

    body = _qc_body(qc_folder)
    assert "self ligation (not removed)\t5\t(1.67)\n" in body
    assert "Hi-C contacts\t200\t(20.35)\n" in body
    assert "Intra short range (< 20kb)\t46\t(23.00)\n" in body
    assert "Read pair type: inward pairs\t47\t(23.50)\n" in body
    shutil.rmtree(qc_folder)


def test_build_matrix_max_distance_overrides_max_library_insert_size(tmp_path):
    """--maxDistance is the obsolete spelling and simply replaces the other."""
    obsolete_file, obsolete_qc = _small_fixture(str(tmp_path),
                                                extra="--maxDistance 1500")
    current_file, current_qc = _small_fixture(str(tmp_path),
                                              extra="--maxLibraryInsertSize 1500")
    assert _matrix_fingerprint(obsolete_file) == _matrix_fingerprint(current_file)
    assert _qc_body(obsolete_qc) == _qc_body(current_qc)
    assert "Max library insert size\t1500\t\t\n" in _qc_body(obsolete_qc)
    shutil.rmtree(obsolete_qc)
    shutil.rmtree(current_qc)


def test_build_matrix_min_mapping_quality(tmp_path):
    """--minMappingQuality, which no test asserted on.

    Raising it from the default 15 to 30 moves 71 pairs out of the mappable
    set. HiCExplorer 4 (Q1) counts a pair as "One mate not unique" when the
    quality of both mates is 0 (buildMatrixMethods.py:523 wrote `mate1.mapq ==
    0 & mate2.mapq == 0`, a chained comparison equal to `mate1.mapq == 0`, and
    3.7.6 charged 13 pairs to it here). No pair of the fixture has two mates at
    quality 0, so the 13 pairs are "Low mapping quality" now and the count is
    248 instead of 235.
    """
    matrix_file, qc_folder = _small_fixture(str(tmp_path),
                                            extra="--minMappingQuality 30")
    fingerprint = _matrix_fingerprint(matrix_file)
    assert fingerprint["nnz"] == 274
    assert fingerprint["sum"] == 279
    assert fingerprint["sha256"] == (
        "135bd57eef7d12d9292d24d9993dd67fd52d9cfe3d6807f0e627b1ee2682d4c2")

    body = _qc_body(qc_folder)
    assert "Pairs mappable, unique and high quality\t229\t(23.30)\n" in body
    assert "Hi-C contacts\t164\t(16.68)\n" in body
    assert "One mate not unique\t0\t(0.00)\n" in body
    assert "Low mapping quality\t248\t(25.23)\n" in body
    shutil.rmtree(qc_folder)


def test_build_matrix_skip_duplication_check(tmp_path):
    """--skipDuplicationCheck, which no test asserted on.

    R1_1000 holds no duplicate pair, so the run is identical to the checked
    one. What that pins is that skipping the check perturbs nothing else, and
    that the duplicate counter really is zero rather than merely unreported.
    """
    checked_file, checked_qc = _small_fixture(str(tmp_path))
    skipped_file, skipped_qc = _small_fixture(str(tmp_path),
                                              extra="--skipDuplicationCheck")
    assert _matrix_fingerprint(checked_file) == _matrix_fingerprint(skipped_file)
    assert _qc_body(checked_qc) == _qc_body(skipped_qc)
    assert "duplicated pairs\t0\t(0.00)\n" in _qc_body(checked_qc)
    shutil.rmtree(checked_qc)
    shutil.rmtree(skipped_qc)


def test_build_matrix_do_test_run_writes_qc_but_no_matrix(tmp_path):
    """--doTestRun and --doTestRunLines, which no test asserted on.

    The quick QC mode fills the counters and writes the whole QC folder but
    never fills the pixel arrays and never saves a matrix. The output file the
    argument parser created is unlinked and not written back, so nothing is
    left at that path at all.
    """
    matrix_file, qc_folder = _small_fixture(
        str(tmp_path), extra="--doTestRun --doTestRunLines 500")
    assert not os.path.exists(matrix_file)
    assert _qc_body(qc_folder)[1:] == EXPECTED_QC_LOG
    assert set(os.listdir(qc_folder)) >= {
        "QC.log", "QC_table.txt", "unmapable_table.txt", "discarded_table.txt",
        "distance_table.txt", "read_orientation_table.txt"}
    shutil.rmtree(qc_folder)


def test_build_matrix_input_buffer_size_crashes_on_an_exact_multiple(tmp_path):
    """PINNED DEFECT, not a feature.

    readBamFiles returns (None, None, True, ...) when a call finds no accepted
    pair left, which happens whenever the number of accepted pairs is an exact
    multiple of --inputBufferSize. createMatrix then evaluates
    `len(buffer_workers1[i])` on that None at buildMatrixMethods.py:1083 and
    dies with a TypeError. This fixture accepts exactly 300 pairs, so a buffer
    of 50 reproduces it and a buffer of 40 does not.
    """
    with pytest.raises(Exception) as failure:
        _small_fixture(str(tmp_path), extra="--inputBufferSize 50")
    assert "NoneType" in str(failure.value)

    # The same run with a buffer that does not divide 300 completes and gives
    # the same answer as the default buffer, so --inputBufferSize is otherwise
    # invisible in the output.
    buffered_file, buffered_qc = _small_fixture(str(tmp_path),
                                                extra="--inputBufferSize 40")
    default_file, default_qc = _small_fixture(str(tmp_path))
    assert _matrix_fingerprint(buffered_file) == _matrix_fingerprint(default_file)
    assert _qc_body(buffered_qc) == _qc_body(default_qc)
    shutil.rmtree(buffered_qc)
    shutil.rmtree(default_qc)


def test_build_matrix_out_bam_contents(tmp_path):
    """The output BAM, which the tests above check only by byte size.

    Every accepted pair is written as two records in input order, the header is
    the first input's header verbatim, and four fields are patched: the paired
    and first/second-in-pair flag bits, the mate reference and the mate
    position.

    The fifth patch, to the insert size, is applied since HiCExplorer 4 (Q5):
    buildMatrixMethods.py:814 patched it in the worker process and wrote the
    master's own copy, so every record of 3.7.6 has isize 0. Now a pair on one
    chromosome carries mate2.pos - mate1.pos and the opposite sign.
    """
    out_bam = str(tmp_path / "valid.bam")
    _, qc_folder = _small_fixture(str(tmp_path), out_bam=out_bam)

    with pysam.Samfile(out_bam, "rb") as handle:
        header = str(handle.header)
        records = list(handle)
    with pysam.Samfile(R1_1000, "rb") as source:
        assert header == str(source.header)

    assert len(records) == 2 * 195  # two per Hi-C contact
    for first, second in zip(records[0::2], records[1::2]):
        assert first.qname == second.qname
        assert first.flag & 0x1 and second.flag & 0x1
        assert first.flag & 0x40 and not first.flag & 0x80
        assert second.flag & 0x80 and not second.flag & 0x40
        assert first.mrnm == second.rname
        assert second.mrnm == first.rname
        assert first.mpos == second.pos
        assert second.mpos == first.pos
        if first.reference_id == second.reference_id:
            assert first.isize == second.pos - first.pos
            assert second.isize == first.pos - second.pos
    assert any(record.isize != 0 for record in records)

    # The record order is the order the pairs appear in the input files, which
    # is what a single worker produces. It is NOT stable against --threads:
    # the master consumes the worker queues in completion order, so at
    # --threads 8 the same 74,642 records of the larger library come out in a
    # different order. That is a defect of the tool, recorded here rather than
    # asserted, because asserting a race is not a test.
    with pysam.Samfile(R1_1000, "rb") as source:
        input_order = [read.qname for read in source]
    written_order = [record.qname for record in records[0::2]]
    assert written_order == [name for name in input_order if name in set(written_order)]
    shutil.rmtree(qc_folder)


def test_build_matrix_region_keeps_the_restriction_sites_inside_it(tmp_path):
    """HiCExplorer 4 (Q6): a bare --region keeps the cut sites on the chromosome.

    bed2interval_list kept a cut site only when
    `chrom == region_chrom and region_start <= site_start and region_end <=
    site_end`, and region_end defaults to the length of the chromosome, so a
    bare --region chrX demanded that the site reach the chromosome end and the
    list came out empty in 3.7.6. rf_positions was then falsy, no self circle
    was ever counted, and every close inward pair fell through to "same
    fragment". Now a site is kept when it lies inside the region.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False, dir=str(tmp_path))
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_", dir=str(tmp_path))
    args = ("-s {} {} --outFileName {} -bs 100000 --QCfolder {} --threads 4 "
            "--restrictionSequence AAGCTT --danglingSequence AGCT -rs {} "
            "--region chr3R").format(R1_1000, R2_1000, outfile.name, qc_folder,
                                     RS_HINDIII).split()
    compute(hicBuildMatrix.main, args, 5)

    body = _qc_body(qc_folder)
    assert "self circle\t3\t(1.00)\n" in body
    assert "same fragment\t20\t(6.67)\n" in body
    # The matrix covers chr3R alone.
    matrix = hm.hiCMatrix(outfile.name)
    assert list(matrix.chrBinBoundaries) == ['chr3R']
    shutil.rmtree(qc_folder)


def _mcool_run(tmp_path, bin_sizes, suffix='.mcool', extra=''):
    outfile = NamedTemporaryFile(suffix=suffix, delete=False, dir=str(tmp_path))
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_", dir=str(tmp_path))
    args = ("-s {} {} --outFileName {} -bs {} --QCfolder {} "
            "--threads 4 --restrictionSequence AAGCTT --danglingSequence AGCT "
            "-rs {} {}").format(R1_1000, R2_1000, outfile.name, bin_sizes, qc_folder,
                                RS_HINDIII, extra)
    return outfile.name, qc_folder, args.split()


def test_build_matrix_two_resolutions_write_an_mcool(tmp_path):
    """HiCExplorer 4 (fix 3): an .mcool name writes one resolution per bin size.

    The mcool branch of buildMatrixMethods.py:1379 is guarded by
    `len(pBinSize) > 2` in 3.7.6, so two bin sizes fell through to
    hiCMatrix.save, which dispatches on `name.endswith('cool')` and wrote a
    single cooler at the file root at the first bin size, ignoring the second.
    Every resolution carries the provenance attributes as plain strings (3.7.6
    stored the repr "b'HiCExplorer-3.7.6'" and only on the first resolution).
    """
    import h5py

    outfile, qc_folder, args = _mcool_run(tmp_path, '100000 200000')
    compute(hicBuildMatrix.main, args, 5)

    with h5py.File(outfile, 'r') as handle:
        assert list(handle.keys()) == ['resolutions']
        assert sorted(handle['resolutions'].keys()) == ['100000', '200000']
        assert handle['resolutions/100000'].attrs['bin-size'] == 100000
        assert handle['resolutions/200000'].attrs['bin-size'] == 200000
        assert int(handle['resolutions/100000'].attrs['nbins']) == 1697
        assert int(handle['resolutions/200000'].attrs['nbins']) < 1697
        for resolution in ['100000', '200000']:
            metadata = json.loads(handle['resolutions/' + resolution].attrs['metadata'])
            assert metadata['matrix-generated-by'] == 'HiCExplorer-4.0.0'
            assert metadata['matrix-generated-by-url'] == 'https://github.com/deeptools/HiCExplorer'
    first = hm.hiCMatrix(outfile + '::/resolutions/100000')
    second = hm.hiCMatrix(outfile + '::/resolutions/200000')
    # the coarser resolution merges pairs of bins of the first one
    assert first.matrix.shape[0] // 2 <= second.matrix.shape[0] <= first.matrix.shape[0] // 2 + len(first.chrBinBoundaries)
    assert second.matrix.sum() > 0
    shutil.rmtree(qc_folder)


def test_build_matrix_one_bin_size_in_an_mcool(tmp_path):
    """HiCExplorer 4 (fix 3): one bin size is enough for an .mcool."""
    import h5py

    outfile, qc_folder, args = _mcool_run(tmp_path, '100000')
    compute(hicBuildMatrix.main, args, 5)
    with h5py.File(outfile, 'r') as handle:
        assert list(handle.keys()) == ['resolutions']
        assert list(handle['resolutions'].keys()) == ['100000']
        metadata = json.loads(handle['resolutions/100000'].attrs['metadata'])
        assert metadata['matrix-generated-by'] == 'HiCExplorer-4.0.0'
    shutil.rmtree(qc_folder)


def test_build_matrix_provenance_attributes_are_plain_strings(tmp_path):
    """HiCExplorer 4 (fix 3): 3.7.6 built the provenance with np.string_ and
    the cool writer stored its repr, "b'HiCExplorer-3.7.6'"."""
    import h5py

    outfile, qc_folder = _small_fixture(str(tmp_path), suffix='.cool')
    with h5py.File(outfile, 'r') as handle:
        assert handle.attrs['matrix-generated-by'] == 'HiCExplorer-4.0.0'
        assert handle.attrs['matrix-generated-by-url'] == 'https://github.com/deeptools/HiCExplorer'
        for value in handle.attrs.values():
            assert not str(value).startswith("b'")
    shutil.rmtree(qc_folder)


def test_build_matrix_several_bin_sizes_need_an_mcool(tmp_path):
    """HiCExplorer 4 (fix 3): two bin sizes with an .h5 or .cool name used to
    write only the first one without a word; now the run stops with status 1
    before it reads the input."""
    for suffix in ['.h5', '.cool']:
        outfile, qc_folder, args = _mcool_run(tmp_path, '100000 200000', suffix=suffix)
        with pytest.raises(SystemExit) as error:
            hicBuildMatrix.main(args)
        assert error.value.code == 1
        shutil.rmtree(qc_folder)


def test_build_matrix_further_bin_sizes_must_be_larger_multiples_of_the_first(tmp_path):
    """HiCExplorer 4 (fix 3): the coarser resolutions are merged from the first
    one, so a bin size that is not a larger multiple of it exits with status 1."""
    for bin_sizes in ['100000 150000', '200000 100000', '100000 100000']:
        outfile, qc_folder, args = _mcool_run(tmp_path, bin_sizes)
        with pytest.raises(SystemExit) as error:
            hicBuildMatrix.main(args)
        assert error.value.code == 1
        shutil.rmtree(qc_folder)


def test_build_matrix_do_test_run_with_an_mcool_name(tmp_path):
    """HiCExplorer 4 (fix 3): --doTestRun with an .mcool name and several bin
    sizes raised the UnboundLocalError of the Python API path in the reproduced
    createMatrix; now it behaves like the other output names."""
    outfile, qc_folder, args = _mcool_run(tmp_path, '100000 200000 300000', extra='--doTestRun --doTestRunLines 500')
    compute(hicBuildMatrix.main, args, 5)
    assert _qc_body(qc_folder)[1:] == EXPECTED_QC_LOG
    shutil.rmtree(qc_folder)
