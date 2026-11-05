import os.path
from tempfile import NamedTemporaryFile, mkdtemp
import shutil

import logging
log = logging.getLogger(__name__)

from hicexplorer import hicQuickQC
from hicexplorer.test.test_compute_function import compute

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")

sam_R1 = ROOT + "small_test_R1_unsorted.bam"
sam_R2 = ROOT + "small_test_R2_unsorted.bam"


def are_files_equal(file1, file2, delta=None):
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


def test_main():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    qc_folder = mkdtemp(prefix="testQC_")
    args = "-s {} {} --QCfolder {} -seq {} --danglingSequence {} --restrictionCutFile {} --lines 1000".format(sam_R1, sam_R2,
                                                                                                              qc_folder, 'GATC', 'GATC', ROOT + "DpnII.bed").split()
    # hicQuickQC.main(args)
    compute(hicQuickQC.main, args, 5)

    # print(set(os.listdir(ROOT + "hicQuickQC/")))
    assert are_files_equal(ROOT + "hicQuickQC/QC.log", qc_folder + "/QC.log")
    assert set(os.listdir(ROOT + "hicQuickQC/")) == set(os.listdir(qc_folder))

    shutil.rmtree(qc_folder)


# ---------------------------------------------------------------------------
# Full comparison, added 2026-09-13 for the C++ port.
#
# are_files_equal above walks both files with zip(), which stops at the end of
# the shorter one, and it skips every line that starts with 'File'. A QC.log
# with extra lines at its end therefore still passes: measured with a mutation
# of buildMatrixMethods.py that writes one more empty line after the last
# counter, which test_main accepts and the line count below rejects. test_main
# also never opens the five tables hicPrepareQCreport derives from QC.log,
# which are what hicQC reads back: a mutation of the discarded table's
# denominator passes test_main and fails
# test_qc_tables_carry_the_counts_of_the_reference.
#
# Two properties of the reference are pinned here:
#
#  1. The second line of QC.log, and the first column of every table, name the
#     NamedTemporaryFile that hicQuickQC.py:100 passes to hicBuildMatrix as
#     --outFileName. The name is random, and buildMatrixMethods.py:1366
#     unlinks the file after the report, so the name points at nothing.
#  2. The checked-in QC_table.txt and discarded_table.txt were written by an
#     older release, whose dangling end label was 'dangling end GATC'. The
#     current QC.log, and so the current tables, say
#     'dangling end GATC (restriction sequence GATC)'. The counts are
#     unchanged, and the other three tables and QC.log match the reference.

import re

TEMPORARY_MATRIX_NAME = re.compile(r'^(/[^\t\n]*/)?tmp[a-z0-9_]{8}\.h5$')
QC_TABLES = ["QC_table.txt", "discarded_table.txt", "distance_table.txt",
             "read_orientation_table.txt", "unmapable_table.txt"]


def _run_quick_qc(pQCfolder):
    args = "-s {} {} --QCfolder {} -seq {} --danglingSequence {} " \
           "--restrictionCutFile {} --lines 1000".format(
               sam_R1, sam_R2, pQCfolder, 'GATC', 'GATC', ROOT + "DpnII.bed").split()
    hicQuickQC.main(args)


def _read_lines(pPath):
    with open(pPath) as handle:
        return handle.read().split("\n")


def test_qc_log_matches_the_reference_line_for_line():
    qc_folder = mkdtemp(prefix="testQC_")
    _run_quick_qc(qc_folder)

    expected = _read_lines(ROOT + "hicQuickQC/QC.log")
    actual = _read_lines(qc_folder + "/QC.log")
    assert len(actual) == len(expected), \
        'zip() in are_files_equal would hide a missing or an extra line'

    file_line = actual[1].split("\t")
    assert file_line[0] == "File"
    assert file_line[2:] == ["", ""]
    assert TEMPORARY_MATRIX_NAME.match(file_line[1]), file_line[1]
    assert not os.path.exists(file_line[1]), \
        'the temporary matrix is unlinked after the QC report'

    for number, (line_actual, line_expected) in enumerate(zip(actual, expected), 1):
        if number == 2:
            continue
        assert line_actual == line_expected, 'QC.log line {}'.format(number)
    shutil.rmtree(qc_folder)


def test_qc_tables_carry_the_counts_of_the_reference():
    qc_folder = mkdtemp(prefix="testQC_")
    _run_quick_qc(qc_folder)
    temporary_name = _read_lines(qc_folder + "/QC.log")[1].split("\t")[1]

    for table in QC_TABLES:
        expected = _read_lines(ROOT + "hicQuickQC/" + table)
        actual = _read_lines(qc_folder + "/" + table)
        assert len(actual) == len(expected) == 3, table  # header, row, ''
        assert actual[2] == ""

        header_expected = expected[0]
        if table in ("QC_table.txt", "discarded_table.txt"):
            # property 2: the reference predates the longer label
            header_expected = header_expected.replace(
                "dangling end GATC", "dangling end GATC (restriction sequence GATC)")
        assert actual[0] == header_expected, table

        row_actual = actual[1].split("\t")
        row_expected = expected[1].split("\t")
        assert row_actual[0] == temporary_name, table
        assert TEMPORARY_MATRIX_NAME.match(row_expected[0])
        assert row_actual[1:] == row_expected[1:], table
    shutil.rmtree(qc_folder)
