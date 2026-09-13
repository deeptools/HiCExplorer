from hicexplorer import chicQualityControl
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure

import matplotlib as mpl
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import warnings

import numpy as np
from hicmatrix import HiCMatrix as hm

from hicexplorer.test.test_compute_function import compute

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")


def are_files_equal(file1, file2, delta=1, skip=0):
    equal = True
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            # if x.startswith('File'):
            #     continue
            if i < skip:
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


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_two_matrices():

    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile_histogram = NamedTemporaryFile(suffix='.png', delete=False)
    outfile_sparsity = NamedTemporaryFile(suffix='.png', delete=False)

    outfile.close()
    args = "--matrices {} {} --referencePoints {} --sparsity {} --outFileName {} --outFileNameHistogram {} --outFileNameSparsity  {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool',
                                                                                                                                                   ROOT + 'MB-E10-5_chr1.cool',
                                                                                                                                                   ROOT + 'referencePoints_qc.bed',
                                                                                                                                                   0.05,
                                                                                                                                                   outfile.name, outfile_histogram.name, outfile_sparsity.name, 1).split()
    # chicQualityControl.main(args)
    compute(chicQualityControl.main, args, 5)
    assert are_files_equal(
        ROOT + "chicQualityControl/new_referencepoints.bed", outfile.name)
    assert are_files_equal(
        ROOT + "chicQualityControl/new_referencepoints.bed_raw_filter", outfile.name + '_raw_filter', skip=4)
    assert are_files_equal(ROOT + "chicQualityControl/new_referencepoints.bed_rejected_filter",
                           outfile.name + '_rejected_filter', skip=2)
    assert are_files_equal(ROOT + "chicQualityControl/new_referencepoints.bed_report",
                           outfile.name + '_report', skip=2)
    assert are_files_equal(ROOT + "chicQualityControl/new_referencepoints.bed_failed_reference_points",
                           outfile.name + '_report', skip=2)

    res = compare_images(
        ROOT + "chicQualityControl/histogram.png", outfile_histogram.name, 50)
    assert res is None, res

    res = compare_images(
        ROOT + "chicQualityControl/sparsity.png", outfile_sparsity.name, 50)
    assert res is None, res


# ---------------------------------------------------------------------------
# Characterization tests added for the C++ port (cpp/AGENTS_CONTRACT.md rule 1).
#
# What test_two_matrices above actually constrains: are_files_equal walks the
# two files with zip(), so an output shorter than the expected file, including
# an empty one, passes; it tolerates one differing line per file; and the last
# assertion compares the failed reference points against the *report* with
# skip=2, which checks nothing, since the failed file has exactly two lines. The
# test is also marked xfail, so it is reported as XPASS rather than PASSED. The
# sparsity values themselves are never compared with anything computed
# independently of lib/viewpoint.py.

QC_OUTPUT_SUFFIXES = ('', '_raw_filter', '_rejected_filter', '_report',
                      '_failed_reference_points')

# The committed expected files were written by version 3.4.3-dev from relative
# matrix paths. The version string and the matrix directory are provenance, and
# are the only two things normalised before the comparison.
VERSION_LINE_PREFIX = '# Created with chicQualityControl version '
MATRIX_LIST_PREFIXES = ('# Used Matrices ', '# QC report for matrices: ')

# Both cHi-C test matrices hold chr1 only, in fixed 1 kb bins, and the last bin
# ends at 197,195,432.
CHR1_LENGTH = 197195432
BIN_SIZE = 1000


def normalised_lines(pPath, pMatrixDirectory=None):
    with open(pPath) as handle:
        lines = handle.read().split('\n')
    result = []
    for line in lines:
        if line.startswith(VERSION_LINE_PREFIX):
            line = VERSION_LINE_PREFIX + '<version>'
        elif pMatrixDirectory and line.startswith(MATRIX_LIST_PREFIXES):
            line = line.replace(pMatrixDirectory, '')
        result.append(line)
    return result


def assert_same_lines(pExpected, pObserved, pMatrixDirectory=None):
    """Line by line equality, including the line count and the final newline."""
    expected = normalised_lines(pExpected)
    observed = normalised_lines(pObserved, pMatrixDirectory)
    assert len(expected) == len(observed), \
        '{}: {} lines expected, {} observed'.format(pObserved, len(expected), len(observed))
    for index, (left, right) in enumerate(zip(expected, observed)):
        assert left == right, '{} line {}:\n  expected {!r}\n  observed {!r}'.format(
            pObserved, index + 1, left, right)


def run_quality_control(pReferencePoints, pSparsity, pThreads=1, pExtraArguments=''):
    output_folder = mkdtemp(prefix='chicQualityControl_')
    outfile = os.path.join(output_folder, 'referencepoints.bed')
    args = "--matrices {} {} --referencePoints {} --sparsity {} --outFileName {} " \
        "--outFileNameHistogram {} --outFileNameSparsity {} -t {} {}".format(
            ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool', pReferencePoints,
            pSparsity, outfile, os.path.join(output_folder, 'histogram.png'),
            os.path.join(output_folder, 'sparsity.png'), pThreads,
            pExtraArguments).split()
    chicQualityControl.main(args)
    return outfile


def read_lines(pPath):
    with open(pPath) as handle:
        return handle.readlines()


def read_raw_filter(pPath):
    return [line.rstrip('\n').split('\t') for line in read_lines(pPath)
            if not line.startswith('#')]


def independent_sparsity(pMatrix, pChromosome, pStart, pEnd, pFixateRange):
    """The sparsity of one viewpoint, recomputed without lib/viewpoint.py.

    The window is fixateRange on either side of the reference point, clipped to
    0 on the left and to one base before the chromosome end on the right, as
    calculateViewpointRange clips it. A reference point spanning several bins
    is summed into a single element. A position outside the chromosome, or a
    chromosome the matrix does not have, is a faulty reference point, -1.0.
    """
    if pChromosome != 'chr1' or not (0 <= pStart < CHR1_LENGTH) \
            or not (0 <= pEnd < CHR1_LENGTH):
        return -1.0
    first_bin = pStart // BIN_SIZE
    last_bin = pEnd // BIN_SIZE
    region_start = max(pStart - pFixateRange, 0)
    region_end = pEnd + pFixateRange
    if region_end > CHR1_LENGTH:
        region_end = CHR1_LENGTH - 1
    low = region_start // BIN_SIZE
    high = region_end // BIN_SIZE
    rows = np.asarray(pMatrix.matrix[first_bin:last_bin + 1, low:high + 1].todense()).sum(axis=0)
    collapsed = np.concatenate([rows[:first_bin - low],
                                [rows[first_bin - low:last_bin - low + 1].sum()],
                                rows[last_bin - low + 1:]])
    return np.count_nonzero(collapsed) / len(collapsed)


def test_text_outputs_are_exact():
    outfile = run_quality_control(ROOT + 'referencePoints_qc.bed', 0.05)
    for suffix in QC_OUTPUT_SUFFIXES:
        assert_same_lines(ROOT + 'chicQualityControl/new_referencepoints.bed' + suffix,
                          outfile + suffix, pMatrixDirectory=ROOT)


def test_sparsity_matches_an_independent_computation():
    """--fixateRange is also untested elsewhere, so a non default value is used."""
    outfile = run_quality_control(ROOT + 'referencePoints_qc.bed', 0.05,
                                  pExtraArguments='--fixateRange 100000')
    matrices = [hm.hiCMatrix(ROOT + 'FL-E13-5_chr1.cool'),
                hm.hiCMatrix(ROOT + 'MB-E10-5_chr1.cool')]
    rows = read_raw_filter(outfile + '_raw_filter')
    assert len(rows) == 45
    faulty = 0
    positive = 0
    for row in rows:
        assert len(row) == 6, row
        chromosome, start, end = row[0], int(row[1]), int(row[2])
        for matrix, observed in zip(matrices, row[4:6]):
            expected = independent_sparsity(matrix, chromosome, start, end, 100000)
            assert float(observed) == expected, '{}: {} observed, {} expected'.format(
                row[:4], observed, expected)
            faulty += expected == -1.0
            positive += expected > 0.0
    assert faulty == 4
    assert positive > 70


def test_a_viewpoint_equal_to_the_threshold_is_rejected():
    """The comparison is strict: sparsity > --sparsity accepts.

    Sox17 has 0.14485514485514486 in FL-E13-5 and 0.16083916083916083 in
    MB-E10-5, so a threshold equal to the larger value rejects it.
    """
    outfile = run_quality_control(ROOT + 'referencePoints_qc.bed', '0.16083916083916083')
    line = 'chr1\t4487435\t4487435\tSox17\n'
    assert line not in read_lines(outfile)
    assert line in read_lines(outfile + '_rejected_filter')


def test_one_passing_matrix_is_enough_to_accept():
    """Pinned as the Python behaves, which is not what its help text says.

    --sparsity is documented as removing a viewpoint 'as soon as it is of bad
    quality in at least one matrix', but chicQualityControl.py:228-240 accepts
    a reference point when any matrix is above the threshold. Tfap2b has
    0.01898101898101898 in FL-E13-5 and 0.016983016983016984 in MB-E10-5, so
    at 0.018 it is below the threshold in one matrix and still accepted.
    """
    outfile = run_quality_control(ROOT + 'referencePoints_qc.bed', 0.018)
    line = 'chr1\t19198995\t19198995\tTfap2b\n'
    assert line in read_lines(outfile)
    assert line not in read_lines(outfile + '_rejected_filter')


def test_thread_count_does_not_change_the_output():
    single = run_quality_control(ROOT + 'referencePoints_qc.bed', 0.05, pThreads=1)
    several = run_quality_control(ROOT + 'referencePoints_qc.bed', 0.05, pThreads=4)
    for suffix in QC_OUTPUT_SUFFIXES:
        assert_same_lines(single + suffix, several + suffix)
