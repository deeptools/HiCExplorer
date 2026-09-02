import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import hashlib
import json
import os
import shutil
from tempfile import mkdtemp

import numpy.testing as nt
import pytest
from hicmatrix import HiCMatrix as hm

from hicexplorer import hicFindTADs
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

# The five text files every run of the tool produces, without the prefix.
OUTPUT_SUFFIXES = ['_boundaries.bed', '_boundaries.gff', '_domains.bed',
                   '_score.bedgraph', '_tad_score.bm']


def are_files_equal(file1, file2, pDifference=10):
    """The original comparison, kept because the four original tests use it.

    It is far weaker than it looks and the tests below exist because of that:
    ``zip`` stops at the shorter of the two files, so a file with a different
    number of lines is not detected at all, and up to ``pDifference``
    characters may differ on every line that is compared.
    """
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


def assert_files_identical(expected, actual):
    """Byte-for-byte equality, with the first difference reported.

    This is what ``are_files_equal`` was meant to be. It compares the number of
    lines as well as their content, so a truncated or over-long output fails.
    """
    with open(expected) as fh:
        expected_lines = fh.readlines()
    with open(actual) as fh:
        actual_lines = fh.readlines()

    for index, (want, got) in enumerate(zip(expected_lines, actual_lines)):
        assert want == got, (
            "{} line {} differs\n  expected: {!r}\n  actual:   {!r}".format(
                actual, index + 1, want, got))
    assert len(expected_lines) == len(actual_lines), (
        "{} has {} lines, expected {}".format(
            actual, len(actual_lines), len(expected_lines)))


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        digest.update(fh.read())
    return digest.hexdigest()


def line_count(path):
    with open(path) as fh:
        return sum(1 for _ in fh)


def assert_file_pinned(path, digest, lines):
    """Full-precision assertion on a text output without a reference file.

    The line count is asserted separately from the digest only so that the
    failure message says which of the two went wrong.
    """
    assert line_count(path) == lines, (
        "{} has {} lines, expected {}".format(path, line_count(path), lines))
    assert sha256_of(path) == digest, (
        "{} content changed; first lines:\n{}".format(
            path, ''.join(open(path).readlines()[:5])))


def run_tool(argument_string):
    compute(hicFindTADs.main, argument_string.split(), 5)


def prepare_precomputed(folder, prefix, reference):
    """Copy a TAD-separation score and its z-score matrix next to ``prefix``.

    This is what the bonferroni and None cases do: it makes the run skip the
    expensive spectrum computation and exercise only the boundary calling.
    """
    shutil.copy(ROOT + "find_TADs/{0}/multi{1}_tad_score.bm".format(reference[0], reference[1]),
                os.path.join(folder, prefix + "_tad_score.bm"))
    shutil.copy(ROOT + "find_TADs/{0}/multi{1}_zscore_matrix.h5".format(reference[0], reference[1]),
                os.path.join(folder, prefix + "_zscore_matrix.h5"))


# ---------------------------------------------------------------------------
# The four original tests, unchanged.
# ---------------------------------------------------------------------------


def test_find_TADs_fdr():
    # full test case with build of the matrix and search for tads
    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_fdr")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 --step 20000 \
    --outPrefix {}/test_multiFDR --minBoundaryDistance 20000 \
    --correctForMultipleTesting fdr --thresholdComparisons 0.1".format(matrix, tad_folder).split()

    # hicFindTADs.main(args)
    compute(hicFindTADs.main, args, 5)
    new = hm.hiCMatrix(tad_folder + "/test_multiFDR_zscore_matrix.h5")
    test = hm.hiCMatrix(ROOT + 'find_TADs/FDR/multiFDR_zscore_matrix.h5')
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(tad_folder + "/test_multiFDR_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR/multiFDR_boundaries.bed", tad_folder + "/test_multiFDR_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR/multiFDR_domains.bed", tad_folder + "/test_multiFDR_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR/multiFDR_tad_score.bm", tad_folder + "/test_multiFDR_tad_score.bm")
    assert are_files_equal(ROOT + "find_TADs/FDR/multiFDR_boundaries.gff", tad_folder + "/test_multiFDR_boundaries.gff")
    # assert are_files_equal
    assert are_files_equal(ROOT + "find_TADs/FDR/multiFDR_score.bedgraph", tad_folder + "/test_multiFDR_score.bedgraph")

    shutil.rmtree(tad_folder)


def test_find_TADs_fdr_chromosomes():
    # full test case with build of the matrix and search for tads
    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_fdr_chromosomes")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 --step 20000 \
    --outPrefix {}/test_multiFDR_chromosomes --minBoundaryDistance 20000 \
    --correctForMultipleTesting fdr --thresholdComparisons 0.5 --chromosomes chr2L chr3R".format(matrix, tad_folder).split()

    # hicFindTADs.main(args)
    compute(hicFindTADs.main, args, 5)

    new = hm.hiCMatrix(tad_folder + "/test_multiFDR_chromosomes_zscore_matrix.h5")
    test = hm.hiCMatrix(ROOT + 'find_TADs/FDR_chromosomes/multiFDR_zscore_matrix.h5')
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    print(tad_folder + "/test_multiFDR_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR_chromosomes/multiFDR_boundaries.bed", tad_folder + "/test_multiFDR_chromosomes_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR_chromosomes/multiFDR_domains.bed", tad_folder + "/test_multiFDR_chromosomes_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/FDR_chromosomes/multiFDR_tad_score.bm", tad_folder + "/test_multiFDR_chromosomes_tad_score.bm")
    assert are_files_equal(ROOT + "find_TADs/FDR_chromosomes/multiFDR_boundaries.gff", tad_folder + "/test_multiFDR_chromosomes_boundaries.gff")
    # assert are_files_equal
    assert are_files_equal(ROOT + "find_TADs/FDR_chromosomes/multiFDR_score.bedgraph", tad_folder + "/test_multiFDR_chromosomes_score.bedgraph")

    shutil.rmtree(tad_folder)


def test_find_TADs_bonferroni():
    # reduced test case, the z-score matrix is given to decrease run time
    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_bonferroni")
    shutil.copy(ROOT + "find_TADs/bonferroni/multiBonferroni_tad_score.bm", tad_folder + "/test_multiBonferroni_tad_score.bm")
    shutil.copy(ROOT + 'find_TADs/bonferroni/multiBonferroni_zscore_matrix.h5', tad_folder + "/test_multiBonferroni_zscore_matrix.h5")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 --step 20000 \
    --outPrefix {}/test_multiBonferroni --minBoundaryDistance 20000 \
    --correctForMultipleTesting bonferroni --thresholdComparisons 0.1".format(matrix, tad_folder).split()

    # hicFindTADs.main(args)
    compute(hicFindTADs.main, args, 5)

    print(tad_folder + "/test_multiBonferroni_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/bonferroni/multiBonferroni_boundaries.bed", tad_folder + "/test_multiBonferroni_boundaries.bed", pDifference=10)
    assert are_files_equal(ROOT + "find_TADs/bonferroni/multiBonferroni_domains.bed", tad_folder + "/test_multiBonferroni_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/bonferroni/multiBonferroni_boundaries.gff", tad_folder + "/test_multiBonferroni_boundaries.gff")
    assert are_files_equal(ROOT + "find_TADs/bonferroni/multiBonferroni_score.bedgraph", tad_folder + "/test_multiBonferroni_score.bedgraph")

    shutil.rmtree(tad_folder)


def test_find_TADs_none():
    # reduced test case, the z-score matrix is given to decrease run time
    matrix = ROOT + "small_test_matrix.h5"
    tad_folder = mkdtemp(prefix="test_case_find_tads_none")
    shutil.copy(ROOT + "find_TADs/None/multiNone_tad_score.bm", tad_folder + "/test_multiNone_tad_score.bm")
    shutil.copy(ROOT + 'find_TADs/None/multiNone_zscore_matrix.h5', tad_folder + "/test_multiNone_zscore_matrix.h5")
    args = "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 --step 20000 \
    --outPrefix {}/test_multiNone --minBoundaryDistance 20000 \
    --correctForMultipleTesting None --thresholdComparisons 1.0".format(matrix, tad_folder).split()

    # hicFindTADs.main(args)
    compute(hicFindTADs.main, args, 5)

    assert are_files_equal(ROOT + "find_TADs/None/multiNone_boundaries.bed", tad_folder + "/test_multiNone_boundaries.bed")
    assert are_files_equal(ROOT + "find_TADs/None/multiNone_domains.bed", tad_folder + "/test_multiNone_domains.bed")
    assert are_files_equal(ROOT + "find_TADs/None/multiNone_boundaries.gff", tad_folder + "/test_multiNone_boundaries.gff")
    assert are_files_equal(ROOT + "find_TADs/None/multiNone_score.bedgraph", tad_folder + "/test_multiNone_score.bedgraph")

    shutil.rmtree(tad_folder)


# ---------------------------------------------------------------------------
# Characterization tests added for the C++ port (cpp/AGENTS_CONTRACT.md rule 1).
#
# Everything below asserts at full precision. The reference files committed
# under test_data/find_TADs/ are byte-identical to what the Python tool
# produces today, verified 2026-09-01, so the exact comparison is a valid
# characterization of the current behaviour and not a tightening of a target
# the reference never met.
# ---------------------------------------------------------------------------


def test_find_TADs_fdr_is_byte_identical_to_the_reference():
    """The FDR case at full precision, in place of a 10-character budget.

    This also pins the two outputs the original FDR test could not constrain:
    ``domains.bed`` is empty in the reference, so any comparison against it
    passes, and ``are_files_equal`` would accept a truncated ``score.bedgraph``.
    """
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_fdr_exact")
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 "
             "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
             "--correctForMultipleTesting fdr --thresholdComparisons 0.1".format(matrix, folder))

    for suffix in OUTPUT_SUFFIXES:
        assert_files_identical(ROOT + "find_TADs/FDR/multiFDR" + suffix, folder + "/t" + suffix)

    # The z-score matrix has to agree in the values, the sparsity pattern and
    # the bin table, not only in the values.
    new = hm.hiCMatrix(folder + "/t_zscore_matrix.h5")
    reference = hm.hiCMatrix(ROOT + 'find_TADs/FDR/multiFDR_zscore_matrix.h5')
    nt.assert_equal(reference.matrix.indptr, new.matrix.indptr)
    nt.assert_equal(reference.matrix.indices, new.matrix.indices)
    nt.assert_equal(reference.matrix.data, new.matrix.data)
    nt.assert_equal(reference.cut_intervals, new.cut_intervals)
    nt.assert_equal(reference.nan_bins, new.nan_bins)

    shutil.rmtree(folder)


def test_find_TADs_none_is_byte_identical_to_the_reference():
    """The None case at full precision.

    This is the richest call set in the corpus, 1,407 boundaries and 1,402
    domains, and therefore the strongest regression detector for the boundary
    calling step. The original test accepted 10 differing characters per line
    and stopped comparing at the end of the shorter file.
    """
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_none_exact")
    prepare_precomputed(folder, "t", ("None", "None"))
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 "
             "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
             "--correctForMultipleTesting None --thresholdComparisons 1.0".format(matrix, folder))

    for suffix in ['_boundaries.bed', '_boundaries.gff', '_domains.bed', '_score.bedgraph']:
        assert_files_identical(ROOT + "find_TADs/None/multiNone" + suffix, folder + "/t" + suffix)
    assert line_count(folder + "/t_boundaries.bed") == 1407
    assert line_count(folder + "/t_domains.bed") == 1402

    shutil.rmtree(folder)


def test_find_TADs_bonferroni_is_byte_identical_to_the_reference():
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_bonferroni_exact")
    prepare_precomputed(folder, "t", ("bonferroni", "Bonferroni"))
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 "
             "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
             "--correctForMultipleTesting bonferroni --thresholdComparisons 0.1".format(matrix, folder))

    for suffix in ['_boundaries.bed', '_boundaries.gff', '_domains.bed', '_score.bedgraph']:
        assert_files_identical(ROOT + "find_TADs/bonferroni/multiBonferroni" + suffix,
                               folder + "/t" + suffix)

    shutil.rmtree(folder)


def test_find_TADs_chromosomes_is_byte_identical_to_the_reference():
    """--chromosomes restricts and reorders the matrix before anything else."""
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_chrom_exact")
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 2 "
             "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
             "--correctForMultipleTesting fdr --thresholdComparisons 0.5 "
             "--chromosomes chr2L chr3R".format(matrix, folder))

    for suffix in OUTPUT_SUFFIXES:
        assert_files_identical(ROOT + "find_TADs/FDR_chromosomes/multiFDR" + suffix,
                               folder + "/t" + suffix)

    # Only the two requested chromosomes survive, in the order they were given.
    with open(folder + "/t_score.bedgraph") as fh:
        chromosomes = []
        for line in fh:
            name = line.split('\t')[0]
            if name not in chromosomes:
                chromosomes.append(name)
    assert chromosomes == ['chr2L', 'chr3R']

    shutil.rmtree(folder)


def test_find_TADs_result_does_not_depend_on_numberOfProcessors():
    """--numberOfProcessors must not change a single value.

    The Python splits the bins into one contiguous range per process and
    concatenates the results in range order, so the result is independent of
    the process count. The C++ port replaces the process pool with threads and
    has to keep exactly this property (cpp/OPTIMIZATION.md 3).
    """
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_processors")
    for processors in (1, 4):
        run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 "
                 "--numberOfProcessors {} --step 20000 --outPrefix {}/p{} "
                 "--minBoundaryDistance 20000 --correctForMultipleTesting fdr "
                 "--thresholdComparisons 0.1".format(matrix, processors, folder, processors))

    for suffix in OUTPUT_SUFFIXES:
        assert_files_identical(folder + "/p1" + suffix, folder + "/p4" + suffix)
        assert_files_identical(ROOT + "find_TADs/FDR/multiFDR" + suffix, folder + "/p1" + suffix)

    shutil.rmtree(folder)


def test_find_TADs_high_numberOfProcessors_is_a_defect():
    """Pinned defect: a large --numberOfProcessors aborts the run.

    ``np.array_split`` gives every worker a contiguous range of bins.
    ``compute_matrix`` drops a bin whose TAD-separation score is None or NaN at
    any window size, and when every bin of one worker's range is dropped the
    unpack at hicFindTADs.py:345 raises

        ValueError: not enough values to unpack (expected 3, got 0)

    On small_test_matrix.h5 that happens from --numberOfProcessors 12 upwards,
    while 1 to 10 all succeed and produce byte-identical output. So the process
    count decides whether the tool runs at all, which is a defect and not a
    numeric difference. The C++ port must not reproduce it: an empty partition
    contributes no rows and the run continues. Reported, not fixed here.
    """
    matrix = ROOT + "small_test_matrix.h5"
    folder = mkdtemp(prefix="test_find_tads_processors_defect")
    with pytest.raises(ValueError, match="not enough values to unpack"):
        hicFindTADs.main(
            "--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 12 "
            "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
            "--correctForMultipleTesting fdr --thresholdComparisons 0.1".format(
                matrix, folder).split())

    shutil.rmtree(folder)


def test_find_TADs_on_a_cool_matrix():
    """cool input, which no existing test covers.

    The same region as small_test_matrix.h5 but without the masked bins the h5
    carries, so it produces a much larger call set: 555 boundaries against the
    2 of the FDR h5 case.
    """
    matrix = ROOT + "small_test_matrix.cool"
    folder = mkdtemp(prefix="test_find_tads_cool")
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --numberOfProcessors 1 "
             "--step 20000 --outPrefix {}/t --minBoundaryDistance 20000 "
             "--correctForMultipleTesting fdr --thresholdComparisons 0.1".format(matrix, folder))

    # The z-score matrix keeps the suffix of the input, so this is a .cool.
    assert os.path.isfile(folder + "/t_zscore_matrix.cool")
    assert_file_pinned(folder + "/t_boundaries.bed",
                       "cdc50fb8e6ed6458c49ff154a6da5ee26e24cf3e54fd45adab5fc47d315329da", 555)
    assert_file_pinned(folder + "/t_boundaries.gff",
                       "de4542b1a382d5163b4b911e71e83ab4c25286c6dc159df929c37c50c6c320dc", 555)
    assert_file_pinned(folder + "/t_domains.bed",
                       "29be12b023a79c2bb5b6ffb5b7670d32a982fce3da91a8dd96b829330d0c610a", 550)
    assert_file_pinned(folder + "/t_score.bedgraph",
                       "b525a83e9bf2bdb3a90dfd845b73ce61d6bad175b16db8ea32909a98dbc1b66b", 16750)
    assert_file_pinned(folder + "/t_tad_score.bm",
                       "305f640136e94b5359c1feeb271490454e3bd907fccbc13ba951fd6b81a0fbdf", 16766)

    shutil.rmtree(folder)


def test_find_TADs_default_depth_parameters():
    """--minDepth, --maxDepth and --step left out, so set_variables decides.

    For a 50 kb matrix that is 5 x binsize for minDepth, 10 x for maxDepth and
    2 x for step, which the .bm header records verbatim. Nothing exercised
    these branches before.
    """
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    folder = mkdtemp(prefix="test_find_tads_defaults")
    run_tool("--matrix {} --outPrefix {}/t --correctForMultipleTesting fdr".format(matrix, folder))

    with open(folder + "/t_tad_score.bm") as fh:
        header = json.loads(fh.readline().lstrip('#'))
        columns = len(fh.readline().rstrip('\n').split('\t'))
    assert header == {"step": 100000, "minDepth": 250000, "maxDepth": 500000, "binsize": 50000}
    # get_incremental_step_size(250000, 500000, 100000) yields two window
    # sizes, so three coordinate columns plus two scores.
    assert columns == 5

    assert_file_pinned(folder + "/t_boundaries.bed",
                       "2b2d56df6bb5a6be9a7be2763ebdfc4ad26d1dc3fd33e1afde013817ae6e45b7", 8)
    assert_file_pinned(folder + "/t_boundaries.gff",
                       "c9bf377841ed765b52d8832666ef7d6d7079e33deab725f9b6925f627ed9e6d7", 8)
    assert_file_pinned(folder + "/t_domains.bed",
                       "4a07c5a795c5627bcd1f9b07a1222d57dff3ea0fae552eb1adca9909f452ae08", 3)
    assert_file_pinned(folder + "/t_score.bedgraph",
                       "a479e4f720b8e2b2634ac9d1b658dca9d8ac404a22255db4472faf44365d9ee4", 2436)
    assert_file_pinned(folder + "/t_tad_score.bm",
                       "63afedb6ecebfdae91d994b4548077cd20a3499c71bbccf7dc03e06f967eed2b", 2449)

    shutil.rmtree(folder)


def _call_boundaries(folder, extra):
    """Run the boundary calling only, on the precomputed None spectrum."""
    prepare_precomputed(folder, "t", ("None", "None"))
    run_tool("--matrix {} --minDepth 60000 --maxDepth 180000 --step 20000 "
             "--outPrefix {}/t {}".format(ROOT + "small_test_matrix.h5", folder, extra))


@pytest.mark.parametrize("delta,digest,lines,domain_digest,domain_lines", [
    ("0.0", "b7749f01b624742378ad471bd1054e7dd6a4a3ca6042091285ef5256c9a31ab0", 1523,
     "266cc3dd822d39f2c616ae8c25befed21d41b6161bcb016e9540240ad3b4f7e1", 1518),
    ("0.01", "a94fa71154f4c43f0c508fcf83d3bbb51a4dd8391968efe3d541b0b45cf314fb", 1407,
     "176773444ec9320fc3e03636407d4dbb64620cfeba1bcc25a562b8192867886f", 1402),
    ("0.06", "f61c73095b9fac7aa00480889bda9d575e9f7a84e35fbb436aa983e3838e7b35", 170,
     "b1c727efca5b6045ff42303cf0bac0031a3a32eea4f637fa650796fa2e708334", 165),
])
def test_find_TADs_delta(delta, digest, lines, domain_digest, domain_lines):
    """--delta was never passed by any test although it filters every call."""
    folder = mkdtemp(prefix="test_find_tads_delta")
    _call_boundaries(folder, "--minBoundaryDistance 20000 --correctForMultipleTesting None "
                             "--thresholdComparisons 1.0 --delta " + delta)
    assert_file_pinned(folder + "/t_boundaries.bed", digest, lines)
    assert_file_pinned(folder + "/t_domains.bed", domain_digest, domain_lines)
    shutil.rmtree(folder)


@pytest.mark.parametrize("method,threshold,digest,lines", [
    ("None", "1.0", "a94fa71154f4c43f0c508fcf83d3bbb51a4dd8391968efe3d541b0b45cf314fb", 1407),
    ("None", "0.001", "9f4d941cb82f542d87c9676c0706bc7e74071f8904b204234594f228de2dbeb9", 3),
    ("fdr", "0.5", "9f4d941cb82f542d87c9676c0706bc7e74071f8904b204234594f228de2dbeb9", 3),
    ("bonferroni", "0.001",
     "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", 0),
])
def test_find_TADs_thresholdComparisons(method, threshold, digest, lines):
    """--thresholdComparisons means a q-value for FDR and a p-value otherwise.

    The FDR row is the interesting one: at a q-value of 0.5 the
    Benjamini-Hochberg cutoff computed at hicFindTADs.py:1236-1241 lands on the
    same three boundaries that an uncorrected p-value threshold of 0.001 keeps.
    """
    folder = mkdtemp(prefix="test_find_tads_threshold")
    _call_boundaries(folder, "--minBoundaryDistance 20000 --correctForMultipleTesting {} "
                             "--thresholdComparisons {}".format(method, threshold))
    assert_file_pinned(folder + "/t_boundaries.bed", digest, lines)
    shutil.rmtree(folder)


def test_find_TADs_minBoundaryDistance_sets_the_lookahead():
    """--minBoundaryDistance divided by the bin size is the peak lookahead."""
    folder = mkdtemp(prefix="test_find_tads_lookahead")
    _call_boundaries(folder, "--minBoundaryDistance 100000 --correctForMultipleTesting None "
                             "--thresholdComparisons 1.0")
    assert_file_pinned(folder + "/t_boundaries.bed",
                       "f22007e48bd38a1045397a710e6abe77f4d3597a6641e8e62e6193ee0bf46f33", 477)
    assert_file_pinned(folder + "/t_boundaries.gff",
                       "3aad2cdbe4d6e464f29ab880ed0743cbaa8b3f0fd598114118416735a77c70b5", 477)
    assert_file_pinned(folder + "/t_domains.bed",
                       "f1cb8590b1f083d5c9fcb47cdb2524f12f366e5824c16ab60d6cbcd9fb1e3266", 472)
    shutil.rmtree(folder)


def test_find_TADs_TAD_sep_score_prefix():
    """--TAD_sep_score_prefix, never exercised before.

    The spectrum and the z-score matrix are read from the given prefix, the
    --minDepth, --maxDepth and --step on the command line are replaced by the
    values in the .bm header, and nothing is written next to --outPrefix except
    the four calling outputs.
    """
    source = mkdtemp(prefix="test_find_tads_prefix_src")
    folder = mkdtemp(prefix="test_find_tads_prefix_out")
    prepare_precomputed(source, "pre", ("None", "None"))

    # Deliberately absurd depths: they must be ignored in favour of the header.
    run_tool("--matrix {} --minDepth 999000 --maxDepth 1000000 --step 500000 "
             "--outPrefix {}/t --TAD_sep_score_prefix {}/pre --minBoundaryDistance 20000 "
             "--correctForMultipleTesting None --thresholdComparisons 1.0".format(
                 ROOT + "small_test_matrix.h5", folder, source))

    for suffix in ['_boundaries.bed', '_boundaries.gff', '_domains.bed', '_score.bedgraph']:
        assert_files_identical(ROOT + "find_TADs/None/multiNone" + suffix, folder + "/t" + suffix)
    assert not os.path.exists(folder + "/t_tad_score.bm")
    assert not os.path.exists(folder + "/t_zscore_matrix.h5")

    shutil.rmtree(source)
    shutil.rmtree(folder)


def test_find_TADs_TAD_sep_score_prefix_missing_files_exit():
    folder = mkdtemp(prefix="test_find_tads_prefix_missing")
    with pytest.raises(SystemExit) as exit_info:
        hicFindTADs.main(
            "--matrix {} --outPrefix {}/t --TAD_sep_score_prefix {}/absent "
            "--correctForMultipleTesting fdr".format(
                ROOT + "small_test_matrix.h5", folder, folder).split())
    assert exit_info.value.code == 1
    shutil.rmtree(folder)


def test_find_TADs_maxDepth_must_exceed_minDepth():
    folder = mkdtemp(prefix="test_find_tads_depth_order")
    with pytest.raises(SystemExit):
        hicFindTADs.main(
            "--matrix {} --outPrefix {}/t --correctForMultipleTesting fdr "
            "--minDepth 100000 --maxDepth 100000".format(
                ROOT + "small_test_matrix.h5", folder).split())
    shutil.rmtree(folder)


def test_find_TADs_minDepth_below_three_bins_exits():
    folder = mkdtemp(prefix="test_find_tads_min_depth")
    with pytest.raises(SystemExit) as exit_info:
        hicFindTADs.main(
            "--matrix {} --outPrefix {}/t --correctForMultipleTesting fdr "
            "--minDepth 10000 --maxDepth 100000".format(
                ROOT + "small_test_matrix.h5", folder).split())
    assert exit_info.value.code == 1
    shutil.rmtree(folder)


def test_find_TADs_step_below_bin_size_exits():
    folder = mkdtemp(prefix="test_find_tads_step")
    with pytest.raises(SystemExit) as exit_info:
        hicFindTADs.main(
            "--matrix {} --outPrefix {}/t --correctForMultipleTesting fdr "
            "--minDepth 60000 --maxDepth 180000 --step 1000".format(
                ROOT + "small_test_matrix.h5", folder).split())
    assert exit_info.value.code == 1
    shutil.rmtree(folder)


def test_find_TADs_unreadable_matrix_never_reaches_the_suffix_check():
    """Pinned defect: the "could not determine file ending" branch is dead.

    hicFindTADs.py:1319-1325 checks the suffix of --matrix and exits 1 when it
    is neither cool nor h5, but the HicFindTads constructor at :1313 has
    already loaded the matrix by then, and hicmatrix treats every non-h5 name
    as a cooler. A file that is neither raises OSError out of cooler before the
    check can run, so the intended diagnostic is unreachable. Pinned rather
    than fixed; the C++ port checks the suffix before it opens anything and
    exits 1 with the message the Python only intends to print.
    """
    folder = mkdtemp(prefix="test_find_tads_suffix")
    with pytest.raises(OSError, match="is not an HDF5 file"):
        hicFindTADs.main(
            "--matrix {} --outPrefix {}/t --correctForMultipleTesting fdr".format(
                ROOT + "small_test_matrix.mtx", folder).split())
    shutil.rmtree(folder)


def test_get_incremental_step_size():
    """The window sizes the spectrum is computed at.

    Pinned defect in the documentation, not the code: the --step help text
    promises 20,000, 30,000, 40,000, 70,000 and 100,000 for step=10,000,
    minDepth=20,000 and maxDepth=150,000, but ``min + int(step * x ** 1.5)``
    produces 20,000, 30,000, 48,284, 71,961, 100,000 and 131,803. The port
    reproduces the code, so this test pins the code.
    """
    assert hicFindTADs.get_incremental_step_size(20000, 150000, 10000) == \
        [20000, 30000, 48284, 71961, 100000, 131803]
    assert hicFindTADs.get_incremental_step_size(60000, 180000, 20000) == \
        [60000, 80000, 116568, 163923]
    assert hicFindTADs.get_incremental_step_size(250000, 500000, 100000) == [250000, 350000]
    assert hicFindTADs.get_incremental_step_size(50000, 200000, 10000) == \
        [50000, 60000, 78284, 101961, 130000, 161803, 196969]
    # min == max yields the single starting window.
    assert hicFindTADs.get_incremental_step_size(60000, 60000, 20000) == [60000]


def test_peakdetect_resets_at_a_chromosome_change():
    """The local minimum search must not carry state across chromosomes.

    The signal falls monotonically to the end of the first chromosome and is
    flat and high on the second. Read as one chromosome that is a local
    minimum at index 4; read as two it is the end of a chromosome and must not
    be reported.
    """
    import numpy as np
    values = np.concatenate([np.array([1.0, 0.8, 0.6, 0.4, 0.2]), np.ones(12)])

    _max, _min = hicFindTADs.HicFindTads.peakdetect(
        values, lookahead=3, chrom=np.array(['a'] * 5 + ['b'] * 12))
    assert [list(peak) for peak in _min] == []

    _max, _min = hicFindTADs.HicFindTads.peakdetect(
        values, lookahead=3, chrom=np.array(['a'] * 17))
    assert [list(peak) for peak in _min] == [[4, 0.2]]


def test_delta_wrt_window_window_is_asymmetric_and_includes_the_minimum():
    """Pinned defect: the local window is not what the docstring describes.

    delta_wrt_window says it averages the ``window_len`` scores to the left and
    to the right "the minimum itself is excluded", but hicFindTADs.py:612-613
    concatenates ``matrix_avg[i - w : i + 3]`` with ``matrix_avg[i + 4 : i + w]``.
    That includes the minimum and the two bins after it, drops the bin at
    ``i + 3``, and reaches only ``w - 1`` bins to the right. Reproduced, not
    fixed, because every boundary the tool calls is filtered on this number.
    """
    import numpy as np
    scores = np.arange(80, dtype=float)
    chrom = np.array(['a'] * 40 + ['b'] * 40)
    result = hicFindTADs.HicFindTads.delta_wrt_window([5, 20, 60, 79], scores, chrom,
                                                      window_len=10)
    # Too close to a chromosome border, and the last index is never inside a
    # range because np.unique's trick makes the last range end at len - 1.
    assert np.isnan(result[5])
    assert np.isnan(result[79])

    window = np.concatenate([scores[10:23], scores[24:30]])
    nt.assert_allclose(result[20], window.mean() - scores[20], rtol=0, atol=0)
    nt.assert_allclose(result[20], -0.6842105263157912, rtol=0, atol=1e-15)
    assert not np.isnan(result[60])
