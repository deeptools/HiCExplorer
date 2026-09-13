"""Characterization tests for hicDetectLoops.

Written before the C++ port of this tool (cpp/AGENTS_CONTRACT.md rule 1), and
deliberately not written the way the four tests it replaces were.

What the previous file constrained, honestly:

  * `test_main_h5` called the tool and asserted nothing at all.
  * The three cool tests all passed `--maxLoopDistance 3000000` against a
    matrix binned at 2,500,000 bp. hicmatrix loads a distance restricted
    chromosome by keeping pixels with `bin2_id - bin1_id < distance //
    binsize`, which is `< 1` here, so only the main diagonal survives; the
    tool then deletes the main diagonal and has nothing left. All three runs
    produce **no loops and no output file**, and the committed reference
    `hicDetectLoops/loops.bedgraph` (one line) is never reproduced.
  * They still passed, because `are_files_equal` compares the two files with
    `zip`, which stops at the shorter one. An empty output file is therefore
    equal to every reference. The old assertions could not fail.

So the whole file was vacuous. It is replaced by:

  * full precision assertions on the exact bedgraph the tool writes, for the
    seven option combinations below, since every value in it is stable to the
    last digit across repeated runs;
  * `test_are_files_equal_cannot_fail_on_a_truncated_file` and
    `test_committed_reference_loops_bedgraph_is_not_reproduced`, which pin the
    two defects above so that they are not rediscovered as regressions;
  * `test_output_order_depends_on_scheduling`, which pins that the tool's line
    order is not reproducible while its call set is;
  * coverage of `--expected` (all three values), `--obsExpThreshold`,
    `--pValuePreselection` as a threshold file, `--peakWidth`/`--windowSize`,
    `--maxLoopDistance` and the two failure exits.

Every run passes `--threads 1 --threadsPerChromosome 1`, except the one test
that is about thread dependence, because that is the only setting under which
the reference's own output order is defined.
"""
import os.path
import subprocess
import sys
from tempfile import NamedTemporaryFile, TemporaryDirectory
import logging

import pytest

from hicexplorer import hicDetectLoops
from hicexplorer.test.test_compute_function import compute

log = logging.getLogger(__name__)

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")

MATRIX_H5 = ROOT + "small_test_matrix.h5"
MATRIX_COOL = ROOT + "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool"
THRESHOLD_FILE = ROOT + "hicCreateThresholdFile/thresholdFile_loose_pValue.txt"

# The bin size of MATRIX_COOL. Named because three of the tests below only make
# sense next to it.
COOL_BIN_SIZE = 2500000


def are_files_equal(file1, file2, delta=None):
    """Kept for the two tests that pin its behaviour. Do not use it to assert.

    It compares the two files with zip, so a file that is a prefix of the
    other compares equal, and an empty file compares equal to anything. Note
    also that `delta=0` is falsy, so the `delta` branch is dead for the value
    the old tests passed.
    """
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


def run_loops(arguments, expect_exit=None):
    """Run the tool and return the lines it wrote, or None when it wrote no file.

    hicDetectLoops writes the output file only when it found at least one loop
    (hicDetectLoops.py:1090), so "no file" is a distinct, observable outcome
    from "an empty file" and the tests keep them apart.
    """
    with TemporaryDirectory() as directory:
        out = os.path.join(directory, "loops.bedgraph")
        args = arguments.format(out=out).split()
        if expect_exit is None:
            compute(hicDetectLoops.main, args, 5)
        else:
            # main() calls exit(1) on failure, which compute() would propagate
            # as SystemExit; catching it here keeps the assertion explicit.
            with pytest.raises(SystemExit) as raised:
                hicDetectLoops.main(args)
            assert raised.value.code == expect_exit
        if not os.path.exists(out):
            return None
        with open(out) as handle:
            return handle.read().splitlines()


def assert_lines(actual, expected_text):
    expected = [line for line in expected_text.strip().split("\n")]
    assert actual is not None, "the tool wrote no output file"
    assert len(actual) == len(expected), (
        "expected %d loops, got %d" % (len(expected), len(actual)))
    for index, (got, want) in enumerate(zip(actual, expected)):
        assert got == want, "line %d:\n  got  %r\n  want %r" % (index + 1, got, want)


# ---------------------------------------------------------------------------
# The two defects in the tests this file replaces


def test_are_files_equal_cannot_fail_on_a_truncated_file():
    """zip() stops at the shorter file, so a prefix compares equal.

    This is why the three cool tests in the previous version of this file
    passed while producing no loops at all.
    """
    with TemporaryDirectory() as directory:
        full = os.path.join(directory, "full.bedgraph")
        empty = os.path.join(directory, "empty.bedgraph")
        with open(full, "w") as handle:
            handle.write("chr1\t1\t2\tchr1\t3\t4\t0.5\n"
                         "chr1\t5\t6\tchr1\t7\t8\t0.25\n")
        with open(empty, "w") as handle:
            handle.write("")
        assert are_files_equal(full, empty, delta=0) is True
        # And a one line prefix is equal to the two line file as well.
        prefix = os.path.join(directory, "prefix.bedgraph")
        with open(prefix, "w") as handle:
            handle.write("chr1\t1\t2\tchr1\t3\t4\t0.5\n")
        assert are_files_equal(full, prefix, delta=0) is True
        # It does detect a difference in a line that is present in both.
        wrong = os.path.join(directory, "wrong.bedgraph")
        with open(wrong, "w") as handle:
            handle.write("chr1\t1\t2\tchr1\t3\t4\t0.6\n")
        assert are_files_equal(full, wrong, delta=0) is False


def test_maxLoopDistance_below_the_bin_size_finds_nothing():
    """--maxLoopDistance 3000000 on a 2.5 Mb matrix leaves only the diagonal.

    hicmatrix keeps `bin2_id - bin1_id < maxLoopDistance // binsize` bins,
    which is `< 1`, so the loaded chromosome is the main diagonal alone.
    hicDetectLoops then deletes the main diagonal (hicDetectLoops.py:861-865)
    and the matrix is empty. This is the configuration the previous version of
    this file used for all three of its cool tests.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --maxLoopDistance 3000000 -pit 1 -w 5 -pw 2 "
        "-p 0.5 -pp 0.55 --chromosomes 1 2 -t 1 -tpc 1".format(MATRIX_COOL))
    assert lines is None


def test_committed_reference_loops_bedgraph_is_not_reproduced():
    """test_data/hicDetectLoops/loops.bedgraph is stale.

    It holds one loop at a distance of exactly one 2.5 Mb bin, which the
    current code cannot produce from this matrix at any --maxLoopDistance:
    below 5000000 the loader keeps nothing off the diagonal, and from 5000000
    upwards the loop is not called. Kept as a test rather than deleted so that
    the file is not mistaken for a valid expectation again.
    """
    with open(ROOT + "hicDetectLoops/loops.bedgraph") as handle:
        reference = handle.read().splitlines()
    assert reference == ["1\t112500000\t115000000\t1\t115000000\t117500000\t0.001"]
    lines = run_loops(
        "--matrix {} -o {{out}} --maxLoopDistance 50000000 -pit 1 -w 5 -pw 2 "
        "-p 0.5 -pp 0.55 --chromosomes 1 -t 1 -tpc 1".format(MATRIX_COOL))
    assert lines is not None
    assert reference[0] not in lines


# ---------------------------------------------------------------------------
# Full precision expectations, one per option path


CHR2L_DEFAULT_EXPECTED = """
chr2L	485000	490000	chr2L	490000	495000	0.4385623894062053
chr2L	960000	965000	chr2L	1045000	1050000	0.4016704795074202
chr2L	1190000	1195000	chr2L	1195000	1200000	0.21901405427555742
chr2L	2950000	2955000	chr2L	2965000	2970000	0.09730180547398455
chr2L	3385000	3390000	chr2L	3390000	3395000	0.3179460186301115
chr2L	4420000	4425000	chr2L	4430000	4435000	0.3565987181527521
chr2L	5200000	5205000	chr2L	5205000	5210000	0.24787229188765736
chr2L	5260000	5265000	chr2L	5265000	5270000	0.14527668564350318
chr2L	7475000	7480000	chr2L	7480000	7485000	0.12442987548470683
chr2L	8675000	8680000	chr2L	8680000	8685000	0.3133150561644056
chr2L	8960000	8965000	chr2L	9185000	9190000	0.3565987181527521
chr2L	8965000	8970000	chr2L	9195000	9200000	0.3565987181527521
chr2L	8970000	8975000	chr2L	9195000	9200000	0.3565987181527521
chr2L	9805000	9810000	chr2L	9820000	9825000	0.44616777885985415
chr2L	13800000	13805000	chr2L	13910000	13915000	0.3565987181527521
chr2L	14755000	14760000	chr2L	14855000	14860000	0.3565987181527521
chr2L	16340000	16345000	chr2L	16345000	16350000	0.3565987181527521
chr2L	17005000	17010000	chr2L	17020000	17025000	0.44616777885985415
chr2L	17045000	17050000	chr2L	17120000	17125000	0.3565987181527521
chr2L	18670000	18675000	chr2L	18675000	18680000	0.23000508888575644
chr2L	18940000	18945000	chr2L	19030000	19035000	0.3565987181527521
chr2L	19585000	19590000	chr2L	19670000	19675000	0.3565987181527521
chr2L	19755000	19760000	chr2L	19875000	19880000	0.3565987181527521
chr2L	20795000	20800000	chr2L	20810000	20815000	0.35827102425110247
chr2L	22030000	22035000	chr2L	22160000	22165000	0.3565987181527521
"""


def test_h5_single_chromosome_full_precision():
    """The default --expected mean on an integer count matrix.

    Note what obs_exp_matrix does to an integer matrix: it records the input
    dtype before casting the data to float32 and casts the quotient back to it
    (hicexplorer/utilities.py:582,588), so the observed over expected ratio is
    truncated to an integer. It only survives that because the expected value
    per distance here is around 0.05, so the ratios are in the hundreds.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.5 "
        "-t 1 -tpc 1".format(MATRIX_H5))
    assert_lines(lines, CHR2L_DEFAULT_EXPECTED)


def test_h5_all_chromosomes_full_precision():
    """The case the previous test_main_h5 ran without asserting anything."""
    lines = run_loops(
        "--matrix {} -o {{out}} -pit 1 -p 0.5 -pp 0.5 -t 1 -tpc 1".format(
            MATRIX_H5))
    assert lines is not None
    assert len(lines) == 52
    # Only two of the 15 chromosomes produce any loop at all, and they come out
    # in bin table order, which puts chr3L before chr2L.
    assert [line.split("\t")[0] for line in lines] == ["chr3L"] * 27 + ["chr2L"] * 25
    assert lines[0] == ("chr3L\t125000\t130000\tchr3L\t135000\t140000\t"
                        "0.3633174801431127")
    assert lines[26] == ("chr3L\t21545000\t21550000\tchr3L\t21550000\t21555000\t"
                         "0.16869892359210836")
    # The chr2L block is exactly the single chromosome run, line for line.
    assert_lines(lines[27:], CHR2L_DEFAULT_EXPECTED)


CHR2L_MEAN_NONZERO_EXPECTED = """
chr2L	485000	490000	chr2L	490000	495000	0.43478781014221823
chr2L	960000	965000	chr2L	1045000	1050000	0.4016704795074202
chr2L	1190000	1195000	chr2L	1195000	1200000	0.21901405427555742
chr2L	2950000	2955000	chr2L	2965000	2970000	0.09730180547398455
chr2L	3385000	3390000	chr2L	3390000	3395000	0.3133150561644056
chr2L	4000000	4005000	chr2L	4005000	4010000	0.28936977075534764
chr2L	4360000	4365000	chr2L	4370000	4375000	0.21901405427555742
chr2L	4420000	4425000	chr2L	4430000	4435000	0.3565987181527521
chr2L	5260000	5265000	chr2L	5265000	5270000	0.14527668564350318
chr2L	6100000	6105000	chr2L	6110000	6115000	0.2531510616226985
chr2L	7475000	7480000	chr2L	7480000	7485000	0.12442987548470683
chr2L	8430000	8435000	chr2L	8435000	8440000	0.17472055277458887
chr2L	8675000	8680000	chr2L	8680000	8685000	0.31178131839409173
chr2L	9805000	9810000	chr2L	9820000	9825000	0.44808079459400096
chr2L	10350000	10355000	chr2L	10355000	10360000	0.3667063700481331
chr2L	12400000	12405000	chr2L	12405000	12410000	0.3752644608113599
chr2L	14755000	14760000	chr2L	14855000	14860000	0.3565987181527521
chr2L	16340000	16345000	chr2L	16345000	16350000	0.3565987181527521
chr2L	17005000	17010000	chr2L	17020000	17025000	0.44616777885985415
chr2L	18670000	18675000	chr2L	18675000	18680000	0.22752936411075952
chr2L	19585000	19590000	chr2L	19670000	19675000	0.3565987181527521
chr2L	20255000	20260000	chr2L	20265000	20270000	0.14794825611385778
chr2L	20830000	20835000	chr2L	20840000	20845000	0.17472055277458887
chr2L	21625000	21630000	chr2L	21630000	21635000	0.24787229188765736
chr2L	21630000	21635000	chr2L	21635000	21640000	0.31178131839409173
chr2L	21635000	21640000	chr2L	21640000	21645000	0.277871806255898
chr2L	22030000	22035000	chr2L	22160000	22165000	0.3565987181527521
"""


def test_h5_expected_mean_nonzero():
    """--expected mean_nonzero, which no test has ever run.

    obs_exp_matrix_non_zero assigns each quotient back into an array it has
    already cast to float32 and never casts back (utilities.py:533,540), so
    unlike the default path its result is float32 rather than the input dtype.
    Its epsilon for a NaN or infinite quotient is 1e-9, not the 1e-6 the
    default path uses.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L --expected mean_nonzero "
        "-pit 1 -p 0.5 -pp 0.5 -t 1 -tpc 1".format(MATRIX_H5))
    assert_lines(lines, CHR2L_MEAN_NONZERO_EXPECTED)


def test_h5_expected_mean_nonzero_ligation():
    """--expected mean_nonzero_ligation, also never run before.

    The Homer style correction multiplies the expected value by
    row_sum(i) * row_sum(j) / total, which is zero whenever either bin has an
    empty row, so the quotient is infinite and is replaced by the 1e-9 epsilon.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L "
        "--expected mean_nonzero_ligation -pit 1 -p 0.5 -pp 0.5 "
        "-t 1 -tpc 1".format(MATRIX_H5))
    assert lines is not None
    assert len(lines) == 62
    assert lines[0] == ("chr2L\t220000\t225000\tchr2L\t240000\t245000\t"
                        "0.4016704795074202")
    assert lines[-1] == ("chr2L\t22035000\t22040000\tchr2L\t22170000\t22175000\t"
                         "0.398086603716326")
    # Every line is a chr2L intra chromosomal pair with seven fields and a
    # p-value at or below the --pValue of 0.5.
    for line in lines:
        fields = line.split("\t")
        assert len(fields) == 7
        assert fields[0] == "chr2L" and fields[3] == "chr2L"
        assert 0.0 <= float(fields[6]) <= 0.5


CHR2L_WINDOW7_EXPECTED = """
chr2L	485000	490000	chr2L	490000	495000	0.30371216516921395
chr2L	950000	955000	chr2L	1105000	1110000	0.38243248357483484
chr2L	2950000	2955000	chr2L	2965000	2970000	0.14387965556669857
chr2L	4360000	4365000	chr2L	4370000	4375000	0.38243248357483484
chr2L	4370000	4375000	chr2L	4375000	4380000	0.4164297080385082
chr2L	7475000	7480000	chr2L	7480000	7485000	0.2749302852924145
chr2L	15825000	15830000	chr2L	15980000	15985000	0.41429834009735467
chr2L	17005000	17010000	chr2L	17020000	17025000	0.45364184883329106
chr2L	17280000	17285000	chr2L	17290000	17295000	0.38243248357483484
chr2L	21185000	21190000	chr2L	21205000	21210000	0.2679017767822581
chr2L	21245000	21250000	chr2L	21275000	21280000	0.23510553107318533
"""


def test_h5_window_size_and_peak_width():
    """-w 7 -pw 3, neither of which the previous file varied on the h5 input.

    A wider peak makes the negative slice bounds of candidate_region_test
    reachable: for a candidate within --peakWidth of the edge of its window,
    `neighborhood[:peak_row - peakWidth, :]` has a negative stop and numpy
    reads that as counting from the end, not as an empty slice.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -w 7 -pw 3 "
        "-pit 1 -p 0.5 -pp 0.5 -t 1 -tpc 1".format(MATRIX_H5))
    assert_lines(lines, CHR2L_WINDOW7_EXPECTED)


def test_h5_tight_preselection_changes_the_call_set():
    """--pValuePreselection is what the negative binomial fit is used for.

    At 0.5 nearly every pixel above --obsExpThreshold survives the tail test
    and the calls are decided downstream; at 0.01 the fitted parameters choose
    the candidates. Pinning both is what makes the fit itself observable.
    """
    loose = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.5 "
        "-t 1 -tpc 1".format(MATRIX_H5))
    tight = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.01 "
        "-t 1 -tpc 1".format(MATRIX_H5))
    assert len(loose) == 25
    assert len(tight) == 14
    assert set(tight).issubset(set(loose))
    assert tight[0] == ("chr2L\t485000\t490000\tchr2L\t490000\t495000\t"
                        "0.4385623894062053")
    assert tight[-1] == ("chr2L\t22030000\t22035000\tchr2L\t22160000\t22165000\t"
                         "0.3565987181527521")


def test_h5_threshold_file_preselection():
    """--pValuePreselection given as a file, the read_threshold_file branch.

    The threshold is looked up by `distance_in_bins * binSize`, so a distance
    the file does not cover raises KeyError inside a worker and the run fails.
    The committed threshold file stops at 200000, which is why
    --maxLoopDistance is capped here.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L --maxLoopDistance 200000 "
        "-pp {} -pit 1 -p 0.5 -t 1 -tpc 1".format(MATRIX_H5, THRESHOLD_FILE))
    assert lines is not None
    assert len(lines) == 22
    assert lines[0] == ("chr2L\t485000\t490000\tchr2L\t490000\t495000\t"
                        "0.4385623894062053")
    assert lines[-1] == ("chr2L\t22030000\t22035000\tchr2L\t22160000\t22165000\t"
                         "0.3565987181527521")


def test_h5_max_loop_distance_filters_the_calls():
    """--maxLoopDistance at 20 bins rather than 400.

    On an h5 input the cut is `distances > maxLoopDistance / binSize`, a float
    division with a non strict comparison (hicDetectLoops.py:820-823), which is
    not the rule the cool loader applies to the same option.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L --maxLoopDistance 100000 "
        "-pit 1 -p 0.5 -pp 0.5 -t 1 -tpc 1".format(MATRIX_H5))
    assert lines is not None
    assert len(lines) == 19
    for line in lines:
        fields = line.split("\t")
        assert abs(int(fields[1]) - int(fields[4])) <= 100000


def test_h5_variable_bin_sizes():
    """Li_et_al_2015.h5: 11,104 restriction fragment bins, not fixed width.

    The largest call set this tool produces anywhere in the corpus, 832 loops,
    and the only input where getBinSize() is a median over bins that genuinely
    differ, so the bin distance cut and the genomic distance filter in
    cluster_to_genome_position_mapping are not the same rule.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --maxLoopDistance 100000 -pit 1 -p 0.5 -pp 0.5 "
        "-t 1 -tpc 1".format(ROOT + "Li_et_al_2015.h5"))
    assert lines is not None
    assert len(lines) == 832
    assert lines[0] == "X\t43223\t44570\tX\t69085\t70854\t2.286311894128384e-08"
    assert lines[1] == "X\t77008\t78499\tX\t78499\t79559\t0.18090119770452806"
    assert lines[-1] == ("X\t22384424\t22386244\tX\t22409588\t22411768\t"
                         "3.591168242233997e-12")
    # Every call is inside --maxLoopDistance measured on the genomic starts,
    # which is what cluster_to_genome_position_mapping enforces. On a fixed
    # width matrix that filter is redundant because the bin distance cut is
    # already stricter; here the bins are 1,843 bp on median but range widely,
    # so it is the filter that decides.
    for line in lines:
        fields = line.split("\t")
        assert abs(int(fields[1]) - int(fields[4])) <= 100000


def test_h5_obs_exp_threshold_and_peak_interactions_threshold():
    """--obsExpThreshold and --peakInteractionsThreshold, both able to empty
    the result.

    --obsExpThreshold gates which pixels get a p-value at all; on this integer
    obs/exp matrix the ratios run into the thousands, so it takes a value in
    the hundreds to bite. --peakInteractionsThreshold gates the raw counts and
    the default of 10 is already enough to find nothing here, which is why
    every other test in this file passes -pit 1.
    """
    # 5 is below every obs/exp value at every distance, so it changes nothing.
    unchanged = run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.5 "
        "-oet 5 -t 1 -tpc 1".format(MATRIX_H5))
    assert_lines(unchanged, CHR2L_DEFAULT_EXPECTED)
    # 500 removes every candidate, and with a single chromosome the tool then
    # exits 1 rather than writing an empty file.
    assert run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.5 "
        "-oet 500 -t 1 -tpc 1".format(MATRIX_H5), expect_exit=1) is None
    # The stock defaults find nothing on this matrix either.
    assert run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -t 1 -tpc 1".format(
            MATRIX_H5), expect_exit=1) is None
    # --peakInteractionsThreshold on its own: 1 gives the 25 loops above and
    # 20 gives none, so this pins that the raw count mask is applied at all.
    # Without it the 20 run would produce the same 25 loops.
    assert run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -pit 20 -p 0.5 -pp 0.5 "
        "-t 1 -tpc 1".format(MATRIX_H5), expect_exit=1) is None


def test_window_size_must_exceed_peak_width():
    """main() rejects --windowSize <= --peakWidth, and does it first.

    The matrix deliberately does not exist: the check runs before anything is
    loaded, so a correct implementation exits 1 without touching it, while one
    that only rejects --windowSize < --peakWidth walks on into the loader and
    raises. Testing it against a real matrix would not distinguish the two,
    because -w 2 -pw 2 finds no loops on any matrix in the corpus and the tool
    exits 1 for that reason instead.
    """
    with TemporaryDirectory() as directory:
        out = os.path.join(directory, "loops.bedgraph")
        args = ("--matrix /nonexistent/matrix.h5 -o {} -w 2 -pw 2 "
                "-t 1 -tpc 1".format(out)).split()
        with pytest.raises(SystemExit) as raised:
            hicDetectLoops.main(args)
        assert raised.value.code == 1
        assert not os.path.exists(out)
    # And --windowSize above --peakWidth is accepted, so the check is not
    # simply rejecting everything.
    assert run_loops(
        "--matrix {} -o {{out}} --chromosomes chr2L -w 5 -pw 2 -pit 1 "
        "-p 0.5 -pp 0.5 -t 1 -tpc 1".format(MATRIX_H5)) is not None


# ---------------------------------------------------------------------------
# The cool loader


COOL_CHR1_EXPECTED = """
1	5000000	7500000	1	17500000	20000000	0.004067531885338184
1	120000000	122500000	1	145000000	147500000	0.0788656354876193
1	207500000	210000000	1	230000000	232500000	8.805294560913203e-07
1	212500000	215000000	1	245000000	247500000	0.031624120627730484
"""

COOL_CHR1_CHR2_EXPECTED = COOL_CHR1_EXPECTED.strip() + """
2	10000000	12500000	2	27500000	30000000	0.1351009131589282
2	25000000	27500000	2	27500000	30000000	0.48530165636498557
2	27500000	30000000	2	45000000	47500000	0.1324698253645566
2	72500000	75000000	2	85000000	87500000	0.008350491904646425
2	120000000	122500000	2	127500000	130000000	0.01622154620376431
"""


def test_cool_single_chromosome_full_precision():
    """The cool loader with a workable --maxLoopDistance.

    50000000 is 20 bins of 2.5 Mb. Below 5000000 this matrix yields nothing at
    all; see test_maxLoopDistance_below_the_bin_size_finds_nothing.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes 1 --maxLoopDistance 50000000 "
        "-pit 1 -w 5 -pw 2 -p 0.5 -pp 0.55 -t 1 -tpc 1".format(MATRIX_COOL))
    assert_lines(lines, COOL_CHR1_EXPECTED)


def test_cool_two_chromosomes_full_precision():
    """Two chromosomes, so the tool leaves its single core branch.

    At --threads 1 only one worker is in flight, so the line order is the
    order of --chromosomes and the comparison can be made line for line.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes 1 2 --maxLoopDistance 50000000 "
        "-pit 1 -w 5 -pw 2 -p 0.5 -pp 0.55 -t 1 -tpc 1".format(MATRIX_COOL))
    assert_lines(lines, COOL_CHR1_CHR2_EXPECTED)


def test_cool_expected_mean_nonzero():
    """The cool loader together with the non zero expected value.

    Both the loaded matrix and the obs/exp matrix are then float32, which is
    what makes fit_nbinom evaluate two of the five terms of its log likelihood
    in single precision.
    """
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes 1 --maxLoopDistance 50000000 "
        "--expected mean_nonzero -pit 1 -w 5 -pw 2 -p 0.5 -pp 0.55 "
        "-t 1 -tpc 1".format(MATRIX_COOL))
    assert lines is not None
    assert len(lines) == 11
    assert lines[0] == ("1\t5000000\t7500000\t1\t17500000\t20000000\t"
                        "0.001529886778055783")
    assert lines[2] == ("1\t52500000\t55000000\t1\t92500000\t95000000\t"
                        "4.178127655870946e-07")
    assert lines[-1] == ("1\t200000000\t202500000\t1\t205000000\t207500000\t"
                         "0.026757454806728666")


def test_cool_tight_preselection():
    lines = run_loops(
        "--matrix {} -o {{out}} --chromosomes 1 --maxLoopDistance 50000000 "
        "-pit 1 -w 5 -pw 2 -p 0.5 -pp 0.1 -t 1 -tpc 1".format(MATRIX_COOL))
    assert lines is not None
    assert len(lines) == 3
    assert lines == [
        "1\t5000000\t7500000\t1\t17500000\t20000000\t0.004067531885338184",
        "1\t120000000\t122500000\t1\t145000000\t147500000\t0.0788656354876193",
        "1\t212500000\t215000000\t1\t245000000\t247500000\t0.031624120627730484",
    ]


# ---------------------------------------------------------------------------
# Thread dependence


def test_threadsPerChromosome_does_not_change_the_result():
    """--threadsPerChromosome only partitions independent work.

    Every stage it parallelises is per genomic distance or per candidate, and
    the partial results are concatenated in worker index order, so the answer
    is the same at any value.
    """
    outputs = []
    for tpc in (1, 4):
        outputs.append(run_loops(
            "--matrix {} -o {{out}} --chromosomes chr2L -pit 1 -p 0.5 -pp 0.5 "
            "-t 1 -tpc {}".format(MATRIX_H5, tpc)))
    assert outputs[0] == outputs[1]
    assert_lines(outputs[0], CHR2L_DEFAULT_EXPECTED)


@pytest.mark.parametrize("threads", [1, 4])
def test_output_order_depends_on_scheduling(threads):
    """The call set is reproducible; the line order is not.

    main() starts one process per chromosome and appends each result to the
    output as that process happens to finish (hicDetectLoops.py:1060-1066), so
    with more than one worker the order of the 52 loops is whatever the
    scheduler produced. Measured on this matrix: three runs at --threads 4
    gave two distinct orderings, and --threads 8 --threadsPerChromosome 4 gave
    a third. This test asserts the invariant that does hold, and records that
    the stronger one does not.
    """
    runs = [run_loops(
        "--matrix {} -o {{out}} -pit 1 -p 0.5 -pp 0.5 -t {} -tpc 4".format(
            MATRIX_H5, threads)) for _ in range(2)]
    assert sorted(runs[0]) == sorted(runs[1])
    assert len(runs[0]) == 52
    if threads == 1:
        # With one worker the order is defined: chromosomes in bin table order.
        assert runs[0] == runs[1]
