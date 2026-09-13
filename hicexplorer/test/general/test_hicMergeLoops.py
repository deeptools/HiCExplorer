from hicexplorer import hicMergeLoops
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicMergeLoops/")


def are_files_equal(file1, file2, delta=1, skip=0):

    lines_file1_dict = {}
    mismatches = 0
    matches = 0
    line_count_file1 = 0
    with open(file1, 'r') as textfile1:
        file_content = textfile1.readlines()

        for i, line in enumerate(file_content):
            if i < skip:
                continue
            lines_file1_dict[line] = True
            line_count_file1 += 1
    with open(file2, 'r') as textfile2:

        file_content = textfile2.readlines()
        for i, line in enumerate(file_content):
            if i < skip:
                continue
            if line in lines_file1_dict:
                matches += 1
            else:
                mismatches += 1
    if mismatches < delta and line_count_file1 - delta <= matches:
        return True
    else:
        return False


def assert_files_identical(pExpected, pObserved):
    """Line by line equality, in order.

    The set comparison this file used before is blind to the row order, and to
    a duplicated line. The output here is a subset of the bedtools-sorted input
    table written in that table's order, so the order is part of the result.
    """
    with open(pExpected) as expected_file:
        expected = expected_file.readlines()
    with open(pObserved) as observed_file:
        observed = observed_file.readlines()
    assert len(expected) == len(observed), \
        'line count differs: {} against {}'.format(len(expected), len(observed))
    for index, (left, right) in enumerate(zip(expected, observed)):
        assert left == right, 'line {} differs:\n  expected {!r}\n  observed {!r}'.format(
            index + 1, left, right)


def read_lines(pPath):
    with open(pPath) as handle:
        return handle.read().splitlines()


def run(pOutFile, pResolution, *pInputs):
    args = "-i {} -o {} -r {}".format(
        ' '.join(pInputs), pOutFile, pResolution).split()
    hicMergeLoops.main(args)


def test_loop_narrow_peak():
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()

    args = "-i {} {} {} -o {} -r {}".format(ROOT + 'gm12878_10kb.bedgraph', ROOT + 'gm12878_25kb.bedgraph',
                                            ROOT + 'gm12878_5kb.bedgraph', outfile.name, 5000).split()
    hicMergeLoops.main(args)

    assert_files_identical(ROOT + 'gm12878_all.bedgraph', outfile.name)


def test_result_is_a_subsequence_of_the_sorted_input():
    """The output rows keep the order and the text of the merged input table.

    Reconstructed here rather than taken from the master, so that the ordering
    assertion does not depend on the checked-in file: the three inputs are
    concatenated, sorted with bedtools, stripped of every duplicated row by
    drop_duplicates(keep=False), and the tool's output must be a subsequence of
    that.
    """
    import pandas as pd
    from pybedtools import BedTool

    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()
    inputs = [ROOT + 'gm12878_10kb.bedgraph', ROOT + 'gm12878_25kb.bedgraph',
              ROOT + 'gm12878_5kb.bedgraph']
    run(outfile.name, 5000, *inputs)

    frame = None
    for path in inputs:
        part = pd.read_csv(path, sep='\t', header=None)
        frame = part if frame is None else pd.concat([frame, part])
    frame = BedTool.from_dataframe(frame).sort().to_dataframe(
        disable_auto_names=True, header=None)
    frame.drop_duplicates(keep=False, inplace=True)
    table = frame.to_csv(sep='\t', header=False, index=False).splitlines()

    merged = read_lines(outfile.name)
    position = 0
    for line in merged:
        while position < len(table) and table[position] != line:
            position += 1
        assert position < len(table), \
            'merged line not found in the sorted table, or out of order: ' + line
        position += 1
    nt.assert_equal(len(merged), 16549)


def test_single_file_drops_only_its_own_overlaps():
    """One input file, which no test covered.

    With a single 10 kb file every loop has the same anchor width, so the merge
    can only drop a loop whose anchors both overlap another loop's. 12,320 loops
    go in and 12,317 come out. It also pins that the tool works at all with one
    --inputFiles argument, where the concat branch is never taken.
    """
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()
    run(outfile.name, 10000, ROOT + 'gm12878_10kb.bedgraph')

    merged = read_lines(outfile.name)
    nt.assert_equal(len(merged), 12317)

    # The coordinates must all come from the input. The score column is
    # deliberately excluded: the frame goes through pandas' read_csv twice,
    # whose default float converter keeps 17 significant digits and then
    # scales, so 0.0019430592210407135 in the input is written back as
    # 0.0019430592210407. That is the reference's behaviour and is pinned in
    # its own test below rather than folded into this one.
    source = set(tuple(line.split('\t')[:6])
                 for line in read_lines(ROOT + 'gm12878_10kb.bedgraph'))
    for line in merged:
        anchors = tuple(line.split('\t')[:6])
        assert anchors in source, \
            'merged line is not one of the input lines: ' + line


def test_score_column_is_reformatted_by_the_pandas_float_parser():
    """A 17+ digit score does not survive the round trip unchanged.

    pandas' C parser converts a float with precise_xstrtod, which accumulates
    at most 17 significant decimal digits and then scales by a power of ten
    rather than rounding correctly, and to_csv writes repr() of the result. The
    10 kb input holds 0.0019430592210407135; the merged output holds
    0.0019430592210407. Pinned because it is the sort of difference a port
    would otherwise 'fix' by using strtod, and because it is the reason the
    coordinate check above ignores the score column.
    """
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()
    run(outfile.name, 10000, ROOT + 'gm12878_10kb.bedgraph')

    assert any(line.endswith('\t0.0019430592210407135')
               for line in read_lines(ROOT + 'gm12878_10kb.bedgraph'))
    merged = read_lines(outfile.name)
    assert any(line.endswith('\t0.0019430592210407') for line in merged)
    assert not any(line.endswith('\t0.0019430592210407135') for line in merged)


def test_lowest_resolution_changes_the_result_monotonically():
    """--lowestResolution is the search radius, and no test varied it.

    The neighbourhood factor is `lowestResolution - anchor width`
    (hicMergeLoops.py:72-75), so a larger value widens the search and can only
    merge more. The counts are pinned exactly, and their monotonicity is
    asserted as well so that a port which inverted the sign of the factor, or
    ignored the option, fails here rather than in one specific number.
    """
    inputs = [ROOT + 'gm12878_10kb.bedgraph', ROOT + 'gm12878_25kb.bedgraph',
              ROOT + 'gm12878_5kb.bedgraph']
    counts = {}
    for resolution in (5000, 10000, 25000, 50000):
        outfile = NamedTemporaryFile(suffix='out', delete=True)
        outfile.close()
        run(outfile.name, resolution, *inputs)
        counts[resolution] = len(read_lines(outfile.name))

    nt.assert_equal(counts, {5000: 16549, 10000: 16242, 25000: 15908, 50000: 15602})
    assert counts[5000] > counts[10000] > counts[25000] > counts[50000]


def test_duplicated_rows_are_dropped_entirely():
    """drop_duplicates(keep=False) removes every copy, not all but one.

    A row that appears in two of the input files disappears from the merge
    completely (hicMergeLoops.py:165). The 5 kb and the 25 kb file share no
    rows, so the case is built here by passing the same file twice: every one of
    its rows is then a duplicate, the frame empties, and the tool writes an
    empty file rather than the file's own content.
    """
    outfile = NamedTemporaryFile(suffix='out', delete=True)
    outfile.close()
    run(outfile.name, 25000, ROOT + 'gm12878_25kb.bedgraph',
        ROOT + 'gm12878_25kb.bedgraph')

    nt.assert_equal(os.path.getsize(outfile.name), 0)
