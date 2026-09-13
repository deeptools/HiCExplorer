import os.path
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import logging
log = logging.getLogger(__name__)

from hicexplorer import hicCreateThresholdFile
from hicexplorer._version import __version__

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")


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


def read_lines(pPath):
    with open(pPath) as handle:
        return handle.read().splitlines()


def assert_body_identical(pExpected, pObserved):
    """Everything below the two header lines, in order and with the count.

    The old comparison walked the files with zip, which stops at the shorter of
    the two, and tolerated one differing line. A file truncated to a single row
    passed it. The first header line carries the HiCExplorer version and the
    master was written by 3.5-dev, so it is compared separately and by shape.
    """
    expected = read_lines(pExpected)
    observed = read_lines(pObserved)
    assert len(expected) == len(observed), \
        'line count differs: {} against {}'.format(len(expected), len(observed))
    assert observed[0] == \
        "# Threshold file of HiCExplorer's hicCreateThresholdFile version " + __version__
    assert expected[1] == observed[1]
    for index, (left, right) in enumerate(zip(expected[2:], observed[2:]), start=3):
        assert left == right, 'line {} differs: {!r} against {!r}'.format(
            index, left, right)


def test_main():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)

    args = "--range {} {} -tv {} -o {}".format(
        200000, 200000, 0.5, outfile.name).split()
    hicCreateThresholdFile.main(args)
    assert_body_identical(
        ROOT + "hicCreateThresholdFile/thresholdFile_loose_pValue.txt", outfile.name)


def test_resolution():
    """--resolution/-r, which no test exercised.

    It is both the step of the generated rows and part of the stop value, since
    the loop runs to `range[1] + resolution` exclusive: with --range -20000
    20000 and -r 10000 the rows run from -20000 to 20000 inclusive, five of
    them. A port that used `range[1]` as the stop, or that stepped without
    adding the resolution to it, drops the last row.
    """
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()

    hicCreateThresholdFile.main("--range -20000 20000 -tv 0.25 -r 10000 -o {}".format(
        outfile.name).split())

    assert read_lines(outfile.name) == [
        "# Threshold file of HiCExplorer's hicCreateThresholdFile version " + __version__,
        '# Standard threshold 0.25',
        '-20000\t0.25',
        '-10000\t0.25',
        '0\t0.25',
        '10000\t0.25',
        '20000\t0.25',
    ]
    os.unlink(outfile.name)


def test_positive_upstream_is_negated():
    """--range takes the upstream value positive or negative, to the same end.

    hicCreateThresholdFile.py:50-51 negates a positive first value in place, so
    `--range 20000 20000` and `--range -20000 20000` must produce the same file.
    Never tested, and the only test in the file happened to pass a positive
    value, so a port that dropped the negation would have passed it as well
    while producing a single row.
    """
    positive = NamedTemporaryFile(suffix='.txt', delete=False)
    positive.close()
    negative = NamedTemporaryFile(suffix='.txt', delete=False)
    negative.close()

    hicCreateThresholdFile.main(
        "--range 20000 20000 -tv 0.25 -r 10000 -o {}".format(positive.name).split())
    hicCreateThresholdFile.main(
        "--range -20000 20000 -tv 0.25 -r 10000 -o {}".format(negative.name).split())

    assert read_lines(positive.name) == read_lines(negative.name)
    assert len(read_lines(positive.name)) == 7

    os.unlink(positive.name)
    os.unlink(negative.name)


def test_threshold_value_is_written_as_a_python_float():
    """The threshold is a float and is formatted with str(), not with %g.

    -tv 1 becomes 1.0 in both the header and every row, and -tv 1e-05 keeps
    Python's exponent spelling. That is `'{}'.format(float)`, which is the
    shortest round-tripping representation, and it is the only piece of
    formatting this tool has.
    """
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()

    hicCreateThresholdFile.main(
        "--range 0 0 -tv 1 -r 10000 -o {}".format(outfile.name).split())
    assert read_lines(outfile.name)[1:] == ['# Standard threshold 1.0', '0\t1.0']

    hicCreateThresholdFile.main(
        "--range 0 0 -tv 0.00001 -r 10000 -o {}".format(outfile.name).split())
    assert read_lines(outfile.name)[1:] == ['# Standard threshold 1e-05', '0\t1e-05']

    os.unlink(outfile.name)
