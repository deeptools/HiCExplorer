from hicexplorer import hicFindRestSite
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/hicFindRestSite/")


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


def assert_files_identical(pExpected, pObserved):
    """Byte equality, and the same number of lines.

    `are_files_equal` above walks the two files with zip, which stops at the
    shorter of the two, and then tolerates one differing line. A port that
    wrote the first half of the sites, or that dropped one site, passed it. The
    output of this tool is a few dozen lines of pure text with no floats in it,
    so there is nothing to be tolerant about.
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


def test_fasta_gz():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} -o {}".format(ROOT + 'dm3_chrM.fasta.gz', 'AAGCTT', outfile.name).split()
    hicFindRestSite.main(args)

    assert_files_identical(ROOT + "hindIII_chrM.bed", outfile.name)


def test_fasta():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} -o {}".format(ROOT + 'dm3_chrM.fasta', 'AAGCTT', outfile.name).split()
    hicFindRestSite.main(args)

    assert_files_identical(ROOT + "hindIII_chrM.bed", outfile.name)


def test_fasta_gz_two_patterns():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} {} -o {}".format(ROOT + 'dm3_chrM.fasta.gz', 'AAGCTT', 'GATC', outfile.name).split()
    hicFindRestSite.main(args)

    assert_files_identical(ROOT + "hindIII_DpnII.bed", outfile.name)


def test_fasta_two_patterns():
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()
    args = "-f {}  -p {} {} -o {}".format(ROOT + 'dm3_chrM.fasta', 'AAGCTT', 'GATC', outfile.name).split()
    hicFindRestSite.main(args)

    assert_files_identical(ROOT + "hindIII_DpnII.bed", outfile.name)


def test_gz_and_plain_agree():
    """The two readers must produce the same file, not merely a similar one.

    mimetypes.guess_type decides which one runs, off the file name alone
    (hicFindRestSite.py:108), so this is also the test that a port keeps that
    decision rule rather than sniffing the magic bytes.
    """
    from_gz = NamedTemporaryFile(suffix='.bed', delete=False)
    from_gz.close()
    from_plain = NamedTemporaryFile(suffix='.bed', delete=False)
    from_plain.close()

    hicFindRestSite.main("-f {} -p AAGCTT GATC -o {}".format(
        ROOT + 'dm3_chrM.fasta.gz', from_gz.name).split())
    hicFindRestSite.main("-f {} -p AAGCTT GATC -o {}".format(
        ROOT + 'dm3_chrM.fasta', from_plain.name).split())

    assert_files_identical(from_gz.name, from_plain.name)


def test_non_palindromic_regex_pattern():
    """A pattern that is a regexp and is not its own reverse complement.

    Neither the reverse-strand branch (hicFindRestSite.py:118-122) nor the
    regexp syntax was exercised by this file: AAGCTT and GATC are both
    palindromic literals, so `rev_compl != pattern` was always false and no
    line with a minus strand was ever produced. CG.AG is the example from the
    tool's own doctest; its reverse complement is CT.CG, and on dm3 chrM the
    two strands together give nine sites, four of them on the minus strand.
    The expected file is written out in full because it is nine lines, and
    because embedding it pins the coordinate convention (zero based,
    half open, end = start + len(pattern)) as well as the strand column.
    """
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()

    hicFindRestSite.main("-f {} -p CG.AG -o {}".format(
        ROOT + 'dm3_chrM.fasta', outfile.name).split())

    assert read_lines(outfile.name) == [
        'chrM\t2959\t2964\t.\t0\t-',
        'chrM\t4803\t4808\t.\t0\t+',
        'chrM\t6862\t6867\t.\t0\t+',
        'chrM\t6969\t6974\t.\t0\t+',
        'chrM\t7707\t7712\t.\t0\t+',
        'chrM\t8518\t8523\t.\t0\t+',
        'chrM\t8558\t8563\t.\t0\t-',
        'chrM\t13314\t13319\t.\t0\t-',
        'chrM\t14809\t14814\t.\t0\t-',
    ]

    os.unlink(outfile.name)


def test_two_patterns_at_the_same_start_keep_the_first_pattern():
    """`sort -u` de-duplicates on the sort key, not on the whole line.

    hicFindRestSite.py:129 sorts with `-k1,1 -k2,2n -u`, so two sites that
    share a chromosome and a start collapse into one even when their end and
    their strand differ, and the survivor is the first of them in the temporary
    file. The temporary file is written pattern by pattern, so the survivor is
    the site of the pattern that was given first on the command line.

    AAGC occurs at every one of the four AAGCTT sites of dm3 chrM, since AAGCTT
    begins with AAGC, which makes the collision easy to see: both orders give
    74 sites, and the four shared starts carry the six base site with AAGCTT
    first and the four base one with AAGC first.
    """
    hindiii_first = NamedTemporaryFile(suffix='.bed', delete=False)
    hindiii_first.close()
    aagc_first = NamedTemporaryFile(suffix='.bed', delete=False)
    aagc_first.close()

    hicFindRestSite.main("-f {} -p AAGCTT AAGC -o {}".format(
        ROOT + 'dm3_chrM.fasta', hindiii_first.name).split())
    hicFindRestSite.main("-f {} -p AAGC AAGCTT -o {}".format(
        ROOT + 'dm3_chrM.fasta', aagc_first.name).split())

    first = read_lines(hindiii_first.name)
    second = read_lines(aagc_first.name)
    assert len(first) == 74
    assert len(second) == 74

    shared_starts = ['402', '5267', '5689', '14117']
    for start in shared_starts:
        line_first = [line for line in first if line.split('\t')[1] == start]
        line_second = [line for line in second if line.split('\t')[1] == start]
        assert len(line_first) == 1
        assert len(line_second) == 1
        assert int(line_first[0].split('\t')[2]) - int(start) == 6
        assert int(line_second[0].split('\t')[2]) - int(start) == 4

    os.unlink(hindiii_first.name)
    os.unlink(aagc_first.name)


def test_output_is_sorted_by_chromosome_then_numeric_start():
    """The sort is `-k1,1 -k2,2n`: bytes on the name, numeric on the start.

    A lexicographic sort of the start column would put 10000 before 402, so
    this fails on any port that sorts the whole line as text.
    """
    outfile = NamedTemporaryFile(suffix='.bed', delete=False)
    outfile.close()

    hicFindRestSite.main("-f {} -p AAGCTT GATC -o {}".format(
        ROOT + 'dm3_chrM.fasta', outfile.name).split())

    keys = []
    for line in read_lines(outfile.name):
        fields = line.split('\t')
        keys.append((fields[0], int(fields[1])))
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys)

    os.unlink(outfile.name)
