import os
import sys
from tempfile import NamedTemporaryFile
from tempfile import mkdtemp
from psutil import virtual_memory
import subprocess
import pytest
import logging
log = logging.getLogger(__name__)

import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from matplotlib.testing.exceptions import ImageComparisonFailure
from hicexplorer import hicInterIntraTAD, hicMergeDomains
from hicexplorer.test.test_compute_function import compute

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def are_files_equal(file1, file2, delta=None, skip=0):
    """Line by line comparison after `skip` header lines.

    The earlier version iterated over zip(file1, file2), which stops at the
    shorter file, and test_main called it without assert, so its result was
    discarded. Both files must now have the same number of lines, and at most
    `delta` lines may differ.
    """
    with open(file1) as textfile1:
        lines1 = textfile1.readlines()[skip:]
    with open(file2) as textfile2:
        lines2 = textfile2.readlines()[skip:]
    if len(lines1) != len(lines2):
        print('line count differs: {} vs {}'.format(len(lines1), len(lines2)))
        return False
    mismatches = 0
    for x, y in zip(lines1, lines2):
        if x.startswith('File'):
            continue
        if x != y:
            mismatches += 1
    return mismatches <= (delta or 0)


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_main():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_plot = NamedTemporaryFile(suffix='.png', delete=True)

    args = "-m {} --tadDomains {} --threads {} --outFileNameRatioPlot {} -o {}".format(
        ROOT + 'hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool',
        ROOT + 'hicInterIntraTAD/untreated_R1_domains_chr1_chr2.bed', 5, outfile_plot.name, outfile.name).split()
    compute(hicInterIntraTAD.main, args, 5)
    # The first line carries the version string, which is not part of the data.
    assert are_files_equal(outfile.name, ROOT + 'hicInterIntraTAD/output_test.txt', delta=0, skip=1)
    res = compare_images(ROOT + 'hicInterIntraTAD/ratio.png', outfile_plot.name, tol=40)
    assert res is None, res


# --------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port (2026-09-13). They pin the
# current Python behaviour, including the defects, on the real matrices.


def _run(pMatrix, pDomains, pThreads):
    outdir = mkdtemp(prefix='interintra_tad_')
    outfile = os.path.join(outdir, 'out.txt')
    args = "-m {} --tadDomains {} --threads {} --outFileNameRatioPlot {} -o {}".format(
        pMatrix, pDomains, pThreads, os.path.join(outdir, 'ratio.png'), outfile).split()
    compute(hicInterIntraTAD.main, args, 5)
    return outfile


def _data_rows(pFile):
    with open(pFile) as handle:
        return [line.rstrip('\n').split('\t') for line in handle if not line.startswith('#')]


def test_text_output_two_chromosomes():
    # The text half of test_main without the xfail, so that a wrong table
    # cannot hide behind the expected image failure.
    outfile = _run(ROOT + 'hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool',
                   ROOT + 'hicInterIntraTAD/untreated_R1_domains_chr1_chr2.bed', 1)
    assert are_files_equal(outfile, ROOT + 'hicInterIntraTAD/output_test.txt', delta=0, skip=1)
    with open(outfile) as handle:
        header = handle.readlines()[:2]
    assert header[0].startswith("# Created with HiCExplorer's hicInterIntraTAD version ")
    assert header[1] == ('# Chromosome\tstart\tend\tname\tscore\tstrand\tinter_left_sum\tinter_right_sum\t'
                         'inter_left_density\tinter_right_density\tinter_left_number_of_contacts\t'
                         'inter_right_number_of_contacts\tinter_left_number_of_contacts_nnz\t'
                         'inter_right_number_of_contacts_nnz\tintra_sum\tintra_number_of_contacts\t'
                         'intra_number_of_contacts_nnz\tintra_density\tinter_left_intra_ratio\t'
                         'inter_right_intra_ratio\tinter_left_inter_right_intra_ratio\n')
    assert os.path.getsize(os.path.join(os.path.dirname(outfile), 'ratio.png')) > 0


def test_text_output_untreated_chr1():
    # A second matrix and a single chromosome, so a different TAD list and a
    # different number of processes than test_main.
    outfile = _run(ROOT + 'hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.cool',
                   ROOT + 'hicDifferentialTAD/untreated_R1_domains.bed', 3)
    assert are_files_equal(outfile, ROOT + 'hicInterIntraTAD/output_untreated_chr1.txt', delta=0, skip=1)
    rows = _data_rows(outfile)
    assert len(rows) == 140
    # Every row satisfies the identities between its own columns.
    for row in rows:
        assert float(row[17]) == int(row[16]) / int(row[15])
        if row[8] != '0':
            assert float(row[8]) == int(row[12]) / int(row[10])
        if row[9] != '0':
            assert float(row[9]) == int(row[13]) / int(row[11])


def test_first_and_last_tad_columns():
    # The first TAD of a chromosome has no left neighbour, so its left columns
    # are the integer 0 and print without a decimal point. The last TAD's left
    # block uses the right boundary index left over from the previous TAD.
    outfile = _run(ROOT + 'hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.cool',
                   ROOT + 'hicDifferentialTAD/untreated_R1_domains.bed', 1)
    rows = _data_rows(outfile)
    first = rows[0]
    assert first[6] == '0' and first[8] == '0' and first[10] == '0' and first[12] == '0'
    last = rows[-1]
    assert last[7] == '0' and last[9] == '0' and last[11] == '0' and last[13] == '0'
    assert last[6] != '0'


def test_h5_input_raises_zero_division():
    # Defect pinned, not endorsed. On the h5 path the left-over right boundary
    # of the last TAD is an absolute index equal to that TAD's first bin, so
    # its left inter-TAD block has zero columns and the density divides by
    # zero. The worker reports the failure and main exits with status 1 before
    # anything is written.
    outdir = mkdtemp(prefix='interintra_tad_')
    outfile = os.path.join(outdir, 'out.txt')
    args = "-m {} --tadDomains {} --threads {} --outFileNameRatioPlot {} -o {}".format(
        ROOT + 'hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5',
        ROOT + 'hicDifferentialTAD/untreated_R1_domains.bed', 1,
        os.path.join(outdir, 'ratio.png'), outfile).split()
    with pytest.raises(SystemExit) as error:
        hicInterIntraTAD.main(args)
    assert error.value.code == 1
    assert not os.path.exists(outfile)


def test_threads_small_chromosome_last_duplicates_rows():
    # Defect pinned, not endorsed; the same partitioning code as
    # hicDifferentialTAD. Slots 1..3 still hold chr1's results when chr2 is
    # processed by a single process, so eight chr1 rows are written twice.
    domains = ROOT + 'hicDifferentialTAD/domains_chr1_10_chr2_3.bed'
    matrix = ROOT + 'hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool'
    assert len(_data_rows(_run(matrix, domains, 1))) == 13
    rows = _data_rows(_run(matrix, domains, 4))
    assert len(rows) == 21
    assert len(set(tuple(row[0:3]) for row in rows)) == 13


def test_threads_small_chromosome_first_drops_last_tad():
    domains = ROOT + 'hicDifferentialTAD/domains_chr1_3_chr2_10.bed'
    matrix = ROOT + 'hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool'
    assert len(_data_rows(_run(matrix, domains, 1))) == 13
    assert len(_data_rows(_run(matrix, domains, 4))) == 12
