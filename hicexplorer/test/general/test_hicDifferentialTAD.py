import os.path
from tempfile import NamedTemporaryFile, mkdtemp
from psutil import virtual_memory
import pytest
import logging
log = logging.getLogger(__name__)
import numpy as np
import pandas as pd
from pybedtools import BedTool

from hicexplorer import hicDifferentialTAD
from hicexplorer.test.test_compute_function import compute

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 4
HIGH_MEMORY = 120

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/hicDifferentialTAD/")


def are_files_equal(file1, file2, delta=1, skip=0, eps=0.1):
    """Line by line comparison after `skip` header lines.

    The earlier version iterated over zip(file1, file2), which stops at the
    shorter file, so a produced file that lost or gained rows at its end still
    compared equal. Both files must now have the same number of lines, and at
    most `delta` of them may differ.
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
    if mismatches <= delta:
        return True
    print('mismatches: {}'.format(mismatches))
    return False


def all_tads_present(pOriginTADs, pAccepted, pRejected):

    original_tads = pd.read_csv(pOriginTADs, sep='\t', header=None)[[0, 1, 2]]
    accepted_tads = pd.read_csv(pAccepted, sep='\t', header=None, skiprows=4)[[0, 1, 2]]
    rejected_tads = pd.read_csv(pRejected, sep='\t', header=None, skiprows=4)[[0, 1, 2]]

    original_tads_bedtool = BedTool.from_dataframe(original_tads)
    accepted_tads_bedtool = BedTool.from_dataframe(accepted_tads)
    rejected_tads_bedtool = BedTool.from_dataframe(rejected_tads)

    x = original_tads_bedtool.intersect(accepted_tads_bedtool, c=True).to_dataframe()
    y = original_tads_bedtool.intersect(rejected_tads_bedtool, c=True).to_dataframe()

    mask_x = x['name'] >= 1
    mask_y = y['name'] >= 1

    # print(x)
    # print(mask_x.sum())
    # print(mask_y.sum())

    return mask_x.sum() + mask_y.sum()
    # selection = (mask_x) & (mask_y)


def test_cool_all_single_core():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'all', 'all'
    ).split()
    # compute(hicDifferentialTAD.main, args, 5)
    hicDifferentialTAD.main(args)
    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_rejected.diff_tad', delta=0, skip=4)


def test_cool_all_multi_core():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_rejected.diff_tad', delta=0, skip=4)


# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_cool_all_single_core_intra_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'intra-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_intra-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_intra-TAD_rejected.diff_tad', delta=0, skip=4)


def test_cool_all_multi_core_intra_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'intra-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_intra-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_intra-TAD_rejected.diff_tad', delta=0, skip=4)


# 'left-inter-TAD'
def test_cool_all_single_core_left_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'left-inter-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_left_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_left_inter-TAD_rejected.diff_tad', delta=0, skip=4)


def test_cool_all_multi_core_left_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        11, outfile_pref.name, 'left-inter-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_left_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_left_inter-TAD_rejected.diff_tad', delta=0, skip=4)


# 'left-inter-TAD'
def test_cool_all_single_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'right-inter-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_right_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_right_inter-TAD_rejected.diff_tad', delta=0, skip=4)


def test_cool_all_multi_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'right-inter-TAD', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_right_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_right_inter-TAD_rejected.diff_tad', delta=0, skip=4)


# 'left-inter-TAD'
def test_cool_one_single_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'right-inter-TAD', 'one'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_right_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_right_inter-TAD_rejected.diff_tad', delta=0, skip=4)


def test_cool_one_multi_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'right-inter-TAD', 'one'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_right_inter-TAD_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_right_inter-TAD_rejected.diff_tad', delta=0, skip=4)


def test_h5_all_single_core_one():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'all', 'one'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_h5_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_h5_rejected.diff_tad', delta=0, skip=4)


def test_h5_all_multi_core_all():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_h5_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_h5_rejected.diff_tad', delta=0, skip=4)


@pytest.mark.xfail(reason='Access of a chromosome which is not in the matrix.')
def test_h5_all_multi_core_all_multichr_chromosome_not_in_matrix():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5",
        ROOT + "untreated_R1_domains_chr1_chr2.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)


@pytest.mark.xfail(reason='Access of a chromosome which is not in the matrix.')
def test_cool_all_multi_core_all_multichr_chromosome_not_in_matrix():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains_chr1_chr2.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)


def test_cool_all_multi_core_all_multichr_chromosome():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1_chr2.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool",
        ROOT + "untreated_R1_domains_chr1_chr2.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 5)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains_chr1_chr2.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains_chr1_chr2.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')

    # delta was 2. The output has been measured to match the master exactly at
    # --threads 1, 4, 11 and 16 (2026-09-13), so the slack only hid changes.
    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'multichromosome_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'multichromosome_rejected.diff_tad', delta=0, skip=4)

    # mode_all_reject_all__t4_cool_multi_chromosomes_rejected.diff_tad


def test_cool_all_one_core_all_multichr_chromosome():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1_chr2.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool",
        ROOT + "untreated_R1_domains_chr1_chr2.bed",
        1, outfile_pref.name, 'all', 'all'
    ).split()
    compute(hicDifferentialTAD.main, args, 7)

    # test only on the line numbers
    with open(ROOT + "untreated_R1_domains_chr1_chr2.bed", 'r') as file:
        number_of_tads = len(file.readlines())

    with open(outfile_pref.name + '_accepted.diff_tad', 'r') as file:
        number_output_tads = len(file.readlines())
    with open(outfile_pref.name + '_rejected.diff_tad', 'r') as file:
        number_output_tads += len(file.readlines())

    number_output_tads -= 8
    assert number_of_tads == number_output_tads

    # test on the intersection to exclude the case of duplicated lines
    assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains_chr1_chr2.bed", outfile_pref.name + '_accepted.diff_tad', outfile_pref.name + '_rejected.diff_tad')
    # assert number_of_tads == all_tads_present(ROOT + "untreated_R1_domains_chr1_chr2.bed", ROOT + 'mode_all_reject_all__t4_cool_multi_chromosomes_accepted.diff_tad', ROOT + 'mode_all_reject_all__t4_cool_multi_chromosomes_rejected.diff_tad')

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'multichromosome_accepted.diff_tad', delta=0, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'multichromosome_rejected.diff_tad', delta=0, skip=4)


# --------------------------------------------------------------------------
# Characterization tests added for the v4 C++ port (2026-09-13). They pin the
# current Python behaviour, including the defects, on the real matrices and on
# subsets of the real TAD calls. Nothing here asserts what the tool should do.


def _run(pTarget, pControl, pDomains, pThreads, pMode, pModeReject, pPValue=None):
    outdir = mkdtemp(prefix='differential_tad_')
    prefix = os.path.join(outdir, 'out')
    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        pTarget, pControl, pDomains, pThreads, prefix, pMode, pModeReject).split()
    if pPValue is not None:
        args += ['--pValue', str(pPValue)]
    compute(hicDifferentialTAD.main, args, 5)
    return prefix


def _data_rows(pFile):
    with open(pFile) as handle:
        return [line.rstrip('\n').split('\t') for line in handle if not line.startswith('#')]


def _assert_matches_master(pPrefix, pMaster):
    for kind in ['accepted', 'rejected']:
        assert are_files_equal(pPrefix + '_' + kind + '.diff_tad',
                               ROOT + pMaster + '_' + kind + '.diff_tad', delta=0, skip=4)


COOL_TARGET = ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool"
COOL_CONTROL = ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool"
H5_TARGET = ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5"
H5_CONTROL = ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5"
COOL2_TARGET = ROOT + "GSM2644945_Untreated-R1.100000_chr1_chr2.cool"
COOL2_CONTROL = ROOT + "GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool"


def test_cool_pvalue_001_all_one():
    # --pValue was never varied from its default.
    prefix = _run(COOL_TARGET, COOL_CONTROL, ROOT + "untreated_R1_domains.bed", 2, 'all', 'one', 0.01)
    _assert_matches_master(prefix, 'pvalue_0.01_mode_all_one')
    with open(prefix + '_rejected.diff_tad') as handle:
        header = handle.readlines()[:4]
    assert header[2] == "# Rejected regions with Wilcoxon rank-sum test to p-value: 0.01  with used mode: all and modeReject: one \n"
    assert header[3].startswith('# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD')

    # The split is a pure function of the printed p-values: rejected when any
    # of the three is <= 0.01, a NaN never rejects.
    for kind, expected in [('accepted', False), ('rejected', True)]:
        for row in _data_rows(prefix + '_' + kind + '.diff_tad'):
            p_values = np.array(row[6:9], dtype=float)
            assert bool(np.any(p_values <= 0.01)) == expected


def test_cool_left_inter_TAD_one():
    # left-inter-TAD was only ever run with --modeReject all.
    prefix = _run(COOL_TARGET, COOL_CONTROL, ROOT + "untreated_R1_domains.bed", 1, 'left-inter-TAD', 'one')
    _assert_matches_master(prefix, 'mode_one_left_inter-TAD')


def test_h5_intra_TAD_all():
    # The h5 path was only ever run with --mode all.
    prefix = _run(H5_TARGET, H5_CONTROL, ROOT + "untreated_R1_domains.bed", 1, 'intra-TAD', 'all')
    _assert_matches_master(prefix, 'mode_intra-TAD_h5')


def test_h5_left_inter_TAD_one():
    prefix = _run(H5_TARGET, H5_CONTROL, ROOT + "untreated_R1_domains.bed", 3, 'left-inter-TAD', 'one')
    _assert_matches_master(prefix, 'mode_one_left_inter-TAD_h5')


def test_h5_last_tad_has_no_left_inter_tad_test():
    # For the last TAD of a chromosome right_boundary_index is not assigned in
    # its own iteration, so the value left over from the previous TAD is used.
    # On the h5 path that is an absolute bin index equal to the last TAD's own
    # first bin, so the left inter-TAD block has zero columns and ranksums of
    # two empty samples is NaN.
    prefix = _run(H5_TARGET, H5_CONTROL, ROOT + "untreated_R1_domains.bed", 1, 'all', 'all')
    rows = _data_rows(prefix + '_accepted.diff_tad') + _data_rows(prefix + '_rejected.diff_tad')
    last = max(rows, key=lambda row: int(row[1]))
    assert last[6] == 'nan' and last[9] == 'nan'
    assert last[7] == 'nan' and last[10] == 'nan'


def test_cool_last_tad_left_inter_tad_uses_stale_index():
    # On the cool path the left-over right boundary is relative to the previous
    # TAD's region, not to the current one, so the left inter-TAD test of the
    # last TAD runs on a block that is not the one the name suggests. It is not
    # NaN, which is what distinguishes it from the h5 path above.
    prefix = _run(COOL_TARGET, COOL_CONTROL, ROOT + "untreated_R1_domains.bed", 1, 'all', 'all')
    rows = _data_rows(prefix + '_accepted.diff_tad') + _data_rows(prefix + '_rejected.diff_tad')
    last = max(rows, key=lambda row: int(row[1]))
    assert last[1:3] == ['193700000', '195200000']
    assert last[6] == '0.9185395091955513'
    assert last[7] == 'nan'


def _row_count(pPrefix):
    return len(_data_rows(pPrefix + '_accepted.diff_tad')) + len(_data_rows(pPrefix + '_rejected.diff_tad'))


def test_threads_small_chromosome_first_drops_last_tad():
    # Defect pinned, not endorsed. A chromosome with fewer TADs than --threads
    # sets args.threads = 1 for good. The next chromosome then computes
    # domainsPerThread with 1 thread but restores the thread count, so the
    # first process gets the whole list and skips its last element, and no
    # other process receives it: the last TAD of that chromosome is lost.
    domains = ROOT + "domains_chr1_3_chr2_10.bed"
    assert _row_count(_run(COOL2_TARGET, COOL2_CONTROL, domains, 1, 'all', 'one')) == 13
    prefix = _run(COOL2_TARGET, COOL2_CONTROL, domains, 4, 'all', 'one')
    assert _row_count(prefix) == 12
    rows = _data_rows(prefix + '_accepted.diff_tad') + _data_rows(prefix + '_rejected.diff_tad')
    assert ['chr2', '17100000', '18100000'] not in [row[0:3] for row in rows]


def test_threads_small_chromosome_last_duplicates_rows():
    # Defect pinned, not endorsed. The per process result slots are allocated
    # once for the original thread count. A chromosome processed with a single
    # process only overwrites slot 0, so slots 1..3 still hold the previous
    # chromosome's results and they are written a second time.
    domains = ROOT + "domains_chr1_10_chr2_3.bed"
    assert _row_count(_run(COOL2_TARGET, COOL2_CONTROL, domains, 1, 'all', 'one')) == 13
    prefix = _run(COOL2_TARGET, COOL2_CONTROL, domains, 4, 'all', 'one')
    assert _row_count(prefix) == 21
    rows = _data_rows(prefix + '_accepted.diff_tad') + _data_rows(prefix + '_rejected.diff_tad')
    keys = [tuple(row[0:3]) for row in rows]
    assert len(set(keys)) == 13


def test_threads_one_tad_per_thread_loses_last_left_test():
    # Defect pinned, not endorsed. With exactly one TAD per process the last
    # process sees [TAD n-2, TAD n-1] and the last TAD sits at local index 1,
    # where the `i - 1 > 0` guard fails, so its left inter-TAD test is skipped.
    domains = ROOT + "domains_chr1_4.bed"
    single = _run(COOL_TARGET, COOL_CONTROL, domains, 1, 'all', 'one')
    multi = _run(COOL_TARGET, COOL_CONTROL, domains, 4, 'all', 'one')

    def last_row(pPrefix):
        rows = _data_rows(pPrefix + '_accepted.diff_tad') + _data_rows(pPrefix + '_rejected.diff_tad')
        return max(rows, key=lambda row: int(row[1]))
    assert last_row(single)[6] == '0.8706579826550646'
    assert last_row(multi)[6] == 'nan'
