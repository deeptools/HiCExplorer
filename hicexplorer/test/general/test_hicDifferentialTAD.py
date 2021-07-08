import os.path
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import pytest
import logging
log = logging.getLogger(__name__)
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

    mismatches = 0
    # if delta:
    #     mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            if i < skip:
                continue
            if x.startswith('File'):
                continue
            if x != y:
                # if delta:
                mismatches += 1
                if mismatches > delta:
                    break
                # else:
                #     break
    if mismatches <= delta:
        return True
    else:
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
