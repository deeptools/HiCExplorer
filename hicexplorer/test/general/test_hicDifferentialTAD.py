import os.path
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import logging
log = logging.getLogger(__name__)

from hicexplorer import hicDifferentialTAD

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
    if delta:
        mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for i, (x, y) in enumerate(zip(textfile1, textfile2)):
            if i < skip:
                continue
            if x.startswith('File'):
                continue
            if x != y:
                if delta:
                    mismatches += 1
                    if mismatches > delta:
                        break
                else:
                    break
    if mismatches < delta:
        return True
    else:
        log.debug('mismatches: {}'.format(mismatches))
        return False


def test_cool_all_single_core():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'all', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_all_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_all_rejected.diff_tad', delta=5, skip=4)


def test_cool_all_multi_core():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_all_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_all_rejected.diff_tad', delta=5, skip=4)


# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_cool_all_single_core_intra_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'intra-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_intra-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_intra-TAD_rejected.diff_tad', delta=5, skip=4)


def test_cool_all_multi_core_intra_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'intra-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_intra-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_intra-TAD_rejected.diff_tad', delta=5, skip=4)


# 'left-inter-TAD'
def test_cool_all_single_core_left_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'left-inter-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_left_inter-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_left_inter-TAD_rejected.diff_tad', delta=5, skip=4)


def test_cool_all_multi_core_left_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'left-inter-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_left_inter-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_left_inter-TAD_rejected.diff_tad', delta=5, skip=4)


# 'left-inter-TAD'
def test_cool_all_single_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'right-inter-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_right_inter-TAD_t1_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_right_inter-TAD_t1_rejected.diff_tad', delta=5, skip=4)


def test_cool_all_multi_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'right-inter-TAD', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_right_inter-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_right_inter-TAD_rejected.diff_tad', delta=5, skip=4)


# 'left-inter-TAD'
def test_cool_one_single_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'right-inter-TAD', 'one'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_reject_right_inter-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_reject_right_inter-TAD_rejected.diff_tad', delta=5, skip=4)


def test_cool_one_multi_core_right_inter_TAD():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.cool",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.cool",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'right-inter-TAD', 'one'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_reject_right_inter-TAD_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_reject_right_inter-TAD_rejected.diff_tad', delta=5, skip=4)


def test_h5_all_single_core_one():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5",
        ROOT + "untreated_R1_domains.bed",
        1, outfile_pref.name, 'all', 'one'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_one_reject_all_h5_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_one_reject_all_h5_rejected.diff_tad', delta=5, skip=4)


def test_h5_all_multi_core_all():
    outfile_pref = NamedTemporaryFile(prefix='differential_tad', delete=True)

    args = "--targetMatrix {} --controlMatrix {} --tadDomains {} -t {} -o {} -m {} -mr {}".format(
        ROOT + "GSM2644945_Untreated-R1.100000_chr1.h5",
        ROOT + "GSM2644947_Auxin2days-R1.100000_chr1.h5",
        ROOT + "untreated_R1_domains.bed",
        4, outfile_pref.name, 'all', 'all'
    ).split()
    hicDifferentialTAD.main(args)

    assert are_files_equal(outfile_pref.name + '_accepted.diff_tad', ROOT + 'mode_all_reject_all__t4_h5_accepted.diff_tad', delta=5, skip=4)
    assert are_files_equal(outfile_pref.name + '_rejected.diff_tad', ROOT + 'mode_all_reject_all__t4_h5_rejected.diff_tad', delta=5, skip=4)
