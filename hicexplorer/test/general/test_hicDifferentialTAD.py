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
MID_MEMORY = 7
HIGH_MEMORY = 200

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


def test_cool_all_single_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_all_multi_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_cool_all_multi_core_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_all_multi_core_left_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_all_multi_core_right_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_one_all_multi_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_cool_one_multi_core_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_one_multi_core_left_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_cool_one_multi_core_right_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_all_single_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_all_multi_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_h5_all_multi_core_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_all_multi_core_left_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_all_multi_core_right_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_one_all_multi_core():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

# intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'
def test_h5_one_multi_core_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_one_multi_core_left_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)

def test_h5_one_multi_core_right_intra_TAD():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5 -t 4 -tpc 4".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)