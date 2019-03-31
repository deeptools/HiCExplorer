import os.path
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import logging
log = logging.getLogger(__name__)

from hicexplorer import hicDetectLoops

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


def test_main_h5():
    outfile_loop_h5 = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} -pit 1 -p 0.5 -pp 0.5".format(
        ROOT + "small_test_matrix.h5", outfile_loop_h5.name).split()
    hicDetectLoops.main(args)



def test_main_cool_chromosomes():
    outfile_loop_cool = NamedTemporaryFile(suffix='.bedgraph', delete=True)

    args = "--matrix {} -o {} --maxLoopDistance 30000000 -w 5 -pw 2 -p 0.5 -pp 0.55 --chromosomes 1 2".format(
        ROOT + "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool", outfile_loop_cool.name).split()
    hicDetectLoops.main(args)
    assert are_files_equal(
        ROOT + "hicDetectLoops/loops.bedgraph", outfile_loop_cool.name, delta=2)
