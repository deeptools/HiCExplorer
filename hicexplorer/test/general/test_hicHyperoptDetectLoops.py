import os.path
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import logging
log = logging.getLogger(__name__)

from hicexplorer import hicHyperoptDetectLoops

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
    # if delta:
    mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                mismatches += 1
                if mismatches > delta:
                    return False
    return equal


def test_main():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile.close()
    args = "--matrix {} -p {} -ml {} -r {} --runs {} -o {}".format(
        ROOT + "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool",
        ROOT + 'hicHyperoptDectedLoops/ctcf_sorted.bed', 3210, 10000, 2, outfile.name).split()
    hicHyperoptDetectLoops.main(args)
    assert are_files_equal(outfile.name, ROOT + 'hicHyperoptDetectLoops/hyperopt_result.txt', delta=2)


def test_main_add():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile.close()
    args = "--matrix {} -p {} -ml {} -r {} --runs {} -o {} -cl {}".format(
        ROOT + "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool",
        ROOT + 'hicHyperoptDectedLoops/ctcf_sorted.bed', 3210, 10000, 2, outfile.name, 'add').split()
    hicHyperoptDetectLoops.main(args)
    assert are_files_equal(outfile.name, ROOT + 'hicHyperoptDetectLoops/hyperopt_result.txt', delta=2)


def test_main_remove():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile.close()
    args = "--matrix {} -p {} -ml {} -r {} --runs {} -o {} -cl {}".format(
        ROOT + "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool",
        ROOT + 'hicHyperoptDectedLoops/ctcf_sorted.bed', 3210, 10000, 2, outfile.name, 'remove').split()
    hicHyperoptDetectLoops.main(args)
    assert are_files_equal(outfile.name, ROOT + 'hicHyperoptDetectLoops/hyperopt_result.txt', delta=2)
