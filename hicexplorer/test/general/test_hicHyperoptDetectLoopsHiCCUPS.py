import os
from tempfile import NamedTemporaryFile
from psutil import virtual_memory
import subprocess
import pytest
import logging
log = logging.getLogger(__name__)

from hicexplorer import hicHyperoptDetectLoopsHiCCUPS

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True

ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")

NVCC = os.system('which nvcc')


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


@pytest.mark.skipif(NVCC != 0,
                    reason="No cuda installed, skipping test case")
def test_main():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)

    bashCommand = 'curl https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.21.01.jar -o juicer.jar'
    subprocess.check_output(['bash', '-c', bashCommand])
    args = "-j {} -m {} -p {} -ml {} -r {} --runs {} -o {}".format(
        'juicer.jar', ROOT + 'hicHyperoptDetectLoopsHiCCUPS/SRR1791297_30.hic', ROOT + 'hicHyperoptDetectLoopsHiCCUPS/ctcf_sorted.bed', 7, 10000, 2, outfile.name).split()
    hicHyperoptDetectLoopsHiCCUPS.main(args)
    are_files_equal(outfile.name, ROOT + 'hicHyperoptDetectLoopsHiCCUPS/hyperoptHiCCUPS_result.txt', delta=2)

    os.remove('juicer.jar')
# hicHyperoptDetectLoopsHiCCUPS -j juicer_tools_1.21.01.jar -m 30qc/GSE63525_HMEC_combined_30.hic -p 30qc/hmec/ctcf_sorted.bed -ml 3210  -r 10000 --runs 1 --cpu
