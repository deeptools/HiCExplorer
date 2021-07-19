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


def are_files_equal(file1, file2, delta=None):
    equal = True
    mismatches = 0
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                mismatches += 1
                if mismatches > delta:
                    equal = False
                    break
    return equal


@pytest.mark.xfail(raises=ImageComparisonFailure, reason='Matplotlib plots for reasons a different image size.')
def test_main():
    outfile = NamedTemporaryFile(suffix='.txt', delete=True)
    outfile_plot = NamedTemporaryFile(suffix='.png', delete=True)

    args = "-m {} --tadDomains {} --threads {} --outFileNameRatioPlot {} -o {}".format(
        ROOT + 'hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool',
        ROOT + 'hicInterIntraTAD/untreated_R1_domains_chr1_chr2.bed', 5, outfile_plot.name, outfile.name).split()
    compute(hicInterIntraTAD.main, args, 5)
    are_files_equal(outfile.name, ROOT + 'hicInterIntraTAD/output_test.txt', delta=2)
    res = compare_images(ROOT + 'hicInterIntraTAD/ratio.png', outfile_plot.name, tol=40)
    assert res is None, res
