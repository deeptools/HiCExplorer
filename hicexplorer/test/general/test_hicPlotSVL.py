import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPlotSVL
from hicmatrix import HiCMatrix as hm

from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
import numpy as np
# import pyBigWig
from matplotlib.testing.compare import compare_images

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

import logging
log = logging.getLogger(__name__)


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


def test_plotSVL():
    plot = NamedTemporaryFile(suffix='.png', delete=False)
    outputFileName = NamedTemporaryFile(suffix='.txt', delete=False)
    outputFileNameData = NamedTemporaryFile(suffix='.txt', delete=False)

    plot.close()
    outputFileName.close()
    outputFileNameData.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    matrix2 = ROOT + "small_test_matrix_50kb_res.h5"

    args = "--matrices {} {} --plotFileName {} --outFileName {} --outFileNameData {} --dpi 300"\
        .format(matrix, matrix2, plot.name, outputFileName.name, outputFileNameData.name).split()
    hicPlotSVL.main(args)

    assert are_files_equal(ROOT + 'hicPlotSVL/data.txt', outputFileNameData.name, delta=2)
    assert are_files_equal(ROOT + 'hicPlotSVL/p_values.txt', outputFileName.name, delta=2)
    res = compare_images(ROOT + 'hicPlotSVL/plot.png', plot.name, tol=50)
    assert res is None, res
    os.unlink(plot.name)
    os.unlink(outputFileName.name)
    os.unlink(outputFileNameData.name)
