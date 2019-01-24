import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicConvertFormat
import pytest
from hicmatrix import HiCMatrix as hm
import numpy.testing as nt
import numpy as np

REMOVE_OUTPUT = True
# DIFF = 60

DELTA_DECIMAL = 0

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicConvertFormat")
original_matrix_h5 = ROOT + "/small_test_matrix.h5"
original_matrix_cool = ROOT + "/small_test_matrix.cool"


@pytest.mark.parametrize("matrices", [original_matrix_h5, original_matrix_cool])  # , original_matrix_cool, original_matrix_hic])  # required
@pytest.mark.parametrize("outputFormat", ['cool', 'h5', 'ginteractions', 'mcool'])
@pytest.mark.parametrize("resolutions", [''])  # TODO: check for possible resolutions
def test_trivial_run(
    matrices,
    outputFormat,
    resolutions,
):
    """
        Test for all commandline arguments.
        Options for cool input format are testet seperately.
    """
    from pathlib import Path
    # get suffix of input matrix without the dot
    inputFormat = Path(matrices).suffix[1:]
    # create file corresponding to output format
    outFileName = NamedTemporaryFile(suffix=".{}".format(outputFormat), delete=True)

    args = "--matrices {} --outFileName {} --inputFormat {} --outputFormat {} {}".format(
        matrices,
        outFileName.name,
        inputFormat,
        outputFormat,
        resolutions,
    ).split()

    hicConvertFormat.main(args)


@pytest.mark.parametrize("matrices", [original_matrix_h5, original_matrix_cool])  # required
@pytest.mark.parametrize("outputFormat", ['cool', 'h5'])  # , 'homer', 'ginteractions', 'mcool'])
@pytest.mark.parametrize("resolutions", [''])  # TODO: Check for resolutions
def test_trivial_functionality(
    matrices,
    outputFormat,
    resolutions,
):
    """
        Test for all commandline arguments.
        Options for cool input format are testet seperately.
    """
    from pathlib import Path
    # get suffix of input matrix without the dot
    inputFormat = Path(matrices).suffix[1:]
    # create file corresponding to output format
    outFileName = NamedTemporaryFile(suffix=".{}".format(outputFormat), delete=True)
    outFileName.close()

    args = "--matrices {} --outFileName {} --inputFormat {} --outputFormat {} {}".format(
        matrices,
        outFileName.name,
        inputFormat,
        outputFormat,
        resolutions,
    ).split()

    hicConvertFormat.main(args)

    test = hm.hiCMatrix(matrices)

    new = hm.hiCMatrix(outFileName.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)

    nt.assert_equal(len(new.cut_intervals), len(test.cut_intervals))

    cut_interval_new_ = []
    cut_interval_test_ = []
    for x in new.cut_intervals:
        cut_interval_new_.append(x[:3])
    for x in test.cut_intervals:
        cut_interval_test_.append(x[:3])

    nt.assert_equal(cut_interval_new_, cut_interval_test_)
    os.unlink(outFileName.name)
