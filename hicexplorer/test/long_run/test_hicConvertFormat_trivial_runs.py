import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicConvertFormat
import pytest


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicConvertFormat")
original_matrix_h5 = ROOT + "/small_test_matrix.h5"
original_matrix_cool = ROOT + "/small_test_matrix.cool"

original_matrix_h5_li = ROOT + "/small_test_matrix.h5"


@pytest.mark.parametrize("matrices", [original_matrix_h5, original_matrix_cool])  # , original_matrix_cool, original_matrix_hic])  # required
@pytest.mark.parametrize("outputFormat", ['cool', 'h5', 'homer', 'ginteractions', 'mcool'])
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
