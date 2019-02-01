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

original_matrix_h5_li = ROOT + "/small_test_matrix.h5"


@pytest.mark.parametrize("matrices", [original_matrix_cool])  # required
@pytest.mark.parametrize("outputFormat", ['cool', 'h5', 'homer', 'ginteractions', 'mcool'])
@pytest.mark.parametrize("correction_name", ['weight'])  # need to check hicInfo for more names
@pytest.mark.parametrize("correction_division", ['', '--correction_division'])
@pytest.mark.parametrize("store_applied_correction", ['', '--store_applied_correction'])
@pytest.mark.parametrize("chromosome", ['chrX'])
@pytest.mark.parametrize("enforce_integer", ['', '--enforce_integer'])
@pytest.mark.parametrize("load_raw_values", ['', '--load_raw_values'])
def test_cool_specific_trivial_run(
    matrices,
    outputFormat,
    correction_name,
    correction_division,
    store_applied_correction,
    chromosome,
    enforce_integer,
    load_raw_values,
):
    """
        Cool input format supports some specific options like correction_name, correction_division...
        Therefore, cool input format is explicitly tested in a single test function.
    """
    from pathlib import Path
    # get suffix of input matrix without the dot
    inputFormat = Path(matrices).suffix[1:]
    # create file corresponding to output format
    outFileName = NamedTemporaryFile(suffix="test_ConvertFormat_trivial_run_cool.{}".format(outputFormat), delete=False)
    outFileName.close()

    args = "--matrices {} --outFileName {} --outputFormat {} --inputFormat {} --correction_name {} {} {} --chromosome {} {} {}".format(
        matrices,
        outFileName.name,
        outputFormat,
        inputFormat,
        correction_name,
        correction_division,
        store_applied_correction,
        chromosome,
        enforce_integer,
        load_raw_values,
    ).split()

    hicConvertFormat.main(args)
