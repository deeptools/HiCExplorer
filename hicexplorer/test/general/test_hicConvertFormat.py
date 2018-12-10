
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicConvertFormat
from hicmatrix import HiCMatrix as hm
import numpy.testing as nt
import numpy as np

REMOVE_OUTPUT = True
# DIFF = 60

DELTA_DECIMAL = 0
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicConvertFormat")
original_matrix_h5 = ROOT + "/small_test_matrix.h5"
original_matrix_cool = ROOT + "/small_test_matrix.cool"

# test cases for:
#   - h5 to cool
#   - h5 to homer
#   - h5 to ginteractions
#   - h5's to mcool
#   - cool to h5
#   - hic to cool
#   - cool to h5:
#       - correction name
#       - correction division
#       - store_applied correction
#       - chromosome
#       - enforce integer
#   - resolutions


def test_hicConvertFormat_h5_to_cool():

    # original_matrix = ''
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat cool".format(original_matrix_h5, outfile.name).split()
    hicConvertFormat.main(args)

    test = hm.hiCMatrix(original_matrix_cool)
    # print(outfile.name)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    # os.unlink(outfile.name)


def test_hicConvertFormat_h5_to_cool_enforce_integer():

    # original_matrix = ''
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat cool ".format(original_matrix_h5, outfile.name).split()
    hicConvertFormat.main(args)

    test = hm.hiCMatrix(original_matrix_cool)
    # print(outfile.name)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=0)
    # print('issubclass(test.matrix.data.dtype.type, np.integer) {}'.format(issubclass(test.matrix.data.dtype.type, np.integer)))
    assert issubclass(test.matrix.data.dtype.type, np.integer)


def test_hicConvertFormat_h5_to_homer():
    pass
    # os.unlink(outfile.name)


def test_hicConvertFormat_h5_to_ginteractions():
    pass

    # os.unlink(outfile.name)


def test_hicConvertFormat_h5_to_mcool():
    pass

    # os.unlink(outfile.name)


def test_hicConvertFormat_cool_to_h5():

    # original_matrix = ''
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat cool --outputFormat h5".format(original_matrix_cool, outfile.name).split()
    hicConvertFormat.main(args)

    test = hm.hiCMatrix(original_matrix_h5)
    # print(outfile.name)
    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
