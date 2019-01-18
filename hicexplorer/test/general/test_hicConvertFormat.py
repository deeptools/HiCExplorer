
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os.path
from tempfile import NamedTemporaryFile
from hicexplorer import hicConvertFormat
from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
import gzip
from scipy.sparse import triu
import numpy.testing as nt
import numpy as np

REMOVE_OUTPUT = True
# DIFF = 60

DELTA_DECIMAL = 0
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/hicConvertFormat")
original_matrix_h5 = ROOT + "/small_test_matrix.h5"
original_matrix_cool = ROOT + "/small_test_matrix.cool"
original_matrix_cool_chr4 = ROOT + "/small_test_matrix_chr4.cool"


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

    outfile = NamedTemporaryFile(suffix='.homer', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat cool --outputFormat homer ".format(original_matrix_cool_chr4, outfile.name).split()
    hicConvertFormat.main(args)

    test = hm.hiCMatrix(original_matrix_cool_chr4)
    f = gzip.open(outfile.name, 'rb')
    file_content = f.read()
    outfile2 = NamedTemporaryFile(suffix='.homer', delete=False)
    outfile2.close()
    with open(outfile2.name, 'wb') as matrix_file:
        matrix_file.write(file_content)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='homer', pMatrixFile=outfile2.name)

    _matrix, cut_intervals, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

    nt.assert_array_almost_equal(test.matrix.data, _matrix.data, decimal=0)


def test_hicConvertFormat_h5_to_ginteractions():
    outfile = NamedTemporaryFile(suffix='.ginteractions', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat ginteractions ".format(original_matrix_h5, outfile.name).split()
    hicConvertFormat.main(args)

    # os.unlink(outfile.name)


def test_hicConvertFormat_h5_to_mcool():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} --outFileName {} --inputFormat h5 --outputFormat mcool -r 10000 100000 200000 ".format(original_matrix_h5, outfile.name).split()
    hicConvertFormat.main(args)

    new1 = hm.hiCMatrix(outfile.name + '::/resolutions/10000')  # noqa: F841
    new2 = hm.hiCMatrix(outfile.name + '::/resolutions/100000')  # noqa: F841
    new3 = hm.hiCMatrix(outfile.name + '::/resolutions/200000')  # noqa: F841

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


def test_hicConvertFormat_hicpro_to_cool():

    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    hicprofile = ROOT + '/test_matrix.hicpro'
    bedfile = ROOT + '/test_matrix.bed'
    args = "--matrices {} --outFileName {} --inputFormat hicpro --outputFormat cool --bedFileHicpro {}".format(hicprofile, outfile.name, bedfile).split()
    hicConvertFormat.main(args)

    # test = hm.hiCMatrix(original_matrix_cool)
    # print(outfile.name)
    new = hm.hiCMatrix(outfile.name)

    matrixFileHandlerInput = MatrixFileHandler(pFileType='hicpro', pMatrixFile=hicprofile,
                                               pBedFileHicPro=bedfile)

    _matrix, cut_intervals, nan_bins, \
        distance_counts, correction_factors = matrixFileHandlerInput.load()

    new.matrix = triu(new.matrix)
    nt.assert_array_almost_equal(new.matrix.data, _matrix.data, decimal=0)
