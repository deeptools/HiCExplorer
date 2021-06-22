import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicCompareMatrices
from hicmatrix import HiCMatrix as hm
from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))), "test_data/")


# I doubled the values in ./hicexplorer/test/test_data/hicConvertFormat/small_test_matrix_chr4.cool
# hicNormalize -m ./hicexplorer/test/test_data/small_test_matrix.cool --normalize multiplicative -mv 2 -o ./hicexplorer/test/test_data/hicCompareMatrices/small_test_matrix_twice.cool

def test_hicCompareMatrices_doubleMinusOneEqual0():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} {} --outFileName {} " \
           "--operation diff" \
           .format(ROOT + "hicCompareMatrices/small_test_matrix_twice.cool",
                   ROOT + "small_test_matrix.cool",
                   outfile.name).split()

    compute(hicCompareMatrices.main, args, 5)
    input = hm.hiCMatrix(
        ROOT + "small_test_matrix.cool")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal([], new.matrix.data)
    nt.assert_equal(input.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)


def test_hicCompareMatrices_noNorm_doubleMinusOneEqualOne():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()

    args = "--matrices {} {} --outFileName {} " \
           "--operation diff --noNorm" \
           .format(ROOT + "hicCompareMatrices/small_test_matrix_twice.cool",
                   ROOT + "small_test_matrix.cool",
                   outfile.name).split()

    compute(hicCompareMatrices.main, args, 5)
    input = hm.hiCMatrix(
        ROOT + "small_test_matrix.cool")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(input.matrix.data, new.matrix.data)
    nt.assert_equal(input.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)
