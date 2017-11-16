from hicexplorer import hicExport
from hicexplorer import HiCMatrix as hm

from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/hicExport/"

NPZ_MATRIX = ROOT + "matrix.npz"
COOL_MATRIX = ROOT + "matrix.cool"
H5_MATRIX = ROOT + "matrix.h5"


def run_compare(pInputFile, pInputFormat, pOutputFormat, pChrNameList=None):

    outfile = NamedTemporaryFile(suffix='.' + pOutputFormat, delete=False)
    outfile.close()
    inputMatrix = pInputFile

    args = "--inFile {} --inputFormat {} "\
        "--outFileName {} " \
        "--outputFormat {}".format(inputMatrix,
                                   pInputFormat,
                                   outfile.name,
                                   pOutputFormat).split()
    hicExport.main(args)

    test = hm.hiCMatrix(inputMatrix)
    new = hm.hiCMatrix(outfile.name)

    nt.assert_equal(new.matrix.data, test.matrix.data)
    nt.assert_equal(new.cut_intervals, test.cut_intervals)
    os.unlink(outfile.name)
    return True


def test_npz_to_cool():
    return run_compare(pInputFile=NPZ_MATRIX, pInputFormat='npz', pOutputFormat='cool')


def test_npz_to_h5():
    return run_compare(pInputFile=NPZ_MATRIX, pInputFormat='npz', pOutputFormat='h5')


def test_cool_to_npz():
    return run_compare(pInputFile=COOL_MATRIX, pInputFormat='cool', pOutputFormat='npz')


def test_cool_to_h5():
    return run_compare(pInputFile=COOL_MATRIX, pInputFormat='cool', pOutputFormat='h5')


def test_h5_to_npz():
    return run_compare(pInputFile=H5_MATRIX, pInputFormat='h5', pOutputFormat='npz')


def test_h5_to_cool():
    return run_compare(pInputFile=H5_MATRIX, pInputFormat='h5', pOutputFormat='cool')
