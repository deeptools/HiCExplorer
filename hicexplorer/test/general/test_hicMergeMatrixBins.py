import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicMergeMatrixBins
from hicmatrix import HiCMatrix as hm
from tempfile import NamedTemporaryFile
import os
import numpy.testing as nt


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_correct_matrix():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()

    args = "--matrix {} --numBins 5 " \
        " --outFileName {}".format(ROOT + "small_test_matrix.h5",
                                   outfile.name).split()
    hicMergeMatrixBins.main(args)

    test = hm.hiCMatrix(ROOT + "hicMergeMatrixBins/result.h5")
    new = hm.hiCMatrix(outfile.name)
    nt.assert_equal(test.matrix.data, new.matrix.data)
    nt.assert_equal(test.cut_intervals, new.cut_intervals)

    os.unlink(outfile.name)
