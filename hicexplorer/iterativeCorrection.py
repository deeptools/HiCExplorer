"""hicexplorer.iterativeCorrection, implemented in C++.

iterativeCorrection keeps its 3.7.6 signature, side effects and results; the
body calls hicexplorer._hicx.iterativeCorrection
(cpp/python/src/api/iterativeCorrection.cpp,
cpp/core/include/hicx/iterative_correction.hpp).
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import iterativeCorrection as _impl


def iterativeCorrection(matrix, v=None, M=50, tolerance=1e-5, verbose=False):
    """
    adapted from cytonised version in mirnylab
    original code from: ultracorrectSymmetricWithVector
    https://bitbucket.org/mirnylab/mirnylib/src/924bfdf5ed344df32743f4c03157b0ce49c675e6/mirnylib/numutils_new.pyx?at=default

    Main method for correcting DS and SS read data.
    Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction
    :param matrix: a scipy sparse matrix
    :param tolerance: Tolerance is the maximum allowed relative
                      deviation of the marginals.
    """
    return _impl.iterativeCorrection(matrix, v, M, tolerance, verbose)
