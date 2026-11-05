"""hicexplorer.lib.cnb, implemented in C++.

The continuous generalisation of the negative binomial distribution; pdf and
cdf keep their 3.7.6 signatures and call hicexplorer._hicx.lib_cnb
(cpp/python/src/api/lib_cnb.cpp, cpp/core/include/hicx/cnb.hpp).
"""
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import lib_cnb as _impl


def pdf(pX, pR, pP):
    """
    PDF for a continuous generalization of NB distribution
    """
    return _impl.pdf(pX, pR, pP)


def cdf(pX, pR, pP):
    """
    Cumulative density function of a continuous generalization of NB distribution
    """
    return _impl.cdf(pX, pR, pP)
