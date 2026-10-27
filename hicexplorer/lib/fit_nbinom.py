"""hicexplorer.lib.fit_nbinom, implemented in C++.

The negative binomial maximum likelihood fit of the fit_nbinom package
(version 1.2, Gokcen Eraslan 2014, maintained by Joachim Wolff 2023), which
HiCExplorer carries itself. ``fit`` keeps its signature and calls
hicexplorer._hicx.lib_fit_nbinom (cpp/python/src/api/lib_fit_nbinom.cpp).
"""
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import lib_fit_nbinom as _impl


def fit(X, initial_params=None):
    # X is a numpy array representing the data
    # initial params is a numpy array representing the initial values of
    # size and prob parameters
    return _impl.fit(X, initial_params)
