"""hicTransform: a thin layer over the C++ executable hicTransform.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from scipy.sparse import csr_matrix, lil_matrix
import numpy as np

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_non_zero, obs_exp_matrix
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros


import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicTransform, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicTransform')


def _obs_exp_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount):

    if len(pSubmatrix.data) == 0:
        return pSubmatrix
    obs_exp_matrix_ = obs_exp_matrix_lieberman(pSubmatrix, pLengthChromosome, pChromosomeCount)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_))
    # if len(obs_exp_matrix_.data) == 0:
    #     return np.array()
    return obs_exp_matrix_  # .todense()


def _pearson(pSubmatrix):
    if len(pSubmatrix.data) == 0:
        return pSubmatrix
    pearson_correlation_matrix = np.corrcoef(pSubmatrix)
    pearson_correlation_matrix = convertNansToZeros(csr_matrix(pearson_correlation_matrix))
    pearson_correlation_matrix = convertInfsToZeros(csr_matrix(pearson_correlation_matrix))
    # if len(pearson_correlation_matrix.data) == 0:
    # return np.array([[]])
    return pearson_correlation_matrix  # .todense()


def _obs_exp(pSubmatrix):

    if len(pSubmatrix.data) == 0:
        return pSubmatrix
    obs_exp_matrix_ = obs_exp_matrix(pSubmatrix)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_))
    # log.error('obs_exp_matrix_.data {}'.format(obs_exp_matrix_.data))
    # if len(obs_exp_matrix_.data) == 0:
    #     log.debug('No data!')
    #     return np.array([[]])
    return obs_exp_matrix_  # .todense()


def _obs_exp_non_zero(pSubmatrix, ligation_factor):
    if len(pSubmatrix.data) == 0:
        return pSubmatrix
    obs_exp_matrix_ = obs_exp_matrix_non_zero(pSubmatrix, ligation_factor)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_))
    # if len(obs_exp_matrix_.data) == 0:
    # return np.array([[]])
    return obs_exp_matrix_  # .todense()


def main(args=None):
    """Run the C++ implementation of hicTransform with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicTransform')(args)
