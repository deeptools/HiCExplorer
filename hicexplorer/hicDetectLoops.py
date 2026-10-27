"""hicDetectLoops: the tool is the C++ executable hicDetectLoops.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names keep their 3.7.6 signatures and behaviour
and call hicexplorer._hicx.hicDetectLoops (cpp/python/src/api/hicDetectLoops.cpp,
cpp/core/include/hicx/detect_loops.hpp); the work 3.7.6 spread over processes
runs in C++ threads of this process.
"""
import argparse  # noqa: F401
from multiprocessing import Process, Queue  # noqa: F401
from multiprocessing.sharedctypes import Array, RawArray  # noqa: F401
from copy import deepcopy  # noqa: F401
import logging
log = logging.getLogger(__name__)
import time  # noqa: F401
import gc  # noqa: F401
import cooler  # noqa: F401
import numpy as np  # noqa: F401
from scipy.sparse import csr_matrix, triu  # noqa: F401
from scipy.stats import anderson_ksamp, ranksums  # noqa: F401
from scipy.stats import nbinom  # noqa: F401
import fit_nbinom  # noqa: F401

from hicmatrix import HiCMatrix as hm  # noqa: F401
from hicmatrix.lib import MatrixFileHandler  # noqa: F401
from hicexplorer._version import __version__  # noqa: F401
from hicexplorer.utilities import check_cooler  # noqa: F401
from hicexplorer.hicPlotMatrix import translate_region  # noqa: F401
from hicexplorer.lib import cnb  # noqa: F401
from inspect import currentframe  # noqa: F401
import traceback  # noqa: F401


# compute_loops calls them through this module
from hicexplorer.utilities import obs_exp_matrix, obs_exp_matrix_non_zero  # noqa: F401

from hicexplorer._hicx import hicDetectLoops as _impl


def get_linenumber():
    return _impl.get_linenumber()


def parse_arguments(args=None):
    """The argparse parser of hicDetectLoops, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDetectLoops')


def create_distance_distribution(pDataObsExp, pDistances, pWindowSize, pMinDistance, pMaxDistance, pQueue):
    return _impl.create_distance_distribution(pDataObsExp, pDistances, pWindowSize, pMinDistance, pMaxDistance, pQueue)


def compute_p_values_mask(pGenomicDistanceDistributionsObsExp, pGenomicDistanceDistributionsKeyList,
                          pPValuePreselection, pGenomicDistanceDistributionPosition, pResolution,
                          pMinimumInteractionsThreshold, pObsExpThreshold, pQueue):
    return _impl.compute_p_values_mask(pGenomicDistanceDistributionsObsExp, pGenomicDistanceDistributionsKeyList,
                                       pPValuePreselection, pGenomicDistanceDistributionPosition, pResolution,
                                       pMinimumInteractionsThreshold, pObsExpThreshold, pQueue)


def compute_long_range_contacts(pHiCMatrix, pObsExpMatrix, pWindowSize,
                                pPValue, pPeakWindowSize,
                                pPValuePreselection,
                                pMinimumInteractionsThreshold,
                                pObsExpThreshold, pThreads):
    """
        This function computes the loops by:
            - decreasing the search space by removing values with p-values > pPValuePreselection
            - decreasing the search space with the pMaximumInteractionPercentageThreshold: removes all values with less interactions as maximum value of their distance  * pMaximumInteractionPercentageThreshold
            - calls neighborhood_merge to merging candidates which share a neighborhood to one candidate
            - calling candidate_region_test to test the neighborhood of a candidate against its peak to detect significant peaks

        Input:
            - pHiCMatrix: original interaction matrix for neighborhood and peak region subset selection
            - pWindowSize: integer, the size of (2*pWindowSize)^2 around a candidate defines its neighborhood. It is used for
                    a) merging candidates and their neighborhoods to one candidate per neighborhood which appear in this region.
                    b) same neighborhood is used to test the peak region against the neighborhood for significant difference (candidate_region_test)
            - pMinimumInteractionsThreshold: float, remove candidates with less interactions
            - pPValue: float, test rejection level for H0 and FDR correction
            - pPValuePreselection: float, p-value for negative binomial
            - pPeakWindowSize: integer, size of the peak region: (2*pPeakWindowSize)^2. Needs to be smaller than pWindowSize

        Returns:
            - A list of detected loops [(x,y)] and x, y are matrix index values
            - An associated list of p-values
    """
    return _impl.compute_long_range_contacts(pHiCMatrix, pObsExpMatrix, pWindowSize, pPValue, pPeakWindowSize,
                                             pPValuePreselection, pMinimumInteractionsThreshold,
                                             pObsExpThreshold, pThreads)


def neighborhood_merge_thread(pCandidateList, pWindowSize, pInteractionCountMatrix, pQueue):
    return _impl.neighborhood_merge_thread(pCandidateList, pWindowSize, pInteractionCountMatrix, pQueue)


def neighborhood_merge(pCandidates, pWindowSize, pInteractionCountMatrix, pThreads):
    """
        Clusters candidates together to one candidate if they share / overlap their neighborhood.
        Implemented in an iterative way, the candidate with the highest interaction count is accepted as candidate for the neighborhood.

        Input:
            - pCandidates: List of candidates
            - pWindowSize: integer, neighborhood size (2*pWindowSize)^2
            - pInteractionCountMatrix: csr_matrix: The interaction count matrix

        Returns:
            - Reduced list of candidates with no more overlapping neighborhoods
    """
    return _impl.neighborhood_merge(pCandidates, pWindowSize, pInteractionCountMatrix, pThreads)


def candidate_region_test_thread(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                                 pPeakWindowSize, pQueue):
    return _impl.candidate_region_test_thread(pHiCMatrix, pCandidates, pWindowSize, pPValue, pPeakWindowSize,
                                              pQueue)


def candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                          pPeakWindowSize, pThreads):
    """
        Tests if a candidate is having a significant peak compared to its neighborhood.
            - smoothes neighborhood in x an y orientation
            - remove candidate if smoothed peak value < pMinimumInteractionsThreshold
            - reject candidate if:
                - mean(peak) < mean(background)
                - max(peak) < max(background)
            - Test background vs peak with Mann-Whitney rank test and reject H0 if pvalue < pPValue
                - Size of background is: (2*pWindowSize)^2 - (2*pPeakWindowSize)^2
            - Apply multi-test Bonferonni based on pPValue

        Input:
            - pHiCMatrix: csr_matrix, interaction matrix to extract candidate neighborhood
            - pCandidates: list of candidates to test for enrichment
            - pWindowSize: integer, size of neighborhood (2*pWindowSize)^2
            - pPValue: float, significance level for Mann-Whitney rank test
            - pMinimumInteractionsThreshold: integer, if smoothed candidate interaction count is less, it will be removed
            - pPeakWindowSize: size of peak region (2*pPeakWindowSize)^2

        Returns:
            - List of accepted candidates
            - List of associated p-values
    """
    return _impl.candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue, pPeakWindowSize, pThreads)


def cluster_to_genome_position_mapping(pHicMatrix, pCandidates, pPValueList, pMaxLoopDistance):
    """
        Maps the computed enriched loops from matrix index values to genomic locations.

        Input:
            - pHicMatrix: hicmatrix object
            - pCandidates: List of detect loops
            - pPValueList: Associated p-values of loops
            - pMaxLoopDistance: integer, exclude detected loops if (x - y) has a larger distance
        Returns:
            List of detect loops in genomic coordinate format
    """
    return _impl.cluster_to_genome_position_mapping(pHicMatrix, pCandidates, pPValueList, pMaxLoopDistance)


def write_bedgraph(pLoops, pOutFileName, pStartRegion=None, pEndRegion=None):
    """
        Writes the detect loops to a bedgraph file.

    """
    return _impl.write_bedgraph(pLoops, pOutFileName, pStartRegion, pEndRegion)


def compute_loops(pHiCMatrix, pRegion, pArgs, pIsCooler, pQueue=None):
    """
        Master function to compute the loops for one chromosome.
            - Removes all regions greater maxLoopSize
            - Calls
            - Writes computed loops to a bedgraph file

        Input:
            - pHiCMatrix: Hi-C interaction matrix object
            - pRegion: Chromosome name
            - pArgs: Argparser object
            - pIsCooler: True / False if matrix is stored in a .cool file
            - pQueue: Queue object for multiprocessing communication with parent process
    """
    return _impl.compute_loops(pHiCMatrix, pRegion, pArgs, pIsCooler, pQueue)


def read_threshold_file(pFile):
    return _impl.read_threshold_file(pFile)


def main(args=None):
    """Run the C++ implementation of hicDetectLoops with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDetectLoops')(args)
