import argparse
from multiprocessing import Process, Queue
from copy import deepcopy
import logging
log = logging.getLogger(__name__)
import time

import cooler
import numpy as np
from scipy.sparse import csr_matrix, triu
from scipy.stats import anderson_ksamp, ranksums
from scipy.stats import nbinom
import fit_nbinom

from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
from hicexplorer.hicPlotMatrix import translate_region

from inspect import currentframe


def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Computes enriched regions (peaks) or long range contacts on the given contact matrix.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to compute the loop detection on.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='Outfile name to store the detected loops. The file will in bedgraph format.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--peakWidth', '-pw',
                           type=int,
                           help='The width of the peak region in bins. The square around the peak will include (2 * peakWidth)^2 bins.')

    parserOpt.add_argument('--windowSize', '-w',
                           type=int,
                           help='The window size for the neighborhood region the peak is located in. All values from this region (exclude the values from the peak '
                           ' region) are tested against the peak region for significant difference. The square will have the size of (2 * windowSize)^2 bins')
    parserOpt.add_argument('--pValuePreselection', '-pp',
                           type=float,
                           help='Only candidates with p-values less the given threshold will be considered as candidates. '
                                'For each genomic distance a negative binomial distribution is fitted and for each pixel a p-value given by the cumulative density function is given. '
                                'This does NOT influence the p-value for the neighborhood testing.',
                           default=0.05)
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=float,
                           help='The minimum number of interactions a detected peaks needs to have to be considered.',
                           default=5)
    parserOpt.add_argument('--maximumInteractionPercentageThreshold', '-mip',
                           type=float,
                           help='For each distance the maximum value is considered and all candidates need to have at least \'max_value * maximumInteractionPercentageThreshold\' interactions.',
                           default=0.01)
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.05,
                           help='Rejection level for Anderson-Darling test for H0. H0 is peak region and background have the same distribution.')

    parserOpt.add_argument('--maxLoopDistance',
                           type=int,
                           default=2000000,
                           help='Maximum genomic distance of a loop, usually loops are within a distance of ~2MB.')
    parserOpt.add_argument('--minLoopDistance',
                           type=int,
                           default=100000,
                           help='Minimum genomic distance of a loop to be considered.')

    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes to include in the analysis. If not set, all chromosomes are included.',
                           nargs='+')

    parserOpt.add_argument('--region',
                           help='The format is chr:start-end.',
                           required=False)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads to use, the parallelization is implemented per chromosome.',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--statisticalTest', '-st',
                           help='Which statistical test should be used.',
                           required=False,
                           type=str,
                           default="anderson-darling",
                           choices=['wilcoxon-rank-sum', 'anderson-darling']
                           )

    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def create_distance_distribution(pData, pDistances):
    pGenomicDistanceDistribution = {}
    pGenomicDistanceDistributionPosition = {}
    pGenomicDistanceDistribution_max_value = {}
    for i, distance in enumerate(pDistances):

        if distance in pGenomicDistanceDistribution:
            pGenomicDistanceDistribution[distance].append(pData[i])
            pGenomicDistanceDistributionPosition[distance].append(i)
        else:
            pGenomicDistanceDistribution[distance] = [pData[i]]
            pGenomicDistanceDistributionPosition[distance] = [i]

    for key in pGenomicDistanceDistribution:
        pGenomicDistanceDistribution_max_value[key] = np.max(pGenomicDistanceDistribution[key])
    return pGenomicDistanceDistribution, pGenomicDistanceDistributionPosition, pGenomicDistanceDistribution_max_value

def compute_p_values_mask(pGenomicDistanceDistributions, pGenomicDistanceDistributionsKeyList, 
                            pPValuePreselection, pMask, pGenomicDistanceDistributionPosition, pQueue): 

    # mask = [False] * len(pGenomicDistanceDistributionsKeyList)
    for i, key in enumerate(pGenomicDistanceDistributionsKeyList):

        nbinom_parameters = fit_nbinom.fit(
                np.array(pGenomicDistanceDistributions[key]))
        nbinom_distance = nbinom(
            nbinom_parameters['size'], nbinom_parameters['prob'])
        less_than = np.array(
            pGenomicDistanceDistributions[key]).astype(int) - 1
        mask_less_than = less_than < 0
        less_than[mask_less_than] = 1
        if len(less_than) <= 0:
            continue
        max_element = np.max(less_than)
        if max_element <= 0:
            continue
        max_element = np.max(less_than)
        sum_of_densities = np.zeros(max_element)

        for j in range(max_element):
            if j >= 1:
                sum_of_densities[j] += sum_of_densities[j - 1]
            sum_of_densities[j] += nbinom_distance.pmf(j)

        # if len(sum_of_densities) > less_than - 1:
        p_value = 1 - sum_of_densities[less_than - 1]
        mask_distance = p_value < pPValuePreselection
        # else:

        for j, value in enumerate(mask_distance):
            if value:
                pMask[pGenomicDistanceDistributionPosition[key][j]] = True
    pQueue.put(pMask)
    return

def compute_long_range_contacts(pHiCMatrix, pWindowSize,
                                pMaximumInteractionPercentageThreshold, pPValue, pPeakWindowSize,
                                pPValuePreselection, pStatisticalTest, pMinimumInteractionsThreshold,
                                pMinLoopDistance, pMaxLoopDistance, pThreads):
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
            - pStatisticalTest: str, which statistical test should be used
            - pPeakWindowSize: integer, size of the peak region: (2*pPeakWindowSize)^2. Needs to be smaller than pWindowSize

        Returns:
            - A list of detected loops [(x,y)] and x, y are matrix index values
            - An associated list of p-values
    """
    log.debug('compute_long_range_contacts initial')

    instances, features = pHiCMatrix.matrix.nonzero()
    distance = np.absolute(instances - features)
    mask = [False] * len(distance)
    log.debug('compute distributions initial')

    genomic_distance_distributions, pGenomicDistanceDistributionPosition, pGenomicDistanceDistribution_max_value = create_distance_distribution(
        pHiCMatrix.matrix.data, distance)
    log.debug('compute distributions DONE')

    nbinom_parameters = {}
    genomic_distance_distributions_thread = len(genomic_distance_distributions) // pThreads

    queue = [None] * pThreads
    process = [None] * pThreads
    mask_threads = [None] * pThreads
    # for i, key in enumerate(genomic_distance_distributions):
    # compute_p_values_mask(pGenomicDistanceDistributions, pGenomicDistanceDistributionsKeyList, pMask, pIndex, pPValuePreselection)
    all_data_collected = False
    thread_done = [False] * pThreads
    genomic_keys_list = list(genomic_distance_distributions.keys())
    index_counter = 0
    for i in range(pThreads):

        if i < pThreads - 1:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:(i + 1) * genomic_distance_distributions_thread]
        else:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:]
        index_pos = index_counter
        # index_counter += len(genomic_distance_keys_thread)
        if len(genomic_distance_keys_thread) == 0:
            process[i] = None
            queue[i] = None
            mask_threads[i] = []
            continue
        else:
            queue[i] = Queue()
            process[i] = Process(target=compute_p_values_mask, kwargs=dict(
                pGenomicDistanceDistributions=genomic_distance_distributions,
                pGenomicDistanceDistributionsKeyList=genomic_distance_keys_thread,
                pPValuePreselection=pPValuePreselection,
                pMask = mask,
                pGenomicDistanceDistributionPosition=pGenomicDistanceDistributionPosition,
                pQueue=queue[i]
            )
            )

            process[i].start()

    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                mask_threads[i] = queue[i].get()
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    for mask_i in mask_threads:
        mask = np.logical_or(mask, mask_i)
       

    
        
    log.debug('compute of nbs DONE')

    peak_interaction_threshold_array = np.zeros(len(distance))
    for i, key in enumerate(distance):
        peak_interaction_threshold_array[i] = pGenomicDistanceDistribution_max_value[key] * pMaximumInteractionPercentageThreshold
    mask_interactions = pHiCMatrix.matrix.data > peak_interaction_threshold_array

    mask = np.logical_and(mask, mask_interactions)

    mask_interactions_hard_threshold = pHiCMatrix.matrix.data > pMinimumInteractionsThreshold
    mask = np.logical_and(mask, mask_interactions_hard_threshold)

    instances = instances[mask]
    features = features[mask]

    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])
    log.debug('compute of neigborhood clustering')

    # Clean neighborhood, results in one candidate per neighborhood
    number_of_candidates = 0
    while number_of_candidates != len(candidates):
        number_of_candidates = len(candidates)

        candidates, mask = neighborhood_merge(
            candidates, pWindowSize, pHiCMatrix.matrix, pMinLoopDistance, pMaxLoopDistance, pThreads)

        if len(candidates) == 0:
            return None, None
    log.debug('call of test')

    candidates, p_value_list = candidate_region_test(
        pHiCMatrix.matrix, candidates, pWindowSize, pPValue,
        pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pStatisticalTest)

    # candidates, p_value_list = candidate_region_test(
    #     pHiCMatrix.matrix, candidates, pWindowSize, pPValue,
    #     peak_interaction_threshold_array, pPeakWindowSize, pStatisticalTest)

    return candidates, p_value_list


def filter_duplicates(pCandidates):
    """
        Removes duplicated candidates in a list.

        Input:
            - pCandidates: List of candidates that may contain duplicates
        Returns:
            - List of candidates without duplicates
    """
    mask = []
    seen_values = {}
    for i, candidate in enumerate(pCandidates):
        if candidate[0] in seen_values:
            if candidate[1] in seen_values[candidate[0]]:
                mask.append(False)
            else:
                seen_values[candidate[0]].append(candidate[1])
                mask.append(True)
        else:
            seen_values[candidate[0]] = [candidate[1]]
            mask.append(True)

    return mask


def neighborhood_merge_thread(pCandidateList, pWindowSize, pInteractionCountMatrix, pMinLoopDistance, pMaxLoopDistance, pQueue):

    x_max = pInteractionCountMatrix.shape[0]
    y_max = pInteractionCountMatrix.shape[1]
    new_candidate_list = []

    for candidate in pCandidateList:
        start_x = candidate[0] - \
            pWindowSize if candidate[0] - pWindowSize > 0 else 0
        end_x = candidate[0] + pWindowSize if candidate[0] + \
            pWindowSize < x_max else x_max
        start_y = candidate[1] - \
            pWindowSize if candidate[1] - pWindowSize > 0 else 0
        end_y = candidate[1] + pWindowSize if candidate[1] + \
            pWindowSize < y_max else y_max

        neighborhood = pInteractionCountMatrix[start_x:end_x,
                                               start_y:end_y].toarray().flatten()
        if len(neighborhood) == 0:
            continue
        argmax = np.argmax(neighborhood)
        x = argmax // (pWindowSize * 2)
        y = argmax % (pWindowSize * 2)

        candidate_x = (candidate[0] - pWindowSize) + \
            x if (candidate[0] - pWindowSize + x) < x_max else x_max - 1
        candidate_y = (candidate[1] - pWindowSize) + \
            y if (candidate[1] - pWindowSize + y) < y_max else y_max - 1
        if candidate_x < 0 or candidate_y < 0:
            continue
        if np.absolute(candidate_x - candidate_y) < pMinLoopDistance or np.absolute(candidate_x - candidate_y) > pMaxLoopDistance:
            continue
        new_candidate_list.append([candidate_x, candidate_y])

    pQueue.put(new_candidate_list)
    return

def neighborhood_merge(pCandidates, pWindowSize, pInteractionCountMatrix, pMinLoopDistance, pMaxLoopDistance, pThreads):
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
    # x_max = pInteractionCountMatrix.shape[0]
    # y_max = pInteractionCountMatrix.shape[1]
    new_candidate_list = []

    new_candidate_list_threads = [None] *  pThreads
    interactionFilesPerThread = len(pCandidates) // pThreads
    all_data_collected = False
    queue = [None] *  pThreads
    process = [None] *  pThreads
    thread_done = [False] *  pThreads
    length_of_threads = 0
    for i in range(pThreads):

        if i <  pThreads - 1:
            candidateThread = pCandidates[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            candidateThread = pCandidates[i * interactionFilesPerThread:]
        # length_of_threads += len(candidateThread)
        queue[i] = Queue()
        process[i] = Process(target=neighborhood_merge_thread, kwargs=dict(
                            pCandidateList=candidateThread,
                            pWindowSize=pWindowSize,
                            pInteractionCountMatrix=pInteractionCountMatrix,
                            pMinLoopDistance=pMinLoopDistance, 
                            pMaxLoopDistance=pMaxLoopDistance,
                            pQueue=queue[i]
        )
        )
        if len(candidateThread) == 0:
            process[i] = None
            queue[i] = None
            new_candidate_list_threads[i] = []
            continue
        else:
            process[i].start()
    # log.debug('length_of_threads {}'.format(length_of_threads))
    while not all_data_collected:
        for i in range( pThreads):
            if queue[i] is not None and not queue[i].empty():
                new_candidate_list_threads[i] = queue[i].get()
                # rejected_file_names[i] = background_data_thread
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    new_candidate_list = [item for sublist in new_candidate_list_threads for item in sublist]


    # for candidate in pCandidates:

        # get neighborhood out of pInteractionCountMatrix matrix
        # start_x = candidate[0] - \
        #     pWindowSize if candidate[0] - pWindowSize > 0 else 0
        # end_x = candidate[0] + pWindowSize if candidate[0] + \
        #     pWindowSize < x_max else x_max
        # start_y = candidate[1] - \
        #     pWindowSize if candidate[1] - pWindowSize > 0 else 0
        # end_y = candidate[1] + pWindowSize if candidate[1] + \
        #     pWindowSize < y_max else y_max

        # neighborhood = pInteractionCountMatrix[start_x:end_x,
        #                                        start_y:end_y].toarray().flatten()
        # if len(neighborhood) == 0:
        #     continue
        # argmax = np.argmax(neighborhood)
        # x = argmax // (pWindowSize * 2)
        # y = argmax % (pWindowSize * 2)

        # candidate_x = (candidate[0] - pWindowSize) + \
        #     x if (candidate[0] - pWindowSize + x) < x_max else x_max - 1
        # candidate_y = (candidate[1] - pWindowSize) + \
        #     y if (candidate[1] - pWindowSize + y) < y_max else y_max - 1
        # if candidate_x < 0 or candidate_y < 0:
        #     continue
        # if np.absolute(candidate_x - candidate_y) < pMinLoopDistance or np.absolute(candidate_x - candidate_y) > pMaxLoopDistance:
        #     continue
        # new_candidate_list.append([candidate_x, candidate_y])


    mask = filter_duplicates(new_candidate_list)

    if mask is not None and len(mask) == 0:
        return [], []
    if mask is None:
        return [], []
    mask = np.array(mask)
    pCandidates = np.array(new_candidate_list)
    # log.debug('type of mask: {}'.format(type(mask)))
    # log.debug('mask: {}'.format(mask))

    pCandidates = pCandidates[mask]
    return pCandidates, mask


def get_test_data(pNeighborhood, pVertical):
    x_len = len(pNeighborhood)
    y_len = len(pNeighborhood[0])

    x_third = x_len // 3
    y_third = y_len // 3

    return_list = []
    if pVertical:
        return_list.append(pNeighborhood.T[0:y_third].flatten())
        return_list.append(pNeighborhood.T[y_third:y_third + y_third].flatten())
        return_list.append(pNeighborhood.T[y_third + y_third:].flatten())
    else:
        return_list.append(pNeighborhood[0:x_third].flatten())
        return_list.append(pNeighborhood[x_third:x_third + x_third].flatten())
        return_list.append(pNeighborhood[x_third + x_third:].flatten())

    return return_list

def candidate_region_test_thread(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                          pMinimumInteractionsThreshold, pPeakWindowSize, pQueue, pStatisticalTest=None):
    mask = []
    pvalues = []

    x_max = pHiCMatrix.shape[0]
    y_max = pHiCMatrix.shape[1]
    for i, candidate in enumerate(pCandidates):
        
        if (candidate[0] - pWindowSize) > 0:
            start_x = candidate[0] - pWindowSize
            peak_x = pWindowSize - 1
        else:
            start_x = 0
            peak_x = pWindowSize - candidate[0] - 1

        if (candidate[1] - pWindowSize) > 0:
            start_y = candidate[1] - pWindowSize
            peak_y = pWindowSize - 1
        else:
            start_y = 0
            peak_y = pWindowSize - candidate[1] - 1

        end_x = candidate[0] + pWindowSize if candidate[0] + \
            pWindowSize < x_max else x_max

        end_y = candidate[1] + pWindowSize if candidate[1] + \
            pWindowSize < y_max else y_max

        neighborhood = pHiCMatrix[start_x:end_x,
                                  start_y:end_y].toarray()

        neighborhood_old = neighborhood
        for j in range(len(neighborhood)):
            neighborhood[j, :] = smoothInteractionValues(neighborhood[j, :], 5)
        for j in range(len(neighborhood_old[0])):
            neighborhood_old[:, j] = smoothInteractionValues(
                neighborhood[:, j], 5)
        neighborhood = (neighborhood + neighborhood_old) / 2

        peak_region = [peak_x, peak_y]

        # if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold[i]:
        # if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold:
        if neighborhood[peak_region[0], peak_region[1]] < pMinimumInteractionsThreshold:

            mask.append(False)
            continue

        if pPeakWindowSize > pWindowSize:
            log.warning('Neighborhood window size ({}) needs to be larger than peak width({}).'.format(
                pWindowSize, pPeakWindowSize))
            return None, None

        peak = neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize,
                            peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize].flatten()

        background = []
        background.extend(
            list(neighborhood[:peak_region[0] - pPeakWindowSize, :].flatten()))
        background.extend(
            list(neighborhood[peak_region[0] + pPeakWindowSize:, :].flatten()))
        background.extend(
            list(neighborhood[:, :peak_region[1] - pPeakWindowSize].flatten()))
        background.extend(
            list(neighborhood[:, peak_region[1] + pPeakWindowSize:].flatten()))
        background = np.array(background)

        if len(background) < pWindowSize:
            mask.append(False)
            continue
        if len(peak) < pWindowSize:
            mask.append(False)
            continue
        if np.mean(peak) < np.mean(background):
            mask.append(False)
            continue
        if np.max(peak) < np.max(background):
            mask.append(False)
            continue
        # if np.min(background) * 5 > np.max(peak):
        #     mask.append(False)
        #     continue
        if pStatisticalTest == 'wilcoxon-rank-sum':
            # test vertical
            test_list = get_test_data(neighborhood, pVertical=True)
            statistic, significance_level_test1 = ranksums(sorted(test_list[0]), sorted(test_list[1]))
            statistic, significance_level_test2 = ranksums(sorted(test_list[1]), sorted(test_list[2]))
            if significance_level_test1 <= pPValue or significance_level_test2 <= significance_level_test2:
                test_list = get_test_data(neighborhood, pVertical=False)
                statistic, significance_level_test1 = ranksums(sorted(test_list[0]), sorted(test_list[1]))
                statistic, significance_level_test2 = ranksums(sorted(test_list[1]), sorted(test_list[2]))
                if significance_level_test1 <= pPValue or significance_level_test2 <= significance_level_test2:

                    statistic, significance_level = ranksums(sorted(peak), sorted(background))
                    mask.append(True)
                    pvalues.append(significance_level)

                    continue
            mask.append(False)
            continue
        else:
            test_list = get_test_data(neighborhood, pVertical=True)
            _, _, significance_level_test1 = anderson_ksamp([sorted(test_list[0]), sorted(test_list[1])])
            _, _, significance_level_test2 = anderson_ksamp([sorted(test_list[1]), sorted(test_list[2])])
            if significance_level_test1 <= pPValue or significance_level_test2 <= significance_level_test2:
                test_list = get_test_data(neighborhood, pVertical=False)
                _, _, significance_level_test1 = anderson_ksamp([sorted(test_list[0]), sorted(test_list[1])])
                _, _, significance_level_test2 = anderson_ksamp([sorted(test_list[1]), sorted(test_list[2])])
                if significance_level_test1 <= pPValue or significance_level_test2 <= significance_level_test2:

                    _, _, significance_level = anderson_ksamp([sorted(peak), sorted(background)])
                    if significance_level <= pPValue:
                        mask.append(True)
                        pvalues.append(significance_level)
                        continue
                    mask.append(False)
                    continue

        mask.append(False)
    # pvalues = []

    pQueue.put([mask, pvalues])
    
def candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                          pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pStatisticalTest=None):
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

    mask = []
    pvalues = []
    x_max = pHiCMatrix.shape[0]
    y_max = pHiCMatrix.shape[1]
    log.debug('candidate_region_test initial: {}'.format(len(pCandidates)))

    if len(pCandidates) == 0:
        return None, None

    pCandidates = np.array(pCandidates)

    mask = []

    mask_thread = [None] *  pThreads
    pvalues_thread = [None] *  pThreads

    interactionFilesPerThread = len(pCandidates) // pThreads
    all_data_collected = False
    queue = [None] *  pThreads
    process = [None] *  pThreads
    thread_done = [False] *  pThreads
    length_of_threads = 0
    for i in range(pThreads):

        if i <  pThreads - 1:
            candidateThread = pCandidates[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            candidateThread = pCandidates[i * interactionFilesPerThread:]
        # length_of_threads += len(candidateThread)
        queue[i] = Queue()
        process[i] = Process(target=candidate_region_test_thread, kwargs=dict(
                            pHiCMatrix=pHiCMatrix, 
                            pCandidates=candidateThread, 
                            pWindowSize=pWindowSize,
                            pPValue=pPValue,
                            pMinimumInteractionsThreshold=pMinimumInteractionsThreshold,
                            pPeakWindowSize=pPeakWindowSize,
                            pQueue=queue[i],
                            pStatisticalTest=pStatisticalTest
        )
        )
# candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                        #   pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pStatisticalTest
        if len(candidateThread) == 0:
            process[i] = None
            queue[i] = None
            mask_threads[i] = []
            continue
        else:
            process[i].start()
    # log.debug('length_of_threads {}'.format(length_of_threads))
    while not all_data_collected:
        for i in range( pThreads):
            if queue[i] is not None and not queue[i].empty():
                result_thread = queue[i].get()
                mask_thread[i], pvalues_thread[i] = result_thread
                # rejected_file_names[i] = background_data_thread
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    mask = [item for sublist in mask_thread for item in sublist]
    pvalues = [item for sublist in pvalues_thread for item in sublist]

   

    # if pStatisticalTest == 'anderson-darling':
    if mask is not None and len(mask) == 0:
        return None, None
    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    log.debug('candidate_region_test done: {}'.format(len(pCandidates)))
    pvalues = np.array(pvalues)

    if pStatisticalTest == 'wilcoxon-rank-sum':
        log.debug('len pCandidates: {}'.format(len(pCandidates)))
        log.debug('len pvalues: {}'.format(len(pvalues)))

        log.debug('Apply fdr')
        pvalues = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
        pvalues_ = sorted(pvalues)
        largest_p_i = 0
        for i, p in enumerate(pvalues_):
            if p <= (pPValue * (i + 1) / len(pvalues_)):
                if p >= largest_p_i:
                    largest_p_i = p
        mask = pvalues < largest_p_i
        pCandidates = pCandidates[mask]
        pvalues = pvalues[mask]

    if len(pCandidates) == 0:
        return None, None

    return pCandidates, pvalues


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
    mapped_cluster = []
    for i, candidate in enumerate(pCandidates):
        chr_x, start_x, end_x, _ = pHicMatrix.getBinPos(candidate[0])
        chr_y, start_y, end_y, _ = pHicMatrix.getBinPos(candidate[1])
        distance = abs(int(start_x) - int(start_y))
        if pMaxLoopDistance is not None and distance > pMaxLoopDistance:
            continue
        mapped_cluster.append(
            (chr_x, start_x, end_x, chr_y, start_y, end_y, pPValueList[i]))
    return mapped_cluster


def write_bedgraph(pLoops, pOutFileName, pStartRegion=None, pEndRegion=None):
    """
        Writes the detect loops to a bedgraph file.

    """
    with open(pOutFileName, 'w') as fh:
        for loop_item in pLoops:
            if pStartRegion and pEndRegion:
                if loop_item[1] >= pStartRegion and loop_item[2] <= pEndRegion \
                        and loop_item[4] >= pStartRegion and loop_item[5] <= pEndRegion:
                    fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % loop_item)
            else:
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % loop_item)


def compute_loops(pHiCMatrix, pRegion, pArgs, pQueue=None):
    """
        Master function to compute the loops for one chromosome.
            - Removes all regions smaller minLoopSize, greater maxLoopSize
            - Calls compute_long_range_contacts
            - Writes computed loops to a bedgraph file

        Input:
            - pHiCMatrix: Hi-C interaction matrix object
            - pRegion: Chromosome name
            - pArgs: Argparser object
            - pQueue: Queue object for multiprocessing communication with parent process
    """
    # log.debug('pRegion {}'.format(pRegion))
    # log.debug('pHiCMatrix.matrix {}'.format(pHiCMatrix.matrix))
    # log.debug('pHiCMatrix.shape {}'.format(pHiCMatrix.matrix.shape))
    if pHiCMatrix.matrix.shape[0] < 5 or pHiCMatrix.matrix.shape[1] < 5:
        log.info('Computed loops for {}: 0'.format(pRegion))

        if pQueue is None:
            return None
        else:
            pQueue.put([None])
            return
    if pArgs.windowSize is None:
        bin_size = pHiCMatrix.getBinSize()
        if 0 < bin_size <= 5000:
            pArgs.windowSize = 12
        elif 5000 < bin_size <= 10000:
            pArgs.windowSize = 10
        elif 10000 < bin_size <= 25000:
            pArgs.windowSize = 8
        elif 25000 < bin_size <= 50000:
            pArgs.windowSize = 7
        else:
            pArgs.windowSize = 6
        log.debug('Setting window size to: {}'.format(pArgs.windowSize))
    if pArgs.peakWidth is None:
        pArgs.peakWidth = pArgs.windowSize - 4
    log.debug('Setting peak width to: {}'.format(pArgs.peakWidth))
    pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
    pHiCMatrix.matrix.eliminate_zeros()
    log.debug('candidates region {} {}'.format(
        pRegion, len(pHiCMatrix.matrix.data)))
    # s
    max_loop_distance = 0
    if pArgs.maxLoopDistance:
        try:
            max_loop_distance = pArgs.maxLoopDistance / pHiCMatrix.getBinSize()
        except Exception:
            log.info('Computed loops for {}: 0'.format(pRegion))

            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        log.debug('pArgs.maxLoopDistance {} max_loop_distance {}'.format(pArgs.maxLoopDistance, max_loop_distance))
        instances, features = pHiCMatrix.matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances > max_loop_distance
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()

        min_loop_distance = 0
    if pArgs.minLoopDistance:
        try:
            min_loop_distance = pArgs.minLoopDistance / pHiCMatrix.getBinSize()
        except Exception:
            log.info('Computed loops for {}: 0'.format(pRegion))

            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        log.debug('pArgs.minLoopDistance {} min_loop_distance {}'.format(pArgs.minLoopDistance, min_loop_distance))

        instances, features = pHiCMatrix.matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances < min_loop_distance
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()

    log.debug('candidates region {} min max boundary {}'.format(
        pRegion, len(pHiCMatrix.matrix.data)))

    candidates, pValueList = compute_long_range_contacts(pHiCMatrix,
                                                         pArgs.windowSize,
                                                         pArgs.maximumInteractionPercentageThreshold,
                                                         pArgs.pValue,
                                                         pArgs.peakWidth,
                                                         pArgs.pValuePreselection,
                                                         pArgs.statisticalTest,
                                                         pArgs.peakInteractionsThreshold,
                                                         min_loop_distance,
                                                         max_loop_distance,
                                                         pArgs.threads)

    if candidates is None:
        log.info('Computed loops for {}: 0'.format(pRegion))

        if pQueue is None:
            return None
        else:
            pQueue.put([None])
            return
    mapped_loops = cluster_to_genome_position_mapping(
        pHiCMatrix, candidates, pValueList, pArgs.maxLoopDistance)
    log.info('Computed loops for {}: {}'.format(pRegion, len(mapped_loops)))

    if pQueue is None:
        return mapped_loops
    else:
        pQueue.put([mapped_loops])
    return


def smoothInteractionValues(pData, pWindowSize):
    '''
        Smoothes pData with a sliding window of pWindowSize
        Adds -pWindowsSize/2 and +pWindowsSize/2 around pData[i] and averages pData[i] by pWindowSize to
        smooth the interaction values.
    '''

    if len(pData) < 2:
        return pData
    window_size = np.int(np.floor(pWindowSize / 2))
    window_size_upstream = window_size
    if pWindowSize % 2 == 0:
        window_size_upstream -= 1

    average_contacts = np.zeros(len(pData))

    # add upstream and downstream, handle regular case
    for i in range(window_size_upstream, len(pData) - window_size):
        start = i - window_size_upstream
        end = i + window_size + 1
        average_contacts[i] = np.mean(pData[start:end])

    # handle border conditions
    for i in range(window_size):
        start = i - window_size_upstream
        if start < 0:
            start = 0
        end = i + window_size + 1

        average_contacts[i] = np.mean(pData[start:end])
        average_contacts[-(i + 1)] = np.mean(pData[-end:])

    return average_contacts


def main(args=None):
    args = parse_arguments().parse_args(args)
    log.info('peak interactions threshold set to {}'.format(
        args.peakInteractionsThreshold))

    if args.region is not None and args.chromosomes is not None:
        log.error('Please choose either --region or --chromosomes.')
        exit(1)
    log.debug('args.matrix {}'.format(args.matrix))
    is_cooler = check_cooler(args.matrix)
    log.debug('is_cooler {}'.format(is_cooler))
    if args.region:
        chrom, region_start, region_end = translate_region(args.region)

        if is_cooler:
            hic_matrix = hm.hiCMatrix(
                pMatrixFile=args.matrix, pChrnameList=[args.region])
        else:
            hic_matrix = hm.hiCMatrix(args.matrix)
            hic_matrix.keepOnlyTheseChr([chrom])
        mapped_loops = compute_loops(hic_matrix, args.region, args)
        write_bedgraph(mapped_loops, args.outFileName,
                       region_start, region_end)

    else:
        mapped_loops = []

        if not is_cooler:
            hic_matrix = hm.hiCMatrix(args.matrix)
            # hic_matrix.keepOnlyTheseChr([chromosome])
            matrix = deepcopy(hic_matrix.matrix)
            cut_intervals = deepcopy(hic_matrix.cut_intervals)

        if args.chromosomes is None:
            # get all chromosomes from cooler file
            if not is_cooler:
                chromosomes_list = list(hic_matrix.chrBinBoundaries)
            else:
                chromosomes_list = cooler.Cooler(args.matrix).chromnames
        else:
            chromosomes_list = args.chromosomes

        if len(chromosomes_list) == 1:
            single_core = True
        else:
            single_core = False

        if single_core:
            for chromosome in chromosomes_list:
                if is_cooler:
                    hic_matrix = hm.hiCMatrix(
                        pMatrixFile=args.matrix, pChrnameList=[chromosome])
                else:
                    hic_matrix.setMatrix(
                        deepcopy(matrix), deepcopy(cut_intervals))
                    hic_matrix.keepOnlyTheseChr([chromosome])
                hic_matrix.maskBins(hic_matrix.nan_bins)
                loops = compute_loops(hic_matrix, chromosome, args)
                if loops is not None:
                    mapped_loops.extend(loops)
        else:
            queue = [None] * args.threads
            process = [None] * args.threads
            all_data_processed = False
            all_threads_done = False
            thread_done = [False] * args.threads
            count_call_of_read_input = 0
            while not all_data_processed or not all_threads_done:
                for i in range(args.threads):
                    if queue[i] is None and not all_data_processed:
                        if count_call_of_read_input >= len(chromosomes_list):
                            all_data_processed = True
                            continue
                        queue[i] = Queue()
                        thread_done[i] = False
                        if is_cooler:
                            hic_matrix = hm.hiCMatrix(pMatrixFile=args.matrix, pChrnameList=[
                                                      chromosomes_list[count_call_of_read_input]])
                        else:
                            hic_matrix.setMatrix(
                                deepcopy(matrix), deepcopy(cut_intervals))
                            hic_matrix.keepOnlyTheseChr(
                                [chromosomes_list[count_call_of_read_input]])
                        if len(hic_matrix.matrix.data) > 0:

                            process[i] = Process(target=compute_loops, kwargs=dict(
                                pHiCMatrix=hic_matrix,
                                pRegion=chromosomes_list[count_call_of_read_input],
                                pArgs=args,
                                pQueue=queue[i]
                            ))
                            process[i].start()

                        else:
                            queue[i] = None
                            thread_done[i] = True
                        if count_call_of_read_input < len(chromosomes_list):
                            count_call_of_read_input += 1
                        else:
                            all_data_processed = True
                    elif queue[i] is not None and not queue[i].empty():
                        result = queue[i].get()
                        if result[0] is not None:
                            mapped_loops.extend(result[0])

                        queue[i] = None
                        process[i].join()
                        process[i].terminate()
                        process[i] = None
                        thread_done[i] = True
                    elif all_data_processed and queue[i] is None:
                        thread_done[i] = True
                    else:
                        time.sleep(1)

                if all_data_processed:
                    all_threads_done = True
                    for thread in thread_done:
                        if not thread:
                            all_threads_done = False
        log.debug('done computing. loops {}'.format(mapped_loops))
        if len(mapped_loops) > 0:
            write_bedgraph(mapped_loops, args.outFileName)

    log.info("Number of detected loops for all regions: {}".format(
        len(mapped_loops)))
