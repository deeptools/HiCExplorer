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
from hicexplorer.lib import cnb
from inspect import currentframe

from hicexplorer.utilities import obs_exp_matrix

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
                        #    type=float,
                           help='Only candidates with p-values less the given threshold will be considered as candidates. '
                                'For each genomic distance a negative binomial distribution is fitted and for each pixel a p-value given by the cumulative density function is given. '
                                'This does NOT influence the p-value for the neighborhood testing. Can a single value or a threshold file created by hicCreateThresholdFile.',
                           default=0.075)
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=float,
                           help='The minimum number of interactions a detected peaks needs to have to be considered.',
                           default=1.0)
    parserOpt.add_argument('--maximumInteractionPercentageThreshold', '-mip',
                           type=float,
                           help='For each distance the maximum value is considered and all candidates need to have at least \'max_value * maximumInteractionPercentageThreshold\' interactions.',
                           default=0.01)
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.025,
                           help='Rejection level for Anderson-Darling or Wilcoxon-rank sum test for H0. H0 is peak region and background have the same distribution.')
    parserOpt.add_argument('--qValue', '-q',
                           type=float,
                           default=0.05,
                           help='Level for FDR per chromosome.')

    parserOpt.add_argument('--maxLoopDistance',
                           type=int,
                           default=2000000,
                           help='Maximum genomic distance of a loop, usually loops are within a distance of ~2MB.')
    parserOpt.add_argument('--minLoopDistance',
                           type=int,
                           default=10000,
                           help='Minimum genomic distance of a loop to be considered.')
    parserOpt.add_argument('--maxVariance',
                           type=float,
                           default=0.1,
                           help='Maximum variance in a normalized neighborhood of a candidate.')
    parserOpt.add_argument('--minVariance',
                           type=float,
                           default=0.01,
                           help='Minimum variance in a normalized neighborhood of a candidate.')
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
    parserOpt.add_argument('--threadsPerChromosome', '-tpc',
                           help='Number of threads to use per parallel thread processing a chromosome. E.g. --threads = 4 and --threadsPerChromosome = 4 makes 4 * 4 = 16 threads in total.',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--statisticalTest', '-st',
                           help='Which statistical test should be used.',
                           required=False,
                           type=str,
                           default="wilcoxon-rank-sum",
                           choices=['wilcoxon-rank-sum', 'anderson-darling']
                           )

    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def create_distance_distribution(pData, pDistances, pWindowSize, pMinimumInteractionsThreshold, pMinDistance, pMaxDistance, pQueue):
    log.debug('distance distribution')
    pGenomicDistanceDistribution = {}
    pGenomicDistanceDistributionNeighborhoodWindow = {}

    pGenomicDistanceDistributionPosition = {}
    # pGenomicDistanceDistribution_max_value = {}

    # min_distance = pDistances.min()
    # max_distance = pDistances.max()
    # mask_interactions = pData >= pMinimumInteractionsThreshold
    for distance in range(pMinDistance, pMaxDistance, 1):
        mask = pDistances == distance
        # data = 
        # if len(data) > 0 and data.max() < pMinimumInteractionsThreshold:
        #     # pGenomicDistanceDistributionPosition[distance] = []
        #     # pGenomicDistanceDistribution[distance] = []
        #     # pGenomicDistanceDistributionNeighborhoodWindow[distance] = []
        #     pGenomicDistanceDistribution_max_value[distance] = 0
        #     continue
        # pGenomicDistanceDistributionNeighborhoodWindow[distance] = data
        # mask = np.logical_and(mask, mask_interactions)
        pGenomicDistanceDistribution[distance] = pData[mask]
        pGenomicDistanceDistributionPosition[distance] = np.argwhere(mask == True).flatten()
        # if len(pGenomicDistanceDistributionPosition[distance]) > 0:
        #     pGenomicDistanceDistribution_max_value[distance] = pGenomicDistanceDistribution[distance].max()
        # else:
        #     pGenomicDistanceDistribution_max_value[distance] = 0
        # mask_upper = pDistances <= distance + pWindowSize
        # mask_lower = distance - pWindowSize < pDistances
        # mask = np.logical_and(mask_lower, mask_upper)

        # mask_truncate = 5 < pData 
        # mask = np.logical_and(mask_truncate, mask)


    log.debug('distance distribution...DNE')

    pQueue.put([pGenomicDistanceDistribution, pGenomicDistanceDistributionPosition, None, None])
    return


def compute_p_values_mask(pGenomicDistanceDistributions, pGenomicDistanceDistributionsKeyList, 
                          pPValuePreselection, pMask, pGenomicDistanceDistributionPosition, pResolution, pMinimumInteractionsThreshold, pQueue):

    # mask = [False] * len(pGenomicDistanceDistributionsKeyList)
    log.debug('p-values...')

    if len(pGenomicDistanceDistributionsKeyList) == 0:
        pQueue.put(pMask)
        return

    for i, key in enumerate(pGenomicDistanceDistributionsKeyList):

        data = np.array(pGenomicDistanceDistributions[key])

        # do not fit and not compute any p-value if all values on this distance are small than the pMinimumInteractionsThreshold
        mask = data >= pMinimumInteractionsThreshold
        if np.sum(mask) == 0:#.max() < pMinimumInteractionsThreshold:
            continue
        nbinom_parameters = fit_nbinom.fit(data)

        # less_than = np.array(pGenomicDistanceDistributions[key])
        # mask = less_than >= pMinimumInteractionsThreshold
        # less_than = less_than[mask]
        p_value = 1 - cnb.cdf(data[mask], nbinom_parameters['size'], nbinom_parameters['prob'])

        if isinstance(pPValuePreselection, float):
            mask_distance = p_value < pPValuePreselection
        elif isinstance(pPValuePreselection, dict):
            key_genomic = int(key * pResolution)
            mask_distance = p_value < pPValuePreselection[key_genomic]
        # else:
        j = 0
        for k, value in enumerate(mask):
            if value:
                if mask_distance[j]:
                    pMask[pGenomicDistanceDistributionPosition[key][k]] = True
                j += 1
        # for j, value in enumerate(mask_distance):
        #     if value:
        #         pMask[pGenomicDistanceDistributionPosition[key][j]] = True
    log.debug('p-values...DONE')

    pQueue.put(pMask)
    return


def compute_long_range_contacts(pHiCMatrix, pWindowSize,
                                pMaximumInteractionPercentageThreshold, pPValue, pQValue, pPeakWindowSize,
                                pPValuePreselection, pStatisticalTest, pMinimumInteractionsThreshold,
                                pMinLoopDistance, pMaxLoopDistance, pThreads,  pMinVariance, pMaxVariance):
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

    instances, features = pHiCMatrix.matrix.nonzero()
    distance = np.absolute(instances - features)
    queue = [None] * pThreads
    process = [None] * pThreads
    genomic_distance_distributions_thread = [None] * pThreads
    pGenomicDistanceDistributionPosition_thread = [None] * pThreads
    pGenomicDistanceDistribution_max_value_thread = [None] * pThreads
    genomicDistanceDistributionsNeighborhoodWindow_thread = [None] * pThreads
    all_data_collected = False
    thread_done = [False] * pThreads
    min_distance = distance.min()
    max_distance = distance.max()

    distances_per_threads = (max_distance - min_distance) // pThreads
    for i in range(pThreads):

        if i < pThreads - 1:
            min_distance_thread = min_distance + (i * distances_per_threads)
            max_distance_thread = min_distance + ((i+1) * distances_per_threads)
            #  = genomic_keys_list[i * genomic_distance_distributions_thread:(i + 1) * genomic_distance_distributions_thread]
        else:
            # genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:]
            min_distance_thread = min_distance + (i * distances_per_threads)
            max_distance_thread = max_distance + 1
        # index_pos = index_counter
        # index_counter += len(genomic_distance_keys_thread)
        # if len(genomic_distance_keys_thread) == 0:

        #     process[i] = None
        #     queue[i] = None
        #     thread_done[i] = True

        #     mask_threads[i] = [False] * len(distance)
        #     continue
        # else:
        queue[i] = Queue()
        process[i] = Process(target=create_distance_distribution, kwargs=dict(
            pData=pHiCMatrix.matrix.data,
            pDistances=distance,
            pWindowSize=pWindowSize, 
            pMinimumInteractionsThreshold=pMinimumInteractionsThreshold, 
            pMinDistance=min_distance_thread, 
            pMaxDistance=max_distance_thread,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                genomic_distance_distributions_thread[i], pGenomicDistanceDistributionPosition_thread[i], \
                    _, _ = queue[i].get()
                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
                # log.debug("276")
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    genomic_distance_distributions = {}
    pGenomicDistanceDistributionPosition = {}
    pGenomicDistanceDistribution_max_value = {}
    genomicDistanceDistributionsNeighborhoodWindow = {}
    for thread_data in genomic_distance_distributions_thread:
        genomic_distance_distributions = {**genomic_distance_distributions, **thread_data}

    for thread_data in pGenomicDistanceDistributionPosition_thread:
        pGenomicDistanceDistributionPosition = {**pGenomicDistanceDistributionPosition, **thread_data}
    
    # for thread_data in pGenomicDistanceDistribution_max_value_thread:
    #     pGenomicDistanceDistribution_max_value = {**pGenomicDistanceDistribution_max_value, **thread_data}
    
    # for thread_data in genomicDistanceDistributionsNeighborhoodWindow_thread:
    #     genomicDistanceDistributionsNeighborhoodWindow = {**genomicDistanceDistributionsNeighborhoodWindow, **thread_data}
    # genomic_distance_distributions, pGenomicDistanceDistributionPosition, pGenomicDistanceDistribution_max_value, genomicDistanceDistributionsNeighborhoodWindow = create_distance_distribution(
    #     pHiCMatrix.matrix.data, distance, pWindowSize, pMinimumInteractionsThreshold)

    # nbinom_parameters = {}
    mask = [False] * len(distance)

    genomic_distance_distributions_thread = len(genomic_distance_distributions) // pThreads

    queue = [None] * pThreads
    process = [None] * pThreads
    mask_threads = [None] * pThreads
    # for i, key in enumerate(genomic_distance_distributions):
    # compute_p_values_mask(pGenomicDistanceDistributions, pGenomicDistanceDistributionsKeyList, pMask, pIndex, pPValuePreselection)
    all_data_collected = False
    thread_done = [False] * pThreads
    genomic_keys_list = list(genomic_distance_distributions.keys())
    # index_counter = 0
    for i in range(pThreads):

        if i < pThreads - 1:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:(i + 1) * genomic_distance_distributions_thread]
        else:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:]
        # index_pos = index_counter
        # index_counter += len(genomic_distance_keys_thread)
        if len(genomic_distance_keys_thread) == 0:

            process[i] = None
            queue[i] = None
            thread_done[i] = True

            mask_threads[i] = [False] * len(distance)
            continue
        else:
            queue[i] = Queue()
            process[i] = Process(target=compute_p_values_mask, kwargs=dict(
                pGenomicDistanceDistributions=genomic_distance_distributions,
                pGenomicDistanceDistributionsKeyList=genomic_distance_keys_thread,
                # pGenomicDistanceDistributionsNeighborhoodWindow=genomicDistanceDistributionsNeighborhoodWindow,
                pPValuePreselection=pPValuePreselection,
                pMask=mask,
                pGenomicDistanceDistributionPosition=pGenomicDistanceDistributionPosition,
                pResolution=pHiCMatrix.getBinSize(),
                pMinimumInteractionsThreshold=pMinimumInteractionsThreshold,
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
                # log.debug("276")
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    for mask_i in mask_threads:
        # log.debug('true for p-value mask_i: {}'.format(np.sum(mask_i)))
        mask = np.logical_or(mask, mask_i)

    # peak_interaction_threshold_array = np.zeros(len(distance))
    # for i, key in enumerate(distance):
    #     peak_interaction_threshold_array[i] = pGenomicDistanceDistribution_max_value[key] * pMaximumInteractionPercentageThreshold
    # mask_interactions = pHiCMatrix.matrix.data > peak_interaction_threshold_array

    # log.debug('true for p-value: {}'.format(np.sum(mask)))
    # log.debug('true for interaction count  {}'.format(np.sum(mask_interactions)))

    # mask = np.logical_and(mask, mask_interactions)

    mask_interactions_hard_threshold = pHiCMatrix.matrix.data > pMinimumInteractionsThreshold
    mask = np.logical_and(mask, mask_interactions_hard_threshold)
    # log.debug('candidates; {}'.format(len(instances)))
    instances = instances[mask]
    features = features[mask]

    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])
    log.debug('compute of neigborhood clustering: {}'.format(len(candidates)))

    del instances
    del features
    del distance
    del mask

    # Clean neighborhood, results in one candidate per neighborhood
    # number_of_candidates = 0
    # while number_of_candidates != len(candidates):
    #     number_of_candidates = len(candidates)

    candidates, _ = neighborhood_merge(
        candidates, pWindowSize, pHiCMatrix.matrix, pMinLoopDistance, pMaxLoopDistance, pThreads)

    if len(candidates) == 0:
        return None, None
    log.debug('call of test')

    log.info('candidates before testing: {}'.format(len(candidates)))
    candidates, p_value_list = candidate_region_test(
        pHiCMatrix.matrix, candidates, pWindowSize, pPValue, pQValue,
        pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pMinVariance, pMaxVariance, pStatisticalTest)

    # candidates, p_value_list = candidate_region_test(
    #     pHiCMatrix.matrix, candidates, pWindowSize, pPValue,
    #     peak_interaction_threshold_array, pPeakWindowSize, pStatisticalTest)

    return candidates, p_value_list


# def filter_duplicates(pCandidates):
#     """
#         Removes duplicated candidates in a list.

#         Input:
#             - pCandidates: List of candidates that may contain duplicates
#         Returns:
#             - List of candidates without duplicates
#     """
#     mask = []
#     seen_values = {}
#     for i, candidate in enumerate(pCandidates):
#         if candidate[0] in seen_values:
#             if candidate[1] in seen_values[candidate[0]]:
#                 mask.append(False)
#             else:
#                 seen_values[candidate[0]].append(candidate[1])
#                 mask.append(True)
#         else:
#             seen_values[candidate[0]] = [candidate[1]]
#             mask.append(True)

#     return mask


def neighborhood_merge_thread(pCandidateList, pWindowSize, pInteractionCountMatrix, pMinLoopDistance, pMaxLoopDistance, pQueue):

    new_candidate_list = []

    if len(pCandidateList) == 0:
        pQueue.put(new_candidate_list)
        return
    x_max = pInteractionCountMatrix.shape[0]
    y_max = pInteractionCountMatrix.shape[1]

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
        
        if len(neighborhood) > 0 and np.max(neighborhood) == pInteractionCountMatrix[candidate[0],candidate[1]]:
            new_candidate_list.append(candidate)
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
    del pCandidateList
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

    log.debug('pCandidates {}'.format(pCandidates[:10]))
    new_candidate_list = []

    new_candidate_list_threads = [None] * pThreads
    interactionFilesPerThread = len(pCandidates) // pThreads
    all_data_collected = False
    queue = [None] * pThreads
    process = [None] * pThreads
    thread_done = [False] * pThreads
    # length_of_threads = 0
    for i in range(pThreads):

        if i < pThreads - 1:
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
            thread_done[i] = True

            new_candidate_list_threads[i] = []
            continue
        else:
            process[i].start()
    # log.debug('length_of_threads {}'.format(length_of_threads))
    while not all_data_collected:
        for i in range(pThreads):
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

    return new_candidate_list, True


# def get_test_data(pNeighborhood, pVertical):
#     x_len = len(pNeighborhood)
#     y_len = len(pNeighborhood[0])

#     x_third = x_len // 3
#     y_third = y_len // 3

#     return_list = []
#     if pVertical:
#         return_list.append(pNeighborhood.T[0:y_third].flatten())
#         return_list.append(pNeighborhood.T[y_third:y_third + y_third].flatten())
#         return_list.append(pNeighborhood.T[y_third + y_third:].flatten())
#     else:
#         return_list.append(pNeighborhood[0:x_third].flatten())
#         return_list.append(pNeighborhood[x_third:x_third + x_third].flatten())
#         return_list.append(pNeighborhood[x_third + x_third:].flatten())

#     return return_list


def candidate_region_test_thread(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                                 pMinimumInteractionsThreshold, pPeakWindowSize, pQueue, pMinVariance, pMaxVariance, pStatisticalTest=None, ):
    # neighborhood_list = []
    mask = []
    pvalues = []
    if len(pCandidates) == 0:
        pQueue.put([mask, pvalues])
        return
    x_max = pHiCMatrix.shape[0]
    y_max = pHiCMatrix.shape[1]
    for i, candidate in enumerate(pCandidates):

        if (candidate[0] - pWindowSize) > 0:
            start_x = candidate[0] - pWindowSize
        else:
            start_x = 0

        if (candidate[1] - pWindowSize) > 0:
            start_y = candidate[1] - pWindowSize
        else:
            start_y = 0

        end_x = candidate[0] + pWindowSize+1 if candidate[0] + \
            pWindowSize+1 < x_max else x_max

        end_y = candidate[1] + pWindowSize+1 if candidate[1] + \
            pWindowSize+1 < y_max else y_max

        

        neighborhood = pHiCMatrix[start_x:end_x,
                                  start_y:end_y].toarray()
        if len(neighborhood) == 0:
            mask.append(False)
            continue
        # get index of original candidate
        peak_region = np.array(np.where(neighborhood == pHiCMatrix[candidate[0], candidate[1]])).flatten()
        
        # neighborhood_mask = 1.0 > neighborhood
        # neighborhood = neighborhood[neighborhood_mask]
        # neighborhood_list.append(neighborhood.flatten())
        # neighborhood_old = neighborhood
        # for j in range(len(neighborhood)):
        #     neighborhood[j, :] = smoothInteractionValues(neighborhood[j, :], 5)
        # for j in range(len(neighborhood_old[0])):
        #     neighborhood_old[:, j] = smoothInteractionValues(
        #         neighborhood[:, j], 5)
        # neighborhood = (neighborhood + neighborhood_old) / 2

        # peak_region = [peak_x, peak_y]

        # if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold[i]:
        # if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold:
        # if neighborhood[peak_region[0], peak_region[1]] < pMinimumInteractionsThreshold:

        #     mask.append(False)
        #     continue

        if pPeakWindowSize > pWindowSize:
            log.warning('Neighborhood window size ({}) needs to be larger than peak width({}).'.format(
                pWindowSize, pPeakWindowSize))
            return None, None

        peak = neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize+1,
                            peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize+1].flatten()

        background = []
        # top to peak 
        background.extend(
            list(neighborhood[:peak_region[0] - pPeakWindowSize, :].flatten()))
        # from peak to bottom
        background.extend(
            list(neighborhood[peak_region[0] + pPeakWindowSize+1:, :].flatten()))
        
        # right middle
        background.extend(
            list(neighborhood[peak_region[0]-pPeakWindowSize:peak_region[0]+pPeakWindowSize+1, peak_region[1]+pPeakWindowSize+1:].flatten()))
        # left middle
        background.extend(
            list(neighborhood[peak_region[0]-pPeakWindowSize:peak_region[0]+pPeakWindowSize+1, :peak_region[1]-pPeakWindowSize].flatten()))
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
        
        # if np.mean(peak) < pMinimumInteractionsThreshold:
        #     mask.append(False)
        #     continue


        
        donut_test_data = []

        horizontal = []

        # top middle
        horizontal.extend(neighborhood[:peak_region[0]-pPeakWindowSize, peak_region[1]-pPeakWindowSize:peak_region[1]+pPeakWindowSize+1].flatten())
        # bottom middle
        horizontal.extend(neighborhood[peak_region[0]+pPeakWindowSize+1:, peak_region[1]-pPeakWindowSize:peak_region[1]+pPeakWindowSize+1].flatten())
        horizontal = np.array(horizontal).flatten()

        vertical = []
        # left
        vertical.extend(neighborhood[peak_region[0]-pPeakWindowSize:peak_region[0]+pPeakWindowSize+1, :peak_region[1]-pPeakWindowSize].flatten())
        # right
        vertical.extend(neighborhood[peak_region[0]-pPeakWindowSize:peak_region[0]+pPeakWindowSize+1, peak_region[1]+pPeakWindowSize+1:].flatten())
        vertical = np.array(vertical).flatten()

        # # top left
        # donut_test_data.append(neighborhood[:peak_region[0] - pPeakWindowSize, :peak_region[1] - pPeakWindowSize].flatten())

        # # top right
        # donut_test_data.append(neighborhood[peak_region[0] + pPeakWindowSize:, :peak_region[1] - pPeakWindowSize].flatten())

        

        # bottom left
        bottom_left_corner = []
        bottom_left_corner.extend(neighborhood[peak_region[0]:, :peak_region[1]-pPeakWindowSize].flatten())
        bottom_left_corner.extend(neighborhood[peak_region[0]+pPeakWindowSize+1:, peak_region[1]-pPeakWindowSize:peak_region[1]+1].flatten())
        # bottom_left_corner = np.array(bottom_left_corner).flatten()
        donut_test_data.append(bottom_left_corner)
        donut_test_data.append(horizontal)
        donut_test_data.append(vertical)


        # # bottom right
        # donut_test_data.append(neighborhood[peak_region[0] + pPeakWindowSize:, peak_region[1] + pPeakWindowSize:].flatten())

        



        
        
        # min_value = np.min(neighborhood)
        # max_value = np.max(neighborhood)
        # min_max_difference = np.float64(max_value - min_value)

        # neighborhood -= min_value
        # neighborhood /= min_max_difference
        # variance = np.var(neighborhood)
        # # log.debug('np.var() {}'.format(variance))
        # neighborhood *= min_max_difference
        # neighborhood += min_value

        # if variance < pMinVariance or variance > pMaxVariance:
        #     mask.append(False)
        #     continue




        # if np.min(background) * 5 > np.max(peak):
        #     mask.append(False)
        #     continue
        del neighborhood
        if pStatisticalTest == 'wilcoxon-rank-sum':

            accept_length = len(donut_test_data)
            accept_count = 0
            for data in donut_test_data:
                statistic, significance_level_test1 = ranksums(sorted(peak), sorted(data))
                if significance_level_test1 <= pPValue:
                    accept_count += 1
                # else:
                #     mask.append(False)
                #     break
            if accept_count >= 3:
                # mask.append(True)
                statistic, significance_level = ranksums(sorted(peak), sorted(background))
                if significance_level <= pPValue:
                #     # log.debug('True')
                    mask.append(True)
                    pvalues.append(significance_level)
                    del peak
                    del background
                    del donut_test_data
                    continue
                else:
                    mask.append(False)
                    del peak
                    del background
                    del donut_test_data
                    continue
            else:
                mask.append(False)
                del peak
                del background
                del donut_test_data
                continue

            # test vertical
            # test_list = get_test_data(neighborhood, pVertical=True)
            # statistic, significance_level_test1 = ranksums(sorted(test_list[0]), sorted(test_list[1]))
            # statistic, significance_level_test2 = ranksums(sorted(test_list[1]), sorted(test_list[2]))
            # if significance_level_test1 <= pPValue or significance_level_test2 <= pPValue:
            #     test_list = get_test_data(neighborhood, pVertical=False)
            #     statistic, significance_level_test1 = ranksums(sorted(test_list[0]), sorted(test_list[1]))
            #     statistic, significance_level_test2 = ranksums(sorted(test_list[1]), sorted(test_list[2]))
            #     if significance_level_test1 <= pPValue or significance_level_test2 <= pPValue:

            #         statistic, significance_level = ranksums(sorted(peak), sorted(background))
            #         if significance_level <= pPValue:
            #             # log.debug('True')
            #             mask.append(True)
            #             pvalues.append(significance_level)
            #             continue
            # log.debug('False')
            
            # mask.append(False)
            # continue
        else:
            # test_list = get_test_data(neighborhood, pVertical=True)
            # _, _, significance_level_test1 = anderson_ksamp([sorted(test_list[0]), sorted(test_list[1])])
            # _, _, significance_level_test2 = anderson_ksamp([sorted(test_list[1]), sorted(test_list[2])])
            # if significance_level_test1 <= pPValue or significance_level_test2 <= pPValue:
            #     test_list = get_test_data(neighborhood, pVertical=False)
            #     _, _, significance_level_test1 = anderson_ksamp([sorted(test_list[0]), sorted(test_list[1])])
            #     _, _, significance_level_test2 = anderson_ksamp([sorted(test_list[1]), sorted(test_list[2])])
            #     if significance_level_test1 <= pPValue or significance_level_test2 <= pPValue:

            #         _, _, significance_level = anderson_ksamp([sorted(peak), sorted(background)])
            #         if significance_level <= pPValue:
            #             # log.debug('True')
            #             mask.append(True)
            #             pvalues.append(significance_level)
            #             continue
            #         # log.debug('False')

            # mask.append(False)
            # continue
            accept_length = len(donut_test_data)
            accept_count = 0
            for data in donut_test_data:
                statistic, significance_level_test1 = anderson_ksamp(sorted(peak), sorted(data))
                if significance_level_test1 <= pPValue:
                    accept_count += 1
                # else:
                #     mask.append(False)
                #     break
            if accept_count >= 6:
                # mask.append(True)
                # statistic, significance_level = ranksums(sorted(peak), sorted(background))
                # if significance_level <= pPValue:
                #     # log.debug('True')
                mask.append(True)
                pvalues.append(significance_level_test1)
                continue
            else:
                mask.append(False)
                continue
        
        # mask.append(False)
    # pvalues = []
    # neighborhood_list = np.array(neighborhood_list)
    # np.save('neighborhood.npy', neighborhood_list)
    # np.save('mask.npy', mask)
    # np.save('candidates.npy', pCandidates)


    pQueue.put([mask, pvalues])
    return


def candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue, pQValue,
                          pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pMinVariance, pMaxVariance, pStatisticalTest=None):
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
    # x_max = pHiCMatrix.shape[0]
    # y_max = pHiCMatrix.shape[1]
    log.debug('candidate_region_test initial: {}'.format(len(pCandidates)))

    if len(pCandidates) == 0:
        return None, None

    pCandidates = np.array(pCandidates)

    mask = []

    mask_thread = [None] * pThreads
    pvalues_thread = [None] * pThreads

    interactionFilesPerThread = len(pCandidates) // pThreads
    all_data_collected = False
    queue = [None] * pThreads
    process = [None] * pThreads
    thread_done = [False] * pThreads
    # length_of_threads = 0
    for i in range(pThreads):

        if i < pThreads - 1:
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
            pMinVariance=pMinVariance, 
            pMaxVariance=pMaxVariance,
            pStatisticalTest=pStatisticalTest

        )
        )
# candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
        #   pMinimumInteractionsThreshold, pPeakWindowSize, pThreads, pStatisticalTest
        if len(candidateThread) == 0:
            process[i] = None
            queue[i] = None
            thread_done[i] = True
            pvalues_thread[i] = []
            mask_thread[i] = []
            continue
        else:
            process[i].start()
    # log.debug('length_of_threads {}'.format(length_of_threads))
    while not all_data_collected:
        for i in range(pThreads):
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
    del mask_thread
    del pvalues_thread
    # if pStatisticalTest == 'anderson-darling':
    if mask is not None and len(mask) == 0:
        return None, None
    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    log.debug('candidate_region_test done: {}'.format(len(pCandidates)))
    pvalues = np.array(pvalues)

    # if pStatisticalTest == 'wilcoxon-rank-sum':
    # log.debug('len pCandidates: {}'.format(len(pCandidates)))
    # log.debug('len pvalues: {}'.format(len(pvalues)))

    # log.debug('Apply fdr')
    # pvalues = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
    # pvalues_ = sorted(pvalues)
    # largest_p_i = 0
    # for i, p in enumerate(pvalues_):
    #     if p <= (pQValue * (i + 1) / len(pvalues_)):
    #         if p >= largest_p_i:
    #             largest_p_i = p
    # mask = pvalues < largest_p_i
    # pCandidates = pCandidates[mask]
    # pvalues = pvalues[mask]

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
    del pCandidates
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


def compute_loops(pHiCMatrix, pRegion, pArgs, pIsCooler, pQueue=None):
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

    # if pIsCooler:
    pHiCMatrix = hm.hiCMatrix(pMatrixFile=pArgs.matrix, pChrnameList=[pRegion])
    if not pIsCooler:
        # pHiCMatrix.setMatrix(
        #     deepcopy(matrix), deepcopy(cut_intervals))
        pHiCMatrix.keepOnlyTheseChr([pRegion])

    if len(pHiCMatrix.matrix.data) == 0:
        pQueue.put([None])
        return

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

    pHiCMatrix.matrix = obs_exp_matrix(pHiCMatrix.matrix)
    pHiCMatrix.matrix.eliminate_zeros()

    # handle pValuePreselection
    try:
        pArgs.pValuePreselection = float(pArgs.pValuePreselection)
    except Exception:
        pArgs.pValuePreselection= read_threshold_file(pArgs.pValuePreselection)

    candidates, pValueList = compute_long_range_contacts(pHiCMatrix,
                                                         pArgs.windowSize,
                                                         pArgs.maximumInteractionPercentageThreshold,
                                                         pArgs.pValue,
                                                         pArgs.qValue,
                                                         pArgs.peakWidth,
                                                         pArgs.pValuePreselection,
                                                         pArgs.statisticalTest,
                                                         pArgs.peakInteractionsThreshold,
                                                         min_loop_distance,
                                                         max_loop_distance,
                                                         pArgs.threadsPerChromosome,
                                                         pMinVariance=pArgs.minVariance,
                                                         pMaxVariance=pArgs.maxVariance)

    if candidates is None:
        log.info('Computed loops for {}: 0'.format(pRegion))

        if pQueue is None:
            return None
        else:
            pQueue.put([None])
            return
    mapped_loops = cluster_to_genome_position_mapping(
        pHiCMatrix, candidates, pValueList, pArgs.maxLoopDistance)
    del pHiCMatrix
    del candidates
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

def read_threshold_file(pFile):
    distance_value_dict = {}
    with open(pFile, 'r') as file:
        file_ = True
        while file_:
            line = file.readline().strip()
            if line.startswith('#'):
                continue
            if line == '':
                break
            relative_distance , value = line.split('\t')
            distance_value_dict[int(relative_distance)] = float(value)
    return distance_value_dict

def main(args=None):
    log.debug('Newc')
    args = parse_arguments().parse_args(args)
    log.info('peak interactions threshold set to {}'.format(
        args.peakInteractionsThreshold))

    if args.region is not None and args.chromosomes is not None:
        log.error('Please choose either --region or --chromosomes.')
        exit(1)
    log.debug('args.matrix {}'.format(args.matrix))
    is_cooler = check_cooler(args.matrix)
    log.debug('is_cooler {}'.format(is_cooler))
    if args.threadsPerChromosome < 1:
        args.threadsPerChromosome = 1
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
                # chromosomes_list = cooler.Cooler(args.matrix).chromnames
                chromosome_sizes = cooler.Cooler(args.matrix).chromsizes
                
                # shuffle the processing order of chromosomes. 
                # with this one large chromosome and 4 smalls are in a row
                # peak memory is reduced and more chromosomes can be processed in parallel on low memory systems.
                sorted_sizes_desc = chromosome_sizes.sort_values(ascending=False)

                size = sorted_sizes_desc.size
                chromosome_names_list = sorted_sizes_desc.index.tolist()
                # threads = 5
                chromosomes_list = []
                i = 0
                j = args.threads # biggest + thread smallest; 2nd biggest chr + 4 - 8 smallest
                k = size - 1
                while i < size:
                    chromosomes_list.append(chromosome_names_list[i])
                    while j > 0 and k > 0:
                        if k == i:
                            break
                        chromosomes_list.append(chromosome_names_list[k])
                        k -= 1
                        j -= 1
                    j = args.threads-1
                    if i == k:
                        break
                    i += 1

                
        else:
            chromosomes_list = args.chromosomes

        if len(chromosomes_list) < args.threads:
            args.threads = len(chromosomes_list)
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
                loops = compute_loops(hic_matrix, chromosome, args, is_cooler)
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

                        # if len(hic_matrix.matrix.data) > 0:

                        process[i] = Process(target=compute_loops, kwargs=dict(
                            pHiCMatrix=args.matrix,
                            pRegion=chromosomes_list[count_call_of_read_input],
                            pArgs=args,
                            pIsCooler=is_cooler,
                            pQueue=queue[i]
                        ))
                        process[i].start()

                        # else:
                        #     queue[i] = None
                        #     thread_done[i] = True
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
