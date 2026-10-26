"""hicDetectLoops: a thin layer over the C++ executable hicDetectLoops.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import argparse
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Array, RawArray
from copy import deepcopy
import logging
log = logging.getLogger(__name__)
import time
import gc
import cooler
import numpy as np
from scipy.sparse import csr_matrix, triu
from scipy.stats import anderson_ksamp, ranksums
from scipy.stats import nbinom
# import scipy.sparse
import fit_nbinom

from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
from hicexplorer.hicPlotMatrix import translate_region
from hicexplorer.lib import cnb
from inspect import currentframe
import traceback


from hicexplorer.utilities import obs_exp_matrix, obs_exp_matrix_non_zero


def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno


def parse_arguments(args=None):
    """The argparse parser of hicDetectLoops, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDetectLoops')


def create_distance_distribution(pDataObsExp, pDistances, pWindowSize, pMinDistance, pMaxDistance, pQueue):
    pGenomicDistanceDistributionObsExp = {}
    pGenomicDistanceDistributionPosition = {}
    try:
        for distance in range(pMinDistance, pMaxDistance, 1):
            mask = pDistances == distance

            pGenomicDistanceDistributionObsExp[distance] = pDataObsExp[mask]
            pGenomicDistanceDistributionPosition[distance] = np.argwhere(mask == True).flatten()
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put([pGenomicDistanceDistributionPosition, pGenomicDistanceDistributionObsExp])
    return


def compute_p_values_mask(pGenomicDistanceDistributionsObsExp, pGenomicDistanceDistributionsKeyList,
                          pPValuePreselection, pGenomicDistanceDistributionPosition, pResolution,
                          pMinimumInteractionsThreshold, pObsExpThreshold, pQueue):

    try:
        true_values = []
        if len(pGenomicDistanceDistributionsKeyList) == 0:
            pQueue.put(true_values)
            return

        float_dict = isinstance(pPValuePreselection, float)
        for i, key in enumerate(pGenomicDistanceDistributionsKeyList):

            data_obs_exp = np.array(pGenomicDistanceDistributionsObsExp[key])
            # do not fit and not compute any p-value if all values on this distance are small than the pMinimumInteractionsThreshold
            mask = data_obs_exp >= pObsExpThreshold
            nbinom_parameters = fit_nbinom.fit(data_obs_exp)

            p_value = 1 - cnb.cdf(data_obs_exp[mask], nbinom_parameters['size'], nbinom_parameters['prob'])

            if float_dict:
                mask_distance = p_value <= pPValuePreselection
            else:
                key_genomic = int(key * pResolution)
                mask_distance = p_value <= pPValuePreselection[key_genomic]
            j = 0
            for k, value in enumerate(mask):
                if value:
                    if mask_distance[j]:
                        true_values.append(pGenomicDistanceDistributionPosition[key][k])
                    j += 1
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put(true_values)
    return


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
    # pObsExpMatrix.eliminate_zeros()
    # pHiCMatrix.matrix.eliminate_zeros()
    instances, features = pObsExpMatrix.nonzero()
    distance = np.absolute(instances - features)
    mask_interactions_hard_threshold = pHiCMatrix.matrix.data >= pMinimumInteractionsThreshold

    del instances
    del features
    queue = [None] * pThreads
    process = [None] * pThreads
    genomic_distance_distributions_thread = None
    genomic_distance_distributions_obs_exp_thread = None

    pGenomicDistanceDistributionPosition_thread = None
    all_data_collected = False
    thread_done = [False] * pThreads
    min_distance = distance.min()
    max_distance = distance.max()
    len_distance = len(distance)

    distances_per_threads = (max_distance - min_distance) // pThreads
    for i in range(pThreads):

        if i < pThreads - 1:
            min_distance_thread = min_distance + (i * distances_per_threads)
            max_distance_thread = min_distance + ((i + 1) * distances_per_threads)
        else:
            min_distance_thread = min_distance + (i * distances_per_threads)
            max_distance_thread = max_distance + 1
        queue[i] = Queue()
        process[i] = Process(target=create_distance_distribution, kwargs=dict(
            pDataObsExp=pObsExpMatrix.data,
            pDistances=distance,
            pWindowSize=pWindowSize,
            pMinDistance=min_distance_thread,
            pMaxDistance=max_distance_thread,
            pQueue=queue[i]
        )
        )

        process[i].start()

    del pHiCMatrix.matrix
    del distance
    # genomic_distance_distributions = {}
    genomic_distance_distributions_obs_exp = {}
    pGenomicDistanceDistributionPosition = {}
    fail_flag = False
    fail_message = ''
    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                queue_data = queue[i].get()
                if 'Fail:' in queue_data:
                    fail_flag = True
                    fail_message = queue_data
                else:
                    pGenomicDistanceDistributionPosition_thread, \
                        genomic_distance_distributions_obs_exp_thread = queue_data

                    pGenomicDistanceDistributionPosition = {**pGenomicDistanceDistributionPosition, **pGenomicDistanceDistributionPosition_thread}
                    del pGenomicDistanceDistributionPosition_thread

                    genomic_distance_distributions_obs_exp = {**genomic_distance_distributions_obs_exp, **genomic_distance_distributions_obs_exp_thread}
                    del genomic_distance_distributions_obs_exp_thread

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

    if fail_flag:
        return fail_message, None
    mask = [False] * len_distance
    genomic_distance_distributions_thread = (len(genomic_distance_distributions_obs_exp) // pThreads) + 1

    queue = [None] * pThreads
    process = [None] * pThreads
    all_data_collected = False
    thread_done = [False] * pThreads
    genomic_keys_list = sorted(list(genomic_distance_distributions_obs_exp.keys()))

    for i in range(pThreads):

        if i < pThreads - 1:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:(i + 1) * genomic_distance_distributions_thread]
        else:
            genomic_distance_keys_thread = genomic_keys_list[i * genomic_distance_distributions_thread:]
        if len(genomic_distance_keys_thread) == 0:

            process[i] = None
            queue[i] = None
            thread_done[i] = True
            continue
        else:
            queue[i] = Queue()
            process[i] = Process(target=compute_p_values_mask, kwargs=dict(
                pGenomicDistanceDistributionsObsExp=genomic_distance_distributions_obs_exp,
                pGenomicDistanceDistributionsKeyList=genomic_distance_keys_thread,
                pPValuePreselection=pPValuePreselection,
                pGenomicDistanceDistributionPosition=pGenomicDistanceDistributionPosition,
                pResolution=pHiCMatrix.getBinSize(),
                pMinimumInteractionsThreshold=pMinimumInteractionsThreshold,
                pObsExpThreshold=pObsExpThreshold,
                pQueue=queue[i]
            )
            )

            process[i].start()
            del genomic_distance_keys_thread

    del genomic_distance_distributions_obs_exp
    del pGenomicDistanceDistributionPosition
    del genomic_keys_list

    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                mask_threads = queue[i].get()
                if isinstance(mask_threads, str) and 'Fail: ' in mask_threads:
                    fail_flag = True
                    fail_message = mask_threads
                else:
                    for index in mask_threads:
                        mask[index] = True
                    del mask_threads
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

    if fail_flag:
        return fail_message, None
    mask = np.logical_and(mask, mask_interactions_hard_threshold)
    instances, features = pObsExpMatrix.nonzero()

    instances = instances[mask]
    features = features[mask]

    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])

    del instances
    del features
    del mask

    candidates, _ = neighborhood_merge(
        candidates, pWindowSize, pObsExpMatrix, pThreads)

    if isinstance(candidates, str) and 'Fail: ' in candidates:
        return candidates, None
    if len(candidates) == 0:
        return None, None
    candidates, p_value_list = candidate_region_test(
        pObsExpMatrix, candidates, pWindowSize, pPValue,
        pPeakWindowSize, pThreads)

    return candidates, p_value_list


def neighborhood_merge_thread(pCandidateList, pWindowSize, pInteractionCountMatrix, pQueue):

    try:
        new_candidate_list = []

        if len(pCandidateList) == 0:
            pQueue.put(new_candidate_list)
            return
        x_max = pInteractionCountMatrix.shape[0]
        y_max = pInteractionCountMatrix.shape[1]

        for candidate in pCandidateList:

            if (candidate[0] - pWindowSize) > 0:
                start_x = candidate[0] - pWindowSize
            else:
                start_x = 0

            if (candidate[1] - pWindowSize) > 0:
                start_y = candidate[1] - pWindowSize
            else:
                start_y = 0

            end_x = candidate[0] + pWindowSize + 1 if candidate[0] + \
                pWindowSize + 1 < x_max else x_max

            end_y = candidate[1] + pWindowSize + 1 if candidate[1] + \
                pWindowSize + 1 < y_max else y_max

            neighborhood = pInteractionCountMatrix[start_x:end_x,
                                                   start_y:end_y].toarray().flatten()

            if len(neighborhood) > 0 and np.max(neighborhood) == pInteractionCountMatrix[candidate[0], candidate[1]]:
                new_candidate_list.append(candidate)

            del neighborhood
        del pCandidateList
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put(new_candidate_list)
    return


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
    log.debug('pCandidates {}'.format(pCandidates[:10]))
    new_candidate_list = []

    new_candidate_list_threads = [None] * pThreads
    interactionFilesPerThread = len(pCandidates) // pThreads
    all_data_collected = False
    queue = [None] * pThreads
    process = [None] * pThreads
    thread_done = [False] * pThreads
    for i in range(pThreads):

        if i < pThreads - 1:
            candidateThread = pCandidates[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            candidateThread = pCandidates[i * interactionFilesPerThread:]
        queue[i] = Queue()
        process[i] = Process(target=neighborhood_merge_thread, kwargs=dict(
            pCandidateList=candidateThread,
            pWindowSize=pWindowSize,
            pInteractionCountMatrix=pInteractionCountMatrix,
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
        del candidateThread
    del pInteractionCountMatrix
    fail_flag = False
    fail_message = ''
    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                new_candidate_list_threads[i] = queue[i].get()
                if isinstance(new_candidate_list_threads[i], str) and 'Fail: ' in new_candidate_list_threads[i]:
                    fail_flag = True
                    fail_message = new_candidate_list_threads[i]
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
    if fail_flag:
        return fail_message
    new_candidate_list = [item for sublist in new_candidate_list_threads for item in sublist]
    del new_candidate_list_threads
    return new_candidate_list, True


def candidate_region_test_thread(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                                 pPeakWindowSize, pQueue):
    try:
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

            end_x = candidate[0] + pWindowSize + 1 if candidate[0] + \
                pWindowSize + 1 < x_max else x_max

            end_y = candidate[1] + pWindowSize + 1 if candidate[1] + \
                pWindowSize + 1 < y_max else y_max

            neighborhood = pHiCMatrix[start_x:end_x,
                                      start_y:end_y].toarray()
            if len(neighborhood) == 0:
                mask.append(False)
                del neighborhood
                continue
            # get index of original candidate
            peak_region = np.array(np.where(neighborhood == pHiCMatrix[candidate[0], candidate[1]])).flatten()

            if pPeakWindowSize > pWindowSize:
                log.warning('Neighborhood window size ({}) needs to be larger than peak width({}).'.format(
                    pWindowSize, pPeakWindowSize))
                return None, None

            peak = neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize + 1,
                                peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize + 1].flatten()

            background = []
            # top to peak
            background.extend(
                list(neighborhood[:peak_region[0] - pPeakWindowSize, :].flatten()))
            # from peak to bottom
            background.extend(
                list(neighborhood[peak_region[0] + pPeakWindowSize + 1:, :].flatten()))

            # right middle
            background.extend(
                list(neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize + 1, peak_region[1] + pPeakWindowSize + 1:].flatten()))
            # left middle
            background.extend(
                list(neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize + 1, :peak_region[1] - pPeakWindowSize].flatten()))
            background = np.array(background)

            if len(background) < pWindowSize:
                mask.append(False)
                del peak
                del background
                continue
            if len(peak) < pWindowSize:
                mask.append(False)
                del peak
                del background
                continue
            if np.mean(peak) < np.mean(background):
                mask.append(False)
                del peak
                del background
                continue
            if np.max(peak) < np.max(background):
                mask.append(False)
                del peak
                del background
                continue

            donut_test_data = []
            horizontal = []
            # top middle
            horizontal.extend(neighborhood[:peak_region[0] - pPeakWindowSize, peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize + 1].flatten())
            # bottom middle
            horizontal.extend(neighborhood[peak_region[0] + pPeakWindowSize + 1:, peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize + 1].flatten())
            horizontal = np.array(horizontal).flatten()

            vertical = []
            # left
            vertical.extend(neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize + 1, :peak_region[1] - pPeakWindowSize].flatten())
            # right
            vertical.extend(neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize + 1, peak_region[1] + pPeakWindowSize + 1:].flatten())
            vertical = np.array(vertical).flatten()

            # bottom left
            bottom_left_corner = []
            bottom_left_corner.extend(neighborhood[peak_region[0]:, :peak_region[1] - pPeakWindowSize].flatten())
            bottom_left_corner.extend(neighborhood[peak_region[0] + pPeakWindowSize + 1:, peak_region[1] - pPeakWindowSize:peak_region[1] + 1].flatten())
            # bottom_left_corner = np.array(bottom_left_corner).flatten()
            donut_test_data.append(bottom_left_corner)
            donut_test_data.append(horizontal)
            donut_test_data.append(vertical)

            del neighborhood

            # test vertical, horizontal, bottom left corner and neighborhood vs peak with wilcoxon-rank-sum test
            accept_count = 0
            for data in donut_test_data:
                statistic, significance_level_test1 = ranksums(sorted(peak), sorted(data))
                if significance_level_test1 <= pPValue:
                    accept_count += 1
            if accept_count >= 3:
                statistic, significance_level = ranksums(sorted(peak), sorted(background))
                if significance_level <= pPValue:
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
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    del pHiCMatrix
    del pCandidates
    pQueue.put([mask, pvalues])
    return


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

    mask = []
    pvalues = []
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
    for i in range(pThreads):

        if i < pThreads - 1:
            candidateThread = pCandidates[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            candidateThread = pCandidates[i * interactionFilesPerThread:]
        queue[i] = Queue()
        process[i] = Process(target=candidate_region_test_thread, kwargs=dict(
            pHiCMatrix=pHiCMatrix,
            pCandidates=candidateThread,
            pWindowSize=pWindowSize,
            pPValue=pPValue,
            pPeakWindowSize=pPeakWindowSize,
            pQueue=queue[i]

        )
        )
        if len(candidateThread) == 0:
            process[i] = None
            queue[i] = None
            thread_done[i] = True
            pvalues_thread[i] = []
            mask_thread[i] = []
            continue
        else:
            process[i].start()

    fail_flag = False
    fail_message = ''
    while not all_data_collected:
        for i in range(pThreads):
            if queue[i] is not None and not queue[i].empty():
                result_thread = queue[i].get()
                if isinstance(result_thread, str) and 'Fail: ' in result_thread:
                    fail_flag = True
                    fail_message = result_thread
                else:
                    mask_thread[i], pvalues_thread[i] = result_thread
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

    if fail_flag:
        return fail_message, None
    mask = [item for sublist in mask_thread for item in sublist]
    pvalues = [item for sublist in pvalues_thread for item in sublist]
    del mask_thread
    del pvalues_thread
    if mask is not None and len(mask) == 0:
        return None, None
    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    pvalues = np.array(pvalues)

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
    try:

        if pQueue is not None:
            if pIsCooler:
                pHiCMatrix = hm.hiCMatrix(pMatrixFile=pArgs.matrix, pChrnameList=[pRegion], pDistance=pArgs.maxLoopDistance, pNoIntervalTree=True, pUpperTriangleOnly=False)
            else:
                pHiCMatrix = hm.hiCMatrix(pMatrixFile=pArgs.matrix, pChrnameList=[pRegion], pDistance=pArgs.maxLoopDistance, pNoIntervalTree=False, pUpperTriangleOnly=False)

        if not pIsCooler:
            # cooler files load only what is necessary.
            pHiCMatrix.keepOnlyTheseChr([pRegion])
            max_loop_distance = pArgs.maxLoopDistance / pHiCMatrix.getBinSize()
            instances, features = pHiCMatrix.matrix.nonzero()
            distances = np.absolute(instances - features)
            mask = distances > max_loop_distance
            pHiCMatrix.matrix.data[mask] = 0
            pHiCMatrix.matrix.eliminate_zeros()

        if len(pHiCMatrix.matrix.data) == 0:
            pQueue.put([None])
            return

        if pHiCMatrix.matrix.shape[0] < 5 or pHiCMatrix.matrix.shape[1] < 5:
            log.debug('Computed loops for {}: 0'.format(pRegion))

            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        if pArgs.windowSize is None:
            bin_size = pHiCMatrix.getBinSize()
            if 0 < bin_size <= 5000:
                pArgs.windowSize = 10
            elif 5000 < bin_size <= 10000:
                pArgs.windowSize = 5
            elif 10000 < bin_size <= 25000:
                pArgs.windowSize = 5
            elif 25000 < bin_size <= 50000:
                pArgs.windowSize = 5
            else:
                pArgs.windowSize = 5
            log.debug('Setting window size to: {}'.format(pArgs.windowSize))
        if pArgs.peakWidth is None:
            pArgs.peakWidth = pArgs.windowSize - 3
        log.debug('Setting peak width to: {}'.format(pArgs.peakWidth))
        pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
        pHiCMatrix.matrix.eliminate_zeros()
        # log.debug('candidates region {} {}'.format(
        #     pRegion, len(pHiCMatrix.matrix.data)))

        # delete main diagonal
        instances, features = pHiCMatrix.matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances == 0
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()

        del instances
        del features
        del mask
        del distances

        if pArgs.expected == 'mean':
            obs_exp_csr_matrix = obs_exp_matrix(pHiCMatrix.matrix, pInplace=False, pToEpsilon=True, pThreads=pArgs.threadsPerChromosome)
        elif pArgs.expected == 'mean_nonzero':
            obs_exp_csr_matrix = obs_exp_matrix_non_zero(pHiCMatrix.matrix, ligation_factor=False, pInplace=False, pToEpsilon=True, pThreads=pArgs.threadsPerChromosome)

        elif pArgs.expected == 'mean_nonzero_ligation':
            obs_exp_csr_matrix = obs_exp_matrix_non_zero(pHiCMatrix.matrix, ligation_factor=True, pInplace=False, pToEpsilon=True, pThreads=pArgs.threadsPerChromosome)

        if not isinstance(obs_exp_csr_matrix, csr_matrix):
            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        pHiCMatrix.matrix.eliminate_zeros()
        obs_exp_csr_matrix.eliminate_zeros()
        if len(pHiCMatrix.matrix.data) != len(obs_exp_csr_matrix.data):
            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        # handle pValuePreselection
        try:
            pArgs.pValuePreselection = float(pArgs.pValuePreselection)
        except Exception:
            pArgs.pValuePreselection = read_threshold_file(pArgs.pValuePreselection)

        candidates, pValueList = compute_long_range_contacts(pHiCMatrix,
                                                             obs_exp_csr_matrix,
                                                             pArgs.windowSize,
                                                             pArgs.pValue,
                                                             pArgs.peakWidth,
                                                             pArgs.pValuePreselection,
                                                             pArgs.peakInteractionsThreshold,
                                                             pArgs.obsExpThreshold,
                                                             pArgs.threadsPerChromosome)

        if candidates is None:
            log.info('Computed loops for {}: 0'.format(pRegion))
            if pQueue is None:
                return None
            else:
                pQueue.put([None])
                return
        elif isinstance(candidates, str) and 'Fail: ' in candidates and pQueue is not None:
            pQueue.put(candidates)
            return
        elif isinstance(candidates, str) and 'Fail: ' in candidates and pQueue is None:
            return candidates
        mapped_loops = cluster_to_genome_position_mapping(
            pHiCMatrix, candidates, pValueList, pArgs.maxLoopDistance)
        del pHiCMatrix
        del candidates
        log.debug('Computed loops for {}: {}'.format(pRegion, len(mapped_loops)))
    except Exception as exp:
        if pQueue is not None:
            pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
            return
        else:
            return 'Fail: ' + str(exp) + traceback.format_exc()
    if pQueue is None:
        return mapped_loops
    else:
        pQueue.put([mapped_loops])
    return


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
            relative_distance, value = line.split('\t')
            distance_value_dict[int(relative_distance)] = float(value)
    return distance_value_dict


def main(args=None):
    """Run the C++ implementation of hicDetectLoops with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDetectLoops')(args)
