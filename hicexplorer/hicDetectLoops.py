from __future__ import division
import argparse
from sklearn.cluster import dbscan
from sklearn.metrics.pairwise import euclidean_distances
from scipy.sparse import csr_matrix, triu
# from scipy import stats
from scipy.stats import normaltest, f_oneway, mannwhitneyu, zscore, kstest, norm, expon, maxwell, chisquare
from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicexplorer.utilities import check_cooler
from hicexplorer.hicPlotMatrix import translate_region
import numpy as np
import logging
log = logging.getLogger(__name__)
from copy import deepcopy
import cooler
from multiprocessing import Process, Queue
import time
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_norm, obs_exp_matrix, obs_exp_matrix_non_zero

from scipy.ndimage.filters import gaussian_filter
from scipy.stats import poisson, multivariate_normal, ttest_ind, mannwhitneyu, anderson_ksamp

import matplotlib.pyplot as plt
import matplotlib.patches as patches
def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Computes long range contacts within the given contact matrix
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to compute the long range contacts on.',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                                help='Outfile name with long range contacts clusters, should be a bedgraph.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--peakWidth', '-pw',
                           type=int,
                           default=6,
                           help='The width the peak region should have in bins. The square around the peak will include (2 * peakWidth)^2 bins.')
    parserOpt.add_argument('--zScoreMatrixName', '-zs',
                           help='Saves the computed z-score matrix')
    parserOpt.add_argument('--windowSize', '-w',
                           type=int,
                           default=10,
                           help='The window size for the neighborhood region the peak is located in. All values from this region (exluded the values from the peak '
                           ' region) are tested against the peak region for significant difference. The square will have the size of (2 * windowSize)^2 bins')
    parserOpt.add_argument('--dynamicZScoreThreshold', '-d',
                           type=float,
                           default=1.75,
                           help='Only candidates with z-scores greater mean(z-scores) * dynamicZScoreThreshold are accepted. With '
                                'the dynamic z-score threshold it is possible to decrease or increase this filter criteria.')
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=int,
                           default=20,
                           help='The minimum number of interactions a detected peaks needs to have to be considered.')
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.001,
                           help='Accepance level for t-test for H1.')

    
    parserOpt.add_argument('--maxLoopDistance',
                           type=int,
                           default=1500000,
                           help='maximum distance of a loop, usually loops are within a distance of ~2MB.')
    parserOpt.add_argument('--minLoopDistance',
                           type=int,
                           default=100000,
                           help='Minimum distance of a loop.')
    
    
    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes and order in which the chromosomes should be saved. If not all chromosomes '
                           'are given, the missing chromosomes are left out. For example, --chromosomeOrder chrX will '
                           'export a matrix only containing chromosome X.',
                           nargs='+')
    
    parserOpt.add_argument('--region',
                           help='The format is chr:start-end.',
                           required=False)
    parserOpt.add_argument('--threads', '-t',
                           help='The chromosomes are processed independent of each other. If the number of to be '
                           'processed chromosomes is greater than one, each chromosome will be computed with on its own core.',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def _sum_per_distance(pSum_per_distance, pData, pDistances):
    distance_count = np.zeros(len(pSum_per_distance))
    for i, distance in enumerate(pDistances):
        pSum_per_distance[distance] += pData[i]
        distance_count[distance] += 1
    return pSum_per_distance, distance_count


def compute_zscore_matrix(pInstances, pFeatures, pData, pMaxDistance):

    pData = pData.astype(float)

    data = deepcopy(pData)
    distances = np.absolute(pInstances - pFeatures)
    sigma_2 = np.zeros(pMaxDistance)

    # shifts all values by one, but 0 / mean is prevented
    sum_per_distance = np.ones(pMaxDistance)
    np.nan_to_num(data, copy=False)

    sum_per_distance, distance_count = _sum_per_distance(sum_per_distance, data, distances)
    # compute mean
    mean_adjusted = sum_per_distance / distance_count
    # compute sigma squared
    data_mean = pData - mean_adjusted[distances]
    data = np.square(data_mean)

    for i in range(len(pInstances)):
        sigma_2[distances[i]] += data[i]

    sigma_2 /= distance_count
    sigma = np.sqrt(sigma_2)
    sigma /= np.sqrt(distance_count)
    data_mean /= sigma[distances]

    return data_mean


def compute_long_range_contacts(pHiCMatrix, pZscoreMatrix, pAdjustedHiCMatrix, pWindowSize,
                                pPeakInteractionsThreshold, pZscoreMeanFactor, pPValue, pPeakWindowSize):

    # keep only z-score values if they are higher than pThreshold
    # keep: z-score value, (x, y) coordinate

    # keep only z-score values (peaks) if they have more interactions than pPeakInteractionsThreshold

    instances, features = pAdjustedHiCMatrix.nonzero()

    # filter by distance
    distance = np.absolute(instances - features)
    # mask_distance =  distance >= 10

    # filter by threshold

    mask = pZscoreMatrix.data >= 0
    # mask = np.logical_and(mask, mask_distance)

    zscore_mean = np.mean(pZscoreMatrix.data[mask])
    mask = pZscoreMatrix.data >= (zscore_mean * pZscoreMeanFactor)

    peak_interaction_threshold_array = pPeakInteractionsThreshold / np.log(distance)
    mask_interactions = pAdjustedHiCMatrix.data > peak_interaction_threshold_array

    mask = np.logical_and(mask, mask_interactions)
    # mask = np.logical_and(mask, mask_distance)

    instances = instances[mask]
    features = features[mask]
    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])

    log.debug('Number of candidates after z-score and height pruning: {}'.format(len(candidates)))

    number_of_candidates = 0
    i = 0
    while number_of_candidates != len(candidates):
        number_of_candidates = len(candidates)

        candidates, mask = window_zscore_cluster(candidates, pWindowSize, pZscoreMatrix)
        i += 1

        if len(candidates) == 0:
            return None, None

    candidates, p_value_list = candidate_region_test(
        pHiCMatrix.matrix, candidates, pWindowSize, pPValue,
        # pZscoreMatrix, candidates, pWindowSize, pPValue,

        pPeakInteractionsThreshold, pPeakWindowSize)
    # if candidates is not None:
    #     log.debug('Candidates: {}'.format(candidates[:20]))
    return candidates, p_value_list


def filter_duplicates(pCandidates):

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


def window_zscore_cluster(pCandidates, pWindowSize, pZScoreMatrix):
    x_max = pZScoreMatrix.shape[0]
    y_max = pZScoreMatrix.shape[1]
    new_candidate_list = []

    for candidate in pCandidates:

        start_x = candidate[0] - \
            pWindowSize if candidate[0] - pWindowSize > 0 else 0
        end_x = candidate[0] + pWindowSize if candidate[0] + \
            pWindowSize < x_max else x_max
        start_y = candidate[1] - \
            pWindowSize if candidate[1] - pWindowSize > 0 else 0
        end_y = candidate[1] + pWindowSize if candidate[1] + \
            pWindowSize < y_max else y_max

        neighborhood = pZScoreMatrix[start_x:end_x, start_y:end_y].toarray().flatten()
        if len(neighborhood) == 0:
            continue
        argmax = np.argmax(neighborhood)
        x = argmax // (pWindowSize * 2)
        y = argmax % (pWindowSize * 2)

        candidate_x = (candidate[0] - pWindowSize) + x if (candidate[0] - pWindowSize + x) < x_max else x_max - 1
        candidate_y = (candidate[1] - pWindowSize) + y if (candidate[1] - pWindowSize + y) < y_max else y_max - 1
        if candidate_x < 0 or candidate_y < 0:
            continue
        if np.absolute(candidate_x - candidate_y) < 4:
            continue
        new_candidate_list.append([candidate_x, candidate_y])
    mask = filter_duplicates(new_candidate_list)

    mask = np.array(mask)
    pCandidates = np.array(new_candidate_list)
    pCandidates = pCandidates[mask]
    return pCandidates, mask


def candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                          pPeakInteractionsThreshold, pPeakWindowSize):
    # this function test if the values in the neighborhood of a
    # function follow a normal distribution, given the significance level pValue.
    # if they do, the candidate is rejected, otherwise accepted
    # The neighborhood is defined as: the square around a candidate i.e. x-pWindowSize,  :
    accepted_index = []
    mask = []
    pvalues = []
    deleted_index = []
    high_impact_values = []
    x_max = pHiCMatrix.shape[0]
    y_max = pHiCMatrix.shape[1]
    maximal_value = 0
    log.debug('candidate_region_test initial: {}'.format(len(pCandidates)))

    if len(pCandidates) == 0:
        return None, None

    pCandidates = np.array(pCandidates)

    mask = []

    neighborhood_list = []

    mask = []
    for i, candidate in enumerate(pCandidates):

        if (candidate[0] - pWindowSize) > 0:
            start_x = candidate[0] - pWindowSize
            peak_x = pWindowSize - 1
        else:
            start_x = 0
            peak_x = pWindowSize - candidate[0] -1 
        
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
            neighborhood_old[:, j] = smoothInteractionValues(neighborhood[:, j], 5)
        neighborhood = (neighborhood + neighborhood_old) / 2

        peak_region = [peak_x, peak_y]

        if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold:
            mask.append(False)
            continue
        
        if pPeakWindowSize > pWindowSize:
            log.warning('Neighborhood window size ({}) needs to be larger than peak width({}).'.format(pWindowSize, pPeakWindowSize))
            return None, None
      
        peak = neighborhood[peak_region[0] - pPeakWindowSize:peak_region[0] + pPeakWindowSize, peak_region[1] - pPeakWindowSize:peak_region[1] + pPeakWindowSize].flatten()

        background = []
        background.extend(list(neighborhood[:peak_region[0] - pPeakWindowSize, :].flatten()))
        background.extend(list(neighborhood[peak_region[0] + pPeakWindowSize:, :].flatten()))
        background.extend(list(neighborhood[:, :peak_region[1] - pPeakWindowSize].flatten()))
        background.extend(list(neighborhood[:, peak_region[1] + pPeakWindowSize:].flatten()))
        background = np.array(background)

        # np.savetxt('candidates/candidates {} neighbor'.format(i), np.sort(background))
        # np.savetxt('candidates/candidates {} peak'.format(i), np.sort(peak))
        # # fig, ax = plt.subplots()

        # fig, axis = plt.subplots()
        # heatmap = axis.pcolor(neighborhood, cmap='hot') 
        # plt.colorbar(heatmap)
        # # plt.imshow(neighborhood, cmap='hot', interpolation='nearest')
        # axis.plot([peak_region[0]], [peak_region[1]], 's', lw=2, markerfacecolor='none', markeredgecolor='red')
        # # axis.plot([m:n], [k:l], 's', lw=2, markerfacecolor='none', markeredgecolor='green')
        # rect = patches.Rectangle((peak_region[0] - m, peak_region[1] - k),6,6,linewidth=1,edgecolor='g',facecolor='none')

        # axis.add_patch(rect)
        # plt.savefig('candidates/candidates {}.png'.format(i))

        # plt.close()
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
   
        statistic, pvalue = mannwhitneyu(peak, background)

        if pvalue <= pPValue:
            mask.append(True)
            pvalues.append(pvalue)
            
            continue

        mask.append(False)

    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    log.debug('candidate_region_test done: {}'.format(len(pCandidates)))
    pvalues = np.array(pvalues)

    ###Bonferroni

    if len(pvalues) > 0:
        adjusted_pvalue = pPValue / len(pvalues)

        mask = pvalues <= adjusted_pvalue
        pvalues = pvalues[mask]
        pCandidates = pCandidates[mask]


    ###FDR
    # pvalues = np.array(pvalues)

    # pvalues_ = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
    # pvalues_ = np.sort(pvalues_)
    # largest_p_i = -np.inf
    # for i, p in enumerate(pvalues_):
    #     if p <= (pQValue * (i + 1) / len(pvalues_)):
    #         if p >= largest_p_i:
    #             largest_p_i = p
    # # pvalues_accepted = []

    # mask = pvalues <= largest_p_i
    # pvalues = pvalues[mask]
    # pCandidates = pCandidates[mask]

    log.debug('pCandidates after Bonferroni: {}'.format(len(pCandidates)))
    if len(pCandidates) == 0:
        return None, None

    return pCandidates, pvalues


def smoothInteractionValues(pData, pWindowSize):
    '''
    Adds -pWindowsSize/2 and +pWindowsSize/2 around pData[i] and averages pData[i] by pWindowSize to 
    smooth the interaction values.
    '''
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


def cluster_to_genome_position_mapping(pHicMatrix, pCandidates, pPValueList, pMaxLoopDistance):
    # mapping: chr_X, start, end, chr_Y, start, end, cluster_id
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

    with open(pOutFileName, 'w') as fh:
        for loop_item in pLoops:
            if pStartRegion and pEndRegion:
                if loop_item[1] >= pStartRegion and loop_item[2] <= pEndRegion \
                        and loop_item[4] >= pStartRegion and loop_item[5] <= pEndRegion:
                    fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % loop_item)
            else:
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % loop_item)

    with open('loops_domains.bed', 'w') as fh:
        for loop_item in pLoops:
            if pStartRegion and pEndRegion:

                if loop_item[1] >= pStartRegion and loop_item[2] <= pEndRegion \
                        and loop_item[4] >= pStartRegion and loop_item[5] <= pEndRegion:
                    fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(loop_item[0], loop_item[1], loop_item[4], 1, loop_item[6],
                                                                           ".", loop_item[1], loop_item[4], "x,x,x"))
            else:
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(loop_item[0], loop_item[1], loop_item[4], 1, loop_item[6],
                                                                       ".", loop_item[1], loop_item[4], "x,x,x"))


def compute_loops(pHiCMatrix, pRegion, pArgs, pQueue=None):

    pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
    pHiCMatrix.matrix.eliminate_zeros()

    # s
    max_loop_distance = None
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
        instances, features = pHiCMatrix.matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances > max_loop_distance
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()
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
        instances, features = pHiCMatrix.matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances < min_loop_distance
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()

    pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
    instances, features = deepcopy(pHiCMatrix.matrix.nonzero())
    data_hic = deepcopy(pHiCMatrix.matrix.data)
    if len(pHiCMatrix.matrix.data) == 0:
        log.info('Computed loops for {}: 0'.format(pRegion))

        if pQueue is None:
            return None
        else:
            pQueue.put([None])
            return
    # if pArgs.method == 'lieberman':
    #     obs_exp_norm_matrix = obs_exp_matrix_lieberman(deepcopy(pHiCMatrix.matrix), pHiCMatrix.matrix.shape[0], 19)
    # elif pArgs.method == 'norm':
    #     obs_exp_norm_matrix = obs_exp_matrix_norm(pHiCMatrix.matrix)

    # mask = obs_exp_norm_matrix.data < pArgs.obsExpMinThreshold
    # obs_exp_norm_matrix.data[mask] = 0

    # mask = obs_exp_norm_matrix.data > 0
    # instances = instances[mask]
    # features = features[mask]
    # data_hic = data_hic[mask]

    # obs_exp_norm_matrix.eliminate_zeros()
    # instances_obs_exp, features_obs_exp = obs_exp_norm_matrix.nonzero()
    # data_obs_exp = obs_exp_norm_matrix.data

    z_score_data = compute_zscore_matrix(instances, features, data_hic, pHiCMatrix.matrix.shape[0])
    # z_score_data = compute_zscore_matrix(instances_obs_exp, features_obs_exp, data_obs_exp, pHiCMatrix.matrix.shape[0])

    # z_score_matrix = csr_matrix((z_score_data, (instances_obs_exp, features_obs_exp)), shape=(pHiCMatrix.matrix.shape[0], pHiCMatrix.matrix.shape[1]))
    z_score_matrix = csr_matrix((z_score_data, (instances, features)), shape=(pHiCMatrix.matrix.shape[0], pHiCMatrix.matrix.shape[1]))
    
    if pArgs.zScoreMatrixName:
        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool')

        matrixFileHandlerOutput.set_matrix_variables(z_score_matrix, pHiCMatrix.cut_intervals, pHiCMatrix.nan_bins,
                                                     None, pHiCMatrix.distance_counts)
        matrixFileHandlerOutput.save(
            pRegion + '_' + pArgs.zScoreMatrixName, pSymmetric=True, pApplyCorrection=False)
    # exit()
    adjusted_hic_matrix = csr_matrix((data_hic, (instances, features)), shape=(pHiCMatrix.matrix.shape[0], pHiCMatrix.matrix.shape[1]))
    candidates, pValueList = compute_long_range_contacts(pHiCMatrix, z_score_matrix, adjusted_hic_matrix,
                                                         pArgs.windowSize,
                                                         pArgs.peakInteractionsThreshold,
                                                         pArgs.dynamicZScoreThreshold,
                                                         pArgs.pValue,
                                                         pArgs.peakWidth)

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

    if args.region is not None and args.chromosomes is not None:
        log.error('Please choose either --region or --chromosomes.')
        exit(1)
    is_cooler = check_cooler(args.matrix)

    if args.region:
        chrom, region_start, region_end = translate_region(args.region)

        if is_cooler:
            hic_matrix = hm.hiCMatrix(
                pMatrixFile=args.matrix, pChrnameList=[args.region])
        else:
            hic_matrix = hm.hiCMatrix(args.matrix)
            hic_matrix.keepOnlyTheseChr([chrom])
        mapped_loops = compute_loops(hic_matrix, args.region, args)
        write_bedgraph(mapped_loops, args.outFileName, region_start, region_end)

    else:
        mapped_loops = []

        if not is_cooler:
            hic_matrix = hm.hiCMatrix(args.matrix)
            hic_matrix.keepOnlyTheseChr([chromosome])
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
                    hic_matrix = hm.hiCMatrix(pMatrixFile=args.matrix, pChrnameList=[chromosome])
                else:
                    hic_matrix.setMatrix(deepcopy(matrix), deepcopy(cut_intervals))
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
                            hic_matrix = hm.hiCMatrix(pMatrixFile=args.matrix, pChrnameList=[chromosomes_list[count_call_of_read_input]])
                        else:
                            hic_matrix.setMatrix(deepcopy(matrix), deepcopy(cut_intervals))
                            hic_matrix.keepOnlyTheseChr([chromosomes_list[count_call_of_read_input]])

                        process[i] = Process(target=compute_loops, kwargs=dict(
                            pHiCMatrix=hic_matrix,
                            pRegion=chromosomes_list[count_call_of_read_input],
                            pArgs=args,
                            pQueue=queue[i]
                        ))
                        process[i].start()
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

        if len(mapped_loops) > 0:
            write_bedgraph(mapped_loops, args.outFileName)

    log.info("Number of detected loops for all regions: {}".format(len(mapped_loops)))
