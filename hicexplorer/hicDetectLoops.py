from __future__ import division
import argparse
from sklearn.cluster import dbscan
from sklearn.metrics.pairwise import euclidean_distances
from scipy.sparse import csr_matrix, triu
# from scipy import stats
from scipy.stats import normaltest, f_oneway, mannwhitneyu, zscore, kstest, norm, expon, maxwell
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
    parserOpt.add_argument('--zScoreThreshold', '-zt',
                           type=float,
                           default=5.0,
                           help='z-score threshold to detect long range interactions')
    parserOpt.add_argument('--zScoreMatrixName', '-zs',
                           help='Saves the computed z-score matrix')
    parserOpt.add_argument('--windowSize', '-w',
                           type=int,
                           default=5,
                           help='Window size')

    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.15,
                           help='d-value')
    parserOpt.add_argument('--qValue', '-q',
                           type=float,
                           default=0.05,
                           help='q value for FDR')

    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=int,
                           default=8,
                           help='The minimum number of interactions a detected peaks needs to have to be considered.')
    parserOpt.add_argument('--maxLoopDistance', '-mld',
                           type=int,
                           default=3000000,
                           help='maximum distance of a loop, usually loops are within a distance of ~2MB.')

    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes and order in which the chromosomes should be saved. If not all chromosomes '
                           'are given, the missing chromosomes are left out. For example, --chromosomeOrder chrX will '
                           'export a matrix only containing chromosome X.',
                           nargs='+')

    parserOpt.add_argument('--region',
                           help='The format is chr:start-end.',
                           required=False)
    parserOpt.add_argument('--threads',
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


def compute_zscore_matrix(pMatrix, pInteractionHeight):


    pMatrix.data = pMatrix.data.astype(float)
    # mask = pMatrix.data < pInteractionHeight
    # pMatrix.data[mask] = 0
    # pMatrix.eliminate_zeros()

    instances, features = pMatrix.nonzero()

    data = deepcopy(pMatrix.data)
    distances = np.absolute(instances - features)
    sigma_2 = np.zeros(pMatrix.shape[0])

    # shifts all values by one, but 0 / mean is prevented
    sum_per_distance = np.ones(pMatrix.shape[0])
    np.nan_to_num(data, copy=False)

    sum_per_distance, distance_count = _sum_per_distance(sum_per_distance, data, distances)
    # compute mean
    mean_adjusted = sum_per_distance / distance_count
    # compute sigma squared
    data_mean = pMatrix.data - mean_adjusted[distances]
    data = np.square(data_mean)

    for i in range(len(instances)):
        sigma_2[distances[i]] += data[i]

    sigma_2 /= distance_count
    sigma = np.sqrt(sigma_2)
    sigma /= np.sqrt(distance_count)
    data_mean /= sigma[distances]

    return data_mean


def compute_long_range_contacts(pHiCMatrix, pZscoreData, pZscoreThreshold, pWindowSize, pPValue, pQValue, pPeakInteractionsThreshold):

    # keep only z-score values if they are higher than pThreshold
    # keep: z-score value, (x, y) coordinate

    # keep only z-score values (peaks) if they have more interactions than pPeakInteractionsThreshold

    instances, features = pHiCMatrix.matrix.nonzero()
    log.debug('Number of candidates 147: {}'.format(len(instances)))
    #  instances, features = pHiCMatrix.matrix.nonzero()
    z_score_matrix = csr_matrix((pZscoreData, (instances, features)), shape=(pHiCMatrix.matrix.shape[0], pHiCMatrix.matrix.shape[1]))
        
    # filter by threshold
    mask = pZscoreData >= pZscoreThreshold
    # pZscoreData = None

    mask_interactions = pHiCMatrix.matrix.data > pPeakInteractionsThreshold

    mask = np.logical_and(mask, mask_interactions)

    instances = instances[mask]
    features = features[mask]
    pZscoreData = pZscoreData[mask]
    interactionHeight = pHiCMatrix.matrix.data[mask]
    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])

    log.debug('Number of candidates after z-score and height pruning: {}'.format(len(candidates)))

    number_of_candidates = 0
    i = 0
    while number_of_candidates != len(candidates):
        number_of_candidates = len(candidates)

        candidates, _ = window_zscore_cluster(candidates, pWindowSize, z_score_matrix)
        i +=1 

        log.debug('iterations: {}'.format(i))
        if len(candidates) == 0:
            return None, None

    candidates, pValueList = candidate_peak_exponential_distribution_test(
        z_score_matrix, candidates, pWindowSize, pPValue, pQValue, None)
    if candidates is not None:
        log.debug('Candidates: {}'.format(candidates[:20]))
    return candidates, pValueList

def precluster(pCandidates, pZscore, pWindowSize):
    cluster_candidates = []
    pCandidates = np.array(pCandidates)

    candidateT = pCandidates.T
    mask = np.absolute(candidateT[0] - candidateT[1]) > 4
    pCandidates = pCandidates[mask]
    pZscore = pZscore[mask]
    candidates_per_cluster = 5
    if len(pCandidates) > 100:
        n_bins = (len(pCandidates) // 100) + 1
        x_values = pCandidates.T[0]

        elements, edges = np.histogram(x_values, bins=n_bins)
        j = 0
        for i in range(0, len(elements)):

            _candidates = pCandidates[j:j+elements[i]]
            j += elements[i]
            if len(_candidates) > 0:
                cluster_candidates.append(_candidates)
    else:
        cluster_candidates.append(pCandidates)
    mask = [False] * len(pCandidates)
    z_score_count = 0
    selector = 0
    for j, candidate in enumerate(cluster_candidates):
        
        distances = euclidean_distances(candidate)
        log.debug('distances: {}'.format(distances))
        # # call DBSCAN
        clusters = dbscan(X=distances, eps=1.5,
                        metric='precomputed', min_samples=2, n_jobs=1)[1]
        # log.debug('Size clusters: {}'.format(clusters))                
        cluster_dict = {}
        for i, cluster in enumerate(clusters):
            if cluster == -1:
                # some random high number to not come in conflict with real clusters
                cluster_dict[i+10000] = [[i, pZscore[z_score_count]]]
                z_score_count += 1
                continue
            if cluster in cluster_dict:
                    cluster_dict[cluster].append([i, pZscore[z_score_count]])
            else:
                cluster_dict[cluster] = [[i, pZscore[z_score_count]]]
            z_score_count += 1
        # print(cluster_dict)
        for key in cluster_dict:

            indices = np.argsort(np.array(cluster_dict[key]).T[1])
            for index in indices[:candidates_per_cluster]:

                mask[selector + cluster_dict[key][index][0]] = True
        if len(cluster_candidates) != 1:
            selector += elements[j]
    # Remove values within the window size of other candidates
    mask = np.array(mask)
    pCandidates = np.array(pCandidates)
    # candidates_region = np.array(candidates_region)
    # candidates_region = candidates_region[mask]
    pCandidates = pCandidates[mask]
    pZscore = pZscore[mask]

    return pCandidates, pZscore

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
        if candidate[0] in [81, 82] and candidate[1] in [82, 83, 88]:
            log.debug('x {} y {} max: {}'.format(candidate[0], candidate[1], np.max(neighborhood)))
        # zscore_neighborhood = zscore(neighborhood)
        if len(neighborhood) == 0:
            continue
        argmax = np.argmax(neighborhood)
        # if np.max(zscore_neighborhood) < 3:
            # continue
        x = argmax // (pWindowSize * 2)
        y = argmax % (pWindowSize * 2)
        if candidate[0] in [81, 82]  and candidate[1] in [82, 83, 88]:
            log.debug('x {} y {}'.format(candidate[0], candidate[1]))
            log.debug('pWindowSize {} x {} y {}'.format(pWindowSize, x, y))
            log.debug('x_new {} y_new {} \n'.format((candidate[0] - pWindowSize) + x, (candidate[1] - pWindowSize) + y))



        candidate_x = (candidate[0] - pWindowSize) + x if (candidate[0] - pWindowSize + x) < x_max else x_max - 1
        candidate_y = (candidate[1] - pWindowSize) + y if  (candidate[1] - pWindowSize + y) < y_max else y_max - 1 
        if candidate_x < 0 or candidate_y < 0:
            continue
        new_candidate_list.append([candidate_x, candidate_y])
    # pCandidates = window_zscore_cluster(pCandidates, pWindowSize, pHiCMatrix)
    # log.debug('candidates: {}'.format(new_candidate_list[:50]))
    mask = filter_duplicates(new_candidate_list)

    mask = np.array(mask)
    pCandidates = np.array(new_candidate_list)
    pCandidates = pCandidates[mask]
    return pCandidates, mask

def candidate_peak_exponential_distribution_test(pHiCMatrix, pCandidates, pWindowSize, pPValue, pQValue, pZscore):
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
    log.debug('pCandidates initial: {}'.format(len(pCandidates)))


    if len(pCandidates) == 0:
        return None, None
   
    pCandidates = np.array(pCandidates)
    # for candidate in pCandidates:
    #     if candidate[0] > candidate[1]:
    #         _tmp = candidate[0]
    #         candidate[0] = candidate[1]
    #         candidate[1] = _tmp
    # mask = filter_duplicates(pCandidates)
   
        
    # mask = np.array(mask)
    # pCandidates = pCandidates[mask]
    # pZscore = pZscore[mask]
    mask = []
    # pre-clustering:
    log.debug('pCandidates after duplicate remove: {}'.format(len(pCandidates)))
    # return pCandidates, np.array([1] * len(pCandidates)   ) 

    # pCandidates, pZscore = precluster(pCandidates, pZscore, pWindowSize)

    
    # return pCandidates, np.array([1] * len(pCandidates)   ) 
    mask = []
    # normalized_average_neighborhood = np.zeros((pWindowSize * 2 ) **2)
    neighborhood_list = []
    # pCandidates = window_zscore_cluster(pCandidates, pWindowSize, pHiCMatrix)

    # mask = filter_duplicates(pCandidates)

    # mask = np.array(mask)

    # # Find other candidates within window size
    # # accept only the candidate with the lowest pvalue

    # # remove candidates which
    # #  - do not full fill neighborhood size requirement
    # #  - have a nan pvalue
    # pCandidates = np.array(pCandidates)

    # pCandidates = pCandidates[mask]
    # pZscore = pZscore[mask]
    # for candidate in pCandidates:
    #     if candidate[0] > candidate[1]:
    #         _tmp = candidate[0]
    #         candidate[0] = candidate[1]
    #         candidate[1] = _tmp
    # mask = filter_duplicates(pCandidates)
   
        
    # mask = np.array(mask)
    # pCandidates = pCandidates[mask]
    # # pZscore = pZscore[mask]
    # log.debug('pCandidates after z-score neighborhood: {}'.format(len(pCandidates)))
    mask = []
    for i, candidate in enumerate(pCandidates):

        start_x = candidate[0] - \
            pWindowSize if candidate[0] - pWindowSize > 0 else 0
        end_x = candidate[0] + pWindowSize if candidate[0] + \
            pWindowSize < x_max else x_max
        start_y = candidate[1] - \
            pWindowSize if candidate[1] - pWindowSize > 0 else 0
        end_y = candidate[1] + pWindowSize if candidate[1] + \
            pWindowSize < y_max else y_max

        neighborhood = pHiCMatrix[start_x:end_x,
                                         start_y:end_y].toarray().flatten()
        # zscore_neighborhood = zscore(neighborhood)

        # argmax = np.argmax(zscore_neighborhood)
        # x = argmax // pWindowSize
        # y = argmax % pWindowSize
        # candidate[0] = (candidate[0] - pWindowSize) + x
        # candidate[1] = (candidate[1] - pWindowSize) + y

        mean = np.mean(neighborhood)
        mask_data = neighborhood >= mean
        peak_region = neighborhood[mask_data]
        # if np.absolute(mean - np.max(neighborhood)) < 10 or np.absolute(mean - np.max(neighborhood)) < mean * 2:
        #     mask.append(False)
        #     continue

        if len(peak_region) == 0:
            mask.append(False)
            continue
        loc, scale =expon.fit(peak_region)
        peak_norm = expon(loc=loc, scale=scale)
        

        test_statistic_peak = kstest(peak_region, peak_norm.cdf)#, alternative = 'greater')
        # test_statistic_peak = kstest(peak_region, peak_norm.cdf)
    

        if np.isnan(test_statistic_peak[1]):
            mask.append(False)
            continue
        if test_statistic_peak[1] <= pPValue and test_statistic_peak[1] > 0:
            mask.append(True)
            pvalues.append(test_statistic_peak[1])
        else:
            mask.append(False)
    
    mask = np.array(mask)
    pCandidates = pCandidates[mask]

    log.debug('pCandidates after testing: {}'.format(len(pCandidates)))

    # pZscore = pZscore[mask]
    pvalues = np.array(pvalues)
    if len(pCandidates) == 0:
        return None, None
    # pCandidates, pvalues = precluster(pCandidates, pvalues, pWindowSize)

    # pvalues = pvalues[mask]
    # log.debug('Number of candidates: {}'.format(len(pCandidates)))
    if len(pCandidates) == 0:
        return None, None
  
    mask = []
    

    # for i, candidate in enumerate(pCandidates):
    #     if pvalues[i] <= largest_p_i:
    #         accepted_index.append(candidate)
    #         pvalues_accepted.append(pvalues[i])
    if len(pCandidates) == 0:
        return None, None
    # remove duplicate elements
    # for candidate in pCandidates:
    #     if candidate[0] > candidate[1]:
    #         _tmp = candidate[0]
    #         candidate[0] = candidate[1]
    #         candidate[1] = _tmp
    # duplicate_set_x = set([])
    # duplicate_set_y = set([])

    # delete_index = []
    # for i, candidate in enumerate(pCandidates):
    #     if candidate[0] in duplicate_set_x and candidate[1] in duplicate_set_y:
    #         delete_index.append(i)
    #     else:
    #         duplicate_set_x.add(candidate[0])
    #         duplicate_set_y.add(candidate[1])

    # delete_index = np.array(delete_index)
    # pCandidates = np.array(pCandidates)
    # pvalues = np.array(pvalues)
    # pCandidates = np.delete(pCandidates, delete_index, axis=0)
    # pvalues = np.delete(pvalues, delete_index, axis=0)

    log.debug('pCandidates after duplicate remove: {}'.format(len(pCandidates)))

    # pQValue = 0.01
    
    
    pvalues_ = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
    pvalues_ = np.sort(pvalues_)
    largest_p_i = -np.inf
    for i, p in enumerate(pvalues_):
        if p <= (pQValue * (i + 1) / len(pvalues_)):
            if p >= largest_p_i:
                largest_p_i = p
    # pvalues_accepted = []

    mask = pvalues <= largest_p_i
    pvalues = pvalues[mask]
    pCandidates = pCandidates[mask]

    log.debug('pCandidates after fdr: {}'.format(len(pCandidates)))
    if len(pCandidates) == 0:
        return None, None
    return pCandidates, pvalues


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
    # log.debug("Compute z-score matrix")
    pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
    log.debug('len(instances) after length pruning I {}'.format(len(pHiCMatrix.matrix.data)))
    max_loop_distance = None
    if pArgs.maxLoopDistance:
        max_loop_distance = pArgs.maxLoopDistance / pHiCMatrix.getBinSize()
        instances, features = pHiCMatrix.matrix.nonzero()
        log.debug('len(instances) {}'.format(len(instances)))
        distances = np.absolute(instances - features)
        mask = distances > max_loop_distance
        pHiCMatrix.matrix.data[mask] = 0
        pHiCMatrix.matrix.eliminate_zeros()
        log.debug('len(instances) after length pruning II{}'.format(len(pHiCMatrix.matrix.data)))

    

    z_score_data = compute_zscore_matrix(pHiCMatrix.matrix, pArgs.peakInteractionsThreshold)
    if pArgs.zScoreMatrixName:

        matrixFileHandlerOutput = MatrixFileHandler(pFileType='cool')
        instances, features = pHiCMatrix.matrix.nonzero()
        z_score_matrix = csr_matrix((z_score_data, (instances, features)), shape=(pHiCMatrix.matrix.shape[0], pHiCMatrix.matrix.shape[1]))
        matrixFileHandlerOutput.set_matrix_variables(z_score_matrix, pHiCMatrix.cut_intervals, pHiCMatrix.nan_bins,
                                                     None, pHiCMatrix.distance_counts)
        matrixFileHandlerOutput.save(
            pRegion + '_' + pArgs.zScoreMatrixName, pSymmetric=True, pApplyCorrection=False)

    candidates, pValueList = compute_long_range_contacts(pHiCMatrix, z_score_data, pArgs.zScoreThreshold,
                                                         pArgs.windowSize, pArgs.pValue, pArgs.qValue,
                                                         pArgs.peakInteractionsThreshold)
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


def main():

    args = parse_arguments().parse_args()

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
