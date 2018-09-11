from __future__ import division
import argparse
from sklearn.cluster import dbscan
from sklearn.metrics.pairwise import euclidean_distances
from scipy.sparse import csr_matrix
# from scipy import stats
from scipy.stats import normaltest, f_oneway, mannwhitneyu, zscore
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicexplorer.utilities import check_cooler
from hicexplorer.hicPlotMatrix import translate_region
import numpy as np
import logging
log = logging.getLogger(__name__)


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
                           default=2.0,
                           help='z-score threshold to detect long range interactions')
    parserOpt.add_argument('--scoreMatrixName', '-zs',
                           help='Saves the computed z-score matrix')
    parserOpt.add_argument('--windowSize', '-w',
                           type=int,
                           default=2,
                           help='Window size')

    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.05,
                           help='p value')
    parserOpt.add_argument('--qValue', '-q',
                           type=float,
                           default=0.05,
                           help='q value for FDR')
    parserOpt.add_argument('--minValueRemove', '-r',
                           type=float,
                           default=2,
                           help='p')
    parserOpt.add_argument('--minNeighborhoodSize', '-mns',
                           type=int,
                           default=20,
                           help='p')
    parserOpt.add_argument('--maxLoopDistance', '-ld',
                           type=int,
                           default=5000000,
                           help='maximum distance of a loop')

    parserOpt.add_argument('--chromosomeOrder',
                           help='Chromosomes and order in which the chromosomes should be saved. If not all chromosomes '
                           'are given, the missing chromosomes are left out. For example, --chromosomeOrder chrX will '
                           'export a matrix only containing chromosome X.',
                           nargs='+')

    parserOpt.add_argument('--region',
                           help='The format is chr:start-end. The region should not be longer than 5MB.',
                           required=True)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def _sum_per_distance(pSum_per_distance, pData, pDistances, pN):
    list_of_zero = []
    for i in range(pN):
        if not np.isnan(pData[i]):
            pSum_per_distance[pDistances[i]] += pData[i]
            if pDistances[i] == 0:
                list_of_zero.append(pData[i])
    return pSum_per_distance


def compute_zscore_matrix(pMatrix):

    instances, features = pMatrix.nonzero()

    pMatrix.data = pMatrix.data.astype(float)
    data = pMatrix.data.tolist()
    distances = np.absolute(instances - features)
    sigma_2 = np.zeros(pMatrix.shape[0])
    sum_per_distance = np.zeros(pMatrix.shape[0])
    log.debug('expected distance')
    sum_per_distance = _sum_per_distance(sum_per_distance, data, distances, len(instances))
    # compute mean
    log.debug('mean')

    max_contacts = np.array(range(pMatrix.shape[0], 0, -1))
    mean = sum_per_distance / max_contacts

    # compute sigma squared
    for i in range(len(instances)):
        if np.isnan(data[i]):
            sigma_2[distances[i]] += np.square(mean[distances[i]])
        else:
            sigma_2[distances[i]] += np.square(data[i] - mean[distances[i]])
    log.debug('sigma_2')

    sigma_2 /= max_contacts
    sigma = np.sqrt(sigma_2)
    log.debug('sigma')

    for i in range(len(instances)):
        if np.isnan(data[i]):
            data[i] = (0 - mean[distances[i]]) / sigma[distances[i]]
        else:
            data[i] = (pMatrix.data[i] - mean[distances[i]]) / sigma[distances[i]]

    log.debug('set data')

    return csr_matrix((data, (instances, features)), shape=(pMatrix.shape[0], pMatrix.shape[1])), mean


def compute_long_range_contacts(pHiCMatrix, pZscoreMatrix, pThreshold, pWindowSize, pPValue, pQValue,
                                pMinValueRemove, pMinNeighborhoodSize, pMaxLoopDistance, pExpectedInteractions):

    # keep only z-score values if they are higher than pThreshold
    # keep: z-score value, (x, y) coordinate

    zscore_matrix = pZscoreMatrix
    zscore_matrix.eliminate_zeros()
    instances, features = zscore_matrix.nonzero()
    data = zscore_matrix.data.astype(float)

    log.debug('len(instances) {}, len(features) {}'.format(len(instances), len(features)))
    log.debug('len(data) {}'.format(len(data)))
    # filter by threshold
    mask = data >= pThreshold
    log.debug('len(data) {}, len(mask) {}'.format(len(data), len(mask)))

    _data = data[mask]
    while len(_data) > 30000:
        pThreshold += 0.05
        mask = data >= pThreshold
        _data = data[mask]
        log.debug('Adjusting z-score threshold to: {}'.format(pThreshold))

    data = _data
    instances = instances[mask]
    features = features[mask]
    if len(features) == 0:
        log.info('No loops detected.')
        exit()
    candidates = [*zip(instances, features)]

    candidates, pValueList = candidate_uniform_distribution_test(pHiCMatrix, candidates, pWindowSize, pPValue, pQValue,
                                                                 pMinValueRemove, pMinNeighborhoodSize, data, pExpectedInteractions)
    return cluster_to_genome_position_mapping(pHiCMatrix, candidates, pValueList, pMaxLoopDistance)


def candidate_uniform_distribution_test(pHiCMatrix, pCandidates, pWindowSize, pPValue, pQValue,
                                        pMinValueRemove, pMinNeighborhoodSize, pZscores, pExpectedInteractions):
    # this function test if the values in the neighborhood of a
    # function follow a normal distribution, given the significance level pValue.
    # if they do, the candidate is rejected, otherwise accepted
    # The neighborhood is defined as: the square around a candidate i.e. x-pWindowSize,  :
    accepted_index = []
    mask = []
    pvalues = []
    deleted_index = []
    high_impact_values = []
    x_max = pHiCMatrix.matrix.shape[0]
    y_max = pHiCMatrix.matrix.shape[1]
    maximal_value = 0
    log.info('Candidates: {}'.format(len(pCandidates)))


    log.debug('Filtering neighborhood')
    mask = []

    # distances = euclidean_distances(pCandidates)
    # # # call DBSCAN
    # clusters = dbscan(X=distances, eps=pWindowSize + 1, metric='precomputed', min_samples=2)[1]
    # cluster_dict = {}
    # for i, cluster in enumerate(clusters):
    #     if cluster == -1:
    #         mask.append(True)
    #         continue
    #     if cluster in cluster_dict:
    #         if pZscores[i] > cluster_dict[cluster][1]:
    #             # log.info('cluster: {}, cluster_dict[cluster] {}'.format(cluster, cluster_dict[cluster]))
    #             mask[cluster_dict[cluster][0]] = False
    #             cluster_dict[cluster] = [i, pZscores[i]]
    #             mask.append(True)
    #         else:
    #             mask.append(False)

    #     else:
    #         cluster_dict[cluster] = [i, pZscores[i]]
    #         mask.append(True)

    # # Remove values within the window size of other candidates
    # pCandidates = np.array(pCandidates)

    # mask = np.array(mask)
    # pCandidates = pCandidates[mask]
    # # pvalues = np.array(pvalues)
    # # pvalues = pvalues[mask]
    # mask = []
    # log.debug('Number of candidates: {}'.format(len(pCandidates)))

    log.debug('start testing for significance ')
    for i, candidate in enumerate(pCandidates):

        start_x = candidate[0] - pWindowSize if candidate[0] - pWindowSize > 0 else 0
        end_x = candidate[0] + pWindowSize if candidate[0] + pWindowSize < x_max else x_max
        start_y = candidate[1] - pWindowSize if candidate[1] - pWindowSize > 0 else 0
        end_y = candidate[1] + pWindowSize if candidate[1] + pWindowSize < y_max else y_max

        neighborhood = pHiCMatrix.matrix[start_x:end_x, start_y:end_y].toarray().flatten()

        zscore_neighborhood = zscore(neighborhood) > 0
        np.sum(zscore_neighborhood)
        # log.debug('zscore neighborhood {}'.format(zscore(neighborhood)))
        if np.sum(zscore_neighborhood) < 10:
            mask.append(False)
            # log.debug('No peak detected')
            continue
        if i < 10:
            log.info('candidate: {}'.format(candidate))
            log.info('x_max {}, y_max {}, start_x {}, end_x {}, start_y {}, end_y {}'.format(x_max, y_max, start_x, end_x, start_y, end_y))

        # neighborhood = neighborhood[~np.isnan(neighborhood)]
        # mask_neighborhood = neighborhood > pMinValueRemove
        # neighborhood = neighborhood[mask_neighborhood]
        # if len(neighborhood) < pMinNeighborhoodSize:
        #     mask.append(False)
        #     continue
        
        expected_neighborhood = np.random.uniform(low=np.min(neighborhood), high=np.max(neighborhood), size=len(neighborhood))
        # expectedNeighborhood(pExpectedInteractions, start_x, end_x, #


        # if i < 10:
        #     log.info('neigborhood {}'.format(neighborhood))
        #     log.info('max(neigborhood) {}'.format(np.max(neighborhood)))

        #     log.info('uniform_distribution {}'.format(uniform_distribution))

        # test_result = f_oneway(neighborhood, uniform_distribution)
        neighborhood = np.nan_to_num(neighborhood)
        expected_neighborhood = np.nan_to_num(expected_neighborhood)
        # test_result = mannwhitneyu(neighborhood, expected_neighborhood)
        try:
            test_result = mannwhitneyu(neighborhood, expected_neighborhood)
        except Exception:
            mask.append(False)
            continue
        pvalue = test_result[1]
        # if i < 10:
        #     log.info('pvalue {}'.format(pvalue))
        if np.isnan(pvalue):
            mask.append(False)
            continue
        if pvalue < pPValue:
            mask.append(True)
            pvalues.append(pvalue)
        else:
            mask.append(False)

    log.debug('MW-Test and p-value filtering...DONE')
    mask = np.array(mask)

    # Find other candidates within window size
    # accept only the candidate with the lowest pvalue

    # remove candidates which
    #  - do not full fill neighborhood size requirement
    #  - have a nan pvalue
    pCandidates = np.array(pCandidates)

    pCandidates = pCandidates[mask]
    log.debug('Number of candidates: {}'.format(len(pCandidates)))

    mask = []

    distances = euclidean_distances(pCandidates)
    # # call DBSCAN
    clusters = dbscan(X=distances, eps=pWindowSize, metric='precomputed', min_samples=2, n_jobs=-1)[1]
    cluster_dict = {}
    for i, cluster in enumerate(clusters):
        if cluster == -1:
            mask.append(True)
            continue
        if cluster in cluster_dict:
            if pvalues[i] < cluster_dict[cluster][1]:
                # log.info('cluster: {}, cluster_dict[cluster] {}'.format(cluster, cluster_dict[cluster]))
                mask[cluster_dict[cluster][0]] = False
                cluster_dict[cluster] = [i, pvalues[i]]
                mask.append(True)
            else:
                mask.append(False)

        else:
            cluster_dict[cluster] = [i, pvalues[i]]
            mask.append(True)

    # Remove values within the window size of other candidates
    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    pvalues = np.array(pvalues)
    pvalues = pvalues[mask]
    log.debug('Number of candidates: {}'.format(len(pCandidates)))
    
    # log.debug('Filtering same neighborhood...DONE')

    # FDR
    pvalues_ = np.array([e if ~np.isnan(e) else 1 for e in pvalues])
    pvalues_ = np.sort(pvalues_)
    log.info('pvalues_ {}'.format(pvalues_))
    largest_p_i = -np.inf
    for i, p in enumerate(pvalues_):
        if p <= (pQValue * (i + 1) / len(pvalues_)):
            if p >= largest_p_i:
                largest_p_i = p
    pvalues_accepted = []
    log.info('largest_p_i {}'.format(largest_p_i))
    for i, candidate in enumerate(pCandidates):
        if pvalues[i] < largest_p_i:
            accepted_index.append(candidate)
            pvalues_accepted.append(pvalues[i])
    if len(accepted_index) == 0:
        log.info('No loops detected.')
        exit()
    # remove duplicate elements
    for i, candidate in enumerate(accepted_index):
        if candidate[0] > candidate[1]:
            _tmp = candidate[0]
            candidate[0] = candidate[1]
            candidate[1] = _tmp
    duplicate_set_x = set([])
    duplicate_set_y = set([])

    delete_index = []
    for i, candidate in enumerate(accepted_index):
        if candidate[0] in duplicate_set_x and candidate[1] in duplicate_set_y:
            delete_index.append(i)
        else:
            duplicate_set_x.add(candidate[0])
            duplicate_set_y.add(candidate[1])

    log.debug('FDR filtering...DONE')

    log.debug('delete_index {}'.format(delete_index))
    delete_index = np.array(delete_index)
    accepted_index = np.array(accepted_index)
    pvalues_accepted = np.array(pvalues_accepted)
    accepted_index = np.delete(accepted_index, delete_index, axis=0)
    pvalues_accepted = np.delete(pvalues_accepted, delete_index, axis=0)

    log.debug('number of accepted_index {}'.format(len(accepted_index)))
    # log.info('Accepted candidates: {}'.format(len(accepted_index)))
    # log.info('Accepted candidates: {}'.format(accepted_index))

    return accepted_index, pvalues_accepted


def cluster_to_genome_position_mapping(pHicMatrix, pCandidates, pPValueList, pMaxLoopDistance):
    # mapping: chr_X, start, end, chr_Y, start, end, cluster_id
    mapped_cluster = []
    for i, candidate in enumerate(pCandidates):
        chr_x, start_x, end_x, _ = pHicMatrix.getBinPos(candidate[0])
        chr_y, start_y, end_y, _ = pHicMatrix.getBinPos(candidate[1])
        distance = abs(int(start_x) - int(start_y))
        if distance > pMaxLoopDistance:
            # log.debug('Distance: {}'.format(distance / 1e6))
            continue
        mapped_cluster.append((chr_x, start_x, end_x, chr_y, start_y, end_y, pPValueList[i]))
    return mapped_cluster


def write_bedgraph(pClusters, pOutFileName, pStartRegion, pEndRegion):

    with open(pOutFileName, 'w') as fh:
        for cluster_item in pClusters:
            if cluster_item[1] >= pStartRegion and cluster_item[2] <= pEndRegion \
                    and cluster_item[4] >= pStartRegion and cluster_item[5] <= pEndRegion:
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % cluster_item)

    with open('loops_domains.bed', 'w') as fh:
        for cluster_item in pClusters:
            if cluster_item[1] >= pStartRegion and cluster_item[2] <= pEndRegion \
                    and cluster_item[4] >= pStartRegion and cluster_item[5] <= pEndRegion:
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cluster_item[0], cluster_item[1], cluster_item[4], 1, cluster_item[6],
                                                                       ".", cluster_item[1], cluster_item[4], "x,x,x"))


def main():

    args = parse_arguments().parse_args()

    is_cooler = check_cooler(args.matrix)
    
    # log.info("Cooler or no cooler: {}".format(is_cooler))
    open_cooler_chromosome_order = True
    if args.chromosomeOrder is not None and len(args.chromosomeOrder) > 1:
        open_cooler_chromosome_order = False
    if args.region:
        chrom, region_start, region_end = translate_region(args.region)
    if is_cooler and open_cooler_chromosome_order:
        # log.info("Retrieve data from cooler format and use its benefits.")
        regionsToRetrieve = None
        if args.region:
            regionsToRetrieve = []
            regionsToRetrieve.append(args.region)
        elif args.chromosomeOrder:
            regionsToRetrieve = args.chromosomeOrder

        hic_matrix = hm.hiCMatrix(pMatrixFile=args.matrix, pChrnameList=regionsToRetrieve)
    else:
        hic_matrix = hm.hiCMatrix(args.matrix)
        if args.region:
            # chrom, region_start, region_end = translate_region(args.region)
            hic_matrix.keepOnlyTheseChr([chrom])
        elif args.chromosomeOrder:
            hic_matrix.keepOnlyTheseChr(args.chromosomeOrder)


    hic_matrix.matrix = hic_matrix.matrix + hic_matrix.matrix.transpose() - csr_matrix(np.diag(hic_matrix.matrix.diagonal()))
    log.debug("Compute z-score matrix")
    z_score_matrix, expected_interactions = compute_zscore_matrix(hic_matrix.matrix)
    if args.scoreMatrixName:
        hic_matrix_save = hm.hiCMatrix(pMatrixFile=args.scoreMatrixName)
        hic_matrix_save.setMatrix(z_score_matrix, hic_matrix.cut_intervals)
        hic_matrix_save.save(args.scoreMatrixName)
    log.debug('Compute loops')
    mapped_clusters = compute_long_range_contacts(hic_matrix, z_score_matrix, args.zScoreThreshold,
                                                  args.windowSize, args.pValue, args.qValue, args.minValueRemove,
                                                  args.minNeighborhoodSize, args.maxLoopDistance, expected_interactions)

    # write it to bedgraph
    write_bedgraph(mapped_clusters, args.outFileName, region_start, region_end)
