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
                                'This does NOT influence the p-value for the neighborhood testing.')
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=int,
                           help='The minimum number of interactions a detected peaks needs to have to be considered. The number of interactions decreases with increasing genomic distances: '
                                ' peakInteractionsThreshold/log(genomic distance)')
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           default=0.01,
                           help='Rejection level for Anderson-Darling test for H0. H0 is peak region and background have the same distribution.')

    parserOpt.add_argument('--maxLoopDistance',
                           type=int,
                           default=2000000,
                           help='Maximum genomic distance of a loop, usually loops are within a distance of ~2MB.')
    parserOpt.add_argument('--minLoopDistance',
                           type=int,
                           default=200000,
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
                            default="Wilcoxon rank-sum",
                            choices=['Wilcoxon rank-sum', 'Anderson-Darling']
                           )
                           
                          
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def create_distance_distribution(pData, pDistances):
    pGenomicDistanceDistribution = {}
    pGenomicDistanceDistributionPosition = {}

    for i, distance in enumerate(pDistances):

        if distance in pGenomicDistanceDistribution:
            pGenomicDistanceDistribution[distance].append(pData[i])
            pGenomicDistanceDistributionPosition[distance].append(i)
        else:
            pGenomicDistanceDistribution[distance] = [pData[i]]
            pGenomicDistanceDistributionPosition[distance] = [i]

    return pGenomicDistanceDistribution, pGenomicDistanceDistributionPosition


def compute_long_range_contacts(pHiCMatrix, pWindowSize,
                                pPeakInteractionsThreshold, pPValue, pPeakWindowSize, pPValuePreselection, pStatisticalTest):
    """
        This function computes the loops by:
            - decreasing the search space by removing zScore values < 0
            - decreasing the search space with the zScoreMeanFactor
            - decreasing the search space by excluding candidates with too less counts
            - merging candidates which share a neighborhood to one candidate
            - calling candidate_region_test to test the neighborhood of a candidate against its peak to detect significant peaks

        Input:
            - pHiCMatrix: original interaction matrix for neighborhood and peak region subset selection
            - pZscoreMatrix: z-score matrix
            - pAdjustedHiCMatrix: interation matrix adjusted to dimensions and sparsity of z-score matrix
            - pWindowSize: integer, the size of (2*pWindowSize)^2 around a candidate defines its neighborhood. It is used for
                    a) merging candidates and their neighborhoods to one candidate per neighborhood which appear in this region.
                    b) same neighborhood is used to test the peak region against the neighborhood for significant difference (candidate_region_test)
            - pPeakInteractionsThreshold: integer, remove candidates with less interactions
            - pZscoreMeanFactor: float, per genomic distance prune all z-scores with: z-score < mean(z-scores) * pZscoreMeanFactor
            - pPValue: float, test rejection level for H0 and Bonferroni correction
            - pPeakWindowSize: integer, size of the peak region: (2*pPeakWindowSize)^2. Needs to be smaller than pWindowSize

        Returns:
            - A list of detected loops [(x,y)] and x, y are matrix index values
            - An associated list of p-values
    """
    instances, features = pHiCMatrix.matrix.nonzero()
    distance = np.absolute(instances - features)
    mask = [False] * len(distance)
    genomic_distance_distributions, pGenomicDistanceDistributionPosition = create_distance_distribution(
        pHiCMatrix.matrix.data, distance)
    nbinom_parameters = {}
    for i, key in enumerate(genomic_distance_distributions):
        nbinom_parameters = fit_nbinom.fit(
            np.array(genomic_distance_distributions[key]))
        nbinom_distance = nbinom(
            nbinom_parameters['size'], nbinom_parameters['prob'])
        less_than = np.array(
            genomic_distance_distributions[key]).astype(int) - 1
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
                mask[pGenomicDistanceDistributionPosition[key][j]] = True

    peak_interaction_threshold_array = pPeakInteractionsThreshold / \
        np.log(distance)
    mask_interactions = pHiCMatrix.matrix.data > peak_interaction_threshold_array

    mask = np.logical_and(mask, mask_interactions)

    instances = instances[mask]
    features = features[mask]
    if len(features) == 0:
        return None, None
    candidates = np.array([*zip(instances, features)])

    # Clean neighborhood, results in one candidate per neighborhood
    number_of_candidates = 0
    while number_of_candidates != len(candidates):
        number_of_candidates = len(candidates)

        candidates, mask = neighborhood_merge(
            candidates, pWindowSize, pHiCMatrix.matrix)

        if len(candidates) == 0:
            return None, None

    candidates, p_value_list = candidate_region_test(
        pHiCMatrix.matrix, candidates, pWindowSize, pPValue,
        pPeakInteractionsThreshold, pPeakWindowSize, pStatisticalTest)

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


def neighborhood_merge(pCandidates, pWindowSize, pInteractionCountMatrix):
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
    x_max = pInteractionCountMatrix.shape[0]
    y_max = pInteractionCountMatrix.shape[1]
    new_candidate_list = []

    for candidate in pCandidates:

        # get neighborhood out of pInteractionCountMatrix matrix
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
        if np.absolute(candidate_x - candidate_y) < 4:
            continue
        new_candidate_list.append([candidate_x, candidate_y])
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


def candidate_region_test(pHiCMatrix, pCandidates, pWindowSize, pPValue,
                          pPeakInteractionsThreshold, pPeakWindowSize, pStatisticalTest=None):
    """
        Tests if a candidate is having a significant peak compared to its neighborhood.
            - smoothes neighborhood in x an y orientation
            - remove candidate if smoothed peak value < pPeakInteractionsThreshold
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
            - pPeakInteractionsThreshold: integer, if smoothed candidate interaction count is less, it will be removed
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

        if neighborhood[peak_region[0], peak_region[1]] < pPeakInteractionsThreshold:
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
        
        if pStatisticalTest == 'Wilcoxon rank-sum':
            statistic, significance_level = ranksums(peak, background)
        else:
            statistic, critical_values, significance_level = anderson_ksamp([peak, background])

        if significance_level <= pPValue:
            mask.append(True)
            pvalues.append(significance_level)

            continue

        mask.append(False)

    if mask is not None and len(mask) == 0:
        return None, None
    mask = np.array(mask)
    pCandidates = pCandidates[mask]
    log.debug('candidate_region_test done: {}'.format(len(pCandidates)))
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
            - Computes z-score for this chromosome per genomic distance
            - Calls compute_long_range_contacts
            - Writes computed loops to a bedgraph file

        Input:
            - pHiCMatrix: Hi-C interaction matrix object
            - pRegion: Chromosome name
            - pArgs: Argparser object
            - pQueue: Queue object for multiprocessing communication with parent process
    """
    if pArgs.peakInteractionsThreshold is None:
        max_value = np.max(pHiCMatrix.matrix.data)
        if 0 <= max_value <= 100:
            pArgs.peakInteractionsThreshold = 2
        elif 100 < max_value <= 500:
            pArgs.peakInteractionsThreshold = max_value * 0.01
        elif 500 < max_value <= 1000:
            pArgs.peakInteractionsThreshold = max_value * 0.004
        elif 1000 < max_value <= 10000:
            pArgs.peakInteractionsThreshold = max_value * 0.005
        elif 10000 < max_value <= 100000:
            pArgs.peakInteractionsThreshold = max_value * 0.001
        else:
            pArgs.peakInteractionsThreshold = max_value * 0.0005

        # pArgs.peakInteractionsThreshold = np.ceil(max_value / int(np.log10(max_value) ** 5))
        log.debug('peak interactions threshold set to {}'.format(
            pArgs.peakInteractionsThreshold))
        log.debug('max values {}'.format(
            max_value))
        log.debug("sum: {}".format(np.sum(pHiCMatrix.matrix.data)))
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
    log.debug('Setting values to: {} because default'.format(pArgs.windowSize))
    if pArgs.peakWidth is None:
        pArgs.peakWidth = pArgs.windowSize - 4
    log.debug('Setting peak width to: {}'.format(pArgs.peakWidth))

    pHiCMatrix.matrix = triu(pHiCMatrix.matrix, format='csr')
    pHiCMatrix.matrix.eliminate_zeros()
    log.debug('candidates region {} {}'.format(
        pRegion, len(pHiCMatrix.matrix.data)))
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

    log.debug('candidates region {} min max boundary {}'.format(
        pRegion, len(pHiCMatrix.matrix.data)))

    candidates, pValueList = compute_long_range_contacts(pHiCMatrix,
                                                         pArgs.windowSize,
                                                         pArgs.peakInteractionsThreshold,
                                                         pArgs.pValue,
                                                         pArgs.peakWidth,
                                                         pArgs.pValuePreselection,
                                                         pArgs.statisticalTest)

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
        log.debug('chromosomes_list {}'.format(chromosomes_list))
        log.debug('cooler.Cooler(args.matrix).chromsizes {}'.format(
            cooler.Cooler(args.matrix).extent('Y')))
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
