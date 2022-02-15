import pandas as pd
import numpy as np
from scipy.stats import ranksums
import argparse
from multiprocessing import Process, Queue
import time
import traceback
from copy import deepcopy
from tempfile import NamedTemporaryFile, mkdtemp

import logging
log = logging.getLogger(__name__)

from pygenometracks import plotTracks
from hicmatrix import HiCMatrix as hm

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Computes differential TADs by comparing the precomputed TAD regions of the target matrix with the same regions of the control matrix.
Please notice that the matrices need to have the same read coverage, this can be achieved with hicNormalize and the 'smallest'-mode.
H0 is the assumption that two regions are identical, the rejected files contain therefore the as differential considered regions.""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--targetMatrix', '-tm',
                                help='The matrix which was used to compute the TADs',
                                required=False)
    parserRequired.add_argument('--controlMatrix', '-cm',
                                help='The control matrix to test the TADs for a differential interaction pattern.',
                                required=False)
    parserRequired.add_argument('--tadDomains', '-td',
                                help='The TADs domain file computed by hicFindTADs.',
                                required=False)
    parserRequired.add_argument('--outFileNamePrefix', '-o',
                                help='Outfile name prefix to store the accepted / rejected H0 TADs.',
                                required=False,
                                default='output_differential_tad'
                                )
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           help='H0 is considered as \'two regions are identical.\' i.e. all regions with a test result of <= p-value are rejected and considered as differential'
                           ' (Default: %(default)s).',
                           default=0.05)
    parserOpt.add_argument('--mode', '-m',
                           help='Consider only intra-TAD interactions, or additional left inter-TAD, right inter-TAD or all'
                           ' (Default: %(default)s).',
                           choices=['intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'],
                           default='all')
    parserOpt.add_argument('--modeReject', '-mr',
                           help='All test of a mode must be rejected (all) or reject region (and accept it is differential) as soon as at least one region is having a p-value <= --pValue'
                           ' (Default: %(default)s).',
                           choices=['all', 'one'],
                           default='one')
    parserOpt.add_argument('--initFilePGT', '-i',
                           type=str,
                           help='Define a pyGenomeTracks ini file to plot the differential detected TADs.'
                           ' (Default: %(default)s).',
                           default=None)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads to use, the parallelization is implemented per chromosome'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2, 3, 4, 5]]

    return domains_df


def computeDifferentialTADs(pMatrixTarget, pMatrixControl, pDomainList, pCoolOrH5, pPValue, pThreadId, pQueue):
    try:
        accepted_inter_left = []
        accepted_inter_right = []
        accepted_intra = []
        p_values_list = []
        stats_list = []
        rows = []

        chromosome_list = pDomainList
        for i, row in enumerate(chromosome_list):

            if pThreadId is None:
                log.debug('first thread')
                if i == len(chromosome_list) - 1:
                    continue
            elif pThreadId == True:
                log.debug('middle thread')

                if i == 0 or i == len(chromosome_list) - 1:
                    log.debug('i: {}'.format(i))
                    log.debug('len(chromosome_list): {}'.format(len(chromosome_list)))

                    continue
            elif pThreadId == False:
                log.debug('last thread')

                if i == 0:
                    continue

            if i - 1 >= 0:
                chromosom = chromosome_list[i - 1][0]
                start = chromosome_list[i - 1][1]
            else:
                chromosom = chromosome_list[i][0]
                start = chromosome_list[i][1]
            if i + 1 < len(chromosome_list):
                end = chromosome_list[i + 1][2]
            else:
                end = chromosome_list[i][2]
            # midpos = row[1] + ((row[2] - row[1]) / 2)

            if pCoolOrH5:

                # # get intra-TAD data
                hic_matrix_target = hm.hiCMatrix(
                    pMatrixFile=pMatrixTarget, pChrnameList=[str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])])
                hic_matrix_control = hm.hiCMatrix(
                    pMatrixFile=pMatrixControl, pChrnameList=[str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])])
                matrix_target = hic_matrix_target.matrix.toarray()
                matrix_control = hic_matrix_control.matrix.toarray()

                hic_matrix_target_inter_tad = hm.hiCMatrix(
                    pMatrixFile=pMatrixTarget, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])
                hic_matrix_control_inter_tad = hm.hiCMatrix(
                    pMatrixFile=pMatrixControl, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])

                matrix_target_inter_tad = hic_matrix_target_inter_tad.matrix
                matrix_control_inter_tad = hic_matrix_control_inter_tad.matrix

            else:
                # in case of h5 pMatrixTarget is already a HiCMatrix object
                hic_matrix_target = pMatrixTarget
                hic_matrix_control = pMatrixControl
                hic_matrix_target_inter_tad = pMatrixTarget
                hic_matrix_control_inter_tad = pMatrixControl
                indices_target = hic_matrix_target.getRegionBinRange(str(row[0]), row[1], row[2])
                indices_control = hic_matrix_control.getRegionBinRange(str(row[0]), row[1], row[2])

                matrix_target = hic_matrix_target.matrix[indices_target[0]:indices_target[1], indices_target[0]:indices_target[1]].toarray()
                matrix_control = hic_matrix_control.matrix[indices_control[0]:indices_control[1], indices_control[0]:indices_control[1]].toarray()
                matrix_target_inter_tad = pMatrixTarget.matrix
                matrix_control_inter_tad = pMatrixControl.matrix

            matrix_target = matrix_target.flatten()
            matrix_control = matrix_control.flatten()
            # tad_midpoint = hic_matrix_target_inter_tad.getRegionBinRange(str(row[0]), midpos, midpos)[0]

            # if i - 1 >= 0:
            # get index position left tad with tad
            left_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])[0]
            left_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])[0]
            if pCoolOrH5:
                outer_left_boundary_index_target = 0
                outer_left_boundary_index_control = 0

                outer_right_boundary_index_control = -1
                outer_right_boundary_index_target = -1

            else:
                outer_left_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]
                outer_left_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

                outer_right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]
                outer_right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]

            if i + 1 < len(chromosome_list) and not pCoolOrH5:
                # get index position right tad with tad
                right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
                right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
            elif i + 1 < len(chromosome_list):
                right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
                right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]

            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
                intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
                intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()
                intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
                intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
                intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()

            significance_level_left = None
            significance_level_right = None
            statistic_left = None
            statistic_right = None

            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                intertad_left_target = intertad_left_target.flatten()
                intertad_left_control = intertad_left_control.flatten()
                intertad_right_target = intertad_right_target.flatten()
                intertad_right_control = intertad_right_control.flatten()

                statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)
                statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                intertad_right_target = intertad_right_target.flatten()
                intertad_right_control = intertad_right_control.flatten()
                statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                intertad_left_target = intertad_left_target.flatten()
                intertad_left_control = intertad_left_control.flatten()
                # log.debug('intertad_left_target {}'.format(intertad_left_target))
                # log.debug('intertad_left_control {}'.format(intertad_left_control))

                statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)

            # log.debug('matrix_target {}'.format(matrix_target))
            # log.debug('matrix_control {}'.format(matrix_control))

            statistic, significance_level = ranksums(matrix_target, matrix_control)
            # log.debug('statistic {}, significance_level {}'.format(statistic, significance_level))
            # log.debug('right statistic {}, significance_level {}'.format(statistic_right, significance_level_right))
            # log.debug('left statistic {}, significance_level {}'.format(statistic_left, significance_level_left))

            p_values = []
            stats = []
            if significance_level_left is None or np.isnan(significance_level_left):
                accepted_inter_left.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level_left <= pPValue:
                accepted_inter_left.append(1)
                p_values.append(significance_level_left)
                stats.append(statistic_left)
            else:
                accepted_inter_left.append(0)
                p_values.append(significance_level_left)
                stats.append(statistic_left)

            if significance_level_right is None or np.isnan(significance_level_right):
                accepted_inter_right.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level_right <= pPValue:
                accepted_inter_right.append(1)
                p_values.append(significance_level_right)
                stats.append(statistic_right)
            else:
                accepted_inter_right.append(0)
                p_values.append(significance_level_right)
                stats.append(statistic_right)

            if significance_level is None or np.isnan(significance_level):
                accepted_intra.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level <= pPValue:
                accepted_intra.append(1)
                p_values.append(significance_level)
                stats.append(statistic)
            else:
                accepted_intra.append(0)
                p_values.append(significance_level)
                stats.append(statistic)

            p_values_list.append(p_values)
            stats_list.append(stats)

            rows.append(row)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    # hic_matrix_target_inter_tad.save('manipulated_target.cool')
    # hic_matrix_control_inter_tad.save('manipulated_control.cool')
    pQueue.put([stats_list, p_values_list, accepted_inter_left, accepted_inter_right, accepted_intra, rows])


def main(args=None):
    args = parse_arguments().parse_args(args)

    # read domains file
    domains_df = readDomainBoundaries(args.tadDomains)
    log.debug('len(domains_df) {}'.format(len(domains_df)))
    domains = domains_df.values.tolist()
    old_chromosome = None

    tads_per_chromosome = []

    for j in range(len(domains)):
        if old_chromosome is None:
            old_chromosome = domains[j][0]
            per_chromosome = []
            per_chromosome.append(domains[j])

        elif old_chromosome == domains[j][0]:
            per_chromosome.append(domains[j])
            continue
        else:
            tads_per_chromosome.append(per_chromosome)
            per_chromosome = []
            per_chromosome.append(domains[j])
            old_chromosome = domains[j][0]
    tads_per_chromosome.append(per_chromosome)
    # log.debug('len(tads_per_chromosome) {}'.format(len(tads_per_chromosome[0]) + len(tads_per_chromosome[1])))

    # read full h5 or only region if cooler
    is_cooler_target = check_cooler(args.targetMatrix)
    is_cooler_control = check_cooler(args.controlMatrix)

    if is_cooler_target != is_cooler_control:
        log.error('Matrices are not given in the same format!')
        exit(1)
    if not is_cooler_control:
        hic_matrix_target = hm.hiCMatrix(args.targetMatrix)
        hic_matrix_control = hm.hiCMatrix(args.controlMatrix)
    else:
        hic_matrix_target = args.targetMatrix
        hic_matrix_control = args.controlMatrix
    # accepted_H0 = []
    # rejected_H0 = []
    # log.debug('domains_df {}'.format(domains_df))

    stats_chromosomes = []
    p_values_chromosomes = []
    accepted_inter_left_chromosomes = []
    accepted_inter_right_chromosomes = []
    accepted_intra_chromosomes = []
    rows_chromosomes = []

    stats_threads = [[]] * args.threads
    p_values_threads = [[]] * args.threads
    accepted_left_inter_threads = [[]] * args.threads
    accepted_right_inter_threads = [[]] * args.threads
    accepted_intra_threads = [[]] * args.threads
    rows_threads = [[]] * args.threads

    threads_save = deepcopy(args.threads)
    for chromosome in tads_per_chromosome:
        log.debug('tads_per_chromosome {}'.format(chromosome))
        domainsPerThread = len(chromosome) // args.threads
        if domainsPerThread == 0 and len(chromosome) > 0:
            domainsPerThread = 1
            args.threads = 1
        elif domainsPerThread > 0:
            args.threads = threads_save

        all_data_collected = False
        queue = [None] * args.threads
        process = [None] * args.threads
        thread_done = [False] * args.threads
        # None --> first thread, process first element in list, ignore last one
        # True --> middle thread: ignore first and last element in tad processing
        # False --> last thread: ignore first element, process last one
        thread_id = None
        for i in range(args.threads):

            if args.threads == 1:
                domainListThread = chromosome

            elif i == 0:
                domainListThread = chromosome[i * domainsPerThread:((i + 1) * domainsPerThread) + 1]
                thread_id = None
            elif i < args.threads - 1:
                domainListThread = chromosome[(i * domainsPerThread) - 1:((i + 1) * domainsPerThread) + 1]
                thread_id = True

            else:
                domainListThread = chromosome[(i * domainsPerThread) - 1:]
                thread_id = False

            if args.threads == 1:
                thread_id = ''

            log.debug('len(domainListThread) {}'.format(len(domainListThread)))
            log.debug('len(thread_id) {}'.format(thread_id))

            queue[i] = Queue()
            process[i] = Process(target=computeDifferentialTADs, kwargs=dict(
                pMatrixTarget=hic_matrix_target,
                pMatrixControl=hic_matrix_control,
                pDomainList=domainListThread,
                pCoolOrH5=is_cooler_control,
                pPValue=args.pValue,
                pThreadId=thread_id,
                pQueue=queue[i]
            )
            )

            process[i].start()
        fail_flag = False
        fail_message = ''
        while not all_data_collected:
            for i in range(args.threads):

                if queue[i] is not None and not queue[i].empty():
                    queue_data = queue[i].get()
                    if 'Fail:' in queue_data:
                        fail_flag = True
                        fail_message = queue_data
                    else:
                        stats_threads[i], p_values_threads[i], accepted_left_inter_threads[i], \
                            accepted_right_inter_threads[i], \
                            accepted_intra_threads[i], rows_threads[i] = queue_data

                    queue[i] = None
                    process[i].join()
                    process[i].terminate()
                    process[i] = None
                    thread_done[i] = True
                # elif queue[i] is None and

            all_data_collected = True
            for thread in thread_done:
                if not thread:
                    all_data_collected = False
            time.sleep(1)

        # outfile_names = [item for sublist in outfile_names for item in sublist]
        # target_list_name = [
        #     item for sublist in target_list_name for item in sublist]
        if fail_flag:
            log.error(fail_message[6:])
            exit(1)
        stats_chromosomes.append([item for sublist in stats_threads for item in sublist])
        p_values_chromosomes.append([item for sublist in p_values_threads for item in sublist])
        accepted_inter_left_chromosomes.append([item for sublist in accepted_left_inter_threads for item in sublist])
        accepted_inter_right_chromosomes.append([item for sublist in accepted_right_inter_threads for item in sublist])
        accepted_intra_chromosomes.append([item for sublist in accepted_intra_threads for item in sublist])
        rows_chromosomes.append([item for sublist in rows_threads for item in sublist])

        log.debug('rows_threads {}'.format(rows_threads))

    stats_list = [item for sublist in stats_chromosomes for item in sublist]
    p_values_list = [item for sublist in p_values_chromosomes for item in sublist]
    accepted_inter_left = [item for sublist in accepted_inter_left_chromosomes for item in sublist]
    accepted_inter_right = [item for sublist in accepted_inter_right_chromosomes for item in sublist]
    accepted_intra = [item for sublist in accepted_intra_chromosomes for item in sublist]
    rows = [item for sublist in rows_chromosomes for item in sublist]

    stats_list = np.array(stats_list)
    p_values_list = np.array(p_values_list)
    accepted_inter_left = np.array(accepted_inter_left)
    accepted_inter_right = np.array(accepted_inter_right)
    accepted_intra = np.array(accepted_intra)
    rows = np.array(rows)

    if args.mode == 'intra-TAD':
        mask = np.array(accepted_intra, dtype=bool)
    elif args.mode == 'left-inter-TAD':
        if args.modeReject == 'all':
            mask = np.logical_and(accepted_inter_left, accepted_intra)
        else:
            mask = np.logical_or(accepted_inter_left, accepted_intra)

    elif args.mode == 'right-inter-TAD':
        if args.modeReject == 'all':
            mask = np.logical_and(accepted_intra, accepted_inter_right)
        else:
            mask = np.logical_or(accepted_intra, accepted_inter_right)

    else:
        if args.modeReject == 'all':
            mask = np.logical_and(accepted_inter_left, accepted_inter_right)
            mask = np.logical_and(mask, accepted_intra)
        else:
            mask = np.logical_or(accepted_inter_left, accepted_inter_right)
            mask = np.logical_or(mask, accepted_intra)

    log.debug('len(mask) {}'.format(len(mask)))
    log.debug('mask.sum() {}'.format(mask.sum()))
    log.debug('mask[:10] {}'.format(mask[:10]))

    accepted_H0 = p_values_list[~mask]
    accepted_H0_s = stats_list[~mask]
    rejected_H0 = p_values_list[mask]
    rejected_H0_s = stats_list[mask]
    accepted_rows = rows[~mask]
    rejected_rows = rows[mask]
    with open(args.outFileNamePrefix + '_accepted.diff_tad', 'w') as file:
        header = '# Created with HiCExplorer\'s hicDifferentialTAD version ' + __version__ + '\n'
        header += '# H0 \'regions are equal\' H0 is accepted for all p-value greater the user given p-value threshold; i.e. regions in this file are not considered as differential.\n'
        header += '# Accepted regions with Wilcoxon rank-sum test to p-value: {}  with used mode: {} and modeReject: {} \n'.format(args.pValue, args.mode, args.modeReject)
        header += '# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD\tp-value right-inter-TAD\tp-value intra-TAD\tW left-inter-TAD\tW right-inter-TAD\tW intra-TAD\n'
        file.write(header)
        for i, row in enumerate(accepted_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\t')
            pvalue_list = list(map(str, accepted_H0[i]))
            file.write('\t'.join(pvalue_list))
            file.write('\t')
            stats_list = list(map(str, accepted_H0_s[i]))
            file.write('\t'.join(stats_list))

            file.write('\n')
    with open(args.outFileNamePrefix + '_rejected.diff_tad', 'w') as file:
        header = '# Created with HiCExplorer\'s hicDifferentialTAD version ' + __version__ + '\n'
        header += '# H0 \'regions are equal\' H0 is rejected for all p-value smaller or equal the user given p-value threshold; i.e. regions in this file are considered as differential.\n'
        header += '# Rejected regions with Wilcoxon rank-sum test to p-value: {}  with used mode: {} and modeReject: {} \n'.format(args.pValue, args.mode, args.modeReject)
        header += '# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD\tp-value right-inter-TAD\tp-value intra-TAD\tW left-inter-TAD\tW right-inter-TAD\tW intra-TAD\n'

        file.write(header)

        for i, row in enumerate(rejected_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\t')
            pvalue_list = list(map(str, rejected_H0[i]))
            file.write('\t'.join(pvalue_list))
            file.write('\t')
            stats_list = list(map(str, rejected_H0_s[i]))
            file.write('\t'.join(stats_list))
            file.write('\n')

    with open(args.outFileNamePrefix + '_accepted_raw.bed', 'w') as file:
        for i, row in enumerate(accepted_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\n')
    with open(args.outFileNamePrefix + '_rejected_raw.bed', 'w') as file:
        for i, row in enumerate(rejected_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\n')

    
    # plot differential regions wiht pygenometracks. 
    # search in 10 MB region ranges if they contain a rejected region
    log.debug('rejected_rows[:][:3] {}'.format(rejected_rows[:, :3]))
    hicmatrix = hm.hiCMatrix()
    plot_regions_intervaltree = hicmatrix.intervalListToIntervalTree(rejected_rows[:, :3])[0]
    
    log.debug('plot_regions_intervaltree {}'.format(plot_regions_intervaltree))

    plot_interval_size = 10000000
    start = 0
    plot_file_format = "pdf"
    # end = plot_interval_size
    # per chromosome

    pyGenomeTracksIniFile = NamedTemporaryFile(suffix='.ini', delete=False)
    pyGenomeTracksIniFile.close()
 
    tracks = args.initFilePGT
    outfile_prefix = args.outFileNamePrefix

    for chromosome in plot_regions_intervaltree:
        end = sorted(list(plot_regions_intervaltree[chromosome]))[-1][1]
        for i in range(start, end , plot_interval_size):
            plot_region_set = plot_regions_intervaltree[chromosome].overlap(i, i +plot_interval_size)
            plot_region_list = sorted(list(plot_region_set))
            

            log.debug("plot_region_list {}".format(plot_region_list))
            if len(plot_region_list) > 0:
                first_tad_size = plot_region_list[0][1] - plot_region_list[0][0]
                last_tad_size = plot_region_list[-1][1] - plot_region_list[-1][0]
                # --tracks all_tads.ini --region chr6:45000000-70000000 -o differential_tad.pdf
                start_point_plot = plot_region_list[0][0] - first_tad_size if plot_region_list[0][0] - first_tad_size < i else i 
                end_point_plot = + plot_region_list[-1][1] + last_tad_size if plot_region_list[-1][1]  + last_tad_size > i+plot_interval_size else i + plot_interval_size
                args = '--tracks {} --region {}:{}-{} -o {}_differential_{}_{}_{}.png'.format(tracks, chromosome, start_point_plot, end_point_plot, outfile_prefix, chromosome, plot_region_list[0][0], plot_region_list[-1][1], plot_file_format).split()
                plotTracks.main(args)
        # end = start + plot_interval_size
            # log.debug("envelop {}".format(plot_regions_intervaltree[chromosome].envelop(i, i +plot_interval_size )))
                # log.debug("overlap {}".format(plot_regions_intervaltree[chromosome].overlap(i, i +plot_interval_size)))

        # start = end
        
        # if start in plot_regions_intervaltree[chromosome]
        # first_plot_region = plot_regions_intervaltree[chromosome][0]
        # last_plot_region = plot_regions_intervaltree[chromosome][-1]

        # log.debug('first_plot_region {}'.format(first_plot_region))
        # log.debug('last_plot_region {}'.format(last_plot_region))
