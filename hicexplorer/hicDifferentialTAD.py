# get tad border (domains) file
# get matrix 1 and matrix 2
# get square out of matrix 1 and matrix 2 as defined by border files
# apply statistical test with H0 that they are the same.

import pandas as pd
import numpy as np
from scipy.stats import ranksums
import argparse
from multiprocessing import Process, Queue
import time

import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__


def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno


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
                                required=False)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           help='H0 is considered as \'two regions are identical.\' i.e. all regions with a test result of <= p-value are rejected and considered as differential.',
                           default=0.05)
    parserOpt.add_argument('--mode', '-m',
                           help='Consider only intra-TAD interactions, or additional left inter-TAD, right inter-TAD or all.',
                           choices=['intra-TAD', 'left-inter-TAD', 'right-inter-TAD', 'all'],
                           default='all')
    parserOpt.add_argument('--modeReject', '-mr',
                           help='All test of a mode must be rejected (all) or reject region (and accept it is differential) as soon as at least one region is having a p-value <= --pValue.',
                           choices=['all', 'one'],
                           default='one')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads to use, the parallelization is implemented per chromosome.',
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
    accepted_inter_left = []
    accepted_inter_right = []
    accepted_intra = []
    p_values_list = []
    rows = []
    row2 = ['chr1']
    # hic_matrix_target = hm.hiCMatrix(
    #             pMatrixFile=pMatrixTarget, pChrnameList=['chr1'])
    # hic_matrix_control = hm.hiCMatrix(
    #     pMatrixFile=pMatrixControl, pChrnameList=['chr1'])
    # matrix_target = hic_matrix_target.matrix.toarray()
    # matrix_control = hic_matrix_control.matrix.toarray()

    # hic_matrix_target_inter_tad = hm.hiCMatrix(
    #     pMatrixFile=pMatrixTarget, pChrnameList=['chr1'])
    # hic_matrix_control_inter_tad = hm.hiCMatrix(
    #     pMatrixFile=pMatrixControl, pChrnameList=['chr1'])
    for i, row in enumerate(pDomainList):
        # for domain in domains_df:

        # log.debug('domain {} {} {}'.format(row[0], row[1], row[2]))
        # get inter-tad data

        # None --> first thread, process first element in list, ignore last one
        # True --> middle thread: ignore first and last element in tad processing
        # False --> last thread: ignore first element, process last one
        if pThreadId == None:
            log.debug('first thread')
            if i == len(pDomainList)-1:
                continue
        elif pThreadId == True:
            log.debug('middle thread')

            if i == 0 or i == len(pDomainList)-1:
                continue
        elif pThreadId == False:
            log.debug('last thread')

            if i == 0:
                continue
        
        if i - 1 >= 0:
            chromosom = pDomainList[i - 1][0]
            start = pDomainList[i - 1][1]
        else:
            chromosom = pDomainList[i][0]
            start = pDomainList[i][1]
        if i + 1 < len(pDomainList):
            end = pDomainList[i + 1][2]
        else:
            end = pDomainList[i][2]
        midpos = row[1] + ((row[2] - row[1]) / 2)
        log.debug('tow {} hoe {}'.format(str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])  , str(chromosom) + ':' + str(start) + '-' + str(end) ))

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

            # log.debug('hoe {}'.format())
            # get intra-TAD data
            # hic_matrix_target = hm.hiCMatrix(
            #     pMatrixFile=pMatrixTarget, pChrnameList=[str(row[0])])
            # hic_matrix_control = hm.hiCMatrix(
            #     pMatrixFile=pMatrixControl, pChrnameList=[str(row[0])])
            # matrix_target = hic_matrix_target.matrix.toarray()
            # matrix_control = hic_matrix_control.matrix.toarray()

            # hic_matrix_target_inter_tad = hm.hiCMatrix(
            #     pMatrixFile=pMatrixTarget, pChrnameList=[str(chromosom)])
            # hic_matrix_control_inter_tad = hm.hiCMatrix(
            #     pMatrixFile=pMatrixControl, pChrnameList=[str(chromosom)])
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
        tad_midpoint = hic_matrix_target_inter_tad.getRegionBinRange(str(row[0]), midpos, midpos)[0]

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
            outer_left_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[0] #####
            outer_left_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

            outer_right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]
            outer_right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]


        if i + 1 < len(pDomainList) and not pCoolOrH5:
        # get index position left tad with tad
            right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0] #####
            right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]

        log.debug('right: {}'.format(hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])))
        log.debug('left: {}'.format(hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])))

        log.debug('outer_left_boundary_index_target {}'.format(outer_left_boundary_index_target))
        log.debug('left_boundary_index_target {}'.format(left_boundary_index_target))
        log.debug('right_boundary_index_target {}'.format(right_boundary_index_target))
        # log.debug('outer_left_boundary_index_target {}'.format(outer_left_boundary_index_target))

        if i - 1 >= 0 and i + 1 < len(pDomainList):
            intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
            intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
            intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()
            intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

        elif i - 1 < 0 and i + 1 < len(pDomainList):
            intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
            intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

        elif i - 1 > 0 and i + 1 >= len(pDomainList):
            intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
            intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()

        significance_level_left = None
        significance_level_right = None
        statistic_left = None
        statistic_right = None

        if i - 1 >= 0 and i + 1 < len(pDomainList):
            intertad_left_target = intertad_left_target.flatten()
            intertad_left_control = intertad_left_control.flatten()
            intertad_right_target = intertad_right_target.flatten()
            intertad_right_control = intertad_right_control.flatten()

            statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)
            statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
        elif i - 1 < 0 and i + 1 < len(pDomainList):
            intertad_right_target = intertad_right_target.flatten()
            intertad_right_control = intertad_right_control.flatten()
            statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
        elif i - 1 > 0 and i + 1 >= len(pDomainList):
            intertad_left_target = intertad_left_target.flatten()
            intertad_left_control = intertad_left_control.flatten()
            log.debug('intertad_left_target {}'.format(intertad_left_target))
            log.debug('intertad_left_control {}'.format(intertad_left_control))

            statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)

        # log.debug('matrix_target {}'.format(matrix_target))
        # log.debug('matrix_control {}'.format(matrix_control))

        statistic, significance_level = ranksums(matrix_target, matrix_control)
        log.debug('statistic {}, significance_level {}'.format(statistic, significance_level))
        log.debug('right statistic {}, significance_level {}'.format(statistic_right, significance_level_right))
        log.debug('left statistic {}, significance_level {}'.format(statistic_left, significance_level_left))


        p_values = []
        if significance_level_left is None or np.isnan(significance_level_left):
            accepted_inter_left.append(0)
            p_values.append(np.nan)
        elif significance_level_left <= pPValue:
            accepted_inter_left.append(1)
            p_values.append(significance_level_left)
        else:
            accepted_inter_left.append(0)
            p_values.append(significance_level_left)

        if significance_level_right is None or np.isnan(significance_level_right):
            accepted_inter_right.append(0)
            p_values.append(np.nan)
        elif significance_level_right <= pPValue:
            accepted_inter_right.append(1)
            p_values.append(significance_level_right)
        else:
            accepted_inter_right.append(0)
            p_values.append(significance_level_right)

        if significance_level is None or np.isnan(significance_level):
            accepted_intra.append(0)
            p_values.append(np.nan)
        elif significance_level <= pPValue:
            accepted_intra.append(1)
            p_values.append(significance_level)
        else:
            accepted_intra.append(0)
            p_values.append(significance_level)

        p_values_list.append(p_values)

        rows.append(row)
    # hic_matrix_target_inter_tad.save('manipulated_target.cool')
    # hic_matrix_control_inter_tad.save('manipulated_control.cool')
    pQueue.put([p_values_list, accepted_inter_left, accepted_inter_right, accepted_intra, rows])


def main(args=None):
    args = parse_arguments().parse_args(args)

    # read domains file
    domains_df = readDomainBoundaries(args.tadDomains)
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
    accepted_H0 = []
    rejected_H0 = []
    # log.debug('domains_df {}'.format(domains_df))
    domains = domains_df.values.tolist()

    p_values_threads = [None] * args.threads
    accepted_left_inter_threads = [None] * args.threads
    accepted_right_inter_threads = [None] * args.threads
    accepted_intra_threads = [None] * args.threads
    rows_threads = [None] * args.threads

    domainsPerThread = len(domains) // args.threads
    all_data_collected = False
    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads

    # None --> first thread, process first element in list, ignore last one
    # True --> middle thread: ignore first and last element in tad processing
    # False --> last thread: ignore first element, process last one
    thread_id = None
    for i in range(args.threads):

        if i == 0:
            domainListThread = domains[i * domainsPerThread:((i + 1) * domainsPerThread)+1]
            thread_id = None
        elif i < args.threads - 1:
            domainListThread = domains[(i * domainsPerThread)-1:((i + 1) * domainsPerThread)+1]
            thread_id = True

        else:
            domainListThread = domains[(i * domainsPerThread)-1:]
            thread_id = False

        if args.threads == 1:
            thread_id = ''
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

    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                p_values_threads[i], accepted_left_inter_threads[i], \
                    accepted_right_inter_threads[i], \
                    accepted_intra_threads[i], rows_threads[i] = queue[i].get()

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

    # outfile_names = [item for sublist in outfile_names for item in sublist]
    # target_list_name = [
    #     item for sublist in target_list_name for item in sublist]

    p_values_list = [item for sublist in p_values_threads for item in sublist]
    accepted_inter_left = [item for sublist in accepted_left_inter_threads for item in sublist]
    accepted_inter_right = [item for sublist in accepted_right_inter_threads for item in sublist]
    accepted_intra = [item for sublist in accepted_intra_threads for item in sublist]
    rows = [item for sublist in rows_threads for item in sublist]

    p_values_list = np.array(p_values_list)
    accepted_inter_left = np.array(accepted_inter_left)
    accepted_inter_right = np.array(accepted_inter_right)
    accepted_intra = np.array(accepted_intra)
    rows = np.array(rows)

    if args.mode == 'intra-TAD':
        mask = accepted_intra
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



    accepted_HO = p_values_list[~mask]
    rejected_H0 = p_values_list[mask]
    accepted_rows = rows[~mask]
    rejected_rows = rows[mask]
    with open(args.outFileNamePrefix + '_accepted.diff_tad', 'w') as file:
        header = '# Created with HiCExplorer\'s hicDifferentialTAD version ' + __version__ + '\n'
        header += '# H0 \'regions are equal\' H0 is accepted for all p-value greater the user given p-value threshold; i.e. regions in this file are not considered as differential.\n'
        header += '# Accepted regions with Wilcoxon rank-sum test to p-value: {}  with used mode: {} and modeReject: {} \n'.format(args.pValue, args.mode, args.modeReject )
        header += '# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD\tp-value right-inter-TAD\tp-value intra-TAD\n'
        file.write(header)
        for i, row in enumerate(accepted_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\t')
            pvalue_list = list(map(str, accepted_HO[i]))
            file.write('\t'.join(pvalue_list))

            file.write('\n')
    with open(args.outFileNamePrefix + '_rejected.diff_tad', 'w') as file:
        header = '# Created with HiCExplorer\'s hicDifferentialTAD version ' + __version__ + '\n'
        header += '# H0 \'regions are equal\' H0 is rejected for all p-value smaller or equal the user given p-value threshold; i.e. regions in this file are considered as differential.\n'
        header += '# Rejected regions with Wilcoxon rank-sum test to p-value: {}  with used mode: {} and modeReject: {} \n'.format(args.pValue, args.mode, args.modeReject )
        header += '# Chromosome\tstart\tend\tname\tscore\tstrand\tp-value left-inter-TAD\tp-value right-inter-TAD\tp-value intra-TAD\n'
       
        file.write(header)

        for i, row in enumerate(rejected_rows):
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\t')
            pvalue_list = list(map(str, rejected_H0[i]))
            file.write('\t'.join(pvalue_list))
            file.write('\n')
