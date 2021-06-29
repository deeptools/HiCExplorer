import argparse
import sys
import errno
import os
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy import stats
import h5py

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
chicDifferentialTest tests if two locations under consideration of the reference point have a different interaction count. For this either Fisher's test or the chi2 contingency test can be used.
The files that are accepted for this test can be created with `chicAggregateStatistic`. H0 assumes the interactions are not different. Therefore the differential interaction counts are all where H0 was rejected.

"""
                                     )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--aggregatedFile', '-af',
                                help='path to the aggregated files which should be used for the differential test.',
                                required=True)

    parserRequired.add_argument('--alpha', '-a',
                                help='define a significance level (alpha) for accepting samples',
                                type=float,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='Output file for the differential results'
                           ' (Default: %(default)s).',
                           required=False,
                           default='differentialResults.hdf5')
    parserOpt.add_argument('--statisticTest',
                           help='Type of test used: fisher\'s exact test or chi2 contingency'
                           ' (Default: %(default)s).',
                           choices=['fisher', 'chi2'],
                           default='fisher')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def readInteractionFile(pInteractionFile):

    line_content = []
    data = []

    with open(pInteractionFile, 'r') as file:
        file.readline()
        header = file.readline()
        sum_of_all_interactions = float(
            header.strip().split('\t')[-1].split(' ')[-1])
        header += file.readline()
        for line in file.readlines():
            if line.startswith('#'):
                continue
            _line = line.strip().split('\t')
            if len(_line) <= 1:
                continue
            line_content.append(_line)
            data.append([sum_of_all_interactions, float(_line[-1])])

    return header, line_content, data


def chisquare_test(pDataFile1, pDataFile2, pAlpha):
    # pair of accepted/unaccepted and pvalue
    # True is rejection of H0
    # False acceptance of H0
    test_result = []
    accepted = []
    rejected = []
    # Find the critical value for alpha confidence level
    critical_value = stats.chi2.ppf(q=1 - pAlpha, df=1)
    zero_values_counter = 0
    for i, (group1, group2) in enumerate(zip(pDataFile1, pDataFile2)):
        try:
            chi2, p_value, dof, ex = stats.chi2_contingency(
                [group1, group2], correction=False)
            if chi2 >= critical_value:
                test_result.append(p_value)
                rejected.append([i, p_value])
            else:
                test_result.append(p_value)
                accepted.append([i, p_value])

        except ValueError:
            zero_values_counter += 1
            test_result.append(np.nan)
            accepted.append([i, 1.0])

    if zero_values_counter > 0:
        log.info('{} samples were not tested because at least one condition contained no data in both groups.'.format(
            zero_values_counter))
    return test_result, accepted, rejected


def fisher_exact_test(pDataFile1, pDataFile2, pAlpha):

    test_result = []
    accepted = []
    rejected = []
    for i, (group1, group2) in enumerate(zip(pDataFile1, pDataFile2)):
        try:
            odds, p_value = stats.fisher_exact(np.ceil([group1, group2]))
            if p_value <= pAlpha:
                test_result.append(p_value)
                rejected.append([i, p_value])
            else:
                test_result.append(p_value)
                accepted.append([i, p_value])
        except ValueError:
            test_result.append(np.nan)
            accepted.append([i, 1.0])
    return test_result, accepted, rejected


def writeResult(pOutFileName, pData, pHeaderOld, pHeaderNew, pAlpha, pTest):

    with open(pOutFileName, 'w') as file:
        header = '# Differential analysis result file of HiCExplorer\'s chicDifferentialTest version '
        header += str(__version__)
        header += '\n'

        header += '# This file contains the p-values computed by {} test\n'.format(
            pTest)
        header += '# To test the smoothed (float) values they were rounded up to the next integer\n'
        header += '#\n'

        header += ' '.join(['# Alpha level', str(pAlpha)])
        header += '\n'
        header += ' '.join(['# Degrees of freedom', '1'])
        header += '\n#\n'

        file.write(header)

        file.write(pHeaderOld.split('\n')[0] + '\n')
        file.write(pHeaderNew.split('\n')[0] + '\n')

        file.write('#Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of interactions 2\ttarget_2 raw\tp-value\n')

        if pData:
            for data in pData:
                line = '\t'.join(data[0][:4])
                line += '\t'

                line += '{}'.format(data[0][5])
                line += '\t'
                line += '\t'.join(format(x, '.5f') for x in data[3])
                line += '\t'

                line += '\t'.join(format(x, '.5f') for x in data[4])
                line += '\t'

                line += '\t{}\n'.format(format(data[2], '.5f'))
                file.write(line)


def writeResultHDF(pOutFileName, pAcceptedData, pRejectedData, pAllResultData, pInputData, pAlpha, pTest):
    resultFileH5Object = h5py.File(pOutFileName, 'w')
    resultFileH5Object.attrs['type'] = "differential"
    resultFileH5Object.attrs['alpha'] = pAlpha
    resultFileH5Object.attrs['test'] = pTest

    all_data_dict = {'accepted': pAcceptedData, 'rejected': pRejectedData, 'all': pAllResultData}
    for i, inputData in enumerate(pInputData):
        matrix1_name = inputData[0][1]
        matrix2_name = inputData[1][1]
        chromosome = inputData[0][2]
        gene_name = inputData[0][3]

        if matrix1_name not in resultFileH5Object:
            matrix1_object = resultFileH5Object.create_group(matrix1_name)
        else:
            matrix1_object = resultFileH5Object[matrix1_name]

        if matrix2_name not in matrix1_object:
            matrix2_object = matrix1_object.create_group(matrix2_name)
        else:
            matrix2_object = matrix1_object[matrix2_name]

        if 'genes' not in matrix2_object:
            geneGroup = matrix2_object.create_group('genes')
        else:
            geneGroup = matrix2_object['genes']

        if chromosome not in matrix2_object:
            chromosome_object = matrix2_object.create_group(chromosome)
        else:
            chromosome_object = matrix2_object[chromosome]

        gene_object = chromosome_object.create_group(gene_name)

        gene_object.create_group('accepted')
        gene_object.create_group('rejected')
        gene_object.create_group('all')

        for category in ['accepted', 'rejected', 'all']:
            write_object = gene_object[category]
            data_object = all_data_dict[category][i]
            if len(data_object) == 0:
                continue
            chromosome = None
            start_list = []
            end_list = []
            sum_of_interactions_1 = None
            sum_of_interactions_2 = None

            relative_distance_list = []
            pvalue_list = []

            raw_target_list_1 = []
            raw_target_list_2 = []

            for data in data_object:

                chromosome = data[0][0]
                start_list.append(data[0][1])
                end_list.append(data[0][2])

                relative_distance_list.append(data[0][5])

                sum_of_interactions_1 = float(data[3][0])
                sum_of_interactions_2 = float(data[4][0])

                raw_target_list_1.append(data[3][1])
                raw_target_list_2.append(data[4][1])
                pvalue_list.append(data[2])

            write_object["chromosome"] = str(chromosome)
            write_object.create_dataset("start_list", data=start_list, compression="gzip", compression_opts=9)
            write_object.create_dataset("end_list", data=end_list, compression="gzip", compression_opts=9)
            write_object["gene"] = str(gene_name)
            write_object.create_dataset("relative_distance_list", data=relative_distance_list, compression="gzip", compression_opts=9)

            write_object["sum_of_interactions_1"] = float(sum_of_interactions_1)
            write_object["sum_of_interactions_2"] = float(sum_of_interactions_2)

            write_object.create_dataset("raw_target_list_1", data=raw_target_list_1, compression="gzip", compression_opts=9)
            write_object.create_dataset("raw_target_list_2", data=raw_target_list_2, compression="gzip", compression_opts=9)
            write_object.create_dataset("pvalue_list", data=pvalue_list, compression="gzip", compression_opts=9)

        try:
            geneGroup[gene_name] = chromosome_object[gene_name]
        except Exception as exp:
            log.debug('exception {}'.format(str(exp)))
    resultFileH5Object.close()


def run_statistical_tests(pInteractionFilesList, pArgs, pViewpointObject, pQueue=None):
    accepted_list = []
    rejected_list = []
    all_results_list = []
    try:
        for interactionFile in pInteractionFilesList:

            line_content1, data1 = pViewpointObject.readAggregatedFileHDF(pArgs.aggregatedFile, interactionFile[0])
            line_content2, data2 = pViewpointObject.readAggregatedFileHDF(pArgs.aggregatedFile, interactionFile[1])

            if len(line_content1) == 0 or len(line_content2) == 0:
                continue
            if pArgs.statisticTest == 'chi2':
                test_result, accepted, rejected = chisquare_test(
                    data1, data2, pArgs.alpha)
            elif pArgs.statisticTest == 'fisher':
                test_result, accepted, rejected = fisher_exact_test(
                    data1, data2, pArgs.alpha)

            write_out_lines = []
            for i, result in enumerate(test_result):
                write_out_lines.append(
                    [line_content1[i], line_content2[i], result, data1[i], data2[i]])

            write_out_lines_accepted = []
            for result in accepted:
                write_out_lines_accepted.append(
                    [line_content1[result[0]], line_content2[result[0]], result[1], data1[result[0]], data2[result[0]]])

            write_out_lines_rejected = []
            for result in rejected:
                write_out_lines_rejected.append(
                    [line_content1[result[0]], line_content2[result[0]], result[1], data1[result[0]], data2[result[0]]])

            accepted_list.append(write_out_lines_accepted)
            rejected_list.append(write_out_lines_rejected)
            all_results_list.append(write_out_lines)

    except Exception as exp:
        pQueue.put('Fail: ' + str(exp))
        return

    if pQueue is None:
        return
    pQueue.put([accepted_list, rejected_list, all_results_list])
    return


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()

    aggregatedList = []

    aggregatedFileHDF5Object = h5py.File(args.aggregatedFile, 'r')
    if aggregatedFileHDF5Object.attrs['type'] != 'aggregate':
        log.error('Please provide a file created by chicAggregateStatistic for the parameter --aggregatedFile.')
        exit(1)

    keys_aggregatedFile = list(aggregatedFileHDF5Object.keys())

    for i, combinationOfMatrix in enumerate(keys_aggregatedFile):
        keys_matrix_intern = list(aggregatedFileHDF5Object[combinationOfMatrix].keys())
        if len(keys_matrix_intern) == 0:
            continue

        matrix1 = keys_matrix_intern[0]
        matrix2 = keys_matrix_intern[1]

        matrix_obj1 = aggregatedFileHDF5Object[combinationOfMatrix + '/' + matrix1]
        matrix_obj2 = aggregatedFileHDF5Object[combinationOfMatrix + '/' + matrix2]

        chromosomeList1 = sorted(list(matrix_obj1.keys()))
        chromosomeList2 = sorted(list(matrix_obj2.keys()))
        chromosomeList1.remove('genes')
        chromosomeList2.remove('genes')
        for chromosome1, chromosome2 in zip(chromosomeList1, chromosomeList2):
            geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
            geneList2 = sorted(list(matrix_obj2[chromosome2].keys()))

            for gene1, gene2 in zip(geneList1, geneList2):
                aggregatedList.append([[combinationOfMatrix, matrix1, chromosome1, gene1], [combinationOfMatrix, matrix2, chromosome2, gene2]])

    aggregatedFileHDF5Object.close()

    fail_flag = False
    fail_message = ''
    all_data = [None] * args.threads
    accepted_data = [None] * args.threads
    rejected_data = [None] * args.threads

    aggregatedListPerThread = len(aggregatedList) // args.threads
    all_data_collected = False
    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    length_of_threads = 0
    for i in range(args.threads):

        if i < args.threads - 1:
            aggregatedListThread = aggregatedList[i * aggregatedListPerThread:(i + 1) * aggregatedListPerThread]
        else:
            aggregatedListThread = aggregatedList[i * aggregatedListPerThread:]
        length_of_threads += len(aggregatedListThread)
        queue[i] = Queue()
        process[i] = Process(target=run_statistical_tests, kwargs=dict(
            pInteractionFilesList=aggregatedListThread,
            pArgs=args,
            pViewpointObject=viewpointObj,
            pQueue=queue[i]
        )
        )

        process[i].start()
    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                background_data_thread = queue[i].get()
                if 'Fail:' in background_data_thread:
                    fail_flag = True
                    fail_message = background_data_thread[6:]
                else:
                    accepted_data[i], rejected_data[i], all_data[i] = background_data_thread
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
        log.error(fail_message)
        exit(1)

    accepted_data = [item for sublist in accepted_data for item in sublist]
    rejected_data = [item for sublist in rejected_data for item in sublist]
    all_data = [item for sublist in all_data for item in sublist]

    writeResultHDF(args.outFileName, accepted_data, rejected_data, all_data, aggregatedList, args.alpha, args.statisticTest)
