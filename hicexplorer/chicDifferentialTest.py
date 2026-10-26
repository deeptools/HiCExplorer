"""chicDifferentialTest: a thin layer over the C++ executable chicDifferentialTest.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
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
    """The argparse parser of chicDifferentialTest, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicDifferentialTest')


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
    resultFileH5Object.attrs['version'] = __version__

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
    """Run the C++ implementation of chicDifferentialTest with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicDifferentialTest')(args)
