import argparse
import sys
import os
import math
import logging
log = logging.getLogger(__name__)

import numpy as np
from scipy import stats

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Test per line if two samples are differential expressed via chi2 contingency test.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for differential test.',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--alpha', '-a',
                                help='Accept all samples to significance level alpha',
                                type=float,
                                default=0.05,
                                required=True)
    
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--outputFolder', '-o',
                           help='Output folder of the files.',
                           required=False,
                           default='differentialResults')
    parserOpt.add_argument('--statisticTest',
                           help='Type of test used for testing: fisher\'s exact test or chi2 contingency',
                           choices=['fisher', 'chi2'],
                           default='fisher')
    parserOpt.add_argument('--batchMode', '-bm',
                           help='The given file for --interactionFile and or --targetFile contain a list of the to be processed files.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def readInteractionFile(pInteractionFile):

    line_content = []
    data = []

    with open(pInteractionFile, 'r') as file:
        header = file.readline()
        sum_of_all_interactions = float(header.strip().split('\t')[-1].split(' ')[-1])
        header += file.readline()
        for line in file.readlines():
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
            chi2, p_value, dof, ex = stats.chi2_contingency([group1, group2], correction=False)
            if chi2 >= critical_value:
                test_result.append(p_value)
                rejected.append([i, p_value])
            else:
                test_result.append(p_value)
                accepted.append([i, p_value])

        except ValueError:
            zero_values_counter += 1
            test_result.append(np.nan)

    if zero_values_counter > 0:
        log.info('{} samples were not tested because at least one condition contained no data in both groups.'.format(zero_values_counter))
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
    return test_result, accepted, rejected


def writeResult(pOutFileName, pData, pHeaderOld, pHeaderNew, pViewpoint1, pViewpoint2, pAlpha, pTest):

    with open(pOutFileName, 'w') as file:
        header = '# Differential analysis result file of HiCExplorer\'s chicDifferentialTest version '
        header += str(__version__)
        header += '\n'

        header += '# This file contains the p-values computed by {} test\n'.format(pTest)
        header += '# To test the smoothed (float) values they were rounded up to the next integer\n'
        header += '#\n'

        header += ' '.join(['# Alpha level', str(pAlpha)])
        header += '\n'
        header += ' '.join(['# Degrees of freedom', '1'])
        header += '\n#\n'

        file.write(header)

        # file.write('#Viewpoint\t\t\ttarget\t\t\tgene\tcondition1\t\tcondition2\t\tp-value\n')
        file.write(pHeaderOld.split('\n')[0] + '\n')
        file.write(pHeaderNew.split('\n')[0] + '\n')

        # pHeaderNew.split('\n')[0]

        # log.debug('pHeaderOld {}'.format(pHeaderOld))
        # log.debug('pHeaderNew {}'.format(pHeaderNew))

        file.write('#Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of interactions 2\ttarget_2 raw\tp-value\n')

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


def main(args=None):
    args = parse_arguments().parse_args(args)
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    interactionFileList = []
    if args.batchMode:
        with open(args.interactionFile[0], 'r') as interactionFile:
            file_ = True
            while file_:
            # for line in fh.readlines():
                file_ = interactionFile.readline().strip()
                file2_ = interactionFile.readline().strip()
                if file_ != '' and file2_ != '':
                    interactionFileList.append((file_, file2_))
            log.debug('interactionFileList {}'.format(interactionFileList))
    else:
        if len(args.interactionFile) % 2 == 0:
            
            i = 0
            while i < len(args.interactionFile):
                interactionFileList.append((args.interactionFile[i], args.interactionFile[i + 1]))
                i += 2

    for interactionFile in interactionFileList:

        header1, line_content1, data1 = readInteractionFile(args.interactionFileFolder + '/' + interactionFile[0])
        header2, line_content2, data2 = readInteractionFile(args.interactionFileFolder + '/' + interactionFile[1])

        if args.statisticTest == 'chi2':
            test_result, accepted, rejected = chisquare_test(data1, data2, args.alpha)
        elif args.statisticTest == 'fisher':
            test_result, accepted, rejected = fisher_exact_test(data1, data2, args.alpha)

        write_out_lines = []
        for i, result in enumerate(test_result):
            write_out_lines.append([line_content1[i], line_content2[i], result, data1[i], data2[i]])

        write_out_lines_accepted = []
        for result in accepted:
            write_out_lines_accepted.append([line_content1[result[0]], line_content2[result[0]], result[1], data1[result[0]], data2[result[0]]] )

        write_out_lines_rejected = []
        for result in rejected:
            log.debug('result[1] {}'.format(result[1]))
            write_out_lines_rejected.append([line_content1[result[0]], line_content2[result[0]], result[1], data1[result[0]], data2[result[0]]] )


        header_new = interactionFile[0]
        header_new += ' '
        header_new += interactionFile[1]

        sample_prefix = interactionFile[0].split('/')[-1].split('_')[0] + '_' + interactionFile[1].split('/')[-1].split('_')[0]
        region_prefix = '_'.join(interactionFile[0].split('/')[-1].split('_')[1:6])
        outFileName = sample_prefix + '_' + region_prefix 
        outFileName_accepted = args.outputFolder + '/' + outFileName + '_H0_accepted.bed'
        outFileName_rejected = args.outputFolder + '/' + outFileName + '_H0_rejected.bed'
        outFileName = args.outputFolder + '/' + outFileName + '_results.bed'

        # outFileName = args.outFileName.split('.')

        # resultsNameFile = outFileName[0] + '_results.bed'

        writeResult(outFileName, write_out_lines, header1, header2, line_content1[0][:4], line_content2[0][:4], args.alpha, args.statisticTest)
        writeResult(outFileName_accepted, write_out_lines_accepted, header1, header2, line_content1[0][:4], line_content2[0][:4], args.alpha, args.statisticTest)
        writeResult(outFileName_rejected, write_out_lines_rejected, header1, header2, line_content1[0][:4], line_content2[0][:4], args.alpha, args.statisticTest)
