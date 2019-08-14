import argparse
import sys
import os
import errno
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import pybedtools
import numpy as np

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=
"""
Per viewpoint the significant interactions are detected based on the background model. Each viewpoint file gets as output a file with all recorded significant interactions and
a target file. The target file is especially in the batch mode context useful, it merges for two consecutive listed control and treatment viewpoint the significant interactions which can then be used
to test for a differential interaction scheme.
 
chicSignificantInteractions supports two modes to detect significant interactions, either by an x-fold over the average background or a loose p-value. In both cases neighboring significant peaks are merged together and an additional 
p-value based on the sum of interactions for this neighborhood is computed. Only interactions with a higher p-value as specified by the threshold `--pValue` are accepted as a significant interaction.

An example usage is for single mode is:

`$ chicSignificantInteractions --interactionFile interactionFilesFolder/Sox17_FL-E13-5_chr1_1000_2000.bed --referencePoints referencePointsFile.bed --range 20000 40000 --backgroundModelFile background_model.bed --loosePValue 0.5 --pValue 0.01`

An example usage is for batch mode is:

`$ chicViewpointBackgroundModel --matrices matrix1.cool matrix2.cool matrix3.cool --referencePoints referencePointsFile.bed --range 20000 40000 --outFileName background_model.bed`

The main difference between single mode and batch mode is that in single mode the parameter `--interactionFile` is interpreted as a list of viewpoint files created with `chicViewpoint`, whereas in batch mode only one file is allowed which contains per line the file names of viewpoint files. 
This file is created by `chicViewpoint` and the parameter `--writeFileNamesToFile`. Please have in mind to specify in batch mode the folder via `--interactionFileFolder` where `chicViewpoint` wrote the files to.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for aggregation of the statistics.',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--pValue', '-p',
                                help='p-value threshold value to filter target regions to include them for differential analysis.',
                                type=float,
                                required=True)
    parserMutuallyExclusiveGroupFilter = parser.add_mutually_exclusive_group(
        required=True)
    parserMutuallyExclusiveGroupFilter.add_argument('--xFoldBackground', '-xf',
                                                    help='Filter x-fold over background. Used to merge neighboring bins with a broader peak but '
                                                    'less significant interactions to one peak with high significance. Used only for pValue option.',
                                                    type=float

                                                    )
    parserMutuallyExclusiveGroupFilter.add_argument('--loosePValue', '-lp',
                                                    help='loose p-value threshold value to filter target regions in a first round. '
                                                    'Used to merge neighboring bins with a broader peak but less significant interactions to one peak with high significance.'
                                                    ' Used only for pValue option.',
                                                    type=float)
    parserRequired.add_argument('--backgroundModelFile', '-bmf',
                                help='path to the background file which is necessary to compute the rbz-score',
                                required=True)
    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream',
                                required=True,
                                type=int,
                                nargs=2)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileNameSuffix', '-suffix',
                           help='File name suffix to save the results, prefix is the input file name.',
                           required=False,
                           default='_significant_interactions.bed')

    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--targetFolder', '-tf',
                           help='Folder where the target files are stored.',
                           required=False,
                           default='targetFolder')
    parserOpt.add_argument('--outputFolder', '-o',
                           help='Output folder of the files significant interaction files.',
                           required=False,
                           default='significantFiles')
    parserOpt.add_argument('--writeFileNamesToFile', '-w',
                           help='',
                           default='significantFilesBatch.txt')
    parserOpt.add_argument('--targetFileList', '-tl',
                           help='The file to store the target file names.',
                           default='targetList.txt')
    parserOpt.add_argument('--batchMode', '-bm',
                           help='The given file for --interactionFile and or --targetFile contain a list of the to be processed files.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate range of backgroundmodel starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin.',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument('--xFoldMaxValueNB', '-xfnb',
                           help='x-fold factor to increase the number of precomputed p-values per relative genomic distance. If set to 1, the maximal distance is used. ',
                           required=False,
                           default=10,
                           type=int
                           )
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=int,
                           default=5,
                           help='The minimum number of interactions a detected peaks needs to have to be considered.')

    parserOpt.add_argument('--resolution', '-r',
                           help='Resolution of the bin in genomic units. Values are usually e.g. 1000 for a 1kb, 5000 for a 5kb or 10000 for a 10kb resolution.'
                           'This value is used to merge neighboring bins.',
                           type=int,
                           default=1000,
                           required=False)
                           
    parserOpt.add_argument('--computeSampleNumber', '-csn',
                           help='Number of samples to compute together. Applies only in batch mode.',
                           required=False,
                           default=2,
                           type=int)

    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackgroundSumOfDensities, pQueue=None):
    outfile_names = []
    target_outfile_names = []

    for interactionFile in pInteractionFilesList:
        target_list = []
        sample_prefix = ''
        for sample in interactionFile:
            # header,
            # interaction_data:rel interaction, p-value, raw, x-fold::{-1000:[0.1, 0.01, 2.3, 5]},
            if pArgs.interactionFileFolder != '.':
                absolute_sample_path = pArgs.interactionFileFolder + '/' + sample
            else:
                absolute_sample_path = sample
            data = pViewpointObj.readInteractionFileForAggregateStatistics(absolute_sample_path)
            sample_prefix += sample.split('/')[-1].split('_')[0]
            sample_prefix += '_'
            # filter by x-fold over background value or loose p-value
            # and merge neighbors. Use center position to compute new p-value.
            if pArgs.xFoldBackground is not None:
                accepted_scores, merged_lines_dict = merge_neighbors_x_fold(
                    pArgs.xFoldBackground, data, pViewpointObj, pResolution=pArgs.resolution)
            else:
                accepted_scores, merged_lines_dict = merge_neighbors_loose_p_value(
                    pArgs.loosePValue, data, pViewpointObj, pResolution=pArgs.resolution)

            # compute new p-values
            accepted_scores, target_lines = compute_new_p_values(
                accepted_scores, pBackgroundSumOfDensities, pArgs.pValue, merged_lines_dict, pArgs.peakInteractionsThreshold)

            # filter by new p-value
            if len(accepted_scores) == 0:
                if pArgs.batchMode:
                    with open('errorLog.txt', 'a+') as errorlog:
                        errorlog.write('Failed for: {} and {}.\n'.format(
                            interactionFile[0], interactionFile[1]))
                else:
                    log.info('No target regions found')
            outFileName = '.'.join(sample.split(
                '/')[-1].split('.')[:-1]) + '_' + pArgs.outFileNameSuffix
            if pArgs.batchMode:
                outfile_names.append(outFileName)
            outFileName = pArgs.outputFolder + '/' + outFileName
            # write only significant lines to file
            write(outFileName, data[0], accepted_scores)
            target_list.append(target_lines)

        target_list = [item for sublist in target_list for item in sublist]
        log.debug('interactionFile {}'.format(interactionFile))
        sample_name = '_'.join(interactionFile[0].split('/')[-1].split('.')[0].split('_')[1:])
        target_name = sample_prefix + sample_name + '_target.bed'
        target_outfile_names.append(target_name)
        target_name = pArgs.targetFolder + '/' + target_name
        writeTargetList(target_list, target_name, pArgs)
    if pQueue is None:
        return target_outfile_names
    pQueue.put([outfile_names, target_outfile_names])
    return


def compute_new_p_values(pData, pBackgroundSumOfDensities, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
    accepted = {}
    accepted_lines = []
    for key in pData:
        if key in pBackgroundSumOfDensities:
            if int(float(pData[key][-1])) - 1 < 0:
                pData[key][-3] = pBackgroundSumOfDensities[key][0]
            else:
                try:
                    if int(float(pData[key][-1])) < len(pBackgroundSumOfDensities[key]):
                        pData[key][-3] = 1 - \
                            pBackgroundSumOfDensities[key][int(
                                float(pData[key][-1]))]
                    else:
                        pData[key][-3] = 1 - pBackgroundSumOfDensities[key][-1]

                except Exception:
                    pData[key][-3] = 1 - pBackgroundSumOfDensities[key][-1]
                    log.error('Not enough p-values precomputed, using highest value instead. Please increase --xFoldMaxValueNB value. Value {}, max value {}'.format(
                        int(float(pData[key][-1])), len(pData[key])))

            if pData[key][-3] <= pPValue:
                if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                    accepted[key] = pData[key]
                    target_content = pMergedLinesDict[key][0][:3]
                    target_content[2] = pMergedLinesDict[key][-1][2]
                    accepted_lines.append(target_content)

    return accepted, accepted_lines


def merge_neighbors_x_fold(pXfold, pData, pViewpointObj, pResolution):
    accepted = {}
    accepted_line = {}
    for key in pData[1]:

        if pData[1][key][-1] < pXfold:
            continue
        accepted[key] = pData[1][key]
        accepted_line[key] = pData[2][key]

    if accepted_line:
        return pViewpointObj.merge_neighbors(accepted_line, pMergeThreshold=pResolution)
    return accepted_line, None


def merge_neighbors_loose_p_value(pLoosePValue, pData, pViewpointObj, pResolution):
    accepted = {}
    accepted_line = {}
    for key in pData[1]:

        if pData[1][key][1] > pLoosePValue:
            continue
        accepted[key] = pData[1][key]
        accepted_line[key] = pData[2][key]
    if accepted_line:
        return pViewpointObj.merge_neighbors(accepted_line, pMergeThreshold=pResolution)
    return accepted_line, None


def write(pOutFileName, pHeader, pInteractionLines):

    # sum_of_interactions = float(pHeader.split('\t')[-1].split(' ')[-1])
    with open(pOutFileName, 'w') as file:
        file.write(pHeader)
        file.write(
            '#Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative interactions\tp-value\tx-fold\tRaw target')
        file.write('\n')

        for data in pInteractionLines:
            new_line = '\t'.join(pInteractionLines[data][:6])
            new_line += '\t' + '\t'.join(format(float(x), "0.20f")
                                         for x in pInteractionLines[data][6:])
            new_line += '\n'
            file.write(new_line)


def call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackgroundSumOfDensities):
    outfile_names = [None] * pArgs.threads
    target_list_name = [None] * pArgs.threads
    interactionFilesPerThread = len(pInteractionFilesList) // pArgs.threads
    all_data_collected = False
    queue = [None] * pArgs.threads
    process = [None] * pArgs.threads
    thread_done = [False] * pArgs.threads
    for i in range(pArgs.threads):

        if i < pArgs.threads - 1:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:(
                i + 1) * interactionFilesPerThread]
        else:
            interactionFileListThread = pInteractionFilesList[i *
                                                              interactionFilesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=compute_interaction_file, kwargs=dict(
            pInteractionFilesList=interactionFileListThread,
            pArgs=pArgs,
            pViewpointObj=pViewpointObj,
            pBackgroundSumOfDensities=pBackgroundSumOfDensities,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(pArgs.threads):
            if queue[i] is not None and not queue[i].empty():
                background_data_thread = queue[i].get()
                outfile_names[i], target_list_name[i] = background_data_thread
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

    outfile_names = [item for sublist in outfile_names for item in sublist]
    target_list_name = [
        item for sublist in target_list_name for item in sublist]

    return outfile_names, target_list_name


def writeTargetList(pTargetList, pOutFileName, pArgs):
    # remove duplicates
    target_list_ = []
    for line in pTargetList:
        target_list_.append('\t'.join(line))
    target_set = set(target_list_)
    pTargetList = sorted(list(target_set))

    a = pybedtools.BedTool(pTargetList)

    header = '# Significant interactions result file of HiCExplorer\'s chicSignificantInteractions version '
    header += str(__version__)
    header += '\n'
    header += '# {}\n'.format(pOutFileName)

    if pArgs.xFoldBackground:
        header += '# Mode: x-fold background with '
        header += str(pArgs.xFoldBackground)
        header += '\n'
    else:
        header += '# Mode: loose p-value with '
        header += str(pArgs.loosePValue)
        header += '\n'
    header += '# Used p-value: '
    header += str(pArgs.pValue)
    header += '\n#\n'
    if len(pTargetList) == 0:
        with open(pOutFileName, 'w') as file:
            file.write(header)
    else:
        a.sort().merge(d=pArgs.resolution).saveas(pOutFileName, trackline=header)


def main(args=None):
    args = parse_arguments().parse_args(args)
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if not os.path.exists(args.targetFolder):
        try:
            os.makedirs(args.targetFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    viewpointObj = Viewpoint()
    outfile_names = []

    interactionFileList = []

    background_model = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range)
    background_sum_of_densities_dict = viewpointObj.computeSumOfDensities(
        background_model, args, pXfoldMaxValue=args.xFoldMaxValueNB)

    if args.batchMode:
        with open(args.interactionFile[0], 'r') as interactionFile:
            file_ = True
            while file_:
                lines = []
                for i in range(0, args.computeSampleNumber):
                    file_ = interactionFile.readline().strip()
                    if file_ != '':
                        lines.append(file_)
                if len(lines) > 0:
                    interactionFileList.append(lines)
        log.debug('interactionFileList {}'.format(interactionFileList))
        outfile_names, target_list_name = call_multi_core(
            interactionFileList, args, viewpointObj, background_sum_of_densities_dict)

    else:
        i = 0
        while i < len(args.interactionFile):
            lines = []
            for j in range(0, args.computeSampleNumber):
                if i < len(args.interactionFile):
                    lines.append(args.interactionFile[i])
                i += 1
            interactionFileList.append(lines)
            
        target_list_name = compute_interaction_file(
            interactionFileList, args, viewpointObj, background_sum_of_densities_dict)

    if args.batchMode:
        with open(args.writeFileNamesToFile, 'w') as nameListFile:
            nameListFile.write('\n'.join(outfile_names))

        with open(args.targetFileList, 'w') as targetNamesFile:
            targetNamesFile.write('\n'.join(target_list_name))
