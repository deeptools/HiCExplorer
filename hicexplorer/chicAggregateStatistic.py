import argparse
import sys
import os
import math
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Aggregates the statistics of interaction files and prepares them for chicDifferentialTest')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for aggregation of the statistics.',
                                required=True,
                                nargs='+')

    parserMutuallyExclusiveGroup = parser.add_mutually_exclusive_group(
        required=True)
    parserMutuallyExclusiveGroup.add_argument('--targetFile', '-tf',
                                              help='path to the target files which contains the target regions to prepare data for differential analysis.'
                                              )
    parserMutuallyExclusiveGroup.add_argument('--pValue', '-p',
                                              help='p-value threshold value to filter target regions to include them for differential analysis.',
                                              type=float)
    # parserMutuallyExclusiveGroupFilter = parser.add_mutually_exclusive_group(
    #     required=True)
    # parserMutuallyExclusiveGroupFilter.add_argument('--xFoldBackground', '-xf',
    #                                           help='Filter x-fold over background. Used to merge neighboring bins with a broader peak but '
    #                                                 'less significant interactions to one peak with high significance. Used only for pValue option.',
    #                                           type=float
    #                                           )
    # parserMutuallyExclusiveGroupFilter.add_argument('--loosePValue', '-lp',
    #                                           help='loose p-value threshold value to filter target regions in a first round. '
    #                                                 'Used to merge neighboring bins with a broader peak but less significant interactions to one peak with high significance.'
    #                                                 ' Used only for pValue option.',
    #                                           type=float)

    parserOpt = parser.add_argument_group('Optional arguments')

    # parserOpt.add_argument('--xFoldBackground', '-xf',
    #                         help='Filter x-fold over background. Used only for pValue option.',
    #                         type=float)
    parserOpt.add_argument('--outFileNameSuffix', '-suffix',
                           help='File name suffix to save the result.',
                           required=False,
                           default='_aggregate_target.bed')

    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--outputFolder', '-o',
                           help='Output folder of the files.',
                           required=False,
                           default='aggregatedFiles')
    parserOpt.add_argument('--writeFileNamesToFile', '-w',
                           help='',
                           default='aggregatedFilesBatch.txt')
    parserOpt.add_argument('--batchMode', '-bm',
                           help='The given file for --interactionFile and or --targetFile contain a list of the to be processed files.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument("--mergeBins", "-mb",
                           type=int,
                           default=0,
                           help="Merge neighboring interactions to one. The value of this parameter defines the maximum distance"
                           " a neighbor can have. The values are averaged.")
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    # parserOpt.add_argument('--backgroundModelFile', '-bmf',
    #                         help='path to the background file which is necessary to compute the rbz-score',
    #                         required=True)
    # parserOpt.add_argument('--range',
    #                         help='Defines the region upstream and downstream of a reference point which should be included. '
    #                         'Format is --region upstream downstream',
    #                         required=True,
    #                         type=int,
    #                         nargs=2)
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

# # def merge_neighbors_x_fold()
# def create_target_regions(pInteraction_file_data, pInteraction_file_data_1, pLoosePValue=None, pXFold=None):
#     # log.debug(pInteraction_file_data)
#     accepted_scores_file_1 = []
#     accepted_scores_file_2 = []

#     if pLoosePValue is not None:
#         # get significant regions
#         for key in pInteraction_file_data:
#             if float(pInteraction_file_data[key][-2]) <= pLoosePValue:
#                 accepted_scores_file_1.append(key)

#         for key in pInteraction_file_data_1:
#             if float(pInteraction_file_data_1[key][-2]) <= pLoosePValue:
#                 accepted_scores_file_2.append(key)
#     else:
#         for key in pInteraction_file_data:
#             if float(pInteraction_file_data[key][-1]) >= pXFold:
#                 accepted_scores_file_1.append(key)

#         for key in pInteraction_file_data_1:
#             if float(pInteraction_file_data_1[key][-1]) >= pXFold:
#                 accepted_scores_file_2.append(key)
#     # merge keys
#     accepted_scores_file_1.extend(accepted_scores_file_2)
#     accepted_scores_file_1 = np.unique(accepted_scores_file_1)

#     target_list = []
#     for key in accepted_scores_file_1:
#         # pInteraction_file_data
#         target_list.append(pInteraction_file_data[key][0:3])

#     # log.debug('target_list {}'.format(target_list))
#     return target_list


def filter_scores_target_list(pScoresDictionary, pTargetRegions):

    accepted_scores = {}
    for target in pTargetRegions:
        start = int(target[1])
        end = int(target[2])
        _accepted_scores = {}
        for key in pScoresDictionary:
            if int(pScoresDictionary[key][1]) >= start and int(pScoresDictionary[key][2]) <= end:
                _accepted_scores[key] = pScoresDictionary[key]
        # log.debug('_accepted_scores {}'.format(_accepted_scores))
        if len(_accepted_scores) > 0:

            values = np.array([0.0, 0.0, 0.0])
            for key in _accepted_scores:
                values += np.array(list(map(float,
                                            _accepted_scores[key][-3:])))
            keys_sorted = sorted(_accepted_scores.keys())
            _accepted_scores[keys_sorted[0]
                             ][-5] = _accepted_scores[keys_sorted[-1]][-5]
            _accepted_scores[keys_sorted[0]][-3] = values[0]
            _accepted_scores[keys_sorted[0]][-2] = values[1]
            _accepted_scores[keys_sorted[0]][-1] = values[2]

            accepted_scores[keys_sorted[0]] = _accepted_scores[keys_sorted[0]]
    # log.debug('accepted_scores {}'.format(len(accepted_scores)))
    return accepted_scores





def write(pOutFileName, pHeader, pNeighborhoods, pInteractionLines):

    # sum_of_interactions = float(pHeader.split('\t')[-1].split(' ')[-1])
    # log.debug('sum_of_interactions {}'.format(sum_of_interactions))
    with open(pOutFileName, 'w') as file:
        file.write('# Aggregated file, created with HiCExplorer\'s chicAggregateStatistic version {}\n'.format(__version__))
        file.write(pHeader)
        file.write(
            '#Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative distance\tRel Inter target\tRaw target')
        file.write('\n')
        
        if pNeighborhoods is not None:
            for data in pNeighborhoods:
                # log.debug('pInteractionLines[data] {}'.format(pInteractionLines[data]))
                new_line = '\t'.join(pInteractionLines[data][:6])
                # new_line += '\t' + format(float(sum_of_interactions), "10.5f")

                # new_line += '\t' + '\t'.join(format(float(x), "10.5f") for x in pInteractionLines[0][8:])
                new_line += '\t' + \
                    format(pNeighborhoods[data][-3], '10.5f') + \
                    '\t' + format(pNeighborhoods[data][-1], '10.5f')
                new_line += '\n'
                file.write(new_line)


def run_target_list_compilation(pInteractionFilesList, pArgs, pViewpointObj, pQueue=None):
    outfile_names = []
    for interactionFile in pInteractionFilesList:
        header, interaction_data, interaction_file_data = pViewpointObj.readInteractionFileForAggregateStatistics(
            pArgs.interactionFileFolder + '/' + interactionFile)

        target_regions = utilities.readBed(pArgs.targetFile)
        accepted_scores = filter_scores_target_list(
            interaction_file_data, target_regions)

        if len(accepted_scores) == 0:
            if pArgs.batchMode:
                with open('errorLog.txt', 'a+') as errorlog:
                    errorlog.write('Failed for: {} and {}.\n'.format(interactionFile[0], interactionFile[1]))
                    # break
            else:
                log.info('No target regions found')
                # break
        # outFileName = '.'.join(interactionFile.split('.')[:-1]) + '_' + pArgs.outFileNameSuffix
        outFileName = '.'.join(interactionFile.split('/')[-1].split('.')[:-1]) + '_' + pArgs.outFileNameSuffix

        if pArgs.batchMode:
            outfile_names.append(outFileName)
        outFileName = pArgs.outputFolder + '/' + outFileName

        if pArgs.mergeBins > 0:
            merged_neighborhood, _ = merge_neighbors(
                accepted_scores, pArgs.mergeBins)
            write(outFileName, header, merged_neighborhood,
                  interaction_file_data)
        else:
            write(outFileName, header, accepted_scores,
                  interaction_file_data)
    if pQueue is None:
        return
    pQueue.put(outfile_names)
    return


# def run_pvalue_compilation(pInteractionFilesList, pArgs, pViewpointObj, pQueue=None):
#     outfile_names = []
#     for interactionFile in pInteractionFilesList:

#         # header, interaction_data, interaction_file_data
#         data = [pViewpointObj.readInteractionFileForAggregateStatistics(
#             pArgs.interactionFileFolder + '/' + interactionFile[0])]
#         data.append(pViewpointObj.readInteractionFileForAggregateStatistics(
#             pArgs.interactionFileFolder + '/' + interactionFile[1]))

#         if pArgs.loosePValue is not None:
#             target_regions = create_target_regions(data[0][2], data[1][2], pLoosePValue=pArgs.pValue)
#         else:
#             target_regions = create_target_regions(data[0][2], data[1][2], pXFold=pArgs.xFold)

#         sample_prefix = interactionFile[0].split('/')[-1].split('_')[0] + '_' + interactionFile[1].split('/')[-1].split('_')[0]
#         for j in range(2):

#             accepted_scores = filter_scores_target_list(
#                 data[j][2], target_regions)
#             # log.debug('length of accepted_scores {}'.format(len(accepted_scores)))
#             if len(accepted_scores) == 0:
#                 if pArgs.batchMode:
#                     with open('errorLog.txt', 'a+') as errorlog:
#                         errorlog.write('Failed for: {} and {}.\n'.format(
#                             interactionFile[0], interactionFile[1]))
#                         # break
#                 else:
#                     log.info('No target regions found')
#                     # break
#             outFileName = '.'.join(interactionFile[j].split('/')[-1].split('.')[:-1]) + '_' + sample_prefix + pArgs.outFileNameSuffix

#             if pArgs.batchMode:
#                 outfile_names.append(outFileName)
#             outFileName = pArgs.outputFolder + '/' + outFileName

#             # for key in accepted_scores:
#             #     log.debug(accepted_scores[key])
#             #     exit()

#             # TODO change condition!
#             if pArgs.mergeBins > 0:
#                 merged_neighborhood = merge_neighbors(accepted_scores, pArgs.mergeBins)

#                 TODO  recompute p-values
#                 TODO redo line matching based on p-value
#                 write(outFileName, data[j][0],
#                       merged_neighborhood, data[j][2])
#             else:
#                 write(outFileName, data[j][0], accepted_scores, data[j][2])
#         # exit(0)
#     if pQueue is None:
#         return
#     pQueue.put(outfile_names)
#     return

def call_multi_core(pInteractionFilesList, pFunctionName, pArgs, pViewpointObj):
    outfile_names = [None] * pArgs.threads
    interactionFilesPerThread = len(pInteractionFilesList) // pArgs.threads
    all_data_collected = False
    queue = [None] * pArgs.threads
    process = [None] * pArgs.threads
    thread_done = [False] * pArgs.threads
    # log.debug('matrix read, starting processing')
    for i in range(pArgs.threads):

        if i < pArgs.threads - 1:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=pFunctionName, kwargs=dict(
            pInteractionFilesList=interactionFileListThread,
            pArgs=pArgs,
            pViewpointObj=pViewpointObj,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(pArgs.threads):
            if queue[i] is not None and not queue[i].empty():
                background_data_thread = queue[i].get()
                outfile_names[i] = background_data_thread
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
    return outfile_names

def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    # if args.pValue:
    #     # if args.backgroundModelFile is None or args.range is None:
    #     log.error('background model file and range need to be defined!')
    #         # exit(1)
    outfile_names = []
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if args.targetFile:
        # read all interaction files.
        if args.batchMode:
            interactionFileList = []
            with open(args.interactionFile[0], 'r') as interactionFile:
                file_ = True
                while file_:
                    file_ = interactionFile.readline().strip()
                    if file_ != '':
                        interactionFileList.append(file_)

            outfile_names = call_multi_core(interactionFileList, run_target_list_compilation, args, viewpointObj)
        else:
            interactionFileList = args.interactionFile
            run_target_list_compilation(interactionFileList, args, viewpointObj)

    # elif args.pValue:
    #     interactionFileList = []
    #     # log.debug('rbz')
    #     if args.batchMode:
    #         # log.debug('args.interactionFile {}'.format(args.interactionFile))
    #         with open(args.interactionFile[0], 'r') as interactionFile:

    #             file_ = True
    #             while file_:
    #                 # for line in fh.readlines():
    #                 file_ = interactionFile.readline().strip()
    #                 file2_ = interactionFile.readline().strip()
    #                 if file_ != '' and file2_ != '':
    #                     interactionFileList.append((file_, file2_))
            
    #         outfile_names = call_multi_core(interactionFileList, run_pvalue_compilation, args, viewpointObj)
            
    #     else:
    #         if len(args.interactionFile) % 2 == 0:
    #             i = 0
    #             while i < len(args.interactionFile):
    #                 interactionFileList.append(
    #                     (args.interactionFile[i], args.interactionFile[i + 1]))
    #                 i += 2
    #         else:
    #             log.error('Number of interaction files needs to be even: {}'.format(
    #                 len(args.interactionFile)))
    #             exit(1)
    #         outfile_names = run_pvalue_compilation(interactionFileList, args, viewpointObj)

    if args.batchMode:
        with open(args.writeFileNamesToFile, 'w') as nameListFile:
            nameListFile.write('\n'.join(outfile_names))
