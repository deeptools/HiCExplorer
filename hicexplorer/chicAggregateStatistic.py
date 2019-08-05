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

from intervaltree import IntervalTree, Interval
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

    parserRequired.add_argument('--targetFile', '-tf',
                                help='path to the target files which contains the target regions to prepare data for differential analysis.'
                                )

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileNameSuffix', '-suffix',
                           help='File name suffix to save the result.',
                           required=False,
                           default='_aggregate_target.bed')

    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--targetFileFolder', '-tff',
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
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )

    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def filter_scores_target_list(pScoresDictionary, pTargetList):

    accepted_scores = {}
    same_target_dict = {}
    target_regions = utilities.readBed(pTargetList)
    if len(target_regions) == 0:
        return accepted_scores
        
    hicmatrix = hm.hiCMatrix()
    target_regions_intervaltree = hicmatrix.intervalListToIntervalTree(target_regions)[0]
    for key in pScoresDictionary:
        # try:
        chromosome = pScoresDictionary[key][0]
        start = int(pScoresDictionary[key][1])
        end = int(pScoresDictionary[key][2])
        target_interval = target_regions_intervaltree[chromosome][start:end]
        if target_interval:
            target_interval = sorted(target_interval)[0]
            if target_interval in same_target_dict:
                same_target_dict[target_interval].append(key)
            else:
                same_target_dict[target_interval] = [key]
 
    for target in same_target_dict:

        values = np.array([0.0, 0.0, 0.0])
        same_target_dict[target] = sorted(same_target_dict[target])
        
        for key in same_target_dict[target]:
            values += np.array(list(map(float, pScoresDictionary[key][-3:])))
        new_data_line = pScoresDictionary[same_target_dict[target][0]]
        new_data_line[2] = pScoresDictionary[same_target_dict[target][-1]][2]
        new_data_line[-5] = pScoresDictionary[same_target_dict[target][-1]][-5]
        new_data_line[-3] = values[0]
        new_data_line[-2] = values[1]
        new_data_line[-1] = values[2]

        accepted_scores[same_target_dict[target][0]] = new_data_line

    return accepted_scores




def write(pOutFileName, pHeader, pNeighborhoods, pInteractionLines):

    with open(pOutFileName, 'w') as file:
        file.write('# Aggregated file, created with HiCExplorer\'s chicAggregateStatistic version {}\n'.format(__version__))
        file.write(pHeader)
        file.write(
            '#Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative distance\tRaw target')
        file.write('\n')
        
        if pNeighborhoods is not None:
            for data in pNeighborhoods:
                new_line = '\t'.join(pInteractionLines[data][:6])
                new_line += '\t' +  format(pNeighborhoods[data][-1], '10.5f')
                new_line += '\n'
                file.write(new_line)


def run_target_list_compilation(pInteractionFilesList, pTargetList, pArgs, pViewpointObj, pQueue=None):
    outfile_names = []
    for i, interactionFile in enumerate(pInteractionFilesList):
        for sample in interactionFile:
            header, interaction_data, interaction_file_data = pViewpointObj.readInteractionFileForAggregateStatistics(
                pArgs.interactionFileFolder + '/' + sample)

            if pArgs.batchMode:
                target_file = pArgs.targetFileFolder+'/'+ pTargetList[i]
            else:
                target_file = pTargetList[i]

            accepted_scores = filter_scores_target_list(interaction_file_data, target_file)

            if len(accepted_scores) == 0:
                # do not call 'break' or 'continue'
                # with this an empty file is written and no track of 'no significant interactions' detected files needs to be recorded.
                if pArgs.batchMode:
                    with open('errorLog.txt', 'a+') as errorlog:
                        errorlog.write('Failed for: {} and {}.\n'.format(interactionFile[0], interactionFile[1]))
                else:
                    log.info('No target regions found')
            outFileName = '.'.join(sample.split('/')[-1].split('.')[:-1]) + '_' + pArgs.outFileNameSuffix

            if pArgs.batchMode:
                outfile_names.append(outFileName)
            outFileName = pArgs.outputFolder + '/' + outFileName

            write(outFileName, header, accepted_scores,
                    interaction_file_data)
    if pQueue is None:
        return
    pQueue.put(outfile_names)
    return


def call_multi_core(pInteractionFilesList, pTargetFileList, pFunctionName, pArgs, pViewpointObj):
    outfile_names = [None] * pArgs.threads
    interactionFilesPerThread = len(pInteractionFilesList) // pArgs.threads
    all_data_collected = False
    queue = [None] * pArgs.threads
    process = [None] * pArgs.threads
    thread_done = [False] * pArgs.threads
    for i in range(pArgs.threads):

        if i < pArgs.threads - 1:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
            targetFileListThread = pTargetFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:]
            targetFileListThread = pTargetFileList[i * interactionFilesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=pFunctionName, kwargs=dict(
            pInteractionFilesList=interactionFileListThread,
            pTargetList=targetFileListThread,
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
   outfile_names = []
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
 
    interactionFileList = []
    targetFileList = []

    if args.batchMode:
        with open(args.interactionFile[0], 'r') as interactionFile:
            file_ = True
            while file_:
                file_ = interactionFile.readline().strip()
                file2_ = interactionFile.readline().strip()
                if file_ != '' and file2_ != '':
                    interactionFileList.append((file_, file2_))
        with open(args.targetFile, 'r') as targetFile:
            file_ = True
            while file_:
                file_ = targetFile.readline().strip()
                if file_ != '':
                    targetFileList.append(file_)
        outfile_names = call_multi_core(interactionFileList, targetFileList, run_target_list_compilation, args, viewpointObj)
        
    else:
        targetFileList = args.targetFile
        if len(args.interactionFile) % 2 == 0:
            i = 0
            while i < len(args.interactionFile):
                interactionFileList.append(
                    (args.interactionFile[i], args.interactionFile[i + 1]))
                i += 2
        else:
            log.error('Number of interaction files needs to be even: {}'.format(
                len(args.interactionFile)))
            exit(1)
        run_target_list_compilation(interactionFileList, targetFileList, args, viewpointObj)


    if args.batchMode:
        with open(args.writeFileNamesToFile, 'w') as nameListFile:
            nameListFile.write('\n'.join(outfile_names))
