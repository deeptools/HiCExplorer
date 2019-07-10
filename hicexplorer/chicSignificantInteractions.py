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

    parserRequired.add_argument('--pValue', '-p',
                                help='p-value threshold value to filter target regions to include them for differential analysis.',
                                type=float)
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
                           help='File name suffix to save the result.',
                           required=False,
                           default='_significant_interactions.bed')

    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--outputFolder', '-o',
                           help='Output folder of the files.',
                           required=False,
                           default='significantFiles')
    parserOpt.add_argument('--writeFileNamesToFile', '-w',
                           help='',
                           default='significantFilesBatch.txt')
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
                           default=7,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackgroundSumOfDensities, pQueue=None):
    outfile_names = []
    for interactionFile in pInteractionFilesList:

        # header, 
        # interaction_data:rel interaction, p-value, raw, x-fold::{-1000:[0.1, 0.01, 2.3, 5]},
        # interaction_file_data: 
        data = pViewpointObj.readInteractionFileForAggregateStatistics(pArgs.interactionFileFolder + '/' + interactionFile)

        sample_prefix = interactionFile.split('/')[-1].split('_')[0]
        
        # filter by x-fold over background value or loose p-value
        # and merge neighbors. Use center position to compute new p-value.
        if pArgs.xFoldBackground is not None:
           accepted_scores = merge_neighbors_x_fold(pArgs.xFoldBackground, data, pViewpointObj, pResolution=1000)
        else:
           accepted_scores = merge_neighbors_loose_p_value(pArgs.loosePValue, data, pViewpointObj, pResolution=1000)
        
        log.debug('data: {}'.format(len(data)))
        # log.debug('data: {}'.format(data))

        # compute new p-values
        accepted_scores = compute_new_p_values(accepted_scores, pBackgroundSumOfDensities, pArgs.pValue)

        # filter by new p-value
        if len(accepted_scores) == 0:
            if pArgs.batchMode:
                with open('errorLog.txt', 'a+') as errorlog:
                    errorlog.write('Failed for: {} and {}.\n'.format(
                        interactionFile[0], interactionFile[1]))
                    # break
            else:
                log.info('No target regions found')
                # break
        outFileName = '.'.join(interactionFile.split('/')[-1].split('.')[:-1]) + '_' + sample_prefix + pArgs.outFileNameSuffix

        if pArgs.batchMode:
            outfile_names.append(outFileName)
        outFileName = pArgs.outputFolder + '/' + outFileName

        write(outFileName, data[0], accepted_scores)
        # exit(0)
    if pQueue is None:
        return
    pQueue.put(outfile_names)
    return

def compute_new_p_values(pData, pBackgroundSumOfDensities, pPValue):
    # log.debug('pBackgroundSumOfDensities {}'.format(pBackgroundSumOfDensities))
    accepted = {}
    for key in pData:
        if int(float(pData[key][-1])) - 1 < 0:
            pData[key][-3] = pBackgroundSumOfDensities[key][0]
        else:
            pData[key][-3] = 1 - pBackgroundSumOfDensities[key][int(float(pData[key][-1]))]
        
        if pData[key][-3] <= pPValue:
            accepted[key] = pData[key]
    # log.debug(pData)
    return accepted
    

def merge_neighbors_x_fold(pXfold, pData, pViewpointObj, pResolution):
    accepted = {}
    accepted_line = {}
    log.debug('len(pData[1]) {}'.format(len(pData[1])))

    # log.debug('pData[1] {}'.format(pData[1]))
    for key in pData[1]:
        # log.debug('pData[1][key] {}'.format(pData[1][key]))
        # log.debug('pData[2][key] {}'.format(pData[2][key]))
        
        if pData[1][key][-1] < pXfold:
            continue
        accepted[key] = pData[1][key]
        accepted_line[key] = pData[2][key]
    log.debug('len(accepted_line) {}'.format(len(accepted_line)))

    return pViewpointObj.merge_neighbors(accepted_line, pMergeThreshold=pResolution)

def merge_neighbors_loose_p_value(pLoosePValue, pData, pViewpointObj, pResolution):
    accepted = {}
    accepted_line = {}
    log.debug('len(pData[1]) {}'.format(len(pData[1])))

    # log.debug('pData[1] {}'.format(pData[1]))
    for key in pData[1]:
        # log.debug('pData[1][key] {}'.format(pData[1][key]))
        # log.debug('pData[2][key] {}'.format(pData[2][key]))
        
        if pData[1][key][1] > pLoosePValue:
            continue
        accepted[key] = pData[1][key]
        accepted_line[key] = pData[2][key]
    log.debug('len(accepted_line) {}'.format(len(accepted_line)))

    return pViewpointObj.merge_neighbors(accepted_line, pMergeThreshold=pResolution)

def write(pOutFileName, pHeader, pInteractionLines):

    # sum_of_interactions = float(pHeader.split('\t')[-1].split(' ')[-1])
    # log.debug('sum_of_interactions {}'.format(sum_of_interactions))
    with open(pOutFileName, 'w') as file:
        file.write(pHeader)
        file.write(
            '#Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative interactions\tp-value\tx-fold\tRaw target')
        file.write('\n')
        
        
        for data in pInteractionLines:
            # log.debug('pInteractionLines[data] {}'.format(pInteractionLines[data]))
            new_line = '\t'.join(pInteractionLines[data][:6])
            # new_line += '\t' + format(float(sum_of_interactions), "10.5f")

            new_line += '\t' + '\t'.join(format(float(x), "0.10f") for x in pInteractionLines[data][6:])
            # new_line += '\t' + \
            #     format(pNeighborhoods[data][-3], '10.5f') + \
            #     '\t' + format(pNeighborhoods[data][-1], '10.5f')
            new_line += '\n'
            file.write(new_line)

def call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackgroundSumOfDensities):
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
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    viewpointObj = Viewpoint()
    outfile_names = []

    interactionFileList = []

    # if args.loosePValue:
    background_model = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range)
    
    # log.debug('background_model {}'.format(background_model))
    # log.debug('compute sum of densities')
    background_sum_of_densities_dict = viewpointObj.computeSumOfDensities(background_model, args, pXfoldMaxValue=args.xFoldMaxValueNB)
    # else:
    background_model_mean_values = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range, pMean=True)
        # log.debug('rbz')
    if args.batchMode:
        # log.debug('args.interactionFile {}'.format(args.interactionFile))
        with open(args.interactionFile[0], 'r') as interactionFile:

            file_ = True
            while file_:
                # for line in fh.readlines():
                file_ = interactionFile.readline().strip()
                # file2_ = interactionFile.readline().strip()
                if file_ != '':
                    interactionFileList.append(file_)
        
        outfile_names = call_multi_core(interactionFileList, args, viewpointObj, background_sum_of_densities_dict)
        
    else:
        # if len(args.interactionFile) % 2 == 0:
        i = 0
        while i < len(args.interactionFile):
            interactionFileList.append(args.interactionFile[i])
            i += 1
        # else:
        #     log.error('Number of interaction files needs to be even: {}'.format(
        #         len(args.interactionFile)))
        #     exit(1)
        # pInteractionFilesList, pArgs, pViewpointObj, pBackgroundSumOfDensities,
        outfile_names = compute_interaction_file(interactionFileList, args, viewpointObj, background_sum_of_densities_dict)

    if args.batchMode:
        with open(args.writeFileNamesToFile, 'w') as nameListFile:
            nameListFile.write('\n'.join(outfile_names))