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
from hicexplorer.lib import cnb
import h5py
import traceback

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
Per viewpoint the significant interactions are detected based on the background model. For each viewpoint file, an output file is created with all recorded significant interactions and
a target file. The target file is especially useful in the batch mode context; for two consecutive listed control and treatment viewpoints it merges the significant interactions which can then be used
to test for a differential interaction scheme.

chicSignificantInteractions supports two modes to detect significant interactions, either by an x-fold over the average background or a loose p-value. In both cases neighboring significant peaks are merged together and an additional
p-value is computed based on the sum of interactions for this neighborhood. Only interactions with a higher p-value (as specified by the threshold `--pValue`) are accepted as a significant interaction.

An example usage is for single mode is:

`$ chicSignificantInteractions --interactionFile interactionFilesFolder/Sox17_FL-E13-5_chr1_1000_2000.bed --referencePoints referencePointsFile.bed --range 20000 40000 --backgroundModelFile background_model.bed --loosePValue 0.5 --pValue 0.01`

An example usage is for batch mode is:

`$ chicViewpointBackgroundModel --matrices matrix1.cool matrix2.cool matrix3.cool --referencePoints referencePointsFile.bed --range 20000 40000 --outFileName background_model.bed`

The main difference between single mode and batch mode is that in single mode the parameter `--interactionFile` is interpreted as a list of viewpoint files created with `chicViewpoint`, whereas in batch mode only one file is allowed, containing the file names of viewpoint files (one per line).
This file is created by `chicViewpoint` and the parameter `--writeFileNamesToFile`. In batch mode, please remember to specify the folder (via `--interactionFileFolder`) where `chicViewpoint` wrote the files.
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction file (HDF5) which should be used for aggregation of the statistics.',
                                required=True)

    parserRequired.add_argument('--pValue', '-p',
                                help='p-value threshold to filter target regions for inclusion in differential analysis.',
                                # type=float,
                                required=True)
    parserMutuallyExclusiveGroupFilter = parser.add_mutually_exclusive_group(
        required=False)
    parserMutuallyExclusiveGroupFilter.add_argument('--xFoldBackground', '-xf',
                                                    help='Filter x-fold over background. Used to merge neighboring bins with a broader peak but '
                                                    'less significant interactions to a single peak with high significance. Used only for pValue option.'
                                                    # type=float

                                                    )
    parserMutuallyExclusiveGroupFilter.add_argument('--loosePValue', '-lp',
                                                    help='loose p-value threshold to filter target regions in a first round. '
                                                    'Used to merge neighboring bins with a broader peak but less significant interactions to a single peak with high significance.'
                                                    ' Used only for pValue option.'
                                                    # type=float
                                                    )
    parserRequired.add_argument('--backgroundModelFile', '-bmf',
                                help='path to the background file.',
                                required=True)
    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream, e.g. --region 500000 500000 plots 500kb up- and 500kb downstream. This value should not exceed the range used in the other chic-tools.',
                                required=True,
                                type=int,
                                nargs=2)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='File name suffix to save the results; prefix is the input file name'
                           ' (Default: %(default)s).',
                           required=False,
                           default='significant_interactions.hdf5')
    parserOpt.add_argument('--targetFileName', '-tl',
                           help='The file to store the target file names'
                           ' (Default: %(default)s).',
                           default='targetFile.hdf5')
    parserOpt.add_argument('--combinationMode',
                        '-a',
                        help='This option defines how the interaction data should be computed and combined: '
                        'dual: Combines as follows: [[matrix1_gene1, matrix2_gene1], [matrix2_gene1, matrix3_gene1],[matrix1_gene2, matrix2_gene2], ...]'
                        'single: Combines as follows: [matrix1_gene1, matrix1_gene2, matrix2_gene1, ...], '
                        'allGenes: Combines as follows: [[matrix1_gene1, matrix2_gene1, matrix2_gene1], [matrix1_gene2, matrix2_gene2, matrix3_gene2], ...]'
                        'oneGene: Computes all data of one gene, please specify \'--\'. If a gene is not unique, each viewpoint is treated independently.'
                        'file: A user specific file (\'--\') with one entry per line and tab seperated: matrixName   geneName. Please define how many lines should be combined to one (\'--computeSampleNumber\').'
                        ' (Default: %(default)s).',
                        default='dual',
                        choices=['dual', 'single', 'allGenes', 'oneGene', 'file']
                        )

    parserOpt.add_argument('--combinationName', '-cn',
                           help='Gene name or file name for modes \'oneGene\' or \'file\' of parameter \'--combinationMode\''
                           ' (Default: %(default)s).',
                           required=False,
                           default=None)
    # parserOpt.add_argument('--targetFolder', '-tf',
    #                        help='Folder where the target files are stored'
    #                        ' (Default: %(default)s).',
    #                        required=False,
    #                        default='targetFolder')
    # parserOpt.add_argument('--outputFolder', '-o',
    #                        help='Output folder of the significant interaction files'
    #                        ' (Default: %(default)s).',
    #                        required=False,
    #                        default='significantFiles')
    # parserOpt.add_argument('--writeFileNamesToFile', '-w',
    #                        help='(Default: %(default)s).',
    #                        default='significantFilesBatch.txt')
    # parserOpt.add_argument('--targetFileList', '-tl',
    #                        help='The file to store the target file names'
    #                        ' (Default: %(default)s).',
    #                        default='targetList.txt')
    # parserOpt.add_argument('--batchMode', '-bm',
    #                        help='Turn on batch mode. The given file for --interactionFile and or --targetFile contain a list of the to be processed files.',
    #                        required=False,
    #                        action='store_true')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--truncateZeroPvalues', '-tzpv',
                           help='Sets all p-values which are equal to zero to one. This has the effect that the associated positions are not part of the significance decision.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate range of backgroundmodel starting at distance x. E.g. all values greater than 500kb are set to the value of the 500kb bin'
                           ' (Default: %(default)s).',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument('--peakInteractionsThreshold', '-pit',
                           type=int,
                           default=5,
                           help='The minimum number of interactions a detected peak needs to have to be considered'
                           ' (Default: %(default)s).')

    # parserOpt.add_argument('--resolution', '-r',
    #                        help='Resolution of the bin in genomic units. Values are set as number of bases, e.g. 1000 for a 1kb, 5000 for a 5kb or 10000 for a 10kb resolution.'
    #                        'This value is used to merge neighboring bins'
    #                        ' (Default: %(default)s).',
    #                        type=int,
    #                        default=1000,
    #                        required=False)

    parserOpt.add_argument('--computeSampleNumber', '-csn',
                           help='Number of samples to compute together. Applies only in batch mode'
                           ' (Default: %(default)s).',
                           required=False,
                           default=2,
                           type=int)
    parserOpt.add_argument('--multipleTesting', '-mt',
                           help='Multiple testing correction per relative distance with Bonferroni or FDR'
                           ' (Default: %(default)s).',
                           type=str,
                           default="None",
                           choices=['fdr', 'bonferroni', 'None'],
                           )
    parserOpt.add_argument('--thresholdMultipleTesting', '-tmt',
                           help='Threshold for Bonferroni / FDR. Either a float value for all or a file with one threshold per relative distance.'
                           )
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution, pQueue=None):
    outfile_names = []
    target_outfile_names = []
    try:
        for interactionFile in pInteractionFilesList:
            target_list = []
            sample_prefix = ''
            log.debug('interactionFile {}'.format(interactionFile))
            for sample in interactionFile:
                log.debug('sample {}'.format(sample))

                # header,
                # interaction_data:rel interaction, p-value, raw, x-fold::{-1000:[0.1, 0.01, 2.3, 5]},
                # if pArgs.interactionFileFolder != '.':
                #     absolute_sample_path = pArgs.interactionFileFolder + '/' + sample
                # else:
                #     absolute_sample_path = sample
                data = pViewpointObj.readInteractionFileForAggregateStatistics(pFilePath, sample)
                log.debug(data)
                # exit(1)
                # sample_prefix += sample.split('/')[-1].split('_')[0]
                # sample_prefix += '_'
                # filter by x-fold over background value or loose p-value
                # and merge neighbors. Use center position to compute new p-value.
                log.debug('209')

                if pArgs.xFoldBackground is not None:
                    accepted_scores, merged_lines_dict = merge_neighbors_x_fold(
                        pArgs.xFoldBackground, data, pViewpointObj, pResolution=pResolution)
                else:
                    accepted_scores, merged_lines_dict = merge_neighbors_loose_p_value(
                        pArgs.loosePValue, data, pViewpointObj, pResolution=pResolution, pTruncateZeroPvalues=pArgs.truncateZeroPvalues)
                log.debug('217')

                # compute new p-values and filter by them
                accepted_scores, target_lines = compute_new_p_values(
                    accepted_scores, pBackground, pArgs.pValue, merged_lines_dict, pArgs.peakInteractionsThreshold, pViewpointObj)

                # filter by new p-value
                if len(accepted_scores) == 0:
                    if pArgs.batchMode:
                        with open('errorLog.txt', 'a+') as errorlog:
                            errorlog.write('Failed for: {} and {}.\n'.format(
                                interactionFile[0], interactionFile[1]))
                    else:
                        log.info('No target regions found')
                log.debug('231')
                
                # outFileName = '.'.join(sample.split(
                #     '/')[-1].split('.')[:-1]) + '_' + pArgs.outFileNameSuffix
                # if pArgs.batchMode:
                #     outfile_names.append(outFileName)
                # outFileName = pArgs.outputFolder + '/' + outFileName

                ##### TODO write significant and target files
                # write only significant lines to file
                # write(outFileName, data[0], accepted_scores)
                # target_list.append(target_lines)
            log.debug('243')

            target_list = [item for sublist in target_list for item in sublist]
            log.debug('interactionFile {}'.format(interactionFile))
            # sample_name = '_'.join(interactionFile[0].split('/')[-1].split('.')[0].split('_')[1:])
            # target_name = sample_prefix + sample_name + '_target.txt'
            # target_outfile_names.append(target_name)
            # target_name = pArgs.targetFolder + '/' + target_name
            # writeTargetList(target_list, target_name, pArgs)
        if pQueue is None:
            return target_outfile_names
    except Exception as exp:
        if pQueue is None:
            return 'Fail: ' + str(exp) + traceback.format_exc()
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc()) 
        return
    pQueue.put([outfile_names, target_outfile_names])
    return


def compute_new_p_values(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold, pViewpointObj):
    accepted = {}
    accepted_lines = []
    if isinstance(pPValue, float):
        for key in pData:
            if key in pBackgroundModel:
                log.debug('Recompute p-values. Old: {}'.format(pData[key][-3]))
                pData[key][-3] = 1 - cnb.cdf(float(pData[key][-1]), float(pBackgroundModel[key][0]), float(pBackgroundModel[key][1]))
                log.debug('new {}\n\n'.format(pData[key][-3]))
                if pData[key][-3] <= pPValue:
                    if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                        accepted[key] = pData[key]
                        target_content = pMergedLinesDict[key][0][:3]
                        target_content[2] = pMergedLinesDict[key][-1][2]
                        accepted_lines.append(target_content)
            else:
                log.debug('key not in background')
    elif isinstance(pPValue, dict):
        for key in pData:
            if key in pBackgroundModel:
                log.debug('Recompute p-values. Old: {}'.format(pData[key][-3]))

                pData[key][-3] = 1 - cnb.cdf(float(pData[key][-1]), float(pBackgroundModel[key][0]), float(pBackgroundModel[key][1]))
                log.debug('new {}\n\n'.format(pData[key][-3]))

                if pData[key][-3] <= pPValue[key]:
                    if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                        accepted[key] = pData[key]
                        target_content = pMergedLinesDict[key][0][:3]
                        target_content[2] = pMergedLinesDict[key][-1][2]
                        accepted_lines.append(target_content)
            else:
                log.debug('key not in background')
    return accepted, accepted_lines


def merge_neighbors_x_fold(pXfold, pData, pViewpointObj, pResolution):
    accepted = {}
    accepted_line = {}
    if isinstance(pXfold, float):
        for key in pData[0]:

            if pData[0][key][-1] < pXfold:
                continue
            accepted[key] = pData[0][key]
            accepted_line[key] = pData[1][key]
    elif isinstance(pXfold, dict):
        for key in pData[0]:
            if pData[0][key][-1] < pXfold[key]:
                continue
            accepted[key] = pData[0][key]
            accepted_line[key] = pData[1][key]
    if accepted_line:
        return pViewpointObj.merge_neighbors(accepted_line, pMergeThreshold=pResolution)
    return accepted_line, None


def merge_neighbors_loose_p_value(pLoosePValue, pData, pViewpointObj, pResolution, pTruncateZeroPvalues):
    accepted = {}
    accepted_line = {}
    if isinstance(pLoosePValue, float):
        for key in pData[0]:
            if pTruncateZeroPvalues:
                if pData[0][key][1] == 0 or pData[0][key][1] > pLoosePValue:
                    continue
            else:
                if pData[0][key][1] > pLoosePValue:
                    continue
            accepted[key] = pData[0][key]
            accepted_line[key] = pData[1][key]
    elif isinstance(pLoosePValue, dict):
        for key in pData[0]:
            if pTruncateZeroPvalues:
                if pData[0][key][1] == 0 or pData[0][key][1] > pLoosePValue[key]:
                    continue
            else:
                if pData[0][key][1] > pLoosePValue[key]:
                    continue
            accepted[key] = pData[0][key]
            accepted_line[key] = pData[1][key]
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


def call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution):
    outfile_names = [None] * pArgs.threads
    target_list_name = [None] * pArgs.threads
    interactionFilesPerThread = len(pInteractionFilesList) // pArgs.threads
    all_data_collected = False
    queue = [None] * pArgs.threads
    process = [None] * pArgs.threads
    thread_done = [False] * pArgs.threads

    fail_flag = False
    fail_message = ''

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
            pBackground=pBackground,
            pFilePath=pFilePath,
            pResolution=pResolution,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(pArgs.threads):
            if queue[i] is not None and not queue[i].empty():
                background_data_thread = queue[i].get()
                if 'Fail:' in background_data_thread:
                    fail_flag = True
                    fail_message = background_data_thread[6:]
                else:
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
    if fail_flag:
        log.error(fail_message)
        exit(1)

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


def read_threshold_file(pFile):
    distance_value_dict = {}
    with open(pFile, 'r') as file:
        file_ = True
        while file_:
            line = file.readline().strip()
            if line.startswith('#'):
                continue
            if line == '':
                break
            relative_distance, value = line.split('\t')
            distance_value_dict[int(relative_distance)] = float(value)
    return distance_value_dict


def main(args=None):

    args = parse_arguments().parse_args(args)
    # args.p_value_dict = None
    # args.p_loose_value_dict = None
    # args.x_fold_dict = None
    # if not os.path.exists(args.outputFolder):
    #     try:
    #         os.makedirs(args.outputFolder)
    #     except OSError as exc:  # Guard against race condition
    #         if exc.errno != errno.EEXIST:
    #             raise
    # if not os.path.exists(args.targetFolder):
    #     try:
    #         os.makedirs(args.targetFolder)
    #     except OSError as exc:  # Guard against race condition
    #         if exc.errno != errno.EEXIST:
    #             raise

    if args.pValue:
        try:
            args.pValue = float(args.pValue)
        except Exception:
            args.pValue = read_threshold_file(args.pValue)
    if args.loosePValue:
        try:
            args.loosePValue = float(args.loosePValue)
        except Exception:
            args.loosePValue = read_threshold_file(args.loosePValue)
    if args.xFoldBackground:
        try:
            args.xFoldBackground = float(args.xFoldBackground)
        except Exception:
            args.xFoldBackground = read_threshold_file(args.xFoldBackground)

    viewpointObj = Viewpoint()
    outfile_names = []

    interactionFileList = []

    background_model = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range, args.fixateRange)
    
    ### read hdf file
    interactionFileHDF5Object = h5py.File(args.interactionFile, 'r')
    keys_interactionFile = list(interactionFileHDF5Object.keys()) 
    # if 'matrices' not in keys_interactionFile or 'referencePoints' not in keys_interactionFile:
    #     log.error('Given interaction file has the wrong structure! Please use an interaction file created with HiCExplorer 3.6, earlier version are not supported!')
    #     exit(1)

    log.debug('list(interactionFileHDF5Object.keys()) {}'.format(list(interactionFileHDF5Object.keys()) ))
    # log.debug('list(interactionFileHDF5Object[keys_interactionFile[0]].keys()) {}'.format(list(interactionFileHDF5Object[keys_interactionFile[0]].keys()) ))
    # log.debug('list(interactionFileHDF5Object[keys_interactionFile[1]].keys()) {}'.format(list(interactionFileHDF5Object[keys_interactionFile[1]].keys()) ))

    # interaction_file_keys = interactionFileHDF5Object[keys_interactionFile[0]].keys()
    # matrices_keys = list(interactionFileHDF5Object[keys_interactionFile[1]].keys())


    # referencePoints_keys = interactionFileHDF5Object.get(keys_interactionFile[2]).values
    # log.debug(referencePoints_keys)
    # log.debug('interactionFileHDF5Object[keys_interactionFile[1]][i] {}'.format(interactionFileHDF5Object[keys_interactionFile[1]][matrices_keys[0]]))
# ['mydataset']



    # 'dual', 'single', 'allGenes', 'oneGene', 'file']

    interactionList = []
    if args.combinationMode == 'dual':
        if len(keys_interactionFile) > 1:
            for i, sample in enumerate(keys_interactionFile):
                for sample2 in keys_interactionFile[i + 1:]:
                    
                    matrix_obj1 = interactionFileHDF5Object[sample]
                    matrix_obj2 = interactionFileHDF5Object[sample]

                    chromosomeList1 = sorted(list(matrix_obj1.keys()))
                    chromosomeList2 = sorted(list(matrix_obj2.keys()))
                    chromosomeList1.remove('genes')
                    chromosomeList2.remove('genes')
                    for chromosome1, chromosome2 in zip(chromosomeList1, chromosomeList2):
                        geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                        geneList2 = sorted(list(matrix_obj2[chromosome2].keys()))

                        for gene1, gene2 in zip(geneList1, geneList2):
                            interactionList.append([[sample,chromosome1, gene1],[sample2,chromosome2, gene2]])

                    # for viewpoint, viewpoint2 in zip(sample, sample2):
                    #     writeFileNamesToList.append(viewpoint.encode("ascii", "ignore"))
                    #     writeFileNamesToList.append(viewpoint2.encode("ascii", "ignore"))
            # log.debug(interactionList)
        else:
            log.error('Dual mode selected but only one matrix is stored')
    elif args.combinationMode == 'single':
        for i, sample in enumerate(keys_interactionFile):
                
                matrix_obj1 = interactionFileHDF5Object[sample]
                chromosomeList1 = sorted(list(matrix_obj1.keys()))
                chromosomeList1.remove('genes')
                for chromosome1 in chromosomeList1:
                    geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                    for gene1 in geneList1:
                        interactionList.append([sample,chromosome1, gene1])

    elif args.combinationMode == 'allGenes':

        if len(keys_interactionFile) > 0:
            matrix_obj1 = interactionFileHDF5Object[keys_interactionFile[0]]
            gene_list = matrix_obj1['genes'].keys()

            for gene in gene_list:
                gene_list = []
                for matrix in keys_interactionFile:
                    gene_list.append([matrix, 'genes', gene])
                interactionList.append(gene_list)

    elif args.combinationMode == 'oneGene': 
        if len(keys_interactionFile) > 0:
            matrix_obj1 = interactionFileHDF5Object[keys_interactionFile[0]]
            all_detected = False
            counter = 0
            gene_list = []
            while not all_detected:
                if counter == 0:
                    check_gene_name = args.combinationName
                else:
                    check_gene_name = args.combinationName + '_' + str(counter)

                if check_gene_name in matrix_obj1['genes']:
                    gene_list.append(check_gene_name)
                else:
                    all_detected = True

                counter += 1

            for gene in gene_list:
                gene_list = []
                for matrix in keys_interactionFile:
                    gene_list.append([matrix, 'genes', gene])
                interactionList.append(gene_list)

            
    elif args.combinationMode == 'file':
        pass
    
    resolution = interactionFileHDF5Object.attrs['resolution'][()]
    interactionFileHDF5Object.close()
    
    outfile_names, target_list_name = call_multi_core(
                            interactionList, args, viewpointObj, background_model, args.interactionFile, resolution)


    # if args.batchMode:
    #     if len(list(matrices_keys)) <= args.computeSampleNumber:
    #         for referencePoint in referencePoints_keys:
    #             lines = []
    #             for i in range(0, args.computeSampleNumber):
    #                 lines.append(interactionFileHDF5Object[keys_interactionFile[1]][i])
    #                 for j in range(i, args.computeSampleNumber):
    #                     lines.append()
    #     else:
    #         log.error('{} samples are requested to considered together, but only {} different samples are available.'.format(args.computeSampleNumber, len(list(matrices_keys))))
    #         exit(1)

    #     with open(args.interactionFile[0], 'r') as interactionFile:
    #         file_ = True
    #         while file_:
    #             lines = []
    #             for i in range(0, args.computeSampleNumber):
    #                 file_ = interactionFile.readline().strip()
    #                 if file_ != '':
    #                     lines.append(file_)
    #             if len(lines) > 0:
    #                 interactionFileList.append(lines)
    #     log.debug('interactionFileList {}'.format(interactionFileList))
    #     outfile_names, target_list_name = call_multi_core(
    #         interactionFileList, args, viewpointObj, background_model)

    # else:
    #     i = 0
    #     while i < len(args.interactionFile):
    #         lines = []
    #         for j in range(0, args.computeSampleNumber):
    #             if i < len(args.interactionFile):
    #                 lines.append(args.interactionFile[i])
    #             i += 1
    #         interactionFileList.append(lines)

    #     target_list_name = compute_interaction_file(
    #         interactionFileList, args, viewpointObj, background_model)
    #     if 'Fail: ' in target_list_name:
    #         log.error(target_list_name[6:])

    # if args.batchMode:
    #     with open(args.writeFileNamesToFile, 'w') as nameListFile:
    #         nameListFile.write('\n'.join(outfile_names))

    #     with open(args.targetFileList, 'w') as targetNamesFile:
    #         targetNamesFile.write('\n'.join(target_list_name))
