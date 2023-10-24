import argparse
import sys
import errno
import os
import math
from multiprocessing import Process, Queue
import time
import traceback
import logging
log = logging.getLogger(__name__)

import numpy as np
import h5py

from intervaltree import IntervalTree, Interval
import hicmatrix.HiCMatrix as hm

from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
chicAggregateStatistic is a preprocessing tool for chicDifferentialTest. It takes two consecutive viewpoint files and one target file and creates one
file containing all locations which should be tested for differential interactions. Either one target file for two consecutive viewpoint files or one
target file for all viewpoints is accepted.
"""
                                     )
    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for aggregation of the statistics.',
                                required=True)

    parserRequired.add_argument('--targetFile', '-tf',
                                help='path to the target files which contains the target regions to prepare data for differential analysis. This is either the target file in the hdf format created by chicSignificantInteractions or a regular, three column bed file.'
                                )

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the result'
                           ' (Default: %(default)s).',
                           required=False,
                           default='aggregate_target.hdf')

    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)ist'
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


def filter_scores_target_list(pScoresDictionary, pTargetList=None, pTargetPosDict=None, pTargetFile=None):
    isBED = True
    accepted_scores = {}
    same_target_dict = {}
    target_regions_intervaltree = None

    if pTargetList is not None:
        if pTargetPosDict is not None:
            chromosome = pTargetPosDict[pTargetList[-1]]['chromosome']
            start_list = pTargetPosDict[pTargetList[-1]]['start_list']
            end_list = pTargetPosDict[pTargetList[-1]]['end_list']
        else:
            # read hdf content for this specific combination
            targetFileHDF5Object = h5py.File(pTargetFile, 'r')
            target_object = targetFileHDF5Object['/'.join(pTargetList)]
            chromosome = target_object.get('chromosome')[()].decode("utf-8")
            start_list = list(target_object['start_list'][:])
            end_list = list(target_object['end_list'][:])
            targetFileHDF5Object.close()

        chromosome = [chromosome] * len(start_list)
        target_regions = list(zip(chromosome, start_list, end_list))

        if len(target_regions) == 0:
            return accepted_scores

        hicmatrix = hm.hiCMatrix()
        target_regions_intervaltree = hicmatrix.intervalListToIntervalTree(target_regions)[0]
    else:
        log.error('No target list given.')
        raise Exception('No target list given.')
    
    for key in pScoresDictionary:
        chromosome = pScoresDictionary[key][0]
        start = int(pScoresDictionary[key][1])
        end = int(pScoresDictionary[key][2])
        if chromosome in target_regions_intervaltree:
            target_interval = target_regions_intervaltree[chromosome][start:end]
        else:
            continue
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


def writeAggregateHDF(pOutFileName, pOutfileNamesList, pAcceptedScoresList, pArgs):

    aggregateFileH5Object = h5py.File(pOutFileName, 'w')

    aggregateFileH5Object.attrs['type'] = "aggregate"
    aggregateFileH5Object.attrs['version'] = __version__

    counter = 0
    for key_outer, data_outer in zip(pOutfileNamesList, pAcceptedScoresList):

        matrix_combination_key = key_outer[0][0] + '_' + key_outer[1][0]
        if matrix_combination_key not in aggregateFileH5Object:
            matrixCombinationGroup = aggregateFileH5Object.create_group(matrix_combination_key)
        else:
            matrixCombinationGroup = aggregateFileH5Object[matrix_combination_key]

        for key, data in zip(key_outer, data_outer):
            if len(data) == 0:
                continue
            else:
                counter += 1
            chromosome = None
            start_list = []
            end_list = []
            gene_name = None
            sum_of_interactions = None
            relative_distance_list = []

            raw_target_list = []

            for key_accepted in data:
                chromosome = data[key_accepted][0]
                start_list.append(data[key_accepted][1])
                end_list.append(data[key_accepted][2])
                gene_name = data[key_accepted][3]
                sum_of_interactions = data[key_accepted][4]
                relative_distance_list.append(data[key_accepted][5])
                raw_target_list.append(data[key_accepted][-1])

            if key[0] not in matrixCombinationGroup:
                matrixGroup = matrixCombinationGroup.create_group(key[0])
            else:
                matrixGroup = matrixCombinationGroup[key[0]]
            if key[1] not in matrixGroup:
                chromosomeObject = matrixGroup.create_group(key[1])
            else:
                chromosomeObject = matrixGroup[chromosome]

            if 'genes' not in matrixGroup:
                geneGroup = matrixGroup.create_group('genes')
            else:
                geneGroup = matrixGroup['genes']

            groupObject = chromosomeObject.create_group(key[-1])
            groupObject["chromosome"] = chromosome
            groupObject.create_dataset("start_list", data=start_list, compression="gzip", compression_opts=9)
            groupObject.create_dataset("end_list", data=end_list, compression="gzip", compression_opts=9)
            groupObject["gene_name"] = gene_name
            groupObject["sum_of_interactions"] = sum_of_interactions
            groupObject.create_dataset("relative_distance_list", data=relative_distance_list, compression="gzip", compression_opts=9)
            groupObject.create_dataset("raw_target_list", data=raw_target_list, compression="gzip", compression_opts=9)

            geneGroup[key[-1]] = chromosomeObject[key[-1]]

    aggregateFileH5Object.close()


def run_target_list_compilation(pInteractionFilesList, pTargetList, pArgs, pViewpointObj, pQueue=None, pOneTarget=False):
    outfile_names_list = []
    accepted_scores_list = []
    target_regions_intervaltree = None
    targetPosDict = {}
    try:
        if not h5py.is_hdf5(pArgs.targetFile):
            targetPosDict = {}
            outer_matrix = 'c_adj_norm'
            inner_matrix = 't_adj_norm'
            with open(pArgs.targetFile, 'r') as file:
                for line in file.readlines():
                    if line.startswith('#'):
                        continue
                    _line = line.strip().split('\t')
                    if len(line) == 0:
                        continue
                    try:
                        chrom, start, end = _line[:4]
                    except ValueError:
                        _line = line.strip().split()
                        chrom, start, end, gene = _line[:4]
                    if gene not in targetPosDict:
                        targetPosDict[gene] = {}
                        targetPosDict[gene]['chromosome'] = {}
                        targetPosDict[gene]['start_list'] = []
                        targetPosDict[gene]['end_list'] = []
                    targetPosDict[gene]['chromosome'] = chrom
                    targetPosDict[gene]['start_list'].append(start)
                    targetPosDict[gene]['end_list'].append(end)
        else:
            targetPosDict = None

        for i, interactionFile in enumerate(pInteractionFilesList):
            outfile_names_list_intern = []
            accepted_scores_list_intern = []
            for sample in interactionFile:
                interaction_data, interaction_file_data, _ = pViewpointObj.readInteractionFile(pArgs.interactionFile, sample)
                target_file = pTargetList[i]
                accepted_scores = filter_scores_target_list(interaction_file_data, pTargetList=target_file, pTargetPosDict=targetPosDict, pTargetFile=pArgs.targetFile)
                outfile_names_list_intern.append(sample)
                accepted_scores_list_intern.append(accepted_scores)
            outfile_names_list.append(outfile_names_list_intern)
            accepted_scores_list.append(accepted_scores_list_intern)

    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    if pQueue is None:
        return
    counter = 0
    for item in accepted_scores_list_intern:
        if len(item) == 0:
            counter += 1
    pQueue.put([outfile_names_list, accepted_scores_list])
    return


def call_multi_core(pInteractionFilesList, pTargetFileList, pFunctionName, pArgs, pViewpointObj):
    if len(pInteractionFilesList) < pArgs.threads:
        pArgs.threads = len(pInteractionFilesList)
    outfile_names_list = [None] * pArgs.threads
    accepted_scores_list = [None] * pArgs.threads

    interactionFilesPerThread = len(pInteractionFilesList) // pArgs.threads

    all_data_collected = False
    queue = [None] * pArgs.threads
    process = [None] * pArgs.threads
    thread_done = [False] * pArgs.threads
    one_target = True if len(pTargetFileList) == 1 else False
    fail_flag = False
    fail_message = ''
    for i in range(pArgs.threads):

        if i < pArgs.threads - 1:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
            if len(pTargetFileList) == 1:
                targetFileListThread = pTargetFileList[0]
            else:
                targetFileListThread = pTargetFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
        else:
            interactionFileListThread = pInteractionFilesList[i * interactionFilesPerThread:]
            if len(pTargetFileList) == 1:
                targetFileListThread = pTargetFileList[0]
            else:
                targetFileListThread = pTargetFileList[i * interactionFilesPerThread:]

        queue[i] = Queue()
        process[i] = Process(target=pFunctionName, kwargs=dict(
            pInteractionFilesList=interactionFileListThread,
            pTargetList=targetFileListThread,
            pArgs=pArgs,
            pViewpointObj=pViewpointObj,
            pQueue=queue[i],
            pOneTarget=one_target
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
                    outfile_names_list[i], accepted_scores_list[i] = background_data_thread
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
    outfile_names_list = [item for sublist in outfile_names_list for item in sublist]
    accepted_scores_list = [item for sublist in accepted_scores_list for item in sublist]

    return outfile_names_list, accepted_scores_list


def readTargetBed(targetList):
    present_genes = {}
    targetDict = {}        
    outer_matrix = 'c_adj_norm'
    inner_matrix = 't_adj_norm'
    present_genes[outer_matrix] = {}
    present_genes[outer_matrix][inner_matrix] = []
    with open(targetList, 'r') as file:
        for line in file.readlines():
            if line.startswith('#'):
                continue
            _line = line.strip().split('\t')
            if len(line) == 0:
                continue
            try:
                chrom, start, end = _line[:4]
            except ValueError:
                _line = line.strip().split()
                chrom, start, end, gene = _line[:4]
            if gene not in present_genes[outer_matrix][inner_matrix]:
                present_genes[outer_matrix][inner_matrix].append(gene)
            targetDict[gene] = [outer_matrix, inner_matrix, 'genes', gene]
    return present_genes, targetDict

def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()

    interactionList = []
    interactionDict = {}

    targetList = []
    present_genes = {}
    # read hdf file
    interactionFileHDF5Object = h5py.File(args.interactionFile, 'r')
    keys_interactionFile = list(interactionFileHDF5Object.keys())

    if h5py.is_hdf5(args.targetFile):
        targetDict, present_genes = viewpointObj.readTargetHDFFile(args.targetFile)

    else:
        present_genes, targetDict = readTargetBed(args.targetFile)

    for i, sample in enumerate(keys_interactionFile):
        for sample2 in keys_interactionFile[i + 1:]:

            matrix_obj1 = interactionFileHDF5Object[sample]
            matrix_obj2 = interactionFileHDF5Object[sample2]

            chromosomeList1 = sorted(list(matrix_obj1.keys()))
            chromosomeList2 = sorted(list(matrix_obj2.keys()))
            chromosomeList1.remove('genes')
            chromosomeList2.remove('genes')
            for chromosome1, chromosome2 in zip(chromosomeList1, chromosomeList2):
                geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                geneList2 = sorted(list(matrix_obj2[chromosome2].keys()))

                for gene1, gene2 in zip(geneList1, geneList2):
                    if gene1 in present_genes[sample][sample2]:
                        interactionDict[gene1] = [[sample, chromosome1, gene1], [sample2, chromosome2, gene2]]

    interactionFileHDF5Object.close()

    key_outer_matrix = present_genes.keys()

    for keys_outer in key_outer_matrix:
        keys_inner_matrix = present_genes[keys_outer].keys()
        for keys_inner in keys_inner_matrix:
            for gene in present_genes[keys_outer][keys_inner]:
                interactionList.append(interactionDict[gene])
                targetList.append(targetDict[gene])
    outfile_names_list, accepted_scores_list = call_multi_core(interactionList, targetList, run_target_list_compilation, args, viewpointObj)
    writeAggregateHDF(args.outFileName, outfile_names_list, accepted_scores_list, args)
