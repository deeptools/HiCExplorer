"""chicSignificantInteractions: a thin layer over the C++ executable chicSignificantInteractions.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
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
    """The argparse parser of chicSignificantInteractions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicSignificantInteractions')


def compute_interaction_file(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution, pQueue=None):
    target_outfile_names = []
    significant_data_list = []
    significant_key_list = []

    target_data_list = []
    target_key_list = []

    reference_points_list_target = []
    reference_points_list_significant = []

    try:
        for interactionFile in pInteractionFilesList:
            target_list = []
            sample_prefix = []

            for sample in interactionFile:
                data = pViewpointObj.readInteractionFile(pFilePath, sample)
                sample_prefix.append(sample[0])
                compute_new_p_values_bool = False
                if pArgs.xFoldBackground is not None:
                    accepted_scores, merged_lines_dict = merge_neighbors_x_fold(
                        pArgs.xFoldBackground, data, pViewpointObj, pResolution=pResolution)
                    compute_new_p_values_bool = True
                elif pArgs.loosePValue is not None:
                    accepted_scores, merged_lines_dict = merge_neighbors_loose_p_value(
                        pArgs.loosePValue, data, pViewpointObj, pResolution=pResolution, pTruncateZeroPvalues=pArgs.truncateZeroPvalues)
                    compute_new_p_values_bool = True

                # compute new p-values and filter by them
                if compute_new_p_values_bool:
                    accepted_scores, target_lines = compute_new_p_values(
                        accepted_scores, pBackground, pArgs.pValue, merged_lines_dict, pArgs.peakInteractionsThreshold)
                else:
                    accepted_scores, target_lines = filter_by_pvalue(
                        data[0], pArgs.pValue, data[1], pArgs.peakInteractionsThreshold)

                # filter by new p-value
                if len(accepted_scores) == 0:
                    with open('errorLog.txt', 'a+') as errorlog:
                        errorlog.write('Failed for: {}.\n'.format(
                            interactionFile))

                target_list.append(target_lines)
                significant_data_list.append(accepted_scores)
                significant_key_list.append(sample)
                reference_points_list_significant.append(data[2])

            target_list = [item for sublist in target_list for item in sublist]
            sample_prefix.append('::')
            sample_prefix.append(interactionFile[0][1])
            sample_prefix.append(interactionFile[0][2])
            target_data_list.append(target_list)
            target_key_list.append(sample_prefix)
            reference_points_list_target.append(data[2])

        if pQueue is None:
            return target_outfile_names
    except Exception as exp:
        if pQueue is None:
            return 'Fail: ' + str(exp) + traceback.format_exc()
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put([significant_data_list, significant_key_list, target_data_list, target_key_list, reference_points_list_target, reference_points_list_significant])
    return


def filter_by_pvalue(pData, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
    accepted = {}
    accepted_lines = []
    if isinstance(pPValue, float):
        for key in pData:
            if pData[key][-3] <= pPValue:
                if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                    accepted[key] = pMergedLinesDict[key]
                    target_content = pMergedLinesDict[key][:3]
                    accepted_lines.append(target_content)
    elif isinstance(pPValue, dict):
        for key in pData:

            if pData[key][-3] <= pPValue[key]:
                if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                    accepted[key] = pMergedLinesDict[key]
                    target_content = pMergedLinesDict[key][:3]
                    accepted_lines.append(target_content)

    return accepted, accepted_lines


def compute_new_p_values(pData, pBackgroundModel, pPValue, pMergedLinesDict, pPeakInteractionsThreshold):
    accepted = {}
    accepted_lines = []
    if isinstance(pPValue, float):
        for key in pData:
            if key in pBackgroundModel:

                pData[key][-3] = 1 - cnb.cdf(float(pData[key][-1]), float(pBackgroundModel[key][0]), float(pBackgroundModel[key][1]))
                if pData[key][-3] <= pPValue:
                    if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                        accepted[key] = pData[key]
                        target_content = pMergedLinesDict[key][0][:3]

                        target_content[2] = pMergedLinesDict[key][-1][2]
                        accepted_lines.append(target_content)
            else:
                log.debug('key not in background {}'.format(key))
    elif isinstance(pPValue, dict):
        for key in pData:
            if key in pBackgroundModel:
                pData[key][-3] = 1 - cnb.cdf(float(pData[key][-1]), float(pBackgroundModel[key][0]), float(pBackgroundModel[key][1]))

                if pData[key][-3] <= pPValue[key]:
                    if float(pData[key][-1]) >= pPeakInteractionsThreshold:
                        accepted[key] = pData[key]
                        target_content = pMergedLinesDict[key][0][:3]
                        target_content[2] = pMergedLinesDict[key][-1][2]
                        accepted_lines.append(target_content)
            else:
                log.debug('key not in background {}'.format(key))
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


def call_multi_core(pInteractionFilesList, pArgs, pViewpointObj, pBackground, pFilePath, pResolution):
    significant_data_list = [None] * pArgs.threads
    significant_key_list = [None] * pArgs.threads

    target_data_list = [None] * pArgs.threads
    target_key_list = [None] * pArgs.threads
    reference_points_list_target = [None] * pArgs.threads
    reference_points_list_significant = [None] * pArgs.threads

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
                    significant_data_list[i], significant_key_list[i], target_data_list[i], target_key_list[i], reference_points_list_target[i], reference_points_list_significant[i] = background_data_thread
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

    significant_data_list = [item for sublist in significant_data_list for item in sublist]
    significant_key_list = np.array([item for sublist in significant_key_list for item in sublist])
    reference_points_list_target = np.array([item for sublist in reference_points_list_target for item in sublist])
    reference_points_list_significant = np.array([item for sublist in reference_points_list_significant for item in sublist])

    significant_key_list, indices = np.unique(significant_key_list, axis=0, return_index=True)

    significant_data_list_new = []
    for x in indices:
        significant_data_list_new.append(significant_data_list[x])
    significant_data_list = significant_data_list_new

    target_data_list = [item for sublist in target_data_list for item in sublist]
    target_key_list = [item for sublist in target_key_list for item in sublist]

    return significant_data_list, significant_key_list, target_data_list, target_key_list, reference_points_list_target, reference_points_list_significant


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


def writeSignificantHDF(pOutFileName, pSignificantDataList, pSignificantKeyList, pViewpointObj, pReferencePointsList, pArgs):

    significantFileH5Object = h5py.File(pOutFileName, 'w')
    significantFileH5Object.attrs['type'] = "significant"
    significantFileH5Object.attrs['version'] = __version__

    significantFileH5Object.attrs['pvalue'] = pArgs.pValue

    if pArgs.xFoldBackground is not None:
        significantFileH5Object.attrs['mode_preselection'] = "xfold"
        significantFileH5Object.attrs['mode_preselection_value'] = pArgs.xFoldBackground

    elif pArgs.loosePValue is not None:
        significantFileH5Object.attrs['mode_preselection'] = "loosePValue"
        significantFileH5Object.attrs['mode_preselection_value'] = pArgs.loosePValue
    else:
        significantFileH5Object.attrs['mode_preselection'] = "None"
        significantFileH5Object.attrs['mode_preselection_value'] = "None"

    significantFileH5Object.attrs['range'] = pArgs.range
    significantFileH5Object.attrs['combinationMode'] = pArgs.combinationMode

    significantFileH5Object.attrs['truncateZeroPvalues'] = pArgs.truncateZeroPvalues
    significantFileH5Object.attrs['fixateRange'] = pArgs.fixateRange
    significantFileH5Object.attrs['peakInteractionsThreshold'] = pArgs.peakInteractionsThreshold

    keys_seen = {}

    for i, (key, data) in enumerate(zip(pSignificantKeyList, pSignificantDataList)):
        if len(data) == 0:
            continue
        chromosome = None
        start_list = []
        end_list = []
        gene_name = None
        sum_of_interactions = None
        relative_distance_list = []
        relative_interactions_list = []
        pvalue_list = []
        xfold_list = []
        raw_target_list = []

        for datum in data.values():
            chromosome = datum[0]
            start_list.append(datum[1])
            end_list.append(datum[2])
            gene_name = datum[3]
            sum_of_interactions = datum[4]
            relative_distance_list.append(datum[5])
            relative_interactions_list.append(datum[6])
            pvalue_list.append(datum[7])
            xfold_list.append(datum[8])
            raw_target_list.append(datum[9])

        if key[0] not in significantFileH5Object:
            matrixGroup = significantFileH5Object.create_group(key[0])
            keys_seen[key[0]] = set()
        else:
            matrixGroup = significantFileH5Object[key[0]]

        if chromosome not in matrixGroup:
            chromosomeObject = matrixGroup.create_group(chromosome)
        else:
            chromosomeObject = matrixGroup[chromosome]

        if 'genes' not in matrixGroup:
            geneGroup = matrixGroup.create_group('genes')
        else:
            geneGroup = matrixGroup['genes']

        success = False
        counter = 0
        while not success:

            if counter != 0:
                gene_name_key = key[2] + '_' + str(counter)
            else:
                gene_name_key = key[2]
            if gene_name_key in keys_seen[key[0]]:
                success = False
            else:
                keys_seen[key[0]].add(gene_name_key)
                success = True
            counter += 1

        group_name = pViewpointObj.writeInteractionFileHDF5(
            chromosomeObject, gene_name_key, [chromosome, start_list, end_list, gene_name, sum_of_interactions, relative_distance_list,
                                              relative_interactions_list, pvalue_list, xfold_list, raw_target_list], pReferencePointsList[i])

        try:
            geneGroup[group_name] = chromosomeObject[group_name]
        except Exception as exp:
            log.debug('exception {}'.format(str(exp)))
            log.debug('Gene group given: {}'.format(key[2]))
            log.debug('Gene group return: {}'.format(group_name))

    significantFileH5Object.close()


def writeTargetHDF(pOutFileName, pTargetDataList, pTargetKeyList, pViewpointObj, pResolution, pReferencePointsList, pArgs):

    targetFileH5Object = h5py.File(pOutFileName, 'w')
    targetFileH5Object.attrs['type'] = "target"
    targetFileH5Object.attrs['version'] = __version__

    targetFileH5Object.attrs['pvalue'] = pArgs.pValue

    if pArgs.xFoldBackground is not None:
        targetFileH5Object.attrs['mode_preselection'] = "xfold"
        targetFileH5Object.attrs['mode_preselection_value'] = pArgs.xFoldBackground

    elif pArgs.loosePValue is not None:
        targetFileH5Object.attrs['mode_preselection'] = "loosePValue"
        targetFileH5Object.attrs['mode_preselection_value'] = pArgs.loosePValue
    else:
        targetFileH5Object.attrs['mode_preselection'] = "None"
        targetFileH5Object.attrs['mode_preselection_calue'] = "None"

    targetFileH5Object.attrs['range'] = pArgs.range
    targetFileH5Object.attrs['combinationMode'] = pArgs.combinationMode

    targetFileH5Object.attrs['truncateZeroPvalues'] = pArgs.truncateZeroPvalues
    targetFileH5Object.attrs['fixateRange'] = pArgs.fixateRange
    targetFileH5Object.attrs['peakInteractionsThreshold'] = pArgs.peakInteractionsThreshold
    keys_seen = {}
    for i, (key, data) in enumerate(zip(pTargetKeyList, pTargetDataList)):
        if len(data) == 0:
            continue

        chromosome = None
        start_list = []
        end_list = []
        a = pybedtools.BedTool(data)
        data_sorted_merged = a.sort().merge(d=pResolution)
        for datum in data_sorted_merged:
            chromosome = datum[0]
            start_list.append(datum[1])
            end_list.append(datum[2])

        if key[0] not in targetFileH5Object:
            matrixGroup = targetFileH5Object.create_group(key[0])
            keys_seen[key[0]] = set()
        else:
            matrixGroup = targetFileH5Object[key[0]]

        for matrix_name in key[1:]:
            if matrix_name == '::':
                break
            if matrix_name not in matrixGroup:
                matrixGroup = matrixGroup.create_group(matrix_name)
                keys_seen[matrix_name] = set()
            else:
                matrixGroup = matrixGroup[matrix_name]

        if 'genes' not in matrixGroup:
            geneGroup = matrixGroup.create_group('genes')
        else:
            geneGroup = matrixGroup['genes']
        if chromosome not in matrixGroup:
            chromosomeObject = matrixGroup.create_group(chromosome)
        else:
            chromosomeObject = matrixGroup[chromosome]

        success = False
        counter = 0
        while not success:

            if counter != 0:
                gene_name_key = key[-1] + '_' + str(counter)
            else:
                gene_name_key = key[-1]
            if gene_name_key in keys_seen[key[0]]:
                success = False
            else:
                keys_seen[key[0]].add(gene_name_key)
                success = True
            counter += 1

        groupObject, groupName = pViewpointObj.createUniqueHDFGroup(chromosomeObject, gene_name_key)

        groupObject["chromosome"] = chromosome
        groupObject.create_dataset("start_list", data=start_list)
        groupObject.create_dataset("end_list", data=end_list)
        groupObject.create_dataset("reference_point_start", data=int(pReferencePointsList[i][0]))
        groupObject.create_dataset("reference_point_end", data=int(pReferencePointsList[i][1]))
        try:
            geneGroup[groupName] = chromosomeObject[groupName]
        except Exception as exp:
            log.debug('exception {}'.format(str(exp)))
            log.debug('Gene group: {}'.format(key[-1]))
            log.debug('Adjusted name {}'.format(groupName))
    targetFileH5Object.close()


def main(args=None):
    """Run the C++ implementation of chicSignificantInteractions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicSignificantInteractions')(args)
