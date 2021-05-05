import argparse
import sys
import os
import errno
from multiprocessing import Process, Queue
import time
import math
import logging
log = logging.getLogger(__name__)

import h5py
import numpy as np

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from .lib import Viewpoint
from hicexplorer._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
Computes per input matrix all viewpoints which are defined in the reference points file.

""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='Path to the Hi-C matrices which store the captured Hi-C data per sample.',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be considered in the analysis. Please remember to use the same fixate range setting as '
                                'for the background model computation and that distances of the range larger than the fixate range use the background model of those.'
                                'Format is --region upstream downstream',
                                required=True,
                                type=int,
                                nargs=2)

    parserRequired.add_argument('--referencePoints', '-rp', help='Reference point file. Needs to be in the format: \'chr 100\' for a '
                                'single reference point or \'chr 100 200\' for a reference region and with a single reference point per line',
                                required=True)
    parserRequired.add_argument('--backgroundModelFile', '-bmf',
                                help='path to the background file computed by chicViewpointBackgroundModel',
                                required=True)
    parserRequired.add_argument('--outFileName', '-o',
                           help='This hdf5 file contains all created viewpoint files.',
                           required=False,
                           default='chic_files.hdf5')
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins via a sliding window approach to smooth the values and be less sensitive for outliers'
                           ' (Default: %(default)s).',
                           type=int,
                           default=5)

    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate range of background model starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin'
                           ' (Default: %(default)s).',
                           required=False,
                           default=500000,
                           type=int
                           )

    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def adjustViewpointData(pViewpointObj, pData, pBackground, pReferencePoint, pRegionStart, pRegionEnd):
    data_viewpoint = {}
    data_background = {}
    view_point_start, _ = pViewpointObj.getReferencePointAsMatrixIndices(
        pReferencePoint)
    view_point_range_start, view_point_range_end = \
        pViewpointObj.getViewpointRangeAsMatrixIndices(
            pReferencePoint[0], pRegionStart, pRegionEnd)

    for i, data in zip(range(view_point_range_start, view_point_range_end, 1), pData):
        relative_position = i - view_point_start
        data_viewpoint[relative_position] = data
    for i, data in zip(range(view_point_range_start, view_point_range_end, 1), pBackground):
        relative_position = i - view_point_start
        data_background[relative_position] = data

    for i in data_background:
        if i in data_viewpoint:
            continue
        else:
            data_viewpoint[i] = 0

    data = np.fromiter(data_viewpoint.values(), dtype=np.float32)
    background = list(data_background.values())

    return data, background


def compute_x_fold(pDataList, pBackgroundList):
    return pDataList / pBackgroundList


def compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel, pBackgroundModelRelativeInteractions):
    file_list = []
    interaction_data_list = []

    try:
        for i, referencePoint in enumerate(pReferencePoints):
            # range of viewpoint with reference point in the middle in genomic units
            # get fixateRange for relative interaction computation denominator
            region_start_fixed, region_end_fixed, range_fixed = pViewpointObj.calculateViewpointRange(
                referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

            intermediate_viewpoint, _ = pViewpointObj.computeViewpoint(
                referencePoint, referencePoint[0], region_start_fixed, region_end_fixed)
            denominator_relative_interactions = np.sum(intermediate_viewpoint)

            # viewpoint data uses full range
            region_start, region_end, _range = pViewpointObj.calculateViewpointRange(
                referencePoint, pArgs.range)

            data_list, index_reference_point = pViewpointObj.computeViewpoint(
                referencePoint, referencePoint[0], region_start, region_end)

            # background uses fixed range, handles fixate range implicitly by same range used in background computation

            background_relative_interaction = pViewpointObj.interactionBackgroundData(
                pBackgroundModelRelativeInteractions, _range).flatten()
            data_list_relative = data_list
            if len(data_list) != len(background_relative_interaction):
                data_list, background_relative_interaction = adjustViewpointData(
                    pViewpointObj, data_list_relative, background_relative_interaction, referencePoint, region_start, region_end)

            if pArgs.averageContactBin > 0 and len(data_list) >= pArgs.averageContactBin:
                data_list = pViewpointObj.smoothInteractionValues(
                    data_list, pArgs.averageContactBin)

            data_list_raw = np.copy(data_list)

            data_list = pViewpointObj.computeRelativeValues(
                data_list, denominator_relative_interactions)

            x_fold_list = compute_x_fold(
                data_list, background_relative_interaction)
            p_value_list = pViewpointObj.pvalues(
                pBackgroundModel, data_list_raw, index_reference_point)

            # add values if range is larger than fixate range

            region_start_range, region_end_range, _ = pViewpointObj.calculateViewpointRange(
                referencePoint, (pArgs.range[0], pArgs.range[1]))

            interaction_data = pViewpointObj.createInteractionFileDataHDF5(referencePoint, referencePoint[0],
                                                                       region_start_range, region_end_range, data_list, data_list_raw,
                                                                       pGeneList[i], denominator_relative_interactions, p_value_list, x_fold_list)

            referencePointString = '_'.join(str(j) for j in referencePoint)

            region_start_in_units = utilities.in_units(region_start)
            region_end_in_units = utilities.in_units(region_end)
            denominator_relative_interactions_str = 'Sum of interactions in fixate range: '
            denominator_relative_interactions_str += str(
                denominator_relative_interactions)
            header_information = '# Interaction file, created with HiCExplorer\'s chicViewpoint version ' + \
                __version__ + '\n# '
            header_information += '\t'.join([pMatrix, referencePointString, str(region_start_in_units), str(
                region_end_in_units), pGeneList[i], denominator_relative_interactions_str])
            header_information += '\n# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative Interactions\tp-value\tx-fold\tRaw\n#'
            matrix_name = '.'.join(pMatrix.split('/')[-1].split('.')[:-1])
            matrix_name = '_'.join(
                [matrix_name, referencePointString, pGeneList[i]])
            referencePointGenename = '_'.join([referencePointString, pGeneList[i]])
            referencePointGenename = str(referencePointGenename)
            interaction_data_list.append([matrix_name, interaction_data, header_information, referencePointGenename])
    except Exception as exp:
        log.debug('Error! {}'.format(str(exp)))
        pQueue.put('Fail: ' + str(exp))
        return
    pQueue.put(interaction_data_list)
    return


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()

    referencePoints, gene_list = viewpointObj.readReferencePointFile(
        args.referencePoints)
    
    log.debug('referencePoints {}'.format(referencePoints[:5]))
    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    file_list = []
    background_model = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range, args.fixateRange)
    background_model_mean_values = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range, args.fixateRange, pMean=True)

    # create hdf5 output file
    interactionFileH5Object = h5py.File(args.outFileName, 'w')
    interactionFileH5Object.attrs['type'] = 'interactions'
    interactionFileH5Object.attrs['range'] = args.range
    interactionFileH5Object.attrs['averageContactBin'] = args.averageContactBin
    interactionFileH5Object.attrs['fixateRange'] = args.fixateRange




    fail_flag = False
    fail_message = ''
    matrix_collection = {}
    resolution = 0
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma
        file_list_sample = [None] * args.threads
        interaction_data_list_sample = [None] * args.threads

        all_data_collected = False

        if resolution == 0:
            resolution = hic_ma.getBinSize()
            log.debug('resolution {}'.format(resolution))
            interactionFileH5Object.attrs['resolution'] = resolution

        for i in range(args.threads):

            if i < args.threads - 1:
                referencePointsThread = referencePoints[i * referencePointsPerThread:(i + 1) * referencePointsPerThread]
                geneListThread = gene_list[i * referencePointsPerThread:(i + 1) * referencePointsPerThread]
            else:
                referencePointsThread = referencePoints[i * referencePointsPerThread:]
                geneListThread = gene_list[i * referencePointsPerThread:]

            if len(referencePointsThread) == 0:
                process[i] = None
                queue[i] = None
                file_list_sample[i] = []
                continue
            queue[i] = Queue()
            process[i] = Process(target=compute_viewpoint, kwargs=dict(
                pViewpointObj=viewpointObj,
                pArgs=args,
                pQueue=queue[i],
                pReferencePoints=referencePointsThread,
                pGeneList=geneListThread,
                pMatrix=matrix,
                pBackgroundModel=background_model,
                pBackgroundModelRelativeInteractions=background_model_mean_values
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    file_list_ = queue[i].get()
                    if 'Fail:' in file_list_:
                        fail_flag = True
                        fail_message = file_list_[6:]
                    interaction_data_list_sample[i] = file_list_
                    process[i].join()
                    process[i].terminate()
                    process[i] = None

            all_data_collected = True

            for i in range(args.threads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)

        if fail_flag:
            log.error(fail_message)
            exit(1)

        interaction_data_list = [item for sublist in interaction_data_list_sample for item in sublist]
        matrix_collection[matrix] = interaction_data_list

    for matrix in matrix_collection:
        matrixGroup = interactionFileH5Object.create_group(matrix.split('.')[0])
        geneGroup = matrixGroup.create_group('genes')

        for i, interaction_data in enumerate(matrix_collection[matrix]):
            if interaction_data[1][0] not in matrixGroup:
                chromosomeObject = matrixGroup.create_group(interaction_data[1][0])

            group_name = viewpointObj.writeInteractionFileHDF5(
                    chromosomeObject, interaction_data[1][3], interaction_data[1], referencePoints[i][1:])

            try:
                geneGroup[group_name] = chromosomeObject[group_name]
            except Exception as e:
                log.debug(str(e))
                log.debug('group_name {}'.format(group_name))
                log.debug('gene name {}'.format(interaction_data[1][3]))
