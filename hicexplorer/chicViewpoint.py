import argparse
import sys
import numpy as np
import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from .lib import Viewpoint
from hicexplorer._version import __version__
from scipy.stats import zscore
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from multiprocessing import Process, Queue
import time
import math
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Plots the number of interactions around a given reference point in a region.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='path of the Hi-C matrices to plot',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream',
                                required=True,
                                type=int,
                                nargs=2)

    parserRequired.add_argument('--referencePoints', '-rp', help='Reference point file. Needs to be in the format: \'chr 100\' for a '
                                'single reference point or \'chr 100 200\' for a reference region and per line one reference point',
                                required=True)
    parserRequired.add_argument('--backgroundModelFile', '-bmf',
                                help='path to the background file which is necessary to compute the rbz-score',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins, written to last column.',
                           type=int,
                           default=0)
    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate range of backgroundmodel starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin.',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel):

    for i, referencePoint in enumerate(pReferencePoints):
        # range of viewpoint with reference point in the middle in genomic units
        region_start, region_end, _range = pViewpointObj.calculateViewpointRange(referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

        data_list = pViewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)
        if pArgs.averageContactBin > 0:
            data_list = pViewpointObj.smoothInteractionValues(data_list, pArgs.averageContactBin)
        data_list_raw = np.copy(data_list)

        data_list = pViewpointObj.computeRelativeValues(data_list)

        # if args.backgroundModelFile:
        # _background_model = pViewpointObj.readBackgroundDataFile(args.backgroundModelFile)
        _backgroundModelData, _backgroundModelSEM = pViewpointObj.interactionBackgroundData(pBackgroundModel, _range)
        rbz_score_data = pViewpointObj.rbz_score(data_list, _backgroundModelData, _backgroundModelSEM)

        # add values if range is larger than fixate range

        max_chromosome = pViewpointObj.hicMatrix.getBinPos(pViewpointObj.hicMatrix.getChrBinRange(referencePoint[0])[1] - 1)[2]
        min_chromosome = pViewpointObj.hicMatrix.getBinPos(pViewpointObj.hicMatrix.getChrBinRange(referencePoint[0])[0] - 1)[1]
        if region_end >= max_chromosome - 1:
            downstream_range = pArgs.fixateRange - ((int(referencePoint[2]) + pArgs.fixateRange) - max_chromosome)
        else:
            downstream_range = pArgs.fixateRange

        # if region_start < min_chromosome:
        #     upstream_range = pArgs.fixateRange - ((int(referencePoint[1]) - pArgs.fixateRange) - min_chromosome)
        # else:
        upstream_range = pArgs.fixateRange

        log.debug('min_chromosome {}'.format(min_chromosome))
        log.debug('max_chromosome {}'.format(max_chromosome))
        log.debug('upstream_range {}'.format(upstream_range))
        log.debug('downstream_range {}'.format(downstream_range))

        difference_upstream = -pArgs.range[0] - (-upstream_range)
        difference_downstream = pArgs.range[1] - downstream_range
        difference_upstream //= pViewpointObj.hicMatrix.getBinSize()
        difference_downstream //= pViewpointObj.hicMatrix.getBinSize()

        log.debug('len rbz_score_data {}'.format(len(rbz_score_data)))
        log.debug('len data_list {}'.format(len(data_list)))
        log.debug('len data_list_raw {}'.format(len(data_list_raw)))

        if difference_upstream < 0:
            # extend with first position
            log.debug('extending upstream')

            rbz_score_data_extend = np.empty(np.absolute(difference_upstream))
            data_list_extend = np.empty(np.absolute(difference_upstream))
            data_list_raw_extend = np.empty(np.absolute(difference_upstream))

            rbz_score_data_extend[:] = rbz_score_data[0]
            data_list_extend[:] = data_list[0]
            data_list_raw_extend[:] = data_list_raw[0]


            rbz_score_data = np.concatenate(rbz_score_data_extend, rbz_score_data)
            data_list = np.concatenate(data_list_extend, data_list)
            data_list_raw = np.concatenate(data_list_raw_extend, data_list_raw)

        elif difference_upstream > 0:
            # clip data
            log.debug('clipping upstream')
            log.debug('difference_upstream {}'.format(difference_upstream))

            # log.debug('rbz_score_data {}'.format(rbz_score_data))

            rbz_score_data = rbz_score_data[difference_upstream:]
            # log.debug('rbz_score_data {}'.format(rbz_score_data))

            data_list = data_list[difference_upstream:]
            data_list_raw = data_list_raw[difference_upstream:]
        if difference_downstream < 0:
            # clip
            log.debug('clipping downstream')
            log.debug('difference_downstream {}'.format(difference_downstream))

            rbz_score_data = rbz_score_data[:difference_downstream]
            data_list = data_list[:difference_downstream]
            data_list_raw = data_list_raw[:difference_downstream]
        elif difference_downstream > 0:
            # extend
            log.debug('extending downstream')

            rbz_score_data_extend = np.empty(difference_downstream)
            data_list_extend = np.empty(difference_downstream)
            data_list_raw_extend = np.empty(difference_downstream)

            rbz_score_data_extend[:] = rbz_score_data[-1]
            data_list_extend[:] = data_list[-1]
            data_list_raw_extend[:] = data_list[-1]

            rbz_score_data = np.concatenate(rbz_score_data, rbz_score_data_extend)
            data_list = np.concatenate(data_list, data_list_extend)
            data_list_raw = np.concatenate(data_list_raw, data_list_raw_extend)

        log.debug('len rbz_score_data {}'.format(len(rbz_score_data)))
        log.debug('len data_list {}'.format(len(data_list)))
        log.debug('len data_list_raw {}'.format(len(data_list_raw)))


        region_start_range, region_end_range, _ = pViewpointObj.calculateViewpointRange(referencePoint, (pArgs.range[0], pArgs.range[1]))


        log.debug('region_start_range {}, region_start_range {}'.format(region_start_range, region_end_range))
        log.debug('diff range {}'.format((region_start_range - region_end_range) // pViewpointObj.hicMatrix.getBinSize()))
        
        # if region_start_range > region_start:

        # if region_end_range < region_end:

        interaction_data = pViewpointObj.createInteractionFileData(referencePoint, referencePoint[0],
                                                                    region_start_range, region_end_range, data_list, data_list_raw,
                                                                    pGeneList[i])

        referencePointString = '_'.join(str(j) for j in referencePoint)

        region_start_in_units = utilities.in_units(region_start)
        region_end_in_units = utilities.in_units(region_end)

        header_information = '\t'.join([pMatrix, referencePointString, str(region_start_in_units), str(region_end_in_units), pGeneList[i]])
        header_information += '\n# ChrViewpoint\tStart\tEnd\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\trbz-score\tRaw\n#'
        matrix_name = '.'.join(pMatrix.split('.')[:-1])
        matrix_name = '_'.join([matrix_name, referencePointString, pGeneList[i]])
        pViewpointObj.writeInteractionFile(matrix_name, interaction_data, header_information, rbz_score_data)

    pQueue.put(['Done'])
    return

def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()

    referencePoints, gene_list = viewpointObj.readReferencePointFile(args.referencePoints)
    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    background_model = viewpointObj.readBackgroundDataFile(args.backgroundModelFile, args.range)
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma

        all_data_collected = False

        # log.debug('len(referencePoints) {}'.format(referencePoints))

        # compute_viewpoint(
        #         pViewpointObj = viewpointObj,
        #         pArgs = args,
        #         pQueue = None,
        #         pReferencePoints = referencePoints,
        #         pGeneList = gene_list,
        #         pMatrix = matrix,
        #         pBackgroundModel = background_model)

        for i in range(args.threads):
            
            if i < args.threads - 1:
                referencePointsThread = referencePoints[i*referencePointsPerThread:(i+1)*referencePointsPerThread]
                geneListThread = gene_list[i*referencePointsPerThread:(i+1)*referencePointsPerThread]
            else:
                referencePointsThread = referencePoints[i*referencePointsPerThread:]
                geneListThread = gene_list[i*referencePointsPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=compute_viewpoint, kwargs=dict(
                pViewpointObj = viewpointObj,
                pArgs = args,
                pQueue =queue[i],
                pReferencePoints = referencePointsThread,
                pGeneList = geneListThread,
                pMatrix = matrix,
                pBackgroundModel = background_model
                )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    _ = queue[i].get()
                    process[i].join()
                    process[i].terminate()
                    process[i] = None
                #     log.debug('Thread {} DONE'.format(i))
                # log.debug('Thread {} WAIT'.format(i))
                
            all_data_collected = True
            
            for i in range(args.threads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)


            # log.debug('referencePoint {}'.format(referencePoint))
            # region_start, region_end, _range = viewpointObj.calculateViewpointRange(referencePoint, args.range)

            # data_list = viewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)
            # if args.averageContactBin > 0:
            #     data_list = viewpointObj.smoothInteractionValues(data_list, args.averageContactBin)
            
            
            # bin_start_viewpoint, bin_end_viewpoint = viewpointObj.hicMatrix.getRegionBinRange(referencePoint[0], region_start, region_end)
            # # log.debug('region_start {}'.format(foo))
            # # log.debug('region_end {}'.format(region_end))

            # # log.debug('len(data_list) {}'.format(len(data_list)))
            # data_list = data_list[bin_start_viewpoint:bin_end_viewpoint]
            # # log.debug('len(data_list) {}'.format(len(data_list)))

            # data_list_raw = np.copy(data_list)

            # data_list = viewpointObj.computeRelativeValues(data_list)

            # if args.backgroundModelFile:
            #     _background_model = viewpointObj.readBackgroundDataFile(args.backgroundModelFile)
            #     _backgroundModelData, _backgroundModelSEM = viewpointObj.interactionBackgroundData(_background_model, _range)
            #     rbz_score_data = viewpointObj.rbz_score(data_list, _backgroundModelData, _backgroundModelSEM)

            # interaction_data = viewpointObj.createInteractionFileData(referencePoint, referencePoint[0],
            #                                                           region_start, region_end, data_list, data_list_raw,
            #                                                           gene_list[i])

            # referencePointString = '_'.join(str(j) for j in referencePoint)

            # region_start_in_units = utilities.in_units(region_start)
            # region_end_in_units = utilities.in_units(region_end)

            # header_information = '\t'.join([matrix, referencePointString, str(region_start_in_units), str(region_end_in_units), gene_list[i]])
            # header_information += '\n# ChrViewpoint\tStart\tEnd\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\trbz-score\tRaw\n#'
            # matrix_name = '.'.join(matrix.split('.')[:-1])
            # matrix_name = '_'.join([matrix_name, referencePointString, gene_list[i]])
            # viewpointObj.writeInteractionFile(matrix_name, interaction_data, header_information, rbz_score_data)
