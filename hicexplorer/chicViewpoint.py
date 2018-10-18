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

    # fixateRange = 
    for i, referencePoint in enumerate(pReferencePoints):
        # range of viewpoint with reference point in the middle in genomic units
        region_start, region_end, _range = pViewpointObj.calculateViewpointRange(referencePoint, pArgs.range)
        
        data_list = pViewpointObj.computeViewpoint(referencePoint, referencePoint[0], region_start, region_end)
        if pArgs.averageContactBin > 0:
            data_list = pViewpointObj.smoothInteractionValues(data_list, pArgs.averageContactBin)


        # these are absolute bin values for the full matrix
        # need to be adjusted for the chromosome
        reference_point_start, reference_point_end = pViewpointObj.getReferencePointAsMatrixIndices(referencePoint)
        # log.debug('reference_point_start {}, reference_point_end {}'.format(reference_point_start, reference_point_end ))
        # index values in bin units
        bin_start_viewpoint, bin_end_viewpoint = pViewpointObj.hicMatrix.getRegionBinRange(referencePoint[0], region_start, region_end)
        # log.debug('x {}, z{}'.format(reference_point_start-bin_start_viewpoint, bin_end_viewpoint-reference_point_end) )

        # start and end index of chromosome in bin units
        start_chromosome, end_chromosome = pViewpointObj.hicMatrix.getChrBinRange(referencePoint[0])
        bin_start_viewpoint = bin_start_viewpoint - start_chromosome
        bin_end_viewpoint = bin_end_viewpoint - start_chromosome 

        bin_end_viewpoint = bin_end_viewpoint - (reference_point_end - reference_point_start) + 1

        
        # if bin_start_viewpoint == bin_end_viewpoint:
        #     bin_end_viewpoint += 1
        # log.debug('len(data_list) {}'.format(len(data_list)))

        # log.debug('bin_start_viewpoint {} bin_end_viewpoint {}'.format(bin_start_viewpoint, bin_end_viewpoint))
        data_list_raw = np.copy(data_list)
        data_list_raw = data_list_raw[bin_start_viewpoint:bin_end_viewpoint]
        # log.debug('bin_start_viewpoint{}, bin_end_viewpoint {}'.format(bin_start_viewpoint, bin_end_viewpoint))
        
        data_list = pViewpointObj.computeRelativeValues(data_list)
        # log.debug('data_list before pruning {}'.format(len(data_list)))
        len_data_list = len(data_list)
        data_list = data_list[bin_start_viewpoint:bin_end_viewpoint]
        # log.debug('len(data_list) {}'.format(len(data_list)))
        # data_list_raw = np.copy(data_list)

        # data_list = pViewpointObj.computeRelativeValues(data_list)

        # if pArgs.backgroundModelFile:
        # log.debug('data_list {}'.format(len(data_list)))
        # log.debug('data_list_raw {}'.format(len(data_list_raw)))
        # log.debug('background: {}'.format(pBackgroundModel))
        # log.debug('len background: {}'.format(len(pBackgroundModel)))

        _backgroundModelData, _backgroundModelSEM = pViewpointObj.interactionBackgroundData(pBackgroundModel, _range)


        if len(data_list) != len(_backgroundModelData):
            # pass
            log.debug('referencePoint {}'.format(referencePoint))
            log.debug('bin_start_viewpoint {}  :bin_end_viewpoint {}'.format(bin_start_viewpoint, bin_end_viewpoint))
            log.debug('len(data_list) {}'.format(len(data_list)))
            log.debug('len(len(_backgroundModelData)) {}'.format(len(_backgroundModelData)))
            log.debug('len_data_list {}'.format(len_data_list))
            log.debug('start_chromosome {}, end_chromosome {}'.format(start_chromosome, end_chromosome))
            log.debug('view_point_start {}, view_point_end {}\n' .format(reference_point_start, reference_point_end))
            # log.debug('datalist[375:425] {}'.format(data_list[395:425]))
            # log.debug('_backgroundModelData[375:425] {}'.format(_backgroundModelData[395:425]))


        # if len(data_list) == len(_backgroundModelData):
        # log.debug('_backgroundModelData {}'.format(_backgroundModelData[:10]))
        # log.debug('_backgroundModelSEM {}'.format(_backgroundModelSEM[:10]))

        # set values which are in a distance larger than fixatedRange to value of index of this range.
        region_start_fixated, region_end_fixated, _ = pViewpointObj.calculateViewpointRange(referencePoint, (pArgs.fixateRange, pArgs.fixateRange))
        
        bin_start_viewpoint_fixated, bin_end_viewpoint_fixated = pViewpointObj.hicMatrix.getRegionBinRange(referencePoint[0], region_start_fixated, region_end_fixated)
        # log.debug('region_start_fixated {}, region_end_fixated {}'.format(region_start_fixated, region_end_fixated))
        # log.debug('bin_start_viewpoint_fixated {}, bin_end_viewpoint_fixated {}'.format(bin_start_viewpoint_fixated, bin_end_viewpoint_fixated))
        # log.debug('start_chromosome {}'.format(start_chromosome))
        bin_start_viewpoint_fixated = bin_start_viewpoint_fixated - start_chromosome
        bin_end_viewpoint_fixated = bin_end_viewpoint_fixated - start_chromosome

        # log.debug('>>>bin_start_viewpoint_fixated {}, bin_end_viewpoint_fixated {}'.format(bin_start_viewpoint_fixated, bin_end_viewpoint_fixated))

        bin_end_viewpoint_fixated = bin_end_viewpoint_fixated - (reference_point_end - reference_point_start)
        # log.debug('bin_end_viewpoint_fixated {}, bin_end_viewpoint_fixated {}'.format(region_start_fixated, region_end_fixated))
        # log.debug('bin_start_viewpoint {}, bin_end_viewpoint {}'.format(bin_start_viewpoint, bin_end_viewpoint))
        bin_start_viewpoint_fixated =  bin_start_viewpoint_fixated - bin_start_viewpoint
        bin_end_viewpoint_fixated = len(data_list) - (bin_end_viewpoint - bin_end_viewpoint_fixated) +1
        # log.debug('<<<<<bin_start_viewpoint_fixated {}, bin_end_viewpoint_fixated {}\n'.format(bin_start_viewpoint_fixated, bin_end_viewpoint_fixated))
        # exit(0)

        if bin_start_viewpoint_fixated > 0:
            data_list[:bin_start_viewpoint_fixated] = data_list[bin_start_viewpoint_fixated]
        if bin_end_viewpoint_fixated < len(data_list):
            data_list[bin_end_viewpoint_fixated:] = data_list[bin_end_viewpoint_fixated]



        rbz_score_data = pViewpointObj.rbz_score(data_list, _backgroundModelData, _backgroundModelSEM)
        if rbz_score_data is None:
            continue
        if bin_start_viewpoint_fixated > 0:
            rbz_score_data[:bin_start_viewpoint_fixated] = rbz_score_data[bin_start_viewpoint_fixated]
        if bin_end_viewpoint_fixated < len(rbz_score_data):
            rbz_score_data[bin_end_viewpoint_fixated:] = rbz_score_data[bin_end_viewpoint_fixated]

        # log.debug('rbz_score_data {}'.format(rbz_score_data[:10]))

        # else:
            # continue
        
        interaction_data = pViewpointObj.createInteractionFileData(referencePoint, referencePoint[0],
                                                                    region_start, region_end, data_list, data_list_raw,
                                                                    pGeneList[i])

        referencePointString = '_'.join(str(j) for j in referencePoint)

        region_start_in_units = utilities.in_units(region_start)
        region_end_in_units = utilities.in_units(region_end)

        header_information = '\t'.join([pMatrix, referencePointString, str(region_start_in_units), str(region_end_in_units), pGeneList[i]])
        header_information += '\n# ChrViewpoint\tStart\tEnd\tGene\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\trbz-score\tRaw\n#'
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
