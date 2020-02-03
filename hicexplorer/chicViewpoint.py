import argparse
import sys
import os
import errno
from multiprocessing import Process, Queue
import time
import math
import logging
log = logging.getLogger(__name__)

import numpy as np

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from .lib import Viewpoint
from hicexplorer._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
Computes per input matrix all viewpoints which are defined in the reference points file. All files are stored in the folder defined by `--outputFolder`, and the files
are named by the name of the reference point, the sample name and the location of the reference point:

gene_matrix_name_chr_start_end.bed

If multiple reference points are used and the processing downstream should be automated via batch processing mode, please activate `--writeFileNamesToFile`. All the file names will be written to this file; in the case of multiple samples two consecutive lines are considered as treatment vs. control for the differential analysis.

An example usage is:

$ chicViewpoint --matrices matrix1.cool matrix2.cool matrix3.cool --referencePoints referencePointsFile.bed --range 20000 40000 --outputFolder interactionFilesFolder

An example usage for batch mode is:

$ chicViewpoint --matrices matrix1.cool matrix2.cool matrix3.cool --referencePoints referencePointsFile.bed --range 20000 40000 --outputFolder interactionFilesFolder --writeFileNamesToFile interactionFile.txt

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
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module).',
                           required=False,
                           default=4,
                           type=int)
    parserOpt.add_argument('--averageContactBin',
                           help='Average the contacts of n bins via a sliding window approach to smooth the values and be less sensitive for outliers.',
                           type=int,
                           default=5)
    parserOpt.add_argument('--writeFileNamesToFile', '-w',
                           help='Set this parameter to have a file with all file names of the viewpoint files (useful only in batch processing mode).')

    parserOpt.add_argument('--fixateRange', '-fs',
                           help='Fixate range of background model starting at distance x. E.g. all values greater 500kb are set to the value of the 500kb bin.',
                           required=False,
                           default=500000,
                           type=int
                           )
    parserOpt.add_argument('--xFoldMaxValueNB', '-xfnb',
                           help='x-fold factor to increase the number of precomputed p-values per relative genomic distance. If set to 1, the maximal distance is used. ',
                           required=False,
                           default=1,
                           type=int
                           )
    parserOpt.add_argument('--allViewpointsList', '-avl',
                           help='Writes a file where all viewpoints all samples are sorted by the viewpoints.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--outputFolder', '-o',
                           help='This folder contains all created viewpoint files.',
                           required=False,
                           default='interactionFiles')
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


def compute_viewpoint(pViewpointObj, pArgs, pQueue, pReferencePoints, pGeneList, pMatrix, pBackgroundModel, pBackgroundModelRelativeInteractions, pOutputFolder):
    file_list = []

    for i, referencePoint in enumerate(pReferencePoints):
        # range of viewpoint with reference point in the middle in genomic units
        # get fixateRange for relative interaction computation denominator
        region_start_fixed, region_end_fixed, range_fixed = pViewpointObj.calculateViewpointRange(
            referencePoint, (pArgs.fixateRange, pArgs.fixateRange))

        intermediate_viewpoint = pViewpointObj.computeViewpoint(
            referencePoint, referencePoint[0], region_start_fixed, region_end_fixed)
        denominator_relative_interactions = np.sum(intermediate_viewpoint)

        # viewpoint data uses full range
        region_start, region_end, _range = pViewpointObj.calculateViewpointRange(
            referencePoint, pArgs.range)

        data_list = pViewpointObj.computeViewpoint(
            referencePoint, referencePoint[0], region_start, region_end)

        # background uses fixed range, handles fixate range implicitly by same range used in background computation

        _backgroundModelNBinom = pViewpointObj.interactionBackgroundData(
            pBackgroundModel, _range)

        background_relative_interaction = pViewpointObj.interactionBackgroundData(
            pBackgroundModelRelativeInteractions, _range).flatten()

        data_list_relative = data_list
        if len(data_list) != len(_backgroundModelNBinom):

            data_list, _backgroundModelNBinom, = adjustViewpointData(
                pViewpointObj, data_list, _backgroundModelNBinom, referencePoint, region_start, region_end)

        if len(data_list) != len(background_relative_interaction):
            _, background_relative_interaction = adjustViewpointData(
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
            _backgroundModelNBinom, data_list_raw)

        # add values if range is larger than fixate range

        region_start_range, region_end_range, _ = pViewpointObj.calculateViewpointRange(
            referencePoint, (pArgs.range[0], pArgs.range[1]))

        interaction_data = pViewpointObj.createInteractionFileData(referencePoint, referencePoint[0],
                                                                   region_start_range, region_end_range, data_list, data_list_raw,
                                                                   pGeneList[i], denominator_relative_interactions)

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
        file_list.append(matrix_name + '.bed')

        matrix_name = pOutputFolder + '/' + matrix_name
        pViewpointObj.writeInteractionFile(
            matrix_name, interaction_data, header_information, p_value_list, x_fold_list)

    pQueue.put(file_list)
    return


def main(args=None):
    args = parse_arguments().parse_args(args)

    viewpointObj = Viewpoint()

    referencePoints, gene_list = viewpointObj.readReferencePointFile(
        args.referencePoints)
    referencePointsPerThread = len(referencePoints) // args.threads
    queue = [None] * args.threads
    process = [None] * args.threads
    file_list = []
    background_model = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range)
    background_model_mean_values = viewpointObj.readBackgroundDataFile(
        args.backgroundModelFile, args.range, pMean=True)
    background_sum_of_densities_dict = viewpointObj.computeSumOfDensities(
        background_model, args, pXfoldMaxValue=args.xFoldMaxValueNB)

    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        viewpointObj.hicMatrix = hic_ma
        file_list_sample = [None] * args.threads
        all_data_collected = False

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
                pBackgroundModel=background_sum_of_densities_dict,
                pBackgroundModelRelativeInteractions=background_model_mean_values,
                pOutputFolder=args.outputFolder
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    file_list_ = queue[i].get()
                    file_list_sample[i] = file_list_
                    process[i].join()
                    process[i].terminate()
                    process[i] = None

            all_data_collected = True

            for i in range(args.threads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)
        file_list_sample = [item for sublist in file_list_sample for item in sublist]
        file_list.append(file_list_sample)

    log.debug('file_list {}'.format(file_list))
    if args.writeFileNamesToFile:
        with open(args.writeFileNamesToFile, 'w') as file:
            log.debug('len(file_list) {}'.format(len(file_list)))
            if len(file_list) > 1:
                for i, sample in enumerate(file_list):
                    for sample2 in file_list[i + 1:]:
                        for viewpoint, viewpoint2 in zip(sample, sample2):
                            file.write(viewpoint + '\n')
                            file.write(viewpoint2 + '\n')
            else:
                for viewpoint in file_list[0]:
                    file.write(viewpoint + '\n')
    if args.allViewpointsList:

        with open(args.writeFileNamesToFile + 'all', 'w') as file:
            if len(file_list) > 1:
                for i, sample in enumerate(file_list[0]):
                    file.write(sample + '\n')
                    for j in range(1, len(file_list)):
                        file.write(file_list[j][i] + '\n')
            else:
                for viewpoint in file_list[0]:
                    file.write(viewpoint + '\n')
