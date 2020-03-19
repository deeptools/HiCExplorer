import argparse
import sys
import os
import errno
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
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
chicPlotViewpoint plots one or many viewpoints with the average background model and the computed p-value per sample. In addition, it can highlight differential interactions of two samples and/or significant regions.

An example usage is:

`$ chicPlotViewpoint --interactionFile viewpoint1.txt viewpoint2.txt --range 500000 500000  --backgroundModelFile background_model.txt --pValue --outFileName viewpoint1_2.png --dpi 300`


In batch mode the list of file names and the folders containing the files need to be given:

`$ chicPlotViewpoint --interactionFile viewpoint_names.txt -interactionFileFolder viewpointFilesFolder --differentialTestResult rejected_H0.txt --differentialTestResultsFolder differentialFolder --range 500000 500000 --backgroundModelFile background_model.txt --pValue --outputFolder plotsFOlder --dpi 300 --threads 20`

"""
                                     )
    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for plotting',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream, e.g.: --region 500000 500000 plots 500kb up- and 500kb downstream. This value should not exceed the range used in the other chic-tools.',
                                required=True,
                                type=int,
                                default=[500000, 500000],
                                nargs=2)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--backgroundModelFile', '-bmf',
                           help='path to the background file which should be used for plotting',
                           required=False)
    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--differentialTestResult', '-dif',
                           help='Path to the H0 rejected files to highlight the regions in the plot.',
                           required=False,
                           nargs='+')
    parserOpt.add_argument('--significantInteractionFileFolder', '-siff',
                           help='Folder where the files with detected significant interactions are stored. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--differentialTestResultsFolder', '-diff',
                           help='Folder where the H0 rejected files are stored. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--significantInteractions', '-si',
                           help='Path to the files with detected significant interactions to highlight the regions in the plot.',
                           required=False,
                           nargs='+')
    parserOpt.add_argument('--plotSignificantInteractions', '-psi',
                           help='Highlights the significant interactions in the plot itself. If not set, only the p-values are updated',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--outputFolder', '-of',
                           help='Output folder of the files.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--outputFormat', '-format',
                           help='Output format of the plot.',
                           required=False,
                           default='png')
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image, if'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300,
                           required=False)
    parserOpt.add_argument('--binResolution', '-r',
                           help='Resolution of the bin in genomic units. Values are set as number of bases, e.g. 1000 for a 1kb, 5000 for a 5kb or 10000 for a 10kb resolution.',
                           type=int,
                           default=1000,
                           required=False)

    parserOpt.add_argument('--colorMapPvalue',
                           help='Color map to use for the p-value. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='RdYlBu')
    parserOpt.add_argument('--maxPValue', '-map',
                           help='Maximal value for p-value. Values above this threshold are set to this value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--minPValue', '-mp',
                           help='Minimal value for p-value. Values below this threshold are set to this value.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--pValue', '-p',
                           help='Plot p-values as a colorbar',
                           action='store_true'
                           )
    parserOpt.add_argument('--pValueSignificanceLevels', '-psl',
                           help='Highlight the p-values by the defined significance levels.',
                           type=float,
                           nargs='+'
                           )
    parserOpt.add_argument('--xFold', '-xf',
                           help='Plot x-fold region for the mean background.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save the image. Not used in batch mode.')
    parserOpt.add_argument('--batchMode', '-bm',
                           help='The given file for --interactionFile and or --targetFile contain a list of the to be processed files.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--plotSampleNumber', '-psn',
                           help='Number of samples per plot. Applies only in batch mode.',
                           required=False,
                           default=2,
                           type=int)
    parserOpt.add_argument('--colorList', '-cl',
                           help='Colorlist for the viewpoint lines.',
                           required=False,
                           default=['g', 'b', 'c', 'm', 'y', 'k'],
                           type=str,
                           nargs='+')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def plot_images(pInteractionFileList, pHighlightDifferentialRegionsFileList, pBackgroundData, pArgs, pViewpointObj, pSignificantRegionsFileList, pQueue=None):
    try:
        for j, interactionFile in enumerate(pInteractionFileList):
            number_of_rows_plot = len(interactionFile)
            matplotlib.rcParams.update({'font.size': 9})
            fig = plt.figure(figsize=(9.4, 4.8))

            z_score_heights = [0.07] * number_of_rows_plot
            viewpoint_height_ratio = 0.95 - (0.07 * number_of_rows_plot)
            if viewpoint_height_ratio < 0.4:
                viewpoint_height_ratio = 0.4
                _ratio = 0.6 / number_of_rows_plot
                z_score_heights = [_ratio] * number_of_rows_plot

            if pArgs.pValue:
                gs = gridspec.GridSpec(1 + len(interactionFile), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *z_score_heights], width_ratios=[0.75, 0.25])
                gs.update(hspace=0.5, wspace=0.05)
                ax1 = plt.subplot(gs[0, 0])
                ax1.margins(x=0)
            else:
                ax1 = plt.gca()
            colors = pArgs.colorList
            background_plot = True
            data_plot_label = None
            gene = ''
            for i, interactionFile_ in enumerate(interactionFile):
                if pArgs.interactionFileFolder != '.':
                    absolute_path_interactionFile_ = pArgs.interactionFileFolder + '/' + interactionFile_
                else:
                    absolute_path_interactionFile_ = interactionFile_

                header, data, background_data_plot, p_values, viewpoint_index_start, viewpoint_index_end = pViewpointObj.getDataForPlotting(absolute_path_interactionFile_, pArgs.range, pBackgroundData, pArgs.binResolution)
                # log.debug('data {}'.format(data))
                if len(data) <= 1 or len(p_values) <= 1:
                    log.warning('Only one data point in given range, no plot is created! Interaction file {} Range {}'.format(interactionFile_, pArgs.range))
                    continue
                matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.strip().split('\t')
                log.debug('Matrix_name {}'.format(matrix_name))
                matrix_name = os.path.basename(matrix_name)

                matrix_name = matrix_name.split('.')[0]
                log.debug('matrix_name {}'.format(matrix_name))
                # number_of_data_points = len(data)
                highlight_differential_regions = None
                significant_p_values = None
                significant_regions = None
                if pArgs.differentialTestResult:
                    if pArgs.differentialTestResultsFolder != '.':
                        differentialFilePath = pArgs.differentialTestResultsFolder + '/' + pHighlightDifferentialRegionsFileList[j]
                    else:
                        differentialFilePath = pHighlightDifferentialRegionsFileList[j]

                    highlight_differential_regions = pViewpointObj.readRejectedFile(differentialFilePath, viewpoint_index_start, viewpoint_index_end, pArgs.binResolution, pArgs.range, viewpoint)
                if pArgs.significantInteractions:
                    if pArgs.significantInteractionFileFolder != '.':
                        significantInteractionsFilePath = pArgs.significantInteractionFileFolder + '/' + pSignificantRegionsFileList[j][i]
                    else:
                        significantInteractionsFilePath = pSignificantRegionsFileList[j][i]
                    significant_regions, significant_p_values = pViewpointObj.readSignificantRegionsFile(significantInteractionsFilePath, viewpoint_index_start, viewpoint_index_end, pArgs.binResolution, pArgs.range, viewpoint)
                if not pArgs.plotSignificantInteractions:
                    significant_regions = None
                if data_plot_label:
                    data_plot_label += pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name, pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)
                else:
                    data_plot_label = pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name, pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)

                if background_plot:
                    # log.debug('background_data_plot {}'.format(len(background_data_plot)))
                    if background_data_plot is not None:
                        data_plot_label += pViewpointObj.plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot, pXFold=pArgs.xFold)
                    background_plot = False

                if pArgs.minPValue is not None or pArgs.maxPValue is not None:

                    p_values = np.array(p_values, dtype=np.float32)
                    if significant_p_values:
                        for location in significant_p_values:
                            for x in range(location[0], location[1]):
                                if x < len(p_values):
                                    p_values[x] = location[2]
                    p_values.clip(pArgs.minPValue, pArgs.maxPValue, p_values)
                # if pArgs.pValueSignificanceLevels:

                if pArgs.pValue:
                    pViewpointObj.plotPValue(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pPValueData=p_values,
                                             pLabelText=gene + ': ' + matrix_name, pCmap=pArgs.colorMapPvalue,
                                             pFigure=fig, pValueSignificanceLevels=pArgs.pValueSignificanceLevels)

            if data_plot_label is not None:

                ticks = []
                x_labels = []

                if pArgs.range[0] + pArgs.range[1] <= 2e6:
                    divisor_legend = 1e3
                    mod_legend = 2e5

                    if pArgs.range[0] + pArgs.range[1] <= 1e4:
                        mod_legend = 5e3
                    elif pArgs.range[0] + pArgs.range[1] <= 5e4:
                        mod_legend = 1e4
                    elif pArgs.range[0] + pArgs.range[1] <= 1e5:
                        mod_legend = 5e4
                    elif pArgs.range[0] + pArgs.range[1] <= 5e5:
                        mod_legend = 1e5
                    log.debug('divisor_legend {}'.format(divisor_legend))

                    unit = 'kb'
                elif pArgs.range[0] + pArgs.range[1] > 2e6:
                    divisor_legend = 1e6
                    mod_legend = 1e6
                    unit = 'Mb'

                for k, j in zip(range((pArgs.range[0])), range(pArgs.range[0], 1, -1)):
                    if j % mod_legend == 0:
                        x_labels.append(str(-int(j) // int(divisor_legend)) + unit)
                        ticks.append(k // pArgs.binResolution)
                x_labels.append('RP')
                ticks.append(pArgs.range[0] // pArgs.binResolution)

                referencepoint_index = ticks[-1]
                for k, j in zip(range(pArgs.range[1]), range(1, pArgs.range[1] + 1, 1)):
                    if j % mod_legend == 0:
                        x_labels.append(str(int(j) // int(divisor_legend)) + unit)
                        ticks.append(referencepoint_index + (k // pArgs.binResolution))

                # log.debug('labels: {}'.format(x_labels))
                ax1.set_ylabel('Number of interactions')
                ax1.set_xticks(ticks)
                ax1.set_xticklabels(x_labels)

                # multiple legends in one figure
                data_legend = [label.get_label() for label in data_plot_label]
                ax1.legend(data_plot_label, data_legend, loc=0)

                sample_prefix = ""
                if pArgs.outFileName:
                    if pArgs.outputFolder != '.':
                        outFileName = pArgs.outputFolder + '/' + pArgs.outFileName
                    else:
                        outFileName = pArgs.outFileName

                else:
                    for interactionFile_ in interactionFile:
                        sample_prefix += interactionFile_.split('/')[-1].split('_')[0] + '_'
                    if sample_prefix.endswith('_'):
                        sample_prefix = sample_prefix[:-1]
                    region_prefix = '_'.join(interactionFile[0].split('/')[-1].split('_')[1:4])
                    outFileName = gene + '_' + sample_prefix + '_' + region_prefix
                    if pArgs.outputFolder != '.':
                        outFileName = pArgs.outputFolder + '/' + outFileName

                if pArgs.outputFormat != outFileName.split('.')[-1]:
                    outFileName = outFileName + '.' + pArgs.outputFormat
                plt.savefig(outFileName, dpi=pArgs.dpi)
            plt.close(fig)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp))
        return
    if pQueue is None:
        return
    pQueue.put('done')
    return


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()
    background_data = None

    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if args.pValueSignificanceLevels:
        old = -100
        for element in args.pValueSignificanceLevels:
            if old < element:
                old = element
                continue
            else:
                log.error('--pValueSignificanceLevels levels need to increase: {}'.format(args.pValueSignificanceLevels))
                exit(1)
    if args.backgroundModelFile:
        background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile, args.range, pMean=True)

    interactionFileList = []
    highlightDifferentialRegionsFileList = []
    highlightSignificantRegionsFileList = []

    if args.batchMode:
        with open(args.interactionFile[0], 'r') as interactionFile:

            file_ = True
            while file_:
                lines = []
                for i in range(0, args.plotSampleNumber):
                    file_ = interactionFile.readline().strip()
                    if file_ != '':
                        lines.append(file_)
                interactionFileList.append(lines)
        if args.differentialTestResult:

            if args.differentialTestResult and args.plotSampleNumber != 2:
                log.warning('Cannot use differential data, only possible for two samples in one plot.')
                args.differentialTestResult = None
            else:
                with open(args.differentialTestResult[0], 'r') as differentialTestFile:

                    file_ = True
                    while file_:
                        file_ = differentialTestFile.readline().strip()
                        if file_ != '':
                            highlightDifferentialRegionsFileList.append(file_)
        if args.significantInteractions:
            with open(args.significantInteractions[0], 'r') as significantRegionsFile:

                file_ = True
                while file_:
                    lines = []
                    for i in range(0, args.plotSampleNumber):
                        file_ = significantRegionsFile.readline().strip()
                        if file_ != '':
                            lines.append(file_)
                    if len(lines) > 0:
                        highlightSignificantRegionsFileList.append(lines)
        interactionFilesPerThread = len(interactionFileList) // args.threads
        highlightSignificantRegionsFileListThread = len(highlightSignificantRegionsFileList) // args.threads

        all_data_collected = False
        queue = [None] * args.threads
        process = [None] * args.threads
        thread_done = [False] * args.threads
        fail_flag = False
        fail_message = ''

        for i in range(args.threads):

            if i < args.threads - 1:
                interactionFileListThread = interactionFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
                highlightDifferentialRegionsFileListThread = highlightDifferentialRegionsFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
                highlightSignificantRegionsFileListThread = highlightSignificantRegionsFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
            else:
                interactionFileListThread = interactionFileList[i * interactionFilesPerThread:]
                highlightDifferentialRegionsFileListThread = highlightDifferentialRegionsFileList[i * interactionFilesPerThread:]
                highlightSignificantRegionsFileListThread = highlightSignificantRegionsFileList[i * interactionFilesPerThread:]
            queue[i] = Queue()
            process[i] = Process(target=plot_images, kwargs=dict(
                pInteractionFileList=interactionFileListThread,
                pHighlightDifferentialRegionsFileList=highlightDifferentialRegionsFileListThread,
                pBackgroundData=background_data,
                pArgs=args,
                pViewpointObj=viewpointObj,
                pSignificantRegionsFileList=highlightSignificantRegionsFileListThread,
                pQueue=queue[i]
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    return_content = queue[i].get()
                    if 'Fail:' in return_content:
                        fail_flag = True
                        fail_message = return_content[6:]
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
    else:
        interactionFileList = [args.interactionFile]
        highlightDifferentialRegionsFileList = args.differentialTestResult
        highlightSignificantRegionsFileList = [args.significantInteractions]
        plot_images(pInteractionFileList=interactionFileList,
                    pHighlightDifferentialRegionsFileList=highlightDifferentialRegionsFileList,
                    pBackgroundData=background_data,
                    pArgs=args,
                    pViewpointObj=viewpointObj,
                    pSignificantRegionsFileList=highlightSignificantRegionsFileList)
