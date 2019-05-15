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
                                     description='Plots a viewpoint computed by chicViewpoint.')

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--interactionFile', '-if',
                                help='path to the interaction files which should be used for plotting',
                                required=True,
                                nargs='+')

    # parserRequired.add_argument('--outFileName', '-o',
    #                             help='File name to save the image.',
    #                             required=True)
    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. '
                                'Format is --region upstream downstream',
                                required=True,
                                type=int,
                                default=[500000, 500000],
                                nargs=2)
    # parserRequired.add_argument('--outFileName', '-o',
    #                             help='File name to save the image. Is not used in batch mode.')
    # parserRequired.add_argument('--range',
    #                        help='Defines the region upstream and downstream of a reference point which should be included. '
    #                        'Format is --region upstream downstream',
    #                        required=True,
    #                        type=int,
    #                        nargs=2)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--backgroundModelFile', '-bmf',
                           help='path to the background file which should be used for plotting',
                           required=False)
    parserOpt.add_argument('--interactionFileFolder', '-iff',
                           help='Folder where the interaction files are stored in. Applies only for batch mode.',
                           required=False,
                           default='.')
    parserOpt.add_argument('--differentialTestResult', '-dif',
                            help='Path to the files which with the H0 rejected files to highlight the regions in the plot.',
                            required=False,
                            nargs='+')
    parserOpt.add_argument('--outputFolder', '-of',
                           help='Output folder of the files.',
                           required=False,
                           default='plots')
    parserOpt.add_argument('--outputFormat', '-format',
                           help='Output format of the plot. Ignored if outFileName is given.',
                           required=False,
                           default='png')
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300,
                           required=False)
    parserOpt.add_argument('--binResolution', '-r',
                           help='Resolution of the bin in genomic units. Values are usually e.g. 1000 for a 1kb, 5000 for a 5kb or 10000 for a 10kb resolution.',
                           type=int,
                           default=1000,
                           required=False)
    
    parserOpt.add_argument('--colorMapPvalue',
                           help='Color map to use for the p-value. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='RdYlBu')
    parserOpt.add_argument('--maxPValue', '-map',
                           help='Maximal value for p-value. Values above are set to this value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--minPValue', '-mp',
                           help='Minimal value for p-value. Values below are set to this value.',
                           type=float,
                           default=None)
    parserOpt.add_argument('--pValue', '-p',
                           help='Plot p-values as a colorbar',
                           choices=['heatmap', ''],
                           default=''
                           )
    parserOpt.add_argument('--useRawBackground', '-raw',
                           help='Use the raw background instead of a smoothed one.',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--outFileName', '-o',
                            help='File name to save the image. It is not used in batch mode.')
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
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def plot_images(pInteractionFileList, pHighlightDifferentialRegionsFileList, pBackgroundData, pArgs, pViewpointObj, pQueue=None):
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

        if pArgs.pValue == 'heatmap':
            gs = gridspec.GridSpec(1 + len(interactionFile), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *z_score_heights], width_ratios=[0.75, 0.25])
            gs.update(hspace=0.5, wspace=0.05)
            ax1 = plt.subplot(gs[0, 0])
            ax1.margins(x=0)
        else:
            ax1 = plt.subplot()
        colors = ['g', 'b', 'c', 'm', 'y', 'k']
        background_plot = True
        data_plot_label = None
        for i, interactionFile_ in enumerate(interactionFile):
            # log.debug('interactionFile_ {}'.format(interactionFile_))
            header, data, background_data_plot, data_background_mean, p_values, interaction_file_data_raw, viewpoint_index = pViewpointObj.getDataForPlotting(pArgs.interactionFileFolder + '/' + interactionFile_, pArgs.range, pBackgroundData)
            if len(data) <= 1 or len(p_values) <= 1:
                log.warning('Only one data point in given range, no plot is created! Interaction file {} Range {}'.format(interactionFile_, pArgs.range))
                continue
            matrix_name, viewpoint, upstream_range, downstream_range, gene, _ = header.strip().split('\t')
            matrix_name = matrix_name[1:].split('.')[0]
            # if number_of_data_points < len(data):
            number_of_data_points = len(data)
            highlight_differential_regions = None
            
            if pArgs.differentialTestResult:
                log.debug('pArgs.binResolution {}'.format(pArgs.binResolution))
                highlight_differential_regions = pViewpointObj.readRejectedFile(pHighlightDifferentialRegionsFileList[j], viewpoint_index, pArgs.binResolution, pArgs.range)
            if data_plot_label:
                data_plot_label += pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name, pHighlightRegion=highlight_differential_regions)
            else:
                data_plot_label = pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name, pHighlightRegion=highlight_differential_regions)

            # if background_plot:
            #     if background_data_plot is not None and data_background_mean is not None:
            #         data_plot_label += pViewpointObj.plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot,
            #                                                             pBackgroundDataMean=data_background_mean)
            #     background_plot = False
            
            if pArgs.minPValue is not None or pArgs.maxPValue is not None:
                p_values = np.array(p_values, dtype=np.float32)
                p_values.clip(pArgs.minPValue, pArgs.maxPValue, p_values)
            if pArgs.pValue == 'heatmap':
                pViewpointObj.plotPValue(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pPValueData=p_values,
                                        pLabelText=gene + ': ' + matrix_name, pCmap=pArgs.colorMapPvalue,
                                        pFigure=fig,)
            elif pArgs.pValue == 'integrated':
                data_plot_label += pViewpointObj.plotViewpoint(pAxis=ax1, pData=p_values, pColor=colors[i % len(colors)], pLabelName=gene + ': ' + matrix_name + ' p-value')

            # pViewpointObj.writePlotData(interaction_file_data_raw, matrix_name + '_' + gene + '_raw_plot', pArgs.backgroundModelFile)

        if data_plot_label is not None:

            step_size = number_of_data_points // 10
            ticks = range(0, number_of_data_points, step_size)

            value_range = (pArgs.range[0]) // (5)
            x_labels = [str(-j // 1000) + 'kb' for j in range(pArgs.range[0], 0, -value_range)]
            x_labels.append('viewpoint')
            x_labels_ = [str(j // 1000) + 'kb' for j in range(value_range, pArgs.range[1] + 1, value_range)]
            x_labels.extend(x_labels_)
            # ax.set
            ax1.set_ylabel('Number of interactions')
            ax1.set_xticks(ticks)
            ax1.set_xticklabels(x_labels)

            # multiple legends in one figure
            data_legend = [label.get_label() for label in data_plot_label]
            ax1.legend(data_plot_label, data_legend, loc=0)

            

            if pArgs.outFileName:
                outFileName = pArgs.outFileName
            else:
                sample_prefix = interactionFile[0].split('/')[-1].split('_')[0] + '_' + interactionFile[1].split('/')[-1].split('_')[0]
                region_prefix = '_'.join(interactionFile[0].split('/')[-1].split('_')[1:6])
                outFileName = sample_prefix + '_' + region_prefix 
                outFileName = pArgs.outputFolder + '/' + outFileName + '.' + pArgs.outputFormat
            plt.savefig(outFileName, dpi=pArgs.dpi)
        plt.close(fig)
    
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
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    # if args.backgroundModelFile:
    #     background_data = viewpointObj.readBackgroundDataFile(args.backgroundModelFile, args.range, args.useRawBackground)
   
    interactionFileList = []
    highlightDifferentialRegionsFileList = []
    if args.batchMode:
        with open(args.interactionFile[0], 'r') as interactionFile:

            file_ = True
            while file_:
            # for line in fh.readlines():
                file_ = interactionFile.readline().strip()
                file2_ = interactionFile.readline().strip()
                if file_ != '' and file2_ != '':
                    interactionFileList.append((file_, file2_))
        if args.differentialTestResult:
            with open(args.differentialTestResult[0], 'r') as differentialTestFile:

                file_ = True
                while file_:
                # for line in fh.readlines():
                    file_ = differentialTestFile.readline().strip()
                    # file2_ = differentialTestFile.readline().strip()
                    if file_ != '':
                        highlightDifferentialRegionsFileList.append(file_)

        interactionFilesPerThread = len(interactionFileList) // args.threads
        all_data_collected = False
        queue = [None] * args.threads
        process = [None] * args.threads
        thread_done = [False] * args.threads
        # log.debug('matrix read, starting processing')
        log.debug('len(interactionFileList) {}'.format(len(interactionFileList)))
        log.debug('len(highlightDifferentialRegionsFileList) {}'.format(len(highlightDifferentialRegionsFileList)))
        # exit()

        for i in range(args.threads):

            if i < args.threads - 1:
                interactionFileListThread = interactionFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]
                highlightDifferentialRegionsFileListThread = highlightDifferentialRegionsFileList[i * interactionFilesPerThread:(i + 1) * interactionFilesPerThread]

            else:
                interactionFileListThread = interactionFileList[i * interactionFilesPerThread:]
                highlightDifferentialRegionsFileListThread = highlightDifferentialRegionsFileList[i * interactionFilesPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=plot_images, kwargs=dict(
                pInteractionFileList=interactionFileListThread,
                pHighlightDifferentialRegionsFileList=highlightDifferentialRegionsFileListThread,
                pBackgroundData=background_data,
                pArgs=args,
                pViewpointObj=viewpointObj,
                pQueue=queue[i]
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    background_data_thread = queue[i].get()
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
    else:
        interactionFileList = [args.interactionFile]
        highlightDifferentialRegionsFileList = args.differentialTestResult
        # plot_images(interactionFileList, highlightDifferentialRegionsFileList, args)
    
        plot_images(pInteractionFileList=interactionFileList,
                pHighlightDifferentialRegionsFileList=highlightDifferentialRegionsFileList,
                pBackgroundData=background_data,
                pArgs=args,
                pViewpointObj=viewpointObj)
