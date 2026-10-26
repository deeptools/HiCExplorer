"""chicPlotViewpoint: a thin layer over the C++ executable chicPlotViewpoint.

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
import traceback

import logging
log = logging.getLogger(__name__)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib as mpl
import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint

import h5py
import io
import tarfile
from contextlib import closing


def parse_arguments(args=None):
    """The argparse parser of chicPlotViewpoint, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicPlotViewpoint')


def plot_images(pInteractionFileList, pHighlightDifferentialRegionsFileList, pBackgroundData, pArgs, pViewpointObj, pSignificantRegionsFileList, pResolution, pQueue=None):
    images_array = []
    file_name_list = []
    try:

        for j, interactionFile in enumerate(pInteractionFileList):
            number_of_rows_plot = len(interactionFile)
            matplotlib.rcParams.update({'font.size': 9})
            fig = plt.figure(figsize=(9.4, 4.8), dpi=pArgs.dpi)
            FigureCanvas(fig)
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

            file_name = []
            for i, interactionFile_ in enumerate(interactionFile):
                file_name.append(interactionFile_[0])
                data, background_data_plot, p_values, viewpoint_index_start, viewpoint_index_end, viewpoint = pViewpointObj.getDataForPlotting(pArgs.interactionFile, interactionFile_, pArgs.range, pBackgroundData, pResolution)
                if len(data) <= 1 or len(p_values) <= 1:
                    log.warning('Only one data point in given range, no plot is created! Interaction file {} Range {}'.format(interactionFile_, pArgs.range))
                    continue
                highlight_differential_regions = None
                significant_p_values = None
                significant_regions = None
                if pArgs.differentialTestResult:
                    highlight_differential_regions = pViewpointObj.readRejectedFile(pArgs.differentialTestResult, pHighlightDifferentialRegionsFileList[j], viewpoint_index_start, viewpoint_index_end, pResolution, pArgs.range, viewpoint)
                if pArgs.significantInteractions:
                    significant_regions = None
                    significant_p_values = None
                    try:
                        if i < len(pSignificantRegionsFileList[j]):
                            #         if len(pSignificantRegionsFileList[j][i]) > 0:
                            significant_regions, significant_p_values = pViewpointObj.readSignificantRegionsFile(pArgs.significantInteractions, pSignificantRegionsFileList[j][i], viewpoint_index_start, viewpoint_index_end, pResolution, pArgs.range, viewpoint)
                    except Exception:
                        log.debug('len(pSignificantRegionsFileList) {}'.format(len(pSignificantRegionsFileList)))
                        log.debug('len(pSignificantRegionsFileList[j]) {}'.format(len(pSignificantRegionsFileList[j])))
                        log.debug('i {}'.format(i))

                if not pArgs.plotSignificantInteractions:
                    significant_regions = None
                if data_plot_label:
                    data_plot_label += pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=':'.join(interactionFile_), pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)
                else:
                    data_plot_label = pViewpointObj.plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=':'.join(interactionFile_), pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)

                if background_plot:
                    if background_data_plot is not None:
                        data_plot_label += pViewpointObj.plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot, pXFold=pArgs.xFold)
                    background_plot = False
                if pArgs.truncateZeroPvalues:
                    p_values = np.array(p_values, dtype=np.float32)
                    mask = p_values == 0.0
                    p_values[mask] = 1.0
                if pArgs.pValue and (pArgs.minPValue is not None or pArgs.maxPValue is not None):

                    p_values = np.array(p_values, dtype=np.float32)
                    if significant_p_values:
                        for location in significant_p_values:
                            for x in range(location[0], location[1]):
                                if x < len(p_values):
                                    p_values[x] = location[2]
                    p_values.clip(pArgs.minPValue, pArgs.maxPValue, p_values)

                if pArgs.pValue:
                    pViewpointObj.plotPValue(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pPValueData=p_values,
                                             pLabelText=':'.join(interactionFile_), pCmap=pArgs.colorMapPvalue,
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

                    unit = 'kb'
                elif pArgs.range[0] + pArgs.range[1] > 2e6:
                    divisor_legend = 1e6
                    mod_legend = 1e6
                    unit = 'Mb'

                for k, j in zip(range((pArgs.range[0])), range(pArgs.range[0], 1, -1)):
                    if j % mod_legend == 0:
                        x_labels.append(str(-int(j) // int(divisor_legend)) + unit)
                        ticks.append(k // pResolution)
                x_labels.append('RP')
                ticks.append(pArgs.range[0] // pResolution)

                referencepoint_index = ticks[-1]
                for k, j in zip(range(pArgs.range[1]), range(1, pArgs.range[1] + 1, 1)):
                    if j % mod_legend == 0:
                        x_labels.append(str(int(j) // int(divisor_legend)) + unit)
                        ticks.append(referencepoint_index + (k // pResolution))

                ax1.set_ylabel('Number of interactions')
                ax1.set_xticks(ticks)
                ax1.set_xticklabels(x_labels)

                # multiple legends in one figure
                data_legend = [label.get_label() for label in data_plot_label]
                ax1.legend(data_plot_label, data_legend, loc=0)

                bufferObject = io.BytesIO()
                plt.savefig(bufferObject, format=pArgs.outputFormat, dpi=300)
                images_array.append(bufferObject)
            plt.close(fig)
            file_name.append(interactionFile[0][2])
            file_name_list.append('_'.join(file_name))
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    if pQueue is None:
        return
    pQueue.put([images_array, file_name_list])
    return


def main(args=None):
    """Run the C++ implementation of chicPlotViewpoint with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicPlotViewpoint')(args)
