"""chicExportData: a thin layer over the C++ executable chicExportData.

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
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import hicmatrix.HiCMatrix as hm
from hicexplorer import utilities
from hicexplorer._version import __version__
from .lib import Viewpoint

import h5py
import io
import tarfile
from contextlib import closing
import pyBigWig
from collections import OrderedDict

from tempfile import mkdtemp
import shutil


def parse_arguments(args=None):
    """The argparse parser of chicExportData, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('chicExportData')


def exportData(pFileList, pArgs, pViewpointObject, pDecimalPlace, pChromosomeSizes, pBackgroundData, pFileType, pQueue):

    file_list = []
    file_content_list = []

    file_ending = '.txt' if pArgs.outputFileType == 'txt' else '.bigwig'
    try:
        if pFileType == 'interactions' or pFileType == 'significant':
            header_information = '# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRelative Interactions\tp-value\tx-fold\tRaw\n'

            for file in pFileList:

                if pArgs.outputFileType == 'bigwig':
                    chromosome_name = []
                    start = []
                    end = []
                    values = []
                    relative_distance = {}
                for sample in file:
                    data = pViewpointObject.readInteractionFile(pArgs.file, sample)
                    key_list = sorted(list(data[1].keys()))

                    if pArgs.outputFileType == 'txt':

                        file_content_string = header_information
                        for key in key_list:
                            file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, float) else str(x) for x in data[1][key]) + '\n'
                        # breakpoint()
                    else:
                        for key in key_list:
                            chromosome_name.append(str(data[1][key][0]))
                            start.append(int(data[1][key][1]))
                            end.append(int(data[1][key][2]))
                            if pArgs.outputValueBigwig == 'relative-interactions':
                                values.append(float(data[1][key][6]))
                            elif pArgs.outputValueBigwig == 'p-value':
                                values.append(float(data[1][key][7]))
                            elif pArgs.outputValueBigwig == 'x-fold':
                                values.append(float(data[1][key][8]))
                            elif pArgs.outputValueBigwig == 'raw':
                                values.append(float(data[1][key][9]))

                            relative_distance[data[1][key][5]] = [str(data[1][key][0]), int(data[1][key][1]), int(data[1][key][2])]
                        header = [(chromosome_name[0], pChromosomeSizes[chromosome_name[0]])]

                        if pArgs.backgroundModelFile is not None:
                            key_list_background = sorted(list(pBackgroundData.keys()))
                            chromosome_name_background = []
                            start_background = []
                            end_background = []
                            value_background = []

                            for key in key_list_background:
                                if key in relative_distance:
                                    chromosome_name_background.append(relative_distance[key][0])
                                    start_background.append(relative_distance[key][1])
                                    end_background.append(relative_distance[key][2])
                                    value_background.append(float(pBackgroundData[key][0]))

                if pArgs.outputFileType == 'txt':
                    file_content_list.append(file_content_string)
                    file_name = '_'.join(sample) + '_' + pFileType + file_ending
                else:
                    if pArgs.backgroundModelFile is not None:
                        file_content_list.append([[header, chromosome_name, start, end, values], [header, chromosome_name_background, start_background, end_background, value_background]])

                        file_name = ['_'.join(sample) + file_ending, 'background_' + '_'.join(sample) + '_' + pFileType + file_ending]

                    else:
                        file_content_list.append([[header, chromosome_name, start, end, values]])
                        file_name = ['_'.join(sample) + '_' + pFileType + file_ending]

                file_list.append(file_name)
        elif pFileType == 'target':
            # targetList, present_genes = pViewpointObject.readTargetHDFFile(pArgs.file)
            # header_information = '# Chromosome\tStart\tEnd\n'

            for targetFile in pFileList:
                targetFileHDF5Object = h5py.File(pArgs.file, 'r')
                target_object = targetFileHDF5Object['/'.join(targetFile)]
                chromosome = target_object.get('chromosome')[()]
                start_list = list(target_object['start_list'][:])
                end_list = list(target_object['end_list'][:])
                targetFileHDF5Object.close()
                chromosome = [chromosome] * len(start_list)

                target_regions = list(zip(chromosome, start_list, end_list))
                file_content_string = ''
                # key_list = sorted(list(data[1].keys()))
                for region in target_regions:
                    file_content_string += '\t'.join(x.decode('utf-8') for x in region) + '\n'
                file_content_list.append(file_content_string)
                file_name = '_'.join(targetFile) + '_' + pFileType + '.txt'
                file_list.append(file_name)

        elif pFileType == 'aggregate':
            header_information = '# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRaw\n'

            for file in pFileList:
                for sample in file:
                    line_content, data = pViewpointObject.readAggregatedFileHDF(pArgs.file, sample)
                    file_content_string = header_information
                    for line in line_content:
                        file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, float) else str(x) for x in line) + '\n'
                    # breakpoint()

                    file_content_list.append(file_content_string)

                    file_name = '_'.join(sample) + '_' + pFileType + '.txt'
                    file_list.append(file_name)

        elif pFileType == 'differential':
            header_information = '# Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of interactions 2\ttarget_2 raw\tp-value\n'

            for file in pFileList:
                # accepted_list, all_list, rejected_list
                item_classification = ['accepted', 'all', 'rejected']
                line_content = pViewpointObject.readDifferentialFile(pArgs.file, file)
                for i, item in enumerate(line_content):
                    file_content_string = header_information

                    for line in item:
                        file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, float) else str(x) for x in line) + '\n'
                    # breakpoint()

                    file_content_list.append(file_content_string)
                    file_name = '_'.join(file) + '_' + item_classification[i] + '_' + pFileType + '.txt'
                    file_list.append(file_name)

    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return

    pQueue.put([file_list, file_content_list])
    return


def main(args=None):
    """Run the C++ implementation of chicExportData with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('chicExportData')(args)
