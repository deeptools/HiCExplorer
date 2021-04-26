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
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
chicExportData exports the data stored in the intermediate hdf5 files to text files per reference point.

"""
                                     )
    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--file', '-f',
                                help='path to the file which should be used for data export',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='Output tar.gz of the files. In case of --outputMode == geneName it is ignored.'
                           ' (Default: %(default)s).',
                           required=False,
                           default='data.tar.gz')

    parserOpt.add_argument('--fileType',
                           '-ft',
                           help=''
                           ' (Default: %(default)s).',
                           default='interaction',
                           choices=['interaction', 'significant', 'target', 'aggregated', 'differential']
                           )

    parserOpt.add_argument('--outputFileType',
                           '-oft',
                           help='Output file type can be set for all --fileTypes to txt; except \'interaction\' supports also bigwig'
                           ' (Default: %(default)s).',
                           default='txt',
                           choices=['txt', 'bigwig']
                           )
    parserOpt.add_argument('--outputMode',
                           '-om',
                           help='Output mode: Either all date is written or a gene name must be specified.'
                           ' (Default: %(default)s).',
                           default='all',
                           choices=['all', 'geneName']
                           )
    parserOpt.add_argument('--outputModeName',
                           '-omn',
                           help='ONLY valid if --outputMode geneName! Define the name of the gene',
                        #    default='',
                           )
    parserOpt.add_argument('--decimalPlaces',
                           help='Decimal places for all output floating numbers in the viewpoint files'
                           ' (Default: %(default)s).',
                           type=int,
                           default=12)
    parserOpt.add_argument('--chromosomeSizes', '-cs',
                           help=('File with the chromosome sizes for your genome. A tab-delimited two column layout \"chr_name size\" is expected'
                                 'Usually the sizes can be determined from the SAM/BAM input files, however, '
                                 'for cHi-C or scHi-C it can be that at the start or end no data is present. '
                                 'Please consider that this option causes that only reads are considered which are on the listed chromosomes.'
                                 'Use this option to guarantee fixed sizes. An example file is available via UCSC: '
                                 'http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes'),
                           type=argparse.FileType('r'),
                           metavar='txt file')
    parserOpt.add_argument('--backgroundModelFile', '-bmf',
                           help='Path to the background model file. Required only for fileType=interactions and outputFileTypeBigwig.',
                           required=False)
    parserOpt.add_argument('--oneTargetFile', '-otf',
                           help='Compile all target files to one. Applies only if --fileType is target',
                           required=False,
                           action='store_true')
    parserOpt.add_argument('--range',
                           help='Defines the region upstream and downstream of a reference point which should be included. '
                           'Format is --range upstream downstream, e.g.: --range 500000 500000 plots 500kb up- and 500kb downstream. '
                           'This value should not exceed the range used in the other chic-tools. Applies only for interaction files in the combination with bigwig and a background model file!',
                           required=False,
                           type=int,
                           nargs=2)
    parserOpt.add_argument('--outputValueBigwig',
                           '-ovb',
                           help='Select which value the bigwig file should contain: \'relative-interactions\', \'p-value\', \'x-fold\', \'raw\''
                           ' (Default: %(default)s).',
                           default='relative-interactions',
                           choices=['relative-interactions', 'p-value', 'x-fold', 'raw']
                           )
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def exportData(pFileList, pArgs, pViewpointObject, pDecimalPlace, pChromosomeSizes, pBackgroundData, pQueue):

    file_list = []
    file_content_list = []

    file_ending = '.txt' if pArgs.outputFileType == 'txt' else '.bigwig'
    log.debug('len(pFileList) {}'.format(len(pFileList)))
    try:
        if pArgs.fileType == 'interaction' or pArgs.fileType == 'significant':
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
                            file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, np.float) else str(x) for x in data[1][key]) + '\n'
                    else:
                        for key in key_list:
                            # log.debug('len(data[1][key]) {}'.format(data[1][key]))
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
                    file_name = '_'.join(sample) + file_ending
                else:
                    if pArgs.backgroundModelFile is not None:
                        file_content_list.append([[header, chromosome_name, start, end, values], [header, chromosome_name_background, start_background, end_background, value_background]])

                        file_name = ['_'.join(sample) + file_ending, 'background_' + '_'.join(sample) + file_ending]

                    else:
                        file_content_list.append([[header, chromosome_name, start, end, values]])
                        file_name = ['_'.join(sample) + file_ending]

                file_list.append(file_name)
        elif pArgs.fileType == 'target':
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
                file_name = '_'.join(targetFile) + '.txt'
                file_list.append(file_name)

        elif pArgs.fileType == 'aggregated':
            header_information = '# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative position\tRaw\n'

            for file in pFileList:
                for sample in file:
                    line_content, data = pViewpointObject.readAggregatedFileHDF(pArgs.file, sample)
                    file_content_string = header_information
                    for line in line_content:
                        file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, np.float) else str(x) for x in line) + '\n'
                    file_content_list.append(file_content_string)

                    file_name = '_'.join(sample) + '.txt'
                    file_list.append(file_name)

        elif pArgs.fileType == 'differential':
            header_information = '# Chromosome\tStart\tEnd\tGene\tRelative distance\tsum of interactions 1\ttarget_1 raw\tsum of interactions 2\ttarget_2 raw\tp-value\n'

            for file in pFileList:
                # accepted_list, all_list, rejected_list
                item_classification = ['accepted', 'all', 'rejected']
                line_content = pViewpointObject.readDifferentialFile(pArgs.file, file)
                for i, item in enumerate(line_content):
                    file_content_string = header_information

                    for line in item:
                        file_content_string += '\t'.join('{:.{decimal_places}f}'.format(x, decimal_places=pDecimalPlace) if isinstance(x, np.float) else str(x) for x in line) + '\n'
                    file_content_list.append(file_content_string)
                    file_name = '_'.join(file) + '_' + item_classification[i] + '.txt'
                    file_list.append(file_name)

    except Exception as exp:
        log.debug("FAIL: {}".format(str(exp) + traceback.format_exc()))
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return

    pQueue.put([file_list, file_content_list])
    log.debug('RETRUN')
    return


def main(args=None):
    args = parse_arguments().parse_args(args)
    viewpointObj = Viewpoint()

    fileList = []
    chromosome_sizes = None
    background_dict = None

    if args.outputMode == 'geneName' and args.outputModeName is None:
        log.error('Output mode is \'geneName\'. Please specify a gene name via --outputModeName too!')
        exit(1)

    if args.outputFileType == 'bigwig':
        if args.fileType != 'interaction':
            log.error('Only file type \'interaction\' supports bigwig. Exiting.')
            exit(1)
        if args.backgroundModelFile is not None:
            if args.backgroundModelFile:
                background_dict = viewpointObj.readBackgroundDataFile(args.backgroundModelFile, args.range, args.range[1], pMean=True)
        else:
            log.error('Please define a background file via --backgroundModelFile')
            exit(1)
        if args.chromosomeSizes is not None:
            chromosome_sizes = OrderedDict()
            with open(args.chromosomeSizes.name, 'r') as file:
                file_ = True
                while file_:
                    file_ = file.readline().strip()
                    if file_ != '':
                        line_split = file_.split('\t')
                        chromosome_sizes[line_split[0]] = int(line_split[1])
        else:
            log.error('Bigwig files require the argument \'--chromosomeSizes\'. Exiting.')
            exit(1)

        if args.range is None:
            log.error('Bigwig files require the argument \'--range upstream downstream\'. Exiting.')
            exit(1)
    # read hdf file
    fileHDF5Object = h5py.File(args.file, 'r')
    keys_file = list(fileHDF5Object.keys())

    if args.fileType == 'interaction' or args.fileType == 'significant':

        # if len(keys_file) > 0:
        if args.outputMode == 'all':
            for i, sample in enumerate(keys_file):

                matrix_obj1 = fileHDF5Object[sample]
                chromosomeList1 = sorted(list(matrix_obj1.keys()))
                chromosomeList1.remove('genes')
                for chromosome1 in chromosomeList1:
                    geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                    for gene1 in geneList1:
                        fileList.append([[sample, chromosome1, gene1]])
        else:
            for i, sample in enumerate(keys_file):
                matrix_obj1 = fileHDF5Object[sample]['genes']
                chromosomeList1 = sorted(list(matrix_obj1.keys()))
                gene_name = args.outputModeName 
                counter = 1
                while gene_name in chromosomeList1:
                    fileList.append([[sample, 'genes', gene_name]])
                    gene_name = args.outputModeName + '_' + str(counter)
                    counter += 1

                # for chromosome1 in chromosomeList1:
                #     geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                #     for gene1 in geneList1:
                #         fileList.append([[sample, chromosome1, gene1]])
    elif args.fileType == 'target':

        if args.outputMode == 'all':
            log.debug('foo')
            for outer_matrix in keys_file:
                inner_matrix_object = fileHDF5Object[outer_matrix]
                keys_inner_matrices = list(inner_matrix_object.keys())
                for inner_matrix in keys_inner_matrices:
                    inner_object = inner_matrix_object[inner_matrix]
                    gene_object = inner_object['genes']
                    keys_genes = list(gene_object.keys())
                    for gen in keys_genes:
                        fileList.append([outer_matrix, inner_matrix, 'genes', gen])
        else:
            log.debug('huh')

            for outer_matrix in keys_file:
                inner_matrix_object = fileHDF5Object[outer_matrix]
                keys_inner_matrices = list(inner_matrix_object.keys())
                for inner_matrix in keys_inner_matrices:
                    inner_object = inner_matrix_object[inner_matrix]['genes']
                    keys_genes = list(inner_object.keys())
                    gene_name = args.outputModeName 
                    counter = 1
                    while gene_name in keys_genes:
                        fileList.append([outer_matrix, inner_matrix, 'genes', gene_name])
                        gene_name = args.outputModeName + '_' + str(counter)
                        counter += 1
                
    elif args.fileType == 'aggregated':

        if args.outputMode == 'all':
            for i, combinationOfMatrix in enumerate(keys_file):
                keys_matrix_intern = list(fileHDF5Object[combinationOfMatrix].keys())
                if len(keys_matrix_intern) == 0:
                    continue
                matrix1 = keys_matrix_intern[0]
                matrix2 = keys_matrix_intern[1]

                matrix_obj1 = fileHDF5Object[combinationOfMatrix + '/' + matrix1]
                matrix_obj2 = fileHDF5Object[combinationOfMatrix + '/' + matrix2]

                chromosomeList1 = sorted(list(matrix_obj1.keys()))
                chromosomeList2 = sorted(list(matrix_obj2.keys()))
                chromosomeList1.remove('genes')
                chromosomeList2.remove('genes')
                for chromosome1, chromosome2 in zip(chromosomeList1, chromosomeList2):
                    geneList1 = sorted(list(matrix_obj1[chromosome1].keys()))
                    geneList2 = sorted(list(matrix_obj2[chromosome2].keys()))

                    for gene1, gene2 in zip(geneList1, geneList2):
                        fileList.append([[combinationOfMatrix, matrix1, chromosome1, gene1], [combinationOfMatrix, matrix2, chromosome2, gene2]])
        else:
            for i, combinationOfMatrix in enumerate(keys_file):
                keys_matrix_intern = list(fileHDF5Object[combinationOfMatrix].keys())
                if len(keys_matrix_intern) == 0:
                    continue
                matrix1 = keys_matrix_intern[0]
                matrix2 = keys_matrix_intern[1]

                matrix_obj1 = fileHDF5Object[combinationOfMatrix + '/' + matrix1]['genes']
                matrix_obj2 = fileHDF5Object[combinationOfMatrix + '/' + matrix2]['genes']

                chromosomeList1 = sorted(list(matrix_obj1.keys()))
                chromosomeList2 = sorted(list(matrix_obj2.keys()))
                # chromosomeList1.remove('genes')
                # chromosomeList2.remove('genes')
                gene_name = args.outputModeName 
                counter = 1
                while gene_name in chromosomeList1 and gene_name in chromosomeList2:
                    fileList.append([[combinationOfMatrix, matrix1, 'genes', gene_name], [combinationOfMatrix, matrix2, 'genes', gene_name]])
                    gene_name = args.outputModeName + '_' + str(counter)
                    counter += 1
                
    elif args.fileType == 'differential':
        if args.outputMode == 'all':
            for outer_matrix in keys_file:
                inner_matrix_object = fileHDF5Object[outer_matrix]
                keys_inner_matrices = list(inner_matrix_object.keys())
                for inner_matrix in keys_inner_matrices:
                    inner_object = inner_matrix_object[inner_matrix]
                    chromosomeList = sorted(list(inner_object.keys()))
                    chromosomeList.remove('genes')
                    for chromosome in chromosomeList:
                        geneList = sorted(list(inner_object[chromosome].keys()))

                        for gene in geneList:
                            fileList.append([outer_matrix, inner_matrix, chromosome, gene])
        else:
            for outer_matrix in keys_file:
                inner_matrix_object = fileHDF5Object[outer_matrix]
                keys_inner_matrices = list(inner_matrix_object.keys())
                for inner_matrix in keys_inner_matrices:
                    inner_object = inner_matrix_object[inner_matrix]['genes']
                    chromosomeList = sorted(list(inner_object.keys()))
                    gene_name = args.outputModeName 
                    counter = 1
                    while gene_name in chromosomeList:
                        fileList.append([outer_matrix, inner_matrix, 'genes', gene_name])
                        # fileList.append([outer_matrix, inner_matrix, 'genes', gene_name])
                        gene_name = args.outputModeName + '_' + str(counter)
                        counter += 1

    log.debug('len(fileList) {}'.format(len(fileList)))
    fileHDF5Object.close()

    filesPerThread = len(fileList) // args.threads

    all_data_collected = False
    thread_data = [None] * args.threads
    file_name_list = [None] * args.threads

    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    fail_flag = False
    fail_message = ''

    for i in range(args.threads):

        if i < args.threads - 1:
            fileListPerThread = fileList[i * filesPerThread:(i + 1) * filesPerThread]
        else:
            fileListPerThread = fileList[i * filesPerThread:]
        queue[i] = Queue()

        process[i] = Process(target=exportData, kwargs=dict(
            pFileList=fileListPerThread,
            pArgs=args,
            pViewpointObject=viewpointObj,
            pDecimalPlace=args.decimalPlaces,
            pChromosomeSizes=chromosome_sizes,
            pBackgroundData=background_dict,
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
                    log.debug('fail flag')
                else:
                    file_name_list[i], thread_data[i] = return_content
                    log.debug('return content')

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

    thread_data = [item for sublist in thread_data for item in sublist]
    file_name_list = [item for sublist in file_name_list for item in sublist]
    log.debug('len(thread_data) {}'.format(len(thread_data)))
    log.debug('len(file_name_list) {}'.format(len(file_name_list)))


    if len(thread_data) == 0:
        log.error('Contains not the requested data!')
        exit(1)
    if args.outputFileType == 'txt':
        if args.outputMode == 'geneName':
            for i, file_content_string in enumerate(thread_data):
                with open(file_name_list[i], "w") as file:
                    file.write(file_content_string)
        else:        
            with tarfile.open(args.outFileName, "w:gz") as tar:

                if not args.oneTargetFile and not args.fileType == 'target':
                    log.debug('correct sub')
                    for i, file_content_string in enumerate(thread_data):

                        tar_info = tarfile.TarInfo(name=file_name_list[i])
                        tar_info.mtime = time.time()
                        file_content_string = file_content_string.encode('utf-8')
                        tar_info.size = len(file_content_string)
                        file = io.BytesIO(file_content_string)
                        tar.addfile(tarinfo=tar_info, fileobj=file)
                else:
                    tar_info = tarfile.TarInfo(name='targets.tsv')
                    tar_info.mtime = time.time()
                    file_content_string_all = ''
                    for i, file_content_string in enumerate(thread_data):

                        # tar_info = tarfile.TarInfo(name=file_name_list[i])
                        file_content_string_all += file_content_string
                    
                    file_content_string_all = file_content_string_all.encode('utf-8')
                    tar_info.size = len(file_content_string_all)
                    file = io.BytesIO(file_content_string_all)
                    tar.addfile(tarinfo=tar_info, fileobj=file)

    elif args.outputFileType == 'bigwig':
        bigwig_folder = mkdtemp(prefix="bigwig_folder")
        for i, file_content in enumerate(thread_data):

            for j, file_list in enumerate(file_content):

                bw = pyBigWig.open(bigwig_folder + '/' + file_name_list[i][j], 'w')
                # # set big wig header
                bw.addHeader(file_list[0])

                bw.addEntries(file_list[1], file_list[2], ends=file_list[3], values=file_list[4])
                bw.close()

        with tarfile.open(args.outFileName, "w:gz") as tar_handle:
            for root, dirs, files in os.walk(bigwig_folder):
                for file in files:
                    tar_handle.add(os.path.join(root, file))

        if os.path.exists(bigwig_folder):
            try:
                shutil.rmtree(bigwig_folder)
            except OSError as e:
                log.error("Error: %s - %s." % (e.filename, e.strerror))
