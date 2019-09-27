import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from io import StringIO
from multiprocessing import Process, Queue
import time
import numpy as np
import cooler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import ranksums
import pandas as pd
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Plots the relation between short and long range interactions as boxplots and if more than one matrix is given, p-values of the distributions are computed. 
An example usage is:
$ hicTadDensity -m hmec_10kb.cool -td domains.bed -o tad_densities_output.txt
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to compute the TAD densities on',
                                required=True)
    parserRequired.add_argument('--tadDomains', '-td',
                                help='The TADs domain file computed by hicFindTADs.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--outFileName', '-o',
                           help='File the densities are written to, p-values are only computed if at least two matrices are given.',
                           default='densities.txt')
    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes to include in the analysis. If not set, all chromosomes are included.',
                           nargs='+')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--mode', '-mo',
                           help='Compute the density via binary: number of non-zero elements / number of all elements;  the raw sum interaction counts; or sum of interactions / number of all elements',
                           choices=['binary', 'raw_interactions', 'interactions_size'],
                           default='interactions_size')
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2]]

    return domains_df


def computeRegionsTADs(pMatrix, pDomainList, pCoolOrH5, pI, pFirst, pLast, pRow):
    length_domains_list = len(pDomainList)
    matrix, intertad_left, intertad_right = [None] * 3
    if not pFirst:
        chromosom = pDomainList[pI - 1][0]
        start = pDomainList[pI - 1][1]
    else:
        chromosom = pDomainList[pI][0]
        start = pDomainList[pI][1]
    if not pLast:
        end = pDomainList[pI + 1][2]
    else:
        end = pDomainList[pI][2]
    midpos = pRow[1] + ((pRow[2] - pRow[1]) / 2)

    if pCoolOrH5:

        # get intra-TAD data
        hic_matrix = hm.hiCMatrix(
            pMatrixFile=pMatrix, pChrnameList=[str(pRow[0]) + ':' + str(pRow[1]) + '-' + str(pRow[2])])
        matrix = hic_matrix.matrix

        hic_matrix_inter_tad = hm.hiCMatrix(
            pMatrixFile=pMatrix, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])

        matrix_inter_tad = hic_matrix_inter_tad.matrix

    else:
        # in case of h5 pMatrixTarget is already a HiCMatrix object
        hic_matrix = pMatrix
        hic_matrix_inter_tad = pMatrix
        indices_target = hic_matrix.getRegionBinRange(str(pRow[0]), pRow[1], pRow[2])

        matrix_target = hic_matrix.matrix[indices_target[0]:indices_target[1], indices_target[0]:indices_target[1]].toarray()
        matrix_inter_tad = pMatrix.matrix

    tad_midpoint = hic_matrix_inter_tad.getRegionBinRange(str(pRow[0]), midpos, midpos)[0]

    if not pFirst:
        # get index position left tad with tad
        left_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), pRow[1], pRow[1])[0]
    if pCoolOrH5:
        outer_left_boundary_index = 0
        outer_right_boundary_index = -1

    else:
        outer_left_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

        outer_right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]

    if not pLast:
        # get index position left tad with tad
        right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), pRow[2], pRow[2])[0]

    if pFirst:
        intertad_right = matrix_inter_tad[tad_midpoint:right_boundary_index, tad_midpoint:outer_right_boundary_index]

    elif pLast:
        intertad_left = matrix_inter_tad[outer_left_boundary_index:tad_midpoint, left_boundary_index:tad_midpoint]

    else:
        intertad_left = matrix_inter_tad[outer_left_boundary_index:tad_midpoint, left_boundary_index:tad_midpoint]
        intertad_right = matrix_inter_tad[tad_midpoint:right_boundary_index, tad_midpoint:outer_right_boundary_index]
    return matrix, intertad_left, intertad_right


def computeDensityTADs(pMatrix, pDomainList, pCoolOrH5, pFirst, pLast, pThreads, pMode, pQueue):
    density_inter_left_list = []
    density_inter_right_list = []
    density_intra_list = []
    rows = []
    length_domains_list = len(pDomainList)

    if pThreads > 1:
        if pFirst:
            startValue = 0
        else:
            startValue = 1
        if pLast:
            endValue = len(pDomainList)
        else:
            endValue = len(pDomainList) - 1
    else:
        startValue = 0
        endValue = len(pDomainList)

    for i, row in enumerate(pDomainList[startValue:endValue]):
        # get intra / inter-tad data
        first = pFirst
        last = pLast
        if startValue == 0 and i == 0:
            first = True
            last = False
        elif i + startValue == len(pDomainList) - 1:
            first = False
            last = True
        else:
            last = False
            first = False

        matrix, intertad_left, intertad_right = computeRegionsTADs(pMatrix, pDomainList, pCoolOrH5, i + startValue, first, last, row)
        if first:

            if pMode == 'binary':
                density_inter_right = intertad_right.count_nonzero() / (intertad_right.shape[0] * intertad_right.shape[1])
            elif pMode == 'raw_interactions':
                density_inter_right = np.sum(intertad_right.data)
            else:
                density_inter_right = np.sum(intertad_right.data) / (intertad_right.shape[0] * intertad_right.shape[1])
            density_inter_left = -1
        elif last:
            if pMode == 'binary':
                density_inter_left = intertad_left.count_nonzero() / (intertad_left.shape[0] * intertad_left.shape[1])
            elif pMode == 'raw_interactions':
                density_inter_left = np.sum(intertad_left.data)
            else:
                density_inter_left = np.sum(intertad_left.data) / (intertad_left.shape[0] * intertad_left.shape[1])
            density_inter_right = -1
        else:
            if pMode == 'binary':
                density_inter_left = intertad_left.count_nonzero() / (intertad_left.shape[0] * intertad_left.shape[1])
                density_inter_right = intertad_right.count_nonzero() / (intertad_right.shape[0] * intertad_right.shape[1])
            elif pMode == 'raw_interactions':
                density_inter_left = np.sum(intertad_left.data)
                density_inter_right = np.sum(intertad_right.data)
            else:
                density_inter_left = np.sum(intertad_left.data) / (intertad_left.shape[0] * intertad_left.shape[1])
                density_inter_right = np.sum(intertad_right.data) / (intertad_right.shape[0] * intertad_right.shape[1])

        if pMode == 'binary':
            density_intra = matrix.count_nonzero() / (matrix.shape[0] * matrix.shape[1])
        elif pMode == 'raw_interactions':
            density_intra = np.sum(matrix.data)
        else:
            density_intra = np.sum(matrix.data) / (matrix.shape[0] * matrix.shape[1])

        density_inter_left_list.append(density_inter_left)
        density_inter_right_list.append(density_inter_right)
        density_intra_list.append(density_intra)

    pQueue.put([density_inter_left_list, density_inter_right_list, density_intra_list])
    return


def main(args=None):

    args = parse_arguments().parse_args(args)
    domains_df = readDomainBoundaries(args.tadDomains)
    domains = domains_df.values.tolist()
    tads_list = []
    matrix = args.matrix

    is_cooler = check_cooler(matrix)
    if not is_cooler:
        hic_matrix = hm.hiCMatrix(matrix)
    else:
        hic_matrix = matrix

    domainsListPerThread = [None] * args.threads
    tadResultListPerThread = [None] * args.threads

    numberOfDomainsPerThread = len(domains) // args.threads
    all_data_collected = False
    queue = [None] * args.threads
    process = [None] * args.threads
    thread_done = [False] * args.threads
    for i in range(args.threads):

        if i == 0 and args.threads > 1:
            domainsListPerThread[i] = domains[(i * numberOfDomainsPerThread):((i + 1) * numberOfDomainsPerThread) + 1]
        elif i < args.threads - 1:
            domainsListPerThread[i] = domains[(i * numberOfDomainsPerThread) - 1:((i + 1) * numberOfDomainsPerThread) + 1]
        else:
            if args.threads == 1:
                domainsListPerThread[i] = domains
            else:
                domainsListPerThread[i] = domains[(i * numberOfDomainsPerThread) - 1:]

        queue[i] = Queue()
        if i == 0:
            first = True
        else:
            first = False
        if i == args.threads - 1:
            last = True
        else:
            last = False
        process[i] = Process(target=computeDensityTADs, kwargs=dict(
            pMatrix=hic_matrix,
            pDomainList=domainsListPerThread[i],
            pCoolOrH5=is_cooler,
            pFirst=first,
            pLast=last,
            pThreads=args.threads,
            pMode=args.mode,
            pQueue=queue[i]
        )
        )

        process[i].start()

    while not all_data_collected:
        for i in range(args.threads):
            if queue[i] is not None and not queue[i].empty():
                tadResultListPerThread[i] = queue[i].get()
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

    inter_left_list = [item for sublist in tadResultListPerThread for item in sublist[0]]
    inter_right_list = [item for sublist in tadResultListPerThread for item in sublist[1]]
    intra_list = [item for sublist in tadResultListPerThread for item in sublist[2]]

    with open(args.outFileName, 'w') as file:
        header = '# Created with HiCExplorer\'s hicTadDensity ' + __version__ + '\n'
        header += "# intra- and inter-tad densities; mode: {}\n#\n".format(args.mode)
        header += "# Chromosome\tstart\tend\tinter-left\tinter-right\tintra\n"

        file.write(header)

        for i, domain in enumerate(domains):
            file.write('{}\t{}\t{}\t{}\n'.format(domain, inter_left_list[i], inter_right_list[i], intra_list[i]))
