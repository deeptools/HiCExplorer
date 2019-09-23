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
$ hicPlotSVL -m hmec_10kb.cool nhek_10kb.cool
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The matrix (or multiple matrices) to use for the comparison',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--tadDomains', '-td',
                                help='The TADs domain file computed by hicFindTADs.',
                                required=True,
                                nargs='+')
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--plotFileName', '-pfn',
                           help='Plot name.',
                           default='plot.png')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File the densities are written to, p-values are only computed if at least two matrices are given.',
                           default='densities.txt')
    parserOpt.add_argument('--distance', '-d',
                           help='Distance which should be considered as short range. Default 2MB.',
                           default=2000000,
                           type=int)
    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes to include in the analysis. If not set, all chromosomes are included.',
                           nargs='+')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads. Using the python multiprocessing module. ',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'output is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300,
                           required=False)
    parserOpt.add_argument('--colorList', '-cl',
                           help='Colorlist for the boxplots.',
                           required=False,
                           default=['g', 'b', 'c', 'm', 'y', 'k'],
                           type=str,
                           nargs='+')
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser



def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2]]

    return domains_df


def computeRegionsTADs(pMatrix, pDomainList, pCoolOrH5, pI, pRow):
    length_domains_list = len(pDomainList)
    matrix, intertad_left, intertad_right = [None] * 3
    if pI - 1 >= 0:
        chromosom = pDomainList[pI - 1][0]
        start = pDomainList[pI - 1][1]
    else:
        chromosom = pDomainList[pI][0]
        start = pDomainList[pI][1]
    if pI + 1 < length_domains_list:
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

    if pI - 1 >= 0:
        # get index position left tad with tad
        left_boundary_index_target = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), pRow[1], pRow[1])[0]
    if pCoolOrH5:
        outer_left_boundary_index = 0
        outer_right_boundary_index = -1

    else:
        outer_left_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

        outer_right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]

    if pI + 1 < length_domains_list:
        # get index position left tad with tad
        right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), pRow[2], pRow[2])[0]

    if pI - 1 >= 0 and pI + 1 < length_domains_list:
        intertad_left = matrix_inter_tad[outer_left_boundary_index:tad_midpoint, left_boundary_index:tad_midpoint]
        intertad_right = matrix_inter_tad[tad_midpoint:right_boundary_index, tad_midpoint:outer_right_boundary_index]
    elif pI - 1 < 0 and pI + 1 < length_domains_list:
        intertad_right = matrix_inter_tad[tad_midpoint:right_boundary_index, tad_midpoint:outer_right_boundary_index]
    elif pI - 1 > 0 and pI + 1 >= length_domains_list:
        intertad_left = matrix_inter_tad[outer_left_boundary_index:tad_midpoint, left_boundary_index:tad_midpoint]

    return matrix, intertad_left, intertad_right

def computeDensityTADs(pMatrix, pDomainList, pCoolOrH5, pQueue):
    density_inter_left_list = []
    density_inter_right_list = []
    density_intra_list = []
    # p_values_list = []
    rows = []
    length_domains_list = len(pDomainList)
    for i, row in enumerate(pDomainList):
        # get intra / inter-tad data
        matrix, intertad_left, intertad_right = computeRegionsTADs(pMatrix, pDomainList, pCoolOrH5, i, row)

        if i - 1 > 0 and i + 1 < length_domains_list:
            density_inter_left = intertad_left.count_nonzero() / (intertad_left.shape[0] * intertad_left.shape[1])
            density_right_left = intertad_right.count_nonzero() / (intertad_right.shape[0] * intertad_right.shape[1])
        elif i - 1 <= 0 and i + 1 < length_domains_list:
            density_right_left = intertad_right.count_nonzero() / (intertad_right.shape[0] * intertad_right.shape[1])
            density_inter_left = -1 
        elif i - 1 > 0 and i + 1 >= length_domains_list:
            density_inter_left = intertad_left.count_nonzero() / (intertad_left.shape[0] * intertad_left.shape[1])
            density_right_left = -1

        density_intra = matrix.count_nonzero() / (matrix.shape[0] * matrix.shape[1])

        density_inter_left_list.append(density_inter_left)
        density_inter_right_list.append(density_right_left)
        density_intra_list.append(density_intra)

    pQueue.put([density_inter_left_list, density_inter_right_list, density_intra_list])
    return


def main(args=None):

    args = parse_arguments().parse_args(args)
    short_v_long_range = []
    for matrix in args.matrices:

        is_cooler = check_cooler(matrix)
        if not is_cooler:
            hic_matrix = hm.hiCMatrix(matrix)
            # hic_matrix.keepOnlyTheseChr([chromosome])
            # matrix = deepcopy(hic_matrix.matrix)
            # cut_intervals = deepcopy(hic_matrix.cut_intervals)
        else:
            hic_matrix = matrix
        if args.chromosomes is None:
            # get all chromosomes from cooler file
            if not is_cooler:
                chromosomes_list = list(hic_matrix.chrBinBoundaries)
            else:
                chromosomes_list = cooler.Cooler(matrix).chromnames
        else:
            chromosomes_list = args.chromosomes

        short_v_long_range_matrix_threads = [None] * args.threads
        chromosomesListPerThread = len(chromosomes_list) // args.threads
        all_data_collected = False
        queue = [None] * args.threads
        process = [None] * args.threads
        thread_done = [False] * args.threads
        for i in range(args.threads):

            if i < args.threads - 1:
                chromosomeListThread = chromosomes_list[i * chromosomesListPerThread:(i + 1) * chromosomesListPerThread]
            else:
                chromosomeListThread = chromosomes_list[i * chromosomesListPerThread:]

            queue[i] = Queue()
            # computeDensityTADs(pMatrix, pDomainList, pCoolOrH5, pQueue):
            process[i] = Process(target=computeDensityTADs, kwargs=dict(
                pMatrix=hic_matrix,
                pDomainList=chromosomeListThread,
                pCoolOrH5=is_cooler,
                pQueue=queue[i]
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(args.threads):
                if queue[i] is not None and not queue[i].empty():
                    short_v_long_range_matrix_threads[i] = queue[i].get()
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

        short_v_long_range_matrix = [item for sublist in short_v_long_range_matrix_threads for item in sublist]
        short_v_long_range.append(short_v_long_range_matrix)

    log.debug(short_v_long_range)
    plt.ylabel('Sum short range / long range')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    box_plot = plt.boxplot(short_v_long_range, patch_artist=True)
    legend_handels_color = []
    for i, patch in enumerate(box_plot['boxes']):
        patch.set_facecolor(args.colorList[i % len(args.colorList)])
        legend_handels_color.append(mpatches.Patch(color=args.colorList[i % len(args.colorList)], label=args.matrices[i]))
    # red_patch = mpatches.Patch(color='red', label='The red data')
    plt.legend(handles=legend_handels_color)
    # plt.legend(args.matrices)
    plt.savefig(args.plotFileName, dpi=args.dpi)

    if len(args.matrices) > 1:
        p_values = []
        for i, sample in enumerate(short_v_long_range):
            for sample2 in short_v_long_range[i + 1:]:
                statistic, significance_level = ranksums(sample, sample2)
                p_values.append(significance_level)
        log.debug('p_values {}'.format(p_values))
        with open(args.outFileName, 'w') as file:
            header = '# Created with HiCExplorer\'s hicPlotSVL ' + __version__ + '\n'
            header += "# Short range vs long range contacts per chromosome, p-values of each distribution against each other distribution with Wilcoxon rank-sum\n"
            header += '# Short range contacts: <= ' + str(args.distance) + '\n'
            file.write(header)
            counter = 0
            for i, matrix_0 in enumerate(args.matrices):
                for j, matrix_1 in enumerate(args.matrices[i + 1:]):
                    file.write(matrix_0 + '\t' + matrix_1 + '\t' + str(p_values[counter]) + '\n')
                    counter += 1
