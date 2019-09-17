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

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--plotFileName', '-pfn',
                           help='Plot name.',
                           default='plot.png')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File the p-values are written to, p-values are only computed if at least two matrices are given.',
                           default='p_values.txt')
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


def compute_relation_short_long_range(pHiCMatrix, pChromosomes, pDistance, pIsCooler, pQueue):
    # hic_matrix = hm.hiCMatrix(
    #             pMatrixFile=args.matrix, pChrnameList=[args.region])
    svl_relations = []
    for chromosome in pChromosomes:
        if pIsCooler:
            hic_matrix_obj = hm.hiCMatrix(
                pMatrixFile=pHiCMatrix, pChrnameList=[chromosome])
            max_distance = pDistance / hic_matrix_obj.getBinSize()
            hic_matrix = hic_matrix_obj.matrix

        else:
            # in case it is h5 pHiCMatrix is an HiCMatrix object
            indices_chromosome = hic_matrix_target.getChrBinRange(chromosome)
            hic_matrix = pHiCMatrix.matrix[indices_chromosome[0]:indices_chromosome[1], indices_chromosome[0]:indices_chromosome[1]]
            max_distance = pDistance / pHiCMatrix.getBinSize()

        instances, features = hic_matrix.nonzero()
        distances = np.absolute(instances - features)
        mask = distances <= max_distance

        sum_smaller_max_distance = np.sum(hic_matrix.data[mask])
        sum_greater_max_distance = np.sum(hic_matrix.data[~mask])
        svl_relation = sum_smaller_max_distance / sum_greater_max_distance
        if np.isinf(svl_relation) or np.isnan(svl_relation):
            continue
        svl_relations.append(svl_relation)

    pQueue.put(svl_relations)
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

        if len(chromosomes_list) == 1:
            single_core = True
        else:
            single_core = False

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
            process[i] = Process(target=compute_relation_short_long_range, kwargs=dict(
                pHiCMatrix=hic_matrix,
                pChromosomes=chromosomeListThread,
                pDistance=args.distance,
                pIsCooler=is_cooler,
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
