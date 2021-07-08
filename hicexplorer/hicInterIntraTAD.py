import pandas as pd
import numpy as np
from scipy.stats import ranksums
import argparse
from multiprocessing import Process, Queue
import time
import traceback
from copy import deepcopy
import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Extracts and computes different inter and intra TAD values and ratios.""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix which was used to compute the TADs',
                                required=False)
    parserRequired.add_argument('--tadDomains', '-td',
                                help='The TADs domain file computed by hicFindTADs.',
                                required=False)
    parserRequired.add_argument('--outFileName', '-o',
                                help='Outfile name',
                                required=False,
                                default='output_interintra_tad.tzt'
                                )
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads to use, the parallelization is implemented per chromosome'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2, 3, 4, 5]]

    return domains_df


def computeInterIntraTADs(pMatrix, pDomainList, pCoolOrH5, pThreadId, pQueue):
    try:

        inter_left_sum_list = []
        inter_right_sum_list = []
        inter_left_densit_list = []
        inter_right_density_list = []
        inter_left_number_of_contacts_list = []
        inter_right_number_of_contacts_list = []
        inter_left_number_of_contacts_nnz_list = []
        inter_right_number_of_contacts_nzz_list = []

        intra_sum_list = []
        intra_number_of_contacts_list = []
        intra_number_of_contacts_nnz_list = []
        intra_density_list = []
        inter_left_intra_ratio_list = []
        inter_right_intra_ratio_list = []
        inter_left_inter_right_intra_ratio_list = []

        rows = []

        chromosome_list = pDomainList
        for i, row in enumerate(chromosome_list):

            if pThreadId is None:
                log.debug('first thread')
                if i == len(chromosome_list) - 1:
                    continue
            elif pThreadId == True:
                log.debug('middle thread')

                if i == 0 or i == len(chromosome_list) - 1:
                    log.debug('i: {}'.format(i))
                    log.debug('len(chromosome_list): {}'.format(len(chromosome_list)))

                    continue
            elif pThreadId == False:
                log.debug('last thread')

                if i == 0:
                    continue

            if i - 1 >= 0:
                chromosom = chromosome_list[i - 1][0]
                start = chromosome_list[i - 1][1]
            else:
                chromosom = chromosome_list[i][0]
                start = chromosome_list[i][1]
            if i + 1 < len(chromosome_list):
                end = chromosome_list[i + 1][2]
            else:
                end = chromosome_list[i][2]
            # midpos = row[1] + ((row[2] - row[1]) / 2)

            if pCoolOrH5:

                # # get intra-TAD data
                hic_matrix = hm.hiCMatrix(
                    pMatrixFile=pMatrix, pChrnameList=[str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])])
                matrix = hic_matrix.matrix

                hic_matrix_inter_tad = hm.hiCMatrix(
                    pMatrixFile=pMatrix, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])

                matrix_inter_tad = hic_matrix_inter_tad.matrix

            else:
                # in case of h5 pMatrixTarget is already a HiCMatrix object
                hic_matrix = pMatrix
                hic_matrix_inter_tad = pMatrix
                indices = hic_matrix.getRegionBinRange(str(row[0]), row[1], row[2])

                matrix = hic_matrix.matrix[indices[0]:indices[1], indices[0]:indices[1]]
                matrix_inter_tad = pMatrix.matrix

            # matrix = matrix.flatten()

            # get index position left tad with tad
            left_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])[0]
            if pCoolOrH5:
                outer_left_boundary_index = 0

                outer_right_boundary_index = -1

            else:
                outer_left_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

                outer_right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]

            if i + 1 < len(chromosome_list) and not pCoolOrH5:
                # get index position right tad with tad
                right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
            elif i + 1 < len(chromosome_list):
                right_boundary_index = hic_matrix_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]

            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                intertad_left = matrix_inter_tad[outer_left_boundary_index:left_boundary_index, left_boundary_index:right_boundary_index]
                intertad_right = matrix_inter_tad[left_boundary_index:right_boundary_index, right_boundary_index:outer_right_boundary_index]

            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                intertad_right = matrix_inter_tad[left_boundary_index:right_boundary_index, right_boundary_index:outer_right_boundary_index]

            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                intertad_left = matrix_inter_tad[outer_left_boundary_index:left_boundary_index, left_boundary_index:right_boundary_index]

            inter_left_sum = 0
            inter_right_sum = 0
            inter_left_density = 0
            inter_right_density = 0
            inter_left_number_of_contacts = 0
            inter_right_number_of_contacts = 0
            inter_left_number_of_contacts_nnz = 0
            inter_right_number_of_contacts_nzz = 0

            intra_sum = matrix.sum()
            intra_number_of_contacts = matrix.shape[0] * matrix.shape[1]
            intra_number_of_contacts_nnz = matrix.nnz
            intra_density = intra_number_of_contacts_nnz / intra_number_of_contacts_nnz
            # both inter, left and right is available
            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                # intertad_left = intertad_left.flatten()
                # intertad_right = intertad_right.flatten()
                inter_left_sum = intertad_left.sum()
                inter_right_sum = intertad_right.sum()

                inter_left_number_of_contacts = intertad_left.shape[0] * intertad_left.shape[1]
                inter_right_number_of_contacts = intertad_right.shape[0] * intertad_right.shape[1]
                inter_left_number_of_contacts_nnz = intertad_left.nnz
                inter_right_number_of_contacts_nzz = intertad_right.nnz

                inter_left_density = inter_left_number_of_contacts_nnz / inter_left_number_of_contacts_nnz
                inter_right_density = inter_right_number_of_contacts_nzz / inter_right_number_of_contacts
                # statistic_left, significance_level_left = ranksums(intertad_left, intertad_left_control)
                # statistic_right, significance_level_right = ranksums(intertad_right, intertad_right_control)
            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                # inter right is available
                # intertad_right = intertad_right.flatten()
                inter_right_sum = intertad_right.sum()
                inter_right_number_of_contacts = intertad_right.shape[0] * intertad_right.shape[1]
                inter_right_number_of_contacts_nzz = intertad_right.nnz
                inter_right_density = inter_right_number_of_contacts_nzz / inter_right_number_of_contacts

                # statistic_right, significance_level_right = ranksums(intertad_right, intertad_right_control)
            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                # inter left is available

                # intertad_left = intertad_left.flatten()
                inter_left_sum = intertad_left.sum()
                inter_left_number_of_contacts = intertad_left.shape[0] * intertad_left.shape[1]
                inter_left_number_of_contacts_nnz = intertad_left.nnz
                inter_left_density = inter_left_number_of_contacts_nnz / inter_left_number_of_contacts_nnz

                # statistic_left, significance_level_left = ranksums(intertad_left, intertad_left_control)

            inter_left_intra_ratio = inter_left_sum / intra_sum
            inter_right_intra_ratio = inter_right_sum / intra_sum
            inter_left_inter_right_intra_ratio = (inter_left_sum + inter_right_sum) / intra_sum

            inter_left_sum_list.append(inter_left_sum)
            inter_right_sum_list.append(inter_right_sum)
            inter_left_densit_list.append(inter_left_density)
            inter_right_density_list.append(inter_right_density)
            inter_left_number_of_contacts_list.append(inter_left_number_of_contacts)
            inter_right_number_of_contacts_list.append(inter_right_number_of_contacts)
            inter_left_number_of_contacts_nnz_list.append(inter_left_number_of_contacts_nnz)
            inter_right_number_of_contacts_nzz_list.append(inter_right_number_of_contacts_nzz)

            intra_sum_list.append(intra_sum)
            intra_number_of_contacts_list.append(intra_number_of_contacts)
            intra_number_of_contacts_nnz_list.append(intra_number_of_contacts_nnz)
            intra_density_list.append(intra_density)
            inter_left_intra_ratio_list.append(inter_left_intra_ratio)
            inter_right_intra_ratio_list.append(inter_right_intra_ratio)
            inter_left_inter_right_intra_ratio_list.append(inter_left_inter_right_intra_ratio)

            rows.append(row)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    pQueue.put([inter_left_sum_list,
                inter_right_sum_list,
                inter_left_densit_list,
                inter_right_density_list,
                inter_left_number_of_contacts_list,
                inter_right_number_of_contacts_list,
                inter_left_number_of_contacts_nnz_list,
                inter_right_number_of_contacts_nzz_list,
                intra_sum_list,
                intra_number_of_contacts_list,
                intra_number_of_contacts_nnz_list,
                intra_density_list,
                inter_left_intra_ratio_list,
                inter_right_intra_ratio_list,
                inter_left_inter_right_intra_ratio_list,
                rows])


def main(args=None):
    args = parse_arguments().parse_args(args)

    # read domains file
    domains_df = readDomainBoundaries(args.tadDomains)
    log.debug('len(domains_df) {}'.format(len(domains_df)))
    domains = domains_df.values.tolist()
    old_chromosome = None

    tads_per_chromosome = []

    for j in range(len(domains)):
        if old_chromosome is None:
            old_chromosome = domains[j][0]
            per_chromosome = []
            per_chromosome.append(domains[j])

        elif old_chromosome == domains[j][0]:
            per_chromosome.append(domains[j])
            continue
        else:
            tads_per_chromosome.append(per_chromosome)
            per_chromosome = []
            per_chromosome.append(domains[j])
            old_chromosome = domains[j][0]
    tads_per_chromosome.append(per_chromosome)

    # read full h5 or only region if cooler
    is_cooler = check_cooler(args.matrix)

    if not is_cooler:
        hic_matrix = hm.hiCMatrix(args.matrix)
    else:
        hic_matrix = args.matrix

    inter_left_sum_list_chromosomes = []
    inter_right_sum_list_chromosomes = []
    inter_left_density_list_chromosomes = []
    inter_right_density_list_chromosomes = []
    inter_left_number_of_contacts_list_chromosomes = []
    inter_right_number_of_contacts_list_chromosomes = []
    inter_left_number_of_contacts_nnz_list_chromosomes = []
    inter_right_number_of_contacts_nzz_list_chromosomes = []

    intra_sum_list_chromosomes = []
    intra_number_of_contacts_list_chromosomes = []
    intra_number_of_contacts_nnz_list_chromosomes = []
    intra_density_list_chromosomes = []
    inter_left_intra_ratio_list_chromosomes = []
    inter_right_intra_ratio_list_chromosomes = []
    inter_left_inter_right_intra_ratio_list_chromosomes = []

    rows_chromosomes = []



    inter_left_sum_list_threads = [[]] * args.threads
    inter_right_sum_list_threads = [[]] * args.threads
    inter_left_density_list_threads = [[]] * args.threads
    inter_right_density_list_threads = [[]] * args.threads
    inter_left_number_of_contacts_list_threads = [[]] * args.threads
    inter_right_number_of_contacts_list_threads = [[]] * args.threads
    inter_left_number_of_contacts_nnz_list_threads = [[]] * args.threads
    inter_right_number_of_contacts_nzz_list_threads = [[]] * args.threads

    intra_sum_list_threads = [[]] * args.threads
    intra_number_of_contacts_list_threads = [[]] * args.threads
    intra_number_of_contacts_nnz_list_threads = [[]] * args.threads
    intra_density_list_threads = [[]] * args.threads
    inter_left_intra_ratio_list_threads = [[]] * args.threads
    inter_right_intra_ratio_list_threads = [[]] * args.threads
    inter_left_inter_right_intra_ratio_list_threads = [[]] * args.threads

    rows_threads = [[]] * args.threads

    threads_save = deepcopy(args.threads)
    for chromosome in tads_per_chromosome:
        log.debug('tads_per_chromosome {}'.format(chromosome))
        domainsPerThread = len(chromosome) // args.threads
        if domainsPerThread == 0 and len(chromosome) > 0:
            domainsPerThread = 1
            args.threads = 1
        elif domainsPerThread > 0:
            args.threads = threads_save

        all_data_collected = False
        queue = [None] * args.threads
        process = [None] * args.threads
        thread_done = [False] * args.threads
        # None --> first thread, process first element in list, ignore last one
        # True --> middle thread: ignore first and last element in tad processing
        # False --> last thread: ignore first element, process last one
        thread_id = None
        for i in range(args.threads):

            if args.threads == 1:
                domainListThread = chromosome

            elif i == 0:
                domainListThread = chromosome[i * domainsPerThread:((i + 1) * domainsPerThread) + 1]
                thread_id = None
            elif i < args.threads - 1:
                domainListThread = chromosome[(i * domainsPerThread) - 1:((i + 1) * domainsPerThread) + 1]
                thread_id = True

            else:
                domainListThread = chromosome[(i * domainsPerThread) - 1:]
                thread_id = False

            if args.threads == 1:
                thread_id = ''

            log.debug('len(domainListThread) {}'.format(len(domainListThread)))
            log.debug('len(thread_id) {}'.format(thread_id))

            queue[i] = Queue()
            process[i] = Process(target=computeInterIntraTADs, kwargs=dict(
                pMatrix=hic_matrix,
                # pMatrixControl=hic_matrix_control,
                pDomainList=domainListThread,
                pCoolOrH5=is_cooler,
                # pPValue=args.pValue,
                pThreadId=thread_id,
                pQueue=queue[i]
            )
            )

            process[i].start()
        fail_flag = False
        fail_message = ''
        while not all_data_collected:
            for i in range(args.threads):

                if queue[i] is not None and not queue[i].empty():
                    queue_data = queue[i].get()
                    if 'Fail:' in queue_data:
                        fail_flag = True
                        fail_message = queue_data
                    else:
                        inter_left_sum_list_threads[i], \
                            inter_right_sum_list_threads[i], \
                            inter_left_density_list_threads[i], \
                            inter_right_density_list_threads[i], \
                            inter_left_number_of_contacts_list_threads[i], \
                            inter_right_number_of_contacts_list_threads[i], \
                            inter_left_number_of_contacts_nnz_list_threads[i], \
                            inter_right_number_of_contacts_nzz_list_threads[i], \
                            intra_sum_list_threads[i], \
                            intra_number_of_contacts_list_threads[i], \
                            intra_number_of_contacts_nnz_list_threads[i], \
                            intra_density_list_threads[i], \
                            inter_left_intra_ratio_list_threads[i], \
                            inter_right_intra_ratio_list_threads[i], \
                            inter_left_inter_right_intra_ratio_list_threads[i], \
                            rows_threads[i] = queue_data

                    queue[i] = None
                    process[i].join()
                    process[i].terminate()
                    process[i] = None
                    thread_done[i] = True
                # elif queue[i] is None and

            all_data_collected = True
            for thread in thread_done:
                if not thread:
                    all_data_collected = False
            time.sleep(1)


        if fail_flag:
            log.error(fail_message[6:])
            exit(1)

        inter_left_sum_list_chromosomes.append([item for sublist in inter_left_sum_list_threads for item in sublist])
        inter_right_sum_list_chromosomes.append([item for sublist in inter_right_sum_list_threads for item in sublist])
        inter_left_density_list_chromosomes.append([item for sublist in inter_left_density_list_threads for item in sublist])
        inter_right_density_list_chromosomes.append([item for sublist in inter_right_density_list_threads for item in sublist])
        inter_left_number_of_contacts_list_chromosomes.append([item for sublist in inter_left_number_of_contacts_list_threads for item in sublist])
        inter_right_number_of_contacts_list_chromosomes.append([item for sublist in inter_right_number_of_contacts_list_threads for item in sublist])
        inter_left_number_of_contacts_nnz_list_chromosomes.append([item for sublist in inter_left_number_of_contacts_nnz_list_threads for item in sublist])
        inter_right_number_of_contacts_nzz_list_chromosomes.append([item for sublist in inter_right_number_of_contacts_nzz_list_threads for item in sublist])

        intra_sum_list_chromosomes.append([item for sublist in intra_sum_list_threads for item in sublist])
        intra_number_of_contacts_list_chromosomes.append([item for sublist in intra_number_of_contacts_list_threads for item in sublist])
        intra_number_of_contacts_nnz_list_chromosomes.append([item for sublist in intra_number_of_contacts_nnz_list_threads for item in sublist])
        intra_density_list_chromosomes.append([item for sublist in intra_density_list_threads for item in sublist])
        inter_left_intra_ratio_list_chromosomes.append([item for sublist in inter_left_intra_ratio_list_threads for item in sublist])
        inter_right_intra_ratio_list_chromosomes.append([item for sublist in inter_right_intra_ratio_list_threads for item in sublist])
        inter_left_inter_right_intra_ratio_list_chromosomes.append([item for sublist in inter_left_inter_right_intra_ratio_list_threads for item in sublist])

        rows_chromosomes.append([item for sublist in rows_threads for item in sublist])


    inter_left_sum_list = [item for sublist in inter_left_sum_list_chromosomes for item in sublist]
    inter_right_sum_list = [item for sublist in inter_right_sum_list_chromosomes for item in sublist]
    inter_left_density_list = [item for sublist in inter_left_density_list_chromosomes for item in sublist]
    inter_right_density_list = [item for sublist in inter_right_density_list_chromosomes for item in sublist]
    inter_left_number_of_contacts_list = [item for sublist in inter_left_number_of_contacts_list_chromosomes for item in sublist]
    inter_right_number_of_contacts_list = [item for sublist in inter_right_number_of_contacts_list_chromosomes for item in sublist]
    inter_left_number_of_contacts_nnz_list = [item for sublist in inter_left_number_of_contacts_nnz_list_chromosomes for item in sublist]
    inter_right_number_of_contacts_nzz_list = [item for sublist in inter_right_number_of_contacts_nzz_list_chromosomes for item in sublist]

    intra_sum_list = [item for sublist in intra_sum_list_chromosomes for item in sublist]
    intra_number_of_contacts_list = [item for sublist in intra_number_of_contacts_list_chromosomes for item in sublist]
    intra_number_of_contacts_nnz_list = [item for sublist in intra_number_of_contacts_nnz_list_chromosomes for item in sublist]
    intra_density_list = [item for sublist in intra_density_list_chromosomes for item in sublist]
    inter_left_intra_ratio_list = [item for sublist in inter_left_intra_ratio_list_chromosomes for item in sublist]
    inter_right_intra_ratio_list = [item for sublist in inter_right_intra_ratio_list_chromosomes for item in sublist]
    inter_left_inter_right_intra_ratio_list = [item for sublist in inter_left_inter_right_intra_ratio_list_chromosomes for item in sublist]

    rows = [item for sublist in rows_chromosomes for item in sublist]


    with open(args.outFileName, 'w') as file:
        header = '# Created with HiCExplorer\'s hicInterIntraTAD version ' + __version__ + '\n'
        header += '# Chromosome\tstart\tend\tname\tscore\tstrand\tinter_left_sum\tinter_right_sum\tinter_left_density\tinter_right_density\tinter_left_number_of_contacts\tinter_right_number_of_contacts\t'  \
            'inter_left_number_of_contacts_nnz\tinter_right_number_of_contacts_nnz\tintra_sum\tintra_number_of_contacts\tintra_number_of_contacts_nnz\tintra_density\tinter_left_intra_ratio\tinter_right_intra_ratio\tinter_left_inter_right_intra_ratio\n'
        file.write(header)
        for i, row in enumerate(rows):
            row_list = list(map(str, row))

            file.write('\t'.join(row_list))

            file.write('\t{}'.format(inter_left_sum_list[i]))
            file.write('\t{}'.format(inter_right_sum_list[i]))
            file.write('\t{}'.format(inter_left_density_list[i]))
            file.write('\t{}'.format(inter_right_density_list[i]))
            file.write('\t{}'.format(inter_left_number_of_contacts_list[i]))
            file.write('\t{}'.format(inter_right_number_of_contacts_list[i]))
            file.write('\t{}'.format(inter_left_number_of_contacts_nnz_list[i]))
            file.write('\t{}'.format(inter_right_number_of_contacts_nzz_list[i]))
            file.write('\t{}'.format(intra_sum_list[i]))
            file.write('\t{}'.format(intra_number_of_contacts_list[i]))
            file.write('\t{}'.format(intra_number_of_contacts_nnz_list[i]))
            file.write('\t{}'.format(intra_density_list[i]))
            file.write('\t{}'.format(inter_left_intra_ratio_list[i]))
            file.write('\t{}'.format(inter_right_intra_ratio_list[i]))
            file.write('\t{}'.format(inter_left_inter_right_intra_ratio_list[i]))

            file.write('\n')
