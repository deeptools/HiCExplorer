"""hicInterIntraTAD: a thin layer over the C++ executable hicInterIntraTAD.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import pandas as pd
import argparse
from multiprocessing import Process, Queue
import time
import traceback
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm

from hicexplorer.utilities import check_cooler
from hicexplorer._version import __version__


def parse_arguments(args=None):
    """The argparse parser of hicInterIntraTAD, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicInterIntraTAD')


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
            intra_density = intra_number_of_contacts_nnz / intra_number_of_contacts
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

                inter_left_density = inter_left_number_of_contacts_nnz / inter_left_number_of_contacts
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
                inter_left_density = inter_left_number_of_contacts_nnz / inter_left_number_of_contacts

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
    """Run the C++ implementation of hicInterIntraTAD with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicInterIntraTAD')(args)
