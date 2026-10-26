"""hicDifferentialTAD: a thin layer over the C++ executable hicDifferentialTAD.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
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
    """The argparse parser of hicDifferentialTAD, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicDifferentialTAD')


def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2, 3, 4, 5]]

    return domains_df


def computeDifferentialTADs(pMatrixTarget, pMatrixControl, pDomainList, pCoolOrH5, pPValue, pThreadId, pQueue):
    try:
        accepted_inter_left = []
        accepted_inter_right = []
        accepted_intra = []
        p_values_list = []
        stats_list = []
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
                hic_matrix_target = hm.hiCMatrix(
                    pMatrixFile=pMatrixTarget, pChrnameList=[str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])])
                hic_matrix_control = hm.hiCMatrix(
                    pMatrixFile=pMatrixControl, pChrnameList=[str(row[0]) + ':' + str(row[1]) + '-' + str(row[2])])
                matrix_target = hic_matrix_target.matrix.toarray()
                matrix_control = hic_matrix_control.matrix.toarray()

                hic_matrix_target_inter_tad = hm.hiCMatrix(
                    pMatrixFile=pMatrixTarget, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])
                hic_matrix_control_inter_tad = hm.hiCMatrix(
                    pMatrixFile=pMatrixControl, pChrnameList=[str(chromosom) + ':' + str(start) + '-' + str(end)])

                matrix_target_inter_tad = hic_matrix_target_inter_tad.matrix
                matrix_control_inter_tad = hic_matrix_control_inter_tad.matrix

            else:
                # in case of h5 pMatrixTarget is already a HiCMatrix object
                hic_matrix_target = pMatrixTarget
                hic_matrix_control = pMatrixControl
                hic_matrix_target_inter_tad = pMatrixTarget
                hic_matrix_control_inter_tad = pMatrixControl
                indices_target = hic_matrix_target.getRegionBinRange(str(row[0]), row[1], row[2])
                indices_control = hic_matrix_control.getRegionBinRange(str(row[0]), row[1], row[2])

                matrix_target = hic_matrix_target.matrix[indices_target[0]:indices_target[1], indices_target[0]:indices_target[1]].toarray()
                matrix_control = hic_matrix_control.matrix[indices_control[0]:indices_control[1], indices_control[0]:indices_control[1]].toarray()
                matrix_target_inter_tad = pMatrixTarget.matrix
                matrix_control_inter_tad = pMatrixControl.matrix

            matrix_target = matrix_target.flatten()
            matrix_control = matrix_control.flatten()
            # tad_midpoint = hic_matrix_target_inter_tad.getRegionBinRange(str(row[0]), midpos, midpos)[0]

            # if i - 1 >= 0:
            # get index position left tad with tad
            left_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])[0]
            left_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[1], row[1])[0]
            if pCoolOrH5:
                outer_left_boundary_index_target = 0
                outer_left_boundary_index_control = 0

                outer_right_boundary_index_control = -1
                outer_right_boundary_index_target = -1

            else:
                outer_left_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]
                outer_left_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[0]

                outer_right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]
                outer_right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), start, end)[1]

            if i + 1 < len(chromosome_list) and not pCoolOrH5:
                # get index position right tad with tad
                right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
                right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
            elif i + 1 < len(chromosome_list):
                right_boundary_index_target = hic_matrix_target_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]
                right_boundary_index_control = hic_matrix_control_inter_tad.getRegionBinRange(str(chromosom), row[2], row[2])[0]

            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
                intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
                intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()
                intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                intertad_right_target = matrix_target_inter_tad[left_boundary_index_target:right_boundary_index_target, right_boundary_index_target:outer_right_boundary_index_target].toarray()
                intertad_right_control = matrix_control_inter_tad[left_boundary_index_control:right_boundary_index_control, right_boundary_index_control:outer_right_boundary_index_control].toarray()

            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                intertad_left_target = matrix_target_inter_tad[outer_left_boundary_index_target:left_boundary_index_target, left_boundary_index_target:right_boundary_index_target].toarray()
                intertad_left_control = matrix_control_inter_tad[outer_left_boundary_index_control:left_boundary_index_control, left_boundary_index_control:right_boundary_index_control].toarray()

            significance_level_left = None
            significance_level_right = None
            statistic_left = None
            statistic_right = None

            if i - 1 >= 0 and i + 1 < len(chromosome_list):
                intertad_left_target = intertad_left_target.flatten()
                intertad_left_control = intertad_left_control.flatten()
                intertad_right_target = intertad_right_target.flatten()
                intertad_right_control = intertad_right_control.flatten()

                statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)
                statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
            elif i - 1 < 0 and i + 1 < len(chromosome_list):
                intertad_right_target = intertad_right_target.flatten()
                intertad_right_control = intertad_right_control.flatten()
                statistic_right, significance_level_right = ranksums(intertad_right_target, intertad_right_control)
            elif i - 1 > 0 and i + 1 >= len(chromosome_list):
                intertad_left_target = intertad_left_target.flatten()
                intertad_left_control = intertad_left_control.flatten()
                # log.debug('intertad_left_target {}'.format(intertad_left_target))
                # log.debug('intertad_left_control {}'.format(intertad_left_control))

                statistic_left, significance_level_left = ranksums(intertad_left_target, intertad_left_control)

            # log.debug('matrix_target {}'.format(matrix_target))
            # log.debug('matrix_control {}'.format(matrix_control))

            statistic, significance_level = ranksums(matrix_target, matrix_control)
            # log.debug('statistic {}, significance_level {}'.format(statistic, significance_level))
            # log.debug('right statistic {}, significance_level {}'.format(statistic_right, significance_level_right))
            # log.debug('left statistic {}, significance_level {}'.format(statistic_left, significance_level_left))

            p_values = []
            stats = []
            if significance_level_left is None or np.isnan(significance_level_left):
                accepted_inter_left.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level_left <= pPValue:
                accepted_inter_left.append(1)
                p_values.append(significance_level_left)
                stats.append(statistic_left)
            else:
                accepted_inter_left.append(0)
                p_values.append(significance_level_left)
                stats.append(statistic_left)

            if significance_level_right is None or np.isnan(significance_level_right):
                accepted_inter_right.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level_right <= pPValue:
                accepted_inter_right.append(1)
                p_values.append(significance_level_right)
                stats.append(statistic_right)
            else:
                accepted_inter_right.append(0)
                p_values.append(significance_level_right)
                stats.append(statistic_right)

            if significance_level is None or np.isnan(significance_level):
                accepted_intra.append(0)
                p_values.append(np.nan)
                stats.append(np.nan)
            elif significance_level <= pPValue:
                accepted_intra.append(1)
                p_values.append(significance_level)
                stats.append(statistic)
            else:
                accepted_intra.append(0)
                p_values.append(significance_level)
                stats.append(statistic)

            p_values_list.append(p_values)
            stats_list.append(stats)

            rows.append(row)
    except Exception as exp:
        pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
        return
    # hic_matrix_target_inter_tad.save('manipulated_target.cool')
    # hic_matrix_control_inter_tad.save('manipulated_control.cool')
    pQueue.put([stats_list, p_values_list, accepted_inter_left, accepted_inter_right, accepted_intra, rows])


def main(args=None):
    """Run the C++ implementation of hicDifferentialTAD with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicDifferentialTAD')(args)
