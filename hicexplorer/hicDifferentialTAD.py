# get tad border (domains) file
# get matrix 1 and matrix 2
# get square out of matrix 1 and matrix 2 as defined by border files
# apply statistical test with H0 that they are the same.

import pandas as pd
import numpy as np
from scipy.stats import ranksums
import argparse

import logging
log = logging.getLogger(__name__)
from hicmatrix import HiCMatrix as hm

from hicexplorer.utilities import check_cooler, getRegion
from hicexplorer._version import __version__

def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Computes differential TADs by comparing the precomputed TAD regions of the target-matrix with the same regions of the control matrix. 
Please notice that the matrices need to have the same read coverage, this can be achieved with hicNormalize and the 'smallest'-mode. """)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--targetMatrix', '-tm',
                                help='The matrix which was used to compute the TADs',
                                required=False)
    parserRequired.add_argument('--controlMatrix', '-cm',
                                help='The control matrix to test the TADs for a differential interaction pattern.',
                                required=False)
    parserRequired.add_argument('--tadDomains', '-td',
                                help='The TADs domain file computed by hicFindTADs.',
                                required=False)
    parserRequired.add_argument('--outFileNamePrefix', '-o',
                                help='Outfile name prefix to store the accepted / rejected H0 TADs.',
                                required=False)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--pValue', '-p',
                           type=float,
                           help='Only candidates with p-values less the given threshold will be accepted as differential.',
                           default=0.05)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser

def readDomainBoundaries(pFile):
    domains_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2]]

    return domains_df

def main(args=None):
    args = parse_arguments().parse_args(args)

    # read domains file
    domains_df = readDomainBoundaries(args.tadDomains)
    # read full h5 or only region if cooler
    is_cooler_target = check_cooler(args.targetMatrix)
    is_cooler_control = check_cooler(args.controlMatrix)

    if is_cooler_target != is_cooler_control:
        log.error('Matrices are not given in the same format!')
        exit(1)
    if not is_cooler_control:
        hic_matrix_target = hm.hiCMatrix(args.targetMatrix)
        hic_matrix_control = hm.hiCMatrix(args.controlMatrix)

    accepted_H0 = []
    rejected_H0 = []
    log.debug('domains_df {}'.format(domains_df))
    for row in domains_df.values.tolist():
    # for domain in domains_df:
        log.debug('domain {} {} {}'.format(row[0], row[1], row[2]))
        if is_cooler_target:
            hic_matrix_target = hm.hiCMatrix(
                pMatrixFile=args.targetMatrix, pChrnameList=[str(row[0]) + ':' + str(row[1])+'-'+str(row[2])])
            hic_matrix_control = hm.hiCMatrix(
                pMatrixFile=args.controlMatrix, pChrnameList=[str(row[0]) + ':' + str(row[1])+'-'+str(row[2])])
            matrix_target = hic_matrix_target.matrix.toarray()
            matrix_control = hic_matrix_control.matrix.toarray()
        else:
            chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(args, hic_matrix_target)
            matrix_target = np.asarray(hic_matrix_target.matrix[idx1, :][:, idx2].todense().astype(float))
            chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(args, hic_matrix_control)
            matrix_control = np.asarray(hic_matrix_control.matrix[idx1, :][:, idx2].todense().astype(float))

        matrix_target = matrix_target.flatten()
        matrix_control = matrix_control.flatten()

        statistic, significance_level = ranksums(sorted(matrix_target), sorted(matrix_control))
        if significance_level <= args.pValue:
            rejected_H0.append(row)
            rejected_H0[-1].append(significance_level)
        else:
            accepted_H0.append(row)
            accepted_H0[-1].append(significance_level)

    with open(args.outFileNamePrefix + '_accepted.bed', 'w') as file:
        for row in accepted_H0:
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\n')
    with open(args.outFileNamePrefix + '_rejected.bed', 'w') as file:
        for row in rejected_H0:
            row_list = list(map(str, row))
            file.write('\t'.join(row_list))
            file.write('\n')

