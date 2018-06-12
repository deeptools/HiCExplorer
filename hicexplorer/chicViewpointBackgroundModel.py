from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
import math
import logging
log = logging.getLogger(__name__)

import numpy as np

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
                    Build a background model for viewpoints.
                    An example usage is:
                    $ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
                    """
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The input matrices to build the background model on.',
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--viewpoints', '-b',
                                help='Bed file contains all viewpoints which should be used to build the background model.',
                                type=str,
                                required=True)
    parserRequired.add_argument('--range', '-r',
                                help='Region around the viewpoint to consider.',
                                type=int,
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser

def getViewpointValues(pHiCMatrix, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionList=None):

    hic = pHiCMatrix

    if len(pReferencePoint) == 2:
        view_point_start, view_point_end = pHiCMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
    elif len(pReferencePoint) == 3:
        view_point_start, view_point_end = pHiCMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
    else:
        log.error("No valid reference point given. {}".format(pReferencePoint))
        exit(1)

    view_point_range = pHiCMatrix.getRegionBinRange(pChromViewpoint, pRegion_start, pRegion_end)
    elements_of_viewpoint = view_point_range[1] - view_point_range[0]
    log.debug('elements_of_viewpoint: {}'.format(elements_of_viewpoint))
    data_list = np.zeros(elements_of_viewpoint)
    view_point_start_ = view_point_start
    interactions_list = None
    if pInteractionList is not None:
        interactions_list = []
    while view_point_start_ <= view_point_end:
        chrom, start, end, _ = hic.getBinPos(view_point_start_)
        for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
            data_list[j] += hic.matrix[view_point_start_, idx]
            if interactions_list is not None:
                chrom_second, start_second, end_second, _ = pHiCMatrix.getBinPos(idx)
                interactions_list.append((chrom, start, end, chrom_second, start_second, end_second, pHiCMatrix.matrix[view_point_start_, idx]))
        view_point_start_ += 1

    return [view_point_start, view_point_end, view_point_range, data_list, interactions_list]

def smooth_viewpoint():
    pass

def main():

    args = parse_arguments().parse_args()

    viewpoints = []
    interactions = []

    with open(args.viewpoints, 'r') as file:
        for line in file.readlines():
            # log.debug("line {}".format(line.split('\t')))
            line_ = line.split('\t')
            if len(line_) == 2:
                chrom, start, end = line_[0], line_[1], line_[1]
            else:
                chrom, start, end = line_
            viewpoints.append((chrom, start, end))
        # parser.read_file(file_h)
    
    # elements_of_viewpoint = args.range * 2 / hic_ma.bin_size
    relative_counts_conditions = []
    relative_positions = set()
    # log.debug('elements_of_viewpoint_max {}'.format(elements_of_viewpoint))
    bin_size = 0

    ####### - compute for each condition (matrix):
    #######  - all viewpoints and smooth them: sliding window approach
    #######  - after smoothing, sum all viewpoints up to one
    #######  - compute the percentage of each position with respect to the total interaction count
    ####### for models of all conditions:
    #######  - compute mean percentage for each bin
    #######  - compute SEM = Standard deviation / (square root of sample size)

    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        background_model_data = {}
        bin_size = hic_ma.getBinSize()
        for viewpoint in viewpoints:
            log.debug('len cutintervals{}'.format(len(hic_ma.cut_intervals)))
            log.debug('foo {}'.format(hic_ma.getChrBinRange(viewpoint[0])))
            log.info('getBinPos{}'.format(hic_ma.getBinPos(hic_ma.getChrBinRange(viewpoint[0])[1]-1)))
            max_length = hic_ma.getBinPos(hic_ma.getChrBinRange(viewpoint[0])[1]-1)[2]
            region_start = int(viewpoint[1]) - args.range
            if region_start < 0:
                region_start = 0
            log.info('region_start {}'.format(region_start))

            region_end = int(viewpoint[2]) + args.range
            if region_end > max_length:
                region_end = max_length
            log.info('region_end {}'.format(region_end))
            
            view_point_start, view_point_end, view_point_range, data_list, interaction_list = \
                getViewpointValues(hic_ma, viewpoint, viewpoint[0], region_start, region_end, True)

            # elements_ = view_point_range[1] - view_point_range[0]
            ### set data in relation to viewpoint, upstream are negative values, downstream positive, zero is viewpoint
            for i, data in zip(range(view_point_range[0], view_point_range[1], 1), data_list):
                relative_position = i - view_point_start
                if relative_position in background_model_data:
                    background_model_data[relative_position] += data
                else:
                    background_model_data[relative_position] = data
                    relative_positions.add(relative_position)
                
            
            log.debug('view_point_start {}'.format(view_point_start))
            log.debug('view_point_end {}'.format(view_point_end))
            log.debug('view_point_range {}'.format(view_point_range))
            log.debug('data_list {}'.format(data_list))
            log.debug('interaction_list {}'.format(len(interaction_list)))

        ## compute the percentage of each position with respect to the total interaction count
        total_count = 0.0
        for key, value in background_model_data.items():
            total_count += value
        
        relative_counts = background_model_data

        for key in relative_counts:
            relative_counts[key] /= total_count
        relative_counts_conditions.append(relative_counts)

    ###### for models of all conditions:
    #######  - compute mean percentage for each bin
    #######  - compute SEM = Standard deviation / (square root of sample size)
    mean_percentage = {}
    sem = {}
    relative_positions = sorted(relative_positions)
    for relative_position in relative_positions:
        i = 0
        count = 0
        for condition in relative_counts_conditions:
            if relative_position in condition:
                i += 1
                count += condition[relative_position]
        mean_percentage[relative_position] = count / i
        sem[relative_position] = (count / i) / math.sqrt(i)

    with open('background_model.bed', 'w') as file:
        for relative_position in relative_positions:
            relative_position_in_genomic_scale = relative_position * bin_size
            file.write("{}\t{:.12f}\t{:.12f}\n".format(relative_position_in_genomic_scale, mean_percentage[relative_position],
                                    sem[relative_position]))