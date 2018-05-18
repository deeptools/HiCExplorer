from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString

import logging
log = logging.getLogger(__name__)


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



def main():

    args = parse_arguments().parse_args()

    viewpoints = []
    interactions = []

    with open(args.viewpoints, 'r') as file:
        for line in file.readlines():
            chrom, start, end = line.split('\t')
            viewpoints.append((chrom, start, end))
        # parser.read_file(file_h)
    background_model_data = []
    
    ####### - compute for each condition (matrix):
    #######  - all viewpoints and smooth them: sliding window approach
    #######  - after smoothing, sum all viewpoints up to one
    #######  - compute the percentage of each position with respect to the total interaction count
    ####### for models of all conditions:
    #######  - compute mean percentage for each bin
    #######  - compute SEM = Standard deviation / (square root of sample size)
    for matrix in args.matrices:
        hic_ma = hm.hiCMatrix(matrix)
        
        
        for viewpoint in viewpoints:
            max_length = hic_ma.getBinPos(hic_ma.getChrBinRange(viewpoint[0])[2])
            region_start = int(viewpoint[1]) - args.range
            if int(viewpoint[1]) - args.range < 0:
                continue
            region_end = int(viewpoint[2]) - args.range
            if int(viewpoint[2]) + args.range >  max_length:
                continue
            view_point_start, view_point_end, view_point_range, data_list = \
                getViewpointValues(hic_ma, viewpoint, viewpoint[0], region_start, region_end)

