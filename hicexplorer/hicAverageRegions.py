"""hicAverageRegions: a thin layer over the C++ executable hicAverageRegions.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
import logging
log = logging.getLogger(__name__)
import numpy as np
from scipy.sparse import csr_matrix, save_npz, lil_matrix, coo_matrix


def parse_arguments(args=None):
    """The argparse parser of hicAverageRegions, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAverageRegions')


def calculateViewpointRange(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    '''
    This function computes the correct start and end position of a viewpoint given the reference and the range.
    '''
    start_out_of_range = False
    end_out_of_range = False
    max_length = pHiCMatrix.getBinPos(pHiCMatrix.getChrBinRange(pViewpoint[0])[1] - 1)[2]
    # bin_size = pHiCMatrix.getBinSize()
    # _range = [pRange[0], pRange[1]]

    if pCoordinatesToBinMapping == 'start':
        region_start = int(pViewpoint[1]) - pRange[0]
        region_end = int(pViewpoint[1]) + pRange[1]
    elif pCoordinatesToBinMapping == 'end':
        region_start = int(pViewpoint[2]) - pRange[0]
        region_end = int(pViewpoint[2]) + pRange[1]
    elif pCoordinatesToBinMapping == 'center':
        viewpoint_center_value = int(float(pViewpoint[1]) + ((float(pViewpoint[2]) - float(pViewpoint[1])) / 2))
        region_start = viewpoint_center_value - pRange[0]
        region_end = viewpoint_center_value + pRange[1]

    if region_start < 0:
        region_start = 0
        start_out_of_range = True

    if region_end > max_length:
        # -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
        region_end = max_length - 1
        end_out_of_range = True
    return region_start, region_end, start_out_of_range, end_out_of_range


def getBinIndices(pHiCMatrix, pViewpoint):

    return pHiCMatrix.getRegionBinRange(pViewpoint[0], pViewpoint[1], pViewpoint[2])


def calculateViewpointRangeBins(pHiCMatrix, pViewpoint, pRange, pCoordinatesToBinMapping):
    # if pCoordinatesToBinMapping == 'start_end':
    #     viewpoint_index_start = getBinIndices(pHiCMatrix, pViewpoint)[0]
    #     viewpoint_index_end = getBinIndices(pHiCMatrix, pViewpoint)[1]
    start_out_of_range = False
    end_out_of_range = False
    if pCoordinatesToBinMapping == 'start':
        viewpoint_index = getBinIndices(pHiCMatrix, pViewpoint)[0]
    elif pCoordinatesToBinMapping == 'end':
        viewpoint_index = getBinIndices(pHiCMatrix, pViewpoint)[1]
    else:
        viewpoint_center_value = int(float(pViewpoint[1]) + ((float(pViewpoint[2]) - float(pViewpoint[1])) / 2))
        viewpoint_center = [pViewpoint[0], viewpoint_center_value, viewpoint_center_value]
        viewpoint_index = getBinIndices(pHiCMatrix, viewpoint_center)[1]

    # if pCoordinatesToBinMapping == 'start_end':
    #     start = viewpoint_index_start - pRange[0]
    #     end = viewpoint_index_end + pRange[1]
    # else:
    first_bin, last_bin = pHiCMatrix.getChrBinRange(pViewpoint[0])
    start = viewpoint_index - pRange[0]
    end = viewpoint_index + pRange[1]
    if start < first_bin:
        start = first_bin
        start_out_of_range = True
    if end > last_bin:
        end = last_bin
        end_out_of_range = True
    return start, end, start_out_of_range, end_out_of_range


def main(args=None):
    """Run the C++ implementation of hicAverageRegions with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAverageRegions')(args)
