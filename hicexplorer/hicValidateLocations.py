"""hicValidateLocations: a thin layer over the C++ executable hicValidateLocations.

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
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)

import pandas as pd
from pybedtools import BedTool
import numpy as np
import cooler


def parse_arguments(args=None):
    """The argparse parser of hicValidateLocations, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicValidateLocations')


def readThreeColumn(pFile, pChrPrefixProtein):
    protein_df = pd.read_csv(pFile, sep='\t', header=None)[[0, 1, 2]]
    log.debug('protein_df {}'.format(protein_df))
    if pChrPrefixProtein == 'add':
        protein_df[0] = 'chr' + protein_df[0].astype(str)
    elif pChrPrefixProtein == 'remove':
        protein_df[0] = protein_df[0].str.lstrip('chr')
    return protein_df


def readLoopFile(pInputFile, pChrPrefixLoops):
    full_loop = pd.read_csv(pInputFile, sep='\t', header=None)
    if pChrPrefixLoops == 'add':
        full_loop[0] = 'chr' + full_loop[0].astype(str)
        full_loop[3] = 'chr' + full_loop[3].astype(str)
    elif pChrPrefixLoops == 'remove':
        full_loop[0] = full_loop[0].str.lstrip('chr')
        full_loop[3] = full_loop[3].str.lstrip('chr')

    return full_loop


def overlapLoop(pDataFrameLoop, pDataFrameProtein):
    loop_bedtool_x = BedTool.from_dataframe(pDataFrameLoop[[0, 1, 2]])
    loop_bedtool_y = BedTool.from_dataframe(pDataFrameLoop[[3, 4, 5]])

    protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
    x = loop_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()
    y = loop_bedtool_y.intersect(protein_bedtool, c=True).to_dataframe()

    mask_x = x['name'] >= 1
    mask_y = y['name'] >= 1

    selection = (mask_x) & (mask_y)

    return selection


def overlapTAD(pDataFrameTAD, pDataFrameProtein):
    tad_bedtool_x = BedTool.from_dataframe(pDataFrameTAD)

    log.debug('tad_bedtool_x {}'.format(tad_bedtool_x))

    protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
    x = tad_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()

    mask_x = x['name'] >= 1
    return pd.Series(mask_x)


def correlateCool(pCoolFile, pDataFrameLoop):
    cool_obj = cooler.Cooler(pCoolFile)
    number_of_peaks = cool_obj.info['nnz']
    match_list = []
    for line in pDataFrameLoop.itertuples():
        try:
            matrix_data = cool_obj.matrix(balance=False, sparse=True).fetch(str(line[1]) + ':' + str(line[2]) + '-' + str(line[3]), str(line[4]) + ':' + str(line[5]) + '-' + str(line[6]))

            if matrix_data.nnz > 0:
                match_list.append(True)
            else:
                match_list.append(False)
        except Exception:
            match_list.append(False)

    selection = pd.Series(match_list)
    log.debug(selection)

    return number_of_peaks, selection


# def correlate2DText(pDataFrameProtein, pDataFrameLoop):
#     loop_bedtool_x = BedTool.from_dataframe(pDataFrameLoop[[0, 1, 2, 3, 4, 5]])

#     protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
#     x = loop_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()

#     mask_x = x['name'] >= 1

#     log.debug(mask_x)
#     return mask_x


def applyBinning(pDataFrame, pBinSize, pMerge=True):
    pDataFrame_out = pDataFrame.copy()
    pDataFrame_out[1] = (pDataFrame[1] / pBinSize).astype(int) * pBinSize
    pDataFrame_out[2] = ((pDataFrame[2] / pBinSize).astype(int) + 1) * pBinSize

    log.debug('pDataFrame_out {}'.format(pDataFrame_out))
    pDataFrame_out = pDataFrame_out.drop_duplicates()
    log.debug('pDataFrame_out {}'.format(pDataFrame_out))

    if pMerge:
        bedtools_data = BedTool.from_dataframe(pDataFrame_out)
        bedtools_data = bedtools_data.merge()
        log.debug('bedtools_data {}'.format(bedtools_data.to_dataframe()))

        bedtools_data = bedtools_data.sort()
        log.debug('bedtools_data {}'.format(bedtools_data.to_dataframe()))

        return bedtools_data.to_dataframe()
    else:
        return pDataFrame_out


def writeLoopFile(pOutFileName, pLoopDataFrame):
    pLoopDataFrame.to_csv(pOutFileName, sep='\t', header=False, index=False)


def main(args=None):
    """Run the C++ implementation of hicValidateLocations with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicValidateLocations')(args)
