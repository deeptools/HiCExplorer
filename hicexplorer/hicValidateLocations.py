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


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
 
""")


    # TAD or Loop data track
    # protein data, compute correlation
    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--data', '-d',
                                help='The loop or TAD domains file',
                                required=True)
    parserRequired.add_argument('--method', '-m', # loop or domain
                                help='The loop or TAD domains file',
                                required=True)
    parserRequired.add_argument('--resolution', '-r', # loop or domain
                                help='The used resolution of the Hi-C interaction matrix.',
                                required=True,
                                type=int)
    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--outFileName', '-o',
                           help=''
                           )
    parserOpt.add_argument('--removeLocations', '-rl',
                           help='Remove locations with no protein overlap.'
                           )
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser

def readProtein(pFile):
    return pd.read_csv(pFile, sep='\t', header=(-1))[[0, 1, 2]]

def readLoopFile(pInputFile, pAddChr):
    full_loop= pd.read_csv(pInputFile, sep='\t', header=(-1))
    if pAddChr:
        full_loop[0] = 'chr' + full_loop[0].astype(str)
        full_loop[3] = 'chr' + full_loop[3].astype(str)
    return full_loop

def overlapLoop(pDataFrameLoop, pDataFrameProtein):
    loop_bedtool_x = BedTool.from_dataframe(pDataFrameLoop[[0,1,2]])
    loop_bedtool_y = BedTool.from_dataframe(pDataFrameLoop[[3,4,5]])

    protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
    x = loop_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()
    y = loop_bedtool_y.intersect(protein_bedtool, c=True).to_dataframe()
    
    mask_x = x['name'] >= 1
    mask_y = y['name'] >= 1
    
    selection  = (mask_x) & (mask_y)
    
    return selection

def overlapTAD(pDataFrameTAD, pDataFrameProtein):
    loop_bedtool_x = BedTool.from_dataframe(pDataFrameTAD[[0,1,2]])
    # loop_bedtool_y = BedTool.from_dataframe(pDataFrameLoop[[3,4,5]])

    protein_bedtool = BedTool.from_dataframe(pDataFrameProtein)
    x = loop_bedtool_x.intersect(protein_bedtool, c=True).to_dataframe()
    # y = loop_bedtool_y.intersect(protein_bedtool, c=True).to_dataframe()
    
    mask_x = x['name'] >= 1
    # mask_y = y['name'] >= 1
    
    # selection  = (mask_x) & (mask_y)
    
    return mask_x

def applyBinning(pDataFrame, pBinSize):
    pDataFrame_out = pDataFrame.copy()
    pDataFrame_out[1] = (pDataFrame[1] / pBinSize).astype(int) * pBinSize
    pDataFrame_out[2] = ((pDataFrame[2] / pBinSize).astype(int)+1) * pBinSize
    pDataFrame_out.drop_duplicates()
    bedtools_data = BedTool.from_dataframe(pDataFrame_out)
    bedtools_data = bedtools_data.merge()
    bedtools_data = bedtools_data.sort()
    return bedtools_data.to_dataframe()



def main(args=None):

    args = parse_arguments().parse_args(args)

    if args.method == 'loops':
        loop_df = readLoopFile(args.data, True)
        protein_df = applyBinning(ctcf_df, args.resolution)

        overlap_mask_df = overlapLoop(loop_df, protein_df)

        if args.removeLocations:
            loop_df = loop_df[overlap_mask_df]

    elif args.method == 'tads':
        tad_df = readLoopFile(args.data, True)
        protein_df = applyBinning(tad_df, args.resolution)

        overlap_mask_df = overlapLoop(tad_df, protein_df)

        if args.removeLocations:
            tad_df = loop_df[overlap_mask_df]

        