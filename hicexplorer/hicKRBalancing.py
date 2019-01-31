from KRBalancing import *
import numpy as np
import argparse

from hicmatrix import HiCMatrix as hm
from scipy import sparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
Computes a balanced matrix using the Knight and Ruiz balancing algorithm.

    $ hicKRBalancing --matrix hic_matrix.h5 -o balanced_matrix.h5

"""
    )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='HiCExplorer matrix in h5 format.',
                                required=True)

    parserRequired.add_argument('--outputFileName', '-o',
                                help='HiCExplorer matrix in h5 format.',
                                required=True)

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    ma = hm.hiCMatrix(args.matrix) #TODO slow for big matrices!
    ma.maskBins(ma.nan_bins)

    print("Load a float sparse matrix for balancing")
    d = kr_balancing(ma.matrix.astype(float))
    d.outer_loop()

    print("get out put")
    ma.matrix = d.get_output()


    ma.save(args.outputFileName, pSymmetric=False, pApplyCorrection=False)

#if __name__ == "__main__":
#    main()
