from KRBalancing import *
import numpy as np
import argparse

from hicmatrix import HiCMatrix as hm
from scipy import sparse

import logging
log = logging.getLogger(__name__)

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
    log.info('Reading matrix')

    ma = hm.hiCMatrix(args.matrix) #TODO slow for big matrices!
    log.debug('Reading matrix... Done')

    ma.maskBins(ma.nan_bins)
    log.debug('masking matrix... Done')

    #TODO remove zero rows and columns
    log.debug("Load a float sparse matrix for balancing")
    d = kr_balancing(ma.matrix.astype(float))
    d.outer_loop()

    print("get out put")
    ma.matrix = d.get_output()
    log.debug("save matrix")

    ma.save(args.outputFileName, pApplyCorrection=False)

#if __name__ == "__main__":
#    main()
