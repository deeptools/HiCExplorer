import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=('Adds Hi-C matrices of the same size. Format '
                                                  'has to be hdf5 or npz'))

    parser.add_argument('--matrices', '-m',
                        help='matrices to add. Must have the same shape.',
                        metavar='.h5 file format',
                        nargs='+',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix. The output is '
                             'also a .h5 file. But don\'t add the suffix',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main():

    args = parse_arguments().parse_args()
    hic = hm.hiCMatrix(args.matrices[0])
    summed_matrix = hic.matrix
    nan_bins = set(hic.nan_bins)
    for matrix in args.matrices[1:]:
        hic_to_append = hm.hiCMatrix(matrix)
        if hic.chrBinBoundaries != hic_to_append.chrBinBoundaries:
            exit("The two matrices have different chromosome order. Use the tool `hicExport` to change the order.\n"
                 "{}: {}\n"
                 "{}: {}".format(args.matrices[0], hic.chrBinBoundaries.keys(),
                                 matrix, hic_to_append.chrBinBoundaries.keys()))

        try:
            summed_matrix = summed_matrix + hic_to_append.matrix
            if len(hic_to_append.nan_bins):
                nan_bins = nan_bins.union(hic_to_append.nan_bins)
        except:
            print "\nMatrix {} seems to be corrupted or of different " \
                  "shape".format(matrix)
            exit(1)

    # save only the upper triangle of the
    # symmetric matrix
    hic.setMatrixValues(summed_matrix)
    hic.maskBins(sorted(nan_bins))
    hic.save(args.outFileName)
