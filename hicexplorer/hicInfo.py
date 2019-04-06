import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
from io import StringIO

import cooler

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicmatrix.HiCMatrix import check_cooler
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Prints information about a matrix or matrices including matrix size,
number of elements, sum of elements, etc.
An example usage is:
$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The matrix (or multiple matrices) to get information about. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                nargs='+',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--outFileName', '-o',
                           help='File name to save information of the matrix instead of writing it to the bash.'
                           )
    parserOpt.add_argument('--no_metadata', '-nm', action='store_false', help='Do not use meta data from cooler file to display information. '
                           'This method is slower and was default until version 2.2 of HiCExplorer. H5 files use always this parameter.')

    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    for matrix in args.matrices:
        # if
        generated_by = None
        genome_assembly = None
        statistics = None
        generated_by_cooler_lib = None
        tool_url = None
        matrix_generated_by = None
        matrix_generated_by_url = None
        creation_date = None
        chromosomes = None
        bin_length = None
        size = None
        nchroms = None
        num_non_zero = None
        min_non_zero = None
        max_non_zero = None
        sum_elements = None
        num_nan_bins = None

        if check_cooler(matrix) and args.no_metadata:
            cooler_file = cooler.Cooler(matrix)

            if cooler_file.info is not None:
                # log.debug('cooler_file.info {}'.format(cooler_file.info))
                if 'bin-size' in cooler_file.info:
                    bin_length = cooler_file.info['bin-size']
                if 'nbins' in cooler_file.info:
                    size = cooler_file.info['nbins']
                if 'nchroms' in cooler_file.info:
                    nchroms = cooler_file.info['nchroms']
                if 'chromosomes' in cooler_file.info:
                    chromosomes = cooler_file.info['chromosomes']
                if 'nnz' in cooler_file.info:
                    num_non_zero = cooler_file.info['nnz']
                if 'min-value' in cooler_file.info:
                    min_non_zero = cooler_file.info['min-value']
                if 'max-value' in cooler_file.info:
                    max_non_zero = cooler_file.info['max-value']
                if 'generated-by' in cooler_file.info:
                    generated_by = toString(cooler_file.info['generated-by'])
                if 'genome-assembly' in cooler_file.info:
                    genome_assembly = toString(
                        cooler_file.info['genome-assembly'])
                if 'metadata' in cooler_file.info:
                    if cooler_file.info['metadata'] is not None:
                        if 'statistics' in cooler_file.info['metadata']:
                            statistics = cooler_file.info['metadata']['statistics']
                if 'generated-by-cooler-lib' in cooler_file.info:
                    generated_by_cooler_lib = toString(
                        cooler_file.info['generated-by-cooler-lib'])
                if 'tool-url' in cooler_file.info:
                    tool_url = toString(cooler_file.info['tool-url'])
                if 'matrix-generated-by' in cooler_file.info:
                    matrix_generated_by = toString(
                        cooler_file.info['matrix-generated-by'])
                if 'matrix-generated-by-url' in cooler_file.info:
                    matrix_generated_by_url = toString(
                        cooler_file.info['matrix-generated-by-url'])
                if 'creation-date' in cooler_file.info:
                    creation_date = cooler_file.info['creation-date']
                if 'sum-elements' in cooler_file.info:
                    sum_elements = cooler_file.info['sum-elements']

        else:
            hic_ma = hm.hiCMatrix(matrix)
            size = hic_ma.matrix.shape[0]
            num_non_zero = hic_ma.matrix.nnz
            sum_elements = hic_ma.matrix.sum() / 2
            bin_length = hic_ma.getBinSize()
            num_nan_bins = len(hic_ma.nan_bins)
            min_non_zero = hic_ma.matrix.data.min()
            max_non_zero = hic_ma.matrix.data.max()

            chromosomes = list(hic_ma.chrBinBoundaries)

        information = StringIO()
        information.write(
            "# Matrix information file. Created with HiCExplorer's hicInfo version {}\n".format(__version__))

        if matrix is not None:
            information.write("File:\t{}\n".format(matrix))
        if creation_date is not None:
            information.write("Date:\t{}\n".format(creation_date))

        if genome_assembly is not None:
            information.write("Genome assembly:\t{}\n".format(genome_assembly))
        if size is not None:
            information.write("Size:\t{:,}\n".format(size))
        if bin_length is not None:
            information.write("Bin_length:\t{}\n".format(bin_length))
        if sum_elements is not None:
            information.write("Sum of matrix:\t{}\n".format(sum_elements))
        if chromosomes is not None:
            information.write("Chromosomes:\t{}\n".format(
                ", ".join(toString(chromosomes))))
        if nchroms is not None:
            information.write("Number of chromosomes:\t{}\n".format(nchroms))
        if num_non_zero is not None:
            information.write(
                "Non-zero elements:\t{:,}\n".format(num_non_zero))
        if min_non_zero is not None:
            information.write("Minimum (non zero):\t{}\n".format(min_non_zero))
        if max_non_zero is not None:
            information.write("Maximum:\t{}\n".format(max_non_zero))
        if num_nan_bins is not None:
            information.write("NaN bins:\t{}\n".format(num_nan_bins))

        if check_cooler(matrix):
            information.write('The following columns are available: {}\n'.format(
                cooler.Cooler(matrix).bins().columns.values))
        if generated_by is not None:
            information.write("\n\nGenerated by:\t{}\n".format(generated_by))

        if generated_by_cooler_lib is not None:
            information.write("Cooler library version:\t{}\n".format(
                generated_by_cooler_lib))
        if tool_url is not None:
            information.write("HiCMatrix url:\t{}\n".format(tool_url))
        if matrix_generated_by is not None:
            information.write(
                "Interaction matrix created with:\t{}\n".format(matrix_generated_by))
        if matrix_generated_by_url is not None:
            information.write("URL:\t{}\n".format(matrix_generated_by_url))

        if statistics is not None:
            information.write("\n\nBuild statistics:\n{}\n".format(statistics))

        if args.outFileName:
            with open(args.outFileName, 'w') as file:
                file.write(information.getvalue())
        else:
            print(information.getvalue())

        information.close()
