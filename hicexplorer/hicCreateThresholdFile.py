import argparse
from hicexplorer._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
                                     """
                                     )

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--thresholdValue', '-tv',
                                help='Standard threshold value for all relative distances. Please manipulate this value manual later to achieve a different threshold per distance where it is necessary.',
                                type=float,
                                required=True)
    parserRequired.add_argument('--range',
                                help='Defines the region upstream and downstream of a reference point which should be included. Smallest upstream value is 0, all positive upstream values are interpreted as negative.'
                                'Format is --region -upstream downstream, e.g. --region -500000 500000 plots 500kb up- and 500kb downstream.',
                                required=True,
                                type=int,
                                nargs=2)
    parserRequired.add_argument('--resolution', '-r',
                                help='Resolution of the bin in genomic units. Values are set as number of bases, e.g. 1000 for a 1kb, 5000 for a 5kb or 10000 for a 10kb resolution.'
                                'This value is used to merge neighboring bins'
                                ' (Default: %(default)s).',
                                type=int,
                                default=1000,
                                required=False)

    parserRequired.add_argument('--outFileName', '-o',
                                help='The name and path of the created threshold file.',
                                required=True)

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)

    with open(args.outFileName, 'w') as file:
        header = '# Threshold file of HiCExplorer\'s hicCreateThresholdFile version '
        header += str(__version__)
        header += '\n'
        header += '# Standard threshold {}\n'.format(args.thresholdValue)
        file.write(header)

        if args.range[0] > 0:
            args.range[0] = -args.range[0]

        for i in range(args.range[0], args.range[1] + args.resolution, args.resolution):
            file.write('{}\t{}\n'.format(i, args.thresholdValue))
