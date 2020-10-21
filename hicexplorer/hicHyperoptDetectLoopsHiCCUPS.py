import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from hyperopt import hp, fmin, tpe, space_eval, STATUS_OK, Trials
import logging
log = logging.getLogger(__name__)
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import hicValidateLocations
from hicexplorer._version import __version__

import subprocess


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to compute the loops on.',
                                required=True)
    parserRequired.add_argument('--proteinFile', '-p',
                                help='The protein file to validate the detected loops',
                                required=True)
    parserRequired.add_argument('--maximumNumberOfLoops', '-ml',
                                help='The maximum number of loops that should be used for optimization computation.',
                                required=True,
                                type=int)
    parserRequired.add_argument('--juicerPath', '-j',
                                help='path to juicer.jar',
                                required=True)
    parserRequired.add_argument('--outputFileName', '-o',
                                help='File names for the result of the optimization'
                                ' (Default: %(default)s).',
                                default='hyperoptHiCCUPS_result.txt',
                                required=False)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--resolution', '-r',
                           type=int,
                           default=10000,
                           help='Resolution of matrix'
                           ' (Default: %(default)s).')
    parserOpt.add_argument('--runs',
                           type=int,
                           default=100,
                           help='Number of runs of hyperopt.')
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--normalization', '-k',
                           help='Normalization table name',
                           required=True)
    parserOpt.add_argument('--cpu',
                           help='use the CPU version',
                           action='store_true'
                           )
    parserOpt.add_argument('--restricted',
                           help='If the GPU version is used, search only within 8 MB.',
                           action='store_true')

    parserOpt.add_argument('--help', '-h', action='help',
                           help='Show this help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def compute_score(pLoopFile, pProteinFile, pMaximumNumberOfLoops, pResolution):
    with open(pLoopFile, 'r') as file:
        lines = file.readlines()
        if len(lines) == 0:
            return 1
    outfile_statistics = NamedTemporaryFile()
    args = "--data {} --protein {} -cl --resolution {} --outFileName {}".format(pLoopFile, pProteinFile, pResolution, outfile_statistics.name).split()
    hicValidateLocations.main(args)
    data_dict = {}

    with open(outfile_statistics.name + '_statistics', 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            line_split = line.split(':')
            data_dict[line_split[0]] = float(line_split[1])

    if data_dict['Matched Loops'] > float(pMaximumNumberOfLoops):
        return 1 - ((data_dict['Loops match protein'] * 2 + 1.0) / 3)
    if data_dict['Matched Loops'] < 500:
        return 1 - (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops))
    return 1 - ((data_dict['Loops match protein'] * 2 + (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops) / 2)) / 3)


def objective(pArgs):
    if pArgs['i'] <= pArgs['p']:
        return 1
    output_folder = mkdtemp(prefix="output_")
    bashCommand = "java -jar {} hiccups {} {} --threads {} -k KR -p {} -i {} -f {} -r {} {} {}".format(pArgs['juicer'], pArgs['restricted'],
                                                                                                       pArgs['cpu'], pArgs['threads'], pArgs['p'], pArgs['i'],
                                                                                                       pArgs['f'], pArgs['resolution'],
                                                                                                       pArgs['matrixFile'], output_folder
                                                                                                       )
    subprocess.check_output(['bash', '-c', bashCommand])
    bashCommand = 'bedtools sort -i {}/merged_loops.bedpe > {}/merged_sorted.bedpe'.format(output_folder, output_folder)
    subprocess.check_output(['bash', '-c', bashCommand])
    error_score = compute_score(output_folder + '/merged_sorted.bedpe', pArgs['proteinFile'], pArgs['maximumNumberOfLoops'], pArgs['resolution'])
    print('Error score: {}'.format(error_score))
    return error_score


def main(args=None):
    args = parse_arguments().parse_args(args)

    space = {
        'p': hp.choice('t', list(range(1, 10))),
        'i': hp.choice('i', list(range(2, 15))),
        'f': hp.uniform('f', 0.0000001, 0.5),
        'matrixFile': args.matrix,
        'proteinFile': args.proteinFile,
        'maximumNumberOfLoops': args.maximumNumberOfLoops,
        'resolution': args.resolution,
        'juicer': args.juicerPath,
        'restricted': '--restricted' if args.restricted else '',
        'cpu': '--cpu' if args.cpu else '',
        'threads': args.threads

    }

    # minimize the objective over the space
    trials = Trials()
    best = fmin(objective, space, algo=tpe.suggest, max_evals=args.runs, trials=trials)

    with open(args.outputFileName, 'w') as file:
        file.write("# Created by HiCExplorer hicHyperoptDetectLoopsHiCCUPS {}\n\n".format(__version__))
        file.write("{}".format(space_eval(space, best)))
