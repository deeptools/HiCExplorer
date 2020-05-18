import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from hyperopt import hp, fmin, tpe, space_eval, STATUS_OK, Trials
import logging
log = logging.getLogger(__name__)
from tempfile import NamedTemporaryFile, mkdtemp

from hicexplorer import hicDetectLoops
from hicexplorer import hicValidateLocations
from hicexplorer._version import __version__

# matrixFile = ''
# proteinFile = ''
# proteinMaximum = 10000


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
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--resolution', '-r',
                           type=int,
                           default=10000,
                           help='Resolution of matrix')

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
    print(args)
    hicValidateLocations.main(args)
    data_dict = {}

    with open(outfile_statistics.name + '_statistics', 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            line_split = line.split(':')
            data_dict[line_split[0]] = float(line_split[1])

    if data_dict['Matched Loops']  > float(pMaximumNumberOfLoops):
        return 1 - ((data_dict['Loops match protein']*2 + 1.0) / 3) 
    if pMaximumNumberOfLoops > 500 and data_dict['Matched Loops'] < 500:
        return 1 - (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops))
    return 1 - ((data_dict['Loops match protein']*2 + (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops) / 2)) / 3)


def objective(pArgs):

    if pArgs['windowSize'] <= pArgs['peakWidth']:
        return 1
    outfile_loop = NamedTemporaryFile()
    args = "--matrix {} -o {} -pit {} -oet {} --windowSize {} --peakWidth {} -pp {} -p {} " \
        "--maxLoopDistance {}  -t 10 -tpc 10".format(
            pArgs['matrixFile'], outfile_loop.name,
            pArgs['pit'], pArgs['oet'], pArgs['windowSize'], pArgs['peakWidth'],pArgs['pp'], pArgs['p'],
            pArgs['maxLoopDistance']).split()
    hicDetectLoops.main(args)

    error_score = compute_score(outfile_loop.name, pArgs['proteinFile'], pArgs['maximumNumberOfLoops'], pArgs['resolution'])
    print('Error score: {}'.format(error_score))
    return error_score

def main(args=None):

    args = parse_arguments().parse_args(args)

    # print(compute_score(args.matrix, args.proteinFile, args.maximumNumberOfLoops, args.resolution))
    space = {

        'pit': hp.uniform('pit', 1, 100),
        'oet': hp.uniform('oet', 0, 5),
        'peakWidth': hp.choice('peakWidth', list(range(1, 10))),

        'windowSize': hp.choice('windowSize', list(range(4, 15))),
        'pp': hp.uniform('pp', 0.0000001, 0.15),
        'p': hp.uniform('p', 0.0000001, 0.1),
        # 'minLoopDistance': hp.choice('minLoopDistance', list(range(50000, 200000, 10000))),
        'maxLoopDistance': hp.choice('maxLoopDistance', list(range(1000000, 3000000, 100000))),
        # 'mip': hp.uniform('mip', 0.0, 0.5),
        # 'st': 'wilcoxon-rank-sum',#, 'anderson-darling']),
        # 'minVariance': hp.uniform('minVariance', 0.000001, 0.09),
        # 'maxVariance': hp.uniform('maxVariance', 0.1, 0.3),
        'matrixFile': args.matrix,
        'proteinFile': args.proteinFile,
        'maximumNumberOfLoops': args.maximumNumberOfLoops,
        'resolution': args.resolution

    }

    # minimize the objective over the space
    # best = fmin(objective, space, algo=tpe.suggest, max_evals=100)

    # print('best {}'.format(best))
    trials = Trials()
    # print(space_eval(space, best))
    best = fmin(objective, space, algo=tpe.suggest, max_evals=100, trials=trials)

    # print(best)
    print('best {}'.format(best))

    print('trails.result {}'.format(trials.results))
    print(space_eval(space, best))

    # -pit 10 --windowSize 8 -pp 0.1 -p 0.1 -o k562.bedgraph
    # --minLoopDistance 100000 -mip 0.1 -st wilcoxon-rank-sum
    # -t 20 -tpc 20 -
    # -minVariance 0.04 --maxVariance 0.4 -q 0.01 --maxLoopDistance 2000000
