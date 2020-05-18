import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from hyperopt import hp, fmin, tpe, space_eval, STATUS_OK, Trials
import logging
log = logging.getLogger(__name__)
from tempfile import NamedTemporaryFile, mkdtemp

# from hicexplorer import hicDetectLoops
from hicexplorer import hicValidateLocations
from hicexplorer._version import __version__

import subprocess
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
    parserRequired.add_argument('--juicerPath', '-j',
                                help='path to juicer.jar',
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
    if data_dict['Matched Loops'] < 500:
        return 1 - (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops))
    return 1 - ((data_dict['Loops match protein']*2 + (data_dict['Matched Loops'] / float(pMaximumNumberOfLoops) / 2)) / 3)


def objective(pArgs):
    print('objective start')
    if pArgs['i'] <= pArgs['p']:
        return 1
    # outfile_loop = NamedTemporaryFile()
    output_folder = mkdtemp(prefix="output_")
    print('juicer')
    bashCommand = "java -jar {} hiccups --restrict --threads 4 -k KR -p {} -i {} -f {} -r {} {} {}".format(pArgs['juicer'], pArgs['p'], pArgs['i'], 
                                                                            pArgs['f'], pArgs['resolution'],
                                                                            pArgs['matrixFile'], output_folder
                                                                            )
    subprocess.check_output(['bash','-c', bashCommand])
    print('juicer done')
    # log.info('{}'.format(output))
    bashCommand = 'bedtools sort -i {}/merged_loops.bedpe > {}/merged_sorted.bedpe'.format(output_folder, output_folder)
    subprocess.check_output(['bash','-c', bashCommand])


    error_score = compute_score(output_folder+'/merged_sorted.bedpe', pArgs['proteinFile'], pArgs['maximumNumberOfLoops'], pArgs['resolution'])
    print('Error score: {}'.format(error_score))
    return error_score

def main(args=None):
    print('foo')
    args = parse_arguments().parse_args(args)

    # print(compute_score(args.matrix, args.proteinFile, args.maximumNumberOfLoops, args.resolution))
    space = {
# [-p peak width] [-i window] [-t thresholds] 
# 		[-d centroid distances]
        'p': hp.choice('t', list(range(1, 10))),
        'i': hp.choice('i', list(range(2, 15))),
        'f': hp.uniform('f', 0.0000001, 0.5),
        
        'matrixFile': args.matrix,
        'proteinFile': args.proteinFile,
        'maximumNumberOfLoops': args.maximumNumberOfLoops,
        'resolution': args.resolution,
        'juicer' : args.juicerPath

    }

    # minimize the objective over the space
    # best = fmin(objective, space, algo=tpe.suggest, max_evals=100)

    # print('best {}'.format(best))
    trials = Trials()
    print("goo")
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
if __name__ == "__main__":
    main()

