from .lib import TADClassifier
import argparse
import logging

from hicexplorer._version import __version__

# taken and altered from hicFindTads
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve',
        description="""
Check out hicTADClassifier for calling TADs using our default classifiers or your trained classifier.

This program can be used to train and test new and existing classifiers for hicTADClassifier. These classifiers can later be run to call boundaries for TADs. There are four modes available: train_new, train_existing, train_test and predict_test. By default, an EasyEnsembleClassifier as described in Liu et al.: “Exploratory Undersampling for Class-Imbalance Learning” will be trained, but you can pass any sklearn classifier that allows for a warm start. You may also vary the resampling method and a range of hyperparameters to fine tune the model. Do mind to set the correct normalization method and resolution for the classifier. The program will check and raise warnings, when resolutions and normalization methods are mixed up. Also, a protein track file in the narrowPeak format with a threshold value may be passed to filter out low quality boundaries.


train_test mode: this is a convenience function, where a single matrix/domains set can be passed to quickly assert the performance of a new classifier. Nothing will be saved from this mode, instead, the classifier will be trained on 80% of the data and tested on the remaining 20%. The output will be a performance report. A quick usage example can be seen here:

$ hicTrainTADClassifier -m 'train_test' -f 'my_test_matrix.cool' -d 'domains.bed' -o 'report.txt' -n 'range' -r 10000

train_new mode: this mode allows the training of a new classifier. Note that range of optional arguments, that can be used to fine tune. The resulting classifier will be pickled at the specified out_file. A quick example can be seen here, where we varied the feature distance:

$ hicTrainTADClassifier -m 'train_new' -f 'my_test_matrix.cool' -d 'domains.bed' -o 'new_classifier.data' -n 'range' -r 10000 --distance 18

train_existing mode: train the classifier specified in saved_classifier on new data. When not setting the saved_classifier, the preset classifiers will be used as preset. The output will be the classifier trained with additional data and more internal estimators.

$ hicTrainTADClassifier -m 'train_existing' -f 'my_test_matrix.cool' -d 'domains.bed' -o 'updated_classifier.data' -n 'range' -r 10000

predict_test mode: predict using an existing classifier and produce a classification report. The difference in using this over hicTADClassifier is, that this version will predict on a balanced test set. Normally, HiC-Matrices contain a lot more non-boundaries than boundaries, which skews the classification report to the point, where it does not contain usefull information anymore. By passing a domain file produced by another TAD Caller, hicTrainClassifier will build a test set using the boundaries of this domain file and will pick at random as many non-boundaries from the passed matrix. Use this over hicTADClassifier to produce a meaningful output, but not for TAD calling.

$ hicTrainTADClassifier -m 'predict_test' -f 'my_test_matrix.cool' -d 'domains.bed' -o 'report.txt' -n 'range' -r 10000
        """)

    def check_leniency(le):
        li = int(le)
        if li < 0:
            raise argparse.ArgumentTypeError(
                "leniency %s is not within valid range 0 -" % le)
        return li

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--mode', '-mo',
                                choices=[
                                    'train_new',
                                    'train_existing',
                                    'train_test',
                                    'predict_test'],
                                help='choice of program',
                                required=True,
                                type=str)

    parserRequired.add_argument('--matrices', '-m',
                                help='HiC-Matrix file or list of files for input. Only COOLER files are supported!',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--domain_file', '-d',
                                help='domain file or list of files containing tad boundaries',
                                required=True,
                                nargs='+')

    parserRequired.add_argument('--out_file', '-o',
                                help='output file for either the classification report or the saved classifier',
                                required=True,
                                type=str)

    parserRequired.add_argument('--normalization_method', '-n',
                                help='set the normalization mode, with which the passed matrices will be normalized',
                                type=str,
                                choices=[
                                    'obs_exp',
                                    'range'
                                ])

    parserRequired.add_argument('--resolution', '-r',
                                help='resolution in bases of the classifier',
                                type=int)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--threshold',
                           help='threshold for protein quality check',
                           type=float,
                           default=None)

    parserOpt.add_argument('--leniency',
                           help='leniency for protein quality check. Widens peaks of protein file by leniency*resolution',
                           type=check_leniency,
                           default=0)

    parserOpt.add_argument('--saved_classifier',
                           help='pickled classifier to be trained or used for prediction',
                           type=str,
                           default=None)

    parserOpt.add_argument('--unselect_border_cases',
                           help='set whether genes at the border of the matrix up to set distance will not be used for training and testing',
                           action='store_false')

    parserOpt.add_argument('--protein_file',
                           help='provide a bed file for TAD quality control',
                           nargs='+',
                           default=None)

    parserOpt.add_argument('--threads', '-t',
                           help='number of threads used',
                           default=4,
                           type=int)
    parserOpt.add_argument('--chromosomes',
                           help='Chromosomes to include in the analysis. If not set, all chromosomes are included.',
                           nargs='+')
    parserOpt.add_argument('--concatenate_before_resample',
                           help='whether features build from matrix list are concatenated and resampled together or resampled separatly per matrix. Not important for random undersampling, but alter for other resampling methods and check if performance increases.',
                           default=False,
                           action='store_true')

    parserOpt.add_argument('--estimators_per_step',
                           help='how many estimators are added in each training step for the classifier (new classifier)',
                           #    choices=range(5, 1001),
                           metavar="[5-1000]",
                           default=20,
                           type=int)

    parserOpt.add_argument('--resampling_method',
                           help='the method used to resample the training set(new classifier)',
                           type=str,
                           choices=[
                               'undersample_cluster_centroids',
                               'undersample_random',
                               'passed_method'
                           ],
                           default='undersample_random')

    parserOpt.add_argument('--alternative_resampling_method',
                           help='pass alternative resampling method from imblearn library (new classifier)',
                           default=None)

    parserOpt.add_argument('--distance',
                           help='max distance between TADs to be used in calculation (new classifier)',
                           #    choices=range(5, 31),
                           metavar="[5-30]",
                           default=15,
                           type=int)

    parserOpt.add_argument('--impute_value',
                           help='non-numerical float values in matrix will be replaced by this value (new classifier)',
                           type=float,
                           default=-1.0)

    parserOpt.add_argument('--alternative_classifier',
                           help='pass custom classifier, needs to implement warm_start (new classifier)',
                           default=None)

    parserOpt.add_argument('--use_cleanlab',
                           help='use Confident Learning with the cleanlab library (new classifier)',
                           action='store_true')
    # parserOpt.add_argument('--chrPrefixProtein', '-cp',
    #                        help='Adding / removing / do nothing a \'chr\'-prefix to chromosome name of the protein.',
    #                        choices=[None, 'add', 'remove'],
    #                        default=None
    #                        )
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser
# def print_args(args):
#     """
#     Print to stderr the parameters used
#     Parameters
#     ----------
#     args

#     Returns
#     -------

#     """
#     for key, value in args._get_kwargs():
#         log.info("{}:\t{}\n".format(key, value))


def main(args=None):

    args = parse_arguments().parse_args(args)
    program = TADClassifier(args.mode,
                            args.out_file,
                            saved_classifier=args.saved_classifier,
                            normalization_method=args.normalization_method,
                            unselect_border_cases=args.unselect_border_cases,
                            threshold=args.threshold,
                            leniency=args.leniency,
                            resolution=args.resolution,
                            threads=args.threads,
                            concatenate_before_resample=args.concatenate_before_resample,
                            distance=args.distance,
                            impute_value=args.impute_value,
                            resampling_method=args.resampling_method,
                            alternative_resampling_method=args.alternative_resampling_method,
                            alternative_classifier=args.alternative_classifier,
                            use_cleanlab=args.use_cleanlab,
                            estimators_per_step=args.estimators_per_step,
                            pAddRemoveChrPrexix=None
                            )
    program.run_hicTrainClassifier(
        args.matrices,
        args.domain_file,
        args.protein_file,
        args.chromosomes)
