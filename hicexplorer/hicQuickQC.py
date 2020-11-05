import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from tempfile import NamedTemporaryFile
import shutil
import argparse
import os
import errno

from hicexplorer import hicBuildMatrix as hicBuildMatrix
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description="""
The tool hicQuickQC considers the first n lines of two bam/sam files to get a first estimate of the quality of the data. It is highly recommended to set the restriction enzyme and dangling end parameter to get a good quality report.
""")
    parserRequired = parser.add_argument_group('Required arguments')
    parserRequired.add_argument('--samFiles', '-s',
                                help='The two PE alignment sam files to process.',
                                metavar='two sam files',
                                nargs=2,
                                required=True)

    parserRequired.add_argument('--QCfolder',
                                help='Path of folder to save the quality control data of the matrix. The log files '
                                'produced this way can be loaded into `hicQC` in order to compare the quality of multiple '
                                'Hi-C libraries.',
                                metavar='FOLDER',
                                required=True)
    parserRequired.add_argument('--restrictionCutFile', '-rs',
                                help=('BED file(s) with all restriction cut places '
                                      '(output of "findRestSite" command). '
                                      'Should contain only  mappable '
                                      'restriction sites. If given, the bins are '
                                      'set to match the restriction fragments (i.e. '
                                      'the region between one restriction site and '
                                      'the next). Alternatively, a fixed binSize can be defined instead. '
                                      'However, either binSize or restrictionCutFile must be defined. '
                                      'To use more than one restriction enzyme, generate for each one a restrictionCutFile and list them space seperated.'),
                                type=argparse.FileType('r'),
                                metavar='BED file',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--restrictionSequence', '-seq',
                                help='Sequence of the restriction site, if multiple are used, please list them space seperated. If a dangling sequence '
                                'is listed at the same time, please preserve the same order.',
                                type=str,
                                nargs='+',
                                required=True)

    parserRequired.add_argument('--danglingSequence',
                                help='Sequence left by the restriction enzyme after cutting, if multiple are used, please list them space '
                                'seperated and preserve the order. Each restriction enzyme '
                                'recognizes a different DNA sequence and, after cutting, they leave behind a specific '
                                '"sticky" end or dangling end sequence.  For example, for HindIII the restriction site '
                                'is AAGCTT and the dangling end is AGCT. For DpnII, the restriction site and dangling '
                                'end sequence are the same: GATC. This information is easily found on the description '
                                'of the restriction enzyme. The dangling sequence is used to classify and report reads '
                                'whose 5\' end starts with such sequence as dangling-end reads. A significant portion '
                                'of dangling-end reads in a sample are indicative of a problem with the re-ligation '
                                'step of the protocol.',
                                type=str,
                                nargs='+',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--lines',
                           help='Number of lines to consider for the QC test run'
                           ' (Default: %(default)s).',
                           required=False,
                           default=1000000,
                           type=int
                           )
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    if not os.path.exists(args.QCfolder):
        try:
            os.makedirs(args.QCfolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    outFile = NamedTemporaryFile(suffix='.h5', delete=False)
    args_hicBuildMatrix = "--samFiles {} {} --outFileName {}  --QCfolder {} --doTestRun --doTestRunLines {} --threads 1 ".format(args.samFiles[0], args.samFiles[1],
                                                                                                                                 outFile.name,
                                                                                                                                 args.QCfolder,
                                                                                                                                 str(args.lines)
                                                                                                                                 ).split()

    args_hicBuildMatrix.append('--binSize')
    args_hicBuildMatrix.append(str(10000))

    if args.restrictionSequence:
        args_hicBuildMatrix.append('--restrictionSequence')
        for restrictionSequence in args.restrictionSequence:
            args_hicBuildMatrix.append(restrictionSequence)

    if args.danglingSequence:
        args_hicBuildMatrix.append('--danglingSequence')
        for danglingSequence in args.danglingSequence:
            args_hicBuildMatrix.append(danglingSequence)

    if args.danglingSequence:
        args_hicBuildMatrix.append('--restrictionCutFile')
        for restrictionCutFile in args.restrictionCutFile:
            args_hicBuildMatrix.append(restrictionCutFile.name)

    log.debug('args_hicBuildMatrix {}'.format(args_hicBuildMatrix))

    hicBuildMatrix.main(args_hicBuildMatrix)
