import argparse
import os

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


# own tools
from hicexplorer.utilities import genomicRegion
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)

from hicexplorer.lib.buildMatrixMethods import *


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description=('Using an alignment from a program that supports '
                     'local alignment (eg. Bowtie2) where both '
                     'PE reads are mapped using  the --local '
                     'option, this program reads such file and '
                     'creates a matrix of interactions.'
                     ))

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--samFiles', '-s',
                                help='The two PE alignment sam files to process',
                                metavar='two sam files',
                                nargs=2,
                                type=argparse.FileType('r'),
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='Output file name for the Hi-C matrix.',
                                metavar='FILENAME',
                                type=argparse.FileType('w'),
                                required=True)

    parserRequired.add_argument('--QCfolder',
                                help='Path of folder to save the quality control data for the matrix. The log files '
                                'produced this way can be loaded into `hicQC` in order to compare the quality of multiple '
                                'Hi-C libraries.',
                                metavar='FOLDER',
                                required=True)
    parserRequired.add_argument('--restrictionCutFile', '-rs',
                                help='BED file(s) with all restriction cut sites '
                                '(output of "hicFindRestSite" command). '
                                'Should only contain the restriction sites of the same genome which has been used '
                                'to generate the input sam files. Using regions of a different genome version can '
                                'generate false results! To use more than one restriction enzyme, generate '
                                'a restrictionCutFile for each enzyne and list them space seperated.',
                                type=argparse.FileType('r'),
                                metavar='BED file',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--restrictionSequence', '-seq',
                                help='Sequence of the restriction site, if multiple are used, '
                                'please list them space seperated. If a dangling sequence '
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

    parserOpt.add_argument('--outBam', '-b',
                           help='Output bam file to process. Optional parameter. '
                           'A bam file containing all valid Hi-C reads can be created '
                           'using this option. This bam file could be useful to inspect '
                           'the distribution of valid Hi-C reads pairs or for other '
                           'downstream analyses, but is not used by any HiCExplorer tool. '
                           'Computation will be significantly longer if this option is set.',
                           metavar='bam file',
                           type=argparse.FileType('w'),
                           required=False)

    # group = parserOpt.add_mutually_exclusive_group(required=True)

    parserOpt.add_argument('--binSize', '-bs',
                           help='Size in bp for the bins. The bin size depends '
                           'on the depth of sequencing. Use a larger bin size for '
                           'libraries sequenced with lower depth. If not given, matrices of restriction site resolution will be built. '
                           'Optionally for mcool file format: Define multiple resolutions which are all a multiple of the first value. '
                           ' Example: --binSize 10000 20000 50000 will create a mcool file formate containing the three defined resolutions.',
                           type=int,
                           nargs='+')

    parserOpt.add_argument('--minDistance',
                           help='Minimum distance between restriction sites. '
                           'Restriction sites that are closer than this '
                           'distance are merged into one. This option only '
                           'applies if --restrictionCutFile is given'
                           ' (Default: %(default)s).',
                           type=int,
                           default=300,
                           required=False)

    parserOpt.add_argument('--maxDistance',
                           help='This parameter is now obsolete. Use --maxLibraryInsertSize instead',
                           type=int)

    parserOpt.add_argument('--maxLibraryInsertSize',
                           help='The maximum library insert size defines different cut offs based on the maximum expected '
                           'library size. *This is not the average fragment size* but the higher end of the '
                           'the fragment size distribution (obtained using for example a Fragment Analyzer or a Bioanalyzer) '
                           'which usually is between 800 to 1500 bp. If this value is not known use the default value.'
                           ' The insert value is used to decide if two mates belong to the same fragment (by '
                           'checking if they are within this max insert size) and to decide if a mate is too far '
                           'away from the nearest restriction site'
                           ' (Default: %(default)s).',
                           type=int,
                           default=1000,
                           required=False)

    parserOpt.add_argument('--genomeAssembly', '-ga',
                           help='The genome the reads were mapped to. Used for metadata of cool file.')

    parserOpt.add_argument('--region', '-r',
                           help='Region of the genome to limit the operation to. '
                           'The format is chr:start-end. It is also possible to just '
                           'specify a chromosome, for example --region chr10',
                           metavar="CHR:START-END",
                           required=False,
                           type=genomicRegion
                           )
    parserOpt.add_argument('--keepSelfLigation',
                           help='If set, inward facing reads less than 1000 bp apart and having a restriction'
                           'site in between are removed. Although this reads do not contribute to '
                           'any distant contact, they are useful to account for bias in the data'
                           ' (for the moment is always True).',
                           #    help=argparse.SUPPRESS,
                           #    default=True
                           action='store_true'
                           )

    parserOpt.add_argument('--keepSelfCircles',
                           help='If set, outward facing reads without any restriction fragment (self circles) are kept. '
                           'They will be counted and shown in the QC plots.',
                           required=False,
                           action='store_true'
                           )

    parserOpt.add_argument('--minMappingQuality',
                           help='minimum mapping quality for reads to be accepted. '
                           'Because the restriction enzyme site could be located '
                           'on top of the read, this may reduce the '
                           'reported quality of the read. Thus, this parameter '
                           'may be adjusted if too many low quality '
                           '(but otherwise perfectly valid Hi-C reads) are found. '
                           'A good strategy is to make a test run (using the --doTestRun), '
                           'then checking the results to see if too many low quality '
                           'reads are present and then using the bam file generated to '
                           'check if those low quality reads are caused by the read '
                           'not being mapped entirely'
                           ' (Default: %(default)s).',
                           required=False,
                           default=15,
                           type=int
                           )
    parserOpt.add_argument('--threads',
                           help='Number of threads. Using the python multiprocessing module. '
                           'One master process which is used to read the input file into the buffer and one process which is merging '
                           'the output bam files of the processes into one output bam file. All other threads do the actual computation. '
                           'Minimum value for the \'--thread\' parameter is 2. '
                           'The usage of 8 threads is optimal if you have an HDD. A higher number of threads is only '
                           'useful if you have a fast SSD. Have in mind that the performance of hicBuildMatrix is influenced by '
                           'the number of threads, the speed of your hard drive and the inputBufferSize. To clarify: the performance '
                           'with a higher thread number is not negative influenced but not positive too. With a slow HDD and a high number of '
                           'threads many threads will do nothing most of the time'
                           ' (Default: %(default)s).',
                           required=False,
                           default=4,
                           type=int
                           )
    parserOpt.add_argument('--inputBufferSize',
                           help='Size of the input buffer of each thread. 400,000 read pairs per input file per thread is the default value. '
                           'Reduce this value to decrease memory usage.',
                           required=False,
                           default=400000,
                           type=int
                           )
    parserOpt.add_argument('--doTestRun',
                           help='A test run is useful to test the quality '
                           'of a Hi-C experiment quickly. It works by '
                           'testing only 1,000,000 reads. This option '
                           'is useful to get an idea of quality control '
                           'values like inter-chromosomal interactions, '
                           'duplication rates etc.',
                           action='store_true'
                           )
    parserOpt.add_argument('--doTestRunLines',
                           help='Number of lines to consider for the qc test run'
                           ' (Default: %(default)s).',
                           required=False,
                           default=1000000,
                           type=int
                           )

    parserOpt.add_argument('--skipDuplicationCheck',
                           help='Identification of duplicated read pairs is '
                           'memory consuming. Thus, in case of memory '
                           'errors this check can be skipped. However, '
                           'consider running a `--doTestRun` first to '
                           'get an estimation of the duplicated reads. ',
                           action='store_true'
                           )
    parserOpt.add_argument('--chromosomeSizes', '-cs',
                           help=('File with the chromosome sizes for your genome. A tab-delimited two column layout \"chr_name size\" is expected'
                                 'Usually the sizes can be determined from the SAM/BAM input files, however, '
                                 'for cHi-C or scHi-C it can be that at the start or end no data is present. '
                                 'Please consider that this option causes that only reads are considered which are on the listed chromosomes.'
                                 'Use this option to guarantee fixed sizes. An example file is available via UCSC: '
                                 'http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes'),
                           type=argparse.FileType('r'),
                           metavar='txt file')
    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    """
    Reads line by line two bam files that are not sorted.
    Each line in the two bam files should correspond
    to the mapped position of the two ends of a Hi-C
    fragment.

    Each mate pair is assessed to determine if it is
    a valid Hi-C pair, in such case a matrix
    reporting the counts of mates is constructed.

    A bam file containing the valid Hi-C reads
    is also constructed.
    """

    args = parse_arguments().parse_args(args)

    createMatrix(pOutFileName=args.outFileName, pMaxDistance=args.maxDistance, pMaxLibraryInsertSize=args.maxLibraryInsertSize, pQCfolder=args.QCfolder,
                 pThreads=args.threads, pDanglingSequence=args.danglingSequence, pRestrictionSequence=args.restrictionSequence, pSamFiles=args.samFiles,
                 pDoTestRun=args.doTestRun, pOutBam=args.outBam, pChromosomeSizes=args.chromosomeSizes, pRestrictionCutFile=args.restrictionCutFile,
                 pRegion=args.region, pBinSize=args.binSize, pInputBufferSize=args.inputBufferSize, pMinDistance=args.minDistance,
                 pDoTestRunLines=args.doTestRunLines, pSkipDuplicationCheck=args.skipDuplicationCheck, pMinMappingQuality=args.minMappingQuality,
                 pKeepSelfCircles=args.keepSelfCircles, pKeepSelfLigation=args.keepSelfLigation, pGenomeAssembly=args.genomeAssembly)
