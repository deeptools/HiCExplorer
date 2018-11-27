import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import re
from tempfile import NamedTemporaryFile
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Identifies the genomic locations of restriction sites. ',
                                     add_help=False,
                                     usage='%(prog)s --fasta mm10.fa '
                                           '--searchPattern AAGCTT -o rest_site_positions.bed')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--fasta', '-f',
                                help='Path to fasta file for the organism genome.',
                                type=argparse.FileType('r'),
                                required=True)

    # define the arguments
    parserRequired.add_argument('--searchPattern', '-p',
                                help='Search pattern. For example, for HindIII this pattern is "AAGCTT". '
                                'Both, forward and reverse strand are searched for a match. The pattern '
                                'is a regexp and can contain regexp specif syntax '
                                '(see https://docs.python.org/2/library/re.html). For example the pattern'
                                'CG..GC will find all occurrence of CG followed by any two bases and then GC.',
                                required=True)

    parserRequired.add_argument('--outFile', '-o',
                                help='Name for the resulting bed file.',
                                type=argparse.FileType('w'),
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser


def find_pattern(pattern, fasta_file, out_file):
    r"""
    Finds the occurrences of the match in the fasta file
    and saves a bed file.

    The coordinate system is zero based

    :param pattern: Sequence to search for
    :param fasta_file:
    :param out_file: file handler

    :return: none

    This is an example from Emboss remap for Csp6I with
    restriction site GTAC


                       Csp6I
               Csp6I   |    Csp6I Csp6I
               \       \    \     \
          CTACGGTACGAACGTACGGTACGcgtaCGNAGTCATG
                   10        20        30
          ----:----|----:----|----:----|----:--
          GATGCCATGCTTGCATGCCATGCgcatGCNTCAGTAC
                 /       /    /     /
                 Csp6I   |    Csp6I Csp6I
                         Csp6I

    >>> fa = open("/tmp/test.fa", 'w')
    >>> foo = fa.write(">chr1\nCTACGGTACGAACGTACGGTACGcgtaCGNAGTCATG\n")
    >>> fa.close()
    >>> find_pattern("GTAC", "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t5\t9\t.\t0\t+\n', 'chr1\t13\t17\t.\t0\t+\n', 'chr1\t18\t22\t.\t0\t+\n', 'chr1\t24\t28\t.\t0\t+\n']

    Test with non palindromic sequence with regexp
    >>> find_pattern("CG.AG", "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t0\t5\t.\t0\t-\n', 'chr1\t27\t32\t.\t0\t+\n']

    """

    # get the reverse complement of the pattern
    rev_compl = str(Seq(pattern, generic_dna).reverse_complement())

    temp = NamedTemporaryFile(suffix=".bed", delete=False, mode='wt')
    for record in SeqIO.parse(fasta_file, 'fasta', generic_dna):
        # find all the occurrences of pattern
        for match in re.finditer(pattern, str(record.seq), re.IGNORECASE):
            _ = temp.write('{}\t{}\t{}\t.\t0\t+\n'.format(record.name,
                                                          match.start(),
                                                          match.end()))
        if rev_compl != pattern:
            # search for the reverse complement only if the pattern is not palindromic
            for match in re.finditer(rev_compl, str(record.seq), re.IGNORECASE):
                _ = temp.write('{}\t{}\t{}\t.\t0\t-\n'.format(record.name,
                                                              match.start(),
                                                              match.end()))
    log.info("Sorting file ...")
    tmpfile_name = temp.name
    temp.close()
    subprocess.check_output(["cat", tmpfile_name])
    # sort bed file using system tools
    cmd = 'sort -k1,1 -k2,2n -u {}'.format(tmpfile_name)
    # LC_ALL=C is to set the appropriate collation order
    proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, env={'LC_ALL': ' C'}, universal_newlines=True)
    stdout, _ = proc.communicate()

    out_file.write(stdout)
    out_file.close()

    temp.close()


def main():
    args = parse_arguments().parse_args()
    find_pattern(args.searchPattern, args.fasta, args.outFile)
