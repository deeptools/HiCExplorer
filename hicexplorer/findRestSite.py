import sys
import argparse
import re
from tempfile import NamedTemporaryFile
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Identifies the genomic locations of restriction sites. ',
                                     usage='An example usage is: %(prog)s --fasta mm10.fa '
                                           '--searchPattern AAGCTT -o rest_site_positions.bed')
    
    # define the arguments
    parser.add_argument('--fasta', '-f',
                        help='Path to fasta file for the organism genome.',
                        type=argparse.FileType('r'),
                        required=True)

    # define the arguments
    parser.add_argument('--searchPattern', '-p',
                        help='Search pattern. For example, for HindII this pattern is "AAGCTT". '
                             'Both, forward and reverse strand are searched for a match. The pattern '
                             'is a regexp and can contain regexp specif syntax '
                             '(see https://docs.python.org/2/library/re.html). For example the pattern'
                             'CG..GC will find all occurrence of CG followed by any two bases and then GC.',
                        required=True)

    parser.add_argument('--outFile', '-o',
                        help='Name for the resulting bed file.',
                        type=argparse.FileType('w'),
                        required=True)

    return parser


def find_pattern(pattern, fasta_file, out_file):
    r"""
    Finds the occurences of the match in the fasta file
    and saves a bed file.

    The coordinate system is zero based

    :param pattern: Sequence to search for
    :param fasta_file:
    :param out_file: file handler

    :return: none

    >>> fa = open("/tmp/test.fa", 'w')
    >>> fa.write(">chr1\nGCATCGTACGAACGTACGGTACGcgtaCGNAGT\n")
    >>> fa.close()
    >>> find_pattern("CGTA", "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t0\t4\t.\t0\t-\n', 'chr1\t4\t8\t.\t0\t+\n', 'chr1\t12\t16\t.\t0\t+\n', 'chr1\t23\t27\t.\t0\t+\n']
    >>> find_pattern("CG.AG", "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t0\t5\t.\t0\t-\n', 'chr1\t27\t32\t.\t0\t+\n']
    """

    rev_compl = str(Seq(pattern, generic_dna).complement())
    temp = NamedTemporaryFile(suffix=".bed", delete=False)
    for record in SeqIO.parse(fasta_file, 'fasta', generic_dna):
        # find all the occurrences of pattern
        for match in re.finditer(pattern, str(record.seq), re.IGNORECASE):
            temp.write('{}\t{}\t{}\t.\t0\t+\n'.format(record.name,
                                                      match.start(),
                                                      match.end()))
        for match in re.finditer(rev_compl, str(record.seq), re.IGNORECASE):
            temp.write('{}\t{}\t{}\t.\t0\t-\n'.format(record.name,
                                                      match.start(),
                                                      match.end()))

    sys.stderr.write("Sorting file ...\n")
    tmpfile_name = temp.name
    temp.close()
    subprocess.check_output(["cat", tmpfile_name])
    # sort bed file using system tools
    cmd = 'sort -k1,1V -k2,2n {}'.format(tmpfile_name)
    # LC_ALL=C is to set the appropriate collation order
    proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, env={'LC_ALL':' C'})
    stdout, _ = proc.communicate()

    out_file.write(stdout)
    out_file.close()

    temp.close()


def main():
    args = parse_arguments().parse_args()
    find_pattern(args.searchPattern, args.fasta, args.outFile)
