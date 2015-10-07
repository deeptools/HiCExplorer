import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='returns a bed file '
                                     'containing the positions '
                                     'of a string from a given fasta file',
                                     usage='An example usage is: %(prog)s --fasta mm10.fa '
                                           '--searchPattern AAGCTT -o rest_site_positions.bed')
    
    # define the arguments
    parser.add_argument('--fasta', '-f',
                        help='Path to fasta file',
                        type=argparse.FileType('r'),
                        required=True)

    # define the arguments
    parser.add_argument('--searchPattern', '-p',
                        help='Search pattern. For example, for HindII this pattern is "AAGCTT". '
                             'Both, forward and reverse strand are searched for a match.',
                        required=True)

    parser.add_argument('--outFile', '-o',
                        help='Name for the resulting bed file',
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
    >>> fa.write(">chr1\nCGTACGAACGTACGGTACGcgtaGTACGGCATT\n")
    >>> fa.close()
    >>> find_pattern("CGTA", "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t0\t4\t.\t0\t+\n', 'chr1\t8\t12\t.\t0\t+\n', 'chr1\t19\t23\t.\t0\t+\n', 'chr1\t28\t32\t.\t0\t-\n']
    """

    rev_compl = str(Seq(pattern, generic_dna).complement())
    for record in SeqIO.parse(fasta_file, 'fasta', generic_dna):
        # find all the occurrences of pattern
        for match in re.finditer(pattern, str(record.seq), re.IGNORECASE):
            out_file.write('{}\t{}\t{}\t.\t0\t+\n'.format(record.name,
                                                     match.start(),
                                                     match.end()))
        for match in re.finditer(rev_compl, str(record.seq), re.IGNORECASE):
            out_file.write('{}\t{}\t{}\t.\t0\t-\n'.format(record.name,
                                                     match.start(),
                                                     match.end()))



def main():
    args = parse_arguments().parse_args()
    find_pattern(args.searchPattern, args.fasta, args.outFile)

