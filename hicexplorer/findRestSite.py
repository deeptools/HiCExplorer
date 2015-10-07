import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='returns a bed file'
                                     'containing the positions '
                                     'of a string, giving a fasta file',
                                     usage='An example usage is: %(prog)s --fasta mm10.fa '
                                           '--searchPattern AAGCTT -o rest_site_positions.bed')
    
    # define the arguments
    parser.add_argument('--fasta', '-f',
                        help='Path to fasta file',
                        type=argparse.FileType('r'),
                        required=True)

    # define the arguments
    parser.add_argument('--searchPattern', '-p',
                        help='Search pattern. For example for HindII this pattern is "AAGCTT". '
                             'Both, forward and reverse strand are searched for a match.',
                        required=True)

    parser.add_argument('--outFile', '-o',
                        help='Output file',
                        type=argparse.FileType('w'),
                        required=True)

    return parser


def main():
    args = parse_arguments().parse_args()

    rev_compl = str(Seq(args.searchPattern, generic_dna).complement())
    for record in SeqIO.parse(args.fasta, 'fasta', generic_dna):
        # find all the occurrences of pattern
        for match in re.finditer(args.searchPattern, str(record.seq), re.IGNORECASE):
            args.outFile.write('{}\t{}\t{}\t.\t0\t+\n'.format(record.name,
                                                     match.start(),
                                                     match.end()))
        for match in re.finditer(rev_compl, str(record.seq), re.IGNORECASE):
            args.outFile.write('{}\t{}\t{}\t.\t0\t.\n'.format(record.name,
                                                     match.start(),
                                                     match.end()))
