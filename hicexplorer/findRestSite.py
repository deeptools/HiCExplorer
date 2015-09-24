#!/usr/bin/env python
#-*- coding: utf-8 -*-

import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parseArguments(args=None):
    parser = argparse.ArgumentParser(description = 'returns a bed file'
                                     'containing the positions '
                                     'of a string, giving a fasta file')
    
    # define the arguments
    parser.add_argument('--fasta', '-f',
                        help='Path to fasta file',
                        type=argparse.FileType('r'),
                        required = True)

    # define the arguments
    parser.add_argument('--searchPattern', '-p',
                        help='Search pattern',
                        required = True)

    parser.add_argument('--outFile', '-o',
                        help = 'Output file',
                        type=argparse.FileType('w'),
                        required=True)



    return(parser.parse_args(args))


def main(args):
    rev_compl = str(Seq(args.searchPattern, generic_dna).complement())
    for record in SeqIO.parse(args.fasta, 'fasta', generic_dna):
        # find all the occurrences of pattern
        for match in re.finditer(args.searchPattern, str(record.seq)):
            args.outFile.write('{}\t{}\t{}\n'.format(record.name, 
                                                     match.start(),
                                                     match.end()))
        for match in re.finditer(rev_compl, str(record.seq)):
            args.outFile.write('{}\t{}\t{}\n'.format(record.name, 
                                                     match.start(),
                                                     match.end()))
if __name__ == "__main__":
    args = parseArguments()
    main(args)
