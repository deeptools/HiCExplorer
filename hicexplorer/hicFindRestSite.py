"""hicFindRestSite: a thin layer over the C++ executable hicFindRestSite.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). find_pattern calls hicexplorer._hicx.hicFindRestSite
(C++, cpp/python/src/api/hicFindRestSite.cpp), the same scan the executable
runs.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import re
from tempfile import NamedTemporaryFile
import subprocess

import gzip
from mimetypes import guess_type
from functools import partial
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicFindRestSite as _impl


def parse_arguments(args=None):
    """The argparse parser of hicFindRestSite, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicFindRestSite')


def find_pattern(pPattern, fasta_file, out_file):
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
    >>> find_pattern(["GTAC"], "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t5\t9\t.\t0\t+\n', 'chr1\t13\t17\t.\t0\t+\n', 'chr1\t18\t22\t.\t0\t+\n', 'chr1\t24\t28\t.\t0\t+\n']

    Test with non palindromic sequence with regexp
    >>> find_pattern(["CG.AG"], "/tmp/test.fa", open("/tmp/test.bed", 'w'))
    >>> open("/tmp/test.bed", 'r').readlines()
    ['chr1\t0\t5\t.\t0\t-\n', 'chr1\t27\t32\t.\t0\t+\n']

    """

    return _impl.find_pattern(pPattern, fasta_file, out_file)


def main(args=None):
    """Run the C++ implementation of hicFindRestSite with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicFindRestSite')(args)
