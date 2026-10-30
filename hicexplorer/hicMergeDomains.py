"""hicexplorer.hicMergeDomains, implemented in C++.

main() runs the C++ executable hicMergeDomains and parse_arguments() is built
from its --help-json (hicexplorer._cpp). Every other function keeps its 3.7.6
signature, docstring and behaviour, and its body calls
hicexplorer._hicx.hicMergeDomains (cpp/python/src/api/hicMergeDomains.cpp,
cpp/core/include/hicx/merge_domains.hpp). The deviations from 3.7.6 are listed
in hicexplorer/test/api/cases_hicMergeDomains.py.
"""
import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicMergeDomains as _impl


def parse_arguments(args=None):
    """The argparse parser of hicMergeDomains, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicMergeDomains')


def create_list_of_file(file, binSizeList):
    """Creates a list from a file and recognizes its BinSize at the same time"""
    return _impl.create_list_of_file(file, binSizeList)


def merge_list(d1, d2, pValue=5000):
    """Combine two lists into one, taking into account the value,
       with which one can decide whether two very similar TADs should
       be seen as one or as two """
    return _impl.merge_list(d1, d2, pValue)


def add_id(domainList):
    """Give each element of a list its individual number """
    return _impl.add_id(domainList)


def create_relationsship_list(domainList, pPercent=0.5):
    """Give each element of a list its individual number """
    return _impl.create_relationsship_list(domainList, pPercent)


def add_relation_to_list(rList, chromosom, parent, child):
    """ Add a discovered relationship to a list """
    return _impl.add_relation_to_list(rList, chromosom, parent, child)


def write_in_file(pList, name, relation=False):
    """ Create and save various lists in their individual files """
    return _impl.write_in_file(pList, name, relation)


def create_tree(rList, dList, pFileNameSuffix, pFormat):
    """ Create a tree diagram from a list and
       thus represent the relationships of individual TADs per chromosome """
    return _impl.create_tree(rList, dList, pFileNameSuffix, pFormat)


def create_small_list(dList):
    """ Creates a nested list from a domain file, in which the TADs
    are sorted according to their chromosomes """
    return _impl.create_small_list(dList)


def read_protein(file):
    """ Reads a protein file and saves it in a list """
    return _impl.read_protein(file)


def merge_protein(pProteinList, binSize, minPeak):
    """ Adapts the protein list to the corresponding BinSize, taking into
       account the minimum number of peaks per section """
    return _impl.merge_protein(pProteinList, binSize, minPeak)


def compare_boundaries_protein(bList, cList, paraScore=0.2):
    """ check the boundaries for a protein-peak """
    return _impl.compare_boundaries_protein(bList, cList, paraScore)


def create_list_with_protein(bList, minPeak, cList=None):
    return _impl.create_list_with_protein(bList, minPeak, cList)


def main(args=None):
    """Run the C++ implementation of hicMergeDomains with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicMergeDomains')(args)
