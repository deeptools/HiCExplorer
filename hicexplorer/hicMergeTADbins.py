"""hicexplorer.hicMergeTADbins, implemented in C++.

main() runs the C++ executable hicMergeTADbins and parse_arguments() is built
from its --help-json (hicexplorer._cpp). merge_tad_bins and
get_boundary_bin_id keep their 3.7.6 signatures, docstrings and behaviour, and
their bodies call hicexplorer._hicx.hicMergeTADbins
(cpp/python/src/api/hicMergeTADbins.cpp, cpp/core/include/hicx/reduce_matrix.hpp).
The deviations from 3.7.6 are listed in
hicexplorer/test/api/cases_hicMergeTADbins.py.
"""
# -*- coding: utf-8 -*-
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicMergeTADbins as _impl


def parse_arguments(args=None):
    """The argparse parser of hicMergeTADbins, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicMergeTADbins')


def merge_tad_bins(hic, boundary_id_list, filename):
    """
    Reduces the HiCMatrix by merging the counts of tad bins.
    :param hic: HiCMatrix object
    :param boundary_id_list list of tad boundary bin ids
    :param filename Name to save the resulting matrix
    :return: HiCMatrix object
    """
    return _impl.merge_tad_bins(hic, boundary_id_list, filename)


def get_boundary_bin_id(hic, bed_fh):
    """
    :param hic: HiCMatrix object
    :param bed_fh: file handle of the bed file
    :return: Sorted list of bin indices.
    """
    return _impl.get_boundary_bin_id(hic, bed_fh)


def main(args=None):
    """Run the C++ implementation of hicMergeTADbins with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicMergeTADbins')(args)
