"""Reference files whose row order is that of the C++ standard library.

bedtools sort orders a chromosome's rows by start alone with the non-stable
std::sort, and hicValidateLocations and hicMergeLoops reproduce it
(cpp/core/include/hicx/bedtools_ops.hpp). Rows that share a start come out in
the order the standard library's std::sort leaves them in, which differs
between libstdc++ (Linux) and libc++ (macOS), for bedtools as for these
tools. The references next to the inputs are those of libstdc++; the ones of
libc++ are in a libc++ directory beside them, made by the same tools built
against libc++ 19.
"""
import os

from hicexplorer._hicx import STANDARD_LIBRARY


def sort_tie_reference(pPath):
    """`pPath`, or its libc++ counterpart when the tools are built against libc++."""
    if STANDARD_LIBRARY == 'libc++':
        directory, name = os.path.split(pPath)
        return os.path.join(directory, 'libc++', name)
    return pPath
