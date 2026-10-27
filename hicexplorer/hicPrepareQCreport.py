"""hicPrepareQCreport: a thin layer over the C++ executable hicPrepareQCreport.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool. save_html and
make_sure_path_exists run in C++ (hicexplorer._hicx.hicPrepareQCreport); the
make_figure_* functions draw with matplotlib and pandas from the tables C++
computes.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)

from hicexplorer._hicx import hicPrepareQCreport as _impl


def parse_arguments(args=None):
    """The argparse parser of hicPrepareQCreport, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPrepareQCreport')


def save_html(filename, unmap_table, discard_table, distance_table, orientation_table, all_table):
    return _impl.save_html(filename, unmap_table, discard_table, distance_table, orientation_table, all_table)


def make_sure_path_exists(path):
    return _impl.make_sure_path_exists(path)


def make_figure_pairs_used(table, filename, dpi):
    prc_table = _impl.pairs_used_data(table)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    prc_table.plot(kind='barh', ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left',
                    bbox_to_anchor=(1, 0.5))
    ax.set_xlabel("Number of reads in millions")
    ax.set_ylabel("")
    plt.savefig(filename, bbox_extra_artists=(
        lgd,), bbox_inches='tight', dpi=dpi)
    plt.close()


def make_figure_umappable_non_unique_reads(table, filename, dpi):
    prc_table, ret_table = _impl.umappable_non_unique_reads_data(table)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    prc_table.plot.bar(ax=ax)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left',
                    bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction w.r.t. pairs sequenced")
    plt.savefig(filename, bbox_extra_artists=(
        lgd,), bbox_inches='tight', dpi=dpi)
    plt.close()

    # the counts table merged with the percentages table (C++)
    return ret_table


def make_figure_pairs_discarded(table, filename, dpi):
    prc_table, ret_table = _impl.pairs_discarded_data(table)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    prc_table.plot.bar(ax=ax)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left',
                    bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction w.r.t. mappable and unique pairs")
    plt.savefig(filename, bbox_extra_artists=(
        lgd,), bbox_inches='tight', dpi=dpi)
    plt.close()

    # the counts table merged with the percentages table (C++)
    return ret_table


def make_figure_distance(table, filename):

    prc_table2, ret_table = _impl.distance_data(table)
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    prc_table2.plot.bar(ax=ax)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left',
                    bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction w.r.t. valid Hi-C contacts")

    plt.savefig(filename, bbox_extra_artists=(
        lgd,), bbox_inches='tight', dpi=200)
    plt.close()

    # the counts table merged with the percentages table (C++)
    return ret_table


def make_figure_read_orientation(table, filename, dpi):
    prc_table3, ret_table = _impl.read_orientation_data(table)
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    prc_table3.plot.bar(ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left',
                    bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction w.r.t. valid Hi-C contacts")
    plt.savefig(filename, bbox_extra_artists=(
        lgd,), bbox_inches='tight', dpi=dpi)
    plt.close()

    # the counts table merged with the percentages table (C++)
    return ret_table


def main(args=None):
    """Run the C++ implementation of hicPrepareQCreport with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPrepareQCreport')(args)
