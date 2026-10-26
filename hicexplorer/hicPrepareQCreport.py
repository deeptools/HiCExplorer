"""hicPrepareQCreport: a thin layer over the C++ executable hicPrepareQCreport.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import os
import errno
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicPrepareQCreport, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPrepareQCreport')


def save_html(filename, unmap_table, discard_table, distance_table, orientation_table, all_table):
    root = os.path.dirname(os.path.abspath(__file__))

    html = open(os.path.join(root, "qc_template.html"), "r")
    html_content = html.read()
    # the html code has a placeholder for the html table
    html_content = html_content.replace("%%TABLE_UNMAP%%", unmap_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).to_html(classes='df'))
    html_content = html_content.replace("%%TABLE_DISCARDED%%", discard_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).to_html(classes='df'))
    html_content = html_content.replace("%%TABLE_DISTANCE%%", distance_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).to_html(classes='df'))
    html_content = html_content.replace("%%TABLE_ORIENTATION%%", orientation_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).to_html(classes='df'))

    if 'Min rest. site distance' in all_table.columns:
        all_table = all_table.drop(['Min rest. site distance'], axis=1)

    if 'Max library insert size' in all_table.columns:
        all_table = all_table.drop(['Max library insert size'], axis=1)
    # all_table = all_table.drop(['Min rest. site distance', 'Max library insert size'], axis=1)

    html_content = html_content.replace("%%TABLE%%", all_table.style.to_html(classes='df'))
    with open(filename, 'w') as fh:
        fh.write(html_content)
    html.close()


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_figure_pairs_used(table, filename, dpi):
    prc_table = table[[
        'Hi-C contacts', 'Pairs mappable, unique and high quality', 'Sequenced reads']] / 1e6

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
    prc_table = table[['Hi-C contacts', 'Low mapping quality', 'One mate not unique',
                       'One mate unmapped']].T / table['Sequenced reads']

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

    # merge the counts table with the percentages table
    ret_table = table[['Hi-C contacts', 'Low mapping quality', 'One mate not unique',
                       'One mate unmapped']].join(prc_table.T, rsuffix='_%')

    return ret_table[[u'Hi-C contacts', u'Hi-C contacts_%', u'Low mapping quality',
                      u'Low mapping quality_%', u'One mate not unique',
                      u'One mate not unique_%',
                      u'One mate unmapped', u'One mate unmapped_%']]


def make_figure_pairs_discarded(table, filename, dpi):
    column_names_prefix = ['One mate not close to rest site', 'dangling end', 'duplicated pairs',
                           'same fragment', 'self circle',
                           'self ligation (removed)']
    column_names = []
    column_names_out = []

    column_names_table = list(table.columns)
    for name in column_names_prefix:
        if name in column_names_table:
            column_names.append(name)
            column_names_out.append(name)
            column_names_out.append(name + ' %')

        else:
            # if 'dangling end' in
            for name_ in column_names_table:
                if name in name_:
                    column_names.append(name_)
                    column_names_out.append(name_ + ' %')

    log.debug('column_names {}'.format(column_names))
    prc_table = table[column_names].T / table['Pairs mappable, unique and high quality']

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

    # merge the counts table with the percentages table
    ret_table = table[column_names].join(prc_table.T, rsuffix=' %')

    return ret_table[column_names_out]


def make_figure_distance(table, filename):

    prc_table2 = table[['inter chromosomal',
                        'Intra short range (< 20kb)', 'Intra long range (>= 20kb)']].T / table['Hi-C contacts']
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

    # merge the counts table with the percentages table
    ret_table = table[['inter chromosomal', 'Intra short range (< 20kb)', 'Intra long range (>= 20kb)']].join(
        prc_table2.T, rsuffix=' %')

    return ret_table[['inter chromosomal', 'inter chromosomal %', 'Intra short range (< 20kb)', 'Intra short range (< 20kb) %', 'Intra long range (>= 20kb)', 'Intra long range (>= 20kb) %']]


def make_figure_read_orientation(table, filename, dpi):
    _t = table[[u'Read pair type: inward pairs', u'Read pair type: outward pairs',
                u'Read pair type: left pairs', u'Read pair type: right pairs']].T
    prc_table3 = _t / _t.sum(axis=0)
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

    # merge the counts table with the percentages table
    ret_table = table[[u'Read pair type: inward pairs', u'Read pair type: outward pairs',
                       u'Read pair type: left pairs', u'Read pair type: right pairs']].join(prc_table3.T, rsuffix=' %')

    return ret_table[[u'Read pair type: inward pairs', u'Read pair type: inward pairs %', u'Read pair type: outward pairs', u'Read pair type: outward pairs %',
                      u'Read pair type: left pairs', u'Read pair type: left pairs %', u'Read pair type: right pairs', u'Read pair type: right pairs %']]


def main(args=None):
    """Run the C++ implementation of hicPrepareQCreport with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPrepareQCreport')(args)
