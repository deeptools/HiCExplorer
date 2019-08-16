#!/usr/bin/env python
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
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Tabulates and plots QC measures from  '
                                                 'hicBuildMatrix log files within an HTML output',
                                     add_help=False,
                                     usage='%(prog)s --logfiles matrix1_QCfolder/QC.log matrix2_QCfolder/QC.log '
                                           '--labels "sample 1" "sample 2" --outputFolder QC_all_samples)')

    parserRequired = parser.add_argument_group('Required arguments')
    # define the arguments
    parserRequired.add_argument('--logfiles', '-l',
                                help='Path to the log files to be processed',
                                type=argparse.FileType('r'),
                                nargs="+",
                                required=True)

    parserRequired.add_argument('--labels',
                                help='Label to assign to each log file. Each label should be separated by a space. Quote '
                                'labels that contain spaces: E.g. --labels label1 "labels 2"',
                                nargs="+")

    parserRequired.add_argument('--outputFolder', '-o',
                                help='Several files with be saved under this folder: A table containing the results and '
                                'a html file with several images.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--dpi',
                           help='Image resolution. By default high resolution png images with a 200 dpi are created.',
                           type=int,
                           default=200)

    parserOpt.add_argument("--help", "-h", action="help",
                           help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def save_html(filename, unmap_table, discard_table, distance_table, orientation_table, all_table):
    root = os.path.dirname(os.path.abspath(__file__))

    html = open(os.path.join(root, "qc_template.html"), "r")
    html_content = html.read()
    # the html code has a placeholder for the html table
    html_content = html_content.replace("%%TABLE_UNMAP%%", unmap_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).render())
    html_content = html_content.replace("%%TABLE_DISCARDED%%", discard_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).render())
    html_content = html_content.replace("%%TABLE_DISTANCE%%", distance_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).render())
    html_content = html_content.replace("%%TABLE_ORIENTATION%%", orientation_table.style
                                        .format(lambda x: '{:,}'.format(x) if x > 1 else '{:.2%}'.format(x)).render())

    all_table = all_table[['Sequenced reads', 'Pairs mappable, unique and high quality', 'Hi-C contacts',
                           'One mate unmapped', 'One mate not unique', 'Low mapping quality', 'dangling end',
                           'self ligation (removed)', 'One mate not close to rest site', 'same fragment',
                           'self circle', 'duplicated pairs', 'inter chromosomal', 'Intra short range (< 20kb)',
                           'Intra long range (>= 20kb)', 'Read pair type: inward pairs', 'Read pair type: outward pairs', 'Read pair type: left pairs', 'Read pair type: right pairs']]

    html_content = html_content.replace("%%TABLE%%", all_table.style.render())
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

    # merge the counts table with the percentages table
    ret_table = table[['Hi-C contacts', 'Low mapping quality', 'One mate not unique',
                       'One mate unmapped']].join(prc_table.T, rsuffix='_%')

    return ret_table[[u'Hi-C contacts', u'Hi-C contacts_%', u'Low mapping quality',
                      u'Low mapping quality_%', u'One mate not unique',
                      u'One mate not unique_%',
                      u'One mate unmapped', u'One mate unmapped_%']]


def make_figure_pairs_discarded(table, filename, dpi):
    prc_table = table[['One mate not close to rest site', 'dangling end', 'duplicated pairs',
                       'same fragment', 'self circle',
                       'self ligation (removed)']].T / table['Pairs mappable, unique and high quality']

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

    # merge the counts table with the percentages table
    ret_table = table[['One mate not close to rest site', 'dangling end', 'duplicated pairs',
                       'same fragment', 'self circle',
                       'self ligation (removed)']].join(prc_table.T, rsuffix=' %')

    return ret_table[['One mate not close to rest site', 'One mate not close to rest site %',
                      'dangling end', 'dangling end %', 'duplicated pairs', 'duplicated pairs %',
                      'same fragment', 'same fragment %',
                      'self circle', 'self circle %', 'self ligation (removed)', 'self ligation (removed) %']]


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

    # merge the counts table with the percentages table
    ret_table = table[[u'Read pair type: inward pairs', u'Read pair type: outward pairs',
                       u'Read pair type: left pairs', u'Read pair type: right pairs']].join(prc_table3.T, rsuffix=' %')

    return ret_table[[u'Read pair type: inward pairs', u'Read pair type: inward pairs %', u'Read pair type: outward pairs', u'Read pair type: outward pairs %',
                      u'Read pair type: left pairs', u'Read pair type: left pairs %', u'Read pair type: right pairs', u'Read pair type: right pairs %']]


def main(args=None):
    """
    The structure of the log file is as follows:
    --------------------------------------------

    File    /tmp/test
    Sequenced reads        99983
    Min rest. site distance 150
    Max rest. site distance 1500


    #       count   (percentage w.r.t. total sequenced reads)
    Pairs mappable, unique and high quality 52726   (52.73)
    Hi-C contacts      36552   (36.56)
    One mate unmapped       8777    (8.78)
    One mate not unique     3603    (3.60)
    Low mapping quality    34877   (34.88)

    #       count   (percentage w.r.t. mappable, unique, high quality pairs)
    dangling end    209     (0.40)
    self ligation (removed) 5056    (9.59)
    One mate not close to rest site 751     (1.42)
    same fragment  10146   (19.24)
    self circle     4274    (8.11)
    duplicated pairs        12      (0.02)

    #       count   (percentage w.r.t. total valid Hi-C contacts)
    inter chromosomal       5849    (16.00)
    Intra short range (< 20kb)      8689    (23.77)
    Intra long range (>= 20kb)      22014   (60.23)
    Read pair type: inward pairs    6977    (19.09)
    Read pair type: outward pairs   9525    (26.06)
    Read pair type: left pairs      7012    (19.18)
    Read pair type: right pairs     7189    (19.67)
    """

    args = parse_arguments().parse_args(args)
    params = dict()
    make_sure_path_exists(args.outputFolder)
    for fh in args.logfiles:
        in_log_part = False
        log.debug('Processing {}\n'.format(fh.name))
        for line in fh.readlines():
            if line.startswith("File"):
                in_log_part = True
            if in_log_part is True:
                if line.strip() == "" or line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) == 1:
                    continue
                if fields[0] not in params:
                    params[fields[0]] = []

                try:
                    params[fields[0]].append(int(fields[1]))
                except ValueError:
                    params[fields[0]].append(fields[1])

    table = pd.DataFrame(params)
    if args.labels and len(args.labels) == len(args.logfiles):
        try:
            table['Labels'] = args.labels
        except ValueError:
            log.error("*ERROR* Some log files may not be valid. Please check that the log files contain "
                      "at the end the summary information.")
            exit()

        table = table.set_index('Labels')
    else:
        table = table.set_index('File')

    if 'Pairs mappable, unique and high quality' not in table.columns:
        table['Pairs mappable, unique and high quality'] = \
            table['Sequenced reads'] - (table['One mate unmapped'] +
                                        table['One mate not unique'] + table['Low mapping quality'])

    if 'same fragment (800 bp)' in table.columns:
        # older versions of the QC used the label 'same fragment (800 bp)'
        table['same fragment'] = table['same fragment (800 bp)']

    make_figure_pairs_used(table, args.outputFolder +
                           "/pairs_sequenced.png", args.dpi)
    unmap_table = make_figure_umappable_non_unique_reads(table, args.outputFolder + "/unmappable_and_non_unique.png",
                                                         args.dpi)

    discarded_table = make_figure_pairs_discarded(
        table, args.outputFolder + "/pairs_discarded.png", args.dpi)
    distance_table = make_figure_distance(
        table, args.outputFolder + "/distance.png")
    read_orientation_table = make_figure_read_orientation(
        table, args.outputFolder + "/read_orientation.png", args.dpi)

    save_html(args.outputFolder + "/hicQC.html", unmap_table, discarded_table, distance_table,
              read_orientation_table, table)

    unmap_table.to_csv(args.outputFolder + "/unmapable_table.txt", sep="\t")
    discarded_table.to_csv(args.outputFolder +
                           "/discarded_table.txt", sep="\t")
    distance_table.to_csv(args.outputFolder + "/distance_table.txt", sep="\t")
    read_orientation_table.to_csv(
        args.outputFolder + "/read_orientation_table.txt", sep="\t")
    table.to_csv(args.outputFolder + "/QC_table.txt", sep="\t")
