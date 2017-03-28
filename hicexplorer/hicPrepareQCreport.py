#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import errno
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hicexplorer._version import __version__


def parse_arguments():
    parser = argparse.ArgumentParser(description='Tabulates and plots QC measures from  '
                                                 'hicBuildMatrix log files.',
                                     usage='%(prog)s --logfiles matrix1.h5 matrix2.h5 '
                                           '--labels "sample 1" "sample 2" -outFile QC.txt)')

    # define the arguments
    parser.add_argument('--logfiles', '-l',
                        help='Path to the log files to be processed',
                        type=argparse.FileType('r'),
                        nargs="+",
                        required=True)

    parser.add_argument('--labels',
                        help='Label to assign to each log file. Each label should be separated by a space. Quote '
                             'labels that contain spaces: E.g. --labels label1 "labels 2"',
                        nargs="+")

    parser.add_argument('--outputFolder', '-o',
                        help='Several files with be saved under this folder: A table containing the results and '
                             'a html file with several images.',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def save_html(filename, html_table):
    root = os.path.dirname(os.path.abspath(__file__))

    html = open(os.path.join(root, "qc_template.html"), "r").read()
    # the html code has a placeholder for the html table
    html = html.replace("%%TABLE%%", html_table)
    with open(filename, 'w') as fh:
        fh.write(html)


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_figure_pairs_considered(table, filename):

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    table['Pairs considered'].plot.barh(ax=ax, color='#999999')
    xticks = ax.get_xticks()
    if xticks[-1] > 1e6:
        labels = ["{:.0f}".format(float(x) / 1e6) for x in xticks]
        labels[-1] += 'M'
    else:
        labels = ["{:,}".format(int(x)) for x in xticks]

    ax.set_xticklabels(labels)
    ax.set_ylabel("")
    ax.set_xlabel("Number of reads")
    plt.tight_layout()
    plt.savefig(filename)


def make_figure_pairs_used(table, filename):
    prc_table = table[['Pairs used', 'One mate low quality', 'One mate not unique',
                       'One mate not close to rest site', 'One mate unmapped', 'dangling end', 'duplicated pairs',
                       'same fragment (800 bp)', 'self circle',
                       'self ligation (removed)']].T / table['Pairs considered']

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    prc_table.plot.bar(ax=ax)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction")
    plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight')


def make_figure_distance(table, filename):

    prc_table2 = table[['inter chromosomal', 'short range < 20kb', 'long range']].T / table['Pairs used']
    fig = plt.figure(figsize=(4, 5))
    ax = fig.add_subplot(111)
    prc_table2.plot.bar(ax=ax)
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction")

    plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight')


def make_figure_read_orientation(table, filename):
    _t = table[[u'inward pairs', u'outward pairs', u'left pairs', u'right pairs']].T
    prc_table3 = _t / _t.sum(axis=0)
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    prc_table3.plot.bar(ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("fraction")
    plt.savefig(filename, bbox_extra_artists=(lgd,), bbox_inches='tight')


def main(args=None):

    """
    The structure of the log file is as follows:
    --------------------------------------------

    File    hic/Beaf32_KD_repA_rf.h5
    Pairs considered        102449560
    Min rest. site distance 150
    Max rest. site distance 1500


    Pairs used      14537174        (14.19) (14.62)
    One mate unmapped       3031462 (2.96)  (3.05)
    One mate not unique     29915787        (29.20) (30.09)
    One mate low quality    6471448 (6.32)  (6.51)
    dangling end    2398169 (2.34)  (2.41)
    self ligation (removed) 14263507        (13.92) (14.35)
    One mate not close to rest site 371977  (0.36)  (0.37)
    same fragment (800 bp)  28543891        (27.86) (28.71)
    self circle     2845649 (2.78)  (2.86)
    duplicated pairs        2916144 (2.85)  (2.93)
    Of pairs used:
    inter chromosomal       951976  (6.55)
    short range < 20kb      9337137 (64.23)
    long range      4248061 (29.22)
    inward pairs    2686016 (18.48)
    outward pairs   3709843 (25.52)
    left pairs      3589334 (24.69)
    right pairs     3600005 (24.76)

    """

    args = parse_arguments().parse_args(args)
    params = dict()
    make_sure_path_exists(args.outputFolder)
    for fh in args.logfiles:
        in_log_part = False
        for line in fh.readlines():
            if line.startswith("File"):
                in_log_part = True
            if in_log_part is True:
                if line.strip() == "":
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

    import pandas as pd
    table = pd.DataFrame(params)
    if args.labels and len(args.labels) == len(args.logfiles):
            table['Labels'] = args.labels
            table = table.set_index('Labels')
    else:
        table = table.set_index('File')

    table.to_csv(args.outputFolder + "/table.txt", sep="\t")

    make_figure_pairs_considered(table, args.outputFolder + "/pairs_considered.png")
    make_figure_pairs_used(table, args.outputFolder + "/pairs_used.png")
    make_figure_distance(table, args.outputFolder + "/distance.png")
    make_figure_read_orientation(table, args.outputFolder + "/read_orientation.png")

    save_html(args.outputFolder + "/hicQC.html", table.T.to_html())
