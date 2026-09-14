"""The figures and hicQC.html of hicPrepareQCreport (and hicQC, hicBuildMatrix
and hicQuickQC, which draw the same report).

The C++ binary writes the five tables. This module reads QC_table.txt back
for the charts and the other four tables for the HTML; the calls are those of
hicexplorer/hicPrepareQCreport.py (make_figure_* and save_html), unchanged,
and qc_template.html is a copy of hicexplorer/qc_template.html. The tables
round trip through to_csv and read_csv with the dtypes they were written
with: integer counts, float64 fractions, text labels.

Data:
    outputFolder, dpi
"""

import os

import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402


def save_html(filename, unmap_table, discard_table, distance_table, orientation_table, all_table):
    root = os.path.dirname(os.path.abspath(__file__))

    html = open(os.path.join(root, "qc_template.html"), "r")
    html_content = html.read()
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

    html_content = html_content.replace("%%TABLE%%", all_table.style.to_html(classes='df'))
    with open(filename, 'w') as fh:
        fh.write(html_content)
    html.close()


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


def make_figure_pairs_discarded(table, filename, dpi):
    column_names_prefix = ['One mate not close to rest site', 'dangling end', 'duplicated pairs',
                           'same fragment', 'self circle',
                           'self ligation (removed)']
    column_names = []

    column_names_table = list(table.columns)
    for name in column_names_prefix:
        if name in column_names_table:
            column_names.append(name)
        else:
            for name_ in column_names_table:
                if name in name_:
                    column_names.append(name_)

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


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    folder = data["outputFolder"]
    dpi = data["dpi"]

    def read(name):
        return pd.read_csv(os.path.join(folder, name), sep="\t", index_col=0)

    table = read("QC_table.txt")
    make_figure_pairs_used(table, folder + "/pairs_sequenced.png", dpi)
    make_figure_umappable_non_unique_reads(table, folder + "/unmappable_and_non_unique.png", dpi)
    make_figure_pairs_discarded(table, folder + "/pairs_discarded.png", dpi)
    make_figure_distance(table, folder + "/distance.png")
    make_figure_read_orientation(table, folder + "/read_orientation.png", dpi)

    save_html(folder + "/hicQC.html", read("unmapable_table.txt"), read("discarded_table.txt"),
              read("distance_table.txt"), read("read_orientation_table.txt"), table)
