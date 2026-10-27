"""The five QC bar charts of hicPrepareQCreport (and hicQC, hicBuildMatrix,
hicBuildMatrixMicroC and hicQuickQC, which draw the same report).

The C++ binary writes the five tables and hicQC.html, and computes the
DataFrame every chart plots (hicx/qc_report.hpp: the column selections and
divisions of make_figure_* in hicexplorer/hicPrepareQCreport.py). This module
only rebuilds each DataFrame from those values and draws it with the
reference's pandas plotting calls, unchanged.

Data:
    outputFolder, dpi
    charts   {name: {"index", "index_name", "columns", "columns_name",
                     "values" (row major)}} for pairs_sequenced,
             unmappable_and_non_unique, pairs_discarded, distance and
             read_orientation
"""

import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402


def _frame(chart):
    index = pd.Index(chart["index"], name=chart["index_name"], dtype=object)
    columns = pd.Index(chart["columns"], name=chart["columns_name"], dtype=object)
    values = np.array(chart["values"], dtype=np.float64).reshape(len(index), len(columns))
    return pd.DataFrame(values, index=index, columns=columns)


def make_figure_pairs_used(prc_table, filename, dpi):
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


def make_figure_umappable_non_unique_reads(prc_table, filename, dpi):
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


def make_figure_pairs_discarded(prc_table, filename, dpi):
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


def make_figure_distance(prc_table2, filename):
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


def make_figure_read_orientation(prc_table3, filename, dpi):
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
    charts = data["charts"]
    make_figure_pairs_used(_frame(charts["pairs_sequenced"]), folder + "/pairs_sequenced.png", dpi)
    make_figure_umappable_non_unique_reads(_frame(charts["unmappable_and_non_unique"]),
                                           folder + "/unmappable_and_non_unique.png", dpi)
    make_figure_pairs_discarded(_frame(charts["pairs_discarded"]), folder + "/pairs_discarded.png", dpi)
    make_figure_distance(_frame(charts["distance"]), folder + "/distance.png")
    make_figure_read_orientation(_frame(charts["read_orientation"]), folder + "/read_orientation.png", dpi)
