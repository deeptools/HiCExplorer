"""hicAggregateContacts, implemented in C++.

main() runs the C++ executable hicAggregateContacts and parse_arguments() is
built from its --help-json (hicexplorer._cpp). The public functions keep the
name, the signature and the behaviour they had in HiCExplorer 3.7.6; their
bodies call hicexplorer._hicx.hicAggregateContacts
(cpp/python/src/api/hicAggregateContacts.cpp). plot_aggregated_contacts and
plot_diagnostic_heatmaps are matplotlib drawing code and stay Python
(cpp/scripts/public_api_allowlist.tsv); the averages, files and sorted rows
they draw and write come from hicexplorer._hicx.hicAggregateContacts.
"""
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colormaps as cm
from mpl_toolkits.mplot3d import Axes3D

import logging
log = logging.getLogger(__name__)
from collections import OrderedDict

from hicexplorer._hicx import hicAggregateContacts as _impl


def parse_arguments(args=None):
    """The argparse parser of hicAggregateContacts, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAggregateContacts')


def read_bed_per_chrom(fh, chrom_list, pConsiderStrandDirection=None):
    """
    Reads the given BED file returning
    a dictionary that contains, per each chromosome
    a list of start, end
    """
    return _impl.read_bed_per_chrom(fh, chrom_list, pConsiderStrandDirection)


def aggregate_contacts(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range=None, transform=None, mode='', pConsiderStrandDirection=None):
    """
    To aggregate the contacts of desired submatrices.
    """
    return _impl.aggregate_contacts(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range, transform, mode,
                                    pConsiderStrandDirection)


def aggregate_contacts_per_row(bed1, bed2, agg_info, ma, chrom_list, M_half, largeRegionsOperation, range=None, transform=None, mode='', perChr=False, pConsiderStrandDirection=None):
    """
    To aggregate the contacts of the desired submatrices , if row-wise.
    """
    return _impl.aggregate_contacts_per_row(bed1, bed2, agg_info, ma, chrom_list, M_half, largeRegionsOperation,
                                            range, transform, mode, perChr, pConsiderStrandDirection)


def count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range=None, transform=None):
    """
    To count the number of contacts for a given pair of intervals
    """
    return _impl.count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range,
                                transform)


def get_outlier_indices(data, max_deviation=200):
    """
    The method is based on the median absolute deviation. See
    Boris Iglewicz and David Hoaglin (1993),
    "Volume 16: How to Detect and Handle Outliers",
    The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

    The max_deviation=200 is like selecting a z-score
    larger than 200, just that it is based on the median
    and the median absolute deviation instead of the
    mean and the standard deviation.

    Returns
    ------

    Boolean of len (data.shape[0])
    """
    return _impl.get_outlier_indices(data, max_deviation)


def compute_clusters(updated_info, k, method="kmeans", how='full', max_deviation=2, keep_outlier=False):
    return _impl.compute_clusters(updated_info, k, method, how, max_deviation, keep_outlier)


def cluster_matrices(agg_info, k, method='kmeans', how='full', perChr=False, max_deviation=2, keep_outlier=False):
    """
    clusters the submatrices .

    Parameters
    ----------
    agg_info["agg_matrix"] key: chrom name, values, a list of submatrices
    k number of clusters
    method either kmeans, hierarchical or spectral
    how how to cluster. Options are 'full', 'center' and 'diagonal'. More info in the argparse options

    Returns
    -------

    indices dict key: chrom_name, value: list of list, with one list per cluster with the ids of the submatrices
                 that belong to that list
    """
    return _impl.cluster_matrices(agg_info, k, method, how, perChr, max_deviation, keep_outlier)


def compute_avg(submatrices, operationType):
    return _impl.compute_avg(submatrices, operationType)


def plot_aggregated_contacts(clustered_info, num_clusters, M_half, args):

    num_figs = len(clustered_info.keys())

    fig = plt.figure(figsize=(5.5 * num_figs, 5.5 * num_clusters + 0.5))
    gs = gridspec.GridSpec(num_clusters + 1, num_figs,
                           width_ratios=[10] * num_figs,
                           height_ratios=[10] * num_clusters + [0.6])

    gs.update(wspace=0.01, hspace=0.2)
    vmin, vmax = (args.vMin, args.vMax)
    cmap = cm.get_cmap(args.colorMap)

    log.debug("vmax: {}, vmin: {}".format(vmax, vmin))
    chrom_avg = OrderedDict()
    for idx, (chrom1, v1) in enumerate(clustered_info.items()):
        assert (v1 != {})
        if chrom1 not in chrom_avg.keys():
            chrom_avg[chrom1] = []

        for cluster_number, cluster_indices in enumerate(clustered_info[chrom1]["clustered_dict"]):
            # the average of the cluster's submatrices and its log records (C++)
            average, empty = _impl._cluster_average(clustered_info[chrom1], cluster_indices,
                                                    args.operationType, cluster_number)
            chrom_avg[chrom1].append(average)
            if empty:
                continue
            title = "cluster_{}".format(cluster_number + 1)
            if args.plotType == '2d':
                ax = plt.subplot(gs[cluster_number, idx])
                if num_clusters != 1:
                    ax.set_title(title)
                if (cluster_number + 1) == num_clusters:
                    ax.set_xlabel("{}".format(chrom1))
                img = ax.imshow(chrom_avg[chrom1][cluster_number], aspect='equal',
                                interpolation='nearest', vmax=vmax, vmin=vmin,
                                cmap=cmap,
                                extent=[-M_half, M_half + 1, -M_half, M_half + 1])
            else:
                # Axes3D is required for projection='3d' to work
                # but since is imported but not used, flake8 will complain
                # thus I add this dummy variable to avoid the error
                Axes3D(fig)
                ax = plt.subplot(gs[cluster_number, idx], projection='3d')
                # ax.set_aspect('equal')
                ax.margins(0)
                X, Y = np.meshgrid(range(-M_half, M_half + 1),
                                   range(-M_half, M_half + 1))
                Z = chrom_avg[chrom1][cluster_number].copy()
                img = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                                      vmax=vmax, vmin=vmin, edgecolor='none')

                ax.set_zticklabels([])
                if vmax is not None and vmax is not None:
                    ax.set_zlim(vmin, vmax)

            # --outFilePrefixMatrix and --outFileContactPairs (C++)
            _impl._cluster_outputs(clustered_info[chrom1], chrom1, cluster_number, cluster_indices,
                                   num_clusters, chrom_avg[chrom1][cluster_number], args)

        cbar_x = plt.subplot(gs[-1, idx])
        fig.colorbar(img, cax=cbar_x, orientation='horizontal')

    if args.disable_bbox_tight:
        plt.savefig(args.outFileName.name, dpi=args.dpi)
    else:
        plt.savefig(args.outFileName.name, dpi=args.dpi, bbox_inches='tight')

    plt.close()


def plot_diagnostic_heatmaps(clustered_info, M_half, args):
    num_chromosomes = len(clustered_info.keys())
    vmax_heat = args.vMax
    if vmax_heat is not None:
        vmax_heat *= 5

    vmin_heat = args.vMin
    if vmin_heat is not None:
        vmin_heat *= 5
    else:
        vmin_heat = 0

    num_plots = len(clustered_info.keys())
    fig = plt.figure(figsize=(num_plots * 4, 20))

    gs0 = gridspec.GridSpec(2, num_plots + 1, width_ratios=[10] * num_plots + [0.5], height_ratios=[1, 5],
                            wspace=0.1, hspace=0.1)

    gs_list = []
    for idx, (chrom_name, values) in enumerate(clustered_info.items()):  # One more for loop should be added as it is or this container needs to be updated too. I would go for the latter
        # the diagonals of every cluster, sorted by their centre value, and their
        # column means (C++); None when they do not stack
        diagnostic = _impl._diagnostic_data(clustered_info[chrom_name], chrom_name, M_half)
        if diagnostic is None:
            continue

        # get size of each cluster for the given chrom
        clust_len, cluster_rows = diagnostic

        # prepare layout
        gs_list.append(gridspec.GridSpecFromSubplotSpec(len(clust_len), 1,
                                                        subplot_spec=gs0[1, idx],
                                                        height_ratios=clust_len,
                                                        hspace=0.03))
        summary_plot_ax = plt.subplot(gs0[0, idx])
        summary_plot_ax.set_title(chrom_name)

        for cluster_number, cluster_indices in enumerate(clustered_info[chrom_name]["clustered_dict"]):
            # sorted by the value at the center of the rows
            heatmap_to_plot, y_values, x_values, order_length = cluster_rows[cluster_number]

            # add line to summary plot ax
            cluster_label = "cluster_{}".format(cluster_number + 1)
            summary_plot_ax.plot(x_values, y_values, label=cluster_label)
            ax = plt.subplot(gs_list[-1][cluster_number, 0])
            ax.set_yticks([])
            if num_chromosomes > 1:
                ax.set_ylabel(cluster_label)

            if cluster_number < num_chromosomes - 1:
                ax.set_xticks([])

            heat_fig = ax.imshow(heatmap_to_plot, aspect='auto',
                                 interpolation='nearest',
                                 cmap=cm.get_cmap(args.colorMap),
                                 origin='upper', vmax=vmax_heat, vmin=vmin_heat,
                                 extent=[-M_half, M_half + 1,
                                         0, order_length])

        summary_plot_ax.legend(ncol=1, frameon=False, markerscale=0.5)

    cbar_x = plt.subplot(gs0[1, -1])
    fig.colorbar(heat_fig, cax=cbar_x, orientation='vertical')

    file_name = args.diagnosticHeatmapFile.name
    log.info('Heatmap file saved under: {}'.format(file_name))
    plt.savefig(file_name, dpi=args.dpi, bbox_inches='tight')
    plt.close()


def main(args=None):
    """Run the C++ implementation of hicAggregateContacts with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAggregateContacts')(args)
