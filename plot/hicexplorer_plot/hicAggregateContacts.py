"""The aggregate and diagnostic heatmap figures of hicAggregateContacts.

The drawing calls are those of hicexplorer/hicAggregateContacts.py
plot_aggregated_contacts and plot_diagnostic_heatmaps, unchanged, and the
rcParams line of main(). The C++ binary clusters the submatrices, computes
the aggregate of every cluster (compute_avg) and writes the tables the
reference writes from inside plot_aggregated_contacts; this module receives
the aggregates, the cluster member indices and, for the heatmap, the
diagonal of every kept submatrix.

Data:
    outFileName          --outFileName
    diagnosticHeatmapFile  --diagnosticHeatmapFile or null
    dpi, vMin, vMax, colorMap, plotType, disable_bbox_tight   the options
    M_half               (numberOfBins - 1) // 2 as the reference computes it
    num_clusters         the requested cluster count (1 without clustering)
    chroms              [{"name", "clusters": [{"average": rows, "indices": ids}],
                           "diagonal": rows}]
"""

from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
from matplotlib import colormaps as cm  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402
import numpy as np  # noqa: E402


class _Args(object):
    pass


def _rebuild(data):
    args = _Args()
    args.dpi = data["dpi"]
    args.vMin = data["vMin"]
    args.vMax = data["vMax"]
    args.colorMap = data["colorMap"]
    args.plotType = data["plotType"]
    args.disable_bbox_tight = data["disable_bbox_tight"]
    args.outFileName = data["outFileName"]
    args.diagnosticHeatmapFile = data["diagnosticHeatmapFile"]

    # Every submatrix is todense().astype(float) in the reference, so the
    # aggregates and the diagonals are float64 whatever the matrix dtype.
    clustered_info = OrderedDict()
    chrom_avg = OrderedDict()
    for chrom in data["chroms"]:
        clustered_info[chrom["name"]] = {
            "clustered_dict": [np.array(c["indices"], dtype=np.int64) for c in chrom["clusters"]],
            "diagonal": [np.array(row, dtype=np.float64) for row in chrom["diagonal"]],
        }
        chrom_avg[chrom["name"]] = [np.array(c["average"], dtype=np.float64)
                                    for c in chrom["clusters"]]
    return args, clustered_info, chrom_avg


def plot_aggregated_contacts(clustered_info, chrom_averages, num_clusters, M_half, args):

    num_figs = len(clustered_info.keys())

    fig = plt.figure(figsize=(5.5 * num_figs, 5.5 * num_clusters + 0.5))
    gs = gridspec.GridSpec(num_clusters + 1, num_figs,
                           width_ratios=[10] * num_figs,
                           height_ratios=[10] * num_clusters + [0.6])

    gs.update(wspace=0.01, hspace=0.2)
    vmin, vmax = (args.vMin, args.vMax)
    cmap = cm.get_cmap(args.colorMap)

    chrom_avg = OrderedDict()
    for idx, (chrom1, v1) in enumerate(clustered_info.items()):
        assert (v1 != {})
        if chrom1 not in chrom_avg.keys():
            chrom_avg[chrom1] = []

        for cluster_number, cluster_indices in enumerate(clustered_info[chrom1]["clustered_dict"]):
            chrom_avg[chrom1].append(chrom_averages[chrom1][cluster_number])

            if chrom_avg[chrom1][cluster_number].shape[0] == 0:
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
                Axes3D(fig)
                ax = plt.subplot(gs[cluster_number, idx], projection='3d')
                ax.margins(0)
                X, Y = np.meshgrid(range(-M_half, M_half + 1),
                                   range(-M_half, M_half + 1))
                Z = chrom_avg[chrom1][cluster_number].copy()
                img = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                                      vmax=vmax, vmin=vmin, edgecolor='none')

                ax.set_zticklabels([])
                if vmax is not None and vmax is not None:
                    ax.set_zlim(vmin, vmax)

        cbar_x = plt.subplot(gs[-1, idx])
        fig.colorbar(img, cax=cbar_x, orientation='horizontal')

    if args.disable_bbox_tight:
        plt.savefig(args.outFileName, dpi=args.dpi)
    else:
        plt.savefig(args.outFileName, dpi=args.dpi, bbox_inches='tight')

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
    for idx, (chrom_name, values) in enumerate(clustered_info.items()):
        try:
            heatmap = np.asarray(np.vstack(clustered_info[chrom_name]['diagonal']))
        except ValueError:
            continue

        # get size of each cluster for the given chrom
        clust_len = [(len(v)) for v in clustered_info[chrom_name]["clustered_dict"]]

        # prepare layout
        gs_list.append(gridspec.GridSpecFromSubplotSpec(len(clust_len), 1,
                                                        subplot_spec=gs0[1, idx],
                                                        height_ratios=clust_len,
                                                        hspace=0.03))
        summary_plot_ax = plt.subplot(gs0[0, idx])
        summary_plot_ax.set_title(chrom_name)

        for cluster_number, cluster_indices in enumerate(clustered_info[chrom_name]["clustered_dict"]):
            # sort by the value at the center of the rows
            heatmap_to_plot = heatmap[cluster_indices, :]

            order = np.argsort(heatmap_to_plot[:, M_half])[::-1]
            heatmap_to_plot = heatmap_to_plot[order, :]

            # add line to summary plot ax
            y_values = heatmap_to_plot.mean(axis=0)
            x_values = np.arange(len(y_values)) - M_half
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
                                         0, len(order)])

        summary_plot_ax.legend(ncol=1, frameon=False, markerscale=0.5)

    cbar_x = plt.subplot(gs0[1, -1])
    fig.colorbar(heat_fig, cax=cbar_x, orientation='vertical')

    plt.savefig(args.diagnosticHeatmapFile, dpi=args.dpi, bbox_inches='tight')
    plt.close()


def draw(data):
    matplotlib.rcParams['pdf.fonttype'] = 42
    args, clustered_info, chrom_avg = _rebuild(data)
    plot_aggregated_contacts(clustered_info, chrom_avg, data["num_clusters"], data["M_half"], args)
    if args.diagnosticHeatmapFile:
        plot_diagnostic_heatmaps(clustered_info, data["M_half"], args)
