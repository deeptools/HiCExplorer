from __future__ import division
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import numpy as np

from future.utils import iteritems

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from hicmatrix import HiCMatrix as hm
import hicexplorer.utilities
from .utilities import toString
from .utilities import check_chrom_str_bytes
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)
from collections import OrderedDict


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Takes a list of positions in the hic-matrix and '
                                                 'makes a pooled image.')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='Path of the Hi-C matrix to plot.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-out',
                                help='File name to save the image. ',
                                type=argparse.FileType('w'),
                                required=True)

    parserRequired.add_argument('--BED',
                                help='Interactions between regions in this BED file are plotted.',
                                type=argparse.FileType('r'),
                                required=True)

    parserRequired.add_argument('--range',
                                help='Range of contacts that will be considered for plotting the aggregate contacts '
                                'in bp with the format low_range:high_range for example 1000000:20000000. The range '
                                'should start at contacts larger than TAD size to reduce background interactions.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--BED2',
                           help='Optional second BED file. Interactions between regions in first '
                           'and second BED file are plotted.',
                           type=argparse.FileType('r'),
                           required=False)

    parserOpt.add_argument('--numberOfBins',
                           help='Number of  bins to include in the submatrix. The bed regions will be centered between '
                           '- half number of bins and the other half number of bins.',
                           default='51',
                           type=int)

    parserOpt.add_argument('--transform',
                           help='Type of transformation for the matrix. The options are "none",  '
                           '"total-counts", "z-score" or "obs/exp". If total counts are selected, '
                           'then the sub-matrix values are divided by the total counts for normalization. '
                           'If z-score or obs/exp are selected, then H-C matrix is converted into a '
                           'z-score or observed / expected matrix.',
                           choices=['total-counts', 'z-score', 'obs/exp', 'none'],
                           default='none')

    parserOpt.add_argument('--avgType',
                           help='Type of average to compute final matrix. Options are mean and median. Default is median.',
                           choices=['mean', 'median'],
                           default='median')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    parserOut = parser.add_argument_group('Output options')

    parserOut.add_argument('--outFilePrefixMatrix',
                           help='If this option is given, then the values underlying the final matrix will be '
                           'saved to tab-delimited tables (one per chromosome) using the indicated prefix, '
                           'for example TSS_to_TSS_chrX.tab. If clustering is performed, then the values are '
                           'saved including the cluster_id a in TSS_to_TSS_chrX_cluster_1.tab',
                           required=False)

    parserOut.add_argument('--outFileContactPairs',
                           help='The value should be a prefix. If this option is given, then the position '
                           'of the contacts positions are saved as (chrom1, start1, end1, chrom2, start2, end2) '
                           'where chrom_n, start_n, end_n correspond to the pair of positions used to compute '
                           'the submatrix. The data is saved per chromosome and per '
                           'cluster separately (one file each). The format is {prefix}_{chrom}_{cluster_id}.tab. '
                           'If no clusters were computed, then only one file per chromosome is produced.',
                           required=False)

    parserOut.add_argument('--diagnosticHeatmapFile',
                           help='If given, a heatmap file (per chromosome) is saved. Each row in the heatmap contains the'
                           'diagonal of each of the submatrices centered on the bed file. This file is useful to '
                           'get an idea of the values that are used for the aggregate matrix and to determine '
                           'the fraction of sub-matrices that are aggregated that may have an enrichment at the '
                           'center.',
                           type=argparse.FileType('w'),
                           required=False)

    parserClust = parser.add_argument_group('Clustering options')

    parserClust.add_argument('--kmeans',
                             help='Number of clusters to compute. When this '
                             'option is set, the submatrices are split into clusters (per chromosome)'
                             'using the k-means algorithm.',
                             type=int)

    parserClust.add_argument('--hclust',
                             help='Number of clusters to compute (per chromosome). When this '
                             'option is set, then the matrix is split into clusters '
                             'using the hierarchical clustering algorithm, using "ward linkage". '
                             ' --hclust could be very slow if you have '
                             '>1000 submatrices per chromosome. In those cases, you might prefer --kmeans',
                             type=int)

    parserClust.add_argument('--howToCluster',
                             help='Options are "full", "center" and "diagonal". The full clustering is the default and '
                             'takes all values of each submatrix for clustering. "center", takes only a square of '
                             'length 3x3 from each submatrix and uses only  this values for clustering. With the '
                             '"diagonal" option the clustering is only carried out based on the submatrix diagonal '
                             '(representing values at the same distance to each other.)',
                             choices=['full', 'center', 'diagonal'],
                             default='full')

    parserPlot = parser.add_argument_group('Plotting options')

    parserPlot.add_argument('--chromosomes', '-C',
                            help='List of chromosomes to plot.',
                            nargs='+')

    parserPlot.add_argument('--colorMap',
                            help='Color map to use for the heatmap. Available '
                            'values can be seen here: '
                            'http://matplotlib.org/examples/color/colormaps_reference.html',
                            default='RdYlBu_r')

    parserPlot.add_argument('--plotType',
                            help='Plot type.',
                            choices=['2d', '3d'],
                            default='2d')

    parserPlot.add_argument('--vMin',
                            help='Minimum value of the plotted score.',
                            type=float,
                            default=None)

    parserPlot.add_argument('--vMax',
                            help='Maximum value of the plotted score.',
                            type=float,
                            default=None)

    #  If set, then the bbox=tight option is disable. Used for automatic testing
    parserPlot.add_argument('--disable_bbox_tight',
                            help=argparse.SUPPRESS,
                            action='store_true')
    parserOpt.add_argument('--dpi',
                           help='Optional parameter: Resolution for the image in case the'
                           'ouput is a raster graphics image (e.g png, jpg)',
                           type=int,
                           default=300)
    return parser


def read_bed_per_chrom(fh):
    """
    Reads the given BED file returning
    a dictionary that contains, per each chromosome
    a list of start, end
    """
    interval = {}
    for line in fh:
        if line[0] == "#":
            continue
        fields = line.strip().split()
        if fields[0] not in interval:
            interval[fields[0]] = []

        interval[fields[0]].append((int(fields[1]), int(fields[2])))

    return interval


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
    data = np.asarray(data)
    median = np.median(data)
    b_value = 1.4826  # value set for a normal distribution
    mad = b_value * np.median(np.abs(data))
    if mad == 0:
        return None
    deviation = abs(data - median) / mad
    # any is used to get, per row, if any of the values is True.
    # as a row will be considered as an outlier if any of the elements
    # is an outlier
    outliers = (deviation > max_deviation).any(axis=1)

    return outliers


def cluster_matrices(submatrices_dict, k, method='kmeans', how='full'):
    """
    clusters the submatrices per chromosome


    Parameters
    ----------
    submatrices_dict key: chrom name, values, a list of submatrices
    k number of clusters
    method either kmeans or hierarchical
    how how to cluster. Options are 'full', 'center' and 'diagonal'. More info in the argparse options

    Returns
    -------

    indices dict key: chrom_name, value: list of list, with one list per cluster with the ids of the submatrices
                 that belong to that list
    """
    clustered_dict = {}
    for chrom in submatrices_dict:
        log.info("Length of entry: {}".format(len(submatrices_dict[chrom])))
        submat_vectors = []
        shape = submatrices_dict[chrom][0].shape
        center_bin = (shape[0] + 1) // 2
        for submatrix in submatrices_dict[chrom]:
            if how == 'diagonal':
                # take from each matrix the diagonal
                submat_vectors.append(submatrix.diagonal())
            elif how == 'center':
                # take the mean of a  smaller submatrix of 3 x 3 centered on the submatrix
                submat_vectors.append(
                    submatrix[center_bin - 2:center_bin + 1, center_bin - 2:center_bin + 1].reshape((1, 9)).mean())
            else:
                # Transform list of submatrices in an array of shape:
                # shape = (num_submatrices, submatrix.shape[0] * submatrix.shape[1]
                # In other words, each submatrix is converted into a row of the matrix
                submat_vectors.append(submatrix.reshape((1, shape[0] * shape[1])))

        matrix = np.vstack(submat_vectors)
        if how == 'diagonal':
            assert matrix.shape == (len(submatrices_dict[chrom]), shape[0])
        elif how == 'center':
            assert matrix.shape == (len(submatrices_dict[chrom]), 1)
        else:
            assert matrix.shape == (len(submatrices_dict[chrom]), shape[0] * shape[1])

        # remove outliers
        out_ind = get_outlier_indices(matrix, max_deviation=2)
        if out_ind is not None and len(np.flatnonzero(out_ind)) > 0:
            log.info("Outliers detected in chrom: {}. Number of outliers: {}".
                     format(chrom, len(np.flatnonzero(out_ind))))

            # keep in matrix all indices that are not outliers
            matrix = matrix[np.logical_not(out_ind), :]

        if np.any(np.isnan(matrix)):
            # replace nans for 0 otherwise kmeans produces a weird behaviour
            log.warning("For clustering nan values have to be replaced by zeros.")
            matrix[np.isnan(matrix)] = 0

        if method == 'kmeans':
            from scipy.cluster.vq import vq, kmeans

            centroids, _ = kmeans(matrix, k)
            # order the centroids in an attempt to
            # get the same cluster order
            cluster_labels, _ = vq(matrix, centroids)

        if method == 'hierarchical':
            # normally too slow for large data sets
            from scipy.cluster.hierarchy import fcluster, linkage
            Z = linkage(matrix, method='ward', metric='euclidean')
            cluster_labels = fcluster(Z, k, criterion='maxclust')
            # hierarchical clustering labels from 1 .. k
            # while k-means labels 0 .. k -1
            # Thus, for consistency, we subtract 1
            cluster_labels -= 1

        # sort clusters
        clustered_dict[chrom] = []
        for cluster in range(k):
            cluster_ids = np.flatnonzero(cluster_labels == cluster)
            clustered_dict[chrom].append(cluster_ids)

    return clustered_dict


def plot_aggregated_contacts(chrom_matrix, chrom_contact_position, cluster_ids, num_clusters, M_half, args):

    num_chromosomes = len(chrom_matrix)

    fig = plt.figure(figsize=(5.5 * num_chromosomes, 5.5 * num_clusters + 0.5))
    gs = gridspec.GridSpec(num_clusters + 1, num_chromosomes,
                           width_ratios=[10] * len(chrom_matrix),
                           height_ratios=[10] * num_clusters + [0.6])

    gs.update(wspace=0.01, hspace=0.2)
    chrom_avg = {}
    chrom_cluster_len = {}
    for idx, chrom in enumerate(chrom_matrix):
        chrom_avg[chrom] = []
        chrom_cluster_len[chrom] = []
        for cluster_number, cluster_indices in enumerate(cluster_ids[chrom]):
            # compute median values
            if num_clusters == 1:
                # this means no clustering
                submatrices = np.array(chrom_matrix[chrom])
            else:
                submatrices = np.array([chrom_matrix[chrom][x] for x in cluster_indices])

            chrom_cluster_len[chrom].append(len(cluster_ids))

            if args.avgType == 'median':
                _median = np.median(submatrices, axis=0)
                if _median.sum() == 0 or np.isnan(_median.sum()):
                    # test if the mean matrix is not zero
                    if np.mean(submatrices, axis=0).sum() != 0:
                        log.info("The median of the matrices is zero. Consider using "
                                 "the mean instead.")
                    else:
                        log.info("Apparently no matrices could be computed. All are "
                                 "zeros or nans.")
                chrom_avg[chrom].append(_median)
            else:
                chrom_avg[chrom].append(np.mean(submatrices, axis=0))

            log.info("Mean aggregate matrix values: {}".format(chrom_avg[chrom][cluster_number].mean()))

    vmin, vmax = (args.vMin, args.vMax)
    cmap = cm.get_cmap(args.colorMap)

    log.debug("vmax: {}, vmin: {}".format(vmax, vmin))
    for idx, chrom in enumerate(chrom_matrix):
        for cluster_number, cluster_indices in enumerate(cluster_ids[chrom]):
            log.info("total pairs considered for {}, cluster_{}: {}".format(chrom, cluster_number + 1,
                                                                            len(cluster_indices)))
            try:
                chrom_avg[chrom][cluster_number].shape[0]
            except IndexError:
                continue
            if chrom_avg[chrom][cluster_number].shape[0] == 0:
                log.debug("matrix for chrom {} is empty".format(chrom))
                continue
            if num_clusters == 1:
                title = chrom
            else:
                title = "{} cluster_{}".format(chrom, cluster_number + 1)
            if args.plotType == '2d':
                ax = plt.subplot(gs[cluster_number, idx])

                ax.set_title(title)
                img = ax.imshow(chrom_avg[chrom][cluster_number], aspect='equal',
                                interpolation='nearest', vmax=vmax, vmin=vmin,
                                cmap=cmap,
                                extent=[-M_half, M_half + 1, -M_half, M_half + 1])
            else:
                from mpl_toolkits.mplot3d import Axes3D
                # Axes3D is required for projection='3d' to work
                # but since is imported but not used, flake8 will complain
                # thus I add this dummy variable to avoid the error
                Axes3D(fig)
                ax = plt.subplot(gs[cluster_number, idx], projection='3d')
                ax.set_aspect('equal')
                ax.margins(0)
                X, Y = np.meshgrid(range(-M_half, M_half + 1),
                                   range(-M_half, M_half + 1))
                Z = chrom_avg[chrom][cluster_number].copy()

                img = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                                      vmax=vmax, vmin=vmin, edgecolor='none')

                ax.set_zticklabels([])
                if vmax is not None and vmax is not None:
                    ax.set_zlim(vmin, vmax)

            if args.outFilePrefixMatrix:
                # save aggregate matrix values
                if num_clusters == 1:
                    output_matrix_name = "{file}_{chrom}.tab".format(file=args.outFilePrefixMatrix, chrom=chrom)
                else:
                    output_matrix_name = "{file}_{chrom}_cluster_{id}.tab".format(file=args.outFilePrefixMatrix,
                                                                                  chrom=chrom, id=cluster_number + 1)
                np.savetxt(output_matrix_name, chrom_avg[chrom][cluster_number], '%0.5f', delimiter='\t')

        cbar_x = plt.subplot(gs[-1, idx])
        fig.colorbar(img, cax=cbar_x, orientation='horizontal')

    if args.disable_bbox_tight:
        plt.savefig(args.outFileName.name, dpi=args.dpi)
    else:
        plt.savefig(args.outFileName.name, dpi=args.dpi, bbox_inches='tight')

    plt.close()


def plot_diagnostic_heatmaps(chrom_diagonals, cluster_ids, M_half, args):

    num_chromosomes = len(chrom_diagonals)
    vmax_heat = args.vMax
    if vmax_heat is not None:
        vmax_heat *= 5

    vmin_heat = args.vMin
    if vmin_heat is not None:
        vmin_heat *= 5
    else:
        vmin_heat = 0

    num_plots = len(chrom_diagonals)
    fig = plt.figure(figsize=(num_plots * 4, 20))

    gs0 = gridspec.GridSpec(2, num_plots + 1, width_ratios=[10] * num_plots + [0.5], height_ratios=[1, 5],
                            wspace=0.1, hspace=0.1)

    gs_list = []
    for idx, (chrom_name, values) in enumerate(iteritems(chrom_diagonals)):
        try:
            heatmap = np.asarray(np.vstack(values))
        except ValueError:
            log.error("Error computing diagnostic heatmap for chrom: {}".format(chrom_name))
            continue

        # get size of each cluster for the given chrom
        clust_len = [(len(v)) for v in cluster_ids[chrom_name]]

        # prepare layout
        gs_list.append(gridspec.GridSpecFromSubplotSpec(len(clust_len), 1,
                                                        subplot_spec=gs0[1, idx],
                                                        height_ratios=clust_len,
                                                        hspace=0.03))
        summary_plot_ax = plt.subplot(gs0[0, idx])
        summary_plot_ax.set_title(chrom_name)

        for cluster_number, cluster_indices in enumerate(cluster_ids[chrom_name]):
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

    file_name = args.diagnosticHeatmapFile.name
    log.info('Heatmap file saved under: {}'.format(file_name))
    plt.savefig(file_name, dpi=args.dpi, bbox_inches='tight')
    plt.close()


def main(args=None):
    args = parse_arguments().parse_args(args)

    ma = hm.hiCMatrix(args.matrix)
    ma.maskBins(ma.nan_bins)
    ma.matrix.data[np.isnan(ma.matrix.data)] = 0

    bin_size = ma.getBinSize()
    ma.maskBins(ma.nan_bins)
    ma.matrix.data = ma.matrix.data
    new_intervals = hicexplorer.utilities.enlarge_bins(ma.cut_intervals)
    ma.setCutIntervals(new_intervals)
    min_dist, max_dist = args.range.split(":")

    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)
    chrom_sizes = ma.get_chromosome_sizes()
    chrom_list = chrom_sizes.keys()
    log.info("checking range {}-{}".format(min_dist, max_dist))
    min_dist = int(min_dist)
    max_dist = int(max_dist)
    assert min_dist < max_dist, "Error lower range larger than upper range"

    if args.transform == "z-score":
        # use zscore matrix
        log.info("Computing z-score matrix. This may take a while.\n")
        ma.convert_to_zscore_matrix(maxdepth=max_dist * 2.5, perchr=True)
    elif args.transform == "obs/exp":
        # use zscore matrix
        log.info("Computing observed vs. expected matrix. This may take a while.\n")
        ma.convert_to_obs_exp_matrix(maxdepth=max_dist * 2.5, perchr=True)

    min_dist_in_bins = int(min_dist) // bin_size
    max_dist_in_bins = int(max_dist) // bin_size

    # read and sort bedgraph.
    bed_intervals = read_bed_per_chrom(args.BED)
    if args.BED2:
        bed_intervals2 = read_bed_per_chrom(args.BED2)
    else:
        bed_intervals2 = bed_intervals

    M = args.numberOfBins if args.numberOfBins % 2 == 1 else args.numberOfBins + 1
    M_half = int((M - 1) // 2)
    # make a new matrix for each chromosome.
    chrom_matrix = OrderedDict()
    chrom_total = {}
    chrom_diagonals = OrderedDict()
    chrom_contact_position = {}
    seen = {}

    center_values = {}

    chrom_list = check_chrom_str_bytes(bed_intervals, chrom_list)

    for chrom in chrom_list:
        if chrom not in bed_intervals:
            continue

        chrom_matrix[chrom] = []
        chrom_total[chrom] = 1
        chrom_diagonals[chrom] = []
        chrom_contact_position[chrom] = []
        center_values[chrom] = []
        seen[chrom] = set()
        over_1_5 = 0
        empty_mat = 0
        chrom_bin_range = ma.getChrBinRange(toString(chrom))

        log.info("processing {}".format(chrom))

        counter = 0
        for start, end in bed_intervals[chrom]:
            # check all other regions that may interact with the
            # current interval at the given depth range
            if end > chrom_sizes[chrom]:
                continue
            bin_id = ma.getRegionBinRange(toString(chrom), start, end)
            if bin_id is None:
                continue
            else:
                bin_id = bin_id[0]

            for start2, end2 in bed_intervals2[chrom]:
                counter += 1
                if counter % 50000 == 0:
                    log.info("Number of contacts considered: {:,}".format(counter))

                if end2 > chrom_sizes[chrom]:
                    continue
                bin_id2 = ma.getRegionBinRange(toString(chrom), start2, end2)
                if bin_id2 is None:
                    continue
                else:
                    bin_id2 = bin_id2[0]
                if bin_id2 in seen[chrom]:
                    continue
                if bin_id == bin_id2:
                    continue
                if min_dist_in_bins <= abs(bin_id2 - bin_id) <= max_dist_in_bins:
                    idx1, idx2 = sorted([bin_id, bin_id2])
                    if (idx1, idx2) in seen[chrom]:
                        continue
                    seen[chrom].add((idx1, idx2))
                    if idx1 - M_half < chrom_bin_range[0] or idx2 + 1 + M_half > chrom_bin_range[1]:
                        continue
                    try:
                        mat_to_append = ma.matrix[idx1 - M_half:idx1 + M_half + 1, :][:, idx2 - M_half:idx2 + M_half + 1].todense().astype(float)
                    except IndexError:
                        log.info("index error for {} {}".format(idx1, idx2))
                        continue
                    counter += 1
                    if counter % 1000 == 0:
                        log.info("Number of contacts within range computed: {:,}".format(counter))
                    if mat_to_append.sum() == 0:
                        empty_mat += 1
                        continue
                    # to account for the fact that submatrices
                    # close to the diagonal have more counts thatn
                    # submatrices far from the diagonal
                    # the submatrices values are normalized using the
                    # total submatrix sum.

                    if args.transform == 'total_counts' and mat_to_append.sum() > 0:
                        mat_to_append = mat_to_append / mat_to_append.sum()

                    chrom_total[chrom] += 1
                    chrom_matrix[chrom].append(mat_to_append)
                    chrom_diagonals[chrom].append(mat_to_append.diagonal())
                    center_values[chrom].append(ma.matrix[idx1, idx2])
                    chrom_contact_position[chrom].append((start, end, start2, end2))
                    if ma.matrix[idx1, idx2] > 1.5:
                        over_1_5 += 1

        if len(chrom_matrix[chrom]) == 0:
            log.warn("No valid submatrices were found for chrom: {}".format(chrom))
            chrom_matrix.pop(chrom, None)

        log.info("Number of matrices with ratio over 1.5 at center {}, fraction w.r.t. non empty submatrices: ({:.2f})".
                 format(over_1_5, float(over_1_5) / len(chrom_matrix[chrom])))

        log.info("Number of discarded empty submatrices  {} ({:.2f})".
                 format(empty_mat, float(empty_mat) / counter))

    if args.kmeans is not None:
        cluster_ids = cluster_matrices(chrom_matrix, args.kmeans, method='kmeans', how=args.howToCluster)
        num_clusters = args.kmeans
    elif args.hclust is not None:
        log.info("Performing hierarchical clustering."
                 "Please note that it might be very slow for large datasets.\n")
        cluster_ids = cluster_matrices(chrom_matrix, args.hclust, method='hierarchical',
                                       how=args.howToCluster)
        num_clusters = args.hclust
    else:
        # make a 'fake' clustering to generalize the plotting of the submatrices
        cluster_ids = {}
        num_clusters = 1
        for chrom in chrom_list:
            cluster_ids[chrom] = [range(len(chrom_matrix[chrom]))]

    plot_aggregated_contacts(chrom_matrix, chrom_contact_position, cluster_ids, num_clusters, M_half, args)

    if args.outFileContactPairs:
        for idx, chrom in enumerate(chrom_matrix):

            for cluster_number, cluster_indices in enumerate(cluster_ids[chrom]):
                center_values_to_order = np.array(center_values[chrom])[cluster_indices]
                center_values_order = np.argsort(center_values_to_order)[::-1]

                output_name = "{file}_{chrom}_cluster_{id}.tab".format(file=args.outFileContactPairs,
                                                                       chrom=chrom, id=cluster_number + 1)
                with open(output_name, 'w') as fh:
                    for cl_idx in center_values_order:
                        value = center_values_to_order[cl_idx]
                        start, end, start2, end2 = chrom_contact_position[chrom][cl_idx]
                        fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, chrom, start2, end2, value))

    # plot the diagonals
    # the diagonals plot is useful to see individual cases and if they had a contact in the center
    if args.diagnosticHeatmapFile:
        plot_diagnostic_heatmaps(chrom_diagonals, cluster_ids, M_half, args)
