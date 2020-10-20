import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
# from scipy.cluster.vq import vq, kmeans
# from scipy.cluster.hierarchy import fcluster, linkage
import sklearn.cluster import sk_clust
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
                                     description='Takes a list of positions in the Hi-C matrix and '
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

    parserRequired.add_argument('--mode',
                                choices=['inter-chr', 'intra-chr', 'all'],
                                required=True)



    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--range',
                            help='Range of contacts that will be considered for plotting the aggregate contacts '
                            'in bp with the format low_range:high_range for example 1000000:20000000. The range '
                            'should start at contacts larger than TAD size to reduce background interactions. '
                            'This will be ignored if inter-chromosomal contacts are of interest.',
                            default=None)

    parserOpt.add_argument('--row_wise',
                            help='If given,the insteractions between each row of the BED file and its '
                            'corresponding row of the BED2 file are computed. If intera-chromosomal '
                            'contacts are computed, the rows with different chromosomes are ignored. '
                            'If inter-chromosomal, the rows with same chromosomes are ignored. '
                            'It keeps all the rows if `all`.',
                            action='store_true',
                            required=False)

    parserOpt.add_argument('--BED2',
                           help='Optional second BED file. Interactions between regions in first '
                           'and second BED file are plotted.',
                           type=argparse.FileType('r'),
                           required=False)

    parserOpt.add_argument('--numberOfBins',
                           help='Number of  bins to include in the submatrix. The bed regions will be centered between '
                           'half number of bins and the other half number of bins.',
                           default='51',
                           type=int)

    parserOpt.add_argument('--transform',
                           help='Type of transformation for the matrix. The options are "none",  '
                           '"total-counts", "z-score" and "obs/exp". If total counts are selected, '
                           'the sub-matrix values are divided by the total counts for normalization. '
                           'If z-score or obs/exp are selected, the Hi-C matrix is converted into a '
                           'z-score or observed / expected matrix.',
                           choices=['total-counts', 'z-score', 'obs/exp', 'none'],
                           default='none')

    parserOpt.add_argument('--operationType',
                           help='Type of the operation to be applied to the output matrix. Options are mean and median. Default is median.',
                           choices=['sum', 'mean', 'median'],
                           default='median')

    parserOpt.add_argument('--perChr',
                           help='if set, it generates a plot per chromosome. It is only affected if '
                           'intera-chromosomal contacts are of interest.',
                           action='store_true',
                           required=False)

    parserOpt.add_argument('--largeRegionsOperation',
                           help='If a given coordinate in the bed file is larger than '
                           'a bin of the inpt matrix, by default only the first bin is taken into account. '
                           'However there are more posibilities to handel such a case. User can '
                           'ask for the last bin, sum of the bins, mean or median of the bins which cover '
                           'this region.',
                           choices=['first', 'last', 'sum', 'mean', 'median'],
                           default='first')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    parserOut = parser.add_argument_group('Output options')

    parserOut.add_argument('--outFilePrefixMatrix',
                           help='If this option is set, the values underlying the output matrix will be '
                           'saved to tab-delimited tables (one per chromosome) using the indicated prefix, '
                           'for example TSS_to_TSS_chrX.tab. If clustering is performed, the values are '
                           'saved including the cluster_id in the file TSS_to_TSS_chrX_cluster_1.tab.',
                           required=False)

    parserOut.add_argument('--outFileContactPairs',
                           help='Output file prefix. If this option is set, the position '
                           'of the contact positions are saved as (chrom1, start1, end1, chrom2, start2, end2) '
                           'where chrom_n, start_n, end_n correspond to the pair of positions used to compute '
                           'the submatrix. The data is saved per chromosome and per '
                           'cluster separately (one file each). The format is {prefix}_{chrom}_{cluster_id}.tab. '
                           'If no clusters were computed, then only one file per chromosome is produced.',
                           required=False)

    parserOut.add_argument('--diagnosticHeatmapFile',
                           help='If given, a heatmap (per chromosome) is saved. Each row in the heatmap contains the'
                           'diagonal of each of the submatrices centered on the bed file. This file is useful to '
                           'get an idea of the values that are used for the aggregate matrix and to determine '
                           'the fraction of sub-matrices that are aggregated that may have an enrichment at the '
                           'center.',
                           type=argparse.FileType('w'),
                           required=False)


    parserClust = parser.add_argument_group('Clustering options')

    parserClust.add_argument('--kmeans',
                             help='Number of clusters to compute. When this '
                             'option is set, the submatrices are split into clusters (per chromosome) '
                             'using the k-means algorithm.',
                             type=int)

    parserClust.add_argument('--hclust',
                             help='Number of clusters to compute (per chromosome). When this '
                             'option is set, the matrix is split into clusters '
                             'using the hierarchical clustering algorithm, using "ward linkage". '
                             ' --hclust could be very slow if you have '
                             '>1000 submatrices per chromosome. In those cases, you might prefer --kmeans',
                             type=int)

    parserClust.add_argument('--spectral',
                             help='Number of clusters to compute (per chromosome). When this '
                             'option is set, the matrix is split into clusters '
                             'using the spectral clustering algorithm.',
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
                           'output is a raster graphics image (e.g png, jpg).',
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

def aggregate_contacts(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range = None, transform = None,  mode='', perChr = False):
    """

    """

    for k1, v1 in bed1.items():
        for k2, v2 in bed2.items():
            if (mode == 'inter-chr') & (k1 == k2):
                continue
            if (mode == 'intra-chr') & (k1 != k2):
                continue
            print("v1: ", v1)
            print("v2: ", v2)
            for coord1 in v1:
                for coord2 in v2:
                    if (k1 == k2) and (coord1 == coord2):
                        continue
                    interval = [(k1, coord1[0], coord1[1]), (k2, coord2[0], coord2[1])]
                    print(interval)
                    count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, range, transform, perChr)

    # over_1_5 = 0
    #

        # if len(chrom_matrix[chrom]) == 0:
        #     log.warn("No valid submatrices were found for chrom: {}".format(chrom))
        #     chrom_matrix.pop(chrom)
        # else:
        #     log.info("Number of matrices with ratio over 1.5 at center {}, fraction w.r.t. non empty submatrices: ({:.2f})".
        #              format(over_1_5, float(over_1_5) / len(chrom_matrix[chrom])))
        #
        # log.info("Number of discarded empty submatrices  {} ({:.2f})".
        #          format(empty_mat, float(empty_mat) / counter))


def aggregate_contacts_per_row(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range = None, transform = None,  mode='', perChr = False):
    intervals = []
    for line1, line2 in zip(bed1, bed2):
        line1 = line1.strip().split()
        line2 = line2.strip().split()
        if mode == 'inter-chr': # skip the lines with same chrom
            if line1[0] == line2[0]:
                continue
        elif mode == 'intra-chr': # skip the lines with different chrom
            if line1[0]!=line2[0]:
                continue
        intervals.append((line1[0:3], line2[0:3]))

    for interval in intervals:
        count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, range, transform, perChr)


def count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, range = None, transform = None, perChr = False):
    chrom1, start1, end1 = interval[0]
    chrom2, start2, end2 = interval[1]
    if (chrom1 not in agg_info["chrom_coord"]) or (chrom2 not in agg_info["chrom_coord"]):
        return
    if (int(end1) > agg_info["chrom_coord"][chrom1][1]) or (int(end2) > agg_info["chrom_coord"][chrom2][1]): # TODO these intervals may still partially be overlapped, shall we keep them?
        return
    if (int(start1) < agg_info["chrom_coord"][chrom1][0]) or (int(start2) < agg_info["chrom_coord"][chrom2][0]):
        return

    bin_id1 = ma.getRegionBinRange(toString(chrom1), start1, end1)
    bin_id2 = ma.getRegionBinRange(toString(chrom2), start2, end2)
    if bin_id1 == bin_id2:
        return
    if (bin_id1 is None) or (bin_id2 is None):
        return
    else: # If the regions size is bigger than a bin then:
        bin_id1 = bin_id1[0]
        bin_id2 = bin_id2[0]
        if largeRegionsOperation == 'last':
            bin_id1 = bin_id1[-1]
            bin_id2 = bin_id2[-1]
        elif largeRegionsOperation == 'center':
            bin_id1 = bin_id1[floor(len(bin_id1)/2)]
            bin_id2 = bin_id2[floor(len(bin_id2)/2)]
        elif largeRegionsOperation == 'sum':
            bin_id1 = np.sum(bin_id1)
            bin_id2 = np.sum(bin_id2)
        elif largeRegionsOperation == 'mean':
            bin_id1 = np.mean(bin_id1)
            bin_id2 = np.mean(bin_id2)
        elif largeRegionsOperation == 'median':
            bin_id1 = np.median(bin_id1)
            bin_id2 = np.median(bin_id2)

    print("bin ids:", bin_id1, bin_id2)
    if bin_id1 > bin_id2:
        if chrom1 == chrom2 :
            bin_id1, bin_id2 = sorted(bin_id1, bin_id2)
        else:
            chr1, str1, en1, bin1 = chrom1, start1, end1, bin_id1
            chrom1, start1, end1, bin_id1 = chrom2, start2, end2, bin_id2
            chrom2, start2, end2, bin_id2 = chr1, str1, en1, bin1


    agg_info["counter"] += 1
    if agg_info["counter"] % 50000 == 0:
        log.info("Number of contacts considered: {:,}".format(agg_info["counter"]))

    bin_size = ma.getBinSize()
    chrom1_bin_range = ma.getChrBinRange(toString(chrom1))
    chrom2_bin_range = ma.getChrBinRange(toString(chrom2))

    if chrom1 == chrom2: # chrom1 == chrom2 can happen in intra or all
        if (chrom1 not in agg_info["agg_total"]) and (mode == "intra-chr") and (perChr == True):
            agg_info["agg_total"][chrom1] = 0
            agg_info["agg_matrix"][chrom1] = []
            agg_info["agg_diagonals"][chrom1] = []
            agg_info["agg_contact_position"][chrom1] = []
            agg_info["agg_center_values"][chrom1] = []

        min_dist, max_dist = range.split(":")
        log.info("checking range {}-{}".format(min_dist, max_dist))
        assert int(min_dist) < int(max_dist), "Error lower range larger than upper range"
        min_dist_in_bins = int(min_dist) // bin_size
        max_dist_in_bins = int(max_dist) // bin_size
        if (min_dist_in_bins > abs(bin_id2 - bin_id1)) or (abs(bin_id2 - bin_id1) > max_dist_in_bins):
            print(min_dist_in_bins, max_dist_in_bins)
            return
    print("sort:")
    if (bin_id1, bin_id2) in agg_info["seen"]:
        print("seen")
        return
    agg_info["seen"].append((bin_id1, bin_id2))
    if bin_id1 - M_half < chrom1_bin_range[0] or bin_id2 + 1 + M_half > chrom2_bin_range[1]:
        print("exceeds chr range")
        return
    try:
        print("append!")
        mat_to_append = ma.matrix[bin_id1 - M_half:bin_id1 + M_half + 1, :][:, bin_id2 - M_half:bin_id2 + M_half + 1].todense().astype(float)
    except IndexError:
        log.info("index error for {} {}".format(bin_id1, bin_id2))
        return
    agg_info["counter"] += 1
    if agg_info["counter"] % 1000 == 0:
        log.info("Number of contacts within range computed: {:,}".format(counter))
    if mat_to_append.sum() == 0:
        agg_info["empty_mat"] += 1
        return
    # to account for the fact that submatrices close to the diagonal have more counts than
    # submatrices far from the diagonal submatrices values are normalized using the
    # total submatrix sum.
    if transform == 'total-counts' and mat_to_append.sum() > 0:
        mat_to_append = mat_to_append / mat_to_append.sum()

    if (mode == "intra-chr") and (perChr == True): # chrom1 == chrom2
        agg_info["agg_total"][chrom1] += 1
        agg_info["agg_matrix"][chrom1].append(mat_to_append)
        agg_info["agg_diagonals"][chrom1].append(mat_to_append.diagonal())
        agg_info["agg_center_values"][chrom1].append(ma.matrix[bin_id1, bin_id2])
        agg_info["agg_contact_position"][chrom1].append((start1, end1, start2, end2))

    else:
        if 'genome' not in agg_info["agg_total"]:
            agg_info["agg_total"]["genome"] = 0
            agg_info["agg_matrix"]["genome"] = []
            agg_info["agg_diagonals"]["genome"] = []
            agg_info["agg_contact_position"]["genome"] = []
            agg_info["agg_center_values"]["genome"] = []

        agg_info["agg_total"]["genome"] += 1
        agg_info["agg_matrix"]["genome"].append(mat_to_append)
        agg_info["agg_diagonals"]["genome"].append(mat_to_append.diagonal())
        agg_info["agg_center_values"]["genome"].append(ma.matrix[bin_id1, bin_id2])
        agg_info["agg_contact_position"]["genome"].append((start1, end1, start2, end2))


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
            clustering =  KMeans(n_clusters=k, random_state=0).fit(matrix)
            cluster_labels = clustering.labels_
        if method == 'hierarchical':
            clustering = AgglomerativeClustering(n_clusters=k).fit(matrix)
            cluster_labels = clustering.labels_
        if method == 'spectral':
            clustering = SpectralClustering(n_clusters=k, assign_labels="discretize",random_state=0).fit(matrix)
            cluster_labels = clustering.labels_

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

            if args.operationType == 'median':
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
                # Axes3D is required for projection='3d' to work
                # but since is imported but not used, flake8 will complain
                # thus I add this dummy variable to avoid the error
                Axes3D(fig)
                ax = plt.subplot(gs[cluster_number, idx], projection='3d')
                # ax.set_aspect('equal')
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
    for idx, (chrom_name, values) in enumerate(chrom_diagonals.items()):
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
    ma.maskBins(ma.nan_bins)
    ma.matrix.data = ma.matrix.data
    new_intervals = hicexplorer.utilities.enlarge_bins(ma.cut_intervals)
    ma.setCutIntervals(new_intervals)

    if args.chromosomes:
        ma.keepOnlyTheseChr(args.chromosomes)

    default_range = '1000000:20000000'
    if args.range == None:
        log.warning("You have not set any range. This is by default set to {} for intra-chr.".format(default_range))
        args.range = default_range
    min_dist, max_dist = args.range.split(":")

    if args.transform == "z-score": # use zscore matrix
        log.info("Computing z-score matrix. This may take a while.\n")
        if args.mode == 'intra-chr':
            ma.convert_to_zscore_matrix(maxdepth = max_dist * 2.5, perchr = True)
        else:
            ma.convert_to_zscore_matrix(maxdepth = None, perchr = True)
    elif args.transform == "obs/exp": # use obs/exp matrix
        log.info("Computing observed vs. expected matrix. This may take a while.\n")
        if args.mode == 'intra-chr':
            ma.convert_to_obs_exp_matrix(maxdepth = max_dist * 2.5, perchr=True)
        else:
            ma.convert_to_obs_exp_matrix(maxdepth = None, perchr=True)


    M = args.numberOfBins if args.numberOfBins % 2 == 1 else args.numberOfBins + 1
    M_half = int((M - 1) // 2)

    chrom_coord = dict()
    chrom_list = ma.getChrNames()
    for chrom in chrom_list:
        first, last = ma.getChrBinRange(chrom)
        first = ma.getBinPos(first)
        last = ma.getBinPos(last-1)
        chrom_coord[chrom] = (first[1],last[2])

    agg_info = dict()
    agg_info["chrom_coord"] = chrom_coord
    agg_info["seen"] =[]
    agg_info["agg_matrix"] = OrderedDict()
    agg_info["agg_total"] = {}
    agg_info["agg_diagonals"] = OrderedDict()
    agg_info["agg_contact_position"] = {}
    agg_info["agg_center_values"] = {}
    agg_info["counter"] = 0
    agg_info["empty_mat"] = 0

    if args.row_wise:
        # read bed files
        bed_intervals = args.BED.readlines()
        if args.BED2:
            bed_intervals2 = args.BED2.readlines()
        else:
            log.error("Error computing row-wise contacts requires two bed files!")
            exit("Error computing row-wise contacts requires two bed files!")
        # agg_matrix could be either per chromosome or genome wide
        aggregate_contacts_per_row(bed_intervals, bed_intervals2, agg_info, ma, M_half, args.largeRegionsOperation, args.range, args.transform, mode = args.mode, perChr=args.perChr)
    else: # not row-wise
        # read and sort bedgraph.
        print("not row-wise")
        bed_intervals = read_bed_per_chrom(args.BED)
        if args.BED2:
            bed_intervals2 = read_bed_per_chrom(args.BED2)
        else:
            bed_intervals2 = bed_intervals
        # agg_matrix could be either per chromosome or genome wide
        aggregate_contacts(bed_intervals, bed_intervals2, agg_info, ma, M_half, args.largeRegionsOperation, args.range, args.transform, mode = args.mode, perChr=args.perChr)


    if args.kmeans is not None:
        cluster_ids = cluster_matrices(agg_info["agg_matrix"], args.kmeans, method='kmeans', how=args.howToCluster)
        num_clusters = args.kmeans
    elif args.hclust is not None:
        log.info("Performing hierarchical clustering."
                 "Please note that it might be very slow for large datasets.\n")
        cluster_ids = cluster_matrices(agg_info["agg_matrix"], args.hclust, method='hierarchical',
                                       how=args.howToCluster)
        num_clusters = args.hclust
    else:
        # make a 'fake' clustering to generalize the plotting of the submatrices
        cluster_ids = {}
        num_clusters = 1
        for k in agg_info["agg_matrix"].keys():
            print(k)
            print(agg_info["agg_matrix"][k])
            cluster_ids[k] = [range(len(agg_info["agg_matrix"][k]))]

    plot_aggregated_contacts(agg_info["agg_matrix"], agg_info["agg_contact_position"], cluster_ids, num_clusters, M_half, args)

    if args.outFileContactPairs:
        for idx, chrom in enumerate(chrom_matrix):
            if chrom not in bed_intervals or chrom not in bed_intervals2:
                continue
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
        plot_diagnostic_heatmaps(agg_info["agg_diagonals"], cluster_ids, M_half, args)
