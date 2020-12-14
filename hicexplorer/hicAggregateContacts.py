import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse
import numpy as np
from itertools import compress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
# from scipy.cluster.vq import vq, kmeans
# from scipy.cluster.hierarchy import fcluster, linkage
import sklearn.cluster as skclust
from hicmatrix import HiCMatrix as hm
from hicmatrix.lib import MatrixFileHandler
import hicexplorer.utilities
from .utilities import check_chrom_str_bytes, change_chrom_names, toString
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)
from collections import OrderedDict


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Takes a list of positions in the Hi-C matrix and '
                                                 'makes a pooled image.')
    # define the arguments
    parserRequired = parser.add_argument_group('Required arguments')

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
                           'corresponding row of the BED2 file are computed. If intra-chromosomal '
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
                           'half number of bins and the other half number of bins'
                           ' (Default: %(default)s).',
                           default='51',
                           type=int)

    parserOpt.add_argument('--transform',
                           help='Type of transformation for the matrix. The options are "none",  '
                           '"total-counts", "z-score" and "obs/exp". If total counts are selected, '
                           'the sub-matrix values are divided by the total counts for normalization. '
                           'If z-score or obs/exp are selected, the Hi-C matrix is converted into a '
                           'z-score or observed / expected matrix'
                           ' (Default: %(default)s).',
                           choices=['total-counts', 'z-score', 'obs/exp', 'none'],
                           default='none')

    parserOpt.add_argument('--operationType',
                           help='Type of the operation to be applied to summerize '
                           'the submatrices into a single matrix. Options are sum, '
                           'mean and median. (Default: %(default)s)',
                           choices=['sum', 'mean', 'median'],
                           default='median')

    parserOpt.add_argument('--perChr',
                           help='if set, it generates a plot per chromosome. It is only affected if '
                           'intra-chromosomal contacts are of interest.',
                           action='store_true',
                           required=False)

    parserOpt.add_argument('--largeRegionsOperation',
                           help='If a given coordinate in the bed file is larger than '
                           'a bin of the input matrix, by default only the first bin '
                           'is taken into account. However there are more posibilities '
                           'to handel such a case. Users can ask for the last bin or '
                           'for center of the region. As an example if a region falls into bins [4,5,6] '
                           'and `--numberOfBins = 2` then if first, bins [3,4,5] are kept. '
                           'If last: [5,6,7] and if center: [4,5,6].',
                           choices=['first', 'last', 'center'],
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

    parserOut.add_argument('--outFileObsExp',
                           help='writes the obs/exp matrix to a file, if --transform=obs/exp.',
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
                             help='Options are "full", "center" and "diagonal". The full clustering '
                             'takes all values of each submatrix for clustering. "center", takes only a square of '
                             'length 3x3 from each submatrix and uses only  this values for clustering. With the '
                             '"diagonal" option the clustering is only carried out based on the submatrix diagonal '
                             '(representing values at the same distance to each other)'
                             ' (Default: %(default)s).',
                             choices=['full', 'center', 'diagonal'],
                             default='full')

    parserPlot = parser.add_argument_group('Plotting options')

    parserPlot.add_argument('--chromosomes', '-C',
                            help='List of chromosomes to plot.',
                            nargs='+')

    parserPlot.add_argument('--colorMap',
                            help='Color map to use for the heatmap. Available '
                            'values can be seen here: '
                            'http://matplotlib.org/examples/color/colormaps_reference.html'
                            ' (Default: %(default)s).',
                            default='RdYlBu_r')

    parserPlot.add_argument('--plotType',
                            help='Plot type'
                            ' (Default: %(default)s).',
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
                           'output is a raster graphics image (e.g png, jpg)'
                           ' (Default: %(default)s).',
                           type=int,
                           default=300)
    return parser


def read_bed_per_chrom(fh, chrom_list):
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
        if fields[0] not in chrom_list:
            if change_chrom_names(fields[0]) in chrom_list:
                fields[0] = change_chrom_names(fields[0])
            else:
                continue
        if fields[0] not in interval:
            interval[fields[0]] = []
        interval[fields[0]].append((int(fields[1]), int(fields[2])))
    return interval


def aggregate_contacts(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range=None, transform=None, mode=''):
    """
    To aggregate the contacts of desired sumatrices.
    """
    seen_chrs = []
    for k1, v1 in bed1.items():
        for k2, v2 in bed2.items():
            if (mode == 'inter-chr') & (k1 == k2):
                if (len(bed1) == 1) and (len(bed2) == 1):
                    exit("Error: 'inter-chr' mode needs at least a pair of coordinates with different chromoses to be available in the bed files.")
                else:
                    continue
            if (mode == 'intra-chr') & (k1 != k2):
                continue
            for coord1 in v1:
                for coord2 in v2:
                    if (k1 == k2) and (coord1 == coord2):
                        continue
                    interval = [(k1, coord1[0], coord1[1]), (k2, coord2[0], coord2[1])]
                    count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range, transform)
    to_del = []
    for k1, v1 in agg_info["agg_matrix"].items():
        print("v1: ", v1)
        if v1 == {}:
            print("no matrix on {}".format(k1))
            to_del.append(k1)
        elif (mode == 'intra-chr') and (v1[k1] == []):
            print("no matrix between {} and {}".format(k1, k1))
            to_del.append(k1)
    agg_info["agg_matrix"] = {key:val for key, val in agg_info["agg_matrix"].items() if key not in to_del}
    agg_info["agg_contact_position"] = {key:val for key, val in agg_info["agg_contact_position"].items() if key not in to_del}


def aggregate_contacts_per_row(bed1, bed2, agg_info, ma, chrom_list, M_half, largeRegionsOperation, range=None, transform=None, mode='', perChr=False):
    """
    To aggregate the contacts of the desired submatrices , if row-wise.
    """
    seen_chrs = []
    for line1, line2 in zip(bed1, bed2):
        line1 = line1.strip().split()
        line2 = line2.strip().split()
        if line1[0] not in chrom_list:
            line1[0] = change_chrom_names(line1[0])
            if line1[0] not in chrom_list:
                continue
        if line2[0] not in chrom_list:
            line2[0] = change_chrom_names(line2[0])
            if line2[0] not in chrom_list:
                continue
        if mode == 'inter-chr':  # skip the lines with same chrom
            if line1[0] == line2[0]:
                continue
        elif mode == 'intra-chr':  # skip the lines with different chrom
            if line1[0] != line2[0]:
                continue
        interval = [(line1[0],line1[1], line1[2]), (line2[0],line2[1], line2[2])]
        count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range, transform)
    to_del = []
    for k1, v1 in agg_info["agg_matrix"].items():
        if v1 == {}:
            print("no matrix on {}".format(k1))
            to_del.append(k1)
        elif (mode == 'intra-chr') and (v1[k1] == []):
            print("no matrix between {} and {}".format(k1, k1))
            to_del.append(k1)
    agg_info["agg_matrix"] = {key:val for key, val in agg_info["agg_matrix"].items() if key not in to_del}
    agg_info["agg_contact_position"] = {key:val for key, val in agg_info["agg_contact_position"].items() if key not in to_del}

def count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range=None, transform=None):
    """
    To count the number of contacts for a given pair of intervals
    """

    chrom1, start1, end1 = interval[0]
    chrom2, start2, end2 = interval[1]

    if chrom1 not in agg_info["chrom_coord"]:
        return
    if chrom2 not in agg_info["chrom_coord"]:
        return
    if (int(end1) > agg_info["chrom_coord"][chrom1][1]) or (int(end2) > agg_info["chrom_coord"][chrom2][1]):  # TODO these intervals may still partially be overlapped, shall we keep them?
        return
    if (int(start1) < agg_info["chrom_coord"][chrom1][0]) or (int(start2) < agg_info["chrom_coord"][chrom2][0]):
        return

    bin_id1 = ma.getRegionBinRange(toString(chrom1), start1, end1)
    bin_id2 = ma.getRegionBinRange(toString(chrom2), start2, end2)
    if bin_id1 == bin_id2: # does not count the contact between a bin and itslef
        return
    if (bin_id1 is None) or (bin_id2 is None):
        return
    else:  # If the regions size is bigger than a bin then:
        if largeRegionsOperation == 'first':
            bin_id1 = bin_id1[0]
            bin_id2 = bin_id2[0]
        elif largeRegionsOperation == 'last':
            bin_id1 = bin_id1[-1]
            bin_id2 = bin_id2[-1]
        elif largeRegionsOperation == 'center':
            bin_id1 = int(np.floor(np.mean(bin_id1)))
            bin_id2 = int(np.floor(np.mean(bin_id2)))
    if bin_id1 > bin_id2:
        if chrom1 == chrom2:
            bin_id1, bin_id2 = sorted([bin_id1, bin_id2])
        else:
            chr1, str1, en1, bin1 = chrom1, start1, end1, bin_id1 # intermediate variables
            chrom1, start1, end1, bin_id1 = chrom2, start2, end2, bin_id2
            chrom2, start2, end2, bin_id2 = chr1, str1, en1, bin1
    agg_info["counter"] += 1
    if agg_info["counter"] % 50000 == 0:
        log.info("Number of contacts considered: {:,}".format(agg_info["counter"]))

    bin_size = ma.getBinSize()
    chrom1_bin_range = ma.getChrBinRange(toString(chrom1))
    chrom2_bin_range = ma.getChrBinRange(toString(chrom2))
    # stable:
    # if chrom1 == chrom2:  # chrom1 == chrom2 can happen in intra or all
    #     if (chrom1 not in agg_info["agg_total"]) and (mode == "intra-chr") and (perChr == True):
    #         agg_info["agg_total"][chrom1] = 0
    #         agg_info["agg_matrix"][chrom1] = []
    #         agg_info["agg_diagonals"][chrom1] = []
    #         agg_info["agg_contact_position"][chrom1] = []
    #         agg_info["agg_center_values"][chrom1] = []
    # alternative:
    if [chrom1, chrom2] not in seen_chrs: # we only want to count one side of the diagonal
        agg_info["agg_total"][chrom1][chrom2] = 0
        agg_info["agg_matrix"][chrom1][chrom2] = []
        agg_info["agg_diagonals"][chrom1][chrom2] = []
        agg_info["agg_contact_position"][chrom1][chrom2] = []
        agg_info["agg_center_values"][chrom1][chrom2] = []
        seen_chrs.append([chrom1,chrom2])
# # #
    if mode == "intra-chr":
        min_dist, max_dist = range.split(":")
        min_dist_in_bins = int(min_dist) // bin_size
        max_dist_in_bins = int(max_dist) // bin_size
        if (min_dist_in_bins > abs(bin_id2 - bin_id1)) or (abs(bin_id2 - bin_id1) > max_dist_in_bins):
            print("out of range")
            return
    if (bin_id1, bin_id2) in agg_info["seen"]:
        return
    agg_info["seen"].append((bin_id1, bin_id2))
    if bin_id1 - M_half < chrom1_bin_range[0] or bin_id1 + M_half >= chrom1_bin_range[1]:
        log.info("The given interval exceeds the chromosome range on {}. It is skipped.".format(chrom1))
        return

    if bin_id2 - M_half < chrom2_bin_range[0] or bin_id2 + M_half >= chrom2_bin_range[1]:
        log.info("The given interval exceeds the chromosome range on {}. It is skipped.".format(chrom2))
        return
    try:
        mat_to_append = ma.matrix[bin_id1 - M_half:bin_id1 + M_half + 1, :][:, bin_id2 - M_half:bin_id2 + M_half + 1].todense().astype(float)
    except IndexError:
        log.info("index error for {} {}".format(bin_id1, bin_id2))
        return

    if mat_to_append.sum() == 0:
        agg_info["empty_mat"] += 1
        return

    agg_info["used_counter"] += 1
    if agg_info["used_counter"] % 50000 == 0:
        log.info("Number of used contacts within the given range: {:,}".format(agg_info["used_counter"]))
    # to account for the fact that submatrices close to the diagonal have more counts than
    # submatrices far from the diagonal submatrices values are normalized using the
    # total submatrix sum.
    if transform == 'total-counts' and mat_to_append.sum() > 0:
        mat_to_append = mat_to_append / mat_to_append.sum()
    # stable:
    # if (mode == "intra-chr") and (perChr == True):  # chrom1 == chrom2
    #     agg_info["agg_total"][chrom1] += 1
    #     agg_info["agg_matrix"][chrom1].append(mat_to_append)
    #     agg_info["agg_diagonals"][chrom1].append(mat_to_append.diagonal())
    #     agg_info["agg_center_values"][chrom1].append(ma.matrix[bin_id1, bin_id2])
    #     agg_info["agg_contact_position"][chrom1].append((start1, end1, start2, end2))
    #
    # else:
    #     print(chrom1, chrom2)
    #     if 'genome' not in agg_info["agg_total"]:
    #         agg_info["agg_total"]["genome"] = 0
    #         agg_info["agg_matrix"]["genome"] = []
    #         agg_info["agg_diagonals"]["genome"] = []
    #         agg_info["agg_contact_position"]["genome"] = []
    #         agg_info["agg_center_values"]["genome"] = []
    #
    #     agg_info["agg_total"]["genome"] += 1
    #     agg_info["agg_matrix"]["genome"].append(mat_to_append)
    #     agg_info["agg_diagonals"]["genome"].append(mat_to_append.diagonal())
    #     agg_info["agg_center_values"]["genome"].append(ma.matrix[bin_id1, bin_id2])
    #     agg_info["agg_contact_position"]["genome"].append((start1, end1, start2, end2))
    # Alternative:
    agg_info["agg_total"][chrom1][chrom2] += 1
    agg_info["agg_matrix"][chrom1][chrom2].append(mat_to_append)
    agg_info["agg_diagonals"][chrom1][chrom2].append(mat_to_append.diagonal())
    agg_info["agg_center_values"][chrom1][chrom2].append(ma.matrix[bin_id1, bin_id2])
    agg_info["agg_contact_position"][chrom1][chrom2].append((start1, end1, start2, end2))
    print(chrom1, chrom2, start1, end1, start2, end2)


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


def compute_clusters(updated_info, k, method="kmeans", how='full'):
    submat_vectors = []
    shape = updated_info["submatrices"][0].shape
    center_bin = (shape[0] + 1) // 2
    for submatrix in updated_info["submatrices"]:
        print(submatrix.shape)
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
    print(submat_vectors)
    matrix = np.vstack(submat_vectors)
    if how == 'diagonal':
        assert matrix.shape == (len(updated_info["submatrices"]), shape[0])
    elif how == 'center':
        assert matrix.shape == (len(updated_info["submatrices"]), 1)
    else:
        assert matrix.shape == (len(updated_info["submatrices"]), shape[0] * shape[1])

    # remove outliers
    updated_coords = updated_info["coords"]
    updated_centers = updated_info["centers"]
    out_ind = get_outlier_indices(matrix, max_deviation=10)
    print(out_ind, "outliers")
    if out_ind is not None and len(np.flatnonzero(out_ind)) > 0:
        if len(np.flatnonzero(out_ind)) == matrix.shape[0]:
            exit("ERROR: all submatrices have been detected as outliers. You can consider changing the threshold")
        log.info("Outliers detected. Number of outliers: {}".
                format(len(np.flatnonzero(out_ind))))
        # keep in matrix all indices that are not outliers
        matrix = matrix[np.logical_not(out_ind), :]
        print("matrix: ")
        print(matrix.shape)
        updated_coords = list(compress(updated_info["coords"], np.logical_not(out_ind)))
        updated_centers = list(compress(updated_info["centers"], np.logical_not(out_ind)))
        print(len(updated_coords))
    if np.any(np.isnan(matrix)):
        # replace nans for 0 otherwise kmeans produces a weird behaviour
        log.warning("For clustering nan values have to be replaced by zeros.")
        matrix[np.isnan(matrix)] = 0

    if (k == 1) and (method == 'no_clust'): # no clustering
        cluster_labels = np.asarray([0] * matrix.shape[0])
    if method == 'kmeans':
        clustering = skclust.KMeans(n_clusters=k, random_state=0).fit(matrix)
        cluster_labels = clustering.labels_
    if method == 'hierarchical':
        clustering = skclust.AgglomerativeClustering(n_clusters=k, distance_threshold=None).fit(matrix)
        cluster_labels = clustering.labels_
    if method == 'spectral':
        clustering = skclust.SpectralClustering(n_clusters=k, assign_labels="discretize", random_state=0).fit(matrix)
        cluster_labels = clustering.labels_

    # sort clusters
    clustered_dict = []
    for cluster in range(k):
        cluster_ids = np.flatnonzero(cluster_labels == cluster)
        clustered_dict.append(cluster_ids)
    print(matrix.shape)
    matrix = np.array(matrix).reshape((matrix.shape[0], int(np.sqrt(matrix.shape[1])), int(np.sqrt(matrix.shape[1]))))
    print(matrix.shape)
    print(matrix)
    updated_info={'coords': updated_coords, 'centers': updated_centers, 'submatrices': matrix, 'clustered_dict': clustered_dict}
    return updated_info


def cluster_matrices(submatrices_dict, positions, center_values, k, method='kmeans', how='full', perChr=False):
    """
    clusters the submatrices .

    Parameters
    ----------
    submatrices_dict key: chrom name, values, a list of submatrices
    k number of clusters
    method either kmeans, hierarchical or spectral
    how how to cluster. Options are 'full', 'center' and 'diagonal'. More info in the argparse options

    Returns
    -------

    indices dict key: chrom_name, value: list of list, with one list per cluster with the ids of the submatrices
                 that belong to that list
    """
    updated_info = dict()
    if not perChr:
        updated_info["genome"]={"coords": [], "centers": [], "submatrices": [],
                                "clustered_dict": []}
    for chrom1 in submatrices_dict.keys():
        for chrom2 , matrices in  submatrices_dict[chrom1].items():
            if perChr:
                updated_info[chrom1]={"coords": positions[chrom1][chrom2], "centers": center_values[chrom1][chrom2],
                                      "submatrices": submatrices_dict[chrom1][chrom2],
                                      "clustered_dict": []}
                assert(chrom1 == chrom2)
                log.info("Length of entry on chr {}: {}".format(chrom1, len(submatrices_dict[chrom1][chrom2])))
                if len(submatrices_dict[chrom1][chrom2]) < k: # TODO test it!
                    log.info("number of the submatrices on chromosome {} is less than {}. Clustering is skipped.".format(chrom1, k))
                    # updated_info[chrom1]["clustered_dict"]  = [range(len(submatrices_dict[chrom1][chrom2]))]
                    k = 1
                updated_info[chrom1] = compute_clusters(updated_info[chrom1], k, method, how)
            else:
                updated_info['genome']["submatrices"] += submatrices_dict[chrom1][chrom2]
                # Add all corrdinates in a new container
                for idx, matrix in enumerate(matrices):
                    start1, end1, start2, end2 = positions[chrom1][chrom2][idx]
                    updated_info['genome']["coords"].append((chrom1,start1,end1,chrom2, start2, end2))
                    updated_info['genome']["centers"].append(center_values[chrom1][chrom2][idx])

    if not perChr:
        print(len(updated_info['genome']["submatrices"]))
        print(updated_info['genome']["submatrices"])
        updated_info["genome"]  = compute_clusters( updated_info["genome"], k, method, how)
        print(len(updated_info['genome']["centers"]), "here", len(updated_info['genome']["submatrices"]), len(updated_info['genome']["coords"]))
        print(updated_info['genome']["submatrices"])
        # TODO Do I need to update submatrices and positions here too?

    return updated_info


def compute_avg(submatrices, operationType):
    print(operationType)
    if operationType == 'median':
        _median = np.median(submatrices, axis=0)
        if _median.sum() == 0 or np.isnan(_median.sum()):
            # test if the mean matrix is not zero
            if np.mean(submatrices, axis=0).sum() != 0:
                log.info("The median of the matrices is zero. Consider using "
                         "the mean instead.")
            else:
                log.info("Apparently no matrices could be computed. All of them "
                         "are zeros or nans.")
        return _median
    elif operationType == 'mean':
        return np.mean(submatrices, axis=0)
    else:
        return np.sum(submatrices, axis=0)
    return


def plot_aggregated_contacts_perChr(clustered_info, num_clusters, M_half, args):

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
        assert(v1 != {})
        if chrom1 not in chrom_avg.keys():
            chrom_avg[chrom1] = []

        for cluster_number, cluster_indices in enumerate(clustered_info[chrom1]["clustered_dict"]):
            print(cluster_number, cluster_indices)
            # compute median values
            submatrices = np.array([clustered_info[chrom1]["submatrices"][x] for x in cluster_indices])

            print("sub:")
            print(submatrices)
            chrom_avg[chrom1].append(compute_avg(submatrices, args.operationType))
            log.info("Mean aggregate matrix values: {}".format(chrom_avg[chrom1][cluster_number].mean()))
            log.info("total pairs considered on cluster_{}: "
                     "{}".format(cluster_number + 1, len(cluster_indices)))

            if chrom_avg[chrom1][cluster_number].shape[0] == 0:
                log.debug("matrix for contacts on cluster_{} is empty".format(cluster_number + 1))
                continue
            title = "cluster_{}".format(cluster_number + 1)
            if args.plotType == '2d':
                ax = plt.subplot(gs[cluster_number, idx])
                ax.set_title(title)
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
                    Z = chrom_avg[chrom][cluster_number].copy()
                    img = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                                          vmax=vmax, vmin=vmin, edgecolor='none')

                    ax.set_zticklabels([])
                    if vmax is not None and vmax is not None:
                        ax.set_zlim(vmin, vmax)

        cbar_x = plt.subplot(gs[-1, idx])
        fig.colorbar(img, cax=cbar_x, orientation='horizontal')
    # for idx, chrom in enumerate(chrom_matrix):
    #     for cluster_number, cluster_indices in enumerate(cluster_ids[chrom]):
    #         log.info("total pairs considered for {}, cluster_{}: {}".format(chrom, cluster_number + 1,
    #                                                                         len(cluster_indices)))
    #         try:
    #             chrom_avg[chrom][cluster_number].shape[0]
    #         except IndexError:
    #             continue
    #         if chrom_avg[chrom][cluster_number].shape[0] == 0:
    #             log.debug("matrix for chrom {} is empty".format(chrom))
    #             continue
    #         if num_clusters == 1:
    #             title = chrom
    #         else:
    #             title = "{} cluster_{}".format(chrom, cluster_number + 1)
    #         if args.plotType == '2d':
    #             ax = plt.subplot(gs[cluster_number, idx])
    #
    #             ax.set_title(title)
    #             img = ax.imshow(chrom_avg[chrom][cluster_number], aspect='equal',
    #                             interpolation='nearest', vmax=vmax, vmin=vmin,
    #                             cmap=cmap,
    #                             extent=[-M_half, M_half + 1, -M_half, M_half + 1])
    #         else:
    #             # Axes3D is required for projection='3d' to work
    #             # but since is imported but not used, flake8 will complain
    #             # thus I add this dummy variable to avoid the error
    #             Axes3D(fig)
    #             ax = plt.subplot(gs[cluster_number, idx], projection='3d')
    #             # ax.set_aspect('equal')
    #             ax.margins(0)
    #             X, Y = np.meshgrid(range(-M_half, M_half + 1),
    #                                range(-M_half, M_half + 1))
    #             Z = chrom_avg[chrom][cluster_number].copy()
    #
    #             img = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
    #                                   vmax=vmax, vmin=vmin, edgecolor='none')
    #
    #             ax.set_zticklabels([])
    #             if vmax is not None and vmax is not None:
    #                 ax.set_zlim(vmin, vmax)
    #
    #         if args.outFilePrefixMatrix:
    #             # save aggregate matrix values
    #             if num_clusters == 1:
    #                 output_matrix_name = "{file}_{chrom}.tab".format(file=args.outFilePrefixMatrix, chrom=chrom)
    #             else:
    #                 output_matrix_name = "{file}_{chrom}_cluster_{id}.tab".format(file=args.outFilePrefixMatrix,
    #                                                                               chrom=chrom, id=cluster_number + 1)
    #             np.savetxt(output_matrix_name, chrom_avg[chrom][cluster_number], '%0.5f', delimiter='\t')
    #             print("sum of matrix", np.sum(chrom_avg[chrom][cluster_number]))
    #
    #     cbar_x = plt.subplot(gs[-1, idx])
    #     fig.colorbar(img, cax=cbar_x, orientation='horizontal')
    #
    if args.disable_bbox_tight:
        plt.savefig(args.outFileName.name, dpi=args.dpi)
    else:
        plt.savefig(args.outFileName.name, dpi=args.dpi, bbox_inches='tight')

    plt.close()

def plot_aggregated_contacts(clustered_info, num_clusters, M_half, args):

    fig = plt.figure(figsize=(5.5 , 5.5 * num_clusters + 0.5))
    gs = gridspec.GridSpec(num_clusters+1, 1,
                           width_ratios=[10],
                           height_ratios=[10] * num_clusters + [0.6])

    gs.update(wspace=0.01, hspace=0.2)
    vmin, vmax = (args.vMin, args.vMax)
    cmap = cm.get_cmap(args.colorMap)
    log.debug("vmax: {}, vmin: {}".format(vmax, vmin))
    chrom_avg = []
    for cluster_number, cluster_indices in enumerate(clustered_info["clustered_dict"]):
        # compute median values
        submatrices = np.array([clustered_info["submatrices"][x] for x in cluster_indices])
        print(compute_avg(submatrices, args.operationType))
        chrom_avg.append(compute_avg(submatrices, args.operationType))
        log.info("Mean aggregate matrix values: {}".format(chrom_avg[cluster_number].mean()))
        log.info("total pairs considered on cluster_{}: "
                 "{}".format(cluster_number + 1, len(cluster_indices)))

        if chrom_avg[cluster_number].shape[0] == 0:
            log.debug("matrix for contacts on cluster_{} is empty".format(cluster_number + 1))
            continue
        title = "cluster_{}".format(cluster_number + 1)
        if args.plotType == '2d':
            ax = plt.subplot(gs[cluster_number, 0])
            ax.set_title(title)
            img = ax.imshow(chrom_avg[cluster_number], aspect='equal',
                            interpolation='nearest', vmax=vmax, vmin=vmin,
                            cmap=cmap,
                            extent=[-M_half, M_half + 1, -M_half, M_half + 1])

        else:
            # Axes3D is required for projection='3d' to work
            # but since is imported but not used, flake8 will complain
            # thus I add this dummy variable to avoid the error
            Axes3D(fig)
            ax = plt.subplot(gs[cluster_number, 1], projection='3d')
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
                output_matrix_name = "{file}_genome.tab".format(file=args.outFilePrefixMatrix)
            else:
                output_matrix_name = "{file}_genome_cluster_{id}.tab".format(file=args.outFilePrefixMatrix,
                                                                             id=cluster_number + 1)
            np.savetxt(output_matrix_name, chrom_avg[cluster_number], '%0.5f', delimiter='\t')

        if args.outFileContactPairs:
            center_values_to_order = np.array(clustered_info["centers"])[cluster_indices]
            center_values_order = np.argsort(center_values_to_order)[::-1]

            output_name = "{file}_cluster_{id}.tab".format(file=args.outFileContactPairs,
                                                           id=cluster_number + 1)
            with open(output_name, 'w') as fh:
                for cl_idx in center_values_order:
                    value = center_values_to_order[cl_idx]
                    chrom1, start1, end1, chrom2, start2, end2 = clustered_info["coords"][cl_idx]
                    fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom1, start1, end1, chrom2, start2, end2, value))


        cbar_x = plt.subplot(gs[-1, 0])
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
    if args.range is None:
        if args.mode == "intra-chr":
            log.warning("You have not set any range. This is by default set to {} for intra-chr.".format(default_range))
        args.range = default_range
    min_dist, max_dist = args.range.split(":")
    log.info("checking range {}-{}".format(min_dist, max_dist))
    assert int(min_dist) < int(max_dist), "Error lower range is larger than upper range!"
    if args.transform == "z-score":  # use zscore matrix
        log.info("Computing z-score matrix. This may take a while.\n")
        if args.mode == 'intra-chr':
            ma.convert_to_zscore_matrix(maxdepth=int(max_dist) * 2.5, perchr=True)
        else:
            ma.convert_to_zscore_matrix(maxdepth=None, perchr=True)
    elif args.transform == "obs/exp":  # use obs/exp matrix
        log.info("Computing observed vs. expected matrix. This may take a while.\n")
        if args.mode == 'intra-chr':
            ma.convert_to_obs_exp_matrix(maxdepth=int(max_dist) * 2.5, perchr=True)
        else:
            ma.convert_to_obs_exp_matrix(maxdepth=None, perchr=True)
        if args.outFileObsExp:
            file_type = 'cool'
            if args.outFileObsExp.endswith('.h5'):
                file_type = 'h5'
            matrixFileHandlerOutput = MatrixFileHandler(pFileType=file_type)
            matrixFileHandlerOutput.set_matrix_variables(ma.matrix,
                                                         ma.cut_intervals,
                                                         ma.nan_bins,
                                                         ma.correction_factors,
                                                         ma.distance_counts)
            matrixFileHandlerOutput.save(args.outFileObsExp, pSymmetric=True, pApplyCorrection=False)

    M = args.numberOfBins if args.numberOfBins % 2 == 1 else args.numberOfBins + 1
    M_half = int((M - 1) // 2)

    chrom_coord = dict()
    chrom_list = ma.getChrNames()
    for chrom in chrom_list:
        first, last = ma.getChrBinRange(chrom)
        first = ma.getBinPos(first)
        last = ma.getBinPos(last - 1)
        chrom_coord[chrom] = (first[1], last[2])

    agg_info = dict()
    agg_info["chrom_coord"] = chrom_coord  # coordinates of each chrom
    agg_info["seen"] = [] # seen bins
    agg_info["agg_matrix"] = {chrom:{} for chrom in chrom_list} # important
    agg_info["agg_total"] = {chrom:{} for chrom in chrom_list}
    agg_info["agg_diagonals"] = {chrom:{} for chrom in chrom_list}
    agg_info["agg_contact_position"] = {chrom:{} for chrom in chrom_list} # important
    agg_info["agg_center_values"] = {chrom:{} for chrom in chrom_list} # important
    agg_info["counter"] = 0
    agg_info["used_counter"] = 0
    agg_info["empty_mat"] = 0
    if (args.mode == 'inter-chr') and (len(agg_info["chrom_coord"]) == 1):
        exit("Error: 'inter-chr' mode can not be applied on matrices of only one chromosme.")
    if (args.mode == 'inter-chr') and (args.perChr):
            exit("Error: 'inter-chr' mode can not be used along with --perChr.")
    if args.row_wise:
        # read bed files
        bed_intervals = args.BED.readlines()
        if args.BED2:
            bed_intervals2 = args.BED2.readlines()
        else:
            log.error("Error computing row-wise contacts requires two bed files!")
            exit("Error computing row-wise contacts requires two bed files!")
        if len(bed_intervals) != len(bed_intervals2):
            log.error("row_wise only works if both bed files have the same length.")
            exit("Error row_wise only works if both bed files have the same length.")
        # agg_matrix could be either per chromosome or genome wide
        aggregate_contacts_per_row(bed_intervals, bed_intervals2, agg_info, ma, chrom_list,
                                   M_half, args.largeRegionsOperation, args.range,
                                   args.transform, mode=args.mode, perChr=args.perChr)
    else:  # not row-wise
        # read and sort bed files.
        bed_intervals = read_bed_per_chrom(args.BED, chrom_list)
        if args.BED2:
            bed_intervals2 = read_bed_per_chrom(args.BED2, chrom_list)
        else:
            bed_intervals2 = bed_intervals
        # agg_matrix could be either per chromosome or genome wide
        aggregate_contacts(bed_intervals, bed_intervals2, agg_info, ma, M_half,
                           args.largeRegionsOperation, args.range, args.transform,
                           mode=args.mode)
    if len(agg_info["agg_matrix"]) == 0:
        exit("No susbmatrix found to be aggregated.")

    if args.kmeans is not None:
        assert(args.kmeans > 1)
        if args.perChr == True:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                           agg_info["agg_center_values"],
                                           k=args.kmeans, method='kmeans', how=args.howToCluster,
                                           perChr=args.perChr)
            print(clustered_info)
        else:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                              agg_info["agg_center_values"],
                                              k=args.kmeans, method='kmeans', how=args.howToCluster,
                                              perChr=False)
        num_clusters = args.kmeans
    elif args.hclust is not None:
        assert(args.hclust > 1)
        log.info("Performing hierarchical clustering."
                 "Please note that it might be very slow for large datasets.\n")
        if args.perChr == True:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                           agg_info["agg_center_values"],
                                           k=args.hclust, method='hierarchical',
                                           how=args.howToCluster, perChr=args.perChr)
        else:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                              agg_info["agg_center_values"],
                                              k=args.hclust, method='hierarchical',
                                              how=args.howToCluster, perChr=False)
        num_clusters = args.hclust
    else:
        # make a 'fake' clustering to generalize the plotting of the submatrices
        k = 1
        if args.perChr == True:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                              agg_info["agg_center_values"], k= k, method='no_clust',
                                              how='full', perChr=args.perChr)
            print(clustered_info)
        else:
            clustered_info = cluster_matrices(agg_info["agg_matrix"], agg_info["agg_contact_position"],
                                              agg_info["agg_center_values"], k= k, method='no_clust',
                                              how='full', perChr=False)


        num_clusters = k
    if args.perChr:
        plot_aggregated_contacts_perChr(clustered_info, num_clusters, M_half, args)
    else:
        plot_aggregated_contacts(clustered_info["genome"], num_clusters, M_half, args)
    # stable

    # # plot the diagonals
    # # the diagonals plot is useful to see individual cases and if they had a contact in the center
    # if args.diagnosticHeatmapFile:
    #     plot_diagnostic_heatmaps(agg_info["agg_diagonals"], cluster_ids, M_half, args)
