"""hicAggregateContacts: a thin layer over the C++ executable hicAggregateContacts.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
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
from matplotlib import colormaps as cm
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
    """The argparse parser of hicAggregateContacts, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicAggregateContacts')


def read_bed_per_chrom(fh, chrom_list, pConsiderStrandDirection=None):
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
        if pConsiderStrandDirection and len(fields) < 6:
            log.error('Strand information should be considered, but BED file has not at least six columns. Exiting!')
            exit(1)
        # log.debug('fields {}'.format(fields))
        if fields[0] not in chrom_list:
            if change_chrom_names(fields[0]) in chrom_list:
                fields[0] = change_chrom_names(fields[0])
            else:
                continue
        if fields[0] not in interval:
            interval[fields[0]] = []

        if pConsiderStrandDirection:
            interval[fields[0]].append((int(fields[1]), int(fields[2]), str(fields[5])))
        else:
            interval[fields[0]].append((int(fields[1]), int(fields[2])))
    return interval


def aggregate_contacts(bed1, bed2, agg_info, ma, M_half, largeRegionsOperation, range=None, transform=None, mode='', pConsiderStrandDirection=None):
    """
    To aggregate the contacts of desired submatrices.
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
                    if pConsiderStrandDirection:
                        interval = [(k1, coord1[0], coord1[1], coord1[2]), (k2, coord2[0], coord2[1], coord2[2])]
                    else:
                        interval = [(k1, coord1[0], coord1[1]), (k2, coord2[0], coord2[1])]

                    count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range, transform)
    to_del = []
    for k1, v1 in agg_info["agg_matrix"].items():
        if v1 == {}:
            log.info("no matrix on {}".format(k1))
            to_del.append(k1)
        else:
            new_v1 = {k: v for k, v in v1.items() if v != []}
            if new_v1 != {}:
                agg_info["agg_matrix"][k1] = new_v1
            else:
                to_del.append(k1)
    agg_info["agg_matrix"] = {key: val for key, val in agg_info["agg_matrix"].items() if key not in to_del}
    agg_info["agg_contact_position"] = {key: val for key, val in agg_info["agg_contact_position"].items() if key not in to_del}


def aggregate_contacts_per_row(bed1, bed2, agg_info, ma, chrom_list, M_half, largeRegionsOperation, range=None, transform=None, mode='', perChr=False, pConsiderStrandDirection=None):
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
        if pConsiderStrandDirection:
            interval = [(line1[0], line1[1], line1[2], line1[5]), (line2[0], line2[1], line2[2], line2[5])]
        else:
            interval = [(line1[0], line1[1], line1[2]), (line2[0], line2[1], line2[2])]

        count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range, transform)
    to_del = []
    for k1, v1 in agg_info["agg_matrix"].items():
        if v1 == {}:
            log.info("no matrix on {}".format(k1))
            to_del.append(k1)
        elif (mode == 'intra-chr') and (v1[k1] == []):
            log.info("no matrix between {} and {}".format(k1, k1))
            to_del.append(k1)
    agg_info["agg_matrix"] = {key: val for key, val in agg_info["agg_matrix"].items() if key not in to_del}
    agg_info["agg_contact_position"] = {key: val for key, val in agg_info["agg_contact_position"].items() if key not in to_del}


def count_contacts(interval, ma, M_half, mode, agg_info, largeRegionsOperation, seen_chrs, range=None, transform=None):
    """
    To count the number of contacts for a given pair of intervals
    """
    log.debug('interval {}'.format(interval))

    orientation1 = None
    orientation2 = None
    if len(interval[0]) == 3 and len(interval[1]) == 3:
        chrom1, start1, end1 = interval[0]
        chrom2, start2, end2 = interval[1]
    elif len(interval[0]) == 4 and len(interval[1]) == 4:
        chrom1, start1, end1, orientation1 = interval[0]
        chrom2, start2, end2, orientation2 = interval[1]

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
    if bin_id1 == bin_id2:  # does not count the contact between a bin and itslef
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
            tmp_chrom, tmp_start, tmp_end, tmp_bin = chrom1, start1, end1, bin_id1  # intermediate variables
            chrom1, start1, end1, bin_id1 = chrom2, start2, end2, bin_id2
            chrom2, start2, end2, bin_id2 = tmp_chrom, tmp_start, tmp_end, tmp_bin
    agg_info["counter"] += 1
    if agg_info["counter"] % 50000 == 0:
        log.info("Number of contacts considered: {:,}".format(agg_info["counter"]))

    bin_size = ma.getBinSize()
    chrom1_bin_range = ma.getChrBinRange(toString(chrom1))
    chrom2_bin_range = ma.getChrBinRange(toString(chrom2))

    if [chrom1, chrom2] not in seen_chrs:  # we only want to count one side of the diagonal
        agg_info["agg_total"][chrom1][chrom2] = 0
        agg_info["agg_matrix"][chrom1][chrom2] = []
        agg_info["agg_diagonals"][chrom1][chrom2] = []
        agg_info["agg_contact_position"][chrom1][chrom2] = []
        agg_info["agg_center_values"][chrom1][chrom2] = []
        seen_chrs.append([chrom1, chrom2])

    if mode == "intra-chr":
        min_dist, max_dist = range.split(":")
        min_dist_in_bins = int(min_dist) // bin_size
        max_dist_in_bins = int(max_dist) // bin_size
        if (min_dist_in_bins > abs(bin_id2 - bin_id1)) or (abs(bin_id2 - bin_id1) > max_dist_in_bins):
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

        if (orientation1 is None and orientation2 is None) or (orientation1 == '+' and orientation2 == '+'):
            pass
        elif orientation1 == '+' and orientation2 == '-':
            # flip values concerning the y axis
            mat_to_append = np.fliplr(mat_to_append)
        elif orientation1 == '-' and orientation2 == '+':
            # flip values concerning the x axis
            mat_to_append = np.flipud(mat_to_append)
        elif orientation1 == '-' and orientation2 == '-':
            # flip values concerning the x and y axis
            mat_to_append = mat_to_append.T

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
    agg_info["agg_total"][chrom1][chrom2] += 1
    agg_info["agg_matrix"][chrom1][chrom2].append(mat_to_append)
    agg_info["agg_diagonals"][chrom1][chrom2].append(mat_to_append.diagonal())
    agg_info["agg_center_values"][chrom1][chrom2].append(ma.matrix[bin_id1, bin_id2])
    agg_info["agg_contact_position"][chrom1][chrom2].append((start1, end1, start2, end2))


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


def compute_clusters(updated_info, k, method="kmeans", how='full', max_deviation=2, keep_outlier=False):

    submat_vectors = []
    shape = updated_info["submatrices"][0].shape
    center_bin = (shape[0] + 1) // 2
    for submatrix in updated_info["submatrices"]:
        if how == 'diagonal':
            # take from each matrix the diagonal
            submat_vectors.append(submatrix.diagonal())
        elif how == 'center':
            # take the mean of a  smaller submatrix of 3 x 3 centered on the submatrix
            submat_vectors.append(submatrix[center_bin - 2:center_bin + 1,
                                            center_bin - 2:center_bin + 1].reshape((1, 9)).mean())
        else:
            # Transform list of submatrices in an array of shape:
            # shape = (num_submatrices, submatrix.shape[0] * submatrix.shape[1]
            # In other words, each submatrix is converted into a row of the matrix
            submat_vectors.append(submatrix.reshape((1, shape[0] * shape[1])))
    matrix = np.asarray(np.vstack(submat_vectors))
    if how == 'diagonal':
        assert matrix.shape == (len(updated_info["submatrices"]), shape[0])
    elif how == 'center':
        assert matrix.shape == (len(updated_info["submatrices"]), 1)
    else:
        assert matrix.shape == (len(updated_info["submatrices"]), shape[0] * shape[1])

    # remove outliers
    updated_coords = updated_info["coords"]
    updated_centers = updated_info["centers"]
    updated_submatrices = updated_info["submatrices"]
    updated_diags = updated_info["diagonal"]
    if not keep_outlier:
        out_ind = get_outlier_indices(matrix, max_deviation=max_deviation)
        if out_ind is not None and len(np.flatnonzero(out_ind)) > 0:
            if len(np.flatnonzero(out_ind)) == matrix.shape[0]:
                exit("ERROR: all submatrices have been detected as outliers. You can consider changing the threshold")
                log.info("Outliers detected. Number of outliers: {}".format(len(np.flatnonzero(out_ind))))
            # keep in matrix all indices that are not outliers
            matrix = matrix[np.logical_not(out_ind), :]

            updated_coords = list(compress(updated_info["coords"], np.logical_not(out_ind)))
            updated_centers = list(compress(updated_info["centers"], np.logical_not(out_ind)))
            updated_submatrices = list(compress(updated_info["submatrices"], np.logical_not(out_ind)))
            updated_diags = list(compress(updated_info["diagonal"], np.logical_not(out_ind)))

    if np.any(np.isnan(matrix)):
        # replace nans for 0 otherwise kmeans produces a weird behaviour
        log.warning("For clustering nan values have to be replaced by zeros.")
        matrix[np.isnan(matrix)] = 0

    if (k == 1) and (method == 'no_clust'):  # no clustering
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

    # matrix = np.array(matrix).reshape((matrix.shape[0], int(np.sqrt(matrix.shape[1])), int(np.sqrt(matrix.shape[1]))))
    updated_info = {'coords': updated_coords, 'centers': updated_centers,
                    'submatrices': updated_submatrices, 'clustered_dict': clustered_dict,
                    'diagonal': updated_diags}
    return updated_info


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
    updated_info = dict()
    if not perChr:
        updated_info["genome"] = {"coords": [], "centers": [], "submatrices": [],
                                  "clustered_dict": [], "diagonal": []}
    for chrom1 in agg_info["agg_matrix"].keys():
        for chrom2, matrices in agg_info["agg_matrix"][chrom1].items():
            if perChr:
                full_coords = [(chrom1, i[0], i[1], chrom2, i[2], i[3]) for i in agg_info["agg_contact_position"][chrom1][chrom2]]
                updated_info[chrom1] = {"coords": full_coords, "centers": agg_info["agg_center_values"][chrom1][chrom2],
                                        "submatrices": agg_info["agg_matrix"][chrom1][chrom2],
                                        "clustered_dict": [], "diagonal": agg_info["agg_diagonals"][chrom1][chrom2]}
                assert (chrom1 == chrom2)
                log.info("Length of entry on chr {}: {}".format(chrom1, len(agg_info["agg_matrix"][chrom1][chrom2])))
                if len(agg_info["agg_matrix"][chrom1][chrom2]) < k:
                    log.info("number of the submatrices on chromosome {} is less than {}. Clustering is skipped.".format(chrom1, k))
                    k = 1
                updated_info[chrom1] = compute_clusters(updated_info[chrom1], k, method, how, max_deviation, keep_outlier)
            else:
                updated_info['genome']["submatrices"] += agg_info["agg_matrix"][chrom1][chrom2]
                # Add all corrdinates in a new container
                for idx, matrix in enumerate(matrices):
                    start1, end1, start2, end2 = agg_info["agg_contact_position"][chrom1][chrom2][idx]
                    updated_info['genome']["coords"].append((chrom1, start1, end1, chrom2, start2, end2))
                    updated_info['genome']["centers"].append(agg_info["agg_center_values"][chrom1][chrom2][idx])
                    updated_info['genome']["diagonal"].append(agg_info["agg_diagonals"][chrom1][chrom2][idx])

    if not perChr:
        this_sum = []
        for i in range(len(updated_info['genome']["submatrices"])):
            this_sum.append(updated_info['genome']["submatrices"][i].sum())
        updated_info["genome"] = compute_clusters(updated_info["genome"], k, method, how, max_deviation, keep_outlier)
        # TODO Do I need to update submatrices and positions here too?

    return updated_info


def compute_avg(submatrices, operationType):
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
            # compute median values
            submatrices = np.array([clustered_info[chrom1]["submatrices"][x] for x in cluster_indices])

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

            if args.outFilePrefixMatrix:
                # save aggregate matrix values
                if num_clusters == 1:
                    output_matrix_name = "{file}_{chrom}.tab".format(file=args.outFilePrefixMatrix, chrom=chrom1)
                else:
                    output_matrix_name = "{file}_{chrom}_cluster_{id}.tab".format(file=args.outFilePrefixMatrix,
                                                                                  chrom=chrom1, id=cluster_number + 1)
                np.savetxt(output_matrix_name, chrom_avg[chrom1][cluster_number], '%0.5f', delimiter='\t')

            if args.outFileContactPairs:
                center_values_to_order = np.array(clustered_info[chrom1]["centers"])[cluster_indices]
                center_values_order = np.argsort(center_values_to_order)[::-1]
                if num_clusters == 1:
                    output_name = "{file}_{chrom}.tab".format(file=args.outFileContactPairs, chrom=chrom1)
                else:
                    output_name = "{file}_{chrom}_cluster_{id}.tab".format(file=args.outFileContactPairs,
                                                                           chrom=chrom1, id=cluster_number + 1)
                with open(output_name, 'w') as fh:
                    for cl_idx in center_values_order:
                        value = center_values_to_order[cl_idx]
                        chrom, start, end, chrom2, start2, end2 = clustered_info[chrom1]["coords"][cl_idx]
                        fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, chrom2, start2, end2, value))

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
        try:
            heatmap = np.asarray(np.vstack(clustered_info[chrom_name]['diagonal']))
        except ValueError:
            log.error("Error computing diagnostic heatmap for chrom: {}".format(chrom_name))
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

    file_name = args.diagnosticHeatmapFile.name
    log.info('Heatmap file saved under: {}'.format(file_name))
    plt.savefig(file_name, dpi=args.dpi, bbox_inches='tight')
    plt.close()


def main(args=None):
    """Run the C++ implementation of hicAggregateContacts with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicAggregateContacts')(args)
