
"""

[hic]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/hic_data/s2/hic-norm/RF/merge.npz
title =
colormap = RdYlBu_r
depth = 100000
#min_value =2.8
#max_value = 3.0
transform = log1p
boundaries_file = conductance_vs_hic/boundaries_all.bed
show_masked_bins = yes
x labels = yes
type = arcplot
type = interaction
#optional in case it can not be guessed by the file ending
file_type = hic_matrix


[x-axis]

[spacer]

[bigwig]
file = /data/manke/group/ramirez/HiC-ThomasLing/data/external/Graveley_mRNA-seq/GSM390060_Kc167-4_spa.bw
title = Kc RNA-seq (Cherbas et al.)
color = black
min_value = 0
#max_value = auto
width = 1.5
number of bins = 500
nans to zeros = True
#optional in case it can not be guessed by the file ending
file_type = bigwig

[simple bed]
file = file.bed
title = peaks
color = read
width = 0.5
# optional. If not given is guessed from the file ending
file_type = bed

[genes]
# example of a genes track
# has the same options as a simple
# bed but if the type=genes is given
# the the file is interpreted as gene
# file. If the bed file contains the exon
# structure then this is plottee. Otherwise
# a region **with direction** is plotted.
file = genes.bed
title = genes
color = darkblue
width = 5
# to turn off/on printing of labels
labels = off
# options are 'genes' or 'domains'.
type = genes
# optional. If not given is guessed from the file ending
file_type = bed

[chrom states]
# this is a case of a bed file that is plotted 'collapsed'
# useful to plot chromatin states if the bed file contains
# the color to plot
file = chromatinStates.bed
title = chromatin states
# color is replaced by the color in the bed file
# in this case
color = black
# default behaviour when plotting intervals from a
# bed file is to 'expand' them such that they
# do not overlap. The display = collapsed
# directive overlaps the intervals.
display = collapsed
width = 0.3

[bedgraph]
file = file.bg
title = bedgraph track
color = green
width = 0.2
# optional, otherwise guseed from file ending
file_type = bedgraph


[bedgraph matrix]
# a bedgraph matrix file is like a bedgraph, except that per bin there
# are more than one value separated by tab: E.g.
# chrX	18279	40131	0.399113	0.364118	0.320857	0.274307
# chrX	40132	54262	0.479340	0.425471	0.366541	0.324736
file = spectra_conductance.bm
title = conductance spectra
width = 1.5
orientation = inverted
min_value = 0.10
max_value = 0.70
# if type is set as lines, then the TAD score lines are drawn instead
# of the matrix
# set to lines if a heatmap representing the matrix
# is not wanted
type = lines
file_type = bedgraph_matrix

[dendrogram]
# a dendrogram is saved in a bed file that contains
# further fields defining the linkage.
file =  linkage.bed
title = dendrogram
width = 2
orientation = inverted
# required on this case. Otherwise
# the bed file will be plotted as regions
file_type = dendrogram
# y positions to draw horizontal lines
hlines = 0.2 0.3

[vlines]
# vertical dotted lines from the top to the bottom of the figure
# can be drawn. For this a bed file is required
# but only the first two columns (chromosome name and start
# are used.
# vlines can also be given at the command line as a list
# of genomic positions. However, sometimes to give a file
# is more convenient in case many lines want to be plotted.
file = regions.bed
type = vlines


"""

from __future__ import division
import sys
from collections import OrderedDict

import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import enlarge_bins
from hicexplorer._version import __version__

from scipy import sparse
from scipy.sparse import triu, dia_matrix

import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import mpl_toolkits.axisartist as axisartist
from bx.intervals.intersection import IntervalTree, Interval


DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 3  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.95, 0.05)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.12, 'top': 0.9}


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description='Plots the diagonal,  and some values close to '
        'the diagonal of a  HiC matrix. The diagonal of the matrix is '
        'plotted horizontally for a region. I will not draw the diagonal '
        'for the whole chromosome',
        usage="%(prog)s --tracks tracks.ini --region chr1:1000000-4000000 -o image.png")

    parser.add_argument('--tracks',
                        help='File containing the instructions to plot the tracks ',
                        type=argparse.FileType('r'),
                        required=True,
                        )

    parser.add_argument('--width',
                        help='figure width in centimeters',
                        type=float,
                        default=DEFAULT_FIGURE_WIDTH)

    parser.add_argument('--height',
                        help='Figure height in centimeters. If not given, the figure height is computed '
                             'based on the widths of the tracks and on the depth '
                             'of the Hi-C tracks, if any. Setting this '
                             'this parameter can cause the Hi-C tracks to look collapsed or expanded.',
                        type=float)

    parser.add_argument('--title', '-t',
                        help='Plot title',
                        required=False)

    parser.add_argument('--scoreName', '-s',
                        help='Score name',
                        required=False)

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image. ',
                        type=argparse.FileType('w'),
                        required=True)

    parser.add_argument('--region',
                        help='Plot only this region. The format is '
                        'chr:start-end ',
                        required=True
                        )

    parser.add_argument('--zMax',
                        help='zMax',
                        type=float,
                        default=None)

    parser.add_argument('--vlines',
                        help='Genomic cooordindates separated by space. E.g. '
                        '--vlines 150000 3000000 124838433 ',
                        type=int,
                        nargs='+'
                        )

    parser.add_argument('--fontSize',
                        help='Font size for the labels of the plot',
                        type=float,
                        )

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             'ouput is a raster graphics image (e.g png, jpg)',
                        type=int,
                        default=72
                        )

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def plot_dendrogram(ax, label_ax, properties, region):
    """
    Uses a linkage stored in a BED format to plot a dendrogram.
    The BED format for the linkage is:
    chrom pos_a pos_b, id, distance, ., id_cluster_a, id_cluster_b, num_clusters


    Each chromosome has a separate linkage

    """
    chrom_region, start_region, end_region = region

    linkage_file = open(properties['file'], 'r')

    # the order of the file is important and must not be sorted.
    # thus, the chromosomes are supposed to be one after the other
    z_value = []
    for line in linkage_file.readlines():
        [chrom, pos_a, pos_b, clust_id, distance, strand, id_cluster_a, id_cluster_b, num_clusters] = \
            line.strip().split('\t')
        if chrom != chrom_region:
            continue

        try:
            pos_a = int(pos_a)
            pos_b = int(pos_b)
            id_cluster_a = int(id_cluster_a)
            id_cluster_b = int(id_cluster_b)
            distance = float(distance)

        except ValueError:
            exit("BED values not valid")

        z_value.append((id_cluster_a, id_cluster_b, distance, num_clusters, pos_a, pos_b))

    boxes = _dendrogram_calculate_info(z_value)
    min_y = 1
    max_y = 0
    for box in boxes:
        if box[1, :].min() < min_y:
            min_y = box[1, :].min()
        if box[1, :].max() > max_y:
            max_y = box[1, :].max()

        ax.plot(box[0, :], box[1, :], 'black')

    if 'hlines' in properties and properties['hlines'] != '':
        # split by space
        for hline in properties['hlines'].split(" "):
            try:
                hline = float(hline)
            except ValueError:
                sys.stderr.write("hlines value: {} in dendrogram is not valid.\n".format(hline))

            ax.hlines(hline, start_region, end_region, 'red', '--')

    if 'orientation' in properties and properties['orientation'] == 'inverted':
        ax.set_ylim(max_y, min_y)
    else:
        ax.set_ylim(min_y, max_y)

    ax.set_xlim(start_region, end_region)
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    label_ax.text(0.15, 0, properties['title'], horizontalalignment='left', size='large',
                  verticalalignment='bottom', transform=label_ax.transAxes)


def _dendrogram_calculate_info(z):
    """
    :param z: list of 6-tuples containing the linkage
    information in :func:`hierarchical_clustering`

    >>> _z_ = [(1, 2, 0.5, 2, 0, 10), (4, 3, 0.6, 3, 5, 15)]
    >>> _dendrogram_calculate_info(_z_)
    [array([[  0. ,   0. ,  10. ,  10. ],
           [  0.5,   0.5,   0.5,   0.5]]), array([[  5. ,   5. ,  15. ,  15. ],
           [  0.6,   0.6,   0.6,   0.5]])]
    """
    boxes = []

    def is_cluster_leaf(cluster_id):
        """
        A cluster is a leaf if the id is less than
        the length of Z
        """

        return False if cluster_id >= num_leafs else True

    def prev_cluster_y(cluster_id):
        """
        The prev_cluster y has as index:
        cluster_id - num_leafs and the
        distance is found at position 2.
        """
        return z[cluster_id - num_leafs][2]

    # at each iteration a sort of box is drawn:
    #
    #           _______
    #   |       |     |
    #   y       |       pos-b, y - prev_cluster_y(id_cluster_b)
    #   |       |
    #           pos_a, y - prev_cluster_y(id_cluster_a)
    #
    # ``y`` is the distance.
    #
    # Four points are required to define such box which
    # are obtained in the following code

    num_leafs = len(z) + 1
    for id_cluster_a, id_cluster_b, distance, num_clusters, pos_a, pos_b in z:
        if is_cluster_leaf(id_cluster_a):
            y_a = 0.5
        else:
            y_a = prev_cluster_y(id_cluster_a)

        if is_cluster_leaf(id_cluster_b):
            y_b = 0.5
        else:
            y_b = prev_cluster_y(id_cluster_b)

        boxes.append(np.array([[pos_a, pos_a, pos_b, pos_b],
                               [y_a, distance, distance, y_b]]))
    return boxes


def plot_boundaries(ax, file_name, region):
    """
    Plots the boundaries as triangles in the given ax.

    :param ax:
    :param file_name: boundaries file in bed format.
                      The file should have either the boundary positions or
                      the continuous genomic intervals. Numerous files, separated by
                      a space are accepted.
    :param region: chromosome region to plot
    :return:
    """
    chrom_region, start_region, end_region = region

    file_names = [x for x in file_name.split(" ") if x != '']

    for file_name in file_names:
        file_h = open(file_name, 'r')

        # get intervals as list
        intervals = []
        line_number = 0
        prev_chrom = None
        prev_start = -1
        prev_line = None
        for line in file_h.readlines():
            line_number += 1
            if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = None
            start = None
            end = None
            try:
                chrom, start, end = fields[0:3]
            except Exception as detail:
                msg = "Error reading line: {}\nError message: {}".format(line_number, detail)
                exit(msg)

            if prev_chrom == chrom:
                assert prev_start <= start, \
                    "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

            try:
                start = int(start)
            except ValueError as detail:
                msg = "Error reading line: {}. The start field is not " \
                      "an integer.\nError message: {}".format(line_number, detail)
                exit(msg)

            try:
                end = int(end)
            except ValueError as detail:
                msg = "Error reading line: {}. The end field is not " \
                      "an integer.\nError message: {}".format(line_number, detail)
                exit(msg)

            # use only intervals in the selected region
            if chrom == chrom_region and \
               start_region - 100000 <= start and \
               end_region + 100000 >= end:

                assert start > 0, \
                    "negative values found in bed file in line {}".format(line)

                intervals.append((start, end))

            prev_line = line
            prev_start = start
            prev_chrom = chrom

        file_h.close()
        if len(intervals) == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(file_name))
            continue
        start, end = zip(*intervals)
        start_list = np.array(start)
        end = np.array(end)
        # check if intervals are consecutive or 1bp positions demarcating the boundaries
        if np.any(end - start_list > 1):
            # intervals are consecutive, but only the boundaries are need.
            start_list = end

        prev_start = -1
        x = []
        y = []
        for start in start_list:
            if prev_start is None:
                # draw only half a triangle
                length = start - prev_start
                x1 = prev_start
                y1 = length
                x2 = start
                y2 = 0
                x.extend([x1, x2])
                y.extend([y1, y2])
            else:
                x1 = prev_start
                x2 = x1 + (start - prev_start) / 2
                x3 = start
                y1 = 0
                y2 = (start - prev_start)
                y3 = 0
                x.extend([x1, x2, x3])
                y.extend([y1, y2, y3])

            prev_start = start
        ax.plot(x, y,  color='black')


def plot_matrix(ax, label_ax, cbar_ax, matrix_properties, region):

    matrix_file = matrix_properties['file']
    hic_ma = HiCMatrix.hiCMatrix(matrix_file)
    chrom, region_start, region_end = region
    hic_ma.keepOnlyTheseChr(chrom)
    if 'show_masked_bins' in matrix_properties and \
            matrix_properties['show_masked_bins'] == 'yes':
        pass
    else:
        hic_ma.maskBins(hic_ma.nan_bins)
    if 'orientation' in matrix_properties and \
            matrix_properties['orientation'] == 'inverted':
        plot_inverted = True
    else:
        plot_inverted = False

    new_intervals = enlarge_bins(hic_ma.cut_intervals)
    hic_ma.interval_trees, hic_ma.chrBinBoundaries = \
        hic_ma.intervalListToIntervalTree(new_intervals)

    hic_ma.cut_intervals = new_intervals

    # expand region to plus depth on both sides
    # to avoid a 45 degree 'cut' on the edges

    # get bin id of start and end of region in given chromosome
    chr_start_id, chr_end_id = hic_ma.getChrBinRange(chrom)
    chr_start = hic_ma.cut_intervals[chr_start_id][1]
    chr_end = hic_ma.cut_intervals[chr_end_id-1][1]
    start_bp = max(chr_start, region_start - matrix_properties['depth'])
    end_bp = min(chr_end, region_end + matrix_properties['depth'])

    idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                           enumerate(hic_ma.cut_intervals)
                           if x[0] == chrom and x[1] >= start_bp and x[2] <= end_bp])

    idx = idx[0:-1]
    # select only relevant matrix part
    hic_ma.matrix = hic_ma.matrix[idx, :][:, idx]
    try:
        new_nan_bins = hic_ma.nan_bins[np.in1d(hic_ma.nan_bins, idx)]
        hic_ma.nan_bins = new_nan_bins - idx[0]
    except:
        pass
    """
    print "filling matrix gaps..."
    hic_ma.matrix, _ = fill_gaps(hic_ma, depth_bins, False)
    """

    # fill the main diagonal, otherwise it looks
    # not so good. The main diagonal is filled
    # with an array containing the max value found
    # in the matrix
    if sum(hic_ma.matrix.diagonal()) == 0:
        print "filling main diagonal because is empty and " \
            "otherwise it looks bad..."
        max_value = hic_ma.matrix.data.max()
        main_diagonal = dia_matrix(([max_value]*hic_ma.matrix.shape[0], [0]),
                                   shape=hic_ma.matrix.shape)
        hic_ma.matrix = hic_ma.matrix + main_diagonal
    """
    # code to select some meaningful max and min values
    min_zscore = -5
    max_zscore = None

    mean = np.mean(hic_ma.matrix.data)
    std = np.std(hic_ma.matrix.data)
    z_score = (hic_ma.matrix.data - mean) / std
    hic_ma.matrix.data[z_score < min_zscore] = 0
    if max_zscore is not None:
        min_ = hic_ma.matrix.data[z_score >= max_zscore].min()
        hic_ma.matrix.data[z_score >= max_zscore] = min_
        print "{}, {}".format(mean, min_)
    hic_ma.matrix.eliminate_zeros()
    """

    # select only the upper triangle of the matrix
    hic_ma.matrix = triu(hic_ma.matrix, k=0, format='csr')

    matrix = np.asarray(hic_ma.matrix.todense().astype(float))

    norm = None

    if 'transform' in matrix_properties:
        if matrix_properties['transform'] == 'log1p':
            matrix += 1
            norm = LogNorm()

        elif matrix_properties['transform'] == '-log':
            mask = matrix == 0
            matrix[mask] = matrix[mask == False].min()
            matrix = -1 * np.log(matrix)

    if 'max_value' in matrix_properties and matrix_properties['max_value'] != 'auto':
        vmax = matrix_properties['max_value']

    else:
        # try to use a 'aesthetically pleasant' max value
        vmax = np.percentile(matrix.diagonal(1), 80)

    if 'min_value' in matrix_properties and matrix_properties['min_value'] != 'auto':
        vmin = matrix_properties['min_value']
    else:
        bin_size = hic_ma.getBinSize()
        depth_bins = int(matrix_properties['depth'] / bin_size)
        vmin = np.median(matrix.diagonal(depth_bins))

    sys.stderr.write("setting min, max values for track {} to: {}, {}\n".format(matrix_properties['section_name'],
                                                                                vmin, vmax))

    if 'colormap' not in matrix_properties:
        matrix_properties['colormap'] = DEFAULT_MATRIX_COLORMAP

    cmap = cm.get_cmap(matrix_properties['colormap'])
    cmap.set_bad('white')

    img = pcolormesh_45deg(matrix, ax, start_pos, vmax=vmax,
                           vmin=vmin, cmap=cmap, norm=norm)

    img.set_rasterized(True)
    if plot_inverted:
        ax.set_ylim(matrix_properties['depth'], 0)
    else:
        ax.set_ylim(0, matrix_properties['depth'])

    # ##plot boundaries
    # if a boundaries file is given, plot the
    # tad boundaries as line delineating the TAD triangles
    if 'boundaries_file' in matrix_properties:
        plot_boundaries(ax, matrix_properties['boundaries_file'], region)

    ax.set_xlim(region_start, region_end)
    if 'x labels' in matrix_properties and matrix_properties['x labels'] != 'no':
        ticks = ax.get_xticks()
        labels = ["{:.2f}".format((x / 1e6))
                  for x in ticks]
        labels[-1] += "Mbp"
        ax.get_xaxis().set_tick_params(
            which='both',
            bottom='on',
            top='off',
            direction='out')

        ax.set_xticklabels(labels)
    else:
        ax.get_xaxis().set_tick_params(
            which='both',
            bottom='off',
            top='off',
            direction='out')
        ax.axes.get_xaxis().set_visible(False)

#    ax.xaxis.tick_top()

    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    cbar_ax.patch.set_alpha(0.0)
    cobar = plt.colorbar(img, ax=cbar_ax, fraction=0.95)
    cobar.solids.set_edgecolor("face")
    label_ax.text(0.3, 0.0, matrix_properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom', transform=label_ax.transAxes)
    # plt.subplots_adjust(wspace=0, hspace=0.1, top=0.95,
    #                    bottom=0.05, left=0.1, right=0.5)


def plot_x_axis(ax, region, properties):
    chrom_region, region_start, region_end = region
    ax.set_xlim(region_start, region_end)
    ticks = ax.get_xticks()
    if ticks[-1] - ticks[1] <= 1e5:
        labels = ["{:,.0f} kb".format((x / 1e3))
                  for x in ticks]

    elif 1e5 < ticks[-1] - ticks[1] < 4e6:
        labels = ["{:,.0f} kb".format((x / 1e3))
                  for x in ticks]
    else:
        labels = ["{:,.0f}Mbp".format((x / 1e6))
                  for x in ticks]
        # labels[-1] += "Mbp"

    ax.axis["x"] = ax.new_floating_axis(0, 0.5)
    
    ax.axis["x"].axis.set_ticklabels(labels)
    ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

    if 'fontsize' in properties:
        ax.axis["x"].major_ticklabels.set(size=int(properties['fontsize']))

    if 'where' in properties and properties['where'] == 'top':
        ax.axis["x"].set_axis_direction("top")


def pcolormesh_45deg(matrix_c, ax, start_pos_vector, vmin=None,
                     vmax=None, cmap=None, norm=None):
    """
    Turns the matrix 45 degrees and adjusts the
    bins to match the actual start end positions.
    """
    import itertools
    # code for rotating the image 45 degrees
    n = matrix_c.shape[0]
    # create rotation/scaling matrix
    t = np.array([[1, 0.5], [-1, 0.5]])
    # create coordinate matrix and transform it
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    # this is to convert the indices into bp ranges
    x = matrix_a[:, 1].reshape(n+1, n+1)
    y = matrix_a[:, 0].reshape(n+1, n+1)
    # plot
    im = ax.pcolormesh(x, y, np.flipud(matrix_c),
                       vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)
    return im


def fill_gaps(hic_ma, depth, fill_contiguous=False):
    """ try to fill in the gaps in the matrix by adding the average values of
    the neighboring rows and cols.

    This method produces best results when the
    missing data is not consecutive.

    """
    from scipy import sparse

    M, N = hic_ma.matrix.shape
    # get only a hic_ma.matrixtrix that contains only
    # the [deph:-depth] diagonals
    hic_ma.matrix = sparse.triu(hic_ma.matrix, k=-depth)
    hic_ma.matrix = sparse.tril(hic_ma.matrix, k=depth, format='csr')

    fill_ma = hic_ma.matrix.copy().tolil()
    if fill_contiguous is True:
        good_nans = hic_ma.nan_bins
        consecutive_nan_idx = np.array([])
    else:
        # find stretches of consecutive nans
        consecutive_nan_idx = np.flatnonzero(np.diff(hic_ma.nan_bins) == 1)
        # the banned list of indices is equal to the actual list
        # and the list plus one, to identify both consecutive nans
        consecutive_nan_idx = np.concatenate([consecutive_nan_idx,
                                              consecutive_nan_idx+1])
        # find the nans that are not consecutive
        good_nans = [x for idx, x in enumerate(hic_ma.nan_bins)
                     if idx not in consecutive_nan_idx]

    for missing_bin in good_nans:
        if 0 < missing_bin < M - 1:
            # the new row value is the mean between the upper
            # and lower rows
            fill_ma[missing_bin, :] = (hic_ma.matrix[missing_bin - 1, :] +
                                       hic_ma.matrix[missing_bin + 1, :]) / 2

            # same for cols
            fill_ma[:, missing_bin] = (hic_ma.matrix[:, missing_bin - 1] +
                                       hic_ma.matrix[:, missing_bin + 1]) / 2

    # identify the intersection points of the failed regions because their
    # neighbors get wrong values
    for bin_a in good_nans:
        for bin_b in good_nans:
            if 0 < bin_a < M and \
                    0 < bin_b < M:
                # the fill value is the average over the
                # neighbors that do have a value

                fill_value = np.mean([
                    hic_ma.matrix[bin_a-1, bin_b-1],
                    hic_ma.matrix[bin_a-1, bin_b+1],
                    hic_ma.matrix[bin_a+1, bin_b-1],
                    hic_ma.matrix[bin_a+1, bin_b+1],
                    ])

                fill_ma[bin_a, bin_b] = fill_value

    # return the matrix and the bins that continue to be nan
    return fill_ma.tocsr(), np.sort(consecutive_nan_idx)


def plot_hic_arcs(ax, label_ax, properties, region):

    """
    Makes and arc connecting two points on a linear scale representing
    interactions between Hi-C bins.
    :param ax: matplotlib axis
    :param label_ax: matplotlib axis for labels
    :param region: tuple containing (chrom, start, end) that will limit the view to
    only the genomic positions between start and end

    :param properties: Dictionary of values that should include:
        at least the file. Optionally: lower threshold, value to filter out low scoring
        enrichments, alpha: (transparency) value for arc lines,  color: default color
        for arcs unless specified in the intervals
        intervals: list of the format chr1:1000000:200000:red
         chr3:13000000:14000000:blue. Only the arcs starting or ending on
         such regions will be plotted. The last field field of the intervals used
         to set the color of the arcs.
    """
    from matplotlib.colors import colorConverter
    from matplotlib.patches import Arc

    hic_matrix = HiCMatrix.hiCMatrix(properties['file'])
    chrom, region_start, region_end = region
    hic_matrix.keepOnlyTheseChr(chrom)
    hic_matrix.diagflat()
    hic_matrix.keepOnlyTheseChr(chrom)

    # filter low scoring enrichments
    hic_matrix.matrix.data[np.isnan(hic_matrix.matrix.data)] = 0
    if 'lower threshold' not in properties:
        sys.stderr.write("Plotting arcs without a lower threshold. This can result "
                         "in too many interactions being plotted. Consider adding a "
                         "'lower threshold' value to the configuration file.")
    else:
        try:
            lower_threshold = float(properties['lower threshold'])
        except ValueError:
            sys.exit("lower threshold value is invalid: {} for {}".format(properties['lower threshold'],
                                                                          properties['file']))
        hic_matrix.matrix.data[hic_matrix.matrix.data < lower_threshold] = 0

    hic_matrix.matrix.eliminate_zeros()
    mat = triu(hic_matrix.matrix, k=0, format='coo')
    max_radius = 0
    count = 0

    if properties['color']:
        color_main = properties['color']
    else:
        color_main = 'blue'

    if properties['alpha']:
        alpha = float(properties['alpha'])
    else:
        alpha = 0.8

    # expand view point
    vp_chrom, vp_start, vp_end = properties['view point'].split(':')
    if vp_chrom == chrom:
        vp_intval = IntervalTree()
        vp_intval.insert_interval(Interval(int(vp_start), int(vp_end)))
    else:
        vp_intval = None
    # process intervals, if any, to plot
    # only interactions from those regions
    intervals = None
    if 'highlights' in properties:
        intervals = {chrom: IntervalTree()}
        if properties['highlights'].endswith(".bed"):
            # the highlights are a bed file
            with open(properties['highlights']) as f:
                for line in f.readlines():
                    if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')

                    try:
                        chrom_vp, start, end = fields[0:3]
                    except Exception as detail:
                        msg = "Error reading line: {}\nError message: {}".format(detail)
                        exit(msg)
                    # skip regions not in the same chromosome
                    if chrom_vp != chrom:
                        continue
                    start = int(start)
                    end = int(end)
                    # insert the interval only if it does not overlap with the view point
                    if vp_intval is not None:
                        if len(vp_intval.find(start, end+1)) > 0:
                            sys.stderr.write("skipping region\n".format(start, end))
                            continue

                    intervals[chrom_vp].insert_interval(Interval(start, end, value='red'))
        else:
            for rangebp in properties['highlights']:
                fields = rangebp.split(":")
                chrom_vp = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if len(fields) == 4:
                    color = fields[3]
                else:
                    color = 'blue'
                print (chrom, start, end, color)
                intervals[chrom_vp].insert_interval(Interval(start, end, value=color))

    for idx, (row, col) in enumerate(np.vstack([mat.row, mat.col]).T):
        chrom1, start1, end1, _ = hic_matrix.cut_intervals[row]
        chrom2, start2, end2, _ = hic_matrix.cut_intervals[col]
        if intervals:
            match1 = intervals[chrom1].find(start1, end1)
            match2 = intervals[chrom2].find(start2, end2)
            if len(match1) == 0 and len(match2) == 0:
                # continue
                color = color_main
            else:
                if len(match1) > 0:
                    color = colorConverter.to_rgba(match1[0].value, alpha=alpha)
                else:
                    color = colorConverter.to_rgba(match2[0].value, alpha=alpha)

        center = start1 + float(start2 - start1) / 2
        radius = start2 - start1
        if radius > max_radius:
            max_radius = radius
        count += 1
        ax.plot([center], [radius])
        if 'line width' in properties:
            line_width = float(properties['line width'])
        else:
            line_width = 0.5*np.sqrt(mat.data[idx])
        ax.add_patch(Arc((center, 0), radius,
                         radius*2, 0, 0, 180, color=color, lw=line_width))

    print "{} arcs plotted".format(count)
    chrom, region_start, region_end = region
    if 'orientation' in properties and properties['orientation'] == 'inverted':
        ax.set_ylim(region_end, -1)
    else:
        ax.set_ylim(-1, region_end)
    ax.set_xlim(region_start, region_end)
    ax.axis('off')

    label_ax.text(0.3, 0.0, properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom', transform=label_ax.transAxes)


def plot_interactions(ax, matrix_properties, region):
    """
    Plots arrows from a viewpoint to nearby
    positions that are in 'regions'
    """
    matrix_file = matrix_properties['file']
    hic_ma = HiCMatrix.hiCMatrix(matrix_file)
    chrom, region_start, region_end = region
    hic_ma.keepOnlyTheseChr(chrom)

    # hard coded threshold
    pval_threshold = float(matrix_properties['extra'][1])
    # hard coded view point (rox2)
    viewpoint = (matrix_properties['extra'][2],
                 int(matrix_properties['extra'][3]),
                 int(matrix_properties['extra'][3])+1)
    viewpoint_bin = hic_ma.getRegionBinRange(*viewpoint)

    # extend view point to include left and right bins
    mat = hic_ma.matrix[viewpoint_bin[0]-1:viewpoint_bin[0]+2, :]
    mat.data[mat.data <= -np.log(pval_threshold)] = 0
    mat.data[np.isnan(mat.data)] = 0
    mat.eliminate_zeros()
    mat = mat.tocoo()

    for target in np.unique(mat.col):
        t_chrom, t_start, t_end, _ = hic_ma.cut_intervals[target]
        x_pos = t_end - (t_end - t_start)/2
        rad = '-0.1' if x_pos < viewpoint[1] else '0.1'
        ax.annotate(" ", xy=(x_pos, 0), xytext=(viewpoint[1], 0),
                    xycoords='data', color='black',
                    arrowprops=dict(arrowstyle="simple", lw=10,
                                    facecolor='black', edgecolor='none',
                                    connectionstyle="arc3,rad={}".format(rad))
                    )
    ax.set_ylim(-5, 0)
    ax.set_xlim(region_start, region_end)

    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)


def get_gene_parts(fields, bed_properties):
    chrom, start, end = fields[0:3]
    start = int(start)
    end = int(end)
    try:
        name = fields[3]
    except IndexError:
        name = ''
    try:
        strand = fields[5]
    except IndexError:
        strand = "."
    try:
        thick_start = int(fields[6])
    except (IndexError, ValueError):
        thick_start = int(fields[1])
    try:
        thick_end = int(fields[7])
    except (IndexError, ValueError):
        thick_end = int(fields[2])

    try:
        rgb = fields[8]
        rgb = rgb.split(',')
        if len(rgb) == 3:
            rgb = [float(x)/255 for x in rgb]
            # edgecolor = 'black'
        else:
            rgb = bed_properties['color']
            # edgecolor = bed_properties['color']
    except IndexError:
        rgb = bed_properties['color']
        # edgecolor = bed_properties['color']
    try:
        block_count = int(fields[9])
        block_sizes = [int(x) for x in fields[10].split(',') if x.strip() != '']
        block_starts = [int(x) for x in fields[11].split(',') if x.strip() != '']
#    except IndexError, ValueError:
    except:
        block_count = 0
        block_sizes = []
        block_starts = []

    return chrom, start, end, name, strand, thick_start, thick_end, rgb, block_count, block_sizes, block_starts


def draw_gene_simple(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor):
    """
    draws a gene using different styles
    """
    from matplotlib.patches import Polygon
    (chrom, start, end, name, strand, thick_start, thick_end, rgb_, block_count,
     block_sizes, block_starts) = get_gene_parts(fields, bed_properties)
    # prepare the vertices for all elements that will be drawn
    # the region length without the tip of the end arrow whose length is 'small_relative'
    if strand == '+':
        x0 = start
        x1 = end - small_relative
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        -----------------\
        ---------------- /

        """

        vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + small_relative, y0 + 50), (x1, y0)]

    elif strand == '-':
        x0 = start + small_relative
        x1 = end
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        /--------___---------_
        \--------   ----------
        """
        vertices = [(x0, y0), (x0 - small_relative, y0 + 50), (x0, y1), (x1, y1), (x1, y0)]

    ax.add_patch(Polygon(vertices, closed=True, fill=True,
                         edgecolor=edgecolor,
                         facecolor=rgb))

    if 'labels' in bed_properties and bed_properties['labels'] != 'off':
        center = start + float(end - start)/2
        ax.text(center, ypos + 125, fields[3], size='small',
                horizontalalignment='center', verticalalignment='top')


def draw_gene_with_introns(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor,
                           fontproperties=None):
    """
    draws a gene using different styles
    """
    from matplotlib.patches import Polygon

    (chrom, start, end, name, strand, thick_start, thick_end, rgb_, block_count,
     block_sizes, block_starts) = get_gene_parts(fields, bed_properties)
    # draw a line from the start until the end of the gene
    if block_count == 0 and thick_start == start and thick_end == end:
        draw_gene_simple(ax, fields, ypos,  bed_properties, small_relative, rgb, edgecolor)
        return

    ax.plot([start, end], [ypos+50, ypos+50], 'black', linewidth=0.5, zorder=-1)
    if strand == '+':
        x1 = thick_start
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        -----------------\
        -----------------/

        """
        start_box = [(start, y0), (start, y1), (x1, y1), (x1, y0)]
        end_box = [(thick_end, y0), (thick_end, y1), (end - small_relative, y1), (end, y0 + 50), (end-small_relative, y0)]

    elif strand == '-':
        x0 = start + min(small_relative, thick_start - start)
        y0 = ypos
        y1 = ypos + 100
        """
        The vertices correspond to 5 points along the path of a form like the following,
        starting in the lower left corner and progressing in a clock wise manner.

        /--------___---------_
        \--------   ----------
        """
        start_box = [(x0, y0), (start, y0 + 50), (x0, y1), (thick_start, y1), (thick_start, y0)]
        end_box = [(thick_end, y0), (thick_end, y1), (end, y1),  (end, y0)]

    for idx in range(0, block_count):
        x0 = start + block_starts[idx]
        x1 = x0 + block_sizes[idx]
        if x1 < thick_start or x0 > thick_end:
            y0 = ypos + 25
            y1 = ypos + 75
        else:
            y0 = ypos
            y1 = ypos + 100

        if x0 < thick_start < x1:
            vertices = ([(x0, ypos+25), (x0, ypos+75), (thick_start, ypos+75), (thick_start, ypos+100),
                         (thick_start, ypos+100), (x1, ypos+100), (x1, ypos),
                         (thick_start, ypos), (thick_start, ypos+25)])

        elif x0 < thick_end < x1:
            vertices = ([(x0, ypos), (x0, ypos+100), (thick_end, ypos+100), (thick_end, ypos+75),
                         (x1, ypos+75), (x1, ypos+25), (thick_end, ypos+25), (thick_end, ypos)])
        else:
            vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             linewidth=0.1,
                             edgecolor='none',
                             facecolor=rgb))

        if idx < block_count - 1:
            intron_length = block_starts[idx+1] - (block_starts[idx] + block_sizes[idx])
            marker = 5 if strand == '+' else 4
            if intron_length > 3*small_relative:
                pos = np.arange(x1 + 1*small_relative, x1 + intron_length + small_relative, int(2 * small_relative))
                ax.plot(pos, np.zeros(len(pos))+ypos+50, '.', marker=marker,
                        fillstyle='none', color='blue', markersize=3)

            elif intron_length > small_relative:
                intron_center = x1 + int(intron_length)/2
                ax.plot([intron_center], [ypos+50], '.', marker=5,
                        fillstyle='none', color='blue', markersize=3)

    if 'labels' in bed_properties and bed_properties['labels'] != 'off':
        center = start + float(end - start)/2
        ax.text(center, ypos + 125, fields[3],
                horizontalalignment='center', verticalalignment='top', fontproperties=fontproperties)


def plot_bed(ax, label_ax, bed_properties, region):

    from matplotlib.patches import Rectangle

    file_h = open(bed_properties['file'], 'r')
    chrom_region, start_region, end_region = region
    counter = 0
    small_relative = 0.003 * (end_region-start_region)
    prev_start = -1
    prev_line = None
    prev_chrom = None
    max_num_row = 1
    region_intervals = IntervalTree()
    colormap = None
    # check if the color given is a color map
    from matplotlib import cm
    color_options = [m for m in cm.datad]

    ax.set_frame_on(False)
    # to improve the visualization of the genes
    # it is good to have an estimation of the label
    # length. In the following code I try to get
    # the length of a 'w'.

    from matplotlib import font_manager
    if 'fontsize' in bed_properties:
        fp = font_manager.FontProperties(size=bed_properties['fontsize'])
    else:
        fp = font_manager.FontProperties()

    # to avoid overlaping gene labels, the size of a 'w' is estimated
    if 'type' in bed_properties and bed_properties['type'] == 'genes' \
            and 'labels' in bed_properties and bed_properties['labels'] != 'off':
        t = matplotlib.textpath.TextPath((0, 0), 'w', prop=fp)
        len_w = t.get_extents().width * 300
    else:
        len_w = 1

    if 'color' not in bed_properties:
        bed_properties['color'] = DEFAULT_BED_COLOR

    if bed_properties['color'] in color_options:
        import matplotlib as mpl
        norm = mpl.colors.Normalize(vmin=bed_properties['min_value'],
                                    vmax=bed_properties['max_value'])
        cmap = cm.get_cmap(bed_properties['color'])
        colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    for line in file_h.readlines():
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')

        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(detail)
            exit(msg)
        start = int(start)
        end = int(end)
        if prev_chrom is None:
            prev_chrom = chrom
        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}\n{} ".format(prev_line, line)

        prev_start = start
        prev_line = line
        assert start >= 0, \
            "negative values found in bed file in line {}".format(line)

        if chrom == chrom_region and end > start_region and \
                start < end_region:
            # Rectangle(xy, width, height, angle=0.0, **kwargs)
            counter += 1
            try:
                strand = fields[5]
            except IndexError:
                strand = "."
            try:
                rgb = fields[8]
                rgb = rgb.split(',')
                rgb[2]  # check that rgb has three elements
                rgb = [float(x)/255 for x in rgb]
                edgecolor = bed_properties['color']
            except IndexError:
                if colormap:
                    rgb = colormap.to_rgba(float(fields[4]))
                    edgecolor = colormap.to_rgba(float(fields[4]))
                else:
                    rgb = bed_properties['color']
                    edgecolor = bed_properties['color']

            except Exception as e:
                exit("Error occurred: {}".format(e))

            if 'type' in bed_properties and \
                    bed_properties['type'] == 'domain':
                ypos = 100 if counter % 2 == 0 else 1
                ax.add_patch(Rectangle(
                        (start, ypos),
                        end-start,
                        100, edgecolor='black',
                        facecolor=bed_properties['color']))

            # check for overlapping features
            match = region_intervals.find(start, end)
            if len(match) == 0:
                min_free_row = 0
            else:
                rows_used = np.zeros(max_num_row + 2)
                for x in match:
                    rows_used[x.value] = 1
                min_free_row = min(np.flatnonzero(rows_used == 0))

            if 'type' in bed_properties and bed_properties['type'] == 'genes' and end - start < len(fields[3]) * len_w:
                region_intervals.add_interval(Interval(start, start + (len(fields[3]) * len_w), min_free_row))
            else:
                region_intervals.add_interval(Interval(start, end+small_relative, min_free_row))
            if min_free_row > max_num_row:
                max_num_row = min_free_row
            if min_free_row > 0:
                # this means we are plotting 
                # on top of previous genes
                ypos = min_free_row * 230
            else:
                ypos = 0

            if 'display' in bed_properties and \
                    bed_properties['display'] == 'collapsed':
                ypos = 0

            # give direction to genes
            if 'type' in bed_properties and \
                    bed_properties['type'] == 'genes' and \
                    strand in ['+', '-']:
                draw_gene_with_introns(ax, fields, ypos, bed_properties, small_relative, rgb, edgecolor,
                                       fontproperties=fp)

            else:
                ax.add_patch(Rectangle(
                        (start, ypos), end-start,
                        100, edgecolor=edgecolor,
                        facecolor=rgb))

    if counter == 0:
        sys.stderr.write("*Warning* No intervals were found for file {} \n"
                         "in section '{}' for the interval plotted ({}:{}-{}).\n"
                         "".format(bed_properties['file'],
                                   bed_properties['section_name'],
                                   chrom_region, start_region, end_region))
    if 'type' in bed_properties and bed_properties['type'] == 'domain':
        ax.set_ylim(-5, 205)
    elif 'display' in bed_properties and bed_properties['display'] == 'collapsed':
        ax.set_ylim(-5, 105)
    else:
        ax.set_ylim((max_num_row+1)*230, -25)

    ax.set_xlim(region[1], region[2])

    label_ax.text(0.15, 1.0, bed_properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='top', transform=label_ax.transAxes)


def plot_bedgraph(ax, label_ax, bedgraph_properties, region):
    file_h = open(bedgraph_properties['file'], 'r')
    chrom_region, start_region, end_region = region
    score_list = []
    pos_list = []

    for line in file_h.readlines():
        if line.startswith('browser') or line.startswith('track'):
            continue
        chrom, start, end, score = line.split('\t')
        start = int(start)
        end = int(end)
        if chrom == chrom_region and start_region - 100000 <= start and \
                end_region + 100000 >= end:
            score_list.append(float(score))
            pos_list.append(start + (end - start)/2)

    if 'color' not in bedgraph_properties:
        bedgraph_properties['color'] = DEFAULT_BEDGRAPH_COLOR

    if 'extra' in bedgraph_properties and \
            bedgraph_properties['extra'][0] == '4C':
        # draw a vertical line for each fragment region center
        ax.fill_between(pos_list, score_list,
                        facecolor=bedgraph_properties['color'],
                        edgecolor='none')
        ax.vlines(pos_list, [0], score_list, color='olive', linewidth=0.5)
        ax.plot(pos_list, score_list, '-', color='slateblue', linewidth=0.7)
    else:
        try:
            ax.fill_between(pos_list, score_list,
                            facecolor=bedgraph_properties['color'])
        except ValueError:
            exit("Invalid color {} for {}".format(bedgraph_properties['color'], bedgraph_properties['file']))
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlim(region[1], region[2])

    ymin, ymax = ax.get_ylim()
    if 'max_value' in bedgraph_properties and bedgraph_properties['max_value'] != 'auto':
        ymax = bedgraph_properties['max_value']
    if 'min_value' in bedgraph_properties and bedgraph_properties['min_value'] != 'auto':
        ymin = bedgraph_properties['min_value']

    if float(ymax) % 1 == 0:
        ymax_print = int(ymax)
    else:
        ymax_print = "{:.1f}".format(ymax)
    ax.set_ylim(ymin, ymax)
    ydelta = ymax - ymin
    small_x = 0.01 * (end_region - start_region)

    if 'show data range' in bedgraph_properties and \
            bedgraph_properties['show data range'] == 'no':
        pass
    else:
        # by default show the data range
        ax.text(start_region-small_x, ymax - ydelta * 0.2,
                "[{}-{}]".format(ymin, ymax_print),
                horizontalalignment='left', size='small',
                verticalalignment='bottom')

    label_ax.text(0.15, 0, bedgraph_properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom', transform=label_ax.transAxes)

    """
    if 'extra' in bedgraph_properties :

        ticks = ax.get_xticks()
        labels = ["{:.1f} Mbp".format((x / 1e6))
                  for x in ticks]

        ax.set_xticklabels(labels, size='large')
        ax.axes.get_xaxis().set_visible(True)
    """


def plot_bedgraph_matrix(ax, label_ax, properties, region):
    """
    Plots a bedgraph matrix file, that instead of having
    a single value per bin, it has several values.

    :param ax:
    :param label_ax:
    :param properties:
    :param region:
    :return:
    """

    fh = open(properties['file'], 'r')
    chrom_region, start_region, end_region = region
    start_pos = []
    matrix_rows = []
    for line in fh:
        line = line.strip()
        region = line.split('\t')
        chrom = region[0]
        start = int(region[1])
        end = int(region[2])
        if chrom == chrom_region and start_region - 100000 <= start and end_region + 100000 >= end:
            start_pos.append(start)
            matrix_rows.append(np.fromiter(region[3:], np.float))

    matrix = np.vstack(matrix_rows).T
    if 'orientation' in properties and properties['orientation'] == 'inverted':
        matrix = np.flipud(matrix)

    vmin = None
    vmax = None
    if 'max_value' in properties and properties['max_value'] != 'auto':
        vmax = properties['max_value']

    if 'min_value' in properties and properties['min_value'] != 'auto':
        vmin = properties['min_value']

    if 'type' in properties and properties['type'] == 'lines':
        for row in matrix:
            ax.plot(start_pos, row)
        ax.plot(start_pos, matrix.mean(axis=0), "--")
    else:
        x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
        img = ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading='gouraud')
        img.set_rasterized(True)
    ax.set_xlim(start_region, end_region)
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    label_ax.text(0.15, 0, properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom', transform=label_ax.transAxes)


def plot_bigwig(ax, label_ax, bigwig_properties, region):

    from bx.bbi.bigwig_file import BigWigFile
    bw = BigWigFile(open(bigwig_properties['file'], 'r'))
    chrom, region_start, region_end = region
    # compute the score in bins of 10000 SLOW
#    score = np.array([bw.query(region[0], x, x+10000,1)[0]['mean']
#                      for x in range(region[1], region[2], 10000)])

    num_bins = 700
    if 'number of bins' in bigwig_properties:
        try:
            num_bins = int(bigwig_properties['number of bins'])
        except TypeError:
            exit("'number of bins' value: {} for bigwig file {} "
                 "is not valid.".format(bigwig_properties['number of bins'],
                                        bigwig_properties['file']))

    scores = bw.get_as_array(chrom, region_start, region_end)
    if scores is None:
        # usually bw returns None when the chromosome
        # is not found. So, we try either
        # removing or appending 'chr'
        if chrom.startswith('chr'):
            scores = bw.get_as_array(chrom[3:] + chrom, region_start, region_end)
        else:
            scores = bw.get_as_array('chr' + chrom, region_start, region_end)

        if scores is None:
            exit("Can not read region {}:{}-{} from bigwig file:\n\n"
                 "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                 "and that the region is valid".format(chrom, region_start, region_end,
                                                       bigwig_properties['file']))

    if 'nans to zeros' in bigwig_properties and bigwig_properties['nans to zeros'] is True:
        scores[np.isnan(scores)] = 0

    scores = np.ma.masked_invalid(scores)

    if 'color' not in bigwig_properties:
        bigwig_properties['color'] = DEFAULT_BIGWIG_COLOR

    if region_end - region_start < 2e6:
        if scores is None:
            chrom = chrom.replace('chr', '')
            scores = np.ma.masked_invalid(bw.get_as_array(chrom, region_start, region_end))
        if scores is None:
            sys.stderr.write("could not find values for region {}\n".format(region))

        else:
            lins = np.linspace(0, len(scores), num_bins).astype(int)
            scores_per_bin = [np.mean(scores[lins[x]:lins[x+1]]) for x in range(len(lins)-1)]
            _x = lins + region_start
            x_values = [float(_x[x] + _x[x+1])/2 for x in range(len(lins)-1)]
            ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                            color=bigwig_properties['color'],
                            facecolor=bigwig_properties['color'])

    else:
        # this method produces shifted regions. It is not clear to me why this happens.
        # Thus I only activate the faster but shifted method for large regions
        # when the previous method would be to slow
        score = bw.query(chrom, region_start, region_end, num_bins)
        if score is None:
            chrom = chrom.replace('chr', '')
            score = bw.query(chrom, region_start, region_end, num_bins)
        if score is None:
            sys.stderr.write("could not find values for region {}\n".format(region))
        else:
            score = [x['mean'] for x in score]
            x_values = np.linspace(region_start, region_end, num_bins)
            ax.fill_between(x_values, score, linewidth=0.1,
                            color=bigwig_properties['color'],
                            facecolor=bigwig_properties['color'], zorder=1)

    ax.set_xlim(region[1], region[2])
    ymin, ymax = ax.get_ylim()
    if 'max_value' in bigwig_properties and ['max_value'] != 'auto':
        ymax = bigwig_properties['max_value']
    if 'min_value' in bigwig_properties and bigwig_properties['min_value'] != 'auto':
        ymin = bigwig_properties['min_value']

    if 'orientation' in bigwig_properties and bigwig_properties['orientation'] == 'inverted':

        ax.set_ylim(ymax, ymin)
    else:
        ax.set_ylim(ymin, ymax)

#    ax.set_yticks([ymax])
    ydelta = ymax - ymin

#    ax.set_yticklabels(["{}-{}".format(int(ymin), int(ymax))], size='large')
    # set min max
    if float(ymax) % 1 == 0:
        ymax_print = int(ymax)
    else:
        ymax_print = "{:.1f}".format(ymax)
    small_x = 0.01 * (region_end - region_start)
    if 'show data range' in bigwig_properties and bigwig_properties['show data range'] == 'no':
        pass
    else:
        # by default show the data range
        ax.text(region_start-small_x, ymax - ydelta * 0.2,
                "[{}-{}]".format(int(ymin), ymax_print),
                horizontalalignment='left', size='small',
                verticalalignment='bottom')

    """
    ax.text(region_end, ymax - ydelta * 0.2, bigwig_properties['title'],
            horizontalalignment='right', size='large',
            verticalalignment='bottom')

    """
    label_ax.text(0.15, 0, bigwig_properties['title'],
                  horizontalalignment='left', size='large',
                  verticalalignment='bottom')


def plot_vlines(vlines_list, vlines_file, axis_list, region):
    """
    Plots dotted lines from the top of the first plot to the bottom
    of the last plot at the specified positions.

    :param vlines_list: list of positions
    :param vlines_file: bed file
    :param axis_list: list of plotted axis
    :param region: tuple containing the region to plot
    :return: None
    """
    chrom_region, start_region, end_region = region
    if vlines_list is None:
        vlines_list = []

    if vlines_file:
        file_h = open(vlines_file, 'r')
        for line in file_h.readlines():
            if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                continue
            # consider only the first two columns (chrom and start)
            chrom, start = line.strip().split('\t')[0:2]
            start = int(start)
            if chrom == chrom_region and \
               start_region <= start and \
               end_region >= start + 1:
                vlines_list.append(start)

    from matplotlib.patches import ConnectionPatch
    a_ymax = axis_list[0].get_ylim()[1]
    b_ymin = axis_list[-1].get_ylim()[0]

    for start_pos in vlines_list:
        con = ConnectionPatch(xyA=(start_pos, a_ymax),
                              xyB=(start_pos, b_ymin),
                              coordsA="data", coordsB="data",
                              axesA=axis_list[0],
                              axesB=axis_list[-1],
                              arrowstyle="-",
                              linestyle='dashed',
                              linewidth=0.5,
                              zorder=100)
        axis_list[0].add_artist(con)


def get_region(region_string):
    """
    splits a region string into
    a chrom, start_region, end_region tuple
    The region_string format is chr:start-end
    """
    if region_string:
        region_string = region_string.translate(
            None, ",.;|!{}()").replace("-", ":")
        region = region_string.split(":")
        chrom = region[0]
        try:
            region_start = int(region[1])
        except IndexError:
            region_start = 0
        try:
            region_end = int(region[2])
        except IndexError:
            region_end = 1e15  # a huge number

        return chrom, region_start, region_end


class MultiDict(OrderedDict):
    """
    Class to allow identically named
    sections in configuration file
    by appending the section number as
    for example:
    1. section name
    """
    _unique = 0

    def __setitem__(self, key, val):
        if isinstance(val, OrderedDict):
            self._unique += 1
            key = "{}. {}".format(str(self._unique), key)
        OrderedDict.__setitem__(self, key, val)


def parse_tracks(tracks_file):
    """
    Parses a configuration file

    :param tracks_file: file path containing the track configuration
    :return: array of dictionaries and vlines_file. One dictionary per track
    """
    from ConfigParser import SafeConfigParser
    from ast import literal_eval
    parser = SafeConfigParser(None, MultiDict)
    parser.readfp(tracks_file)

    track_list = []
    vlines_file = None
    for section_name in parser.sections():
        track_options = dict({"section_name": section_name})
        if section_name.endswith('spacer'):
            track_options['spacer'] = True
        elif section_name.endswith('x-axis'):
            track_options['x-axis'] = True
        for name, value in parser.items(section_name):
            if name in ['max_value', 'min_value', 'depth', 'width'] and value != 'auto':
                track_options[name] = literal_eval(value)
            else:
                track_options[name] = value

        if 'type' in track_options and track_options['type'] == 'vlines':
            vlines_file = track_options['file']
        else:
            track_list.append(track_options)

    return track_list, vlines_file


def check_file_exists(track_dict):
    """
    Checks if a file or list of files exists
    :param track_dict: dictionary of track values. Should contain
                        a 'file' key containing the path of the file
                        or files to be checked separated by space
                        For example: file1 file2 file3
    :return: None
    """
    file_names = [x for x in track_dict['file'].split(" ") if x != '']
    for file_name in file_names:
        try:
            open(file_name, 'r').close()
        except IOError:
            sys.stderr.write("\n*ERROR*\nFile in section [{}] "
                             "not found:\n{}\n\n".format(track_dict['section_name'],
                                                         file_name))
            exit(1)


def guess_filetype(track_dict):
    """

    :param track_dict: dictionary of track values with the 'file' key
                containing a string path of the file or files. Only the ending
                 of the last file is used in case when there are more files
    :return: string file type detected
    """
    file = track_dict['file'].strip()
    file_type = None

    if file.endswith(".bed"):
        file_type = 'bed'
    elif file.endswith(".npz"):
        file_type = 'hic_matrix'
    elif file.endswith(".bw"):
        file_type = 'bigwig'
    elif file.endswith(".bg"):
        file_type = 'bedgraph'
    elif file.endswith(".bm"):
        file_type = 'bedgraph_matrix'
    else:
        exit("Section [{}]: can not identify file type. Please specify "
             "the file_type for {}".format(track_dict['section_name'],
                                           file))

    return file_type


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def main(args=None):

    args = parse_arguments().parse_args(args)

    region = get_region(args.region)
    chrom, region_start, region_end = region
    if region_end <= region_start:
        exit("Please check that the region end is larger than the region start.\n"
             "Values given:\nstart: {}\nend: {}\n".format(region_start, region_end))

    track_properties, vlines_file = parse_tracks(args.tracks)

    # prepare layout based on the tracks given.
    # The main purpose of the following loop is
    # to get the height of each of the tracks
    track_height = []
    for track_dict in track_properties:
        warn = None
        height = DEFAULT_TRACK_HEIGHT
        if 'file' in track_dict and track_dict['file'] != '':
            if 'file_type' not in track_dict:
                track_dict['file_type'] = guess_filetype(track_dict)
            check_file_exists(track_dict)

            #  set some default values
            if 'title' not in track_dict:
                warn = "\ntitle not set for 'section {}'\n".format(track_dict['section_name'])
                track_dict['title'] = ''
            if warn:
                sys.stderr.write(warn)

        if 'width' in track_dict:
            height = track_dict['width']

        # compute the height of a Hi-C track
        # based on the depth such that the
        # resulting plot appears proportional
        #
        #      /|\
        #     / | \
        #    /  |d \   d is the depth that we want to be proportional
        #   /   |   \  when plotted in the figure
        # ------------------
        #   region len
        #
        # d (in cm) =  depth (in bp) * width (in cm) / region len (in bp)

        elif 'depth' in track_dict and track_dict['file_type'] == 'hic_matrix':
            # to compute the actual width of the figure the margins and the region
            # set for the legends have to be considered
            # DEFAULT_MARGINS[1] - DEFAULT_MARGINS[0] is the proportion of plotting area

            hic_width = args.width * (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) * DEFAULT_WIDTH_RATIOS[0]
            scale_factor = 0.6  # the scale factor is to obtain a pleasing result.
            height = scale_factor * track_dict['depth'] * hic_width / (region_end - region_start)

#            height = (track_dict['depth'] * args.width /
#                      (1.8*(region_end - region_start)))

        track_height.append(height)

    if args.height:
        fig_height = args.height
    else:
        fig_height = sum(track_height)

    sys.stderr.write("Figure size in cm is {} x {}. Dpi is set to {}\n".format(args.width, fig_height, args.dpi))
    fig = plt.figure(figsize=cm2inch(args.width, fig_height))
    fig.suptitle(args.title)

    if args.fontSize:
        fontsize = args.fontSize
    else:
        fontsize = float(args.width) * 0.4

    font = {'size': fontsize}
    matplotlib.rc('font', **font)

    grids = gridspec.GridSpec(len(track_height), 2,
                              height_ratios=track_height,
                              width_ratios=[1, 0.05])

    # iterate again to plot each track
    axis_list = []
    for idx, properties in enumerate(track_properties):
        if 'spacer' in properties:
            continue
        axis = axisartist.Subplot(fig, grids[idx, 0])
        fig.add_subplot(axis)
        axis.axis[:].set_visible(False)
        # to make the background transparent
        axis.patch.set_visible(False)

        if 'x-axis' in properties:
            # ideally the axisartis would allow
            # to have a floating axis but this is not
            # working
            plot_x_axis(axis, region, properties)
            continue
        else:
            label_axis = plt.subplot(grids[idx, 1])
            label_axis.set_axis_off()
        if properties['file_type'] == 'bed':
            plot_bed(axis, label_axis, properties, region)
        elif properties and properties['file_type'] == 'dendrogram':
            plot_dendrogram(axis, label_axis, properties, region)
        elif properties['file_type'] == 'bedgraph':
            plot_bedgraph(axis, label_axis, properties, region)
        elif properties['file_type'] == 'bigwig':
            plot_bigwig(axis, label_axis, properties, region)
        elif properties['file_type'] == 'bedgraph_matrix':
            plot_bedgraph_matrix(axis, label_axis, properties, region)
        elif properties['file_type'] == 'hic_matrix':
            if 'type' in properties:
                if properties['type'] == 'interaction':
                    plot_interactions(axis, properties, region)
                elif properties['type'] == 'arcplot':
                    plot_hic_arcs(axis, label_axis, properties, region)
                else:
                    exit("Hi-C plot type invalid ({}) for {}".format(properties['type'],
                                                                     properties['file']))
            else:
                # to avoid the color bar to span all the
                # width of the axis I pass two axes
                # to plot_matrix
                cbar_axis = label_axis
                label_axis = plt.subplot(grids[idx, 1])
                label_axis.set_axis_off()
                plot_matrix(axis, label_axis, cbar_axis,
                            properties, region)
        axis_list.append(axis)
    if args.vlines or vlines_file:
        plot_vlines(args.vlines, vlines_file, axis_list, region)


    plt.subplots_adjust(wspace=0, hspace=0.1,
                        left=DEFAULT_MARGINS['left'],
                        right=DEFAULT_MARGINS['right'],
                        bottom=DEFAULT_MARGINS['bottom'],
                        top=DEFAULT_MARGINS['top'])

    plt.savefig(args.outFileName.name, dpi=args.dpi)
