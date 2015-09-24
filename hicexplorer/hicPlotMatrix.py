import sys
import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import writableFile
from hicexplorer._version import __version__
from scipy.sparse import tril
import numpy as np

import argparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

def parseArguments(args=None):
    parser=argparse.ArgumentParser(description='Plots a heatmap of a HiC matrix')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='path of the Hi-C matrix to plot',
                        required=True)

    parser.add_argument('--title', '-t',
                        help='Plot title')

    parser.add_argument('--scoreName', '-s',
                        help='Score name')

    parser.add_argument('--outFileName', '-out',
                        help='File name to save the image. ',
                        type=writableFile,
                        required=True)

    parser.add_argument('--perChromosome',
                        help='Instead of plotting the whole matrix, '
                        'each chromosome is plotted next to the other. '
                        'This parameter is not compatible with --region',
                        action='store_true')

    parser.add_argument('--clearMaskedBins',
                        help='if set, masked bins are removed from the matrix',
                        action='store_true')

    parser.add_argument('--whatToShow',
                        help='Options are: "heatmap", "3D", and "both". '
                        'Default is both',
                        default="both",
                        choices=["heatmap", "3D", "both"])

    parser.add_argument('--chromosomeOrder',
                        help='Chromosomes and order in which the '
                        'chromosomes should be plotted. This option '
                        'overrides --region and --region2 ',
                        nargs='+')

    parser.add_argument('--region',
                        help='Plot only this region. The format is '
                        'chr:start-end The plotted region contains '
                        'the main diagonal and is symmetric unless '
                        ' --region2 is given'
                        )

    parser.add_argument('--region2',
                        help='If given, then only the region defined by '
                        '--region and --region2 is given. The format '
                        'is the same as --region1'
                        )

    parser.add_argument('--log1p',
                        help='Plot the log1p of the matrix values.',
                        action='store_true')

    parser.add_argument('--log',
                        help='Plot the *MINUS* log of the matrix values.',
                        action='store_true')


    color_options = "', '".join([m for m in cm.datad
                                 if not m.endswith('_r')])

    parser.add_argument('--colorMap',
                        help='Color map to use for the heatmap. Available '
                        'values can be seen here: '
                        'http://www.astro.lsa.umich.edu/~msshin/science/code/'
                        'matplotlib_cm/ The available options are: \'' +
                        color_options + '\'',
                        default='jet')

    parser.add_argument('--vMin',
                        help='vMin',
                        type=float,
                        default=None)

    parser.add_argument('--vMax',
                        help='vMax',
                        type=float,
                        default=None)

    parser.add_argument('--zMax',
                        help='zMax for 3D plot',
                        type=float,
                        default=None)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def plotHeatmap3D(ma, fig, position, args, cmap):
    axHeat1 = fig.add_axes( position, projection='3d' )
    axHeat1.view_init(50, -45)
    axHeat1.margins(0)

    X, Y = np.meshgrid(range(ma.shape[0]), range(ma.shape[0]))
    ma = np.ma.array( ma, mask=np.isnan(ma) )
    Z = ma.copy()

    Z[np.where(np.tril(Z)==0)] = np.nan

    axHeat1.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                         vmax=args.vMax, vmin=args.vMin)
    if args.zMax:
        axHeat1.set_zlim(0, args.zMax)


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, figWidth, cmap):

    axHeat2 = fig.add_axes(position, title=args.title)
    norm = None
    if args.log1p:
        ma += 1
        norm = LogNorm()

    img3 = axHeat2.imshow(ma, 
                          interpolation='nearest',
#                          interpolation='spline16',
                          vmax=args.vMax, vmin=args.vMin, cmap=cmap,
                          norm=norm
                          )

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axHeat2)
    cax = divider.append_axes("right", size="2.5%", pad=0.09)
    cbar = fig.colorbar(img3, cax=cax)
    cbar.solids.set_edgecolor("face") # to avoid white lines in the color bar in pdf plots
    cbar.ax.set_ylabel(args.scoreName, rotation=270, size=8)

    if args.region:
        """
        ticks = axHeat2.get_xticks()
        labels = ["{:.3f} Mbp".format((x / 1e6))
                  for x in ticks]
        axHeat2.get_xaxis().set_tick_params(
            which='both',
            bottom='on')
        axHeat2.set_xticklabels(labels, size='small')
        """

        axHeat2.set_xticks([0, ma.shape[0]])
        axHeat2.set_xticklabels([args.region[1], args.region[2]], size=4, rotation=90)
        axHeat2.set_axis_off()
    else:
        ticks = [ int(pos[0] + (pos[1]-pos[0]) /2 ) for pos in chrBinBoundaries.values()]
        labels = chrBinBoundaries.keys()
        axHeat2.set_xticks(ticks)
        if len(labels)> 20:
            axHeat2.set_xticklabels(labels, size=4, rotation=90)
        else:
            axHeat2.set_xticklabels(labels, size=8)


    axHeat2.get_yaxis().set_visible(False)

def translate_region(region_string):
    """
    Takes an string and returns a list
    of chrom, start, end.
    If the region string only contains
    the chrom, then start and end
    are set to a 0 and 1e15
    """

    region_string = region_string.translate(None, ",.;|!{}()").replace("-", ":")
    fields = region_string.split(":")
    chrom = fields[0]
    try:
        region_start = int(fields[1])
    except IndexError:
        region_start = 0
    try:
        region_end = int(fields[2])
    except IndexError:
        region_end = 1e15 # vert large number

    return chrom, region_start, region_end


def plotPerChr(hic_matrix, cmap, args):
    """
    plots each chromosome individually, one after the other
    in one row. scale bar is added at the end
    """
    chromosomes = hic_matrix.getChrNames()
    width_ratios = [1] * len(chromosomes) + [0.05]
    grids = gridspec.GridSpec(1, len(chromosomes) + 1,
                              width_ratios=width_ratios)

    fig_height = 6
    fig_width = sum((np.array(width_ratios)+0.95) * fig_height)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    norm = None
    if args.log1p:
        norm = LogNorm()
    chrom, start, end, _ = zip(*hic_matrix.cut_intervals)

    for idx, chrname in enumerate(chromosomes):
        axis = plt.subplot(grids[idx])
        axis.set_title(chrname)
        chrom_range = hic_matrix.getChrBinRange(chrname)
        mat = hic_matrix.matrix[chrom_range[0]:chrom_range[1],
                                chrom_range[0]:chrom_range[1]].todense().astype(float)

        if args.log1p:
            mat += 1
        img = axis.imshow(mat, aspect='auto',
                          interpolation='spline16',
#                          interpolation='nearest',
                          vmax=args.vMax, vmin=args.vMin, cmap=cmap,
                          norm=norm,
                          extent=[start[chrom_range[0]], end[chrom_range[1] - 1],
                                  end[chrom_range[1]-1], start[chrom_range[0]]])

        xticks = axis.get_xticks()
        xlabels = ["{:.0f}".format(int(x) / 1e6)
                  for x in xticks]
        print xlabels

        axis.set_xticklabels(xlabels, size='small')
        yticks = axis.get_yticks()

        ylabels = ["{:.0f}".format(int(x) / 1e6)
                  for x in yticks]

        axis.get_xaxis().set_tick_params(
            which='both',
            bottom='on',
            direction='out')

        axis.set_yticklabels(ylabels, size='small')

    cbar3 = plt.subplot(grids[-1])
    cbar = fig.colorbar(img, cax=cbar3)
    cbar.solids.set_edgecolor("face") # to avoid white lines in
                                      # the color bar in pdf plots
    cbar.ax.set_ylabel(args.scoreName, rotation=270, labelpad=20)

def main():
    args = parseArguments().parse_args()

    ma = HiCMatrix.hiCMatrix(args.matrix)
    if args.perChromosome and args.region:
        sys.stderr.write('ERROR, choose from the option '
                         '--perChromosome or --region, the two '
                         'options at the same time are not '
                         'compatible.')
    if args.chromosomeOrder:
        args.region = None
        args.region2 = None
        ma.reorderChromosomes(args.chromosomeOrder)

    print ma.interval_trees.keys()
    ma.restoreMaskedBins()

    if args.clearMaskedBins:
        ma.maskBins(ma.nan_bins)

    #ma.matrix.data[np.isnan(ma.matrix.data)] = 0
    """
    #temp
    sys.stderr.write("\n\n**Warning, normalizing values from "
                     "0 to 1 for comparizon\n\n")

    ma.diagflat(0)
    ma.matrix.data = ma.matrix.data / ma.matrix.data.max()
    """
    sys.stderr.write("min: {}, max: {}\n".format(ma.matrix.data.min(), ma.matrix.data.max()))
    if args.region:
        chrom, region_start, region_end = translate_region(args.region)
        args.region = [chrom, region_start, region_end]
        idx1 = [idx for idx, x in enumerate(ma.cut_intervals) if \
                   x[0] == chrom and x[1] >= region_start and \
                   x[2] < region_end]
        if args.region2:
            chrom2, region_start2, region_end2 = \
                translate_region(args.region2)
            idx2 = [idx for idx, x in enumerate(ma.cut_intervals) if \
                       x[0] == chrom2 and x[1] >= region_start2 and \
                       x[2] < region_end2]

            # select only relevant part
            matrix = ma.matrix[idx1, :][:, idx2].todense().astype(float)
        else:
            matrix = ma.matrix[idx1, :][:, idx1].todense().astype(float)

    else:
        matrix = ma.getMatrix().astype(float)

    cmap = cm.get_cmap(args.colorMap)
    sys.stderr.write("Nan values set to black\n")
    cmap.set_bad('black')

    if args.perChromosome:
        plotPerChr(ma, cmap, args)

    else:
        # set as nan regions the nan bins

        fig_height = 6
        height = 5.0/fig_height
        if args.whatToShow == 'both':
            fig_width = 11
            width = 4.9/fig_width
            left_margin = (1.0-2*width) * 0.35
        else:
            fig_width = 7
            width = 4.4/fig_width
            left_margin = (1.0-width) * 0.35


        fig = plt.figure(figsize=(fig_width, fig_height), dpi=720)
        bottom = 0.6 /fig_height

        if args.log:
            mask = matrix == 0
            matrix[mask] = matrix[mask == False].min()
    #        matrix = -1 *  np.log(matrix)
            matrix = np.log(matrix)

        if args.whatToShow == '3D':
            position = [left_margin, bottom, width*1.2, height*1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

        if args.whatToShow == 'both':
            position = [left_margin, bottom, width*1.2, height*1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

            left_margin2 = 5.5/fig_width + left_margin
            position = [left_margin2, bottom, width, height]
            plotHeatmap(matrix, ma.chrBinBoundaries, fig, position, args,
                        fig_width, cmap)

        if args.whatToShow == 'heatmap':
            position = [left_margin, bottom, width, height]
            plotHeatmap(matrix, ma.chrBinBoundaries, fig, position,
                        args, fig_width, cmap)

    plt.savefig(args.outFileName, dpi=320)
