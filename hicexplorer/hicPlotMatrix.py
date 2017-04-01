import sys
import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import writableFile
from hicexplorer._version import __version__
import numpy as np

import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Creates a Heatmap of a HiC matrix')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='Path of the Hi-C matrix to plot',
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
                        'Default is heatmap',
                        default="heatmap",
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

    parser.add_argument('--colorMap',
                        help='Color map to use for the heatmap. Available '
                        'values can be seen here: '
                        'http://matplotlib.org/examples/color/colormaps_reference.html',
                        default='RdYlBu_r')

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

    parser.add_argument('--dpi',
                        help='Resolution for the image in case the'
                             'ouput is a raster graphics image (e.g png, jpg)',
                        type=int,
                        default=72)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def plotHeatmap3D(ma, fig, position, args, cmap):
    axHeat1 = fig.add_axes(position, projection='3d')
    axHeat1.view_init(50, -45)
    axHeat1.margins(0)

    X, Y = np.meshgrid(range(ma.shape[0]), range(ma.shape[0]))
    ma = np.ma.array(ma, mask=np.isnan(ma))
    Z = ma.copy()

    Z[np.where(np.tril(Z) == 0)] = np.nan

    axHeat1.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, cmap=cmap,
                         vmax=args.vMax, vmin=args.vMin)
    if args.zMax:
        axHeat1.set_zlim(0, args.zMax)


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, figWidth, cmap):

    axHeat2 = fig.add_axes(position)
    if args.title:
        axHeat2.set_title(args.title)
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

    img3.set_rasterized(True)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axHeat2)
    cax = divider.append_axes("right", size="2.5%", pad=0.09)
    cbar = fig.colorbar(img3, cax=cax)
    cbar.solids.set_edgecolor("face")  # to avoid white lines in the color bar in pdf plots
    if args.scoreName:
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
        ticks = [int(pos[0] + (pos[1] - pos[0]) / 2) for pos in chrBinBoundaries.values()]
        labels = chrBinBoundaries.keys()
        axHeat2.set_xticks(ticks)
        if len(labels) > 20:
            axHeat2.set_xticklabels(labels, size=4, rotation=90)
        else:
            axHeat2.set_xticklabels(labels, size=8)

    axHeat2.get_yaxis().set_visible(False)


def plotHeatmap_region(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                       ylabel=None, start_pos=None, start_pos2=None):

    axHeat2 = fig.add_axes(position)
    if args.title:
        axHeat2.set_title(args.title)
    norm = None
    if args.log1p:
        ma += 1
        norm = LogNorm()

    if start_pos is None:
        start_pos = np.arange(ma.shape[0])
    if start_pos2 is None:
        start_pos2 = start_pos

    xmesh, ymesh = np.meshgrid(start_pos, start_pos2)
    img3 = axHeat2.pcolormesh(xmesh.T, ymesh.T, ma, vmin=args.vMin, vmax=args.vMax, cmap=cmap, norm=norm)

    img3.set_rasterized(True)

    if args.region:

        # relabel xticks
        def relabel_ticks(ticks):
            if ticks[-1] - ticks[0] > 100000:
                labels = ["{:.2f} ".format((x / 1e6))
                          for x in ticks]
                labels[-1] += "Mbp"
            else:
                labels = ["{:,} ".format((x))
                          for x in ticks]
                labels[-1] += "bp"
            return labels
        xtick_lables = relabel_ticks(axHeat2.get_xticks())
        axHeat2.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_xticklabels(xtick_lables, size='small', rotation=45)

        ytick_lables = relabel_ticks(axHeat2.get_yticks())
        axHeat2.get_yaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_yticklabels(ytick_lables, size='small')

        """
        axHeat2.set_xticks([0, ma.shape[0]])
        axHeat2.set_xticklabels([args.region[1], args.region[2]], size=4, rotation=90)
        axHeat2.set_axis_off()
        """
    else:
        ticks = [int(pos[0] + (pos[1] - pos[0]) / 2) for pos in chrBinBoundaries.values()]
        labels = chrBinBoundaries.keys()
        axHeat2.set_xticks(ticks)
        if len(labels) > 20:
            axHeat2.set_xticklabels(labels, size=4, rotation=90)
        else:
            axHeat2.set_xticklabels(labels, size=8)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axHeat2)
    cax = divider.append_axes("right", size="2.5%", pad=0.09)
    if args.log1p:
        from matplotlib.ticker import LogFormatter
        formatter = LogFormatter(10, labelOnlyBase=False)
        # get a useful log scale
        # that looks like [1, 2, 5, 10, 20, 50, 100, ... etc]
        aa = np.array([1, 2, 5])
        tick_values = np.concatenate([aa * 10 ** x for x in range(10)])
        cbar = fig.colorbar(img3, ticks=tick_values, format=formatter, cax=cax)
    else:
        cbar = fig.colorbar(img3, cax=cax)

    cbar.solids.set_edgecolor("face")  # to avoid white lines in the color bar in pdf plots
    if args.scoreName:
        cbar.ax.set_ylabel(args.scoreName, rotation=270, size=8)
    axHeat2.set_ylim(start_pos2[0], start_pos2[-1])
    axHeat2.set_xlim(start_pos[0], start_pos[-1])

    if ylabel is not None:
        axHeat2.set_ylabel(ylabel)

    if xlabel is not None:
        axHeat2.set_xlabel(xlabel)


def translate_region(region_string):
    """
    Takes an string and returns a list
    of chrom, start, end.
    If the region string only contains
    the chrom, then start and end
    are set to a 0 and 1e15
    """

    region_string = region_string.translate(None, ",;!").replace("-", ":")
    fields = region_string.split(":")
    chrom = fields[0]
    try:
        region_start = int(fields[1])
    except IndexError:
        region_start = 0
    try:
        region_end = int(fields[2])
    except IndexError:
        region_end = 1e15  # vert large number

    return chrom, region_start, region_end


def plotPerChr(hic_matrix, cmap, args):
    """
    plots each chromosome individually, one after the other
    in one row. scale bar is added at the end
    """
    from math import ceil
    chromosomes = hic_matrix.getChrNames()
    chrom_per_row = 5
    num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
    num_cols = min(chrom_per_row, len(chromosomes))
    width_ratios = [1.0] * num_cols + [0.05]
    grids = gridspec.GridSpec(num_rows, num_cols + 1,
                              width_ratios=width_ratios,
                              height_ratios=[1] * num_rows)

    fig_height = 6 * num_rows
    fig_width = sum((np.array(width_ratios) + 0.05) * 6)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)

    norm = None
    if args.log1p:
        norm = LogNorm()
    chrom, start, end, _ = zip(*hic_matrix.cut_intervals)

    for idx, chrname in enumerate(chromosomes):
        row = idx / chrom_per_row
        col = idx % chrom_per_row

        axis = plt.subplot(grids[row, col])
        axis.set_title(chrname)
        chrom_range = hic_matrix.getChrBinRange(chrname)
        mat = hic_matrix.matrix[chrom_range[0]:chrom_range[1],
                                chrom_range[0]:chrom_range[1]].todense().astype(float)

        if args.log1p:
            mat += 1
        img = axis.imshow(mat, aspect='auto',
                          interpolation='spline16',
                          vmax=args.vMax, vmin=args.vMin, cmap=cmap,
                          norm=norm,
                          extent=[start[chrom_range[0]], end[chrom_range[1] - 1],
                                  end[chrom_range[1] - 1], start[chrom_range[0]]])

        img.set_rasterized(True)

        xticks = axis.get_xticks()
        if xticks[1] < 1e6:
            xlabels = ["{:.0f}".format(int(x) / 1e3)
                       for x in xticks]
            xlabels[-2] = "{} Kb".format(xlabels[-2])
        else:
            xlabels = ["{:.0f}".format(int(x) / 1e6)
                       for x in xticks]
            xlabels[-2] = "{} Mb".format(xlabels[-2])

        axis.set_xticklabels(xlabels, size='small')
        # yticks = axis.get_yticks()

        # ylabels = ["{:.0f}".format(int(x) / 1e6)
        #          for x in yticks]

        axis.get_xaxis().set_tick_params(
            which='both',
            bottom='on',
            direction='out')

        axis.set_yticklabels(xlabels, size='small')

    cbar3 = plt.subplot(grids[-1])
    cbar = fig.colorbar(img, cax=cbar3)
    cbar.solids.set_edgecolor("face")  # to avoid white lines in
    # the color bar in pdf plots
    cbar.ax.set_ylabel(args.scoreName, rotation=270, labelpad=20)


def main(args=None):
    args = parse_arguments().parse_args(args)

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

    sys.stderr.write("min: {}, max: {}\n".format(ma.matrix.data.min(), ma.matrix.data.max()))
    if args.region:
        chrom, region_start, region_end = translate_region(args.region)
        if chrom not in ma.interval_trees.keys():
            chrom = change_chrom_names(chrom)
            if chrom not in ma.interval_trees.keys():
                exit("Chromosome name {} in --region not in matrix".format(change_chrom_names(chrom)))

        args.region = [chrom, region_start, region_end]
        idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and
                                 x[1] >= region_start and x[2] < region_end])
        if args.region2:
            chrom2, region_start2, region_end2 = translate_region(args.region2)
            if chrom2 not in ma.interval_trees.keys():
                chrom2 = change_chrom_names(chrom2)
                if chrom2 not in ma.interval_trees.keys():
                    exit("Chromosome name {} in --region2 not in matrix".format(change_chrom_names(chrom2)))
            idx2, start_pos2 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom2 and
                                     x[1] >= region_start2 and x[2] < region_end2])
        else:
            idx2 = idx1
            chrom2 = chrom
            start_pos2 = start_pos1
        # select only relevant part
        matrix = np.asarray(ma.matrix[idx1, :][:, idx2].todense().astype(float))

    else:
        # TODO make start_pos1
        start_pos1 = None
        start_pos2 = None
        chrom = None
        chrom2 = None
        matrix = np.asanyarray(ma.getMatrix().astype(float))

    cmap = cm.get_cmap(args.colorMap)
    sys.stderr.write("Nan values set to black\n")
    cmap.set_bad('black')

    if args.perChromosome:
        plotPerChr(ma, cmap, args)

    else:
        fig_height = 6
        height = 4.8 / fig_height
        if args.whatToShow == 'both':
            fig_width = 11
            width = 4.9 / fig_width
            left_margin = (1.0 - 2 * width) * 0.35
        else:
            fig_width = 7
            width = 5.0 / fig_width
            left_margin = (1.0 - width) * 0.4

        fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)
        bottom = 0.8 / fig_height

        if args.log:
            mask = matrix == 0
            matrix[mask] = matrix[mask is False].min()
    #        matrix = -1 *  np.log(matrix)
            matrix = np.log(matrix)

        if args.whatToShow == '3D':
            position = [left_margin, bottom, width * 1.2, height * 1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

        if args.whatToShow == 'both':
            position = [left_margin, bottom, width * 1.2, height * 1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

            left_margin2 = 5.5 / fig_width + left_margin
            position = [left_margin2, bottom, width, height]
            plotHeatmap(matrix, ma.chrBinBoundaries, fig, position, args,
                        fig_width, cmap)

        if args.whatToShow == 'heatmap':
            position = [left_margin, bottom, width, height]
            if args.region:
                plotHeatmap_region(matrix, ma.chrBinBoundaries, fig, position,
                                   args, cmap, xlabel=chrom, ylabel=chrom2,
                                   start_pos=start_pos1, start_pos2=start_pos2)

            else:
                plotHeatmap(matrix, ma.chrBinBoundaries, fig, position,
                            args, fig_width, cmap)

    plt.savefig(args.outFileName, dpi=args.dpi)
