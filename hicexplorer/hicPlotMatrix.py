from __future__ import division

import sys
import hicexplorer.HiCMatrix as HiCMatrix
from hicexplorer.utilities import writableFile
from hicexplorer.utilities import toString, toBytes
from hicexplorer._version import __version__
from hicexplorer.trackPlot import file_to_intervaltree
import numpy as np

from builtins import range
from past.builtins import zip
from future.utils import itervalues

import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from collections import OrderedDict


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

    parser.add_argument('--pca',
                        help='List of eigenvector from pca analysis as bigwig or bedgraph files.',
                        type=str,
                        default=None,
                        nargs='+')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser

# relabel xticks


def relabel_ticks(pXTicks):
    if pXTicks[-1] > 1.5e6:
        labels = ["{:.2f} ".format(x / 1e6)
                  for x in pXTicks]
        labels[-2] += " Mbp"
    elif pXTicks[-1] > 1500:
        labels = ["{:.0f}".format(x / 1e3)
                  for x in pXTicks]
        labels[-2] += " Kbp"
    else:
        labels = ["{:,} ".format((x))
                  for x in pXTicks]
        labels[-2] += " bp"
    return labels


def plotHeatmap3D(ma, fig, position, args, cmap):
    axHeat1 = fig.add_axes(position, projection='3d')
    axHeat1.view_init(50, -45)
    axHeat1.margins(0)

    X, Y = np.meshgrid(list(range(ma.shape[0])), list(range(ma.shape[0])))
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
    chrom = toString(chrom)
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                ylabel=None, start_pos=None, start_pos2=None, pNorm=None, pAxis=None, pPca=None):
    print("PLOT heatmap region")
    if pAxis:
        axHeat2 = pAxis
    else:
        axHeat2 = fig.add_axes(position)

    if args.title:
        axHeat2.set_title(args.title)

    if start_pos is None:
        start_pos = np.arange(ma.shape[0])
    if start_pos2 is None:
        start_pos2 = start_pos

    xmesh, ymesh = np.meshgrid(start_pos, start_pos2)
    # print(xmesh)
    # print(xmesh.T)

    img3 = axHeat2.pcolormesh(ymesh, xmesh, ma, vmin=args.vMin, vmax=args.vMax, cmap=cmap, norm=pNorm)
    # img3 = axHeat2.imshow(ma, cmap=cmap, vmin=args.vMin, vmax=args.vMax, interpolation='nearest', norm=pNorm)
    axHeat2.invert_yaxis()
    img3.set_rasterized(True)
    xticks = None
    if args.region:
        xtick_lables = relabel_ticks(axHeat2.get_xticks())
        axHeat2.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_xticklabels(xtick_lables, size='small', rotation=45)

        ytick_lables = relabel_ticks(axHeat2.get_yticks())
        axHeat2.get_yaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_yticklabels(ytick_lables, size='small')
        xticks = [xtick_lables]
        """
        axHeat2.set_xticks([0, ma.shape[0]])
        axHeat2.set_xticklabels([args.region[1], args.region[2]], size=4, rotation=90)
        axHeat2.set_axis_off()
        """
    else:

        ticks = [int(pos[0] + (pos[1] - pos[0]) / 2) for pos in itervalues(chrBinBoundaries)]
        labels = list(chrBinBoundaries)
        axHeat2.set_xticks(ticks)
        axHeat2.set_yticks(ticks)
        xticks = [labels, ticks]

        if len(labels) > 20:
            axHeat2.set_xticklabels(labels, size=4, rotation=90)
            axHeat2.set_yticklabels(labels, size=4)

        else:
            axHeat2.set_xticklabels(labels, size=8)
            axHeat2.set_yticklabels(labels, size=8)

    if pPca is None:
        divider = make_axes_locatable(axHeat2)
        cax = divider.append_axes("right", size="2.5%", pad=0.09)
    else:
        cax = pPca['axis_colorbar']
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

    if ylabel is not None:
        axHeat2.set_ylabel(ylabel)

    if xlabel is not None:
        axHeat2.set_xlabel(xlabel)

    if pPca:
        axHeat2.xaxis.set_label_position("top")
        axHeat2.xaxis.tick_top()
        if args.region:
            plotEigenvector(pPca['axis'], pPca['args'].pca, pRegion=pPca['args'].region, pXticks=xticks)
        else:
            plotEigenvector(pPca['axis'], pPca['args'].pca, pXticks=xticks, pChromosomeList=labels)


def translate_region(region_string):
    """
    Takes an string and returns a list
    of chrom, start, end.
    If the region string only contains
    the chrom, then start and end
    are set to a 0 and 1e15
    """

    if sys.version_info[0] == 2:
        region_string = region_string.translate(None, ",;!").replace("-", ":")
    if sys.version_info[0] == 3:
        region_string = region_string.replace(",", "")
        region_string = region_string.replace(";", "")
        region_string = region_string.replace("!", "")
        region_string = region_string.replace("-", ":")

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


def plotPerChr(hic_matrix, cmap, args, pPca):
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

    chrom, start, end, _ = zip(*hic_matrix.cut_intervals)

    for idx, chrname in enumerate(chromosomes):
        row = idx // chrom_per_row
        col = idx % chrom_per_row
        if pPca:
            inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, height_ratios=[0.85, 0.15], width_ratios=[0.93, 0.07],
                                                          subplot_spec=grids[row, col], wspace=0.0, hspace=0.1)
            axis = plt.subplot(inner_grid[0, 0])
            axis_eigenvector = plt.subplot(inner_grid[1, 0])
            axis_scale = plt.subplot(inner_grid[0, 1])

        else:
            axis = plt.subplot(grids[row, col])
            axis.set_title(chrname)
        chrom_range = hic_matrix.getChrBinRange(chrname)
        matrix = np.asarray(hic_matrix.matrix[chrom_range[0]:chrom_range[1],
                                              chrom_range[0]:chrom_range[1]].todense().astype(float))

        norm = None
        if args.log or args.log1p:
            mask = matrix == 0
            mask_nan = matrix == np.nan
            try:
                matrix[mask] = matrix[mask == False].min()
                matrix[mask_nan] = matrix[mask_nan == False].min()
                matrix = np.log(matrix)
            except:
                pass
        if args.log1p:
            matrix += 1
            norm = LogNorm()

        pca = None
        if pPca:
            pca = {'args': args, 'axis': None, 'axis_colorbar': None, 'nan_bins': hic_matrix.nan_bins}
            pca['axis'] = axis_eigenvector
            pca['axis_colorbar'] = axis_scale

        chr_bin_boundary = OrderedDict()
        chr_bin_boundary[chrname] = hic_matrix.chrBinBoundaries[chrname]
        plotHeatmap(matrix, chr_bin_boundary, fig, None,
                    args, cmap, xlabel=chrname, ylabel=chrname,
                    start_pos=None, start_pos2=None, pNorm=norm, pAxis=axis, pPca=pca)


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

        valid_chromosomes = []
        invalid_chromosomes = []
        args.chromosomeOrder = toBytes(args.chromosomeOrder)
        for chrom in args.chromosomeOrder:
            if chrom in ma.chrBinBoundaries:
                valid_chromosomes.append(chrom)
            else:
                invalid_chromosomes.append(chrom)

        ma.reorderChromosomes(valid_chromosomes)
        if len(invalid_chromosomes) > 0:
            sys.stderr.write("WARNING: The following chromosome/scaffold names were not found. Please check"
                             "the correct spelling of the chromosome names. \n")
            sys.stderr.write("\n".join(invalid_chromosomes))

    ma.restoreMaskedBins()
    if args.clearMaskedBins:
        ma.maskBins(ma.nan_bins)

    sys.stderr.write("min: {}, max: {}\n".format(ma.matrix.data.min(), ma.matrix.data.max()))
    start_pos1 = None
    start_pos2 = None
    if args.region:
        chrom, region_start, region_end = translate_region(args.region)
        if type(next(iter(ma.interval_trees))) is np.bytes_:
            chrom = toBytes(chrom)

        if chrom not in list(ma.interval_trees):
            chrom = change_chrom_names(chrom)
            if type(next(iter(ma.interval_trees))) is np.bytes_:
                chrom = toBytes(chrom)
            if chrom not in list(ma.interval_trees):
                exit("Chromosome name {} in --region not in matrix".format(change_chrom_names(chrom)))

        args.region = [chrom, region_start, region_end]
        idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and
                                 x[1] >= region_start and x[2] < region_end])

        if args.region2:
            chrom2, region_start2, region_end2 = translate_region(args.region2)

            if type(next(iter(ma.interval_trees))) is np.bytes_:
                chrom2 = toBytes(chrom)
            if chrom2 not in list(ma.interval_trees):
                chrom2 = change_chrom_names(chrom2)
                if type(next(iter(ma.interval_trees))) is np.bytes_:
                    chrom2 = toBytes(chrom)
                if chrom2 not in list(ma.interval_trees):
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
        chrom = None
        chrom2 = None
        matrix = np.asarray(ma.getMatrix().astype(float))

    cmap = cm.get_cmap(args.colorMap)
    sys.stderr.write("Nan values set to black\n")
    cmap.set_bad('black')

    pca = None
    if args.pca:
        pca = {'args': args, 'axis': None, 'axis_colorbar': None, 'nan_bins': ma.nan_bins}

    if args.perChromosome:
        plotPerChr(ma, cmap, args, pPca=pca)

    else:
        norm = None

        if args.log or args.log1p:
            mask = matrix == 0
            mask_nan = matrix == np.nan
            try:
                matrix[mask] = matrix[mask == False].min()
                matrix[mask_nan] = matrix[mask_nan == False].min()
                matrix = np.log(matrix)
            except:
                pass
        if args.log1p:
            matrix += 1
            norm = LogNorm()

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

        if args.pca:
            gs = gridspec.GridSpec(2, 2, height_ratios=[0.85, 0.15], width_ratios=[0.93, 0.07])
            gs.update(hspace=0.1)
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[1, 0])
            ax3 = plt.subplot(gs[0, 1])
            pca['axis'] = ax2
            pca['axis_colorbar'] = ax3

        else:
            ax1 = None
        bottom = 0.8 / fig_height

        if args.whatToShow == '3D':
            position = [left_margin, bottom, width * 1.2, height * 1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

        if args.whatToShow == 'both':
            position = [left_margin, bottom, width * 1.2, height * 1.2]
            plotHeatmap3D(matrix, fig, position, args, cmap)

            left_margin2 = 5.5 / fig_width + left_margin
            position = [left_margin2, bottom, width, height]
            plotHeatmap(matrix, ma.chrBinBoundaries, fig, position, args,
                        fig_width, cmap, pNorm=norm)

        if args.whatToShow == 'heatmap':
            position = [left_margin, bottom, width, height]
            plotHeatmap(matrix, ma.chrBinBoundaries, fig, position,
                        args, cmap, xlabel=chrom, ylabel=chrom2,
                        start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=ax1, pPca=pca)

    plt.savefig(args.outFileName, dpi=args.dpi)


def plotEigenvector(pAxis, pNameOfEigenvectorsList, pChromosomeList=None, pRegion=None, pXticks=None):

    pAxis.set_frame_on(False)

    file_format = pNameOfEigenvectorsList[0].split(".")[-1]
    if file_format != 'bedgraph' and file_format != 'bigwig':
        exit("Given eigenvector files are not bedgraph or bigwig")

    for eigenvector in pNameOfEigenvectorsList:
        if eigenvector.split('.')[-1] != file_format:
            exit("Eigenvector input files have different formats.")

    if pRegion:
        chrom, region_start, region_end = pRegion

    if file_format == "bigwig":
        raise NotImplementedError
    else:
        for i, eigenvectorFile in enumerate(pNameOfEigenvectorsList):
            interval_tree, min_value, max_value = file_to_intervaltree(eigenvectorFile)
            eigenvector = []
            if pChromosomeList:
                for chrom in pChromosomeList:
                    if chrom not in interval_tree:
                        print("Chromosome with no entry in the eigenvector found. Please exclude it from the matrix: {}. The eigenvector is left empty.".format(chrom))
                        return
                    for i, region in enumerate(sorted(interval_tree[chrom])):
                        if i == 0:
                            region_start = region[0]
                        region_end = region[1]
                        eigenvector.append(complex(region.data[0]).real)
                x = np.arange(0, len(eigenvector), 1)

            elif pRegion:
                if chrom not in interval_tree:
                    print("Chromosome with no entry in the eigenvector found. Please exclude it from the matrix: {}. The eigenvector is left empty.".format(chrom))
                    return
                for region in sorted(interval_tree[chrom][region_start:region_end]):
                    eigenvector.append(float(region.data[0]))
                step = (region_end*2 - region_start) // len(eigenvector)

                x = np.arange(region_start, region_end*2, int(step))
                while len(x) < len(eigenvector):
                    x = np.append(x[-1] + int(step))
                while len(eigenvector) < len(x):
                    x = x[:-1]

                pAxis.set_xlim(region_start, region_end*2)
            pAxis.fill_between(x, 0, eigenvector, edgecolor='none')

        if pXticks:
            if len(pXticks) == 1:
                pAxis.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')
                pAxis.set_xticklabels(pXticks[0], size='small', rotation=45)
            else:
                pAxis.set_xlim(0, len(eigenvector))
                pAxis.set_xticks(pXticks[1])
                pAxis.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')

                if len(pXticks[0]) > 20:
                    pAxis.set_xticklabels(pXticks[0], size=4, rotation=90)
                else:
                    pAxis.set_xticklabels(pXticks[0], size=8)
