"""hicPlotMatrix, implemented in C++.

main() runs the C++ executable hicPlotMatrix and parse_arguments() is built
from its --help-json (hicexplorer._cpp). The public functions keep the name,
the signature and the behaviour they had in HiCExplorer 3.7.6. relabel_ticks,
translate_region, getRegion, chromosome_start_end and make_start_pos_array call
hicexplorer._hicx.hicPlotMatrix (cpp/python/src/api/hicPlotMatrix.cpp).
plotHeatmap, plotPerChr, bigwig_axes_config, plotBigwig, plotLongRangeContacts
and plotTADs are matplotlib drawing code and stay Python
(cpp/scripts/public_api_allowlist.tsv); plotPerChr and plotBigwig take the data
they draw from hicexplorer._hicx.hicPlotMatrix.
"""
import logging
from collections import OrderedDict
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import matplotlib as mpl
from past.builtins import zip
import numpy as np
from hicexplorer.utilities import toString
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


matplotlib.use('Agg')


log = logging.getLogger(__name__)

from hicexplorer._hicx import hicPlotMatrix as _impl


def parse_arguments(args=None):
    """The argparse parser of hicPlotMatrix, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicPlotMatrix')


def relabel_ticks(pXTicks):
    return _impl.relabel_ticks(pXTicks)


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                ylabel=None, start_pos=None, start_pos2=None, pNorm=None, pAxis=None, pBigwig=None,
                pLoops=None, pLoopLargeRegionsOperation=None, pHiCMatrix=None, pChromsomeStartEndDict=None, pResolution=None, pTads=None):
    log.debug("plotting heatmap")
    if ma.shape[0] < 5:
        # This happens when a tiny matrix wants to be plotted, or by using per chromosome and
        # a small chromosome (eg. contig) is present.
        # Otherwise, pcolormesh will throw an error if the matrix size is 1.
        chr_names = " ".join([toString(x) for x in chrBinBoundaries.keys()])
        log.info("Matrix for {} too small to plot. Matrix size: {}".format(
            chr_names, ma.shape))
        return
    if pAxis is not None:
        axHeat2 = pAxis
    else:
        axHeat2 = fig.add_axes(position)

    if args.title:
        axHeat2.set_title(toString(args.title))

    if start_pos2 is None:
        start_pos2 = start_pos

    xmesh, ymesh = np.meshgrid(start_pos, start_pos2)

    # patch for newer matplotlib versions where vmin and vmax is for norm not allowed anymore
    vmin = None
    vmax = None
    if pNorm is None:
        vmin = args.vMin
        vmax = args.vMax

    img3 = axHeat2.pcolormesh(
        xmesh.T, ymesh.T, ma, vmin=vmin, vmax=vmax, cmap=cmap, norm=pNorm)
    img3.set_rasterized(True)

    if args.region:
        xtick_lables = relabel_ticks(axHeat2.get_xticks())
        axHeat2.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_xticklabels(
            xtick_lables, size='small', rotation=args.rotationX)

        ytick_lables = relabel_ticks(axHeat2.get_yticks())
        axHeat2.get_yaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_yticklabels(ytick_lables, size='small')
        xticks = [xtick_lables]
        log.debug("223")

        """
        axHeat2.set_xticks([0, ma.shape[0]])
        axHeat2.set_xticklabels([args.region[1], args.region[2]], size=4, rotation=90)
        axHeat2.set_axis_off()
        """
    else:

        pos = 0
        ticks = []

        for index, (k, v) in enumerate(pChromsomeStartEndDict.items()):
            ticks.append(pos)
            pos += v[1] - v[0]

        labels = list(chrBinBoundaries)
        axHeat2.set_xticks(ticks)
        axHeat2.set_yticks(ticks)
        labels = toString(labels)
        xticks = [labels, ticks]

        if len(labels) > 20:
            log.debug("243")

            axHeat2.set_xticklabels(
                labels, rotation=args.rotationX, fontsize=args.fontsize)
            # axHeat2.set_xticklabels(labels, size=4)

            axHeat2.set_yticklabels(
                labels, rotation=args.rotationY, fontsize=args.fontsize)

        else:
            log.debug("249")
            axHeat2.set_xticklabels(
                labels, rotation=args.rotationX, fontsize=args.fontsize)
            axHeat2.set_yticklabels(
                labels, rotation=args.rotationY, fontsize=args.fontsize)

    if pBigwig is None:
        divider = make_axes_locatable(axHeat2)
        cax = divider.append_axes("right", size="2.5%", pad=0.09)
    else:
        cax = pBigwig['axis_colorbar']

    cbar = fig.colorbar(img3, cax=cax)

    # to avoid white lines in the color bar in pdf plots
    cbar.solids.set_edgecolor("face")
    if args.scoreName:
        # cbar.ax.set_ylabel(args.scoreName, rotation=270, size=8)
        cbar.ax.set_ylabel(args.scoreName, size=8)

    if ylabel is not None:
        ylabel = toString(ylabel)
        axHeat2.set_ylabel(ylabel, fontsize=args.fontsize)

    if xlabel is not None:
        xlabel = toString(xlabel)
        axHeat2.set_xlabel(xlabel, fontsize=args.fontsize)
    if pLoops:
        log.debug('pLoops called')

        plotLongRangeContacts(axHeat2, pLoops, pLoopLargeRegionsOperation,
                              args.region, args.chromosomeOrder)
        # pLongRangeContacts=None, pHiCMatrix=None
        # plotLongRangeContacts(pAxis, pNameOfLongRangeContactsFile, pHiCMatrix)
    if pTads:
        plotTADs(axHeat2, pTads, pHiCMatrix,
                 args.region, args.chromosomeOrder)
    axHeat2.invert_yaxis()

    if pBigwig:

        axHeat2.xaxis.set_label_position("top")
        axHeat2.xaxis.tick_top()

        if args.region:
            plotBigwig(pBigwig['axis'], pBigwig['args'].bigwig, pChromosomeSizes=pChromsomeStartEndDict,
                       pRegion=pBigwig['args'].region, pXticks=xticks, pFlipBigwigSign=args.flipBigwigSign,
                       pScaleFactorBigwig=args.scaleFactorBigwig, pVertical=False,
                       pValueMin=args.vMinBigwig, pValueMax=args.vMaxBigwig, pResolution=pResolution)
        else:
            plotBigwig(pBigwig['axis'], pBigwig['args'].bigwig, pXticks=xticks, pChromosomeSizes=pChromsomeStartEndDict,
                       pFlipBigwigSign=args.flipBigwigSign, pScaleFactorBigwig=args.scaleFactorBigwig, pVertical=False,
                       pValueMin=args.vMinBigwig, pValueMax=args.vMaxBigwig, pResolution=pResolution)

        if args.bigwigAdditionalVerticalAxis:
            if args.region:
                plotBigwig(pBigwig['axis_vertical'], pBigwig['args'].bigwig[::-1], pChromosomeSizes=pChromsomeStartEndDict,
                           pRegion=pBigwig['args'].region, pXticks=xticks, pFlipBigwigSign=args.flipBigwigSign,
                           pScaleFactorBigwig=args.scaleFactorBigwig, pVertical=True,
                           pValueMin=args.vMinBigwig, pValueMax=args.vMaxBigwig, pResolution=pResolution)
            else:
                plotBigwig(pBigwig['axis_vertical'], pBigwig['args'].bigwig[::-1], pXticks=xticks, pChromosomeSizes=pChromsomeStartEndDict,
                           pFlipBigwigSign=args.flipBigwigSign, pScaleFactorBigwig=args.scaleFactorBigwig, pVertical=True,
                           pValueMin=args.vMinBigwig, pValueMax=args.vMaxBigwig, pResolution=pResolution)


def translate_region(region_string, ma):
    """
    Takes an string and returns a list
    of chrom, start, end.
    If the region string only contains
    the chrom, then start and end
    are set to matrix start and end on this chromosome
    """
    return _impl.translate_region(region_string, ma)


def plotPerChr(hic_matrix, cmap, args, pBigwig, pResolution):
    """
    plots each chromosome individually, one after the other
    in one row. scale bar is added at the end
    """
    from math import ceil
    chromosomes = _impl._plotPerChr_chromosomes(hic_matrix)
    chrom_per_row = 5
    num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
    num_cols = min(chrom_per_row, len(chromosomes))
    width_ratios = [1.0] * num_cols + [0.05]
    grids = gridspec.GridSpec(num_rows, num_cols + 1,
                              width_ratios=width_ratios,
                              height_ratios=[1] * num_rows)

    fig_height = 6 * num_rows
    fig_width = sum((ratio + 0.05) * 6 for ratio in width_ratios)
    if pBigwig:
        for i in range(len(args.bigwig)):
            fig_height += args.increaseFigureHeight
            # if args.bigwigAdditionalVerticalAxis:
            fig_width += args.increaseFigureWidth

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)

    _impl._plotPerChr_unpack(hic_matrix)
    for idx, chrname in enumerate(chromosomes):
        log.debug('chrom: {}'.format(chrname))
        bigwig_info = None
        # if pBigwig:
        # bigwig_info['axis'] = axis_eigenvector
        # bigwig_info['axis_colorbar'] = axis_scale
        row = idx // chrom_per_row
        col = idx % chrom_per_row
        if pBigwig:
            bigwig_info = {'args': args, 'axis': None,
                           'axis_colorbar': None, 'nan_bins': hic_matrix.nan_bins}
            # bigwig_info, axis = bigwig_axes_config(args, bigwig_info)
            # bigwig_info['nan_bins'] = hic_matrix.nan_bins
            # bigwig_info['args'] = args

            # inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, height_ratios=[0.85, 0.15], width_ratios=[0.93, 0.07],
            #                                               subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
            # axis = plt.subplot(inner_grid[0, 0])
            # axis_eigenvector = plt.subplot(inner_grid[1, 0])
            # axis_scale = plt.subplot(inner_grid[0, 1])
            number_of_rows_plot = len(args.bigwig)
            bigwig_heights = [0.07] * number_of_rows_plot
            bigwig_height_ratio = 0.95 - (0.07 * number_of_rows_plot)
            if bigwig_height_ratio < 0.4:
                bigwig_height_ratio = 0.4
                _ratio = 0.6 / len(number_of_rows_plot)
                bigwig_heights = [_ratio] * number_of_rows_plot

            if args.bigwigAdditionalVerticalAxis:
                # gs = gridspec.GridSpecFromSubplotSpec(1 + len(args.bigwig), 3, height_ratios=[0.90, 0.1], width_ratios=[0.15, 0.82, 0.03],
                #                                       subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
                # # gs = gridspec.GridSpec(1 + len(args.bigwig), 3, height_ratios=[0.90, 0.1], width_ratios=[0.15, 0.82, 0.03])
                # # gs.update(hspace=0.05, wspace=0.05)
                # bigwig_vertical_axis = plt.subplot(gs[0, 0])
                # axis = plt.subplot(gs[0, 1])
                # ax2 = plt.subplot(gs[1, 1])
                # ax3 = plt.subplot(gs[0, 2])

                # bigwig_info['axis'] = ax2
                # bigwig_info['axis_colorbar'] = ax3
                # bigwig_info['axis_vertical'] = bigwig_vertical_axis

                gs = gridspec.GridSpecFromSubplotSpec(1 + len(args.bigwig), 2 + len(args.bigwig), height_ratios=[0.95 - (0.07 * number_of_rows_plot), *bigwig_heights], width_ratios=[*bigwig_heights, 0.97 - (0.07 * number_of_rows_plot), 0.03],
                                                      subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
                # gs.update(hspace=0.05, wspace=0.05)
                # gs.update(hspace=0.05, wspace=0.05)
                axis = plt.subplot(gs[0, len(args.bigwig)])
                ax2_list = []
                for i in range(len(args.bigwig)):
                    ax2_list.append(plt.subplot(gs[1 + i, len(args.bigwig)]))

                bigwig_vertical_axis_list = []
                for i in range(len(args.bigwig)):
                    bigwig_vertical_axis_list.append(plt.subplot(gs[0, i]))
                # ax2 = plt.subplot(gs[1, 0])
                ax3 = plt.subplot(gs[0, len(args.bigwig) + 1])
                bigwig_info['axis'] = ax2_list
                bigwig_info['axis_colorbar'] = ax3
                bigwig_info['axis_vertical'] = bigwig_vertical_axis_list

            else:
                # [0.95 - (0.07 * number_of_rows_plot), *z_score_heights], width_ratios=[0.75, 0.25])
                gs = gridspec.GridSpecFromSubplotSpec(1 + len(args.bigwig), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *bigwig_heights], width_ratios=[0.97, 0.03],
                                                      subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
                # gs.update(hspace=0.05, wspace=0.05)
                axis = plt.subplot(gs[0, 0])
                ax2_list = []
                for i in range(len(args.bigwig)):
                    ax2_list.append(plt.subplot(gs[1 + i, 0]))
                # ax2 = plt.subplot(gs[1, 0])
                ax3 = plt.subplot(gs[0, 1])
                bigwig_info['axis'] = ax2_list
                bigwig_info['axis_colorbar'] = ax3
        else:
            axis = plt.subplot(grids[row, col])
            axis.set_title(toString(chrname))
        matrix, chr_bin_boundary = _impl._plotPerChr_matrix(hic_matrix, chrname, args)

        norm = None
        if args.log1p:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)

        elif args.log:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)

        args.region = toString(chrname)
        chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(
            args, hic_matrix)
        plotHeatmap(matrix, chr_bin_boundary, fig, None,
                    args, cmap, xlabel=chrname, ylabel=chrname,
                    start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=axis, pBigwig=bigwig_info,
                    pChromsomeStartEndDict=chromosome_start_end(hic_matrix), pResolution=pResolution)
    return fig


def getRegion(args, ma):
    return _impl.getRegion(args, ma)


def bigwig_axes_config(pArgs, pBigWigInfo):

    # increase figure height to accommodate bigwig track

    # for i in range(len(pArgs.bigwig)):
    #     pFigureHeight += pArgs.increaseFigureHeight
    #     # if args.bigwigAdditionalVerticalAxis:
    #     pFigureWidth += pArgs.increaseFigureWidth

    number_of_rows_plot = len(pArgs.bigwig) * 2

    if pArgs.bigwigAdditionalVerticalAxis:

        number_of_columns = number_of_rows_plot + 11

        # vertical bigwig plot axes
        bigwig_vertical_axis = []
        for i in range(0, number_of_rows_plot, 2):
            bigwig_vertical_axis.append(plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, 0 + i), colspan=2, rowspan=8))

        # heatmap axis
        ax1 = plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, number_of_rows_plot), colspan=10, rowspan=8)
        # colorbar for the heatmap axis
        ax3 = plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, number_of_rows_plot + 10), colspan=1, rowspan=8)

        # horizontal bigwig axes
        ax2_list = []
        for i in range(0, number_of_rows_plot, 2):
            ax2_list.append(plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (10 + i, number_of_rows_plot), colspan=10, rowspan=2))

        pBigWigInfo['axis'] = ax2_list
        pBigWigInfo['axis_colorbar'] = ax3
        pBigWigInfo['axis_vertical'] = bigwig_vertical_axis
    else:

        # number_of_rows_plot = number_of_rows_plot *2
        ax2_list = []

        # heatmap axis
        ax1 = plt.subplot2grid((10 + number_of_rows_plot, 11), (0, 0), colspan=10, rowspan=10)
        # colorbar for the heatmap axis
        ax3 = plt.subplot2grid((10 + number_of_rows_plot, 11), (0, 10), colspan=1, rowspan=10)

        for i in range(0, number_of_rows_plot, 2):
            ax2_list.append(plt.subplot2grid((10 + number_of_rows_plot, 11), (10 + i, 0), colspan=10, rowspan=2))

        pBigWigInfo['axis'] = ax2_list
        pBigWigInfo['axis_colorbar'] = ax3

    return pBigWigInfo, ax1


def chromosome_start_end(pHiCMatrix):
    return _impl.chromosome_start_end(pHiCMatrix)


def make_start_pos_array(ma):
    return _impl.make_start_pos_array(ma)


def plotBigwig(pAxis, pNameOfBigwigList, pChromosomeSizes=None, pRegion=None, pXticks=None,
               pFlipBigwigSign=None, pScaleFactorBigwig=None, pVertical=False,
               pValueMin=None, pValueMax=None, pResolution=None):
    log.debug('plotting bigwig file')

    # the file names, the scores read from the bigwig files, their positions,
    # the axis limits and the flipped, scaled and clipped values are computed
    # in C++ (hicexplorer._hicx.hicPlotMatrix); this function draws them.
    _impl._bigwig_check(pNameOfBigwigList)
    for i, bigwigFile in enumerate(pNameOfBigwigList):
        pAxis[i].set_frame_on(False)

        if pVertical:

            pAxis[i].yaxis.set_visible(False)

        else:
            pAxis[i].xaxis.set_visible(False)
        track = _impl._bigwig_track(bigwigFile, pChromosomeSizes, pRegion, pVertical, pFlipBigwigSign,
                                    pScaleFactorBigwig, pValueMin, pValueMax, pResolution)
        if track is None:
            return
        x_values, bigwig_scores, limits = track
        if limits is not None:
            if pVertical:
                pAxis[i].set_ylim(limits[0], limits[1])
            else:
                pAxis[i].set_xlim(limits[0], limits[1])
        if pVertical:
            pAxis[i].fill_between(bigwig_scores, x_values, edgecolor='none')
        else:
            pAxis[i].fill_between(x_values, 0, bigwig_scores, edgecolor='none')


def plotLongRangeContacts(pAxis, pNameOfLongRangeContactsFile, pLoopLargeRegionsOperation, pRegion, pChromosomeOrder):

    x_list = []
    y_list = []
    log.debug('pRegion {}'.format(pRegion))
    with open(pNameOfLongRangeContactsFile, 'rb') as file:
        for line in file.readlines():
            line = toString(line)
            fields = line.strip().split('\t')
            try:
                chrom_X, start_X, end_X = fields[0:3]
                chrom_Y, start_Y, end_Y = fields[3:6]

                if pRegion is not None and (chrom_X != pRegion[0] or chrom_Y != pRegion[0]):
                    continue
                elif pChromosomeOrder is not None and (chrom_X not in pChromosomeOrder or chrom_Y not in pChromosomeOrder):
                    continue

                if pLoopLargeRegionsOperation == 'first':
                    x = int(start_X)
                    y = int(start_Y)
                elif pLoopLargeRegionsOperation == 'last':
                    x = int(end_X)
                    y = int(end_Y)
                elif pLoopLargeRegionsOperation == 'center':
                    x = (int(start_X) + int(end_X)) // 2
                    y = (int(start_Y) + int(end_Y)) // 2

                log.debug('x {} y {}'.format(x, y))
                if x >= int(pRegion[1]) and x <= int(pRegion[2]):
                    if y >= int(pRegion[1]) and y <= int(pRegion[2]):
                        x_list.append(x)
                        y_list.append(y)
            except Exception:
                pass

        if pRegion is not None and (int(pRegion[1]) != 0 and int(pRegion[2]) != 1e15):
            pAxis.set_xlim(int(pRegion[1]), int(pRegion[2]))
            pAxis.set_ylim(int(pRegion[1]), int(pRegion[2]))

        pAxis.plot(x_list, y_list, 's', lw=2,
                   markerfacecolor='none', markeredgecolor='red')


def plotTADs(pAxis, pNameOfLongRangeContactsFile, pHiCMatrix, pRegion, pChromosomeOrder):
    x_list = []
    y_list = []
    with open(pNameOfLongRangeContactsFile, 'rb') as file:
        for line in file.readlines():
            line = toString(line)
            fields = line.strip().split('\t')
            try:
                chrom_X, start_X, end_X = fields[0:3]

                if pRegion is not None and chrom_X != pRegion[0]:
                    continue
                elif pChromosomeOrder is not None and chrom_X not in pChromosomeOrder:
                    continue

                x = int(start_X)
                y = int(end_X)
                if x >= int(pRegion[1]) and x <= int(pRegion[2]):
                    if y >= int(pRegion[1]) and y <= int(pRegion[2]):
                        x_list.append(x)
                        y_list.append(y)
            except Exception as exp:
                log.debug('Exception! {}'.format(str(exp)))

        if pRegion is not None and (int(pRegion[1]) != 0 and int(pRegion[2]) != 1e15):
            pAxis.set_xlim(int(pRegion[1]), int(pRegion[2]))
            pAxis.set_ylim(int(pRegion[1]), int(pRegion[2]))

        for x_id, y_id in zip(x_list, y_list):
            pAxis.plot([x_id, x_id], [y_id, x_id], 'k')
            pAxis.plot([x_id, y_id], [y_id, y_id], 'k')


def main(args=None):
    """Run the C++ implementation of hicPlotMatrix with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicPlotMatrix')(args)
