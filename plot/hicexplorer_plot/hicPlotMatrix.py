"""The figure of hicPlotMatrix.

The C++ binary loads the matrix, applies --clearMaskedBins, --chromosomeOrder,
--region and --region2, extracts the dense matrices pcolormesh draws (one per
chromosome with --perChromosome), replaces zeros, NaN and infinite values
for --log and --log1p and adds 1 for --log1p, and computes the bin start
positions, the chromosome extents and the bin size. This module runs the
plotting part of hicexplorer/hicPlotMatrix.py: plotHeatmap, plotBigwig,
bigwig_axes_config, plotLongRangeContacts and plotTADs are the reference's
functions, and draw() and plot_per_chromosome() are its main() and
plotPerChr() from the figure on, fed with the computed data.

Data:
    options      the command line values the plotting reads
    resolution, chrom_names, chromosome_start_end: [[chrom, start, end]]
    chromosome_order_as_bytes   the reference converts --chromosomeOrder to
                                bytes on the whole-matrix path
    perChromosome = false:
        matrix (.npy path), start_pos, start_pos2, xlabel, ylabel, region
    perChromosome = true:
        chromosomes: [{name, matrix, start_pos, start_pos2, region}]
    temporary_files
"""

import logging
from collections import OrderedDict
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
from matplotlib import colormaps as cm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import matplotlib as mpl  # noqa: E402
import pyBigWig  # noqa: E402

log = logging.getLogger(__name__)


def toString(s):
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        return s.decode('ascii')
    if isinstance(s, list):
        return [toString(x) for x in s]
    if isinstance(s, np.ndarray):
        return s.astype(str)
    return s


def toBytes(s):
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [toBytes(x) for x in s]
    return s


def check_chrom_str_bytes(pIteratableObj, pObj):
    if isinstance(pObj, list) and len(pObj) > 0:
        type_ = type(pObj[0])
    else:
        type_ = type(pObj)
    if not isinstance(type(next(iter(pIteratableObj))), type_):
        if type(next(iter(pIteratableObj))) is str:
            pObj = toString(pObj)
        elif type(next(iter(pIteratableObj))) in [bytes, np.bytes_]:
            pObj = toBytes(pObj)
    return pObj


def change_chrom_names(chrom):
    chrom = toString(chrom)
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    else:
        chrom = 'chr' + chrom
    return chrom


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
        labels = ["{:.2f} ".format((x))
                  for x in pXTicks]
        labels[-2] += " bp"
    return labels


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                ylabel=None, start_pos=None, start_pos2=None, pNorm=None, pAxis=None, pBigwig=None,
                pLoops=None, pLoopLargeRegionsOperation=None, pHiCMatrix=None, pChromsomeStartEndDict=None, pResolution=None, pTads=None):
    log.debug("plotting heatmap")
    if ma.shape[0] < 5:
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
            axHeat2.set_xticklabels(
                labels, rotation=args.rotationX, fontsize=args.fontsize)
            axHeat2.set_yticklabels(
                labels, rotation=args.rotationY, fontsize=args.fontsize)

        else:
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

    cbar.solids.set_edgecolor("face")
    if args.scoreName:
        cbar.ax.set_ylabel(args.scoreName, size=8)

    if ylabel is not None:
        ylabel = toString(ylabel)
        axHeat2.set_ylabel(ylabel, fontsize=args.fontsize)

    if xlabel is not None:
        xlabel = toString(xlabel)
        axHeat2.set_xlabel(xlabel, fontsize=args.fontsize)
    if pLoops:
        plotLongRangeContacts(axHeat2, pLoops, pLoopLargeRegionsOperation,
                              args.region, args.chromosomeOrder)
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


def bigwig_axes_config(pArgs, pBigWigInfo):
    number_of_rows_plot = len(pArgs.bigwig) * 2

    if pArgs.bigwigAdditionalVerticalAxis:

        number_of_columns = number_of_rows_plot + 11

        bigwig_vertical_axis = []
        for i in range(0, number_of_rows_plot, 2):
            bigwig_vertical_axis.append(plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, 0 + i), colspan=2, rowspan=8))

        ax1 = plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, number_of_rows_plot), colspan=10, rowspan=8)
        ax3 = plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (0, number_of_rows_plot + 10), colspan=1, rowspan=8)

        ax2_list = []
        for i in range(0, number_of_rows_plot, 2):
            ax2_list.append(plt.subplot2grid((10 + number_of_rows_plot, number_of_columns), (10 + i, number_of_rows_plot), colspan=10, rowspan=2))

        pBigWigInfo['axis'] = ax2_list
        pBigWigInfo['axis_colorbar'] = ax3
        pBigWigInfo['axis_vertical'] = bigwig_vertical_axis
    else:
        ax2_list = []

        ax1 = plt.subplot2grid((10 + number_of_rows_plot, 11), (0, 0), colspan=10, rowspan=10)
        ax3 = plt.subplot2grid((10 + number_of_rows_plot, 11), (0, 10), colspan=1, rowspan=10)

        for i in range(0, number_of_rows_plot, 2):
            ax2_list.append(plt.subplot2grid((10 + number_of_rows_plot, 11), (10 + i, 0), colspan=10, rowspan=2))

        pBigWigInfo['axis'] = ax2_list
        pBigWigInfo['axis_colorbar'] = ax3

    return pBigWigInfo, ax1


def plotBigwig(pAxis, pNameOfBigwigList, pChromosomeSizes=None, pRegion=None, pXticks=None,
               pFlipBigwigSign=None, pScaleFactorBigwig=None, pVertical=False,
               pValueMin=None, pValueMax=None, pResolution=None):
    for file in pNameOfBigwigList:
        file_format = file.split(".")[-1]
        if file_format != 'bigwig' and file_format != 'bw':
            log.error("Given files are not bigwig")
            exit(1)

    if file_format == "bigwig" or file_format == 'bw':
        for i, bigwigFile in enumerate(pNameOfBigwigList):
            x_values = []
            bigwig_scores = []
            pAxis[i].set_frame_on(False)

            if pVertical:
                pAxis[i].yaxis.set_visible(False)
            else:
                pAxis[i].xaxis.set_visible(False)
            bw = pyBigWig.open(bigwigFile)
            bigwig_scores = []
            if pRegion:
                chrom, region_start, region_end = pRegion
                region_end = min(region_end, pChromosomeSizes[chrom][1])
                chrom = check_chrom_str_bytes(bw.chroms(), chrom)
                if chrom not in list(bw.chroms().keys()):
                    chrom = change_chrom_names(chrom)
                    if chrom not in list(bw.chroms().keys()):
                        log.info(
                            "bigwig file has no chromosome named: {}.".format(chrom))
                        return

                bigwig_end = min(bw.chroms()[chrom], region_end)

                num_bins = int(bigwig_end - region_start) // pResolution
                scores_per_bin = np.array(
                    bw.stats(chrom, region_start, bigwig_end, nBins=num_bins)).astype(float)
                if scores_per_bin is None:
                    log.info(
                        "Chromosome {} has no entries in bigwig file.".format(chrom))
                    return

                _x_vals = np.linspace(region_start, region_end, num_bins)
                assert len(_x_vals) == len(scores_per_bin)
                x_values.extend(_x_vals)
                bigwig_scores.extend(scores_per_bin)
                if pVertical:
                    pAxis[i].set_ylim(region_start, region_end)
                else:
                    pAxis[i].set_xlim(region_start, region_end)

            elif pChromosomeSizes:
                chrom_length_sum = 0
                min_start = None
                for chrom in pChromosomeSizes:
                    chrom_ = check_chrom_str_bytes(bw.chroms(), chrom)

                    if chrom_ not in list(bw.chroms().keys()):
                        chrom_ = 'chr' + chrom_
                        if chrom_ not in list(bw.chroms().keys()):
                            log.info(
                                "bigwig file as no chromosome named: {}.".format(chrom))
                            return
                    bigwig_end = min(bw.chroms()[chrom_], pChromosomeSizes[chrom][1])
                    num_bins = int(bigwig_end - pChromosomeSizes[chrom][0]) // pResolution

                    scores_per_bin = np.array(
                        bw.stats(chrom_, pChromosomeSizes[chrom][0], bigwig_end, nBins=num_bins)).astype(float)

                    if scores_per_bin is None:
                        log.info(
                            "Chromosome {} has no entries in bigwig file.".format(chrom))
                        return

                    _x_vals = np.linspace(
                        chrom_length_sum + pChromosomeSizes[chrom][0], chrom_length_sum + pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0], num_bins)

                    assert len(_x_vals) == len(scores_per_bin)
                    x_values.extend(_x_vals)

                    bigwig_scores.extend(scores_per_bin)

                    chrom_length_sum += pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0]
                    if min_start is None:
                        min_start = pChromosomeSizes[chrom][0]
                    elif min_start > pChromosomeSizes[chrom][0]:
                        min_start = pChromosomeSizes[chrom][0]
                if pVertical:
                    pAxis[i].set_ylim(min_start, chrom_length_sum)
                else:
                    pAxis[i].set_xlim(min_start, chrom_length_sum)

            bigwig_scores = np.array(bigwig_scores)
            if pFlipBigwigSign:
                bigwig_scores *= -1
            if pScaleFactorBigwig is not None and pScaleFactorBigwig != 1.0:
                bigwig_scores *= pScaleFactorBigwig
            if pValueMin is not None or pValueMax is not None:
                bigwig_scores = bigwig_scores.clip(pValueMin, pValueMax)

            if x_values is not None and bigwig_scores is not None:
                if pVertical:
                    pAxis[i].fill_between(
                        np.flip(bigwig_scores, 0), x_values, edgecolor='none')
                else:
                    pAxis[i].fill_between(
                        x_values, 0, bigwig_scores, edgecolor='none')


def plotLongRangeContacts(pAxis, pNameOfLongRangeContactsFile, pLoopLargeRegionsOperation, pRegion, pChromosomeOrder):

    x_list = []
    y_list = []
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


def plot_per_chromosome(data, cmap, args, pBigwig, pResolution, chromosome_start_end):
    """plotPerChr from the figure on, with the per chromosome matrices,
    start positions and regions computed by the C++ binary."""
    from math import ceil
    chromosomes = [entry["name"] for entry in data["chromosomes"]]
    chrom_per_row = 5
    num_rows = int(ceil(float(len(chromosomes)) / chrom_per_row))
    num_cols = min(chrom_per_row, len(chromosomes))
    width_ratios = [1.0] * num_cols + [0.05]
    grids = gridspec.GridSpec(num_rows, num_cols + 1,
                              width_ratios=width_ratios,
                              height_ratios=[1] * num_rows)

    fig_height = 6 * num_rows
    fig_width = sum((np.array(width_ratios) + 0.05) * 6)
    if pBigwig:
        for i in range(len(args.bigwig)):
            fig_height += args.increaseFigureHeight
            fig_width += args.increaseFigureWidth

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)

    for idx, entry in enumerate(data["chromosomes"]):
        chrname = entry["name"]
        bigwig_info = None
        row = idx // chrom_per_row
        col = idx % chrom_per_row
        if pBigwig:
            bigwig_info = {'args': args, 'axis': None,
                           'axis_colorbar': None, 'nan_bins': None}
            number_of_rows_plot = len(args.bigwig)
            bigwig_heights = [0.07] * number_of_rows_plot
            bigwig_height_ratio = 0.95 - (0.07 * number_of_rows_plot)
            if bigwig_height_ratio < 0.4:
                bigwig_height_ratio = 0.4
                _ratio = 0.6 / len(number_of_rows_plot)
                bigwig_heights = [_ratio] * number_of_rows_plot

            if args.bigwigAdditionalVerticalAxis:
                gs = gridspec.GridSpecFromSubplotSpec(1 + len(args.bigwig), 2 + len(args.bigwig), height_ratios=[0.95 - (0.07 * number_of_rows_plot), *bigwig_heights], width_ratios=[*bigwig_heights, 0.97 - (0.07 * number_of_rows_plot), 0.03],
                                                      subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
                axis = plt.subplot(gs[0, len(args.bigwig)])
                ax2_list = []
                for i in range(len(args.bigwig)):
                    ax2_list.append(plt.subplot(gs[1 + i, len(args.bigwig)]))

                bigwig_vertical_axis_list = []
                for i in range(len(args.bigwig)):
                    bigwig_vertical_axis_list.append(plt.subplot(gs[0, i]))
                ax3 = plt.subplot(gs[0, len(args.bigwig) + 1])
                bigwig_info['axis'] = ax2_list
                bigwig_info['axis_colorbar'] = ax3
                bigwig_info['axis_vertical'] = bigwig_vertical_axis_list

            else:
                gs = gridspec.GridSpecFromSubplotSpec(1 + len(args.bigwig), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *bigwig_heights], width_ratios=[0.97, 0.03],
                                                      subplot_spec=grids[row, col], wspace=0.1, hspace=0.1)
                axis = plt.subplot(gs[0, 0])
                ax2_list = []
                for i in range(len(args.bigwig)):
                    ax2_list.append(plt.subplot(gs[1 + i, 0]))
                ax3 = plt.subplot(gs[0, 1])
                bigwig_info['axis'] = ax2_list
                bigwig_info['axis_colorbar'] = ax3
        else:
            axis = plt.subplot(grids[row, col])
            axis.set_title(toString(chrname))
        matrix = np.load(entry["matrix"])

        norm = None
        if args.log1p:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)
        elif args.log:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)

        chr_bin_boundary = OrderedDict()
        chr_bin_boundary[chrname] = None

        args.region = entry["region"]
        start_pos1 = entry["start_pos"]
        start_pos2 = entry["start_pos2"]
        plotHeatmap(matrix, chr_bin_boundary, fig, None,
                    args, cmap, xlabel=chrname, ylabel=chrname,
                    start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=axis, pBigwig=bigwig_info,
                    pChromsomeStartEndDict=chromosome_start_end, pResolution=pResolution)
    return fig


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    args = SimpleNamespace(**data["options"])
    if args.title:
        from unidecode import unidecode
        args.title = unidecode(args.title)
    if data["chromosome_order_as_bytes"] and args.chromosomeOrder is not None:
        args.chromosomeOrder = toBytes(args.chromosomeOrder)

    chromosome_start_end = {}
    for chrom, start, end in data["chromosome_start_end"]:
        chromosome_start_end[chrom] = (start, end)
    chrom_sizes = OrderedDict((name, None) for name in data["chrom_names"])
    resolution = data["resolution"]

    cmap = cm.get_cmap(args.colorMap)
    cmap.set_bad('black')

    bigwig_info = None
    if args.bigwig:
        bigwig_info = {'args': args, 'axis': None,
                       'axis_colorbar': None, 'nan_bins': None}

    if args.perChromosome:
        fig = plot_per_chromosome(data, cmap, args, pBigwig=bigwig_info, pResolution=resolution,
                                  chromosome_start_end=chromosome_start_end)

    else:
        matrix = np.load(data["matrix"])
        args.region = data["region"]
        norm = None
        if args.log1p:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)
        elif args.log:
            norm = LogNorm(vmin=args.vMin, vmax=args.vMax)

        fig_height = 7
        fig_width = 8

        if args.bigwig:
            for i in range(len(args.bigwig)):
                fig_height += args.increaseFigureHeight
                fig_width += args.increaseFigureWidth
        height = 4.8 / fig_height

        width = 5.0 / fig_width
        left_margin = (1.0 - width) * 0.5

        fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)

        if args.bigwig:
            bigwig_info, ax1 = bigwig_axes_config(args, bigwig_info)

        else:
            ax1 = None
        bottom = 1.3 / fig_height

        position = [left_margin, bottom, width, height]
        plotHeatmap(matrix, chrom_sizes, fig, position,
                    args, cmap, xlabel=data["xlabel"], ylabel=data["ylabel"],
                    start_pos=data["start_pos"], start_pos2=data["start_pos2"], pNorm=norm, pAxis=ax1, pBigwig=bigwig_info,
                    pLoops=args.loops, pLoopLargeRegionsOperation=args.loopLargeRegionsOperation, pHiCMatrix=None, pChromsomeStartEndDict=chromosome_start_end, pResolution=resolution, pTads=args.tads)

    if not args.disable_tight_layout:
        if args.perChromosome or args.bigwig:
            try:
                plt.tight_layout()
            except UserWarning:
                log.info("Failed to tight layout. Using regular plot.")
            except ValueError:
                log.info("Failed to tight layout. Using regular plot.")
    plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close(fig)
