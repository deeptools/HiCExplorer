"""The figure of hicPlotMatrix.

The C++ binary loads the matrix, applies --clearMaskedBins, --chromosomeOrder,
--region and --region2, extracts the matrices pcolormesh draws (one per
chromosome with --perChromosome) as tiles, aggregated in blocks of bins where
an axis has more bins than the figure has pixels, replaces zeros, NaN and
infinite values for --log and --log1p and adds 1 for --log1p, and computes the bin start
positions, the chromosome extents and the bin size, and it reads the tracks:
the bigwig scores of every plotBigwig call with their positions and axis
limits, and the loop and TAD positions of plotLongRangeContacts and plotTADs.
This module runs the plotting part of hicexplorer/hicPlotMatrix.py:
plotHeatmap, plotBigwig, bigwig_axes_config, plotLongRangeContacts and
plotTADs are the reference's functions with their reading and arithmetic
replaced by that data, and draw() and plot_per_chromosome() are its main()
and plotPerChr() from the figure on.

Data:
    options      the command line values the plotting reads
    A matrix is given either as `matrix`, the path of one .npy file with the
    whole matrix (rows by columns, drawn with the reference's single pcolormesh
    call from `start_pos` and `start_pos2`), or as `tile_spec` with `matrix`
    null, when the matrix is reduced to blocks or does not fit one tile:
        tile_spec = {x_edges, y_edges (cell edges along the first and the
                     second matrix axis, one more than blocks),
                     bins, blocks, bins_per_block, blocks_by_chromosome,
                     aggregate ("none", "sum" or "mean"), dtype,
                     vmin, vmax (of all tiles; null for a single tile),
                     tiles: [{file, r0, r1, c0, c1}]}
        Each tile file holds the blocks r0..r1 of the first axis and c0..c1 of
        the second as an array of shape (c1 - c0, r1 - r0), pcolormesh's
        orientation. The tiles are drawn one after the other on the same axes
        with one colour scale, and each is freed after it is drawn.
    The arrays keep the type the C++ binary chose (int32 counts, int64 where
    int32 does not hold them, float32 otherwise).
    resolution, chrom_names, chromosome_start_end: [[chrom, start, end]]
    chromosome_order_as_bytes   the reference converts --chromosomeOrder to
                                bytes on the whole-matrix path (applied by the
                                C++ to the loop and TAD positions)
    perChromosome = false:
        matrix (.npy path), start_pos, start_pos2, xlabel, ylabel, region,
        bigwig_tracks, loop_positions, tad_positions
    perChromosome = true:
        chromosomes: [{name, matrix, start_pos, start_pos2, region, bigwig_tracks}]
    bigwig_tracks: null or {"horizontal": TRACKS, "vertical": TRACKS or null},
        TRACKS = [{"missing": false, "x", "y", "lim": [low, high]}, ...], ending
        with {"missing": true} where plotBigwig stops at a missing chromosome
    loop_positions, tad_positions: null or {"x", "y"}
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


def as_colour_input(array):
    """float32 tiles are widened to float64 for the colour mapping. matplotlib
    normalises float32 data in float32, which moves a few values across a
    colour boundary compared with the float64 path (up to an RMS of 4 on a
    figure of the test data); the copy is one tile (at most 4096 x 4096) and is
    freed with it. Integer tiles are widened by matplotlib itself and
    need no copy here."""
    if array.dtype == np.float32:
        return array.astype(np.float64)
    return array


def make_norm(args, spec=None):
    """The LogNorm of --log and --log1p; with several tiles the limits the
    C++ binary found over all of them fill the ones the user did not give."""
    if not (args.log1p or args.log):
        return None
    vmin = args.vMin
    vmax = args.vMax
    if spec is not None:
        if vmin is None:
            vmin = spec["vmin"]
        if vmax is None:
            vmax = spec["vmax"]
    return LogNorm(vmin=vmin, vmax=vmax)


def draw_tiles(axis, spec, vmin, vmax, cmap, norm, temporary):
    """One pcolormesh per tile at the block edges of the whole matrix, so the
    tiles meet exactly; the tile files the C++ binary wrote as temporary files
    are removed as they are read."""
    import os
    x_edges = np.asarray(spec["x_edges"], dtype=np.float64)
    y_edges = np.asarray(spec["y_edges"], dtype=np.float64)
    if len(spec["tiles"]) > 1 and norm is None:
        if vmin is None:
            vmin = spec["vmin"]
        if vmax is None:
            vmax = spec["vmax"]
    image = None
    for tile_info in spec["tiles"]:
        tile = np.load(tile_info["file"])
        if tile_info["file"] in temporary:
            try:
                os.remove(tile_info["file"])
            except OSError:
                pass
        image = axis.pcolormesh(
            x_edges[tile_info["r0"]:tile_info["r1"] + 1],
            y_edges[tile_info["c0"]:tile_info["c1"] + 1],
            as_colour_input(tile), vmin=vmin, vmax=vmax, cmap=cmap, norm=norm,
            shading='flat')
        image.set_rasterized(True)
        del tile
    return image


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                ylabel=None, start_pos=None, start_pos2=None, pNorm=None, pAxis=None, pBigwig=None,
                pLoops=None, pLoopLargeRegionsOperation=None, pHiCMatrix=None, pChromsomeStartEndDict=None, pResolution=None, pTads=None,
                pBigwigTracks=None, pLoopPositions=None, pTadPositions=None, pTileSpec=None,
                pTemporary=()):
    log.debug("plotting heatmap")
    shape = ma.shape if pTileSpec is None else tuple(pTileSpec["bins"])
    if shape[0] < 5:
        chr_names = " ".join([toString(x) for x in chrBinBoundaries.keys()])
        log.info("Matrix for {} too small to plot. Matrix size: {}".format(
            chr_names, shape))
        return
    if pAxis is not None:
        axHeat2 = pAxis
    else:
        axHeat2 = fig.add_axes(position)

    if args.title:
        axHeat2.set_title(toString(args.title))

    vmin = None
    vmax = None
    if pNorm is None:
        vmin = args.vMin
        vmax = args.vMax

    if pTileSpec is None:
        if start_pos2 is None:
            start_pos2 = start_pos

        xmesh, ymesh = np.meshgrid(start_pos, start_pos2)

        img3 = axHeat2.pcolormesh(
            xmesh.T, ymesh.T, as_colour_input(ma), vmin=vmin, vmax=vmax, cmap=cmap,
            norm=pNorm)
        img3.set_rasterized(True)
    else:
        img3 = draw_tiles(axHeat2, pTileSpec, vmin, vmax, cmap, pNorm, pTemporary)

    if args.region:
        xtick_lables = relabel_ticks(axHeat2.get_xticks())
        axHeat2.get_xaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_xticklabels(
            xtick_lables, size='small', rotation=args.rotationX)

        ytick_lables = relabel_ticks(axHeat2.get_yticks())
        axHeat2.get_yaxis().set_tick_params(which='both', bottom='on', direction='out')
        axHeat2.set_yticklabels(ytick_lables, size='small')
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
        plotLongRangeContacts(axHeat2, pLoopPositions, args.region)
    if pTads:
        plotTADs(axHeat2, pTadPositions, args.region)
    axHeat2.invert_yaxis()

    if pBigwig:

        axHeat2.xaxis.set_label_position("top")
        axHeat2.xaxis.tick_top()

        plotBigwig(pBigwig['axis'], pBigwigTracks['horizontal'], pVertical=False)

        if args.bigwigAdditionalVerticalAxis:
            plotBigwig(pBigwig['axis_vertical'], pBigwigTracks['vertical'], pVertical=True)


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


def plotBigwig(pAxis, pTracks, pVertical=False):
    """The drawing of the reference's plotBigwig, for the tracks the C++
    binary read (one per bigwig file, in the order of the axes)."""
    for i, track in enumerate(pTracks):
        pAxis[i].set_frame_on(False)

        if pVertical:
            pAxis[i].yaxis.set_visible(False)
        else:
            pAxis[i].xaxis.set_visible(False)
        if track["missing"]:
            # the reference logs that the file lacks the chromosome and returns
            return
        low, high = track["lim"]
        if pVertical:
            pAxis[i].set_ylim(low, high)
        else:
            pAxis[i].set_xlim(low, high)
        x_values = np.array(track["x"], dtype=np.float64)
        bigwig_scores = np.array(track["y"], dtype=np.float64)
        if pVertical:
            pAxis[i].fill_between(
                np.flip(bigwig_scores, 0), x_values, edgecolor='none')
        else:
            pAxis[i].fill_between(
                x_values, 0, bigwig_scores, edgecolor='none')


def plotLongRangeContacts(pAxis, pPositions, pRegion):
    x_list = pPositions["x"]
    y_list = pPositions["y"]

    if pRegion is not None and (int(pRegion[1]) != 0 and int(pRegion[2]) != 1e15):
        pAxis.set_xlim(int(pRegion[1]), int(pRegion[2]))
        pAxis.set_ylim(int(pRegion[1]), int(pRegion[2]))

    pAxis.plot(x_list, y_list, 's', lw=2,
               markerfacecolor='none', markeredgecolor='red')


def plotTADs(pAxis, pPositions, pRegion):
    x_list = pPositions["x"]
    y_list = pPositions["y"]

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
        spec = entry.get("tile_spec")
        matrix = None if spec is not None else np.load(entry["matrix"])
        norm = make_norm(args, spec)

        chr_bin_boundary = OrderedDict()
        chr_bin_boundary[chrname] = None

        args.region = entry["region"]
        start_pos1 = entry["start_pos"]
        start_pos2 = entry["start_pos2"]
        plotHeatmap(matrix, chr_bin_boundary, fig, None,
                    args, cmap, xlabel=chrname, ylabel=chrname,
                    start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=axis, pBigwig=bigwig_info,
                    pChromsomeStartEndDict=chromosome_start_end, pResolution=pResolution,
                    pBigwigTracks=entry["bigwig_tracks"], pTileSpec=spec,
                    pTemporary=set(data.get("temporary_files", [])))
    return fig


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    args = SimpleNamespace(**data["options"])
    if args.title:
        from unidecode import unidecode
        args.title = unidecode(args.title)

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
        spec = data.get("tile_spec")
        matrix = None if spec is not None else np.load(data["matrix"])
        args.region = data["region"]
        norm = make_norm(args, spec)

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
                    pLoops=args.loops, pLoopLargeRegionsOperation=args.loopLargeRegionsOperation, pHiCMatrix=None, pChromsomeStartEndDict=chromosome_start_end, pResolution=resolution, pTads=args.tads,
                    pBigwigTracks=data["bigwig_tracks"], pLoopPositions=data["loop_positions"],
                    pTadPositions=data["tad_positions"], pTileSpec=spec,
                    pTemporary=set(data.get("temporary_files", [])))

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
