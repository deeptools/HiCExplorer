import logging
from collections import OrderedDict
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
import argparse
from past.builtins import zip
import pyBigWig
import numpy as np
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
from hicexplorer.utilities import check_chrom_str_bytes
from hicexplorer.utilities import remove_non_ascii
from hicexplorer.utilities import change_chrom_names
from hicexplorer.utilities import enlarge_bins
from hicexplorer.utilities import toString, toBytes
from hicexplorer.utilities import writableFile
from hicmatrix import HiCMatrix
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


matplotlib.use('Agg')


log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(add_help=False,
                                     description='Creates a heatmap of a Hi-C matrix.')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--matrix', '-m',
                                help='Path of the Hi-C matrix to plot.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-out',
                                help='File name to save the image.',
                                type=writableFile,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--title', '-t',
                           help='Plot title.')

    parserOpt.add_argument('--scoreName', '-s',
                           help='Score name label for the heatmap legend.')

    parserOpt.add_argument('--perChromosome',
                           help='Instead of plotting the whole matrix, '
                           'each chromosome is plotted next to the other. '
                           'This parameter is not compatible with --region.',
                           action='store_true')

    parserOpt.add_argument('--clearMaskedBins',
                           help='If set, masked bins are removed from the matrix '
                           'and the nearest bins are extended to cover the empty space '
                           'instead of plotting black lines.',
                           action='store_true')

    parserOpt.add_argument('--chromosomeOrder',
                           help='Chromosomes and order in which the '
                           'chromosomes should be plotted. This option '
                           'overrides --region and --region2.',
                           nargs='+')

    parserOpt.add_argument('--region',
                           help='Plot only this region. The format is '
                           'chr:start-end The plotted region contains '
                           'the main diagonal and is symmetric unless '
                           '--region2 is given.'
                           )

    parserOpt.add_argument('--region2',
                           help='If given, then only the region defined by '
                           '--region and --region2 is given. The format '
                           'is the same as --region1.'
                           )

    parserOpt.add_argument('--log1p',
                           help='Plot the log1p of the matrix values.',
                           action='store_true')

    parserOpt.add_argument('--log',
                           help='Plot the *MINUS* log of the matrix values.',
                           action='store_true')

    parserOpt.add_argument('--colorMap',
                           help='Color map to use for the heatmap. Available '
                           'values can be seen here: '
                           'http://matplotlib.org/examples/color/colormaps_reference.html',
                           default='RdYlBu_r')

    parserOpt.add_argument('--vMin',
                           help='Minimum score value.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--vMax',
                           help='Maximum score value.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--dpi',
                           help='Resolution for the image in case the'
                           'output is a raster graphics image (e.g png, jpg).',
                           type=int,
                           default=72)

    parserOpt.add_argument('--bigwig',
                           help='Bigwig file to plot below the matrix. This can for '
                           'example be used to visualize A/B compartments or '
                           'ChIP-seq data.',
                           type=str,
                           default=None,
                           nargs='+')
    parserOpt.add_argument('--bigwigAdditionalVerticalAxis',
                           help='Add an additional axis to determine the values of a bigwig file in 2D better.',
                           action='store_true')
    parserOpt.add_argument('--vMinBigwig',
                           help='Minimum score value for bigwig.',
                           type=float,
                           default=None)

    parserOpt.add_argument('--vMaxBigwig',
                           help='Maximum score value for bigwig',
                           type=float,
                           default=None)
    parserOpt.add_argument('--flipBigwigSign',
                           help='The sign of the bigwig values are flipped. Useful if hicPCA gives inverted values.',
                           action='store_true')
    parserOpt.add_argument('--scaleFactorBigwig',
                           help='Scale the values of a bigwig file by the given factor.',
                           type=float,
                           default=1.0)
    parserOpt.add_argument('--fontsize',
                           help='Fontsize in the plot for x and y axis.',
                           type=float,
                           default=10)
    parserOpt.add_argument('--rotationX',
                           help='Rotation in degrees for the labels of x axis.',
                           type=float,
                           default=0)
    parserOpt.add_argument('--rotationY',
                           help='Rotation in degrees for the labels of y axis.',
                           type=float,
                           default=0)
    parserOpt.add_argument('--increaseFigureWidth',
                           help='Plotting additional bigwig tracks can cause compression in x or y direction of the heatmap. Set this factor to a value unequal to 0 to correct this.',
                           type=float,
                           default=0.5)
    parserOpt.add_argument('--increaseFigureHeight',
                           help='Plotting additional bigwig tracks can cause compression in x or y direction of the heatmap. Set this factor to a value unequal to 0 to correct this.',
                           type=float,
                           default=0.5)
    parserOpt.add_argument('--loops',
                           help='Bedgraph file to plot detected long range contacts '
                           'from hicDetectLoops.',
                           type=str,
                           default=None)
    parserOpt.add_argument('--help', '-h', action='help',
                           help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    #  Used for automatic testing
    parserOpt.add_argument('--disable_tight_layout',
                           help=argparse.SUPPRESS,
                           action='store_true')

    return parser


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


# def getRegion(args, ma):
#     chrom = region_start = region_end = idx1 = start_pos1 = chrom2 = region_start2 = region_end2 = idx2 = start_pos2 = None
#     chrom, region_start, region_end = translate_region(args.region)

#     chrom = check_chrom_str_bytes(ma.interval_trees, chrom)
#     # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
#     #     chrom = toBytes(chrom)

#     if chrom not in list(ma.interval_trees):

#         chrom = change_chrom_names(chrom)

#         chrom = check_chrom_str_bytes(ma.interval_trees, chrom)

#         # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
#         #     chrom = toBytes(chrom)

#         if chrom not in list(ma.interval_trees):
#             exit("Chromosome name {} in --region not in matrix".format(change_chrom_names(chrom)))

#     args.region = [chrom, region_start, region_end]
#     is_cooler = check_cooler(args.matrix)
#     if is_cooler:
#         idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and  # noqa: W504
#                                  ((x[1] >= region_start and x[2] < region_end) or                            # noqa: W504
#                                   (x[1] < region_end and x[2] < region_end and x[2] > region_start) or       # noqa: W504
#                                   (x[1] > region_start and x[1] < region_end))])
#     else:
#         idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and  # noqa: W504
#                                  x[1] >= region_start and x[2] < region_end])
#     if hasattr(args, 'region2') and args.region2:
#         chrom2, region_start2, region_end2 = translate_region(args.region2)
#         chrom2 = check_chrom_str_bytes(ma.interval_trees, chrom2)

#         # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
#         #     chrom2 = toBytes(chrom)
#         if chrom2 not in list(ma.interval_trees):
#             chrom2 = change_chrom_names(chrom2)
#             chrom2 = check_chrom_str_bytes(ma.interval_trees, chrom2)

#             # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
#             #     chrom2 = toBytes(chrom)
#             if chrom2 not in list(ma.interval_trees):
#                 exit("Chromosome name {} in --region2 not in matrix".format(change_chrom_names(chrom2)))
#         if is_cooler:
#             idx2, start_pos2 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom2 and  # noqa: W504
#                                      ((x[1] >= region_start2 and x[2] < region_end2) or                           # noqa: W504
#                                       (x[1] < region_end2 and x[2] < region_end2 and x[2] > region_start2) or     # noqa: W504
#                                       (x[1] > region_start2 and x[1] < region_end2))])
#         else:
#             idx2, start_pos2 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom2 and  # noqa: W504
#                                      x[1] >= region_start2 and x[2] < region_end2])
#     else:
#         idx2 = idx1
#         chrom2 = chrom
#         start_pos2 = start_pos1

#     return chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2


def plotHeatmap(ma, chrBinBoundaries, fig, position, args, cmap, xlabel=None,
                ylabel=None, start_pos=None, start_pos2=None, pNorm=None, pAxis=None, pBigwig=None,
                pLoops=None, pHiCMatrix=None, pChromsomeStartEndDict=None, pResolution=None):
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
    img3 = axHeat2.pcolormesh(
        xmesh.T, ymesh.T, ma, vmin=args.vMin, vmax=args.vMax, cmap=cmap, norm=pNorm)
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

        plotLongRangeContacts(axHeat2, pLoops, pHiCMatrix,
                              args.region, args.chromosomeOrder)
        # pLongRangeContacts=None, pHiCMatrix=None
        # plotLongRangeContacts(pAxis, pNameOfLongRangeContactsFile, pHiCMatrix)
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
    # region_string = toBytes(region_string)
    region_string = region_string.replace(",", "")
    region_string = region_string.replace(";", "")
    region_string = region_string.replace("!", "")
    region_string = region_string.replace("-", ":")

    fields = region_string.split(":")
    chrom = fields[0]

    chrom = check_chrom_str_bytes(ma.interval_trees, chrom)
    if chrom not in list(ma.interval_trees):
        chrom = change_chrom_names(chrom)
        chrom = check_chrom_str_bytes(ma.interval_trees, chrom)
        if chrom not in list(ma.interval_trees):
            exit(
                "Chromosome name {} in --region not in matrix".format(change_chrom_names(chrom)))

    first_bin, last_bin = ma.getChrBinRange(chrom)
    try:
        region_start = int(fields[1])
    except IndexError:
        region_start = ma.getBinPos(first_bin)[1]
    try:
        region_end = int(fields[2])
    except IndexError:
        region_end = ma.getBinPos(last_bin - 1)[2]

    return chrom, region_start, region_end


def plotPerChr(hic_matrix, cmap, args, pBigwig, pResolution):
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
    if pBigwig:
        for i in range(len(args.bigwig)):
            fig_height += args.increaseFigureHeight
            # if args.bigwigAdditionalVerticalAxis:
            fig_width += args.increaseFigureWidth

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=args.dpi)

    chrom, start, end, _ = zip(*hic_matrix.cut_intervals)
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
        chrom_range = hic_matrix.getChrBinRange(chrname)
        matrix = np.asarray(hic_matrix.matrix[chrom_range[0]:chrom_range[1],
                                              chrom_range[0]:chrom_range[1]].todense().astype(float))

        norm = None
        if args.log or args.log1p:
            mask = matrix == 0
            mask_nan = np.isnan(matrix)
            mask_inf = np.isinf(matrix)
            log.debug("any nan {}".format(np.isnan(matrix).any()))
            log.debug("any inf {}".format(np.isinf(matrix).any()))

            try:
                matrix[mask] = np.nanmin(matrix[mask == False])
                matrix[mask_nan] = np.nanmin(matrix[mask_nan == False])
                matrix[mask_inf] = np.nanmin(matrix[mask_inf == False])

            except Exception:
                log.debug("Clearing of matrix failed.")
            log.debug("any nanafter remove of nan: {}".format(
                np.isnan(matrix).any()))
            log.debug("any inf after remove of inf: {}".format(
                np.isinf(matrix).any()))
        if args.log1p:
            matrix += 1
            norm = LogNorm()

        elif args.log:
            norm = LogNorm()

        chr_bin_boundary = OrderedDict()
        chr_bin_boundary[chrname] = hic_matrix.get_chromosome_sizes()[chrname]

        args.region = toString(chrname)
        chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(
            args, hic_matrix)
        plotHeatmap(matrix, chr_bin_boundary, fig, None,
                    args, cmap, xlabel=chrname, ylabel=chrname,
                    start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=axis, pBigwig=bigwig_info,
                    pChromsomeStartEndDict=chromosome_start_end(hic_matrix), pResolution=pResolution)
    return fig


def getRegion(args, ma):
    chrom = region_start = region_end = idx1 = start_pos1 = chrom2 = region_start2 = region_end2 = idx2 = start_pos2 = None
    chrom, region_start, region_end = translate_region(args.region, ma)

    args.region = [chrom, region_start, region_end]
    is_cooler = check_cooler(args.matrix)
    if is_cooler:
        idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and  # noqa: W504
                                 ((x[1] >= region_start and x[2] <= region_end) or                            # noqa: W504
                                  (x[1] < region_end and x[2] < region_end and x[2] > region_start) or       # noqa: W504
                                  (x[1] > region_start and x[1] < region_end))])
    else:
        idx1, start_pos1 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom and  # noqa: W504
                                 x[1] >= region_start and x[2] <= region_end])
    if hasattr(args, 'region2') and args.region2:
        chrom2, region_start2, region_end2 = translate_region(args.region2, ma)
        chrom2 = check_chrom_str_bytes(ma.interval_trees, chrom2)

        # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
        #     chrom2 = toBytes(chrom)
        if chrom2 not in list(ma.interval_trees):
            chrom2 = change_chrom_names(chrom2)
            chrom2 = check_chrom_str_bytes(ma.interval_trees, chrom2)

            # if type(next(iter(ma.interval_trees))) in [np.bytes_, bytes]:
            #     chrom2 = toBytes(chrom)
            if chrom2 not in list(ma.interval_trees):
                exit(
                    "Chromosome name {} in --region2 not in matrix".format(change_chrom_names(chrom2)))
        if is_cooler:
            idx2, start_pos2 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom2 and  # noqa: W504
                                     ((x[1] >= region_start2 and x[2] < region_end2) or                           # noqa: W504
                                      (x[1] < region_end2 and x[2] < region_end2 and x[2] > region_start2) or     # noqa: W504
                                      (x[1] > region_start2 and x[1] < region_end2))])
        else:
            idx2, start_pos2 = zip(*[(idx, x[1]) for idx, x in enumerate(ma.cut_intervals) if x[0] == chrom2 and  # noqa: W504
                                     x[1] >= region_start2 and x[2] < region_end2])
    else:
        idx2 = idx1
        chrom2 = chrom
        start_pos2 = start_pos1

    return chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2


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


def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.title:
        args.title = remove_non_ascii(args.title)

    chrom = None
    start_pos1 = None
    chrom2 = None
    start_pos2 = None

    if args.perChromosome and args.region:
        log.error('ERROR, choose from the option '
                  '--perChromosome or --region, the two '
                  'options at the same time are not '
                  'compatible.')
        exit(1)
    if args.bigwig is not None and len(args.bigwig) > 1 and args.bigwigAdditionalVerticalAxis:
        log.error(
            'Either multiple bigwig files on x axis or additional vertical axis are supported')
    # if args.region and args.region2 and args.bigwig:
    #     log.error("Inter-chromosomal pca is not supported.")
    #     exit(1)
    # is_cooler = False
    # if args.matrix.endswith('.cool') or cooler.io.is_cooler(args.matrix) or'.mcool' in args.matrix:
    is_cooler = check_cooler(args.matrix)
    log.info("Cooler or no cooler: {}".format(is_cooler))
    open_cooler_chromosome_order = True
    if args.chromosomeOrder is not None and len(args.chromosomeOrder) > 1:
        open_cooler_chromosome_order = False

    if is_cooler and not args.region2 and open_cooler_chromosome_order:
        log.debug("Retrieve data from cooler format and use its benefits.")
        regionsToRetrieve = None
        if args.region:
            regionsToRetrieve = []
            regionsToRetrieve.append(args.region)
            # if args.region2:
            #     chrom2, region_start2, region_end2 = translate_region(args.region2)
            #     regionsToRetrieve.append(args.region2)
        if args.chromosomeOrder:
            args.region = None
            args.region2 = None
            regionsToRetrieve = args.chromosomeOrder

        ma = HiCMatrix.hiCMatrix(args.matrix, pChrnameList=regionsToRetrieve)
        log.debug('Shape {}'.format(ma.matrix.shape))
        if args.clearMaskedBins:
            ma.maskBins(ma.nan_bins)
            # to avoid gaps in the plot, bins flanking the masked bins
            # are enlarged
            new_intervals = enlarge_bins(ma.cut_intervals)
            ma.setCutIntervals(new_intervals)

        if args.region:
            chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(
                args, ma)

        matrix = np.asarray(ma.matrix.todense().astype(float))
        matrix_length = len(matrix[0])
        log.debug("Number of data points matrix_cool: {}".format(matrix_length))
    else:
        ma = HiCMatrix.hiCMatrix(args.matrix)
        if args.clearMaskedBins:
            ma.maskBins(ma.nan_bins)
            new_intervals = enlarge_bins(ma.cut_intervals)
            ma.setCutIntervals(new_intervals)
        if args.chromosomeOrder:
            args.region = None
            args.region2 = None

            valid_chromosomes = []
            invalid_chromosomes = []
            log.debug('args.chromosomeOrder: {}'.format(args.chromosomeOrder))
            log.debug("ma.chrBinBoundaries {}".format(ma.chrBinBoundaries))

            args.chromosomeOrder = toBytes(args.chromosomeOrder)
            for chrom in toString(args.chromosomeOrder):
                if chrom in ma.chrBinBoundaries:
                    valid_chromosomes.append(chrom)
                else:
                    invalid_chromosomes.append(chrom)

            if len(invalid_chromosomes) > 0:
                log.warning("WARNING: The following chromosome/scaffold names were not found. Please check"
                            "the correct spelling of the chromosome names. \n")
                log.warning("\n".join(invalid_chromosomes))
            ma.reorderChromosomes(valid_chromosomes)
            chrom = None
        log.info("min: {}, max: {}\n".format(
            ma.matrix.data.min(), ma.matrix.data.max()))

        if args.region:
            chrom, region_start, region_end, idx1, start_pos1, chrom2, region_start2, region_end2, idx2, start_pos2 = getRegion(
                args, ma)

            matrix = np.asarray(
                ma.matrix[idx1, :][:, idx2].todense().astype(float))

        else:
            log.debug("Else branch")
            matrix = np.asarray(ma.getMatrix().astype(float))

    resolution = ma.getBinSize()
    matrix_length = len(matrix[0])
    log.debug("Number of data points matrix: {}".format(matrix_length))

    for matrix_ in matrix:
        if not matrix_length == len(matrix_):
            log.error("Matrices do not have the same length: {} , {}".format(
                matrix_length, len(matrix_)))

    cmap = cm.get_cmap(args.colorMap)
    log.debug("Nan values set to black\n")
    cmap.set_bad('black')

    bigwig_info = None
    if args.bigwig:
        bigwig_info = {'args': args, 'axis': None,
                       'axis_colorbar': None, 'nan_bins': ma.nan_bins}

    if args.perChromosome:
        log.debug('583')
        fig = plotPerChr(ma, cmap, args, pBigwig=bigwig_info, pResolution=resolution)

    else:
        norm = None

        if args.log or args.log1p:
            mask = matrix == 0
            try:
                matrix[mask] = np.nanmin(matrix[mask == False])
            except ValueError:
                log.info('Matrix contains only 0. Set all values to {}'.format(
                    np.finfo(float).tiny))
                matrix[mask] = np.finfo(float).tiny
            if np.isnan(matrix).any() or np.isinf(matrix).any():
                log.debug("any nan {}".format(np.isnan(matrix).any()))
                log.debug("any inf {}".format(np.isinf(matrix).any()))
                mask_nan = np.isnan(matrix)
                mask_inf = np.isinf(matrix)
                matrix[mask_nan] = np.nanmin(matrix[mask_nan == False])
                matrix[mask_inf] = np.nanmin(matrix[mask_inf == False])

        log.debug("any nan after remove of nan: {}".format(
            np.isnan(matrix).any()))
        log.debug("any inf after remove of inf: {}".format(
            np.isinf(matrix).any()))
        if args.log1p:
            matrix += 1
            norm = LogNorm()
        elif args.log:
            norm = LogNorm()

        fig_height = 7
        fig_width = 8

        if args.bigwig:
            # increase figure height to accommodate bigwig track

            for i in range(len(args.bigwig)):
                fig_height += args.increaseFigureHeight
                # if args.bigwigAdditionalVerticalAxis:
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

        if start_pos1 is None:
            start_pos1 = make_start_pos_array(ma)

        position = [left_margin, bottom, width, height]
        plotHeatmap(matrix, ma.get_chromosome_sizes(), fig, position,
                    args, cmap, xlabel=chrom, ylabel=chrom2,
                    start_pos=start_pos1, start_pos2=start_pos2, pNorm=norm, pAxis=ax1, pBigwig=bigwig_info,
                    pLoops=args.loops, pHiCMatrix=ma, pChromsomeStartEndDict=chromosome_start_end(ma), pResolution=resolution)

    if not args.disable_tight_layout:
        if args.perChromosome or args.bigwig:
            try:
                plt.tight_layout()
            except UserWarning:
                log.info("Failed to tight layout. Using regular plot.")
            except ValueError:
                log.info("Failed to tight layout. Using regular plot.")
    # plt.setp(bigwig_vertical_axis.get_xticklabels(), rotation=180)
    plt.savefig(args.outFileName, dpi=args.dpi)
    plt.close(fig)


def chromosome_start_end(pHiCMatrix):
    chromosome_start_end_dict = {}

    for chromosome, (start_bin, end_bin) in pHiCMatrix.chrBinBoundaries.items():
        chromosome_start_end_dict[chromosome] = (pHiCMatrix.cut_intervals[start_bin][1], pHiCMatrix.cut_intervals[end_bin - 1][2])

    return chromosome_start_end_dict


def make_start_pos_array(ma):
    # makes an start_pos array that can be used
    # to plot the bins of the matrix using their real length
    # When the whole matrix wants to be plotted, the start_pos needs to be modified
    # such that at each chromosome start, the start_pos does not go back to zero and instead
    # is added
    chrom_sizes = dict()
    for chr, v in chromosome_start_end(ma).items():
        start, end = v
        chrom_sizes[chr] = end - start
    # chrom_sizes = ma.get_chromosome_sizes() ## TODO this is incorrect!!
    prev_chrom = ma.cut_intervals[0][0]
    prev_chroms_sum = 0
    start_pos = []
    shift = 0
    for index, (chrom, start, end, _) in enumerate(ma.cut_intervals): # Shifts all the coords to start from 0
        if index == 0:
            if start != 0:
                shift = start

        if chrom != prev_chrom:
            prev_chroms_sum += chrom_sizes[prev_chrom]
            prev_chrom = chrom
            if start != 0: # shift all the other chrosmomes also back to be started just after the previous chromosome
                shift = start
            else:
                shift = 0

        start_pos.append(start - shift + prev_chroms_sum)
    return start_pos


def plotBigwig(pAxis, pNameOfBigwigList, pChromosomeSizes=None, pRegion=None, pXticks=None,
               pFlipBigwigSign=None, pScaleFactorBigwig=None, pVertical=False,
               pValueMin=None, pValueMax=None, pResolution=None):
    log.debug('plotting bigwig file')

    # pNameOfBigwigList is not a list, but to make room for future options
    # requiring more than one bigwig file I set this to a list intentionally.
    # pNameOfBigwigList = [pNameOfBigwigList]
    for file in pNameOfBigwigList:
        file_format = file.split(".")[-1]
        if file_format != 'bigwig' and file_format != 'bw':
            log.error("Given files are not bigwig")
            exit(1)

    # for bigwig_file in pNameOfBigwigList:
    #     if bigwig_file.split('.')[-1] != file_format:
    #         log.error("Eigenvector input files have different formats.")
    #         exit()

    if file_format == "bigwig" or file_format == 'bw':
        for i, bigwigFile in enumerate(pNameOfBigwigList):
            x_values = []
            bigwig_scores = []
            pAxis[i].set_frame_on(False)

            if pVertical:

                pAxis[i].yaxis.set_visible(False)

            else:
                # pAxis[i].set_frame_on(False)

                pAxis[i].xaxis.set_visible(False)
            bw = pyBigWig.open(bigwigFile)
            bigwig_scores = []
            if pRegion:
                chrom, region_start, region_end = pRegion
                # region_end could be a very large number returned by translate_region
                region_end = min(region_end, pChromosomeSizes[chrom][1])
                # log.info("chromosomes bigwig: {}".format(bw.chroms()))
                chrom = check_chrom_str_bytes(bw.chroms(), chrom)
                if chrom not in list(bw.chroms().keys()):
                    chrom = change_chrom_names(chrom)
                    if chrom not in list(bw.chroms().keys()):
                        log.info(
                            "bigwig file has no chromosome named: {}.".format(chrom))
                        return

                # the bigwig file may end before the region end, to avoid and error
                # the bigwig_end is set for the pyBigwig query
                bigwig_end = min(bw.chroms()[chrom], region_end)

                # TODO, this could be a parameter
                num_bins = int(bigwig_end - region_start) // pResolution
                log.debug('chrom {}, region_start {}, bigwig_end {}, num_bins {}'.format(chrom, region_start, bigwig_end, num_bins))
                scores_per_bin = np.array(
                    bw.stats(chrom, region_start, bigwig_end, nBins=num_bins)).astype(float)
                if scores_per_bin is None:
                    log.info(
                        "Chromosome {} has no entries in bigwig file.".format(chrom))
                    return

                _x_vals = np.linspace(region_start, region_end, num_bins)
                log.debug('_x_vals[:10] {}'.format(_x_vals[:10]))
                assert len(_x_vals) == len(scores_per_bin)
                x_values.extend(_x_vals)
                bigwig_scores.extend(scores_per_bin)
                if pVertical:
                    pAxis[i].set_ylim(region_start, region_end)
                else:
                    log.debug('limits region: {}, {}'.format(region_start, region_end))

                    pAxis[i].set_xlim(region_start, region_end)

            elif pChromosomeSizes:
                log.debug('pChromosomeSizes: {}'.format(pChromosomeSizes))
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
                    # chrom = check_chrom_str_bytes(pChromosomeSizes, chrom)
                    # set the bin size to approximately 100kb
                    # or to the chromosome size if this happens to be less than 100kb
                    # chunk_size = min(1e6, pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0])
                    # num_bins = int(pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0] / chunk_size)
                    bigwig_end = min(bw.chroms()[chrom_], pChromosomeSizes[chrom][1])
                    num_bins = int(bigwig_end - pChromosomeSizes[chrom][0]) // pResolution
                    log.debug('chrom {}, region_start {}, bigwig_end {}, num_bins {}'.format(chrom_, pChromosomeSizes[chrom][0], pChromosomeSizes[chrom][1], num_bins))

                    scores_per_bin = np.array(
                        bw.stats(chrom_, pChromosomeSizes[chrom][0], bigwig_end, nBins=num_bins)).astype(float)

                    if scores_per_bin is None:
                        log.info(
                            "Chromosome {} has no entries in bigwig file.".format(chrom))
                        return

                    _x_vals = np.linspace(
                        chrom_length_sum + pChromosomeSizes[chrom][0], chrom_length_sum + pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0], num_bins)
                    log.debug('_x_vals[:10] {}'.format(_x_vals[:10]))

                    assert len(_x_vals) == len(scores_per_bin)
                    x_values.extend(_x_vals)

                    # log.debug('x_values {}'.format(x_values))
                    bigwig_scores.extend(scores_per_bin)

                    log.debug('len bigwig score {}, x_values {}'.format(len(bigwig_scores), len(x_values)))

                    chrom_length_sum += pChromosomeSizes[chrom][1] - pChromosomeSizes[chrom][0]
                    if min_start is None:
                        min_start = pChromosomeSizes[chrom][0]
                    elif min_start > pChromosomeSizes[chrom][0]:
                        min_start = pChromosomeSizes[chrom][0]
                # min_start = 0
                if pVertical:
                    pAxis[i].set_ylim(min_start, chrom_length_sum)
                else:
                    log.debug('limits: {}, {}'.format(min_start, chrom_length_sum))
                    pAxis[i].set_xlim(min_start, chrom_length_sum)

            log.debug("Number of data points: {}".format(len(bigwig_scores)))
            bigwig_scores = np.array(bigwig_scores)
            if pFlipBigwigSign:
                log.info("Flipping sign of bigwig values.")
                bigwig_scores *= -1
            if pScaleFactorBigwig is not None and pScaleFactorBigwig != 1.0:
                log.info("Scaling bigwig values.")
                bigwig_scores *= pScaleFactorBigwig
            if pValueMin is not None or pValueMax is not None:
                bigwig_scores = bigwig_scores.clip(pValueMin, pValueMax)

            if x_values is not None and bigwig_scores is not None:
                log.debug('Set values')
                if pVertical:
                    log.debug('Set values vertical')

                    pAxis[i].fill_between(
                        np.flip(bigwig_scores, 0), x_values, edgecolor='none')
                else:
                    log.debug('Set values horizontal')

                    pAxis[i].fill_between(
                        x_values, 0, bigwig_scores, edgecolor='none')


def plotLongRangeContacts(pAxis, pNameOfLongRangeContactsFile, pHiCMatrix, pRegion, pChromosomeOrder):

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

                x = int(start_X)
                y = int(start_Y)
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
