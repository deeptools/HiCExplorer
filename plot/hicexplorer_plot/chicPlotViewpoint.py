"""The figures of chicPlotViewpoint.

The C++ binary reads the interaction, background, differential and
significant interaction files and computes, per plot, what
Viewpoint.getDataForPlotting, readRejectedFile and readSignificantRegionsFile
return in hicexplorer/chicPlotViewpoint.py plot_images(), the --xFold band
of the background model, and the p-value row
plot_images() draws for --pValue (float32, --truncateZeroPvalues, the
significant regions' p-values, the clip) with its colorbar ticks. This module
runs the drawing part of plot_images() with the reference's calls
(plotViewpoint, plotBackgroundModel and plotPValue of
hicexplorer/lib/viewpoint.py are reproduced below), once per worker chunk in
the reference's order, and writes image i of the flattened chunks to
imageDirectory/<i>.<outputFormat>; the C++ binary writes the archive from
them as main() does.

Data:
    outputFormat, dpi, range, resolution, colorList, xFold, pValue,
    colorMapPvalue, pValueSignificanceLevels, imageDirectory
    chunks: [[group]], one list per worker process
    group: {"file_name", "items": [item]}
    item: {"label", "skip", "data", "background", "background_upper" (the
           background times --xFold, null without a non-zero --xFold),
           "p_values", "p_value_ticks", "highlight", "significant_regions",
           "significant_p_values"}
"""

import io

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas  # noqa: E402
import matplotlib as mpl  # noqa: E402
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402


def plotViewpoint(pAxis, pData, pColor, pLabelName, pHighlightRegion=None, pHighlightSignificantRegion=None):
    data_plot_label = pAxis.plot(
        range(len(pData)), pData, '-' + pColor, alpha=0.9, label=pLabelName, linewidth=1)
    if pHighlightRegion:
        for region in pHighlightRegion:
            pAxis.axvspan(region[0], region[1], color='red', alpha=0.3)
    if pHighlightSignificantRegion:
        for region in pHighlightSignificantRegion:
            pAxis.axvspan(region[0], region[1], color=pColor, alpha=0.3)
    return data_plot_label


def plotBackgroundModel(pAxis, pBackgroundData, pUpperValues=None):
    pBackgroundData = np.array(pBackgroundData)
    data_plot_label = pAxis.plot(range(len(pBackgroundData)), pBackgroundData, '-r', alpha=0.5, label='background model', linewidth=1)
    if pUpperValues is not None:
        upper_values = np.array(pUpperValues)
        lower_values = pBackgroundData
        pAxis.fill_between(range(len(pBackgroundData)), upper_values, lower_values, facecolor='r', alpha=0.5)
    return data_plot_label


def plotPValue(pAxis, pAxisLabel, pPValueData, pTicks, pLabelText, pCmap, pFigure, pValueSignificanceLevels):

    _z_score = np.empty([2, len(pPValueData)])
    _z_score[:, :] = pPValueData
    pAxis.xaxis.set_visible(False)
    pAxis.yaxis.set_visible(False)
    divider = make_axes_locatable(pAxisLabel)
    cax = divider.append_axes("left", size="20%", pad=0.09)

    if pPValueData is not None:
        img = pAxis.contourf(_z_score, cmap=pCmap)
        colorbar = pFigure.colorbar(
            img, cax=cax, ticks=pTicks)
        colorbar.ax.set_ylabel('p-value', size=6)

    elif pValueSignificanceLevels:
        pValueSignificanceLevels.insert(0, -1)
        pValueSignificanceLevels.append(1)

        img = pAxis.contourf(_z_score, levels=pValueSignificanceLevels, colors=['#CC0000', '#FFD43B', '#306998', '#FFFFFF'])
        colorbar = pFigure.colorbar(img, cax=cax, ticks=[pValueSignificanceLevels[1], pValueSignificanceLevels[2], pValueSignificanceLevels[3]])
        colorbar.ax.tick_params(labelsize=6)
        colorbar.ax.set_ylabel('p-value', size=6)

    pAxisLabel.text(0.45, 0, pLabelText, size=7)
    pAxisLabel.xaxis.set_visible(False)
    pAxisLabel.yaxis.set_visible(False)
    pAxisLabel.set_frame_on(False)


def plot_images(groups, args):
    images_array = []
    file_name_list = []
    pRange = args["range"]
    pResolution = args["resolution"]
    for group in groups:
        interactionFile = group["items"]
        number_of_rows_plot = len(interactionFile)
        matplotlib.rcParams.update({'font.size': 9})
        fig = plt.figure(figsize=(9.4, 4.8), dpi=args["dpi"])
        FigureCanvas(fig)
        z_score_heights = [0.07] * number_of_rows_plot
        viewpoint_height_ratio = 0.95 - (0.07 * number_of_rows_plot)
        if viewpoint_height_ratio < 0.4:
            viewpoint_height_ratio = 0.4
            _ratio = 0.6 / number_of_rows_plot
            z_score_heights = [_ratio] * number_of_rows_plot

        if args["pValue"]:
            gs = gridspec.GridSpec(1 + len(interactionFile), 2, height_ratios=[0.95 - (0.07 * number_of_rows_plot), *z_score_heights], width_ratios=[0.75, 0.25])
            gs.update(hspace=0.5, wspace=0.05)
            ax1 = plt.subplot(gs[0, 0])
            ax1.margins(x=0)
        else:
            ax1 = plt.gca()
        colors = args["colorList"]
        background_plot = True
        data_plot_label = None

        for i, item in enumerate(interactionFile):
            if item["skip"]:
                continue
            data = item["data"]
            background_data_plot = item["background"]
            p_values = item["p_values"]
            p_value_ticks = item["p_value_ticks"]
            highlight_differential_regions = item["highlight"]
            significant_regions = item["significant_regions"]
            if data_plot_label:
                data_plot_label += plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=item["label"], pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)
            else:
                data_plot_label = plotViewpoint(pAxis=ax1, pData=data, pColor=colors[i % len(colors)], pLabelName=item["label"], pHighlightRegion=highlight_differential_regions, pHighlightSignificantRegion=significant_regions)

            if background_plot:
                if background_data_plot is not None:
                    data_plot_label += plotBackgroundModel(pAxis=ax1, pBackgroundData=background_data_plot,
                                                           pUpperValues=item["background_upper"])
                background_plot = False
            if args["pValue"]:
                plotPValue(pAxis=plt.subplot(gs[1 + i, 0]), pAxisLabel=plt.subplot(gs[1 + i, 1]), pPValueData=p_values,
                           pTicks=p_value_ticks, pLabelText=item["label"], pCmap=args["colorMapPvalue"],
                           pFigure=fig, pValueSignificanceLevels=args["pValueSignificanceLevels"])

        if data_plot_label is not None:

            ticks = []
            x_labels = []

            if pRange[0] + pRange[1] <= 2e6:
                divisor_legend = 1e3
                mod_legend = 2e5

                if pRange[0] + pRange[1] <= 1e4:
                    mod_legend = 5e3
                elif pRange[0] + pRange[1] <= 5e4:
                    mod_legend = 1e4
                elif pRange[0] + pRange[1] <= 1e5:
                    mod_legend = 5e4
                elif pRange[0] + pRange[1] <= 5e5:
                    mod_legend = 1e5

                unit = 'kb'
            elif pRange[0] + pRange[1] > 2e6:
                divisor_legend = 1e6
                mod_legend = 1e6
                unit = 'Mb'

            for k, j in zip(range((pRange[0])), range(pRange[0], 1, -1)):
                if j % mod_legend == 0:
                    x_labels.append(str(-int(j) // int(divisor_legend)) + unit)
                    ticks.append(k // pResolution)
            x_labels.append('RP')
            ticks.append(pRange[0] // pResolution)

            referencepoint_index = ticks[-1]
            for k, j in zip(range(pRange[1]), range(1, pRange[1] + 1, 1)):
                if j % mod_legend == 0:
                    x_labels.append(str(int(j) // int(divisor_legend)) + unit)
                    ticks.append(referencepoint_index + (k // pResolution))

            ax1.set_ylabel('Number of interactions')
            ax1.set_xticks(ticks)
            ax1.set_xticklabels(x_labels)

            data_legend = [label.get_label() for label in data_plot_label]
            ax1.legend(data_plot_label, data_legend, loc=0)

            bufferObject = io.BytesIO()
            plt.savefig(bufferObject, format=args["outputFormat"], dpi=300)
            images_array.append(bufferObject)
        plt.close(fig)
        file_name_list.append(group["file_name"])
    return images_array, file_name_list


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    images_array = []
    for chunk in data["chunks"]:
        images, _names = plot_images(chunk, data)
        images_array.append(images)

    images_array = [item for sublist in images_array for item in sublist]
    for i, bufferObject in enumerate(images_array):
        path = "{}/{}.{}".format(data["imageDirectory"], i, data["outputFormat"])
        with open(path, "wb") as handle:
            handle.write(bufferObject.getvalue())
        bufferObject.close()
