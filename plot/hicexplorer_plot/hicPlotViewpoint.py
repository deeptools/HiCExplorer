"""The figure of hicPlotViewpoint.

The calls are those of hicexplorer/hicPlotViewpoint.py main() from
`fig = plt.figure(figsize=(6.4, 4.8))` to `plt.close(fig)`, unchanged. The
C++ binary computes getViewpointValues for every matrix; the tick positions
use the bin ranges of the last matrix, as the reference's loop variables do.

Data:
    outFileName        the figure file
    dpi                --dpi
    referencePoint     the reference point split at ':' after the reference's
                       character replacements, 2 or 3 strings
    region_start, region_end
    view_point_start, view_point_end, view_point_range   of the last matrix
    data               one list of summed interactions per matrix
    legend             os.path.basename of every matrix
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402
import numpy as np  # noqa: E402


def relabelTicks(pTick):
    if pTick < 1e6:
        xlabels = "{:.2f} Kb".format(int(pTick) / 1e3)
    else:
        xlabels = "{:.2f} Mb".format(int(pTick) / 1e6)
    return xlabels


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    referencePoint = data["referencePoint"]
    region_start = data["region_start"]
    region_end = data["region_end"]
    view_point_start = data["view_point_start"]
    view_point_end = data["view_point_end"]
    view_point_range = data["view_point_range"]
    data_list = [np.array(values, dtype=np.float64) for values in data["data"]]
    matrix_name_legend = data["legend"]

    fig = plt.figure(figsize=(6.4, 4.8))
    ax = plt.subplot(111)
    matrices_plot_legend = []
    for i, values in enumerate(data_list):
        matrices_plot_legend.append(ax.plot(range(len(values)), values, alpha=0.7, label=matrix_name_legend[i])[0])
    if len(referencePoint) == 2:
        ax.set_xticks([0, view_point_start - view_point_range[0], view_point_range[1] - view_point_range[0]])
        xticklabels = [None] * 3
        xticklabels[0] = relabelTicks((int(referencePoint[1]) - region_start) * (-1))
        xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
        xticklabels[2] = relabelTicks(region_end - int(referencePoint[1]))

    elif len(referencePoint) == 3:
        ax.set_xticks([0, view_point_start - view_point_range[0], view_point_end - view_point_range[0], view_point_range[1] - view_point_range[0]])
        xticklabels = [None] * 4
        xticklabels[0] = relabelTicks((int(referencePoint[1]) - region_start) * (-1))
        xticklabels[1] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[1]))
        xticklabels[2] = referencePoint[0] + ":" + relabelTicks(int(referencePoint[2]))
        xticklabels[3] = relabelTicks(region_end - int(referencePoint[1]))

    ax.set_xticklabels(xticklabels)
    ax.set_ylabel('Number of interactions')

    plt.legend(handles=matrices_plot_legend)
    plt.savefig(data["outFileName"], dpi=data["dpi"])
    plt.close(fig)
