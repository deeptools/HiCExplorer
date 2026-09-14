"""The box plot of hicPlotSVL.

The calls are those of hicexplorer/hicPlotSVL.py main() from
plt.ylabel to plt.savefig, unchanged. The C++ binary passes, per matrix, the
short range / long range ratio of every chromosome the reference keeps. A
ratio of two float32 sums is a numpy float32 there, so boxplot's statistics
are computed on a float32 array; those samples are marked and rebuilt as
float32.

Data:
    plotFileName   --plotFileName
    dpi            --dpi
    colorList      --colorList
    matrices       --matrices, as given
    samples        [{"values": ratios, "float32": bool}] per matrix
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402


def _sample(sample):
    dtype = np.float32 if sample["float32"] else np.float64
    return [dtype(value) for value in sample["values"]]


def draw(data):
    short_v_long_range = [_sample(sample) for sample in data["samples"]]
    color_list = data["colorList"]
    matrices = data["matrices"]

    plt.ylabel('Sum short range / long range')
    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    box_plot = plt.boxplot(short_v_long_range, patch_artist=True)
    legend_handels_color = []
    for i, patch in enumerate(box_plot['boxes']):
        patch.set_facecolor(color_list[i % len(color_list)])
        legend_handels_color.append(mpatches.Patch(color=color_list[i % len(color_list)], label=matrices[i].split('/')[-1]))
    plt.legend(handles=legend_handels_color)
    plt.xlabel('Boxplot shows svl-ratio per chromosome.')
    plt.savefig(data["plotFileName"], dpi=data["dpi"])
