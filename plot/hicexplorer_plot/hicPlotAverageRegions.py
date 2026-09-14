"""The figure of hicPlotAverageRegions.

The data is the npz file hicAverageRegions writes; the C++ binary checks it
and passes its path. The calls are those of
hicexplorer/hicPlotAverageRegions.py main() from `load_npz` to
`plt.savefig`, unchanged, including the 45 degree scipy.ndimage.rotate.

Data:
    matrix, outputFile, log1p, log, colorMap, vMin, vMax, dpi
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import rotate  # noqa: E402
from scipy.sparse import load_npz  # noqa: E402
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42

    matrix = load_npz(data["matrix"])

    matrix = matrix.toarray()
    matrix = np.triu(matrix)
    matrix = rotate(matrix, 45, cval=np.nan)
    matrix_shapes = matrix.shape
    matrix = matrix[:matrix_shapes[0] // 2, :]
    if data["log1p"]:
        matrix += 1

    fig = plt.figure()
    axis = plt.gca()
    if data["log"]:
        norm = LogNorm(vmin=data["vMin"], vmax=data["vMax"])
    elif data["log1p"]:
        if data["vMin"] is not None:
            vMin = data["vMin"] + 1
        else:
            vMin = None
        if data["vMax"] is not None:
            vMax = data["vMax"] + 1
        else:
            vMax = None
        norm = LogNorm(vmin=vMin, vmax=vMax)
    else:
        norm = matplotlib.colors.Normalize(vmin=data["vMin"], vmax=data["vMax"])

    matrix_axis = axis.matshow(matrix, cmap=data["colorMap"], norm=norm)
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)

    fig.colorbar(matrix_axis, cax=cax)
    plt.tight_layout()
    plt.savefig(data["outputFile"], dpi=data["dpi"])
