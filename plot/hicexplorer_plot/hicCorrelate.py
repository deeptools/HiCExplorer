"""The figures of hicCorrelate: the scatter grid and the clustered heatmap.

The C++ binary computes the per matrix vectors (as big_mat.T in an npz), the
pairwise correlations and the complete linkage of the correlation matrix.
The calls are those of hicexplorer/hicCorrelate.py main() from the GridSpec
to plot_correlation, unchanged, except that the correlation values and the
linkage come from the data instead of scipy.stats and
scipy.cluster.hierarchy.linkage; the dendrogram is drawn by scipy from that
linkage.

Data:
    outFileNameHeatmap, outFileNameScatter, labels, method, log1p, zMin, zMax,
    colorMap, plotNumbers
    results   the correlation matrix, upper triangle and diagonal filled
    linkage   the (n - 1) by 4 linkage matrix, row major
    leaves    the dendrogram leaf order (for checks; scipy derives the same)
    big_mat   path of an npz holding big_mat.T
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
from matplotlib import colormaps as cm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
from matplotlib.ticker import FixedLocator  # noqa: E402
import matplotlib.colors as pltcolors  # noqa: E402


def plot_correlation(corr_matrix, labels, plot_filename, linkage, vmax=None,
                     vmin=None, colormap='Reds', pPlotNumbers=None):
    import scipy.cluster.hierarchy as sch
    num_rows = corr_matrix.shape[0]

    if vmax is None:
        vmax = 1
    if vmin is None:
        vmin = 0 if corr_matrix.min() >= 0 else -1

    fig = plt.figure(figsize=(10.5, 9.5))
    axdendro = fig.add_axes([0.02, 0.1, 0.1, 0.7])
    axdendro.set_axis_off()
    y_var = linkage
    z_var = sch.dendrogram(y_var, orientation='left',
                           link_color_func=lambda k: 'black')
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    cmap = cm.get_cmap(colormap)

    if pPlotNumbers:
        cmap = pltcolors.LinearSegmentedColormap.from_list(colormap + "clipped",
                                                           cmap(np.linspace(0, 0.9, 10)))
    axmatrix = fig.add_axes([0.13, 0.1, 0.6, 0.7])
    index = z_var['leaves']
    corr_matrix = corr_matrix[index, :]
    corr_matrix = corr_matrix[:, index]
    img_mat = axmatrix.matshow(corr_matrix, aspect='equal', origin='lower',
                               cmap=cmap, extent=(0, num_rows, 0, num_rows),
                               vmax=vmax, vmin=vmin)
    axmatrix.yaxis.tick_right()
    axmatrix.set_yticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_yticklabels(np.array(labels).astype('str')[index],
                             fontsize=14)

    axmatrix.set_xticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_xticklabels(np.array(labels).astype('str')[index],
                             fontsize=14,
                             rotation=45,
                             ha='left')

    axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
    plt.colorbar(img_mat, cax=axcolor, orientation='horizontal')

    if pPlotNumbers:
        for row in range(num_rows):
            for col in range(num_rows):
                axmatrix.text(row + 0.5, col + 0.5,
                              "{:.2f}".format(corr_matrix[row, col]),
                              ha='center', va='center')

    fig.savefig(plot_filename)


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    labels = data["labels"]
    num_files = len(labels)
    results = np.array(data["results"], dtype='float')
    big_mat = np.ma.masked_invalid(np.load(data["big_mat"]))

    rows, cols = np.triu_indices(num_files)

    grids = gridspec.GridSpec(num_files, num_files)
    grids.update(wspace=0, hspace=0)
    fig = plt.figure(figsize=(2 * num_files, 2 * num_files))
    plt.rcParams['font.size'] = 8.0

    min_value = int(big_mat.min())
    max_value = int(big_mat.max())
    if (min_value % 2 == 0 and max_value % 2 == 0) or \
            (min_value % 1 == 0 and max_value % 2 == 1):
        max_value += 1

    if data["log1p"]:
        major_locator = FixedLocator(list(range(min_value, max_value, 2)))
        minor_locator = FixedLocator(list(range(min_value, max_value, 1)))

    for index in range(len(rows)):
        row = rows[index]
        col = cols[index]
        if row == col:
            ax = fig.add_subplot(grids[row, col])
            ax.text(0.6, 0.6, labels[row],
                    verticalalignment='center',
                    horizontalalignment='center',
                    fontsize=10, fontweight='bold',
                    transform=ax.transAxes)
            ax.set_axis_off()
            continue

        _mat = big_mat[:, [row, col]]
        _mat = _mat[_mat.sum(axis=1) > 1, :]
        vector1 = _mat[:, 0]
        vector2 = _mat[:, 1]

        ax = fig.add_subplot(grids[row, col])
        if data["log1p"]:
            ax.xaxis.set_major_locator(major_locator)
            ax.xaxis.set_minor_locator(minor_locator)
            ax.yaxis.set_major_locator(major_locator)
            ax.yaxis.set_minor_locator(minor_locator)

        ax.text(0.2, 0.8, "{}={:.2f}".format(data["method"],
                                             results[row, col]),
                horizontalalignment='left',
                transform=ax.transAxes)
        ax.get_yaxis().set_tick_params(
            which='both',
            left='off',
            right='off',
            direction='out')

        ax.get_xaxis().set_tick_params(
            which='both',
            top='off',
            bottom='off',
            direction='out')

        if col != num_files - 1:
            ax.set_yticklabels([])
        else:
            ax.yaxis.tick_right()
            ax.get_yaxis().set_tick_params(
                which='both',
                left='off',
                right='on',
                direction='out')
        if col - row == 1:
            ax.xaxis.tick_bottom()
            ax.get_xaxis().set_tick_params(
                which='both',
                top='off',
                bottom='on',
                direction='out')
        else:
            ax.set_xticklabels([])

        ax.hist2d(vector1, vector2, bins=150, cmin=0.1)
    fig.tight_layout()
    fig.savefig(data["outFileNameScatter"], bbox_inches='tight')

    results = results + np.triu(results, 1).T
    linkage = np.array(data["linkage"], dtype=np.float64).reshape(-1, 4)
    plot_correlation(results, labels,
                     data["outFileNameHeatmap"],
                     linkage,
                     data["zMax"],
                     data["zMin"],
                     data["colorMap"],
                     pPlotNumbers=data["plotNumbers"])
