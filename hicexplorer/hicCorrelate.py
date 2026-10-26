"""hicCorrelate: a thin layer over the C++ executable hicCorrelate.

main() runs the executable and parse_arguments() is built from its --help-json
(hicexplorer._cpp). The other names of this module are kept for API
compatibility and are not used by the tool.
"""
import warnings
warnings.simplefilter(action="ignore", category=UserWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import argparse
import os
import numpy as np
from past.builtins import map
from scipy.sparse import triu
from scipy.stats import pearsonr, spearmanr

from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import check_cooler
# for plotting
from matplotlib import use as mplt_use
import matplotlib as mpl
from matplotlib import colormaps as cm

mplt_use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator
import matplotlib.colors as pltcolors
import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """The argparse parser of hicCorrelate, built from the C++ executable's --help-json."""
    from hicexplorer._cpp import build_parser
    return build_parser('hicCorrelate')


def plot_correlation(corr_matrix, labels, plot_filename, vmax=None,
                     vmin=None, colormap='Reds', pPlotNumbers=None):
    import scipy.cluster.hierarchy as sch
    num_rows = corr_matrix.shape[0]

    # set the minimum and maximum values
    if vmax is None:
        vmax = 1
    if vmin is None:
        vmin = 0 if corr_matrix.min() >= 0 else -1

    # Compute and plot dendrogram.
    fig = plt.figure(figsize=(10.5, 9.5))
    axdendro = fig.add_axes([0.02, 0.1, 0.1, 0.7])
    axdendro.set_axis_off()
    y_var = sch.linkage(corr_matrix, method='complete')
    z_var = sch.dendrogram(y_var, orientation='left',
                           link_color_func=lambda k: 'black')
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    cmap = cm.get_cmap(colormap)

    # this line simply makes a new cmap, based on the original
    # colormap that goes from 0.0 to 0.9
    # This is done to avoid colors that
    # are too dark at the end of the range that do not offer
    # a good contrast between the correlation numbers that are
    # plotted on black.
    # cmap = cmap.from_list(colormap + "clipped", cmap([0.0, 0.8]))
    if pPlotNumbers:
        cmap = pltcolors.LinearSegmentedColormap.from_list(colormap + "clipped",
                                                           cmap(np.linspace(0, 0.9, 10)))
    # Plot distance matrix.
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

#    axmatrix.set_xticks([])
    # Plot colorbar.
    axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
    plt.colorbar(img_mat, cax=axcolor, orientation='horizontal')

    if pPlotNumbers:
        for row in range(num_rows):
            for col in range(num_rows):
                axmatrix.text(row + 0.5, col + 0.5,
                              "{:.2f}".format(corr_matrix[row, col]),
                              ha='center', va='center')

    fig.savefig(plot_filename)


def get_vectors(mat1, mat2):
    """
    Uses sparse matrix tricks to convert
    into a vector the matrix values such
    that zero values that appear in only
    one of the matrices is kept. But
    zeros in two matrices are removed

    Requires two sparse matrices as input
    """
    assert mat1.shape == mat2.shape, "Matrices have different shapes. "\
        "Computation of correlation is not possible."

    # create a new matrix that is the sum of the two
    # matrices to compare. The goal is to have
    # a matrix that contains all the positions
    # that are non-zero in both matrices
    _mat = mat1 + mat2

    # add one to each element in the new matrix
    _mat.data += 1

    # get a vector of the values in mat1 from
    # _mat
    values1 = (_mat - mat1).data - 1

    # get a vector of the values in mat2 from
    # _mat
    values2 = (_mat - mat2).data - 1

    return values1, values2


def main(args=None):
    """Run the C++ implementation of hicCorrelate with `args` (default sys.argv[1:])."""
    from hicexplorer._cpp import entry_point
    return entry_point('hicCorrelate')(args)
