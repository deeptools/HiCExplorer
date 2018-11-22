import os
from tempfile import NamedTemporaryFile
from hicexplorer import hicPlotAverageRegions
import numpy.testing as nt
from matplotlib.testing.compare import compare_images

import numpy as np
import pytest
from scipy.sparse import load_npz
import logging
log = logging.getLogger(__name__)

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_average_regions():

    outfile = NamedTemporaryFile(suffix='.png', prefix='average_region', delete=False)
    matrix = ROOT + 'hicAverageRegions/result_range_100000.npz'
    args = "--matrix {} -o {}".format(matrix, outfile.name).split()

    hicPlotAverageRegions.main(args)

    res = compare_images(ROOT + '/hicPlotAverageRegions/defaults.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)


def test_average_regions_plot_log():

    outfile = NamedTemporaryFile(suffix='.png', prefix='average_region_log', delete=False)
    matrix = ROOT + 'hicAverageRegions/result_range_100000.npz'
    args = "--matrix {} -o {} --log".format(matrix, outfile.name).split()
    hicPlotAverageRegions.main(args)

    res = compare_images(ROOT + '/hicPlotAverageRegions/defaults_log.png', outfile.name, tol=40)
    assert res is None, res
    os.remove(outfile.name)


def test_average_regions_plot_log1p():

    outfile = NamedTemporaryFile(suffix='.png', prefix='average_region_log1p', delete=False)
    matrix = ROOT + 'hicAverageRegions/result_range_100000.npz'
    args = "--matrix {} -o {} --log1p".format(matrix, outfile.name).split()
    hicPlotAverageRegions.main(args)

    res = compare_images(ROOT + '/hicPlotAverageRegions/defaults_log1p.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)


def test_average_regions_plot_log_vmin_vmax_colormap():

    outfile = NamedTemporaryFile(suffix='.png', prefix='average_region_log1p', delete=False)
    matrix = ROOT + 'hicAverageRegions/result_range_100000.npz'
    args = "--matrix {} -o {} --vMax 20 --vMin 10 --colorMap plasma".format(matrix, outfile.name).split()
    hicPlotAverageRegions.main(args)

    res = compare_images(ROOT + '/hicPlotAverageRegions/defaults_log_vmin_vmax.png', outfile.name, tol=40)
    assert res is None, res

    os.remove(outfile.name)
