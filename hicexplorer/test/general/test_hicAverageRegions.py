
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import os
from tempfile import NamedTemporaryFile
from hicexplorer import hicAverageRegions
import numpy.testing as nt

from scipy.sparse import load_npz
from hicexplorer.test.test_compute_function import compute

import logging
log = logging.getLogger(__name__)

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")


def test_average_regions_start():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions {} -o {} --range 100000 100000 -cb {}".format(matrix, bed_file, outfile.name, 'start').split()
    log.debug('path: {}'.format(matrix))

    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)
    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_start.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_start():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions  {} -o {} --rangeInBins 100 100 -cb {}".format(matrix, bed_file, outfile.name, 'start').split()
    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_range_in_bins_start.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_center():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions {} -o {} --range 100000 100000 -cb {}".format(matrix, bed_file, outfile.name, 'center').split()
    log.debug('path: {}'.format(matrix))

    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_center.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_center():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions  {} -o {} --rangeInBins 100 100 -cb {}".format(matrix, bed_file, outfile.name, 'center').split()
    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_range_in_bins_center.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_end():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions {} -o {} --range 100000 100000 -cb {}".format(matrix, bed_file, outfile.name, 'end').split()
    log.debug('path: {}'.format(matrix))

    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_end.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_end():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions_multi.bed'
    args = "--matrix {} --regions  {} -o {} --rangeInBins 100 100 -cb {}".format(matrix, bed_file, outfile.name, 'end').split()
    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/regions_multi_range_in_bins_end.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_single():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions.bed'
    args = "--matrix {} --regions {} -o {} --range 100000 100000".format(matrix, bed_file, outfile.name).split()
    log.debug('path: {}'.format(matrix))

    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/result_range_100000.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)


def test_average_regions_range_in_bins_single():

    outfile = NamedTemporaryFile(suffix='.npz', prefix='average_region', delete=False)
    matrix = ROOT + 'small_test_matrix.cool'
    bed_file = ROOT + 'hicAverageRegions/regions.bed'
    args = "--matrix {} --regions  {} -o {} --rangeInBins 100 100".format(matrix, bed_file, outfile.name).split()
    # hicAverageRegions.main(args)
    compute(hicAverageRegions.main, args, 5)

    test_file = load_npz(ROOT + 'hicAverageRegions/result_rangeInBins_100.npz')
    new_file = load_npz(outfile.name)

    nt.assert_almost_equal(test_file.data, new_file.data, decimal=0)

    os.remove(outfile.name)
