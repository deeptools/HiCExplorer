"""Equivalence cases for hicexplorer.hicPlotDistVsCounts (apieq.py).

Both functions take a hiCMatrix, which the drivers load in each interpreter
from the repository's matrices: small_test_matrix_50kb_res (h5 int64 and cool
int32, 50 kb bins), Li_et_al_2015.h5 (float64, restriction fragment bins, so
fit_cut_intervals snaps the starts) and domains.bed for the --domains path.
compute_distance_mean replaces hicmat.matrix, so the drivers return the
matrix after the call and the cached bin size too. exit() raises SystemExit,
which apieq does not catch, so the drivers return it as a value.

Deviations from 3.7.6 (cpp/python/src/api/hicPlotDistVsCounts.cpp):
duplicate matrix entries are summed in storage order; a region overlapping
several bins that were overlapped before names the first in bin order; fh
must yield str.
"""
import os

import numpy as np

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicPlotDistVsCounts.'
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
LI = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')
DOMAINS = os.path.join(TEST_DATA, 'domains.bed')


def load(path):
    from hicmatrix import HiCMatrix as hm
    return hm.hiCMatrix(path)


def matrix_state(hic):
    m = hic.matrix
    return (type(m).__name__, m.format, m.shape, m.dtype.str, m.data, m.indices, m.indptr)


def distance_mean(target, path, *args, **kwargs):
    hic = load(path)
    keep = kwargs.pop('keep', None)
    if keep is not None:
        hic.keepOnlyTheseChr(keep)
    inject_nan = kwargs.pop('inject_nan', False)
    if inject_nan:
        hic.matrix.data = hic.matrix.data.astype(np.float64)
        # far from the diagonal, in both triangles, and one infinite value
        rows, cols = hic.matrix.nonzero()
        far = np.flatnonzero(cols - rows > 100)[:3]
        hic.matrix.data[far[0]] = np.nan
        hic.matrix.data[far[1]] = np.inf
        near = np.flatnonzero(cols - rows == 2)[:2]
        hic.matrix.data[near[0]] = np.nan
        hic.matrix.data[near[1]] = 0
    try:
        result = target(hic, *args, **kwargs)
    except SystemExit as error:
        return ('SystemExit', error.code), hic.bin_size
    return result, matrix_state(hic), hic.bin_size


def with_domains(target, path, domains, *args):
    import importlib
    module = importlib.import_module('hicexplorer.hicPlotDistVsCounts')
    hic = load(path)
    with open(domains) as fh:
        custom = module.from_bed_to_cut_interval(hic, fh)
    result = target(hic, *args, custom_cut_intervals=custom)
    return result, matrix_state(hic)


def wrong_custom(target):
    hic = load(H5)
    return target(hic, None, False, hic.cut_intervals[:-1])


def bed_file(target, path, bed):
    hic = load(path)
    with open(bed) as fh:
        return target(hic, fh)


def bed_lines(target, path, lines):
    return target(load(path), lines)


BED_SMALL = b"# comment\nchr2L\t0\t120000\nchr2L\t150000\t260000\tname\n2R\t100\t200\nchrZ\t1\t2\n" \
            b"chrX\t500000\t1000000\n"
with open(DOMAINS, 'rb') as _handle:
    # every second domain, so that no bin is overlapped by two of them
    DOMAINS_APART = b''.join(_handle.readlines()[::2])
BED_TWICE = b"chr2L\t0\t120000\nchr2L\t100000\t300000\n"


CASES = [
    # compute_distance_mean
    Case('mean/h5', M + 'compute_distance_mean', driver='distance_mean', args=(H5,)),
    Case('mean/h5_maxdepth', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 3000000)),
    Case('mean/h5_perchr', M + 'compute_distance_mean', driver='distance_mean', args=(H5, None, True)),
    Case('mean/h5_perchr_maxdepth', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 1000000, True)),
    Case('mean/h5_float_maxdepth', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 777777.5)),
    Case('mean/h5_maxdepth_zero', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 0)),
    Case('mean/h5_maxdepth_below_binsize', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 1000)),
    Case('mean/cool_int32', M + 'compute_distance_mean', driver='distance_mean', args=(COOL, 2000000, True)),
    Case('mean/cool_all', M + 'compute_distance_mean', driver='distance_mean', args=(COOL,)),
    Case('mean/li_fragments', M + 'compute_distance_mean', driver='distance_mean', args=(LI,)),
    Case('mean/li_maxdepth', M + 'compute_distance_mean', driver='distance_mean', args=(LI, 100000)),
    Case('mean/li_perchr', M + 'compute_distance_mean', driver='distance_mean', args=(LI, 50000, True)),
    Case('mean/some_chromosomes', M + 'compute_distance_mean', driver='distance_mean',
         args=(H5, None, True), kwargs={'keep': ['chr2L', 'chr4', 'chrX']}),
    Case('mean/nan_and_inf', M + 'compute_distance_mean', driver='distance_mean', args=(H5, 500000),
         kwargs={'inject_nan': True}),
    Case('mean/nan_and_inf_no_depth', M + 'compute_distance_mean', driver='distance_mean', args=(H5,),
         kwargs={'inject_nan': True}),
    Case('mean/keywords', M + 'compute_distance_mean', driver='distance_mean', args=(H5,),
         kwargs={'maxdepth': 2000000, 'perchr': True}),
    Case('mean/domains', M + 'compute_distance_mean', driver='with_domains', args=(LI, DOMAINS)),
    Case('mean/domains_perchr', M + 'compute_distance_mean', driver='with_domains', args=(LI, DOMAINS, 200000, True)),
    Case('mean/domains_small', M + 'compute_distance_mean', driver='with_domains', args=(H5, DOMAINS, None, True)),
    Case('mean/domains_apart', M + 'compute_distance_mean', driver='with_domains', args=(LI, File('d.bed')),
         files={'d.bed': DOMAINS_APART}),
    Case('mean/domains_apart_perchr', M + 'compute_distance_mean', driver='with_domains',
         args=(LI, File('d.bed'), 300000, True), files={'d.bed': DOMAINS_APART}),
    Case('mean/domains_apart_small', M + 'compute_distance_mean', driver='with_domains',
         args=(H5, File('d.bed'), None, True), files={'d.bed': DOMAINS_APART}),
    Case('mean/wrong_custom', M + 'compute_distance_mean', driver='wrong_custom'),

    # from_bed_to_cut_interval
    Case('bed/li_domains', M + 'from_bed_to_cut_interval', driver='bed_file', args=(LI, DOMAINS)),
    Case('bed/small_domains', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, DOMAINS)),
    Case('bed/cool_domains', M + 'from_bed_to_cut_interval', driver='bed_file', args=(COOL, DOMAINS)),
    Case('bed/li_domains_apart', M + 'from_bed_to_cut_interval', driver='bed_file', args=(LI, File('d.bed')),
         files={'d.bed': DOMAINS_APART}),
    Case('bed/small', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': BED_SMALL}),
    Case('bed/small_cool', M + 'from_bed_to_cut_interval', driver='bed_file', args=(COOL, File('a.bed')),
         files={'a.bed': BED_SMALL}),
    Case('bed/overlapped_twice', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': BED_TWICE}),
    Case('bed/no_overlap', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chrZ\t1\t2\n# only\n"}),
    Case('bed/empty_region', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chr2L\t200\t100\n"}),
    Case('bed/outside_bins', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chr2L\t900000000\t900000100\n"}),
    Case('bed/blank_line', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chr2L\t0\t100\n\nchr2L\t200000\t300000\n"}),
    Case('bed/two_fields', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chr2L\t0\n"}),
    Case('bed/not_a_number', M + 'from_bed_to_cut_interval', driver='bed_file', args=(H5, File('a.bed')),
         files={'a.bed': b"chr2L\t0\t1e5\n"}),
    Case('bed/list_of_lines', M + 'from_bed_to_cut_interval', driver='bed_lines',
         args=(COOL, ['chrX 0 100000', '  chr3R\t1000000\t1100000  \n', 'chr3R 2000000 2000001'])),
    Case('bed/empty_string', M + 'from_bed_to_cut_interval', driver='bed_lines', args=(COOL, ['chrX 0 100000', ''])),
]
