"""Equivalence cases for hicexplorer.hicMergeTADbins (apieq.py).

The merge cases run main() of 3.7.6 step by step on the Li et al. 2015
matrices of the test data (cool and h5): load the matrix with hicmatrix (the
input), restore the masked bins, read the domain boundaries with
get_boundary_bin_id and merge them with merge_tad_bins. The state of the
hiCMatrix object afterwards and the matrix of the written file (read back
with hicmatrix) are compared. The boundary cases call get_boundary_bin_id on
real and constructed BED files, including every exit of 3.7.6.

Deviations from 3.7.6, kept on purpose (cpp/python/src/api/hicMergeTADbins.cpp):

* merge_tad_bins writes the file with the C++ cool and h5 writers of the
  tools, so the file holds the same matrix, bins and attributes as the one
  hicmatrix writes but is not byte identical (creation date, HDF5 chunking);
  the cases compare what hicmatrix reads back.
* hiCMatrix matrices are compared with sorted column indices.
* The debug records of reduce_matrix (elapsed times, hicexplorer.reduceMatrix)
  are not compared; the merge_tad_bins of HiCExplorer 4 does not log them.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
M = 'hicexplorer.hicMergeTADbins.'
DOMAINS = os.path.join(TEST_DATA, 'domains.bed')
LI_H5 = os.path.join(TEST_DATA, 'Li_et_al_2015.h5')
LI_COOL = os.path.join(TEST_DATA, 'Li_et_al_2015.cool')
LI_CUT = os.path.join(TEST_DATA, 'Li_cut.h5')


def sorted_matrix(matrix):
    copy = matrix.copy()
    copy.sort_indices()
    return copy


def hic_state(h):
    return {
        'matrix': sorted_matrix(h.matrix),
        'cut_intervals': list(h.cut_intervals),
        'nan_bins': h.nan_bins,
        'correction_factors': h.correction_factors,
        'distance_counts': h.distance_counts,
        'orig_bin_ids': h.orig_bin_ids,
        'chrBinBoundaries': h.chrBinBoundaries,
        'interval_trees': [(k, sorted((iv.begin, iv.end, iv.data) for iv in v)) for k, v in h.interval_trees.items()],
    }


def exits(function):
    def driver(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except SystemExit as error:
            return ('SystemExit', error.code)
    driver.__name__ = function.__name__
    return driver


@exits
def boundaries(target, matrix, bed, mode='r'):
    from hicmatrix import HiCMatrix as hm
    hic = hm.hiCMatrix(matrix)
    with open(bed, mode) as handle:
        return target(hic, handle)


@exits
def merge(target, directory, matrix, bed, out_name, mask=None, ids=None):
    from hicmatrix import HiCMatrix as hm
    import hicexplorer.hicMergeTADbins as module
    hic = hm.hiCMatrix(matrix)
    hic.restoreMaskedBins()
    if ids is None:
        with open(bed) as handle:
            ids = module.get_boundary_bin_id(hic, handle)
    if mask is not None:
        hic.maskBins(mask)
    # reduce_matrix logs elapsed times at the debug level; not comparable
    import logging
    timing = logging.getLogger('hicexplorer.reduceMatrix')
    level = timing.level
    timing.setLevel(logging.INFO)
    try:
        result = target(hic, ids, os.path.join(directory, out_name))
    finally:
        timing.setLevel(level)
    written = sorted(n for n in os.listdir(directory) if n != os.path.basename(bed))
    saved = None
    for name in written:
        if name.endswith('.h5') or name.endswith('.cool'):
            saved = hic_state(hm.hiCMatrix(os.path.join(directory, name)))
    return {'result': result, 'hic': hic_state(hic), 'written': written, 'saved': saved}


BED = (b"browser position X:1-100\ntrack name=domains\n# comment\n"
       b"X\t100000\t500000\tID_1\t0.5\t.\nX\t500000\t1200000\tID_2\t0.1\t.\nX\t1200000\t1500000\n")

CASES = [
    Case('boundaries_h5', M + 'get_boundary_bin_id', args=(LI_H5, DOMAINS), driver='boundaries'),
    Case('boundaries_cool_binary', M + 'get_boundary_bin_id', args=(LI_COOL, DOMAINS, 'rb'), driver='boundaries'),
    Case('boundaries_comments', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')), files={'d.bed': BED},
         driver='boundaries'),
    Case('boundaries_empty', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')), files={'d.bed': b'# x\n'},
         driver='boundaries'),
    Case('boundaries_short_line', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'X\t100000\t500000\nX 1 2\n'}, driver='boundaries'),
    Case('boundaries_not_an_integer', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'X\t100000\t500000\nX\t1x\t2\n'}, driver='boundaries'),
    Case('boundaries_python_int_syntax', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'X\t 1_000\t+20000\n'}, driver='boundaries'),
    Case('boundaries_start_after_end', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'X\t600000\t500000\n'}, driver='boundaries'),
    Case('boundaries_outside', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'X\t100000\t900000000\n'}, driver='boundaries'),
    Case('boundaries_unknown_chromosome', M + 'get_boundary_bin_id', args=(LI_H5, File('d.bed')),
         files={'d.bed': b'chrZ\t100000\t500000\n'}, driver='boundaries'),
    Case('merge_h5', M + 'merge_tad_bins', args=(File(), LI_H5, DOMAINS, 'out.h5'), driver='merge'),
    Case('merge_h5_named_cool', M + 'merge_tad_bins', args=(File(), LI_CUT, DOMAINS, 'out.cool'), driver='merge'),
    Case('merge_cool', M + 'merge_tad_bins', args=(File(), LI_COOL, DOMAINS, 'out.cool'), driver='merge'),
    Case('merge_cool_named_h5', M + 'merge_tad_bins', args=(File(), LI_COOL, DOMAINS, 'out.h5'), driver='merge'),
    Case('merge_masked', M + 'merge_tad_bins', args=(File(), LI_CUT, DOMAINS, 'out.h5'),
         kwargs={'mask': [3, 4, 50]}, driver='merge'),
    Case('merge_nothing', M + 'merge_tad_bins', args=(File(), LI_CUT, DOMAINS, 'out.h5'), kwargs={'ids': []},
         driver='merge'),
    Case('merge_other_name', M + 'merge_tad_bins', args=(File(), LI_CUT, DOMAINS, 'out.txt'), driver='merge'),
]
