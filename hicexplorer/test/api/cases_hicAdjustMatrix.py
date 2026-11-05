"""Equivalence cases for hicexplorer.hicAdjustMatrix (apieq.py).

adjustMatrix gets the argparse namespace the tool builds, on the matrices
and region files of the hicAdjustMatrix tests (small_test_matrix_50kb_res
as h5 and cool, gm12878_1_2_3.cool with weights). It returns a hiCMatrix,
which has no value representation; the drivers return its state (matrix,
bin table, nan_bins, correction factors, the masked bin bookkeeping,
chrBinBoundaries, interval trees and the file handler's state that save()
reads) and, for 'saved', the state of the matrix written by save() and
loaded again. exit(1) is compared as the SystemExit code.

Deviation (hic_object.hpp): the cool metadata of the file handler holds
str() values, so only the three entries save() reads are compared.
"""
import os

from apieq import Case, File

HERE = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = os.path.join(HERE, '..', 'test_data')
A = 'hicexplorer.hicAdjustMatrix.adjustMatrix'
H5 = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.h5')
COOL = os.path.join(TEST_DATA, 'small_test_matrix_50kb_res.cool')
GM = os.path.join(TEST_DATA, 'hicAdjustMatrix', 'gm12878_1_2_3.cool')
REGIONS = os.path.join(TEST_DATA, 'regions.bed')
REGIONS_XFAIL = os.path.join(TEST_DATA, 'regions_xfail.bed')


def namespace(matrix, chromosomes=None, regions=None, action='keep', maskBadRegions=None, interIntraHandling=None):
    def make():
        import argparse
        return (argparse.Namespace(matrix=matrix, chromosomes=chromosomes, regions=regions, action=action,
                                   maskBadRegions=maskBadRegions, interIntraHandling=interIntraHandling,
                                   outFileName='unused.cool'),)
    return make


def state(hic):
    if hic is None:
        return None
    trees = getattr(hic, 'interval_trees', None)
    boundaries = getattr(hic, 'chrBinBoundaries', None)
    handler = hic.matrixFileHandler
    handler_state = None
    if handler is not None:
        f = handler.matrixFile
        handler_state = (type(f).__name__, getattr(f, 'matrixFileName', None), getattr(f, 'chrnameList', None),
                         getattr(f, 'correctionOperator', None), getattr(f, 'hic2cool_version', None),
                         getattr(f, 'hicmatrix_version', None),
                         [str((getattr(f, 'hic_metadata', None) or {}).get(k)) for k in
                          ('genome-assembly', 'matrix-generated-by', 'matrix-generated-by-url')])
    return {
        'matrix': hic.matrix,
        'cut_intervals': hic.cut_intervals,
        'nan_bins': hic.nan_bins,
        'correction_factors': hic.correction_factors,
        'distance_counts': hic.distance_counts,
        'orig_bin_ids': hic.orig_bin_ids,
        'orig_cut_intervals': hic.orig_cut_intervals,
        'chrBinBoundaries': None if boundaries is None else list(boundaries.items()),
        'interval_trees': None if trees is None else
        [(name, sorted((iv.begin, iv.end, iv.data) for iv in tree)) for name, tree in trees.items()],
        'handler': handler_state,
        'bin_size': hic.bin_size,
    }


def run(target, args):
    try:
        return state(target(args))
    except SystemExit as error:
        return 'SystemExit', error.code


def run_with_regions(target, regions, matrix, action):
    import argparse
    return run(target, argparse.Namespace(matrix=matrix, chromosomes=None, regions=regions, action=action,
                                          maskBadRegions=None, interIntraHandling=None))


def saved(target, args, out):
    from hicmatrix import HiCMatrix as hm
    hic = target(args)
    hic.save(out)
    written = out if os.path.exists(out) else out + '.h5'
    reloaded = hm.hiCMatrix(written)
    result = state(reloaded)
    result.pop('handler')
    return os.path.basename(written), result


R = dict(driver='run', compare_args=False)
CASES = []
for name, path in (('h5', H5), ('cool', COOL)):
    for action in ('keep', 'remove', 'mask'):
        CASES += [
            Case('{}/chromosomes_{}'.format(name, action), A, args=namespace(path, ['chrX', 'chr3R'], action=action), **R),
            Case('{}/one_chromosome_{}'.format(name, action), A, args=namespace(path, ['chr3R'], action=action), **R),
            Case('{}/regions_{}'.format(name, action), A, args=namespace(path, regions=REGIONS, action=action), **R),
            Case('{}/regions_xfail_{}'.format(name, action), A, args=namespace(path, regions=REGIONS_XFAIL,
                                                                               action=action), **R),
            Case('{}/unknown_chromosome_{}'.format(name, action), A, args=namespace(path, ['chr10'], action=action),
                 **R),
            Case('{}/unknown_and_known_{}'.format(name, action), A,
                 args=namespace(path, ['chr10', 'chrX'], action=action), **R),
        ]
CASES += [
    Case('gm/keep_1', A, args=namespace(GM, ['1']), **R),
    Case('gm/keep_2_3', A, args=namespace(GM, ['3', '2']), **R),
    Case('gm/mask_2', A, args=namespace(GM, ['2'], action='mask'), **R),
    Case('gm/inter', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='inter'), **R),
    Case('gm/intra', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='intra'), **R),
    Case('gm/invalid_handling', A, args=namespace(GM, ['1', '2', '3'], interIntraHandling='both'), **R),
    Case('gm/mask_bad_regions', A, args=namespace(GM, maskBadRegions=True, chromosomes=[]), **R),
    Case('gm/mask_bad_regions_no_chromosomes', A, args=namespace(GM, maskBadRegions=True), **R),
    Case('h5/mask_bad_regions', A, args=namespace(H5, maskBadRegions=True), **R),
    Case('h5/inter', A, args=namespace(H5, ['chrX', 'chr2L'], interIntraHandling='inter'), **R),
    Case('cool/keep_one_unknown', A, args=namespace(COOL, ['chrZ']), **R),
    Case('both_given', A, args=namespace(H5, ['chrX'], regions=REGIONS), **R),
    Case('nothing_given', A, args=namespace(H5), **R),
    Case('nothing_given_handling', A, args=namespace(H5, interIntraHandling='inter'), **R),
    Case('short_region_lines', A, args=(File('r.bed'), H5, 'remove'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\t10000000\t12000000\n\nab\nchr3R\t0\t1000000\tx\n'}),
    Case('region_lines_mask', A, args=(File('r.bed'), COOL, 'mask'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\t10000000\t12000000\nchr2L\t500000\t900000\n'}),
    Case('bad_region_int', A, args=(File('r.bed'), H5, 'keep'), driver='run_with_regions', compare_args=False,
         files={'r.bed': b'chrX\tten\t12000000\n'}),
    Case('missing_region_file', A, args=(File('missing.bed'), H5, 'keep'), driver='run_with_regions',
         compare_args=False),
    Case('saved/h5_to_cool_name', A, args=lambda: namespace(H5, ['chrX'])() + (File('out.cool'),),
         driver='saved', compare_args=False),
    Case('saved/cool_to_cool', A, args=lambda: namespace(COOL, ['chrX', 'chr2L'], action='remove')() +
         (File('out.cool'),), driver='saved', compare_args=False),
    Case('saved/gm_weights', A, args=lambda: namespace(GM, ['2'])() + (File('out.cool'),), driver='saved',
         compare_args=False),
    Case('saved/gm_to_h5', A, args=lambda: namespace(GM, ['1', '3'])() + (File('out.h5'),), driver='saved',
         compare_args=False),
]
