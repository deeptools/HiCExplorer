import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicTransform
from hicmatrix import HiCMatrix as hm
import numpy.testing as nt

from tempfile import NamedTemporaryFile
import os
from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
original_matrix = ROOT + "small_test_matrix_50kb_res.h5"
original_matrix_cool = ROOT + "small_test_matrix.cool"

DELTA_DECIMAL = 0


def test_hic_transfer_obs_exp():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp".format(original_matrix_cool, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)
    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_perChromosome():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp --perChromosome".format(original_matrix_cool, outfile.name).split()
    hicTransform.main(args)
    # compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_per_chromosome.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero".format(original_matrix_cool, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_non_zero.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero_perChromosome():

    outfile = NamedTemporaryFile(suffix='obs_exp_.cool', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero --perChromosome".format(original_matrix_cool, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_non_zero_per_chromosome.cool")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_lieberman():
    outfile = NamedTemporaryFile(suffix='obs_exp_lieberman_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_lieberman".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_lieberman.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero_with_ligation_factor():
    outfile = NamedTemporaryFile(suffix='obs_exp_norm_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero --ligation_factor".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_norm.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data,
                                 new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_obs_exp_non_zero_with_ligation_factor_perChromosome():
    outfile = NamedTemporaryFile(suffix='obs_exp_norm_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method obs_exp_non_zero --ligation_factor --perChromosome".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/obs_exp_norm_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    os.unlink(outfile.name)


def test_hic_transfer_pearson():
    outfile = NamedTemporaryFile(suffix='pearson_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method pearson".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/pearson.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    assert_reference_values_tight(new.matrix, test.matrix, 'pearson.h5')
    os.unlink(outfile.name)


def test_hic_transfer_pearson_perChromosome():
    outfile = NamedTemporaryFile(suffix='pearson_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method pearson --perChromosome".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)

    test = hm.hiCMatrix(ROOT + "hicTransform/pearson_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    assert_reference_values_tight(new.matrix, test.matrix, 'pearson_perChromosome.h5')
    os.unlink(outfile.name)


def test_hic_transfer_covariance():
    outfile = NamedTemporaryFile(suffix='covariance_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method covariance".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    assert_reference_values_tight(new.matrix, test.matrix, 'covariance.h5')
    os.unlink(outfile.name)


def test_hic_transfer_covariance_perChromosome():
    outfile = NamedTemporaryFile(suffix='covariance_.h5', delete=False)
    outfile.close()

    args = "--matrix {} --outFileName {} --method covariance --perChromosome".format(original_matrix, outfile.name).split()
    # hicTransform.main(args)
    compute(hicTransform.main, args, 5)
    test = hm.hiCMatrix(ROOT + "hicTransform/covariance_perChromosome.h5")

    new = hm.hiCMatrix(outfile.name)
    nt.assert_array_almost_equal(test.matrix.data, new.matrix.data, decimal=DELTA_DECIMAL)
    assert_structure_matches_reference(new.matrix, test.matrix)
    assert_reference_values_tight(new.matrix, test.matrix, 'covariance_perChromosome.h5')
    os.unlink(outfile.name)


# ---------------------------------------------------------------------------
# Full precision characterization tests, added 2026-09-02 for the C++ port.
#
# Every test above this line asserts with assert_array_almost_equal at
# decimal=0, that is agreement to the nearest integer, and compares only
# matrix.data. Two consequences, both measured on this tree:
#
#  1. No numeric regression short of a whole unit is detectable. On
#     obs_exp_norm.h5 the checked-in reference holds int64 values while this
#     tree produces float32 ones, differing by up to 1.0 absolute and 9.9e-01
#     relative, and the test passes. The obs/exp references for the .cool
#     input are stale in the same way: obs_exp.cool differs by exactly 1 on
#     some entries because utilities.obs_exp_matrix casts its float result
#     back to the input dtype (int32 here) by truncation, so a last-bit
#     difference in the quotient moves the stored integer by one.
#  2. The sparsity pattern is never asserted. matrix.data is compared without
#     matrix.indices or matrix.indptr, so a transform that put the right
#     values in the wrong places would pass.
#
# The reference files therefore cannot be re-asserted at full precision: they
# were produced by an older environment. What is pinned below instead is the
# pipeline itself, recomputed from numpy and from hicexplorer.utilities inside
# the test, at float64 precision and with the sparsity pattern compared
# exactly. That is a specification of the tool rather than a snapshot of one
# old run of it, and it is what the C++ port is written against.
#
# _reference_pearson_from_first_principles adds an independent check that does
# not go through np.corrcoef at all, so the pearson tests are not merely
# asserting that np.corrcoef equals np.corrcoef.
#
# --chromosomes had no test at all. It has three distinct code paths
# (hicTransform.py:151-156): a cool input with exactly one chromosome goes
# through the cool loader's pChrnameList, a cool input with more than one and
# every h5 input go through keepOnlyTheseChr, and no --chromosomes at all
# keeps the matrix as loaded. All three are covered below.

import numpy as np
import pytest
from scipy.sparse import csr_matrix, lil_matrix, triu
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.utilities import obs_exp_matrix, obs_exp_matrix_non_zero
from hicexplorer.utilities import obs_exp_matrix_lieberman

original_matrix_h5 = original_matrix


def _run(args):
    """hicTransform.main on a split argument string."""
    hicTransform.main(args.split())


def _clean(pDense):
    """convertNansToZeros(csr_matrix(x)) then convertInfsToZeros(x).

    The three-line tail every branch of hicTransform.py shares
    (:98-99, :109-110, :121-122, :134-135). csr_matrix() of a dense block
    stores every entry that is not exactly zero, NaN and infinity included,
    and the two converters then overwrite those with zero without removing
    them, so the result carries explicit zeros. Both file writers call
    eliminate_zeros(), so they do not survive to disk; the pattern before that
    point still matters, because it is what the LIL accumulator of the
    --perChromosome branch is built from.
    """
    matrix = csr_matrix(pDense)
    matrix = convertNansToZeros(matrix)
    matrix = convertInfsToZeros(matrix)
    return matrix


def _guarded(pSubmatrix, pKernel):
    """The `if len(pSubmatrix.data) == 0: return pSubmatrix` guard that opens
    each of hicTransform's four private helpers (:95, :106, :118, :131).

    It is not cosmetic. utilities.obs_exp_matrix returns None for a block with
    no stored entry, because utilities.expected_interactions returns None when
    the distance array is empty, and small_test_matrix.cool has two such
    blocks (chrYHet and chrM). Without the guard the --perChromosome branch
    would crash on that input rather than pass the empty block through.
    """
    if len(pSubmatrix.data) == 0:
        return pSubmatrix
    return _clean(pKernel(pSubmatrix))


def _saved_and_reloaded(pMatrix):
    """What hicmatrix stores with pSymmetric=True and reads back.

    save() keeps triu(matrix, k=0) and calls eliminate_zeros() on it
    (hicmatrix/lib/h5.py:116-121, cool.py does the same); the loader then
    calls fillLowerTriangle, which is matrix + triu(matrix, 1).T
    (HiCMatrix.py:106-118). Applying both to the expected matrix is what makes
    an exact comparison against a reloaded output file possible.
    """
    saved = triu(pMatrix, k=0, format='csr')
    saved.eliminate_zeros()
    return (saved + triu(saved, 1).T).tocsr()


def assert_matrix_identical(pActual, pExpected, pMessage=''):
    """Sparsity pattern exactly equal and every value bit-identical."""
    assert pActual.shape == pExpected.shape, \
        '{}shape {} != {}'.format(pMessage, pActual.shape, pExpected.shape)
    nt.assert_array_equal(pActual.indptr, pExpected.indptr,
                          err_msg=pMessage + 'indptr differs')
    nt.assert_array_equal(pActual.indices, pExpected.indices,
                          err_msg=pMessage + 'indices differs')
    assert pActual.nnz == pExpected.nnz, \
        '{}nnz {} != {}'.format(pMessage, pActual.nnz, pExpected.nnz)
    nt.assert_array_equal(np.asarray(pActual.data),
                          np.asarray(pExpected.data),
                          err_msg=pMessage + 'data differs')


def assert_structure_matches_reference(pActual, pReference):
    """The sparsity pattern of a checked-in reference, which is not stale.

    Only the *values* of the references drifted; every one of the eleven has
    the same shape, indptr and indices as this tree produces, measured
    2026-09-02. Asserting that here turns the eleven decimal=0 tests above
    from value-only into value-and-pattern checks without touching their
    tolerance.
    """
    assert pActual.shape == pReference.shape
    nt.assert_array_equal(pActual.indptr, pReference.indptr)
    nt.assert_array_equal(pActual.indices, pReference.indices)


def _per_chromosome(pHiCMatrix, pBlock):
    """The LIL accumulation of hicTransform.py's --perChromosome branches.

    pBlock(submatrix) returns the transformed csr block. Note what the LIL
    does and the whole-matrix branch does not: assigning a block into a
    lil_matrix drops every explicitly stored zero (verified: a csr with five
    stored entries one of which is an explicit zero becomes four after
    `lil[a:b, a:b] = block`), so the two branches do not produce the same
    sparsity pattern even where they compute the same numbers. Inter
    chromosomal entries are dropped entirely, because the accumulator starts
    empty.
    """
    result = lil_matrix(pHiCMatrix.matrix.shape)
    for chrname in pHiCMatrix.getChrNames():
        chr_range = pHiCMatrix.getChrBinRange(chrname)
        submatrix = pHiCMatrix.matrix[chr_range[0]:chr_range[1],
                                      chr_range[0]:chr_range[1]]
        block = pBlock(submatrix)
        if len(block.data) == 0:
            block = lil_matrix(block.shape)
        else:
            block = lil_matrix(block)
        result[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = block
    return result.tocsr()


# --- pearson ---------------------------------------------------------------


def assert_reference_values_tight(pActual, pReference, pName):
    """The checked-in reference at an E3-style bound instead of decimal=0.

    Only for the four correlation references. Measured on this tree,
    2026-09-02, against test_data/hicTransform/:

        pearson.h5                  max abs 2.02e-14 over 7,038,409 values
        pearson_perChromosome.h5    max abs 1.57e-14 over 1,142,710 values
        covariance.h5               max abs 2.44e-14 over 7,038,409 values
        covariance_perChromosome.h5 max abs 1.42e-14 over 1,142,710 values

    so these four files are *not* stale: they still describe what this tree
    computes, to fourteen decimal places. The obs/exp references are, by up to
    1.0 absolute, because utilities.obs_exp_matrix truncates its result back to
    the input integer dtype and a last-bit difference in the quotient moves the
    stored integer by one; those four keep decimal=0 and are covered instead by
    the full precision tests further down.

    The bound is absolute rather than relative because a correlation matrix is
    full of entries whose true value is zero and which both implementations
    produce as cancellation noise around 1e-17. A pure relative comparison
    fails on those however small they are, which is also why the C++ harness
    declares the per chromosome cases E3 rather than ED.
    """
    nt.assert_allclose(np.asarray(pActual.data, dtype=np.float64),
                       np.asarray(pReference.data, dtype=np.float64),
                       rtol=1e-12, atol=1e-12,
                       err_msg=pName + ': values drifted from the reference')


def test_pearson_matches_numpy_corrcoef_at_full_precision():
    """The whole-matrix branch, hicTransform.py:229.

    np.corrcoef is taken of the *whole genome* densified matrix, not of a
    chromosome block, although the --method help text says the transformation
    is computed per chromosome. Pinned as it stands.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method pearson'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    expected = _saved_and_reloaded(
        _clean(np.corrcoef(reference.matrix.todense())))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_pearson_perChromosome_matches_numpy_corrcoef_at_full_precision():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method pearson --perChromosome'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    expected = _saved_and_reloaded(_per_chromosome(
        reference, lambda block: _clean(np.corrcoef(block.todense()))))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def _reference_pearson_from_first_principles(pDense):
    """Pearson correlation of the rows, written out rather than called.

    Independent of np.corrcoef so that the two tests above are not asserting
    that numpy equals numpy. The covariance uses ddof=1, which is np.cov's
    default and therefore what np.corrcoef divides through; the normalisation
    cancels it, so the result is the plain centred correlation.
    """
    values = np.asarray(pDense, dtype=np.float64)
    centred = values - values.mean(axis=1, keepdims=True)
    norms = np.sqrt((centred * centred).sum(axis=1))
    with np.errstate(divide='ignore', invalid='ignore'):
        result = (centred @ centred.T) / np.outer(norms, norms)
    return result


def test_pearson_agrees_with_the_textbook_correlation_on_a_real_chromosome():
    """chrXHet of the 50 kb matrix, 5 bins, --perChromosome.

    The smallest real chromosome block in this matrix, small enough that the
    correlation can be written out in full and compared against the tool
    without going through np.corrcoef. 1e-12 relative rather than exact,
    because the two expressions accumulate in different orders.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method pearson --perChromosome'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    first, last = reference.getChrBinRange('chrXHet')
    block = reference.matrix[first:last, first:last].todense()
    expected = _reference_pearson_from_first_principles(block)
    expected[~np.isfinite(expected)] = 0.0

    actual = hm.hiCMatrix(outfile.name).matrix[first:last, first:last].todense()
    nt.assert_allclose(np.asarray(actual), expected, rtol=1e-12, atol=1e-12)
    os.unlink(outfile.name)


# --- covariance ------------------------------------------------------------


def test_covariance_matches_numpy_cov_at_full_precision():
    """hicTransform.py:249.

    Note the asymmetry with the pearson branch: covariance never calls
    convertNansToZeros, so a NaN produced by np.cov would reach the file. It
    cannot on this input, because np.cov performs no division by a per-row
    quantity, but the branch is different code and is pinned as such.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method covariance'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    expected = _saved_and_reloaded(
        csr_matrix(np.cov(reference.matrix.todense())))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_covariance_perChromosome_matches_numpy_cov_at_full_precision():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method covariance --perChromosome'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    expected = _saved_and_reloaded(_per_chromosome(
        reference, lambda block: csr_matrix(np.cov(block.todense()))))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_covariance_of_a_real_chromosome_is_the_centred_second_moment():
    """chrXHet again, np.cov written out: (X - mean) (X - mean)^T / (n - 1)."""
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method covariance --perChromosome'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    first, last = reference.getChrBinRange('chrXHet')
    block = np.asarray(reference.matrix[first:last, first:last].todense(),
                       dtype=np.float64)
    centred = block - block.mean(axis=1, keepdims=True)
    expected = (centred @ centred.T) / (block.shape[1] - 1)

    actual = hm.hiCMatrix(outfile.name).matrix[first:last, first:last].todense()
    nt.assert_allclose(np.asarray(actual), expected, rtol=1e-12, atol=1e-12)
    os.unlink(outfile.name)


# --- obs/exp ---------------------------------------------------------------


def test_obs_exp_matches_the_utilities_kernel_at_full_precision():
    """The whole-matrix branch on the cool input, hicTransform.py:174.

    Also pins the dtype: utilities.obs_exp_matrix casts its float64 quotient
    back to `type(pSubmatrix.data[0])`, the dtype the matrix had before the
    float32 step, so an integer input matrix comes back integer and the
    fractional part of every obs/exp value is discarded by truncation.
    """
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp'.format(
        original_matrix_cool, outfile.name))

    reference = hm.hiCMatrix(original_matrix_cool)
    assert reference.matrix.dtype == np.int32
    expected = _saved_and_reloaded(_clean(obs_exp_matrix(reference.matrix)))

    actual = hm.hiCMatrix(outfile.name).matrix
    assert actual.dtype == np.int32, \
        'the obs/exp result is cast back to the input dtype'
    assert_matrix_identical(actual, expected)
    os.unlink(outfile.name)


def test_obs_exp_perChromosome_matches_the_utilities_kernel_at_full_precision():
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp --perChromosome'.format(
        original_matrix_cool, outfile.name))

    reference = hm.hiCMatrix(original_matrix_cool)
    expected = _saved_and_reloaded(_per_chromosome(
        reference, lambda block: _guarded(block, obs_exp_matrix)))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_obs_exp_non_zero_matches_the_utilities_kernel_at_full_precision():
    """utilities.obs_exp_matrix_non_zero, which is *not* obs_exp_matrix.

    It divides element by element and assigns each quotient back into a
    float32 array, so the output stays float32 rounded value by value, where
    obs_exp_matrix divides the whole array at once, comes back float64 and is
    then cast to the input dtype. The two are different computations.
    """
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero'.format(
        original_matrix_cool, outfile.name))

    reference = hm.hiCMatrix(original_matrix_cool)
    expected = _saved_and_reloaded(
        _clean(obs_exp_matrix_non_zero(reference.matrix, False)))

    actual = hm.hiCMatrix(outfile.name).matrix
    assert actual.dtype == np.float32, \
        'obs_exp_non_zero never casts back, so its output is float32'
    assert_matrix_identical(actual, expected)
    os.unlink(outfile.name)


def test_obs_exp_non_zero_with_ligation_factor_at_full_precision():
    """The Homer scaling branch, on the h5 input."""
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero '
         '--ligation_factor'.format(original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    expected = _saved_and_reloaded(
        _clean(obs_exp_matrix_non_zero(reference.matrix, True)))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_obs_exp_lieberman_matches_the_utilities_kernel_at_full_precision():
    """hicTransform.py:194-211.

    obs_exp_lieberman ignores --perChromosome: it always runs the per
    chromosome loop, and it feeds every block the same two genome-wide
    numbers, the total bin count and the chromosome count.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_lieberman'.format(
        original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    chromosome_count = len(reference.getChrNames())
    length_chromosome = 0
    for chrname in reference.getChrNames():
        chr_range = reference.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]
    assert length_chromosome == reference.matrix.shape[0]

    expected = _saved_and_reloaded(_per_chromosome(
        reference,
        lambda block: _guarded(
            block,
            lambda inner: obs_exp_matrix_lieberman(
                inner, length_chromosome, chromosome_count))))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_obs_exp_lieberman_ignores_perChromosome():
    """--perChromosome is accepted and changes nothing for this method."""
    plain = NamedTemporaryFile(suffix='.h5', delete=False)
    plain.close()
    per_chromosome = NamedTemporaryFile(suffix='.h5', delete=False)
    per_chromosome.close()

    _run('--matrix {} --outFileName {} --method obs_exp_lieberman'.format(
        original_matrix_h5, plain.name))
    _run('--matrix {} --outFileName {} --method obs_exp_lieberman '
         '--perChromosome'.format(original_matrix_h5, per_chromosome.name))

    assert_matrix_identical(hm.hiCMatrix(per_chromosome.name).matrix,
                           hm.hiCMatrix(plain.name).matrix)
    os.unlink(plain.name)
    os.unlink(per_chromosome.name)


# --- --chromosomes ---------------------------------------------------------


def test_chromosomes_h5_keeps_only_those_bins():
    """The keepOnlyTheseChr path, hicTransform.py:156. Never tested before."""
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero '
         '--chromosomes chrX chr4'.format(original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    reference.keepOnlyTheseChr(['chrX', 'chr4'])
    expected = _saved_and_reloaded(
        _clean(obs_exp_matrix_non_zero(reference.matrix, False)))

    actual = hm.hiCMatrix(outfile.name)
    assert actual.getChrNames() == reference.getChrNames()
    assert actual.matrix.shape == reference.matrix.shape
    assert_matrix_identical(actual.matrix, expected)
    os.unlink(outfile.name)


def test_chromosomes_changes_the_pearson_result_not_just_its_size():
    """--chromosomes before a whole-matrix pearson.

    The correlation is taken over the retained columns only, so restricting
    the matrix changes every value, not only the shape. Without this the port
    could apply --chromosomes after the transform and still produce a file of
    the right shape.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method pearson '
         '--chromosomes chrX chr4'.format(original_matrix_h5, outfile.name))

    reference = hm.hiCMatrix(original_matrix_h5)
    reference.keepOnlyTheseChr(['chrX', 'chr4'])
    expected = _saved_and_reloaded(
        _clean(np.corrcoef(reference.matrix.todense())))

    actual = hm.hiCMatrix(outfile.name).matrix
    assert_matrix_identical(actual, expected)

    whole = hm.hiCMatrix(original_matrix_h5)
    first, last = whole.getChrBinRange('chr4')
    restricted_first, restricted_last = \
        hm.hiCMatrix(outfile.name).getChrBinRange('chr4')
    assert last - first == restricted_last - restricted_first
    os.unlink(outfile.name)


def test_chromosomes_cool_single_chromosome_uses_the_loader_path():
    """hicTransform.py:151-152: a cool input plus exactly one chromosome.

    The chromosome is passed to the cool loader as pChrnameList instead of
    being removed afterwards with keepOnlyTheseChr. The two are not obviously
    the same operation, so both the bin table and the values are compared.
    """
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero '
         '--chromosomes chrX'.format(original_matrix_cool, outfile.name))

    reference = hm.hiCMatrix(pMatrixFile=original_matrix_cool,
                             pChrnameList=['chrX'])
    expected = _saved_and_reloaded(
        _clean(obs_exp_matrix_non_zero(reference.matrix, False)))

    actual = hm.hiCMatrix(outfile.name)
    assert actual.getChrNames() == ['chrX']
    assert actual.matrix.shape == reference.matrix.shape
    assert_matrix_identical(actual.matrix, expected)
    os.unlink(outfile.name)


def test_chromosomes_cool_two_chromosomes_uses_keepOnlyTheseChr():
    """Two chromosomes take the other branch even for a cool input."""
    outfile = NamedTemporaryFile(suffix='.cool', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero '
         '--chromosomes chrX chr4'.format(original_matrix_cool, outfile.name))

    reference = hm.hiCMatrix(original_matrix_cool)
    reference.keepOnlyTheseChr(['chrX', 'chr4'])
    expected = _saved_and_reloaded(
        _clean(obs_exp_matrix_non_zero(reference.matrix, False)))

    actual = hm.hiCMatrix(outfile.name)
    assert actual.getChrNames() == reference.getChrNames()
    assert_matrix_identical(actual.matrix, expected)
    os.unlink(outfile.name)


def test_unknown_chromosome_is_an_error():
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    with pytest.raises(Exception):
        _run('--matrix {} --outFileName {} --method obs_exp '
             '--chromosomes notAChromosome'.format(original_matrix_h5,
                                                   outfile.name))
    os.unlink(outfile.name)


# --- output naming ---------------------------------------------------------


def test_output_format_follows_the_input_not_the_file_name():
    """hiCMatrix.save reuses the handler built during the load.

    So `-m x.h5 -o out.cool` writes a HiCExplorer h5 file, and the h5 writer
    appends its own suffix, so the file that appears is `out.cool.h5` and
    `out.cool` is never created. hicTransform's own suffix check
    (hicTransform.py:145) only decides whether the run is refused.
    """
    directory = os.path.dirname(NamedTemporaryFile(delete=True).name)
    target = os.path.join(directory, 'hicTransform_naming_test.cool')
    for stale in (target, target + '.h5'):
        if os.path.exists(stale):
            os.unlink(stale)

    _run('--matrix {} --outFileName {} --method obs_exp_lieberman'.format(
        original_matrix_h5, target))

    assert not os.path.exists(target)
    assert os.path.exists(target + '.h5')
    assert hm.hiCMatrix(target + '.h5').matrix.shape == (2794, 2794)
    os.unlink(target + '.h5')


def test_unknown_output_suffix_exits_one():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    with pytest.raises(SystemExit) as error:
        _run('--matrix {} --outFileName {} --method obs_exp'.format(
            original_matrix_h5, outfile.name))
    assert error.value.code == 1
    os.unlink(outfile.name)


# --- the obs/exp kernels from first principles -----------------------------
#
# The obs/exp tests above call hicexplorer.utilities themselves, so they pin
# which submatrix reaches which kernel, the dtype the result is cast back to
# and the LIL accumulation, but not the kernel's own arithmetic: a mutation
# inside utilities.obs_exp_matrix moves both sides equally. The two tests
# below close that by writing the expected value per genomic distance out in
# the test. Measured: a floor() where utilities has ceil() in the distance
# index is caught by them and by the decimal=0 reference comparison, and by
# nothing else in this file.
#
# Li_cut.h5 is the input rather than the 50 kb matrix because its values are
# float64. utilities.obs_exp_matrix casts its result back to the input dtype,
# so on an integer matrix the fractional part is truncated away and a wrong
# expected value usually survives the truncation unchanged.

original_matrix_float = ROOT + "Li_cut.h5"


def _expected_per_distance_written_out(pMatrix):
    """utilities.expected_interactions(pSubmatrix, pThreads=None), written out.

    Note the divisor: np.arange(n + 1, 1, -1), that is n + 1 - d, which is one
    more than the number of cells diagonal d actually has. Reproduced as it
    stands, not corrected.

    Note also that `distance` is derived from nonzero(), which skips an
    explicitly stored zero, while the values are taken from .data, which does
    not. The two are only the same length when the matrix carries no explicit
    zero. Every matrix reaching this point does, because hicmatrix's loader
    drops them, so the mismatch is latent rather than live.
    """
    size = pMatrix.shape[0]
    row, col = pMatrix.nonzero()
    distance = np.absolute(row - col)
    expected = np.zeros(size)
    for step in range(distance.min(), distance.max() + 1):
        expected[step] = np.sum(pMatrix.data[distance == step])
    occurrences = np.arange(size + 1, 1, -1)
    expected = expected / occurrences
    expected[~np.isfinite(expected)] = 0.0
    return expected


def test_obs_exp_expected_value_per_distance_written_out():
    """obs_exp on a float64 matrix, with the kernel spelled out in the test.

    Two properties of utilities.obs_exp_matrix are pinned here and nowhere
    else: the values are rounded through float32 before the division, and the
    expected value is read at index ceil(d / 2) although it was accumulated at
    index d, so consecutive distances share an expected value.
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp'.format(
        original_matrix_float, outfile.name))

    reference = hm.hiCMatrix(original_matrix_float)
    assert reference.matrix.dtype == np.float64
    per_distance = _expected_per_distance_written_out(reference.matrix)
    row, col = reference.matrix.nonzero()
    index = np.ceil(np.absolute(row - col) / 2).astype(np.int32)
    values = np.divide(reference.matrix.data.astype(np.float32),
                       per_distance[index])
    values[~np.isfinite(values)] = 0.0
    expected = _saved_and_reloaded(csr_matrix(
        (values, reference.matrix.indices.copy(),
         reference.matrix.indptr.copy()), shape=reference.matrix.shape))

    assert_matrix_identical(hm.hiCMatrix(outfile.name).matrix, expected)
    os.unlink(outfile.name)


def test_obs_exp_non_zero_expected_value_per_distance_written_out():
    """obs_exp_non_zero, whose divisor counts only the stored entries.

    The other difference from obs_exp is where the rounding lands: this kernel
    assigns each quotient back into a float32 array one element at a time, so
    every value is rounded to float32 *after* the division, and it reads the
    expected value at index d rather than at ceil(d / 2).
    """
    outfile = NamedTemporaryFile(suffix='.h5', delete=False)
    outfile.close()
    _run('--matrix {} --outFileName {} --method obs_exp_non_zero'.format(
        original_matrix_float, outfile.name))

    reference = hm.hiCMatrix(original_matrix_float)
    size = reference.matrix.shape[0]
    row, col = reference.matrix.nonzero()
    distance = np.absolute(row - col)
    totals = np.zeros(size)
    counts = np.zeros(size)
    for position, step in enumerate(distance):
        totals[step] += reference.matrix.data[position]
        counts[step] += 1
    per_distance = np.divide(totals, counts,
                             out=np.zeros(size), where=counts != 0)
    per_distance[~np.isfinite(per_distance)] = 0.0

    values = reference.matrix.data.astype(np.float32)
    for position in range(len(values)):
        values[position] = np.divide(values[position],
                                     per_distance[distance[position]])
    values[~np.isfinite(values)] = 0.0
    expected = _saved_and_reloaded(csr_matrix(
        (values, reference.matrix.indices.copy(),
         reference.matrix.indptr.copy()), shape=reference.matrix.shape))

    actual = hm.hiCMatrix(outfile.name).matrix
    assert actual.dtype == np.float32
    assert_matrix_identical(actual, expected)
    os.unlink(outfile.name)
