import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
from hicexplorer import hicPCA
from hicmatrix import HiCMatrix as hm

from tempfile import NamedTemporaryFile
import os
import pytest
import numpy.testing as nt
import numpy as np
import pyBigWig
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
from hicexplorer.utilities import obs_exp_matrix_lieberman, obs_exp_matrix_non_zero
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.utilities import enlarge_bins, opener
from hicexplorer.readBed import ReadBed
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")
from hicexplorer.test.test_compute_function import compute

import logging
log = logging.getLogger(__name__)

DELTA_DECIMAL = 0

# Full precision characterization tests, added 2026-09-01 for the C++ port.
#
# The tests above this line assert at decimal=0, agreement to the nearest
# integer, which no realistic numeric regression can violate, and
# are_files_equal_bigwig compares np.absolute() of only the first interval
# before falling through to a whole-array comparison, so half of its
# assertions cannot fail. Everything below asserts at float64 precision and
# covers the option combinations the tests above never reach:
# --method dist_norm without --ligation_factor, --whichEigenvectors with three
# vectors, --histonMarkType inactive, the bigwig writer's payload, and the
# argument-count error path.
#
# Two properties of hicPCA make a naive full-precision assertion impossible,
# and both are pinned here rather than tolerated:
#
#  1. The sign of an eigenvector is arbitrary. hicPCA.py:305 calls
#     scipy.linalg.eig and writes out whatever LAPACK's dgeev produced.
#     Measured on this tree: running hicPCA on small_test_matrix_50kb_res.h5
#     with OPENBLAS_NUM_THREADS=1 and with OPENBLAS_NUM_THREADS=16 gives
#     bedgraph files whose chr2L block is sign-inverted, with every magnitude
#     equal to all twelve printed digits. The sign is therefore compared as a
#     single global flip per chromosome, never per value, except with
#     --extraTrack where the flip is decided by the correlation with the track
#     and is deterministic (verified at 1 and at 16 BLAS threads).
#  2. scipy.linalg.eig is the general non-symmetric solver and does not sort.
#     Its column order coincides with descending eigenvalue order only when
#     the eigenvalues are well separated. On small_test_matrix.h5 chrX under
#     --method lieberman the two leading eigenvalues differ by 8.6e-04
#     relative and dgeev returns them swapped, so hicPCA's "first
#     eigenvector" is the eigenvector of the *second* largest eigenvalue.
#     test_pca_eigenvector_order_is_not_sorted_by_eigenvalue pins that.
#
# The reference values are not taken from the .bedgraph and .bw files checked
# into test_data/hicPCA: those were produced by an older environment and
# differ from what this tree produces by up to 2.5e-02 relative, which is why
# the tests above had to be written at decimal=0 in the first place. The
# expectation is recomputed from the documented pipeline instead
# (_reference_eigenvectors), which pins the algorithm rather than a stale
# artefact, and a small table of literal values pins the arithmetic itself.

# hicPCA writes its bedgraph with '{:.12f}', twelve decimal places, so the
# file quantises every value to 5e-13 absolute. That is the resolution of the
# artefact and therefore the floor of any assertion made against it; below it
# there is nothing to compare. The relative bound stays at 1e-12, so a value of
# order 1 is pinned to all twelve of its printed decimals.
RELATIVE_TOLERANCE = 1e-12
ABSOLUTE_TOLERANCE = 1e-12


def _read_bedgraph(path):
    """[(chrom, start, end, value)] from a hicPCA bedgraph."""
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            rows.append((fields[0], int(fields[1]), int(fields[2]),
                         float(fields[3])))
    return rows


def _blocks_by_chromosome(rows):
    """Row index ranges of each chromosome, in the order they appear."""
    blocks = []
    for index, row in enumerate(rows):
        if blocks and blocks[-1][0] == row[0]:
            blocks[-1][2] = index + 1
        else:
            blocks.append([row[0], index, index + 1])
    return [(name, start, end) for name, start, end in blocks]


def _canonical_sign(values):
    """+1 or -1, chosen so that the largest magnitude entry is positive.

    The eigenvector sign is arbitrary, so a comparison has to fix a
    convention. The largest-magnitude component is the stable choice: it is
    the entry furthest from a sign flip caused by rounding.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 1.0
    extreme = values[int(np.argmax(np.abs(values)))]
    return -1.0 if extreme < 0 else 1.0


def assert_bedgraph_equal(actual_path, expected_rows, pAllowSignFlip=True,
                          pRelTolerance=RELATIVE_TOLERANCE):
    """Coordinates exactly, values at float64 precision.

    With pAllowSignFlip the values of one chromosome may be globally negated,
    which is the only freedom LAPACK leaves. A per-value sign difference is a
    failure, and so is a chromosome whose sign flips halfway through.
    """
    actual = _read_bedgraph(actual_path)
    assert len(actual) == len(expected_rows), \
        'line count {} != {}'.format(len(actual), len(expected_rows))
    for got, want in zip(actual, expected_rows):
        assert got[:3] == tuple(want[:3]), \
            'coordinates differ: {} != {}'.format(got[:3], want[:3])

    got_values = np.array([row[3] for row in actual], dtype=np.float64)
    want_values = np.array([row[3] for row in expected_rows], dtype=np.float64)
    for name, start, end in _blocks_by_chromosome(actual):
        got_block = got_values[start:end]
        want_block = want_values[start:end]
        if pAllowSignFlip:
            got_block = got_block * _canonical_sign(got_block)
            want_block = want_block * _canonical_sign(want_block)
        nt.assert_allclose(got_block, want_block, rtol=pRelTolerance,
                           atol=ABSOLUTE_TOLERANCE,
                           err_msg='chromosome {}'.format(name))


def _covariance_per_chromosome(pMatrixFile, pChromosomes, pMethod,
                               pLigationFactor=False, pIgnoreMaskedBins=False):
    """hicPCA.main's per-chromosome pipeline up to the covariance matrix.

    Reimplemented from hicPCA.py:250-304 rather than imported from it, so that
    a change to the tool's pipeline shows up as a test failure.
    """
    ma = hm.hiCMatrix(pMatrixFile)
    if pIgnoreMaskedBins:
        ma.maskBins(ma.nan_bins)
        ma.setCutIntervals(enlarge_bins(ma.cut_intervals))
    if pChromosomes:
        ma.keepOnlyTheseChr(pChromosomes)

    chromosome_count = len(ma.getChrNames())
    length_chromosome = 0
    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        length_chromosome += chr_range[1] - chr_range[0]

    for chrname in ma.getChrNames():
        chr_range = ma.getChrBinRange(chrname)
        submatrix = ma.matrix[chr_range[0]:chr_range[1],
                              chr_range[0]:chr_range[1]]
        if pMethod == 'lieberman':
            obs_exp = obs_exp_matrix_lieberman(submatrix, length_chromosome,
                                               chromosome_count)
        else:
            obs_exp = obs_exp_matrix_non_zero(submatrix, pLigationFactor)
        obs_exp = csr_matrix(obs_exp).todense()

        pearson = np.corrcoef(obs_exp)
        pearson = convertNansToZeros(csr_matrix(pearson)).todense()
        pearson = convertInfsToZeros(csr_matrix(pearson)).todense()

        covariance = np.cov(obs_exp)
        covariance = convertNansToZeros(csr_matrix(covariance)).todense()
        covariance = convertInfsToZeros(csr_matrix(covariance)).todense()

        intervals = ma.cut_intervals[chr_range[0]:chr_range[1]]
        yield (chrname, intervals, np.asarray(obs_exp), np.asarray(pearson),
               np.asarray(covariance))


def _gene_density_per_bin(pMatrixFile, pChromosomes, pGeneTrack):
    """The gene occurrence array hicPCA.py:148-173 builds from a bed track.

    One count per bin, added at the *end* bin of the interval
    (`gene_occurrence[bin_id[1]] += 1`), which is what the tool does; an
    interval whose chromosome is not in the matrix, whose start runs past the
    chromosome, or whose bin range comes back None is skipped.
    """
    ma = hm.hiCMatrix(pMatrixFile)
    if pChromosomes:
        ma.keepOnlyTheseChr(pChromosomes)
    occurrence = np.zeros(len(ma.cut_intervals))
    chromosomes = ma.getChrNames()
    sizes = ma.get_chromosome_sizes()
    for interval in ReadBed(opener(pGeneTrack)):
        if interval.chromosome not in chromosomes:
            continue
        if interval.start > sizes[interval.chromosome]:
            continue
        bin_id = ma.getRegionBinRange(interval.chromosome, interval.start,
                                      interval.end)
        if bin_id is None:
            continue
        occurrence[bin_id[1]] += 1
    return occurrence


def _reference_eigenvectors(pMatrixFile, pChromosomes, pMethod, pWhich,
                            pLigationFactor=False, pIgnoreMaskedBins=False):
    """The bedgraph rows hicPCA must produce, one list per output file."""
    per_output = [[] for _ in pWhich]
    for chrname, intervals, _obs_exp, _pearson, covariance in \
            _covariance_per_chromosome(pMatrixFile, pChromosomes, pMethod,
                                       pLigationFactor, pIgnoreMaskedBins):
        _evals, eigs = linalg.eig(covariance)
        for output_index, which in enumerate(pWhich):
            column = eigs[:, int(which) - 1].real
            for bin_index, interval in enumerate(intervals):
                per_output[output_index].append(
                    (interval[0], interval[1], interval[2],
                     float(column[bin_index])))
    return per_output


def _run_pca(pArgs):
    hicPCA.main(pArgs)


def _temporary_names(pCount, pSuffix):
    names = []
    for _ in range(pCount):
        handle = NamedTemporaryFile(suffix=pSuffix, delete=False)
        handle.close()
        names.append(handle.name)
    return names


def test_pca_bedgraph_lieberman_full_precision():
    """The default two-eigenvector lieberman run, pinned at float64.

    The magnitudes have to agree to twelve significant digits; the only
    freedom is one global sign per chromosome.
    """
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method lieberman"
            .format(matrix, outputs[0], outputs[1])).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, None, 'lieberman', ['1', '2'])
    for path, rows in zip(outputs, expected):
        assert_bedgraph_equal(path, rows)
    for path in outputs:
        os.unlink(path)


def test_pca_bedgraph_dist_norm_without_ligation_factor():
    """--method dist_norm with no --ligation_factor, which nothing tests.

    Every existing dist_norm test passes --ligation_factor, so the branch at
    utilities.obs_exp_matrix_non_zero that skips the Homer scaling has no
    coverage at all.
    """
    matrix = ROOT + "small_test_matrix.h5"
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method dist_norm --chromosomes chrX chrXHet"
            .format(matrix, outputs[0], outputs[1])).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, ['chrX', 'chrXHet'],
                                       'dist_norm', ['1', '2'],
                                       pLigationFactor=False)
    for path, rows in zip(outputs, expected):
        assert_bedgraph_equal(path, rows)
    for path in outputs:
        os.unlink(path)


def test_pca_bedgraph_dist_norm_with_ligation_factor():
    """--method dist_norm --ligation_factor at full precision."""
    matrix = ROOT + "small_test_matrix.h5"
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method dist_norm --ligation_factor "
            "--chromosomes chrX chrXHet"
            .format(matrix, outputs[0], outputs[1])).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, ['chrX', 'chrXHet'],
                                       'dist_norm', ['1', '2'],
                                       pLigationFactor=True)
    for path, rows in zip(outputs, expected):
        assert_bedgraph_equal(path, rows)
    for path in outputs:
        os.unlink(path)


def test_pca_three_eigenvectors_and_non_contiguous_selection():
    """--whichEigenvectors 1 3 5, which no existing test exercises.

    Every test above asks for '1 2'. The loop at hicPCA.py:314-318 hstacks an
    arbitrary set of columns, and nothing checks that column k-1 is the one
    that ends up in output file k.
    """
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    outputs = _temporary_names(3, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} {} -f bedgraph "
            "--whichEigenvectors 1 3 5 --method lieberman"
            .format(matrix, outputs[0], outputs[1], outputs[2])).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, None, 'lieberman',
                                       ['1', '3', '5'])
    for path, rows in zip(outputs, expected):
        assert_bedgraph_equal(path, rows)
    for path in outputs:
        os.unlink(path)


def test_pca_eigenvector_order_is_not_sorted_by_eigenvalue():
    """hicPCA takes dgeev's column k-1, not the k-th largest eigenvector.

    hicPCA.py:305 calls scipy.linalg.eig, the general non-symmetric solver, on
    a symmetric covariance matrix and slices columns without sorting. For a
    well separated spectrum dgeev happens to emerge in descending order, so
    the defect is invisible; for a near degenerate pair it does not.

    Measured on small_test_matrix.h5 chrX with --method lieberman: the two
    leading eigenvalues are 2795.01342998 and 2792.60820894, a relative gap of
    8.6e-04, and dgeev returns them swapped. hicPCA's first eigenvector is
    therefore the eigenvector of the second largest eigenvalue and its second
    eigenvector is the eigenvector of the largest, and the two are orthogonal,
    so no numeric tolerance can absorb the difference.

    This is the property a reimplementation has to match. A port that sorts by
    descending eigenvalue, which is what a symmetric solver returns naturally,
    produces a different vector here, not a slightly different one.
    """
    blocks = list(_covariance_per_chromosome(ROOT + "small_test_matrix.h5",
                                             ['chrX', 'chrXHet'], 'lieberman'))
    covariance = dict((name, cov) for name, _iv, _oe, _pe, cov in blocks)['chrX']

    eigenvalues, eigenvectors = linalg.eig(covariance)
    eigenvalues = eigenvalues.real
    descending = np.argsort(-eigenvalues, kind='stable')

    # The spectrum really is near degenerate at the top, which is what makes
    # the ordering unstable.
    gap = abs(eigenvalues[descending[0]] - eigenvalues[descending[1]]) \
        / abs(eigenvalues[descending[0]])
    assert gap < 1e-3, 'expected a near degenerate leading pair, gap {}'.format(gap)

    # dgeev's column 0 is the second largest eigenvalue and column 1 is the
    # largest. If a future scipy or LAPACK changes that, this assertion is the
    # place it must be reconsidered, not silently absorbed by a tolerance.
    assert descending[0] == 1, \
        'expected the largest eigenvalue in dgeev column 1, found column {}' \
        .format(descending[0])
    assert descending[1] == 0, \
        'expected the second largest eigenvalue in dgeev column 0, found column {}' \
        .format(descending[1])

    # The two columns span different one dimensional subspaces: this is a
    # selection difference, not a rounding difference.
    first = eigenvectors[:, 0].real
    second = eigenvectors[:, 1].real
    cosine = abs(float(first @ second)) / (np.linalg.norm(first) *
                                           np.linalg.norm(second))
    assert cosine < 1e-8, 'the two leading columns are not orthogonal: {}'.format(cosine)

    # And the symmetric solver, which sorts ascending, agrees on the values
    # while disagreeing on which column they sit in.
    symmetric_values, symmetric_vectors = linalg.eigh(covariance)
    nt.assert_allclose(symmetric_values[-1], eigenvalues[1], rtol=1e-10)
    nt.assert_allclose(symmetric_values[-2], eigenvalues[0], rtol=1e-10)
    largest = symmetric_vectors[:, -1]
    overlap = abs(float(largest @ second)) / (np.linalg.norm(largest) *
                                              np.linalg.norm(second))
    nt.assert_allclose(overlap, 1.0, rtol=1e-8)


def test_pca_bedgraph_values_are_pinned_literals():
    """A table of literal values, so the arithmetic is pinned and not only the shape.

    Recomputing the expectation from the same numpy calls the tool makes would
    not notice a change inside hicexplorer.utilities. These twelve numbers
    would. They are the first three bins of four chromosomes of the default
    lieberman run on small_test_matrix_50kb_res.h5, normalised so that the
    largest magnitude entry of each chromosome is positive, and printed at the
    twelve decimals hicPCA writes. Verified identical at OPENBLAS_NUM_THREADS=1
    and 16 after that normalisation, so they pin the arithmetic and not the
    BLAS threading.
    """
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method lieberman"
            .format(matrix, outputs[0], outputs[1])).split()
    _run_pca(args)

    rows = _read_bedgraph(outputs[0])
    values = np.array([row[3] for row in rows], dtype=np.float64)
    got = {}
    for name, start, end in _blocks_by_chromosome(rows):
        block = values[start:end] * _canonical_sign(values[start:end])
        got[name] = block[:3]

    expected = {
        'chr2RHet': [0.491006273550, 0.001908229902, -0.013181219741],
        'chr3RHet': [-0.064516879370, -0.039350633223, -0.015256798888],
        'chr2L': [0.001937021304, -0.039947959528, 0.041819909727],
        'chr3R': [-0.018165184479, -0.005371327000, -0.005847547014],
    }
    for name, wanted in expected.items():
        nt.assert_allclose(got[name], wanted, rtol=1e-9, atol=1e-12,
                           err_msg='chromosome {}'.format(name))
    for path in outputs:
        os.unlink(path)


def test_pca_extratrack_gene_density_fixes_the_sign_exactly():
    """With --extraTrack the sign is decided, not arbitrary, so compare it.

    correlateEigenvectorWithGeneTrack negates a chromosome's eigenvector when
    its Pearson correlation with the gene density is negative
    (hicPCA.py:190-193). That makes the output sign reproducible where the
    correlation is real: verified byte identical at OPENBLAS_NUM_THREADS=1 and
    16, where the same run without --extraTrack flips a whole chromosome. The
    assertion therefore reproduces the decision and forbids any sign
    difference, rather than allowing a global flip; allowing one would make
    the test blind to the rule being inverted.

    chrX only. On chrXHet of this matrix the eigenvector correlates with the
    gene density at -4.3e-16, which is rounding noise rather than a
    correlation, so --extraTrack does not decide that chromosome's sign
    either. The test asserts that the correlation it relies on is bounded away
    from zero, so a future input that has the same problem fails loudly here
    instead of passing by accident.
    """
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method lieberman --extraTrack {} "
            "--chromosomes chrX"
            .format(matrix, outputs[0], outputs[1], gene_track)).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, ['chrX'], 'lieberman', ['1', '2'])
    density = _gene_density_per_bin(matrix, ['chrX'], gene_track)

    for path, rows in zip(outputs, expected):
        actual = _read_bedgraph(path)
        assert len(actual) == len(rows)
        for name, start, end in _blocks_by_chromosome(actual):
            got = np.array([row[3] for row in actual[start:end]])
            want = np.array([row[3] for row in rows[start:end]])
            # The rule the test asserts, and the reason it is asserted rather
            # than a per chromosome global sign being allowed: hicPCA.py:190-193
            # negates a chromosome whose eigenvector correlates *negatively*
            # with the gene density. Accepting any global sign would make this
            # test blind to the rule being inverted, which is exactly the
            # defect the tests above it have.
            correlation = pearsonr(want, density[start:end])[0]
            if correlation < 0:
                want = -want
            nt.assert_allclose(got, want, rtol=RELATIVE_TOLERANCE,
                               atol=ABSOLUTE_TOLERANCE,
                               err_msg='chromosome {} (correlation with the gene '
                                       'density {})'.format(name, correlation))
            # The premise: the correlation is not near zero on either
            # chromosome, so the decision is not itself a coin toss.
            assert abs(correlation) > 1e-3, \
                'chromosome {} correlates with the gene density at {}, too close ' \
                'to zero for the sign rule to be well determined'.format(
                    name, correlation)
    for path in outputs:
        os.unlink(path)


def test_pca_extratrack_gene_density_is_reproducible():
    """Two runs with --extraTrack produce byte identical bedgraphs."""
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    first = _temporary_names(2, '.bedgraph')
    second = _temporary_names(2, '.bedgraph')
    for outputs in (first, second):
        args = ("--matrix {} --outputFileName {} {} -f bedgraph "
                "--whichEigenvectors 1 2 --method lieberman --extraTrack {} "
                "--chromosomes chrX chrXHet"
                .format(matrix, outputs[0], outputs[1], gene_track)).split()
        _run_pca(args)
    for a, b in zip(first, second):
        with open(a, 'rb') as fa, open(b, 'rb') as fb:
            assert fa.read() == fb.read()
    for path in first + second:
        os.unlink(path)


def test_pca_histone_mark_track_active_and_inactive_are_both_no_ops_here():
    """--histonMarkType inactive, which no test reaches, plus why it is inert.

    The active branch flips a chromosome when the positive compartment's mean
    bigwig coverage is below the negative one's, the inactive branch flips on
    the opposite comparison, and both are guarded by
    `(neg_mean != 0) and (pos_mean != 0)` (hicPCA.py:226-236).

    Measured on the only bigwig in the corpus, bigwig_chrx_2e6_5e6.bw against
    chrX of small_test_matrix.h5: the first eigenvector has 10 positive and
    2,032 negative bins and the second 14 and 2,028, and in both cases every
    positive bin falls outside the track's 2 to 5 Mb coverage, so pos_mean is
    exactly 0 and the guard blocks the flip in *both* modes. The corpus cannot
    tell the two branches apart.

    The test therefore pins that fact rather than pretending otherwise: the two
    modes agree with each other and with the unflipped eigenvectors, value for
    value, and pos_mean is asserted to be zero so that a future change to the
    track or to the guard makes this test's premise fail loudly instead of
    silently passing for the wrong reason.
    """
    matrix = ROOT + "small_test_matrix.h5"
    extra_track = ROOT + 'bigwig_chrx_2e6_5e6.bw'
    active = _temporary_names(2, '.bedgraph')
    inactive = _temporary_names(2, '.bedgraph')

    for outputs, mark_type in ((active, 'active'), (inactive, 'inactive')):
        args = ("--matrix {} --outputFileName {} {} -f bedgraph "
                "--whichEigenvectors 1 2 --method lieberman --extraTrack {} "
                "--chromosomes chrX --histonMarkType {}"
                .format(matrix, outputs[0], outputs[1], extra_track,
                        mark_type)).split()
        _run_pca(args)

    expected = _reference_eigenvectors(matrix, ['chrX'], 'lieberman', ['1', '2'])
    for path in active + inactive:
        rows = _read_bedgraph(path)
        assert len(rows) == len(expected[0])
    for outputs in (active, inactive):
        for path, rows in zip(outputs, expected):
            # No sign freedom: the flip is the only thing --extraTrack could
            # have changed and it did not fire, so the values must be the raw
            # eigenvector components.
            actual = _read_bedgraph(path)
            got = np.array([row[3] for row in actual])
            want = np.array([row[3] for row in rows])
            nt.assert_allclose(got, want, rtol=RELATIVE_TOLERANCE,
                               atol=ABSOLUTE_TOLERANCE)

    # And the reason the flip did not fire.
    track = pyBigWig.open(extra_track)
    try:
        rows = _read_bedgraph(active[0])
        for output_index in (0, 1):
            values = np.array([row[3] for row in _read_bedgraph(active[output_index])])
            positive = np.flatnonzero(values > 0)
            negative = np.flatnonzero(values < 0)
            assert positive.size > 0 and negative.size > 0
            positive_sum = 0.0
            for index in positive:
                statistic = track.stats(rows[index][0], rows[index][1],
                                        rows[index][2])[0]
                if statistic:
                    positive_sum += statistic
            assert positive_sum == 0.0, \
                'eigenvector {} now has coverage in its positive bins, so the ' \
                'guard no longer blocks the flip and this test needs ' \
                'rewriting'.format(output_index + 1)
    finally:
        track.close()
    for path in active + inactive:
        os.unlink(path)


def test_pca_bigwig_payload_matches_the_bedgraph_at_full_precision():
    """The bigwig writer carries the same float64 values as the bedgraph.

    are_files_equal_bigwig above compares np.absolute(...) at decimal=0 and
    returns True unconditionally when the first interval of the two files
    already agrees, so the bigwig payload is effectively unchecked. This
    compares every interval of every chromosome against the bedgraph the same
    invocation would write.
    """
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    bigwigs = _temporary_names(2, '.bw')
    bedgraphs = _temporary_names(2, '.bedgraph')

    base = ("--matrix {} --outputFileName {} {} -f {} "
            "--whichEigenvectors 1 2 --method lieberman")
    _run_pca(base.format(matrix, bigwigs[0], bigwigs[1], 'bigwig').split())
    _run_pca(base.format(matrix, bedgraphs[0], bedgraphs[1], 'bedgraph').split())

    for bigwig_path, bedgraph_path in zip(bigwigs, bedgraphs):
        rows = _read_bedgraph(bedgraph_path)
        handle = pyBigWig.open(bigwig_path)
        try:
            expected_by_chrom = {}
            for chrom, start, end, value in rows:
                expected_by_chrom.setdefault(chrom, []).append((start, end, value))
            assert set(handle.chroms().keys()) == set(expected_by_chrom.keys())
            for chrom, wanted in expected_by_chrom.items():
                intervals = handle.intervals(chrom)
                assert len(intervals) == len(wanted), \
                    'chromosome {}: {} intervals, expected {}' \
                    .format(chrom, len(intervals), len(wanted))
                for got, want in zip(intervals, wanted):
                    assert (got[0], got[1]) == (want[0], want[1])
                    # bigWig stores float32 and the bedgraph is quantised to
                    # twelve decimals, so the comparison is bounded by the
                    # coarser of the two, not by float64.
                    nt.assert_allclose(got[2], np.float32(want[2]),
                                       rtol=1e-6, atol=1e-12,
                                       err_msg='chromosome {}'.format(chrom))
        finally:
            handle.close()
    for path in bigwigs + bedgraphs:
        os.unlink(path)


def test_pca_intermediate_matrices_full_precision():
    """--pearsonMatrix and --obsexpMatrix, compared value by value.

    The existing test compares matrix.data at decimal=0 against files checked
    into test_data, and compares the two data arrays without checking that
    they describe the same sparsity pattern.
    """
    matrix = ROOT + "small_test_matrix.h5"
    outputs = _temporary_names(2, '.bedgraph')
    pearson_file = _temporary_names(1, '.h5')[0]
    obsexp_file = _temporary_names(1, '.h5')[0]
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --method lieberman "
            "--chromosomes chrX chrXHet --pearsonMatrix {} --obsexpMatrix {}"
            .format(matrix, outputs[0], outputs[1], pearson_file,
                    obsexp_file)).split()
    _run_pca(args)

    written_pearson = hm.hiCMatrix(pearson_file).matrix.todense()
    written_obsexp = hm.hiCMatrix(obsexp_file).matrix.todense()

    offset = 0
    for _name, intervals, obs_exp, pearson, _cov in \
            _covariance_per_chromosome(matrix, ['chrX', 'chrXHet'], 'lieberman'):
        size = len(intervals)
        block = slice(offset, offset + size)
        nt.assert_allclose(np.asarray(written_obsexp[block, block]), obs_exp,
                           rtol=RELATIVE_TOLERANCE, atol=ABSOLUTE_TOLERANCE)
        nt.assert_allclose(np.asarray(written_pearson[block, block]), pearson,
                           rtol=RELATIVE_TOLERANCE, atol=ABSOLUTE_TOLERANCE)
        offset += size
    assert offset == written_pearson.shape[0]

    for path in outputs + [pearson_file, obsexp_file]:
        os.unlink(path)


def test_pca_ignore_masked_bins_full_precision():
    """--ignoreMaskedBins, at float64 rather than at decimal=0."""
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    outputs = _temporary_names(2, '.bedgraph')
    args = ("--matrix {} --outputFileName {} {} -f bedgraph "
            "--whichEigenvectors 1 2 --ignoreMaskedBins --method lieberman"
            .format(matrix, outputs[0], outputs[1])).split()
    _run_pca(args)

    expected = _reference_eigenvectors(matrix, None, 'lieberman', ['1', '2'],
                                       pIgnoreMaskedBins=True)
    for path, rows in zip(outputs, expected):
        assert_bedgraph_equal(path, rows)
    for path in outputs:
        os.unlink(path)


def test_pca_output_file_count_must_match_eigenvector_count():
    """The argument check at hicPCA.py:242-248 exits 1. Nothing tests it."""
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    outputs = _temporary_names(1, '.bedgraph')
    args = ("--matrix {} --outputFileName {} -f bedgraph "
            "--whichEigenvectors 1 2 --method lieberman"
            .format(matrix, outputs[0])).split()
    with pytest.raises(SystemExit) as raised:
        _run_pca(args)
    assert raised.value.code == 1
    for path in outputs:
        os.unlink(path)


def are_files_equal(file1, file2):
    equal = True
    with open(file1) as textfile1, open(file2) as textfile2:
        for x, y in zip(textfile1, textfile2):
            if x.startswith('File'):
                continue
            if x != y:
                # handle the case of flipped values
                split_x = x.split('\t')
                split_y = y.split('\t')
                if split_x[0] == split_y[0] and split_x[1] == split_y[1] and split_x[2] == split_y[2]:
                    # to ignore rounding errors after 2th digit
                    if 0 <= abs(abs(np.complex128(split_x[3].strip()).real) - abs(np.complex128(split_y[3].strip()).real)) <= 0.01:
                        continue
                    else:
                        log.debug('split_x {} split_y {}'.format(split_x, split_y))
                equal = False
                break
    return equal


def are_files_equal_bigwig(pFile1, pFile2, pChromosomeList):

    bw_file1 = pyBigWig.open(pFile1)
    bw_file2 = pyBigWig.open(pFile2)

    for chrom in pChromosomeList:
        try:
            bins_list_file1 = bw_file1.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)
        try:
            bins_list_file2 = bw_file2.intervals(chrom)
        except Exception:
            log.debug("Chrom not found: {}", chrom)
        # sometimes the values are + / - flipped

        if bins_list_file1 is not None and bins_list_file1[0][2] != bins_list_file2[0][2]:
            bins_list_file1 = np.array(bins_list_file1)
            bins_list_file2 = np.array(bins_list_file2)
            bins_list_file1[:][2] *= -1
        if bins_list_file1 is None and bins_list_file2 is None:
            return True
        nt.assert_array_almost_equal(np.absolute(bins_list_file1),
                                     np.absolute(bins_list_file2),
                                     decimal=DELTA_DECIMAL)
    return True


def test_pca_bedgraph_lieberman():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    args = "--matrix {} --outputFileName {} {} -f bedgraph --whichEigenvectors  1 2 --method lieberman"\
        .format(matrix, pca1.name, pca2.name).split()
    hicPCA.main(args)
    # compute(hicPCA.main, args, 5)
    assert are_files_equal(ROOT + "hicPCA/pca1.bedgraph", pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2.bedgraph", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bedgraph_lieberman_ignore_masked_bins():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    args = "--matrix {} --outputFileName {} {} -f bedgraph --whichEigenvectors  1 2 \
           --ignoreMaskedBins --method lieberman".format(matrix, pca1.name, pca2.name).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    assert are_files_equal(ROOT + "hicPCA/pca1_ignoredMaskedBins.bedgraph", pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2_ignoredMaskedBins.bedgraph", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bigwig_lieberman():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix_50kb_res.h5"
    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 --method lieberman"\
        .format(matrix, pca1.name, pca2.name).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    chrom_list = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr2RHet',
                  'chr3RHet', 'chr2LHet', 'chr4', 'chrU', 'chrX',
                  'chrXHet', 'chr3LHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2.bw",
                                  pca2.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bedgraph_lieberman_gene_density():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bedgraph --whichEigenvectors 1 2 --method lieberman \
    --extraTrack {} --chromosomes {}"\
    .format(matrix, pca1.name, pca2.name, gene_track, chromosomes).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    assert are_files_equal(ROOT + "hicPCA/pca1_gene_track.bedgraph", pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2_gene_track.bedgraph", pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bigwig_lieberman_gene_density():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 --method lieberman \
    --extraTrack {} --chromosomes {}"\
    .format(matrix, pca1.name, pca2.name, gene_track, chromosomes).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    chrom_list = ['chrX', 'chrXHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1_gene_track.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2_gene_track.bw",
                                  pca2.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bigwig_lieberman_gene_density_intermediate_matrices():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)
    pearson_matrix = NamedTemporaryFile(suffix='.h5', delete=False)
    obs_exp_matrix = NamedTemporaryFile(suffix='.h5', delete=False)
    pca1.close()
    pca2.close()
    pearson_matrix.close()
    obs_exp_matrix.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 --method lieberman \
    --extraTrack {} --chromosomes {} --pearsonMatrix {} --obsexpMatrix {}"\
    .format(matrix, pca1.name, pca2.name, gene_track, chromosomes,
            pearson_matrix.name, obs_exp_matrix.name).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    chrom_list = ['chrX', 'chrXHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1_gene_track.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2_gene_track.bw",
                                  pca2.name, chrom_list)

    test_pearson = hm.hiCMatrix(ROOT + "hicPCA/pearson_intermediate.h5")

    new_pearson = hm.hiCMatrix(pearson_matrix.name)

    test_obs_exp = hm.hiCMatrix(ROOT + "hicPCA/obs_exp_intermediate.h5")

    new_obs_exp = hm.hiCMatrix(obs_exp_matrix.name)
    nt.assert_array_almost_equal(test_pearson.matrix.data,
                                 new_pearson.matrix.data,
                                 decimal=DELTA_DECIMAL)
    nt.assert_array_almost_equal(test_obs_exp.matrix.data,
                                 new_obs_exp.matrix.data,
                                 decimal=DELTA_DECIMAL)

    # assert are_files_equal_bigwig(ROOT + "hicPCA/pearson_intermediate.h5", pearson_matrix.name, chrom_list)
    # assert are_files_equal_bigwig(ROOT + "hicPCA/obs_exp_intermediate.h5", obs_exp_matrix.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)
    os.unlink(obs_exp_matrix.name)
    os.unlink(pearson_matrix.name)


def test_pca_bigwig_gene_density_intermediate_matrices_norm():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)
    pearson_matrix = NamedTemporaryFile(suffix='.cool', delete=False)
    obs_exp_matrix = NamedTemporaryFile(suffix='.cool', delete=False)
    pca1.close()
    pca2.close()
    pearson_matrix.close()
    obs_exp_matrix.close()
    matrix = ROOT + "small_test_matrix.h5"
    gene_track = ROOT + 'dm3_genes.bed.gz'
    chromosomes = 'chrX chrXHet'
    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 \
    --extraTrack {} --chromosomes {} --pearsonMatrix {} --obsexpMatrix {} \
    --method dist_norm --ligation_factor"\
    .format(matrix, pca1.name, pca2.name, gene_track, chromosomes,
            pearson_matrix.name, obs_exp_matrix.name).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)

    chrom_list = ['chrX', 'chrXHet']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1_gene_track.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2_gene_track.bw",
                                  pca2.name, chrom_list)

    test_pearson = hm.hiCMatrix(ROOT + "hicPCA/pearson_norm.cool")

    new_pearson = hm.hiCMatrix(pearson_matrix.name)

    test_obs_exp = hm.hiCMatrix(ROOT + "hicPCA/obsexp_norm.cool")

    new_obs_exp = hm.hiCMatrix(obs_exp_matrix.name)

    # load h5 matrices to compare if they are store the same data
    test_obs_exp_h5 = hm.hiCMatrix(ROOT + "hicPCA/obsexp_norm.h5")

    test_pearson_h5 = hm.hiCMatrix(ROOT + "hicPCA/pearson_norm.h5")

    nt.assert_array_almost_equal(test_pearson.matrix.data,
                                 new_pearson.matrix.data,
                                 decimal=DELTA_DECIMAL)
    nt.assert_array_almost_equal(test_obs_exp.matrix.data,
                                 new_obs_exp.matrix.data,
                                 decimal=DELTA_DECIMAL)

    nt.assert_array_almost_equal(test_pearson_h5.matrix.data,
                                 new_pearson.matrix.data,
                                 decimal=DELTA_DECIMAL)
    nt.assert_array_almost_equal(test_obs_exp_h5.matrix.data,
                                 new_obs_exp.matrix.data,
                                 decimal=DELTA_DECIMAL)

    os.unlink(pca1.name)
    os.unlink(pca2.name)
    os.unlink(obs_exp_matrix.name)
    os.unlink(pearson_matrix.name)


def test_pca_bigwig_lieberman_histoneMark_track():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    extra_track = ROOT + 'bigwig_chrx_2e6_5e6.bw'
    chromosomes = 'chrX '
    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 --method lieberman \
    --extraTrack {} --chromosomes {}"\
    .format(matrix, pca1.name, pca2.name, extra_track, chromosomes).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)
    chrom_list = ['chrX']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca1_chip_track.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca2_chip_track.bw",
                                  pca2.name, chrom_list)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_bedgraph_lieberman_histoneMark_track():
    pca1 = NamedTemporaryFile(suffix='.bedgraph', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bedgraph', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "small_test_matrix.h5"
    extra_track = ROOT + 'bigwig_chrx_2e6_5e6.bw'
    chromosomes = 'chrX '
    args = "--matrix {} --outputFileName {} {} -f bedgraph --whichEigenvectors  1 2 --method lieberman \
    --extraTrack {} --chromosomes {}"\
    .format(matrix, pca1.name, pca2.name, extra_track, chromosomes).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)
    assert are_files_equal(ROOT + "hicPCA/pca1_chip_track.bedgraph",
                           pca1.name)
    assert are_files_equal(ROOT + "hicPCA/pca2_chip_track.bedgraph",
                           pca2.name)

    os.unlink(pca1.name)
    os.unlink(pca2.name)


def test_pca_extratrack_extends_chromosomesize():
    pca1 = NamedTemporaryFile(suffix='.bw', delete=False)
    pca2 = NamedTemporaryFile(suffix='.bw', delete=False)

    pca1.close()
    pca2.close()
    matrix = ROOT + "hicPCA/mm9_reduced_chr1.cool"
    extra_track = ROOT + 'hicPCA/mm9_genes_sorted.bed12'

    args = "--matrix {} --outputFileName {} {} -f bigwig --whichEigenvectors  1 2 --method lieberman \
    --extraTrack {}"\
    .format(matrix, pca1.name, pca2.name, extra_track).split()
    # hicPCA.main(args)
    compute(hicPCA.main, args, 5)
    chrom_list = ['chr1']
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca_reduced_mm9_1.bw",
                                  pca1.name, chrom_list)
    assert are_files_equal_bigwig(ROOT + "hicPCA/pca_reduced_mm9_2.bw",
                                  pca2.name, chrom_list)
    os.unlink(pca1.name)
    os.unlink(pca2.name)
