import logging
from hicexplorer import chicViewpointBackgroundModel
import numpy as np
from sys import platform
from tempfile import NamedTemporaryFile, mkdtemp
import os
import pytest
import warnings

from hicmatrix import HiCMatrix as hm
from scipy.special import factorial, gammaln

from hicexplorer.lib import Viewpoint
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from hicexplorer.test.test_compute_function import compute


ROOT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "test_data/cHi-C/")
log = logging.getLogger(__name__)


def are_files_equal(file1, file2, delta=1, skip=0, eps=0.1):

    mismatches = 0
    with open(file1, 'r') as textfile1:
        with open(file2, 'r') as textfile2:

            file1_content = textfile1.readlines()
            file2_content = textfile2.readlines()

            for i, (line1, line2) in enumerate(zip(file1_content, file2_content)):
                if i < skip:
                    continue
                line1_list = np.array(line1.split('\t'))
                line2_list = np.array(line2.split('\t'))

                line1_list = line1_list.astype(np.float64)
                line2_list = line2_list.astype(np.float64)

                for value1, value2 in zip(line1_list, line2_list):
                    if np.abs(value1 - value2) < eps:
                        continue
                    else:
                        log.debug('{}'.format(line1_list))
                        log.debug('{}'.format(line2_list))

                        mismatches += 1

    if mismatches < delta:
        return True
    else:
        log.debug('mismatches: {}'.format(mismatches))
        return False


def test_compute_background_functional():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    # chicViewpointBackgroundModel.main(args)
    compute(chicViewpointBackgroundModel.main, args, 5)
    assert are_files_equal(ROOT + 'background.txt',
                           outfile.name, delta=700, skip=1)


def test_compute_background_functional_truncate_zeros():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {} --truncateZeros".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    # chicViewpointBackgroundModel.main(args)
    compute(chicViewpointBackgroundModel.main, args, 5)

    assert are_files_equal(ROOT + 'background_truncateZeros.txt',
                           outfile.name, delta=1000, skip=1)


def test_compute_background_number_of_lines():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {}".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                      ROOT + 'referencePoints.bed', outfile.name, 1).split()
    # chicViewpointBackgroundModel.main(args)
    compute(chicViewpointBackgroundModel.main, args, 5)

    length_background = 0
    length_background_outfile = 0

    with open(ROOT + 'background.txt') as textfile:
        file_content = textfile.readlines()
        length_background = len(file_content)
    with open(outfile.name) as textfile:
        file_content = textfile.readlines()
        length_background_outfile = len(file_content)

    assert np.abs(length_background - length_background_outfile) < 1


@pytest.mark.xfail
def test_compute_background_functional_fail():
    outfile = NamedTemporaryFile(suffix='.txt', delete=False)
    outfile.close()
    args = "--matrices {} {} --referencePoints {} -o {} -t {} --truncateZeros".format(ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
                                                                                      ROOT + 'referencePoints_qc.bed', outfile.name, 1).split()
    # chicViewpointBackgroundModel.main(args)
    compute(chicViewpointBackgroundModel.main, args, 5)

    assert are_files_equal(ROOT + 'background.txt',
                           outfile.name, delta=700, skip=1)


# ---------------------------------------------------------------------------
# Characterization tests added for the C++ port (cpp/AGENTS_CONTRACT.md rule 1).
#
# What the tests above actually constrain: are_files_equal compares every value
# with an absolute eps of 0.1, tolerates 700 (1,000 with --truncateZeros)
# mismatching values out of 5,005, and walks the files with zip(), so a
# truncated output passes. The mean value column holds numbers of order 1e-4,
# which an absolute eps of 0.1 cannot see, so that column is unconstrained. The
# only exact assertion is the line count. test_compute_background_functional_fail
# runs reference points that make the tool exit(1) and is marked xfail, so it
# cannot fail. --averageContactBin and --fixateRange are never varied.
#
# The fitted size and prob cannot be pinned to digits, and this is measured,
# not assumed: the scipy L-BFGS-B fit sits on a flat likelihood ridge. Refitting
# the 1,001 distributions of this data set with only the order of the values
# permuted (which changes nothing but the last bits of numpy's pairwise sums,
# exactly what the tool's own thread merge order does) moves `size` by a median
# of 55 percent over 11 orderings and `prob` by 6.7e-4; the unmodified tool at
# --threads 8 already writes different digits from --threads 1 on 987 of the
# 1,001 lines. What does not move is the likelihood the fit reaches: over the
# same 11 orderings, the negative log likelihood at the fitted parameters spans
# at most 0.095 nats at any position. The committed background.txt, written by
# an older scipy, is within 0.036 nats of today's fit everywhere. So the fit is
# pinned by the likelihood it reaches, with a tolerance of twice the measured
# ordering spread.
NLL_TOLERANCE = 0.2


class CollectingQueue:
    def put(self, pItem):
        self.item = pItem


def run_background_model(pExtraArguments='', pReferencePoints=None):
    outfile = os.path.join(mkdtemp(prefix='chicViewpointBackgroundModel_'), 'background.txt')
    args = "--matrices {} {} --referencePoints {} -o {} -t {} {}".format(
        ROOT + 'FL-E13-5_chr1.cool', ROOT + 'MB-E10-5_chr1.cool',
        pReferencePoints or ROOT + 'referencePoints.bed', outfile, 1,
        pExtraArguments).split()
    chicViewpointBackgroundModel.main(args)
    return outfile


def read_table(pPath):
    with open(pPath) as handle:
        header = handle.readline()
        rows = [line.rstrip('\n').split('\t') for line in handle]
    return header, rows


def background_distributions(pFixateRange=500000, pAverageContactBin=5, pTruncateZeros=False):
    """The value list the tool fits per relative position, in --threads 1 order."""
    args = chicViewpointBackgroundModel.parse_arguments().parse_args(
        ['--matrices', 'unused', '--referencePoints', 'unused',
         '--fixateRange', str(pFixateRange), '--averageContactBin', str(pAverageContactBin)])
    viewpoint = Viewpoint()
    reference_points, _ = viewpoint.readReferencePointFile(ROOT + 'referencePoints.bed')
    collected = {}
    for matrix in ('FL-E13-5_chr1.cool', 'MB-E10-5_chr1.cool'):
        viewpoint.hicMatrix = hm.hiCMatrix(ROOT + matrix)
        queue = CollectingQueue()
        chicViewpointBackgroundModel.compute_background(reference_points, viewpoint, args, queue)
        assert not isinstance(queue.item, str), queue.item
        for position, values in queue.item[0].items():
            collected.setdefault(position, []).extend(values)
    distributions = {}
    for position in sorted(collected):
        values = np.array(collected[position])
        if pTruncateZeros:
            values = values[values > 0.0]
        distributions[position] = values
    return distributions


def negative_log_likelihood(pX, pSize, pProb):
    """The objective of fit_nbinom.fit, evaluated at the given parameters."""
    infinitesimal = np.finfo(np.float64).eps
    count = pX.size
    with np.errstate(all='ignore'):
        result = np.sum(gammaln(pX + pSize)) \
            - np.sum(np.log(factorial(pX))) \
            - count * gammaln(pSize) \
            + count * pSize * np.log(pProb) \
            + np.sum(pX * np.log(1 - (pProb if pProb < 1 else 1 - infinitesimal)))
    return -result


def moment_start(pX):
    """fit_nbinom's starting point, R's fitdistr moment estimator."""
    mean = np.mean(pX)
    variance = np.var(pX)
    size = (mean ** 2) / (variance - mean) if variance > mean else 10
    prob = size / ((size + mean) if size + mean != 0 else 1)
    return size, prob


@pytest.mark.parametrize('pArguments, pExpected', [
    ('', 'background.txt'),
    ('--truncateZeros', 'background_truncateZeros.txt')])
def test_positions_max_and_mean_columns_are_exact(pArguments, pExpected):
    header_expected, rows_expected = read_table(ROOT + pExpected)
    header_observed, rows_observed = read_table(run_background_model(pArguments))
    assert header_observed == header_expected
    assert len(rows_observed) == len(rows_expected) == 1001
    for expected, observed in zip(rows_expected, rows_observed):
        assert len(observed) == 5
        assert (observed[0], observed[3], observed[4]) == \
            (expected[0], expected[3], expected[4]), (expected, observed)


@pytest.mark.parametrize('pArguments, pExpected, pTruncateZeros', [
    ('', 'background.txt', False),
    ('--truncateZeros', 'background_truncateZeros.txt', True)])
def test_fitted_parameters_reach_the_reference_likelihood(pArguments, pExpected, pTruncateZeros):
    _, rows_expected = read_table(ROOT + pExpected)
    _, rows_observed = read_table(run_background_model(pArguments))
    distributions = background_distributions(pTruncateZeros=pTruncateZeros)
    assert len(rows_observed) == len(rows_expected) == len(distributions)
    finite = 0
    degenerate = 0
    empty = 0
    for expected, observed in zip(rows_expected, rows_observed):
        values = distributions[int(observed[0]) // 1000]
        if values.size == 0:
            # --truncateZeros empties the distribution at +500 kb. fit_nbinom
            # then starts from mean nan: size 10, prob 10 / (10 + nan).
            assert observed[1:] == ['10.000000000000', 'nan', '0.000000000000', '0.000000000000']
            empty += 1
            continue
        reference = negative_log_likelihood(values, float(expected[1]), float(expected[2]))
        fitted = negative_log_likelihood(values, float(observed[1]), float(observed[2]))
        if np.isfinite(reference):
            assert abs(fitted - reference) <= NLL_TOLERANCE, \
                'position {}: NLL {} against the reference {}'.format(observed[0], fitted, reference)
            finite += 1
        else:
            # factorial overflows above 170, the objective is +inf for every
            # parameter pair, and fmin_l_bfgs_b returns its starting point.
            size, prob = moment_start(values)
            assert (observed[1], observed[2]) == ('{:.12f}'.format(size), '{:.12f}'.format(prob))
            degenerate += 1
    assert finite >= 990
    assert empty == (1 if pTruncateZeros else 0)
    assert finite + degenerate + empty == 1001


def test_fixate_range_and_average_contact_bin_are_honoured():
    """An even window, whose smoothing is asymmetric, and a smaller fixateRange."""
    header, rows = read_table(run_background_model('--fixateRange 200000 --averageContactBin 4'))
    assert header == 'Relative position\tsize nbinom\tprob nbinom\tmax value\tmean value\n'
    assert [int(row[0]) for row in rows] == list(range(-200000, 201000, 1000))
    distributions = background_distributions(pFixateRange=200000, pAverageContactBin=4)
    averages = {position: np.average(values) for position, values in distributions.items()}
    total = 0
    for position in sorted(averages):
        total += averages[position]
    improved = 0
    for row in rows:
        position = int(row[0]) // 1000
        values = distributions[position]
        assert row[3] == '{:.12f}'.format(np.max(values))
        assert row[4] == '{:.12f}'.format(averages[position] / total)
        start = moment_start(values)
        start_nll = negative_log_likelihood(values, *start)
        fitted_nll = negative_log_likelihood(values, float(row[1]), float(row[2]))
        if np.isfinite(start_nll):
            assert fitted_nll <= start_nll + 1e-6
            improved += fitted_nll < start_nll - 1e-3
    assert improved > 300


def expected_smoothing(pData, pWindow):
    """smoothInteractionValues written out element by element.

    For an odd window the average runs over pWindow // 2 values on either side,
    clipped at both ends. For an even window it takes one value fewer upstream
    than downstream, except at the right border, where viewpoint.py:561
    averages pData[-(i + window_size + 1):] and so takes the full half window
    upstream. That asymmetry is pinned here as the Python behaves.
    """
    half = pWindow // 2
    upstream = half - 1 if pWindow % 2 == 0 else half
    size = len(pData)
    result = np.empty(size)
    for i in range(size):
        if i >= size - half:
            begin, end = i - half, size
        else:
            begin, end = max(i - upstream, 0), i + half + 1
        result[i] = np.mean(pData[begin:end])
    return result


@pytest.mark.parametrize('pWindow', [1, 2, 3, 4, 5, 6])
def test_smoothing_window_semantics(pWindow):
    data = np.arange(1, 22, dtype=np.float64) ** 1.5
    observed = Viewpoint().smoothInteractionValues(data, pWindow)
    np.testing.assert_array_equal(observed, expected_smoothing(data, pWindow))
