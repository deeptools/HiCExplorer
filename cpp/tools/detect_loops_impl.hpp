// The computational core of hicDetectLoops, separated from the tool so that
// the unit tests can reach it (cpp/tests/test_detect_loops.cpp).
//
// Everything here is a transcription of hicexplorer/hicDetectLoops.py and of
// the two hicexplorer/utilities.py helpers it calls. The Python is followed
// literally, including several behaviours that are defects; each is marked
// QUIRK with the line it comes from, and none is silently improved
// (cpp/AGENTS_CONTRACT.md rule 4).
//
// The matrix this operates on is the state hicDetectLoops.py:855-865 leaves
// behind: triu(matrix) with the main diagonal removed, so every stored entry
// has column > row. It is carried in a hicx::CsrMatrix marked Symmetry::Full,
// because at that point the matrix really is the bare upper triangle and not a
// compressed symmetric one: the neighbourhood windows the tool cuts out of it
// read zeros below the diagonal and that is what the Python sees too.
//
// Determinism (cpp/OPTIMIZATION.md 3). Every parallel loop here runs over
// *individual independent items*: one distance, or one candidate. The item
// index comes from the data alone, each item writes only into its own
// preallocated slot, and the slots are combined in index order afterwards. No
// reduction crosses an item boundary, so the result cannot depend on the
// worker count or on the completion order. That is strictly stronger than the
// Python, which cuts the work into exactly `--threadsPerChromosome` ranges and
// therefore has a partition that moves with the thread count; it happens to
// produce the same answer because its items are independent too, which is the
// property this file relies on and the unit tests pin.

#ifndef HICX_TOOLS_DETECT_LOOPS_IMPL_HPP
#define HICX_TOOLS_DETECT_LOOPS_IMPL_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/parallel.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::loops {

// --------------------------------------------------------------------------
// numpy semantics the Python leans on without meaning to

// numpy basic slicing on one axis of length `length`. A negative bound counts
// from the end and an out of range bound is clamped, which is what makes
// several of the window slices in candidate_region_test_thread behave the way
// they do when the peak sits within --peakWidth of the edge of its
// neighbourhood: `neighborhood[:pr - pw, :]` with `pr < pw` is not empty, it
// is everything except the last `pw - pr` rows.
struct AxisRange {
    std::int64_t begin = 0;
    std::int64_t end = 0;
    [[nodiscard]] std::int64_t size() const noexcept {
        return end > begin ? end - begin : 0;
    }
};

[[nodiscard]] inline AxisRange numpy_slice(std::int64_t length, std::int64_t start,
                                           std::int64_t stop) {
    const auto normalise = [length](std::int64_t value) {
        if (value < 0) {
            value += length;
            if (value < 0) {
                value = 0;
            }
        }
        if (value > length) {
            value = length;
        }
        return value;
    };
    AxisRange range{normalise(start), normalise(stop)};
    if (range.end < range.begin) {
        range.end = range.begin;
    }
    return range;
}

// A dense row major block, the `.toarray()` of a slice of the sparse matrix.
struct DenseBlock {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<double> values;

    [[nodiscard]] double at(std::int64_t row, std::int64_t col) const {
        return values[static_cast<std::size_t>(row * cols + col)];
    }
    [[nodiscard]] std::int64_t size() const noexcept { return rows * cols; }
};

// matrix[row0:row1, col0:col1].toarray(), with the bounds already valid.
[[nodiscard]] DenseBlock dense_block(const CsrMatrix& matrix, std::int64_t row0,
                                     std::int64_t row1, std::int64_t col0,
                                     std::int64_t col1);

// block[rows, cols].flatten() with numpy slicing on both axes.
[[nodiscard]] std::vector<double> flatten_slice(const DenseBlock& block,
                                                std::int64_t row_start,
                                                std::int64_t row_stop,
                                                std::int64_t col_start,
                                                std::int64_t col_stop);

// matrix[row, col] for a single element.
[[nodiscard]] double element(const CsrMatrix& matrix, std::int64_t row,
                             std::int64_t col);

// np.mean, which is np.add.reduce divided by the count, so the reduction order
// is numpy's pairwise one and not a plain loop.
[[nodiscard]] inline double numpy_mean(std::span<const double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return simd::pairwise_sum(values.data(), values.size()) /
           static_cast<double>(values.size());
}

// --------------------------------------------------------------------------
// The stored dtype, because obs/exp casts its result back to it

enum class ValueKind { Integer, Float32, Float64 };

[[nodiscard]] ValueKind value_kind(const std::string& dtype);

// numpy's astype(<the input dtype>).
//
// QUIRK (hicexplorer/utilities.py:582,588). obs_exp_matrix records
// `data_type = type(pSubmatrix.data[0])` *before* it casts the data to
// float32, and casts the quotient back to it at the end. For an integer matrix
// -- which is what an h5 count matrix and an unrestricted cool load both are --
// the observed over expected ratio is therefore **truncated to an integer**.
// It is not a rounding artefact at the edge: on small_test_matrix.h5 the
// expected value per distance is around 0.05, so the ratios are in the
// hundreds and truncation merely discards their fractional part, which is why
// the tool still finds loops instead of dividing everything to zero. Change
// this and the call set changes.
[[nodiscard]] inline double cast_to_kind(double value, ValueKind kind) {
    switch (kind) {
        case ValueKind::Integer:
            return std::trunc(value);
        case ValueKind::Float32:
            return static_cast<double>(static_cast<float>(value));
        case ValueKind::Float64:
        default:
            return value;
    }
}

// --------------------------------------------------------------------------
// Grouping the stored entries by genomic distance

// The stored positions of the matrix, bucketed by |row - column|, each bucket
// in CSR storage order. That is the order np.sum(data[distance == d]) reduces
// in, and the order the per distance obs/exp distribution is built in, so both
// the expected value and the negative binomial fit depend on it.
struct DistanceGroups {
    std::int64_t min_distance = 0;
    std::int64_t max_distance = -1;
    // Bucket d occupies positions[offsets[d - min] .. offsets[d - min + 1]).
    std::vector<std::int64_t> offsets;
    std::vector<std::int64_t> positions;

    [[nodiscard]] bool empty() const noexcept { return max_distance < min_distance; }
    [[nodiscard]] std::int64_t bucket_count() const noexcept {
        return empty() ? 0 : max_distance - min_distance + 1;
    }
    [[nodiscard]] std::span<const std::int64_t> bucket(std::int64_t index) const {
        const std::size_t begin = static_cast<std::size_t>(offsets[static_cast<std::size_t>(index)]);
        const std::size_t end =
            static_cast<std::size_t>(offsets[static_cast<std::size_t>(index) + 1]);
        return std::span<const std::int64_t>(positions.data() + begin, end - begin);
    }
};

[[nodiscard]] DistanceGroups group_by_distance(const CsrMatrix& matrix);

// --------------------------------------------------------------------------
// Observed over expected

// hicexplorer.utilities.expected_interactions: the sum of the values at each
// genomic distance divided by `occurrences`.
//
// QUIRK (utilities.py:362). `occurrences = np.arange(shape[0] + 1, 1, -1)`, so
// the divisor for distance d is n + 1 - d. The number of positions a diagonal
// of an n by n matrix actually has is n - d, so every expected value is a
// factor (n - d) / (n + 1 - d) too small. It is a constant factor per distance
// and the tool only ever compares ratios within one distance, so it shifts
// every obs/exp value at a given distance by the same amount, but the value
// itself is not the mean it is documented to be.
//
// Threaded over distances. Each distance is reduced by exactly one worker with
// numpy's pairwise order, so the result is independent of `threads`.
[[nodiscard]] std::vector<double> expected_interactions(const CsrMatrix& matrix,
                                                        const DistanceGroups& groups,
                                                        unsigned int threads);

// hicexplorer.utilities.expected_interactions_non_zero: the same sums divided
// by the number of stored entries at that distance rather than by the diagonal
// length. Accumulated **sequentially** in a Python loop
// (utilities.py:327-329), not with np.sum, so this one is a plain left to
// right sum and not the pairwise one.
[[nodiscard]] std::vector<double> expected_interactions_non_zero(
    const CsrMatrix& matrix);

// hicexplorer.utilities.obs_exp_matrix(pInplace=False, pToEpsilon=True).
//
// QUIRK (utilities.py:578). The expected value for a pair at distance d is
// looked up at index ceil(d / 2), not at d. Distances 1 and 2 share the
// expected value of distance 1, 3 and 4 share that of distance 2, and so on.
// The docstring says "expected contacts for loci at that genomic distance".
//
// The value is computed as float64(float32(count)) / expected, NaN and
// infinity are replaced by 1e-6 (utilities.py:137,144, pToEpsilon), and the
// result is cast back to the input dtype.
[[nodiscard]] CsrMatrix obs_exp_matrix(const CsrMatrix& matrix,
                                       const DistanceGroups& groups,
                                       unsigned int threads);

// hicexplorer.utilities.obs_exp_matrix_non_zero(pToEpsilon=True), for
// --expected mean_nonzero and mean_nonzero_ligation.
//
// Unlike obs_exp_matrix this one leaves the result in float32, because it
// assigns each quotient back into an array it has already cast to float32
// (utilities.py:533,540) and never casts back. Its epsilon is 1e-9, not the
// 1e-6 of obs_exp_matrix (utilities.py:543).
[[nodiscard]] CsrMatrix obs_exp_matrix_non_zero(const CsrMatrix& matrix,
                                                bool ligation_factor);

// --------------------------------------------------------------------------
// The per distance negative binomial preselection

// One fitted distribution, kept so that the harness and the tests can look at
// the parameters rather than only at the calls they produce.
struct DistanceFit {
    std::int64_t distance = 0;
    std::size_t count = 0;
    double size = 0.0;
    double prob = 0.0;
    // True when np.sum(np.log(factorial(X))) overflowed to infinity, which
    // makes the objective constant and leaves fit_nbinom at its starting
    // point. See the note on preselect_candidates.
    bool degenerate = false;
};

struct PreselectionResult {
    // One flag per stored entry of the obs/exp matrix, in CSR order.
    std::vector<char> mask;
    std::vector<DistanceFit> fits;
};

// hicDetectLoops.compute_p_values_mask (:145-180).
//
// For every genomic distance present in the matrix, fit a negative binomial to
// *all* obs/exp values at that distance, then compute 1 - cnb.cdf for the
// values that reach pObsExpThreshold and keep those whose p-value is at most
// the preselection threshold.
//
// A note on reproducibility that decides this tool's equivalence class. The
// objective fit_nbinom minimises contains the constant term
// np.sum(np.log(factorial(X))). scipy's factorial overflows to infinity above
// 170, so as soon as one obs/exp value at a distance exceeds that, the
// objective is +inf for every parameter pair, the forward difference gradient
// is inf - inf = NaN, and fmin_l_bfgs_b stops at iteration zero and returns
// its starting point unchanged. The fit is then exactly R's fitdistr moment
// estimator and is bit reproducible. That is the case for essentially every
// distance of an integer obs/exp matrix, whose values run into the thousands.
// Where it does not hold the optimiser really runs, and there the reference is
// not reproducible against itself: see the EN measurement in the report.
//
// `threshold_by_genomic_distance` empty means the single float threshold in
// `threshold` applies; otherwise the value for `distance * resolution` is
// looked up, and a missing key is a hard error, as the Python's KeyError is.
[[nodiscard]] PreselectionResult preselect_candidates(
    const CsrMatrix& obs_exp, const DistanceGroups& groups, double threshold,
    const std::map<std::int64_t, double>& threshold_by_genomic_distance,
    std::int64_t resolution, double obs_exp_threshold, unsigned int threads);

// --------------------------------------------------------------------------
// Candidate selection

struct Candidate {
    std::int64_t row = 0;
    std::int64_t col = 0;
};

// hicDetectLoops.neighborhood_merge (:382-496): keep a candidate only when it
// is the maximum of its own (2 * window + 1)^2 neighbourhood.
[[nodiscard]] std::vector<Candidate> neighborhood_merge(
    const std::vector<Candidate>& candidates, std::int64_t window_size,
    const CsrMatrix& obs_exp, unsigned int threads);

struct RegionTestResult {
    std::vector<Candidate> candidates;
    std::vector<double> pvalues;
};

// hicDetectLoops.candidate_region_test (:499-752): the donut test.
//
// QUIRK (:534). `np.array(np.where(neighborhood == value)).flatten()` produces
// [row_0, .., row_{k-1}, col_0, .., col_{k-1}] for k matches, and the code then
// reads peak_region[0] and peak_region[1] as a (row, column) pair. That is only
// the position of the peak when the value occurs exactly once in the window.
// With two or more matches peak_region[1] is the *second row* index, so the
// window is cut around a position that is not the candidate. Reproduced.
[[nodiscard]] RegionTestResult candidate_region_test(
    const CsrMatrix& obs_exp, const std::vector<Candidate>& candidates,
    std::int64_t window_size, double p_value, std::int64_t peak_window_size,
    unsigned int threads);

}  // namespace hicx::loops

#endif  // HICX_TOOLS_DETECT_LOOPS_IMPL_HPP
