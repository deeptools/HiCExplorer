// See hicx/find_tads.hpp. Line references are to hicexplorer/hicFindTADs.py
// of HiCExplorer 3.7.6.

#include "hicx/find_tads.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"
#include "hicx/python_error.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/utilities.hpp"

namespace hicx::find_tads {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// The (start, stop) of slice(start, stop) on a sequence of `length` items, step 1.
std::pair<std::int64_t, std::int64_t> python_slice(std::int64_t start, std::int64_t stop,
                                                   std::int64_t length) {
    const auto clamp = [length](std::int64_t value) {
        if (value < 0) {
            value += length;
            return value < 0 ? std::int64_t{0} : value;
        }
        return value > length ? length : value;
    };
    const std::int64_t first = clamp(start);
    const std::int64_t last = clamp(stop);
    return {first, std::max(first, last)};
}

}  // namespace

// get_incremental_step_size (hicFindTADs.py:283-301).
std::vector<std::int64_t> incremental_step_size(std::int64_t min_window,
                                                std::int64_t max_window,
                                                std::int64_t start_step) {
    std::vector<std::int64_t> steps;
    std::int64_t step = -1;
    while (true) {
        ++step;
        const std::int64_t increment =
            min_window + static_cast<std::int64_t>(static_cast<double>(start_step) *
                                                   std::pow(static_cast<double>(step), 1.5));
        if (step > 1 && !steps.empty() && increment == steps.back()) {
            continue;
        }
        if (increment > max_window) {
            break;
        }
        steps.push_back(increment);
    }
    return steps;
}

// --------------------------------------------------------------------------
// The matrix operations of hiCMatrix that hicFindTADs performs.

// hiCMatrix.diagflat(value=0): the diagonal is replaced by zeros, and scipy's
// sparse addition drops results that are exactly zero, so the entries simply
// disappear.
void drop_diagonal(hicx::CsrMatrix& matrix) {
    const std::int64_t n = matrix.rows();
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    indices.reserve(matrix.stored_nnz());
    values.reserve(matrix.stored_nnz());
    for (std::int64_t row = 0; row < n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (static_cast<std::int64_t>(matrix.indices()[k]) == row) {
                continue;
            }
            indices.push_back(matrix.indices()[k]);
            values.push_back(matrix.data()[k]);
        }
        indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(values.size());
    }
    const hicx::Symmetry symmetry = matrix.symmetry();
    matrix = hicx::CsrMatrix(n, n, std::move(indptr), std::move(indices),
                             std::move(values), matrix.dtype());
    matrix.set_symmetry(symmetry);
}

// The bookkeeping half of hiCMatrix.maskBins / restoreMaskedBins as the tool
// needs it. The Python keeps the removed bins so that save() can put them back
// as empty rows, and hicFindTADs relies on that: when enlarge_bins does not
// change the bin table, the masked bins reappear in the saved z-score matrix.

MaskState mask_bins(hicx::MatrixData& data, const std::vector<std::int64_t>& bin_ids) {
    MaskState state;
    if (bin_ids.empty()) {
        return state;  // maskBins returns before touching anything
    }
    std::set<std::int64_t> masked(bin_ids.begin(), bin_ids.end());
    // The union with the NaN bins the matrix already carries
    // (HiCMatrix.py:794-799).
    masked.insert(data.nan_bins.begin(), data.nan_bins.end());
    data.nan_bins.clear();

    std::vector<std::int64_t> rows;
    rows.reserve(static_cast<std::size_t>(data.matrix.rows()));
    for (std::int64_t bin = 0; bin < data.matrix.rows(); ++bin) {
        if (masked.count(bin) == 0) {
            rows.push_back(bin);
        }
    }
    data.matrix = hicx::select_bins(data.matrix, rows);

    state.orig_bin_ids = rows;
    state.orig_bin_ids.insert(state.orig_bin_ids.end(), masked.begin(), masked.end());
    state.orig_cut_intervals.reserve(state.orig_bin_ids.size());
    for (const std::int64_t bin : state.orig_bin_ids) {
        state.orig_cut_intervals.push_back(
            data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    std::vector<hicx::CutInterval> kept;
    kept.reserve(rows.size());
    for (const std::int64_t bin : rows) {
        kept.push_back(data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    data.cut_intervals = std::move(kept);
    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(rows.size());
        for (const std::int64_t bin : rows) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(bin)]);
        }
        data.correction_factors = std::move(factors);
    }
    return state;
}

void restore_masked_bins(hicx::MatrixData& data, MaskState& state) {
    if (state.orig_bin_ids.empty()) {
        return;
    }
    const std::int64_t m = data.matrix.rows();
    const std::int64_t total = static_cast<std::int64_t>(state.orig_bin_ids.size());
    const std::int64_t added = total - m;

    // argsort(orig_bin_ids), then matrix[rows, :][:, rows] on the matrix padded
    // with `added` empty rows and columns.
    std::vector<std::int64_t> order(static_cast<std::size_t>(total));
    for (std::int64_t i = 0; i < total; ++i) {
        order[static_cast<std::size_t>(i)] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&state](std::int64_t a, std::int64_t b) {
                         return state.orig_bin_ids[static_cast<std::size_t>(a)] <
                                state.orig_bin_ids[static_cast<std::size_t>(b)];
                     });

    if (added > 0) {
        // Pad: the extra rows and columns hold nothing, so only the shape and
        // the row offsets change.
        std::vector<std::int64_t> indptr(static_cast<std::size_t>(total) + 1, 0);
        for (std::int64_t row = 0; row <= m; ++row) {
            indptr[static_cast<std::size_t>(row)] =
                data.matrix.indptr()[static_cast<std::size_t>(row)];
        }
        for (std::int64_t row = m + 1; row <= total; ++row) {
            indptr[static_cast<std::size_t>(row)] =
                static_cast<std::int64_t>(data.matrix.stored_nnz());
        }
        std::vector<std::int32_t> indices(data.matrix.indices());
        std::vector<double> values(data.matrix.data());
        const hicx::Symmetry symmetry = data.matrix.symmetry();
        data.matrix = hicx::CsrMatrix(total, total, std::move(indptr), std::move(indices),
                                      std::move(values), "float64");
        data.matrix.set_symmetry(symmetry);
    }
    data.matrix = hicx::select_bins(data.matrix, order);

    std::vector<hicx::CutInterval> restored;
    restored.reserve(order.size());
    for (const std::int64_t index : order) {
        restored.push_back(state.orig_cut_intervals[static_cast<std::size_t>(index)]);
    }
    data.cut_intervals = std::move(restored);
    data.nan_bins.assign(state.orig_bin_ids.begin() + static_cast<std::ptrdiff_t>(m),
                         state.orig_bin_ids.end());
    std::sort(data.nan_bins.begin(), data.nan_bins.end());

    if (data.correction_factors.has_value()) {
        std::vector<double> padded = *data.correction_factors;
        padded.resize(static_cast<std::size_t>(total), kNaN);
        std::vector<double> reordered;
        reordered.reserve(padded.size());
        for (const std::int64_t index : order) {
            reordered.push_back(padded[static_cast<std::size_t>(index)]);
        }
        data.correction_factors = std::move(reordered);
    }
    state.orig_bin_ids.clear();
    state.orig_cut_intervals.clear();
}

// triu(m, 0) - triu(m, limit): entries at or beyond `limit` cancel, except
// NaN, which cannot cancel and stays stored.
void truncate_to_distance(hicx::CsrMatrix& matrix, std::int64_t limit) {
    const std::int64_t n = matrix.rows();
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    indices.reserve(matrix.stored_nnz());
    values.reserve(matrix.stored_nnz());
    for (std::int64_t row = 0; row < n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(matrix.indices()[k]);
            if (column < row) {
                continue;
            }
            const double value = matrix.data()[k];
            if (column - row >= limit && !std::isnan(value)) {
                continue;
            }
            if (value == 0.0) {
                continue;  // eliminate_zeros
            }
            indices.push_back(static_cast<std::int32_t>(column));
            values.push_back(value);
        }
        indptr[static_cast<std::size_t>(row) + 1] =
            static_cast<std::int64_t>(values.size());
    }
    const std::string dtype = matrix.dtype();
    matrix = hicx::CsrMatrix(n, n, std::move(indptr), std::move(indices),
                             std::move(values), dtype);
    matrix.set_symmetry(hicx::Symmetry::Full);
}

// --------------------------------------------------------------------------
// The TAD-separation score


void extract_block(const hicx::CsrMatrix& matrix, std::int64_t r0, std::int64_t r1,
                   std::int64_t c0, std::int64_t c1, DenseBlock& block) {
    block.rows = std::max<std::int64_t>(0, r1 - r0);
    block.columns = std::max<std::int64_t>(0, c1 - c0);
    const std::size_t count =
        static_cast<std::size_t>(block.rows) * static_cast<std::size_t>(block.columns);
    block.values.assign(count, 0.0);
    block.stored = 0;
    for (std::int64_t row = r0; row < r1; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t column = static_cast<std::int64_t>(matrix.indices()[k]);
            if (column < c0) {
                continue;
            }
            if (column >= c1) {
                break;  // column indices are sorted
            }
            block.values[static_cast<std::size_t>(row - r0) *
                             static_cast<std::size_t>(block.columns) +
                         static_cast<std::size_t>(column - c0)] = matrix.data()[k];
            ++block.stored;
        }
    }
}


// get_idx_of_bins_at_given_distance. Returns false where the Python's
// getRegionBinRange returns None and the caller ends up with a TypeError.
bool bins_at_distance(const ScoreContext& context, std::int64_t bin,
                      std::int64_t window_length, std::int64_t& left,
                      std::int64_t& right) {
    const hicx::CutInterval& interval =
        (*context.cut_intervals)[static_cast<std::size_t>(bin)];
    const std::int64_t left_start = std::max<std::int64_t>(0, interval.start - window_length);
    const auto left_range =
        context.bins->region_bin_range(interval.chrom, left_start, left_start + 1);
    if (!left_range.has_value()) {
        return false;
    }
    const auto size = context.chromosome_sizes.find(interval.chrom);
    if (size == context.chromosome_sizes.end()) {
        return false;
    }
    const std::int64_t right_end =
        std::min(size->second, interval.end + window_length) - 1;
    const auto right_range =
        context.bins->region_bin_range(interval.chrom, right_end, right_end);
    if (!right_range.has_value()) {
        return false;
    }
    left = left_range->first;
    right = right_range->first;
    return true;
}

// get_cut_weight(..., return_mean=True). `defined` is false where the Python
// returns None.
double cut_weight_mean(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch, bool& defined) {
    defined = true;
    if (cut < 0 || cut > context.matrix->rows()) {
        defined = false;
        return kNaN;
    }
    std::int64_t left = 0;
    std::int64_t right = 0;
    if (!bins_at_distance(context, cut, window_length, left, right)) {
        defined = false;  // the TypeError branch of get_cut_weight
        return kNaN;
    }
    extract_block(*context.matrix, left, cut, cut, right, scratch);
    if (scratch.stored == 0) {
        return 0.0;  // the Python returns the integer 0 here
    }
    // The hot reduction of the whole tool: one dense block per bin and window
    // size. hicx::simd::pairwise_sum is bit-identical to the numpy compatible
    // scalar reduction, so this stays byte-identical to the Python.
    const std::size_t count = scratch.values.size();
    return hicx::simd::pairwise_sum(scratch.values.data(), count) /
           static_cast<double>(count);
}

// get_cut_weight(..., return_mean=False): the flattened dense block.
bool cut_weight_values(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch) {
    if (cut < 0 || cut > context.matrix->rows()) {
        return false;
    }
    std::int64_t left = 0;
    std::int64_t right = 0;
    if (!bins_at_distance(context, cut, window_length, left, right)) {
        return false;
    }
    extract_block(*context.matrix, left, cut, cut, right, scratch);
    return true;
}

// --------------------------------------------------------------------------
// peakdetect (hicFindTADs.py:466-572)


void peakdetect(const std::vector<double>& y, const std::vector<std::string>& chrom,
                std::int64_t lookahead, std::vector<Peak>& maxima,
                std::vector<Peak>& minima) {
    maxima.clear();
    minima.clear();
    std::vector<bool> dump;
    const std::int64_t n = static_cast<std::int64_t>(y.size());
    const std::int64_t limit = n - lookahead;

    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    std::int64_t max_pos = 0;
    std::int64_t min_pos = 0;
    int search_for = 0;  // 0 none, 1 min, 2 max
    bool have_previous_chrom = false;
    std::string previous_chrom;

    for (std::int64_t index = 0; index < limit; ++index) {
        const double value = y[static_cast<std::size_t>(index)];
        if (!(value > -std::numeric_limits<double>::infinity() &&
              value < std::numeric_limits<double>::infinity())) {
            throw PythonError("AssertionError",
                              "Error, infinity value detected for value at position " +
                                  std::to_string(index));
        }
        if (!have_previous_chrom) {
            previous_chrom = chrom[static_cast<std::size_t>(index)];
            have_previous_chrom = true;
        }
        if (previous_chrom != chrom[static_cast<std::size_t>(index)]) {
            min_y = std::numeric_limits<double>::infinity();
            max_y = -std::numeric_limits<double>::infinity();
            search_for = 0;
        }
        previous_chrom = chrom[static_cast<std::size_t>(index)];

        if (value > max_y) {
            max_y = value;
            max_pos = index;
        }
        if (value < min_y) {
            min_y = value;
            min_pos = index;
        }

        if (value < max_y - 0.0 && max_y != std::numeric_limits<double>::infinity() &&
            search_for != 1) {
            double ahead = -std::numeric_limits<double>::infinity();
            for (std::int64_t k = index; k < index + lookahead && k < n; ++k) {
                ahead = std::max(ahead, y[static_cast<std::size_t>(k)]);
            }
            if (ahead < max_y) {
                maxima.push_back(Peak{max_pos, max_y});
                dump.push_back(true);
                max_y = value;
                min_y = value;
                min_pos = index;
                search_for = 1;
                continue;
            }
        }

        if (value > min_y + 0.0 && min_y != -std::numeric_limits<double>::infinity() &&
            search_for != 2) {
            double ahead = std::numeric_limits<double>::infinity();
            for (std::int64_t k = index; k < index + lookahead && k < n; ++k) {
                ahead = std::min(ahead, y[static_cast<std::size_t>(k)]);
            }
            if (ahead > min_y) {
                minima.push_back(Peak{min_pos, min_y});
                dump.push_back(false);
                min_y = value;
                max_y = value;
                max_pos = index;
                search_for = 2;
            }
        }
    }

    if (!dump.empty()) {
        if (dump[0]) {
            if (!maxima.empty()) {
                maxima.erase(maxima.begin());
            }
        } else if (!minima.empty()) {
            minima.erase(minima.begin());
        }
    }
}

// delta_wrt_window (hicFindTADs.py:574-616).
std::map<std::int64_t, double> delta_wrt_window(const std::vector<std::int64_t>& minima,
                                                const std::vector<double>& scores,
                                                const std::vector<std::string>& chrom,
                                                std::int64_t window_length) {
    // np.unique(chrom, return_index=True) gives the first index of every
    // distinct name; the ranges are then built from the sorted indices plus
    // len(chrom) - 1.
    std::map<std::string, std::int64_t> first_index;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(chrom.size()); ++i) {
        first_index.emplace(chrom[static_cast<std::size_t>(i)], i);
    }
    std::vector<std::int64_t> boundaries;
    boundaries.reserve(first_index.size() + 1);
    for (const auto& [name, index] : first_index) {
        boundaries.push_back(index);
    }
    boundaries.push_back(static_cast<std::int64_t>(chrom.size()) - 1);
    std::sort(boundaries.begin(), boundaries.end());

    std::map<std::int64_t, double> delta;
    for (const std::int64_t minimum : minima) {
        bool close_to_border = true;
        for (std::size_t r = 0; r + 1 < boundaries.size(); ++r) {
            const std::int64_t start = boundaries[r];
            const std::int64_t end = boundaries[r + 1];
            if (start < minimum && minimum < end) {
                if (minimum - window_length >= start && minimum + window_length < end) {
                    close_to_border = false;
                }
            }
        }
        if (close_to_border) {
            delta[minimum] = kNaN;
            continue;
        }
        // matrix_avg[i - w : i + 3] concatenated with matrix_avg[i + 4 : i + w].
        std::vector<double> window;
        const std::int64_t n = static_cast<std::int64_t>(scores.size());
        const auto append = [&](std::int64_t from, std::int64_t to) {
            from = std::max<std::int64_t>(0, from);
            to = std::min(n, to);
            for (std::int64_t k = from; k < to; ++k) {
                window.push_back(scores[static_cast<std::size_t>(k)]);
            }
        };
        append(minimum - window_length, minimum + 3);
        append(minimum + 4, minimum + window_length);
        const double mean = hicx::npy::pairwise_sum(window.data(), window.size()) /
                            static_cast<double>(window.size());
        delta[minimum] = mean - scores[static_cast<std::size_t>(minimum)];
    }
    return delta;
}


std::vector<double> delta_wrt_window_values(const std::vector<std::int64_t>& minima,
                                            const std::vector<double>& scores,
                                            const std::vector<std::string>& chrom,
                                            std::int64_t window_length) {
    std::map<std::string, std::int64_t> first_index;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(chrom.size()); ++i) {
        first_index.emplace(chrom[static_cast<std::size_t>(i)], i);
    }
    std::vector<std::int64_t> boundaries;
    boundaries.reserve(first_index.size() + 1);
    for (const auto& [name, index] : first_index) {
        boundaries.push_back(index);
    }
    boundaries.push_back(static_cast<std::int64_t>(chrom.size()) - 1);
    std::sort(boundaries.begin(), boundaries.end());

    const std::int64_t n = static_cast<std::int64_t>(scores.size());
    std::vector<double> out;
    out.reserve(minima.size());
    for (const std::int64_t minimum : minima) {
        bool close_to_border = true;
        for (std::size_t r = 0; r + 1 < boundaries.size(); ++r) {
            const std::int64_t start = boundaries[r];
            const std::int64_t end = boundaries[r + 1];
            if (start < minimum && minimum < end && minimum - window_length >= start &&
                minimum + window_length < end) {
                close_to_border = false;
            }
        }
        if (close_to_border) {
            out.push_back(kNaN);
            continue;
        }
        std::vector<double> window;
        const auto append = [&](std::int64_t from, std::int64_t to) {
            const auto [first, last] = python_slice(from, to, n);
            for (std::int64_t k = first; k < last; ++k) {
                window.push_back(scores[static_cast<std::size_t>(k)]);
            }
        };
        append(minimum - window_length, minimum + 3);
        append(minimum + 4, minimum + window_length);
        const double mean = hicx::npy::pairwise_sum(window.data(), window.size()) /
                            static_cast<double>(window.size());
        const std::int64_t at = minimum < 0 ? minimum + n : minimum;
        out.push_back(mean - scores[static_cast<std::size_t>(at)]);
    }
    return out;
}

std::vector<double> incremental_step_size_real(double min_window, double max_window,
                                               double start_step) {
    if (!std::isfinite(start_step)) {
        // k = 0: start_step * 0.0 is NaN, and int(NaN) raises.
        throw PythonError("ValueError", "cannot convert float NaN to integer");
    }
    if (std::isnan(min_window) || std::isnan(max_window) || max_window == kInf ||
        (start_step <= 0 && min_window <= max_window)) {
        throw PythonError("ValueError",
                          "get_incremental_step_size never terminates for min_win_size " +
                              hicx::npy::float_repr(min_window) + ", max_win_size " +
                              hicx::npy::float_repr(max_window) + " and start_step_len " +
                              hicx::npy::float_repr(start_step) +
                              " (HiCExplorer 3.7.6 loops forever; HiCExplorer 4 refuses)");
    }
    std::vector<double> steps;
    std::int64_t step = -1;
    while (true) {
        ++step;
        const double grown = start_step * std::pow(static_cast<double>(step), 1.5);
        const double increment = min_window + std::trunc(grown);
        if (step > 1 && !steps.empty() && increment == steps.back()) {
            continue;
        }
        if (increment > max_window) {
            break;
        }
        steps.push_back(increment);
    }
    return steps;
}

CsrMatrix diagflat_zero(const CsrMatrix& matrix) {
    const std::string& dtype = matrix.dtype();
    std::string result_dtype = "int64";
    if (dtype == "float16" || dtype == "float32" || dtype == "float64" || dtype == "uint64") {
        result_dtype = "float64";
    }
    const bool integral = result_dtype == "int64";
    const std::int64_t n = matrix.rows();
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    indices.reserve(matrix.stored_nnz());
    values.reserve(matrix.stored_nnz());
    std::vector<std::pair<std::int32_t, double>> row_entries;
    for (std::int64_t row = 0; row < n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        row_entries.clear();
        for (std::size_t k = begin; k < end; ++k) {
            row_entries.emplace_back(matrix.indices()[k], matrix.data()[k]);
        }
        // scipy's binary operations sum duplicates and sort the columns.
        std::stable_sort(row_entries.begin(), row_entries.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        std::size_t k = 0;
        while (k < row_entries.size()) {
            const std::int32_t column = row_entries[k].first;
            double value = row_entries[k].second;
            ++k;
            while (k < row_entries.size() && row_entries[k].first == column) {
                value += row_entries[k].second;
                ++k;
            }
            if (column == row) {
                // (x - x) + 0: zero for a finite x, NaN for NaN and infinity.
                value = std::isfinite(value) || integral ? 0.0 : kNaN;
            }
            if (value == 0.0) {
                continue;
            }
            indices.push_back(column);
            values.push_back(value);
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(values.size());
    }
    CsrMatrix out(n, matrix.cols(), std::move(indptr), std::move(indices), std::move(values),
                  result_dtype);
    out.set_symmetry(Symmetry::Full);
    return out;
}

std::vector<double> row_means(const std::vector<double>& values, std::size_t rows,
                              std::size_t columns) {
    std::vector<double> means(rows, 0.0);
    for (std::size_t r = 0; r < rows; ++r) {
        means[r] = hicx::simd::pairwise_sum(values.data() + r * columns, columns) /
                   static_cast<double>(columns);
    }
    return means;
}

PeakResult peakdetect_general(const std::vector<double>& y,
                              const std::vector<std::string>* chrom, std::int64_t lookahead,
                              double delta, bool y_is_array, bool assertions,
                              const std::string& y_type_name) {
    PeakResult result;
    std::vector<bool> dump;
    const std::int64_t n = static_cast<std::int64_t>(y.size());
    const std::int64_t limit = python_slice(0, -lookahead, n).second;

    double min_y = kInf;
    double max_y = -kInf;
    std::int64_t max_pos = -1;
    std::int64_t min_pos = -1;
    std::int64_t max_at = -1;  // the index whose value max_y holds
    std::int64_t min_at = -1;
    int search_for = 0;  // 0 None, 1 'min', 2 'max'
    bool have_previous = false;
    std::string previous;

    const auto ahead = [&](std::int64_t index, bool maximum) {
        if (!y_is_array) {
            throw PythonError("AttributeError", "'" + y_type_name + "' object has no attribute '" +
                                                    std::string(maximum ? "max" : "min") + "'");
        }
        const auto [first, last] = python_slice(index, index + lookahead, n);
        if (first == last) {
            throw PythonError("ValueError", std::string("zero-size array to reduction operation ") +
                                                (maximum ? "maximum" : "minimum") +
                                                " which has no identity");
        }
        double value = y[static_cast<std::size_t>(first)];
        for (std::int64_t k = first; k < last; ++k) {
            const double v = y[static_cast<std::size_t>(k)];
            if (std::isnan(v)) {
                return kNaN;  // np.max and np.min propagate NaN
            }
            value = maximum ? std::max(value, v) : std::min(value, v);
        }
        return value;
    };

    for (std::int64_t index = 0; index < limit; ++index) {
        const double value = y[static_cast<std::size_t>(index)];
        if (assertions && !(-kInf < value && value < kInf)) {
            throw PythonError("AssertionError",
                              "Error, infinity value detected for value at position " +
                                  std::to_string(index));
        }
        if (chrom == nullptr) {
            throw PythonError("TypeError", "'NoneType' object is not subscriptable");
        }
        if (index >= static_cast<std::int64_t>(chrom->size())) {
            throw PythonError("IndexError", "index " + std::to_string(index) +
                                                " is out of bounds for axis 0 with size " +
                                                std::to_string(chrom->size()));
        }
        const std::string& name = (*chrom)[static_cast<std::size_t>(index)];
        if (!have_previous) {
            previous = name;
            have_previous = true;
        }
        if (previous != name) {
            min_y = kInf;
            max_y = -kInf;
            max_pos = -1;
            min_pos = -1;
            search_for = 0;
        }
        previous = name;

        if (value > max_y) {
            max_y = value;
            max_at = index;
            max_pos = index;
        }
        if (value < min_y) {
            min_y = value;
            min_at = index;
            min_pos = index;
        }

        if (value < max_y - delta && max_y != kInf && search_for != 1) {
            if (ahead(index, true) < max_y) {
                result.maxima.push_back(PeakAt{max_pos, max_at});
                dump.push_back(true);
                max_y = value;
                max_at = index;
                min_y = value;
                min_at = index;
                min_pos = index;
                search_for = 1;
                continue;
            }
        }

        if (value > min_y + delta && min_y != -kInf && search_for != 2) {
            if (ahead(index, false) > min_y) {
                result.minima.push_back(PeakAt{min_pos, min_at});
                dump.push_back(false);
                min_y = value;
                min_at = index;
                max_y = value;
                max_at = index;
                max_pos = index;
                search_for = 2;
            }
        }
    }

    if (!dump.empty()) {
        if (dump[0]) {
            result.maxima.erase(result.maxima.begin());
        } else {
            result.minima.erase(result.minima.begin());
        }
    }
    return result;
}

std::optional<double> fdr_cutoff(std::vector<double> pvalues, double q) {
    std::sort(pvalues.begin(), pvalues.end());
    std::optional<double> largest;
    const double n = static_cast<double>(pvalues.size());
    for (std::size_t i = 0; i < pvalues.size(); ++i) {
        const double p = pvalues[i];
        if (p <= q * static_cast<double>(i + 1) / n) {
            if (p >= (largest.has_value() ? *largest : 0.0)) {
                largest = p;
            }
        }
    }
    return largest;
}

double numpy_median(std::vector<double> values) {
    return hicx::utilities::median(values);
}

double numpy_percentile(std::vector<double> values, double q) {
    if (values.empty()) {
        return kNaN;
    }
    for (const double v : values) {
        if (std::isnan(v)) {
            return kNaN;
        }
    }
    std::sort(values.begin(), values.end());
    const double virtual_index = static_cast<double>(values.size() - 1) * (q / 100.0);
    const double previous = std::floor(virtual_index);
    std::size_t below = static_cast<std::size_t>(std::max(0.0, previous));
    below = std::min(below, values.size() - 1);
    const std::size_t above = std::min(below + 1, values.size() - 1);
    const double gamma = virtual_index - previous;
    const double a = values[below];
    const double b = values[above];
    const double diff = b - a;
    // numpy's _lerp
    return gamma >= 0.5 ? b - diff * (1.0 - gamma) : a + diff * gamma;
}

std::vector<std::pair<std::size_t, std::size_t>> numpy_split_bounds(
    std::size_t length, const std::vector<std::int64_t>& indices) {
    std::vector<std::int64_t> points;
    points.reserve(indices.size() + 2);
    points.push_back(0);
    points.insert(points.end(), indices.begin(), indices.end());
    points.push_back(static_cast<std::int64_t>(length));
    std::vector<std::pair<std::size_t, std::size_t>> out;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        const auto [first, last] =
            python_slice(points[i], points[i + 1], static_cast<std::int64_t>(length));
        out.emplace_back(static_cast<std::size_t>(first), static_cast<std::size_t>(last));
    }
    return out;
}

namespace {

std::vector<std::int64_t> argsort_as(const std::vector<double>& values, bool integer) {
    if (integer) {
        std::vector<std::int64_t> ints(values.begin(), values.end());
        return hicx::npy::argsort(ints);
    }
    return hicx::npy::argsort(values);
}

}  // namespace

ClusterResult cluster_chromosome(const ClusterInput& input, const std::vector<double>& cutoffs) {
    ClusterResult result;
    const std::size_t count = input.start.size();
    const std::size_t domain_count = count == 0 ? 0 : count - 1;
    if (static_cast<std::int64_t>(domain_count) != static_cast<std::int64_t>(input.value.size()) - 1) {
        throw PythonError("AssertionError", "error");
    }
    std::vector<std::set<std::size_t>> clusters;
    std::vector<std::int64_t> cluster_x;
    for (std::size_t d = 0; d < domain_count; ++d) {
        clusters.push_back({d});
        const double d_start = input.start[d];
        const double d_end = input.start[d + 1];
        cluster_x.push_back(static_cast<std::int64_t>(d_start + (d_end - d_start) / 2.0));
    }
    // values = value[1:-1]
    std::vector<double> values;
    for (std::size_t i = 1; i + 1 < input.value.size(); ++i) {
        values.push_back(input.value[i]);
    }
    std::vector<std::int64_t> order = argsort_as(values, input.value_integer);
    std::reverse(order.begin(), order.end());

    std::size_t cutoff_at = 0;
    const auto find_in_clusters = [&clusters](std::size_t id) -> std::int64_t {
        for (std::size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].count(id) != 0) {
                return static_cast<std::int64_t>(i);
            }
        }
        throw PythonError("TypeError", "list indices must be integers or slices, not NoneType");
    };
    for (std::size_t idx = 0; idx < order.size(); ++idx) {
        const std::size_t order_idx = static_cast<std::size_t>(order[idx]);
        if (cutoff_at < cutoffs.size() && idx + 1 < order.size() &&
            values[order_idx] >= cutoffs[cutoff_at] &&
            cutoffs[cutoff_at] > values[static_cast<std::size_t>(order[idx + 1])]) {
            ++cutoff_at;
            std::vector<ClusterRegion> regions;
            std::vector<double> starts;
            for (const std::set<std::size_t>& set : clusters) {
                if (set.empty()) {
                    continue;
                }
                regions.push_back(ClusterRegion{*set.begin(), *set.rbegin()});
                starts.push_back(input.start[*set.begin()]);
            }
            const std::vector<std::int64_t> by_start = argsort_as(starts, input.start_integer);
            std::vector<ClusterRegion> ordered;
            ordered.reserve(regions.size());
            for (const std::int64_t k : by_start) {
                ordered.push_back(regions[static_cast<std::size_t>(k)]);
            }
            result.levels.push_back(std::move(ordered));
        }
        const std::int64_t left = find_in_clusters(order_idx);
        const std::int64_t right = find_in_clusters(order_idx + 1);
        const std::size_t l = static_cast<std::size_t>(left);
        const std::size_t r = static_cast<std::size_t>(right);
        result.merges.push_back(ClusterMerge{
            left, right, order_idx + 1,
            static_cast<std::int64_t>(clusters[l].size() + clusters[r].size()), cluster_x[l],
            cluster_x[r]});
        const std::int64_t gen_dist = static_cast<std::int64_t>(
            static_cast<double>(std::llabs(cluster_x[l] - cluster_x[r])) / 2.0);
        cluster_x.push_back(std::min(cluster_x[l], cluster_x[r]) + gen_dist);
        std::set<std::size_t> merged = clusters[l];
        merged.insert(clusters[r].begin(), clusters[r].end());
        clusters.push_back(std::move(merged));
        clusters[l].clear();
        clusters[r].clear();
    }
    return result;
}

std::string format_fixed(double value, int decimals) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return buffer;
}

}  // namespace hicx::find_tads
