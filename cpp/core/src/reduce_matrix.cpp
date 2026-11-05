#include "hicx/reduce_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "hicx/adjust_ops.hpp"
#include "hicx/numpy_compat.hpp"

namespace hicx {

namespace {

// numpy's astype, which is a C cast and therefore truncates toward zero on the
// way to an integer type. reduce_matrix sums with np.bincount, which always
// produces float64, and then hands the result to coo_matrix(dtype=ma.dtype).
double cast_to_dtype(double value, const std::string& dtype) {
    if (!std::isfinite(value)) {
        // A cast of NaN or an infinity to an integer is undefined in C++ and
        // implementation defined in numpy. Leave it alone rather than invent a
        // value; the only matrices that reach this are float ones anyway.
        return value;
    }
    if (dtype == "float32") {
        return static_cast<double>(static_cast<float>(value));
    }
    if (dtype_from_name(dtype) != DType::Integer) {
        return value;
    }
    // Truncate toward zero, then narrow. A double outside the range of the
    // integer type is undefined behaviour to cast directly, so the value is
    // clamped first; numpy wraps instead, but a merged count that overflows
    // its own dtype is a broken input either way and no matrix in the corpus
    // comes close.
    const double truncated = std::trunc(value);
    const double low = -9.2233720368547758e18;
    const double high = 9.2233720368547758e18;
    if (truncated <= low || truncated >= high) {
        return truncated;
    }
    const std::int64_t wide = static_cast<std::int64_t>(truncated);
    if (dtype == "int8") {
        return static_cast<double>(static_cast<std::int8_t>(wide));
    }
    if (dtype == "int16") {
        return static_cast<double>(static_cast<std::int16_t>(wide));
    }
    if (dtype == "int32") {
        return static_cast<double>(static_cast<std::int32_t>(wide));
    }
    if (dtype == "uint8") {
        return static_cast<double>(static_cast<std::uint8_t>(wide));
    }
    if (dtype == "uint16") {
        return static_cast<double>(static_cast<std::uint16_t>(wide));
    }
    if (dtype == "uint32") {
        return static_cast<double>(static_cast<std::uint32_t>(wide));
    }
    if (dtype == "uint64") {
        return static_cast<double>(static_cast<std::uint64_t>(wide));
    }
    return static_cast<double>(wide);
}

}  // namespace

CsrMatrix reduce_matrix(const CsrMatrix& matrix,
                        const std::vector<std::vector<std::int64_t>>& bins_to_merge,
                        bool use_triu, bool diagonal) {
    const std::int64_t rows = matrix.rows();
    const std::int64_t groups = static_cast<std::int64_t>(bins_to_merge.size());
    if (groups == rows) {
        // reduceMatrix.py:154-155 short circuits here and hands back the input
        // itself, neither triangle selected nor symmetrised.
        return matrix;
    }

    // map_ in the Python: the group every input index belongs to, -1 for the
    // indices that are dropped.
    std::vector<std::int64_t> group_of(static_cast<std::size_t>(rows), -1);
    for (std::size_t k = 0; k < bins_to_merge.size(); ++k) {
        for (const std::int64_t index : bins_to_merge[k]) {
            if (index < 0 || index >= rows) {
                throw std::out_of_range("bins_to_merge refers to row " +
                                        std::to_string(index) + " of a matrix with " +
                                        std::to_string(rows) + " rows");
            }
            group_of[static_cast<std::size_t>(index)] = static_cast<std::int64_t>(k);
        }
    }

    // The inverse of map_: the input rows of every group, ascending. Derived
    // from group_of rather than from bins_to_merge, so that a row named by two
    // groups belongs to the later one, which is what the Python's overwriting
    // assignment produces.
    std::vector<std::int64_t> member_start(static_cast<std::size_t>(groups) + 1, 0);
    for (const std::int64_t group : group_of) {
        if (group >= 0) {
            ++member_start[static_cast<std::size_t>(group) + 1];
        }
    }
    for (std::size_t i = 1; i < member_start.size(); ++i) {
        member_start[i] += member_start[i - 1];
    }
    std::vector<std::int64_t> member(static_cast<std::size_t>(member_start.back()));
    {
        std::vector<std::int64_t> cursor(member_start.begin(), member_start.end() - 1);
        for (std::int64_t row = 0; row < rows; ++row) {
            const std::int64_t group = group_of[static_cast<std::size_t>(row)];
            if (group >= 0) {
                member[static_cast<std::size_t>(cursor[static_cast<std::size_t>(group)]++)] =
                    row;
            }
        }
    }

    // np.unique(new_row + 1j * new_col, return_inverse=True) followed by
    // np.bincount(inverse, weights=data) gives one accumulator per distinct
    // (row, col) of the reduced matrix, filled in the order the entries of
    // triu(matrix, k=0, format='coo') appear, which is the row major order of
    // the CSR. That order is reproduced here without staging the entries: one
    // output row is fed only by the input rows of its own group, and those are
    // visited ascending, so a sparse accumulator over the output columns sums
    // every cell in exactly the order np.bincount would. Nothing that scales
    // with the entry count is allocated beyond the result itself.
    //
    // Explicit zeros take part, as scipy's triu keeps them, and the result is
    // pruned only at the very end.
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    std::vector<double> accumulator(static_cast<std::size_t>(groups), 0.0);
    std::vector<char> touched(static_cast<std::size_t>(groups), 0);
    std::vector<std::int32_t> touched_columns;

    std::vector<std::int64_t> result_indptr(static_cast<std::size_t>(groups) + 1, 0);
    std::vector<std::int32_t> result_indices;
    std::vector<double> result_values;
    bool result_is_upper = true;

    for (std::int64_t group = 0; group < groups; ++group) {
        touched_columns.clear();
        const std::size_t first = static_cast<std::size_t>(member_start[static_cast<std::size_t>(group)]);
        const std::size_t last = static_cast<std::size_t>(member_start[static_cast<std::size_t>(group) + 1]);
        for (std::size_t m = first; m < last; ++m) {
            const std::int64_t row = member[m];
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                const std::int64_t col = indices[k];
                if (use_triu && col < row) {
                    continue;
                }
                const std::int64_t new_col = group_of[static_cast<std::size_t>(col)];
                if (new_col < 0) {
                    continue;
                }
                const std::size_t slot = static_cast<std::size_t>(new_col);
                if (touched[slot] == 0) {
                    touched[slot] = 1;
                    touched_columns.push_back(static_cast<std::int32_t>(new_col));
                }
                accumulator[slot] += values[k];
            }
        }
        std::sort(touched_columns.begin(), touched_columns.end());
        for (const std::int32_t column : touched_columns) {
            const std::size_t slot = static_cast<std::size_t>(column);
            result_indices.push_back(column);
            result_values.push_back(cast_to_dtype(accumulator[slot], matrix.dtype()));
            accumulator[slot] = 0.0;
            touched[slot] = 0;
            if (column < group) {
                result_is_upper = false;
            }
        }
        result_indptr[static_cast<std::size_t>(group) + 1] =
            static_cast<std::int64_t>(result_indices.size());
    }

    CsrMatrix result(groups, groups, std::move(result_indptr), std::move(result_indices),
                     std::move(result_values), matrix.dtype());

    if (!use_triu) {
        // result, without the symmetrisation.
        result.eliminate_zeros();
        return result;
    }

    if (result_is_upper) {
        // result + result.T - diagmatrix over an upper triangular result is
        // exactly what an upper triangle with symmetric access represents, so
        // the mirror is never built. Without pDiagonal the Python subtracts
        // twice the diagonal, which leaves every diagonal entry at exactly
        // zero, and eliminate_zeros then prunes it.
        if (!diagonal) {
            std::vector<double>& stored = result.mutable_data();
            for (std::int64_t row = 0; row < groups; ++row) {
                const std::size_t begin =
                    static_cast<std::size_t>(result.indptr()[static_cast<std::size_t>(row)]);
                const std::size_t end = static_cast<std::size_t>(
                    result.indptr()[static_cast<std::size_t>(row) + 1]);
                for (std::size_t k = begin; k < end; ++k) {
                    if (result.indices()[k] == row) {
                        stored[k] = 0.0;
                    }
                }
            }
        }
        result.set_symmetry(Symmetry::UpperTriangle);
        result.eliminate_zeros();
        return result;
    }

    // A group list that is not increasing can put an entry below the diagonal,
    // and then result + result.T genuinely adds two different cells. That path
    // does need both triangles, and no caller in HiCExplorer reaches it: every
    // grouping the tools build is a partition into ascending ranges.
    std::vector<std::int32_t> out_rows;
    std::vector<std::int32_t> out_cols;
    std::vector<double> out_values;
    out_rows.reserve(result.stored_nnz() * 2);
    out_cols.reserve(result.stored_nnz() * 2);
    out_values.reserve(result.stored_nnz() * 2);
    const auto emit = [&](std::int64_t row, std::int64_t col, double value) {
        out_rows.push_back(static_cast<std::int32_t>(row));
        out_cols.push_back(static_cast<std::int32_t>(col));
        out_values.push_back(value);
    };
    // The entries of `result` first and its transpose afterwards, so that a
    // cell fed from both sides accumulates in the order csr_plus_csr uses.
    for (std::int64_t row = 0; row < groups; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(result.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(result.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = result.indices()[k];
            if (col == row && !diagonal) {
                continue;
            }
            emit(row, col, result.data()[k]);
        }
    }
    for (std::int64_t row = 0; row < groups; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(result.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(result.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = result.indices()[k];
            if (col == row) {
                continue;
            }
            emit(col, row, result.data()[k]);
        }
    }

    CsrMatrix symmetric = CsrMatrix::from_coo(groups, groups, out_rows, out_cols,
                                              std::move(out_values), matrix.dtype());
    symmetric.eliminate_zeros();
    return symmetric;
}

BinMergePlan plan_bin_merge(const std::vector<CutInterval>& intervals,
                            std::int64_t num_bins) {
    if (intervals.empty()) {
        throw std::invalid_argument("cannot merge the bins of an empty bin table");
    }
    if (num_bins < 1) {
        throw std::invalid_argument("--numBins must be at least 1");
    }

    BinMergePlan plan;
    const auto coverage_mean = [&](std::size_t first, std::size_t last) {
        std::vector<double> window;
        window.reserve(last - first);
        for (std::size_t i = first; i < last; ++i) {
            window.push_back(intervals[i].extra);
        }
        if (window.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return npy::pairwise_sum(window) / static_cast<double>(window.size());
    };

    std::string previous = intervals[0].chrom;
    std::size_t start_index = 0;
    std::int64_t new_start = intervals[0].start;
    std::int64_t count = 0;
    std::size_t index = 0;
    for (index = 0; index < intervals.size(); ++index) {
        const std::string& chrom = intervals[index].chrom;
        if ((count > 0 && count % num_bins == 0) || chrom != previous) {
            // A trailing piece shorter than half a merged bin is dropped, with
            // its rows and columns, because it never enters bins_to_merge.
            if (static_cast<double>(count) >=
                static_cast<double>(num_bins) / 2.0) {
                std::vector<std::int64_t> group;
                group.reserve(index - start_index);
                for (std::size_t i = start_index; i < index; ++i) {
                    group.push_back(static_cast<std::int64_t>(i));
                }
                plan.intervals.push_back(CutInterval{intervals[start_index].chrom,
                                                     new_start, intervals[index - 1].end,
                                                     coverage_mean(start_index, index),
                                                     ""});
                plan.bins_to_merge.push_back(std::move(group));
            }
            start_index = index;
            new_start = intervals[index].start;
            count = 0;
        }
        previous = chrom;
        ++count;
    }
    // The final group is appended unconditionally, so it survives even when it
    // is shorter than half a merged bin.
    index = intervals.size() - 1;
    std::vector<std::int64_t> group;
    group.reserve(intervals.size() - start_index);
    for (std::size_t i = start_index; i <= index; ++i) {
        group.push_back(static_cast<std::int64_t>(i));
    }
    plan.intervals.push_back(CutInterval{intervals[index].chrom, new_start,
                                         intervals[index].end,
                                         coverage_mean(start_index, intervals.size()),
                                         ""});
    plan.bins_to_merge.push_back(std::move(group));
    return plan;
}

BinMergePlan plan_tad_merge(const std::vector<CutInterval>& intervals,
                            const std::vector<std::int64_t>& boundaries) {
    BinMergePlan plan;
    if (intervals.empty()) {
        return plan;
    }
    std::vector<char> is_boundary(intervals.size(), 0);
    for (const std::int64_t bin : boundaries) {
        if (bin >= 0 && bin < static_cast<std::int64_t>(intervals.size())) {
            is_boundary[static_cast<std::size_t>(bin)] = 1;
        }
    }

    const auto coverage_mean = [&](std::size_t first, std::size_t last) {
        std::vector<double> window;
        window.reserve(last - first);
        for (std::size_t i = first; i < last; ++i) {
            window.push_back(intervals[i].extra);
        }
        if (window.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return npy::pairwise_sum(window) / static_cast<double>(window.size());
    };

    std::string previous = intervals[0].chrom;
    std::size_t start_index = 0;
    std::int64_t new_start = intervals[0].start;
    std::int64_t count = 0;
    std::size_t index = 0;
    for (index = 0; index < intervals.size(); ++index) {
        const std::string& chrom = intervals[index].chrom;
        if ((count > 0 && is_boundary[index] != 0) || chrom != previous) {
            std::vector<std::int64_t> group;
            group.reserve(index - start_index);
            for (std::size_t i = start_index; i < index; ++i) {
                group.push_back(static_cast<std::int64_t>(i));
            }
            plan.intervals.push_back(CutInterval{intervals[start_index].chrom, new_start,
                                                 intervals[index - 1].end,
                                                 coverage_mean(start_index, index), ""});
            plan.bins_to_merge.push_back(std::move(group));
            start_index = index;
            new_start = intervals[index].start;
            count = 0;
        }
        previous = chrom;
        ++count;
    }
    // The tail group is appended only when the loop produced something, which
    // is the len(bins_to_merge) > 0 guard at hicMergeTADbins.py:78.
    if (plan.bins_to_merge.empty()) {
        plan.intervals.clear();
        return plan;
    }
    index = intervals.size() - 1;
    std::vector<std::int64_t> group;
    group.reserve(intervals.size() - start_index);
    for (std::size_t i = start_index; i <= index; ++i) {
        group.push_back(static_cast<std::int64_t>(i));
    }
    plan.intervals.push_back(CutInterval{intervals[index].chrom, new_start,
                                         intervals[index].end,
                                         coverage_mean(start_index, intervals.size()), ""});
    plan.bins_to_merge.push_back(std::move(group));
    return plan;
}

CsrMatrix running_window(const CsrMatrix& matrix, std::int64_t num_bins) {
    if (num_bins % 2 == 0) {
        throw std::invalid_argument("num_bins has to be an odd number");
    }
    const std::int64_t size = matrix.rows();
    const std::int64_t half = (num_bins - 1) / 2;

    // Every entry of triu(matrix, k=0) is spread over the window around it and
    // only the cells at or above the diagonal are kept, because folding to the
    // upper triangle before the duplicates are summed is equivalent to the
    // Python's order: a cell's value depends only on the entries that land on
    // it.
    //
    // The result is produced one output row at a time rather than by staging
    // every shifted entry and sorting them. Output row r is fed only by the
    // input rows r - half .. r + half, so a sparse accumulator over the output
    // columns is enough, and the resident cost is the input, the output and
    // two arrays of one entry per bin. Staging the shifted entries would cost
    // num_bins^2 / 2 times the entry count: on Li_et_al_2015.h5 with
    // --numBins 3 that was 473 MB against a 20 MB matrix.
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();

    std::vector<double> accumulator(static_cast<std::size_t>(size), 0.0);
    std::vector<char> touched(static_cast<std::size_t>(size), 0);
    std::vector<std::int32_t> touched_columns;

    std::vector<std::int64_t> out_indptr(static_cast<std::size_t>(size) + 1, 0);
    std::vector<std::int32_t> out_indices;
    std::vector<double> out_data;

    // The row is walked twice: once to count the cells it will hold and once
    // to fill them. That makes the output arrays exactly the right size, and a
    // vector that grows by doubling would otherwise hold its old and its new
    // buffer at the same time, which on this input is a third of the peak.
    // `accumulate` false is the counting pass.
    const auto walk_row = [&](std::int64_t out_row, bool accumulate) {
        touched_columns.clear();
        const std::int64_t first_row = std::max<std::int64_t>(0, out_row - half);
        const std::int64_t last_row = std::min<std::int64_t>(size - 1, out_row + half);
        for (std::int64_t row = first_row; row <= last_row; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                const std::int64_t column = indices[k];
                if (column < row) {
                    continue;  // triu(matrix, k=0)
                }
                const std::int64_t first_column =
                    std::max<std::int64_t>(out_row, column - half);
                const std::int64_t last_column =
                    std::min<std::int64_t>(size - 1, column + half);
                for (std::int64_t out_column = first_column; out_column <= last_column;
                     ++out_column) {
                    const std::size_t slot = static_cast<std::size_t>(out_column);
                    if (touched[slot] == 0) {
                        touched[slot] = 1;
                        touched_columns.push_back(static_cast<std::int32_t>(out_column));
                    }
                    if (accumulate) {
                        accumulator[slot] += values[k];
                    }
                }
            }
        }
    };

    for (std::int64_t out_row = 0; out_row < size; ++out_row) {
        walk_row(out_row, false);
        out_indptr[static_cast<std::size_t>(out_row) + 1] =
            out_indptr[static_cast<std::size_t>(out_row)] +
            static_cast<std::int64_t>(touched_columns.size());
        for (const std::int32_t column : touched_columns) {
            touched[static_cast<std::size_t>(column)] = 0;
        }
    }
    out_indices.reserve(static_cast<std::size_t>(out_indptr.back()));
    out_data.reserve(static_cast<std::size_t>(out_indptr.back()));

    for (std::int64_t out_row = 0; out_row < size; ++out_row) {
        walk_row(out_row, true);
        std::sort(touched_columns.begin(), touched_columns.end());
        for (const std::int32_t column : touched_columns) {
            const std::size_t slot = static_cast<std::size_t>(column);
            out_indices.push_back(column);
            out_data.push_back(cast_to_dtype(accumulator[slot], matrix.dtype()));
            accumulator[slot] = 0.0;
            touched[slot] = 0;
        }
    }

    CsrMatrix result(size, size, std::move(out_indptr), std::move(out_indices),
                     std::move(out_data), matrix.dtype());
    // R + R.T - diag(R) is exactly what an upper triangle with symmetric access
    // represents, so the mirror is never materialised.
    result.set_symmetry(Symmetry::UpperTriangle);
    result.eliminate_zeros();
    return result;
}

MatrixData merge_bins(const MatrixData& input, std::int64_t num_bins) {
    BinMergePlan plan = plan_bin_merge(input.cut_intervals, num_bins);

    MatrixData merged;
    merged.matrix = reduce_matrix(input.matrix, plan.bins_to_merge, true, true);
    merged.matrix.eliminate_zeros();
    merged.cut_intervals = std::move(plan.intervals);
    merged.correction_factors = input.correction_factors;
    merged.distance_counts = input.distance_counts;
    merged.correction_factors_are_column = input.correction_factors_are_column;

    if (static_cast<std::int64_t>(merged.cut_intervals.size()) !=
        merged.matrix.rows()) {
        throw std::runtime_error("merged bin table of " +
                                 std::to_string(merged.cut_intervals.size()) +
                                 " entries does not match the merged matrix of " +
                                 std::to_string(merged.matrix.rows()) + " rows");
    }

    // hic.nan_bins = np.flatnonzero(hic.matrix.sum(0).A == 0): the columns
    // whose sum is exactly zero, which for a merged matrix means the empty
    // ones. The traversal is row major, as scipy's sum over axis 0 is.
    //
    // empty_column_bins rather than a local loop over the stored entries,
    // because reduce_matrix hands the input straight back when the group count
    // equals the row count (reduceMatrix.py:154-155), and that input may still
    // be held as an upper triangle. --numBins 1 is exactly that case.
    merged.nan_bins = empty_column_bins(merged.matrix);
    return merged;
}

BinMergePlan plan_bin_merge_genome(const std::vector<CutInterval>& intervals,
                                   std::int64_t num_bins,
                                   const ChromosomeLengths& chromosome_lengths,
                                   std::int64_t resolution) {
    if (intervals.empty()) {
        throw std::invalid_argument("cannot merge the bins of an empty bin table");
    }
    if (num_bins < 1) {
        throw std::invalid_argument("--numBins must be at least 1");
    }
    if (resolution < 1) {
        throw std::invalid_argument(
            "cannot merge by chromosome sizes: the matrix resolution is not positive");
    }
    if (chromosome_lengths.empty()) {
        throw std::invalid_argument("the chromosome sizes file is empty");
    }

    // The full genome bin table, independent of which of its bins `intervals`
    // actually holds: every chromosome tiled from 0 in steps of `resolution`,
    // exactly as read_two_dimensional_text builds it.
    std::vector<CutInterval> genome;
    for (const auto& [chrom, length] : chromosome_lengths) {
        if (length < 1) {
            throw std::invalid_argument("chromosome '" + chrom +
                                        "' has a non-positive length in the chromosome "
                                        "sizes file");
        }
        for (std::int64_t start = 0; start < length; start += resolution) {
            genome.push_back(
                CutInterval{chrom, start, std::min(length, start + resolution), 1.0, ""});
        }
    }
    if (genome.empty()) {
        throw std::invalid_argument(
            "the chromosome sizes file produced no bins at this resolution");
    }

    // (chromosome, start) -> row of `intervals`, so a bin the genome layout
    // expects can be looked up in whatever the input actually has. A repeated
    // (chromosome, start) in the input, which a well formed bin table never
    // has, keeps the first occurrence.
    std::unordered_map<std::string, std::unordered_map<std::int64_t, std::size_t>> present;
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        present[intervals[i].chrom].emplace(intervals[i].start, i);
    }

    const auto coverage_mean = [](const std::vector<double>& window) {
        if (window.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return npy::pairwise_sum(window) / static_cast<double>(window.size());
    };

    BinMergePlan plan;
    const auto flush_group = [&](std::size_t first, std::size_t last) {
        std::vector<std::int64_t> group;
        std::vector<double> coverage_values;
        for (std::size_t i = first; i < last; ++i) {
            const auto chrom_it = present.find(genome[i].chrom);
            if (chrom_it == present.end()) {
                continue;
            }
            const auto bin_it = chrom_it->second.find(genome[i].start);
            if (bin_it == chrom_it->second.end()) {
                continue;  // The genome layout has this bin, the input does not.
            }
            group.push_back(static_cast<std::int64_t>(bin_it->second));
            coverage_values.push_back(intervals[bin_it->second].extra);
        }
        plan.intervals.push_back(CutInterval{genome[first].chrom, genome[first].start,
                                             genome[last - 1].end,
                                             coverage_mean(coverage_values), ""});
        plan.bins_to_merge.push_back(std::move(group));
    };

    std::string previous = genome[0].chrom;
    std::size_t start_index = 0;
    std::int64_t count = 0;
    for (std::size_t index = 0; index < genome.size(); ++index) {
        const std::string& chrom = genome[index].chrom;
        if ((count > 0 && count % num_bins == 0) || chrom != previous) {
            // Unlike plan_bin_merge, a short trailing group of the genome
            // layout is kept rather than dropped: it is derived from the true
            // chromosome length, not from whatever the input happens to
            // contain, so there is no reason to reproduce the numBins/2 quirk
            // here. Keeping every group is also what makes the layout only
            // depend on the genome and the resolution, never on the count of
            // bins an input happens to have at its tail.
            flush_group(start_index, index);
            start_index = index;
            count = 0;
        }
        previous = chrom;
        ++count;
    }
    flush_group(start_index, genome.size());
    return plan;
}

MatrixData merge_bins_genome(const MatrixData& input, std::int64_t num_bins,
                             const ChromosomeLengths& chromosome_lengths,
                             std::int64_t resolution) {
    BinMergePlan plan =
        plan_bin_merge_genome(input.cut_intervals, num_bins, chromosome_lengths, resolution);

    MatrixData merged;
    merged.matrix = reduce_matrix(input.matrix, plan.bins_to_merge, true, true);
    merged.matrix.eliminate_zeros();
    merged.cut_intervals = std::move(plan.intervals);
    merged.correction_factors = input.correction_factors;
    merged.distance_counts = input.distance_counts;
    merged.correction_factors_are_column = input.correction_factors_are_column;

    if (static_cast<std::int64_t>(merged.cut_intervals.size()) != merged.matrix.rows()) {
        throw std::runtime_error("merged bin table of " +
                                 std::to_string(merged.cut_intervals.size()) +
                                 " entries does not match the merged matrix of " +
                                 std::to_string(merged.matrix.rows()) + " rows");
    }

    merged.nan_bins = empty_column_bins(merged.matrix);
    return merged;
}

}  // namespace hicx
