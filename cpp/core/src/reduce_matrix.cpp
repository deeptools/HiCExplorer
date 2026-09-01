#include "hicx/reduce_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

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

std::int64_t pack(std::int64_t row, std::int64_t col) {
    return (row << 32) | static_cast<std::int64_t>(static_cast<std::uint32_t>(col));
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

    // The kept entries in the order triu(matrix, k=0, format='coo') produces,
    // which is the row major order of the CSR. Explicit zeros are kept, as
    // scipy's triu keeps them, because np.bincount sums them and the result is
    // pruned only at the very end.
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& values = matrix.data();
    std::vector<std::int64_t> keys;
    std::vector<double> kept;
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = indices[k];
            if (use_triu && col < row) {
                continue;
            }
            const std::int64_t new_row = group_of[static_cast<std::size_t>(row)];
            const std::int64_t new_col = group_of[static_cast<std::size_t>(col)];
            if (new_row < 0 || new_col < 0) {
                continue;
            }
            keys.push_back(pack(new_row, new_col));
            kept.push_back(values[k]);
        }
    }

    // np.unique(new_row + 1j * new_col, return_inverse=True) followed by
    // np.bincount(inverse, weights=data): one accumulator per distinct
    // (row, col), filled in the order the entries appear. Sorting an index
    // permutation rather than the data keeps that order available.
    std::vector<std::size_t> order(keys.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });

    std::vector<std::int64_t> group_id(keys.size(), 0);
    std::vector<std::int64_t> group_key;
    for (std::size_t position = 0; position < order.size();) {
        const std::int64_t key = keys[order[position]];
        const std::int64_t id = static_cast<std::int64_t>(group_key.size());
        group_key.push_back(key);
        while (position < order.size() && keys[order[position]] == key) {
            group_id[order[position]] = id;
            ++position;
        }
    }
    std::vector<double> sums(group_key.size(), 0.0);
    for (std::size_t i = 0; i < kept.size(); ++i) {
        sums[static_cast<std::size_t>(group_id[i])] += kept[i];
    }

    if (!use_triu) {
        // result, without the symmetrisation.
        std::vector<std::int32_t> out_rows;
        std::vector<std::int32_t> out_cols;
        std::vector<double> out_values;
        out_rows.reserve(group_key.size());
        out_cols.reserve(group_key.size());
        out_values.reserve(group_key.size());
        for (std::size_t g = 0; g < group_key.size(); ++g) {
            out_rows.push_back(static_cast<std::int32_t>(group_key[g] >> 32));
            out_cols.push_back(static_cast<std::int32_t>(
                static_cast<std::uint32_t>(group_key[g] & 0xFFFFFFFF)));
            out_values.push_back(cast_to_dtype(sums[g], matrix.dtype()));
        }
        CsrMatrix result = CsrMatrix::from_coo(groups, groups, out_rows, out_cols,
                                               std::move(out_values), matrix.dtype());
        result.eliminate_zeros();
        return result;
    }

    // result + result.T - diagmatrix. With pDiagonal the diagonal keeps its
    // single contribution; without it the Python subtracts twice the diagonal,
    // so every diagonal entry becomes exactly zero and is pruned.
    std::vector<std::int32_t> out_rows;
    std::vector<std::int32_t> out_cols;
    std::vector<double> out_values;
    out_rows.reserve(group_key.size() * 2);
    out_cols.reserve(group_key.size() * 2);
    out_values.reserve(group_key.size() * 2);
    const auto emit = [&](std::int64_t row, std::int64_t col, double value) {
        out_rows.push_back(static_cast<std::int32_t>(row));
        out_cols.push_back(static_cast<std::int32_t>(col));
        out_values.push_back(value);
    };
    // The entries of `result` first and its transpose afterwards, so that a
    // cell fed from both sides accumulates in the order csr_plus_csr uses.
    for (std::size_t g = 0; g < group_key.size(); ++g) {
        const std::int64_t row = group_key[g] >> 32;
        const std::int64_t col =
            static_cast<std::int64_t>(static_cast<std::uint32_t>(group_key[g] & 0xFFFFFFFF));
        const double value = cast_to_dtype(sums[g], matrix.dtype());
        if (row == col && !diagonal) {
            continue;
        }
        emit(row, col, value);
    }
    for (std::size_t g = 0; g < group_key.size(); ++g) {
        const std::int64_t row = group_key[g] >> 32;
        const std::int64_t col =
            static_cast<std::int64_t>(static_cast<std::uint32_t>(group_key[g] & 0xFFFFFFFF));
        if (row == col) {
            continue;
        }
        emit(col, row, cast_to_dtype(sums[g], matrix.dtype()));
    }

    CsrMatrix result = CsrMatrix::from_coo(groups, groups, out_rows, out_cols,
                                           std::move(out_values), matrix.dtype());
    result.eliminate_zeros();
    return result;
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
    const std::int64_t rows = merged.matrix.rows();
    std::vector<double> column_sum(static_cast<std::size_t>(merged.matrix.cols()), 0.0);
    const std::vector<std::int64_t>& indptr = merged.matrix.indptr();
    const std::vector<std::int32_t>& indices = merged.matrix.indices();
    const std::vector<double>& values = merged.matrix.data();
    for (std::int64_t row = 0; row < rows; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            column_sum[static_cast<std::size_t>(indices[k])] += values[k];
        }
    }
    for (std::size_t column = 0; column < column_sum.size(); ++column) {
        if (column_sum[column] == 0.0) {
            merged.nan_bins.push_back(static_cast<std::int64_t>(column));
        }
    }
    return merged;
}

}  // namespace hicx
