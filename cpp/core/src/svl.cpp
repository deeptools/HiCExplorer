#include "hicx/svl.hpp"

#include <cmath>
#include <vector>

#include "hicx/numpy_compat.hpp"

namespace hicx::svl {

namespace {

// The entries of matrix[first:last, first:last] with the lower triangle
// filled in, in scipy's row major order, exact zeros skipped, passed to
// visit(row, column, value).
//
// For Symmetry::UpperTriangle a row's entries below the diagonal are the
// column's entries above it, so they are collected per column first. That
// transposed copy of the block is the one allocation beyond the matrix
// itself, and it is only needed when the order of summation matters, which is
// for float matrices.
template <class Visit>
void visit_block_row_major(const CsrMatrix& matrix, std::int64_t first, std::int64_t last,
                           Visit&& visit) {
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    const auto n = static_cast<std::size_t>(last - first);
    if (matrix.symmetry() != Symmetry::UpperTriangle) {
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                const std::int64_t column = indices[k];
                if (column >= first && column < last && data[k] != 0.0) {
                    visit(row, column, data[k]);
                }
            }
        }
        return;
    }
    std::vector<std::size_t> column_start(n + 1, 0);
    for (std::int64_t row = first; row < last; ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column > row && column < last && data[k] != 0.0) {
                ++column_start[static_cast<std::size_t>(column - first) + 1];
            }
        }
    }
    for (std::size_t c = 0; c < n; ++c) {
        column_start[c + 1] += column_start[c];
    }
    std::vector<std::int32_t> lower_column(column_start[n]);
    std::vector<double> lower_value(column_start[n]);
    std::vector<std::size_t> cursor(column_start.begin(), column_start.end() - 1);
    for (std::int64_t row = first; row < last; ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column > row && column < last && data[k] != 0.0) {
                const std::size_t slot = cursor[static_cast<std::size_t>(column - first)]++;
                lower_column[slot] = static_cast<std::int32_t>(row);
                lower_value[slot] = data[k];
            }
        }
    }
    for (std::int64_t row = first; row < last; ++row) {
        const auto local = static_cast<std::size_t>(row - first);
        for (std::size_t p = column_start[local]; p < column_start[local + 1]; ++p) {
            visit(row, static_cast<std::int64_t>(lower_column[p]), lower_value[p]);
        }
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column >= row && column < last && data[k] != 0.0) {
                visit(row, column, data[k]);
            }
        }
    }
}

}  // namespace

BlockSums block_sums(const CsrMatrix& matrix, std::int64_t first, std::int64_t last,
                     double max_distance) {
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    const bool upper = matrix.symmetry() == Symmetry::UpperTriangle;

    if (!upper) {
        std::size_t stored = 0;
        std::size_t nonzero = 0;
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                if (indices[k] >= first && indices[k] < last) {
                    ++stored;
                    if (data[k] != 0.0) {
                        ++nonzero;
                    }
                }
            }
        }
        if (stored != nonzero) {
            throw ExplicitZeros(stored, nonzero);
        }
    }
    const auto is_short = [max_distance](std::int64_t row, std::int64_t column) {
        const std::int64_t distance = row > column ? row - column : column - row;
        return static_cast<double>(distance) <= max_distance;
    };

    BlockSums result;
    result.kind = matrix.dtype_kind();
    if (result.kind == DType::Integer) {
        // Integer sums do not depend on the order of the additions, so the
        // stored triangle is enough. numpy's int64 accumulator wraps, hence
        // the unsigned arithmetic.
        std::uint64_t smaller = 0;
        std::uint64_t greater = 0;
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                const std::int64_t column = indices[k];
                if (column < first || column >= last || data[k] == 0.0) {
                    continue;
                }
                if (upper && column < row) {
                    continue;
                }
                const auto value = static_cast<std::uint64_t>(static_cast<std::int64_t>(data[k]));
                const std::uint64_t times = (upper && column != row) ? 2U : 1U;
                (is_short(row, column) ? smaller : greater) += times * value;
            }
        }
        result.smaller_int = static_cast<std::int64_t>(smaller);
        result.greater_int = static_cast<std::int64_t>(greater);
        result.ratio = static_cast<double>(result.smaller_int) / static_cast<double>(result.greater_int);
    } else if (result.kind == DType::Float32) {
        npy::PairwiseSumStream<float> smaller;
        npy::PairwiseSumStream<float> greater;
        visit_block_row_major(matrix, first, last, [&](std::int64_t row, std::int64_t column, double value) {
            (is_short(row, column) ? smaller : greater).add(static_cast<float>(value));
        });
        const float small_sum = smaller.result();
        const float large_sum = greater.result();
        result.smaller = static_cast<double>(small_sum);
        result.greater = static_cast<double>(large_sum);
        result.ratio = static_cast<double>(small_sum / large_sum);
    } else {
        npy::PairwiseSumStream<double> smaller;
        npy::PairwiseSumStream<double> greater;
        visit_block_row_major(matrix, first, last, [&](std::int64_t row, std::int64_t column, double value) {
            (is_short(row, column) ? smaller : greater).add(value);
        });
        result.smaller = smaller.result();
        result.greater = greater.result();
        result.ratio = result.smaller / result.greater;
    }
    result.kept = !(std::isinf(result.ratio) || std::isnan(result.ratio));
    return result;
}

}  // namespace hicx::svl
