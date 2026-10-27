// scipy.sparse operations whose rounding the public Python API reproduces,
// in the element type of the matrix (int32, int64, float32, float64):
//
//   * coo.tocsr() of a COO with duplicates: coo_tocsr scatters the entries
//     into their rows in storage order, csr_sort_indices sorts every row with
//     std::sort on (column, value) pairs compared by column only, and
//     csr_sum_duplicates adds equal columns left to right. std::sort is not
//     stable, so the order in which duplicates are added is the one
//     libstdc++'s introsort leaves; the same algorithm runs here.
//   * coo.diagonal() of a COO without canonical format: np.lexsort (stable)
//     by row and column, then np.add.reduceat, which adds the first value of
//     a run to numpy's pairwise sum of the rest.
//   * csr + csr and csr - csr of canonical matrices (csr_binop_csr): one
//     result per cell, op(a, 0) or op(0, b) where only one side stores the
//     cell, zero results dropped.
//   * matrix.sum(axis=0) of a CSR: a sequential accumulation per column in
//     storage order.
//
// Built from them, hicexplorer.hicMergeMatrixBins.running_window_merge.

#ifndef HICX_SCIPY_COO_HPP
#define HICX_SCIPY_COO_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include "hicx/numpy_compat.hpp"

namespace hicx::scipy_coo {

template <class T>
struct Csr {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> indptr;
    std::vector<std::int64_t> indices;
    std::vector<T> data;
};

template <class T>
struct Coo {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<T> data;
};

// Addition and subtraction in T as numpy does them (integers wrap).
template <class T>
T add(T a, T b) {
    if constexpr (std::is_integral_v<T>) {
        using U = std::make_unsigned_t<T>;
        return static_cast<T>(static_cast<U>(a) + static_cast<U>(b));
    } else {
        return a + b;
    }
}

template <class T>
T subtract(T a, T b) {
    if constexpr (std::is_integral_v<T>) {
        using U = std::make_unsigned_t<T>;
        return static_cast<T>(static_cast<U>(a) - static_cast<U>(b));
    } else {
        return a - b;
    }
}

template <class T>
bool nonzero(T value) {
    if constexpr (std::is_floating_point_v<T>) {
        return value != T(0) || std::isnan(value);
    } else {
        return value != T(0);
    }
}

// coo.tocsr() for a COO that is not in canonical format.
template <class T>
Csr<T> tocsr(const Coo<T>& coo) {
    Csr<T> out;
    out.rows = coo.rows;
    out.cols = coo.cols;
    out.indptr.assign(static_cast<std::size_t>(coo.rows) + 1, 0);
    if (coo.data.empty()) {
        return out;
    }
    for (const std::int64_t r : coo.row) {
        ++out.indptr[static_cast<std::size_t>(r) + 1];
    }
    for (std::size_t i = 1; i < out.indptr.size(); ++i) {
        out.indptr[i] += out.indptr[i - 1];
    }
    std::vector<std::int64_t> cursor(out.indptr.begin(), out.indptr.end() - 1);
    std::vector<std::int64_t> indices(coo.data.size());
    std::vector<T> data(coo.data.size());
    for (std::size_t k = 0; k < coo.data.size(); ++k) {
        const std::size_t dest = static_cast<std::size_t>(cursor[static_cast<std::size_t>(coo.row[k])]++);
        indices[dest] = coo.col[k];
        data[dest] = coo.data[k];
    }
    // csr_sort_indices: std::sort of (column, value) pairs, compared by column
    std::vector<std::pair<std::int64_t, T>> temp;
    std::vector<std::int64_t> new_indptr(out.indptr.size(), 0);
    for (std::int64_t i = 0; i < coo.rows; ++i) {
        const std::size_t begin = static_cast<std::size_t>(out.indptr[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(out.indptr[static_cast<std::size_t>(i) + 1]);
        temp.resize(end - begin);
        for (std::size_t k = begin; k < end; ++k) {
            temp[k - begin] = {indices[k], data[k]};
        }
        std::sort(temp.begin(), temp.end(),
                  [](const std::pair<std::int64_t, T>& a, const std::pair<std::int64_t, T>& b) {
                      return a.first < b.first;
                  });
        // csr_sum_duplicates
        std::size_t k = 0;
        while (k < temp.size()) {
            const std::int64_t column = temp[k].first;
            T x = temp[k].second;
            ++k;
            while (k < temp.size() && temp[k].first == column) {
                x = add(x, temp[k].second);
                ++k;
            }
            out.indices.push_back(column);
            out.data.push_back(x);
        }
        new_indptr[static_cast<std::size_t>(i) + 1] = static_cast<std::int64_t>(out.data.size());
    }
    out.indptr = std::move(new_indptr);
    return out;
}

// numpy's add.reduceat over one run: the first value plus the pairwise sum
// of the others (integers simply in order).
template <class T>
T reduceat_run(const std::vector<T>& values) {
    if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
        if (values.size() == 1) {
            return values[0];
        }
        return values[0] + hicx::npy::pairwise_sum_unbuffered(values.data() + 1, values.size() - 1);
    } else {
        T x = values[0];
        for (std::size_t k = 1; k < values.size(); ++k) {
            x = add(x, values[k]);
        }
        return x;
    }
}

// coo.diagonal() for a COO that is not in canonical format.
template <class T>
std::vector<T> diagonal(const Coo<T>& coo) {
    const std::int64_t n = std::min(coo.rows, coo.cols);
    std::vector<T> out(static_cast<std::size_t>(std::max<std::int64_t>(n, 0)), T(0));
    std::vector<std::size_t> positions;
    for (std::size_t k = 0; k < coo.data.size(); ++k) {
        if (coo.row[k] == coo.col[k]) {
            positions.push_back(k);
        }
    }
    std::stable_sort(positions.begin(), positions.end(),
                     [&](std::size_t a, std::size_t b) { return coo.row[a] < coo.row[b]; });
    std::size_t k = 0;
    std::vector<T> run;
    while (k < positions.size()) {
        const std::int64_t r = coo.row[positions[k]];
        run.clear();
        while (k < positions.size() && coo.row[positions[k]] == r) {
            run.push_back(coo.data[positions[k]]);
            ++k;
        }
        out[static_cast<std::size_t>(r)] = reduceat_run(run);
    }
    return out;
}

// csr_binop_csr for canonical inputs.
template <class T, class Op>
Csr<T> binop(const Csr<T>& a, const Csr<T>& b, Op op) {
    Csr<T> out;
    out.rows = a.rows;
    out.cols = a.cols;
    out.indptr.assign(static_cast<std::size_t>(a.rows) + 1, 0);
    for (std::int64_t i = 0; i < a.rows; ++i) {
        std::size_t p = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(i)]);
        const std::size_t p_end = static_cast<std::size_t>(a.indptr[static_cast<std::size_t>(i) + 1]);
        std::size_t q = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(i)]);
        const std::size_t q_end = static_cast<std::size_t>(b.indptr[static_cast<std::size_t>(i) + 1]);
        auto emit = [&](std::int64_t column, T value) {
            if (nonzero(value)) {
                out.indices.push_back(column);
                out.data.push_back(value);
            }
        };
        while (p < p_end && q < q_end) {
            if (a.indices[p] == b.indices[q]) {
                emit(a.indices[p], op(a.data[p], b.data[q]));
                ++p;
                ++q;
            } else if (a.indices[p] < b.indices[q]) {
                emit(a.indices[p], op(a.data[p], T(0)));
                ++p;
            } else {
                emit(b.indices[q], op(T(0), b.data[q]));
                ++q;
            }
        }
        for (; p < p_end; ++p) {
            emit(a.indices[p], op(a.data[p], T(0)));
        }
        for (; q < q_end; ++q) {
            emit(b.indices[q], op(T(0), b.data[q]));
        }
        out.indptr[static_cast<std::size_t>(i) + 1] = static_cast<std::int64_t>(out.data.size());
    }
    return out;
}

// The CSR of a diagonal (dia_matrix.tocsr(): zeros dropped).
template <class T>
Csr<T> diagonal_matrix(const std::vector<T>& diag, std::int64_t n) {
    Csr<T> out;
    out.rows = out.cols = n;
    out.indptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (std::int64_t i = 0; i < n; ++i) {
        if (static_cast<std::size_t>(i) < diag.size() && nonzero(diag[static_cast<std::size_t>(i)])) {
            out.indices.push_back(i);
            out.data.push_back(diag[static_cast<std::size_t>(i)]);
        }
        out.indptr[static_cast<std::size_t>(i) + 1] = static_cast<std::int64_t>(out.data.size());
    }
    return out;
}

// matrix.sum(0) == 0 for every column, the nan_bins hicMergeMatrixBins
// recomputes (np.flatnonzero). Integer sums are exact; float sums are
// accumulated in storage order in T.
template <class T>
std::vector<std::int64_t> zero_columns(const Csr<T>& matrix) {
    using Acc = std::conditional_t<std::is_integral_v<T>, std::int64_t, T>;
    std::vector<Acc> sums(static_cast<std::size_t>(matrix.cols), Acc(0));
    for (std::size_t k = 0; k < matrix.data.size(); ++k) {
        sums[static_cast<std::size_t>(matrix.indices[k])] += static_cast<Acc>(matrix.data[k]);
    }
    std::vector<std::int64_t> out;
    for (std::size_t c = 0; c < sums.size(); ++c) {
        if (sums[c] == Acc(0)) {
            out.push_back(static_cast<std::int64_t>(c));
        }
    }
    return out;
}

// hicMergeMatrixBins.running_window_merge on `matrix` (any CSR) for an odd
// num_bins > 1: every entry of triu(matrix) copied into the num_bins by
// num_bins window around it in the Python's order, positions outside the
// matrix dropped, the COO folded to its upper triangle and mirrored as
// new_ma + new_ma.T - diag.
template <class T>
Csr<T> running_window(const Csr<T>& matrix, std::int64_t num_bins) {
    const std::int64_t M = matrix.rows;
    const std::int64_t half = (num_bins - 1) / 2;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<T> data;
    for (std::int64_t i = 0; i < matrix.rows; ++i) {
        for (std::int64_t k = matrix.indptr[static_cast<std::size_t>(i)];
             k < matrix.indptr[static_cast<std::size_t>(i) + 1]; ++k) {
            if (matrix.indices[static_cast<std::size_t>(k)] >= i) {
                row.push_back(i);
                col.push_back(matrix.indices[static_cast<std::size_t>(k)]);
                data.push_back(matrix.data[static_cast<std::size_t>(k)]);
            }
        }
    }
    Coo<T> window;
    window.rows = window.cols = M;
    auto append = [&](std::int64_t dr, std::int64_t dc) {
        for (std::size_t k = 0; k < data.size(); ++k) {
            const std::int64_t r = row[k] + dr;
            const std::int64_t c = col[k] + dc;
            // the keep mask, then triu(new_ma, k=0)
            if (r > -1 && c > -1 && r < M && c < M && r <= c) {
                window.row.push_back(r);
                window.col.push_back(c);
                window.data.push_back(data[k]);
            }
        }
    };
    append(0, 0);
    for (std::int64_t i = 0; i < num_bins; ++i) {
        for (std::int64_t j = 0; j < num_bins; ++j) {
            const std::int64_t dr = j - half;
            const std::int64_t dc = i - half;
            if (dr == 0 && dc == 0) {
                continue;
            }
            append(dr, dc);
        }
    }
    const std::vector<T> diag = diagonal(window);
    Coo<T> transposed;
    transposed.rows = transposed.cols = M;
    transposed.row = window.col;
    transposed.col = window.row;
    transposed.data = window.data;
    const Csr<T> sum = binop(tocsr(window), tocsr(transposed), [](T a, T b) { return add(a, b); });
    return binop(sum, diagonal_matrix(diag, M), [](T a, T b) { return subtract(a, b); });
}

// The result of assign_zero_block: `inserted` false when every cell of the
// block was stored, so that scipy only overwrote values.
template <class T>
struct BlockAssignment {
    Csr<T> matrix;
    bool inserted = false;
};

template <class T>
BlockAssignment<T> assign_zero_block(const Csr<T>& matrix, bool sorted, std::int64_t row_first,
                                     std::int64_t row_last, std::int64_t col_first, std::int64_t col_last) {
    BlockAssignment<T> out;
    out.matrix.rows = matrix.rows;
    out.matrix.cols = matrix.cols;
    out.matrix.indptr.push_back(0);
    std::vector<char> present;
    std::vector<std::pair<std::int64_t, T>> row;
    for (std::int64_t r = 0; r < matrix.rows; ++r) {
        const std::size_t begin = static_cast<std::size_t>(matrix.indptr[static_cast<std::size_t>(r)]);
        const std::size_t end = static_cast<std::size_t>(matrix.indptr[static_cast<std::size_t>(r) + 1]);
        const bool in_rows = r >= row_first && r < row_last;
        if (in_rows) {
            present.assign(static_cast<std::size_t>(col_last - col_first), 0);
        }
        row.clear();
        for (std::size_t k = begin; k < end; ++k) {
            T value = matrix.data[k];
            const std::int64_t c = matrix.indices[k];
            if (in_rows && c >= col_first && c < col_last) {
                value = T(0);
                present[static_cast<std::size_t>(c - col_first)] = 1;
            }
            row.emplace_back(c, value);
        }
        if (in_rows) {
            std::vector<std::pair<std::int64_t, T>> inserted;
            for (std::int64_t c = col_first; c < col_last; ++c) {
                if (!present[static_cast<std::size_t>(c - col_first)]) {
                    inserted.emplace_back(c, T(0));
                }
            }
            out.inserted = out.inserted || !inserted.empty();
            row.insert(row.begin(), inserted.begin(), inserted.end());
        }
        for (const auto& [c, v] : row) {
            out.matrix.indices.push_back(c);
            out.matrix.data.push_back(v);
        }
        out.matrix.indptr.push_back(static_cast<std::int64_t>(out.matrix.data.size()));
    }
    if (out.inserted && sorted) {
        for (std::int64_t r = 0; r < matrix.rows; ++r) {
            const std::size_t begin = static_cast<std::size_t>(out.matrix.indptr[static_cast<std::size_t>(r)]);
            const std::size_t end = static_cast<std::size_t>(out.matrix.indptr[static_cast<std::size_t>(r) + 1]);
            row.clear();
            for (std::size_t k = begin; k < end; ++k) {
                row.emplace_back(out.matrix.indices[k], out.matrix.data[k]);
            }
            std::stable_sort(row.begin(), row.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
            for (std::size_t k = begin; k < end; ++k) {
                out.matrix.indices[k] = row[k - begin].first;
                out.matrix.data[k] = row[k - begin].second;
            }
        }
    }
    return out;
}

}  // namespace hicx::scipy_coo

#endif  // HICX_SCIPY_COO_HPP
