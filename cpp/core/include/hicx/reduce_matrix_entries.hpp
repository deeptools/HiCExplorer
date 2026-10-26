// hicexplorer.reduceMatrix.reduce_matrix over the stored entries of any
// scipy matrix, for the public Python API.
//
// hicx::reduce_matrix (hicx/reduce_matrix.hpp) is the tools' version on a
// CSR. The Python sums with np.bincount in the order the entries of
// triu(matrix, k=0, format='coo') (or matrix.tocoo()) appear, which for a CSR
// is row major but for a CSC or a COO is their own storage order, and that
// order decides the rounding of every float sum. This version takes the
// entries in that order, so it reproduces the Python for every input format;
// on a CSR it computes what hicx::reduce_matrix computes
// (cpp/tests/test_api_loops_correct.cpp).
//
// The arithmetic of the dtype is the caller's: values are summed in float64
// (np.bincount), cast back with `cast` (coo_matrix(dtype=ma.dtype)), and the
// symmetrisation result + result.T - diag adds and subtracts with `add` and
// `subtract` in the dtype.

#ifndef HICX_REDUCE_MATRIX_ENTRIES_HPP
#define HICX_REDUCE_MATRIX_ENTRIES_HPP

#include <cstdint>
#include <functional>
#include <vector>

namespace hicx::reduce_entries {

struct Entries {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<double> data;  // the values as float64 (bincount's weights)
};

struct Arithmetic {
    std::function<double(double)> cast;              // float64 -> dtype -> float64
    std::function<double(double, double)> add;       // a + b in the dtype
    std::function<double(double, double)> subtract;  // a - b in the dtype
    std::function<double(double)> twice;             // a * 2 in the dtype
};

// The entries of triu(entries, k=0): those with row <= col, order kept.
[[nodiscard]] Entries upper_triangle(const Entries& entries);

// The number of NaN values among the entries.
[[nodiscard]] std::size_t count_nan(const Entries& entries);

// The reduced matrix as sorted (row, col) cells, before or after the
// symmetrisation. `group_of[i]` is the group of input index i, -1 for none.
// With symmetrise false: the coo_matrix of the sums, with its zeros
// eliminated (the non-triu result). With symmetrise true: the CSR of
// result + result.T - diagmatrix (diagmatrix doubled unless `diagonal`),
// without zeros.
struct Cells {
    std::int64_t size = 0;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<double> data;
};

[[nodiscard]] Cells reduce(const Entries& entries, const std::vector<std::int64_t>& group_of,
                           std::int64_t groups, bool symmetrise, bool diagonal,
                           const Arithmetic& arithmetic);

}  // namespace hicx::reduce_entries

#endif  // HICX_REDUCE_MATRIX_ENTRIES_HPP
