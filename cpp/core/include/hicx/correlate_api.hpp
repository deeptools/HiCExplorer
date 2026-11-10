// The C++ behind hicexplorer.hicCorrelate's public functions (HiCExplorer
// 3.7.6, upstream c2ad8630): get_vectors and the numbers plot_correlation
// draws from. The hicCorrelate executable computes its own vectors over the
// union of many matrices (cpp/tools/hicCorrelate.cpp); get_vectors is the
// two matrix function of the Python module, reproducing scipy.sparse's
// arithmetic on the stored entries, including the order the entries come out
// in for matrices whose indices are not sorted.
//
// get_vectors(mat1, mat2) of 3.7.6:
//     _mat = mat1 + mat2;  _mat.data += 1
//     values1 = (_mat - mat1).data - 1
//     values2 = (_mat - mat2).data - 1
// scipy adds and subtracts two compressed matrices of one format with
// csr_binop_csr (the canonical merge when both have sorted, unique indices,
// else the general scan), drops the results equal to 0, and converts the
// second operand to the format of the first (coo is converted to csr, summing
// duplicates). The type of the arithmetic is the upcast of the operand types.

#ifndef HICX_CORRELATE_API_HPP
#define HICX_CORRELATE_API_HPP

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace hicx::correlate {

enum class SparseFormat { Csr, Csc, Coo };

// One scipy.sparse matrix, values already converted to the arithmetic type T.
// Csr/Csc: `indptr` (major + 1 entries) and `minor` (the indices); Coo: `major`
// holds the row and `minor` the column of every stored entry.
template <class T>
struct SparseInput {
    SparseFormat format = SparseFormat::Csr;
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> indptr;
    std::vector<std::int64_t> major;
    std::vector<std::int64_t> minor;
    std::vector<T> data;
};

// (values1, values2) of get_vectors for matrices of equal shape.
template <class T>
[[nodiscard]] std::pair<std::vector<T>, std::vector<T>> get_vectors(const SparseInput<T>& mat1,
                                                                    const SparseInput<T>& mat2);

extern template std::pair<std::vector<std::int32_t>, std::vector<std::int32_t>>
get_vectors(const SparseInput<std::int32_t>&, const SparseInput<std::int32_t>&);
extern template std::pair<std::vector<std::int64_t>, std::vector<std::int64_t>>
get_vectors(const SparseInput<std::int64_t>&, const SparseInput<std::int64_t>&);
extern template std::pair<std::vector<float>, std::vector<float>> get_vectors(
    const SparseInput<float>&, const SparseInput<float>&);
extern template std::pair<std::vector<double>, std::vector<double>> get_vectors(
    const SparseInput<double>&, const SparseInput<double>&);

// plot_correlation's numbers for a rows x cols float64 matrix (row major).
//
// The minimum of the values as ndarray.min() returns it (NaN if any value is
// NaN); throws the ValueError numpy raises for zero elements.
[[nodiscard]] double matrix_minimum(const std::vector<double>& values);

// The condition under which scipy.cluster.hierarchy.linkage warns that the
// observation matrix "looks suspiciously like an uncondensed distance matrix":
// square, a diagonal within 1e-8 of 0, no negative value, symmetric (allclose).
[[nodiscard]] bool looks_like_distance_matrix(const std::vector<double>& values, std::int64_t rows,
                                              std::int64_t cols);

// scipy.cluster.hierarchy.linkage(matrix, method='complete') of the rows: the
// (rows - 1) x 4 linkage matrix, row major. Throws the ValueErrors scipy
// raises (fewer than two observations, a non-finite value or distance).
[[nodiscard]] std::vector<double> complete_linkage_of_rows(const std::vector<double>& values,
                                                           std::int64_t rows, std::int64_t cols);

}  // namespace hicx::correlate

#endif  // HICX_CORRELATE_API_HPP
