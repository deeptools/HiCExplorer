// The three sparse kernels that ICE and Knight-Ruiz balancing spend their time
// in, threaded and vectorised under the determinism rules of
// cpp/OPTIMIZATION.md section 3.
//
// All three work on a *diagonal block* of an upper triangle CSR matrix, seen as
// the symmetric matrix that triangle represents. The whole genome case is the
// block [0, nbins); hicCorrectMatrix --perchr uses one block per chromosome.
// A block never copies anything: for a diagonal block every stored column of a
// row inside the block is already inside the block, because the storage is
// upper triangular and the column of a stored entry is at least its row, so the
// block is described by two per-row offsets into the shared arrays.
//
// Determinism, which is the constraint the whole file is shaped by:
//
//   * The number of row partitions is a function of the matrix alone
//     (kBlockSplitThreshold below), never of the thread count. A partition is
//     reduced sequentially and the partitions are combined in index order, so
//     the result at --threads 1 and at --threads 16 is the same bit pattern.
//   * No accumulator is shared between threads. Each partition owns a private
//     scatter buffer of one double per row of the block, which is 200 kB for
//     the largest matrix in the corpus and therefore stays in L2.
//   * The vectorised paths are used only where vectorising cannot change a
//     result: the two scaling kernels are elementwise, so the AVX2 path
//     computes exactly the same products as the scalar path in the same
//     pairing, and the unit tests assert bit equality rather than a tolerance.
//     The reductions are left scalar, both because their accumulation order is
//     observable against the Python reference and because they are bound by
//     memory bandwidth, not by arithmetic.
//
// Marginal order, and why it is reproduced rather than approximated:
// scipy computes W.sum(axis=1) with coo_matvec, which walks the COO in storage
// order and adds into the row accumulator. The COO comes from .tocoo() on the
// symmetric CSR, so for row i the terms arrive in ascending column: first the
// entries mirrored from the strict upper triangle of earlier rows, then row i's
// own stored entries from the diagonal on. The two phase scatter and gather
// below reproduce exactly that sequence, so with one partition the marginals
// are bit identical to the reference.

#ifndef HICX_MATH_SPARSE_KERNELS_HPP
#define HICX_MATH_SPARSE_KERNELS_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "hicx/sparse_matrix.hpp"

namespace hicx::kernels {

// Which vectorised path the process selected at start up. Fixed for the whole
// run, so two runs on the same machine cannot differ.
enum class SimdPath { Scalar, Avx2 };

[[nodiscard]] SimdPath active_simd_path();
[[nodiscard]] const char* simd_path_name(SimdPath path);
// Test hook. Forcing the scalar path is what lets the unit tests compare the
// two implementations against each other on the same input.
void force_simd_path(SimdPath path);

// A square diagonal block of an upper triangle CSR, addressed in local
// coordinates 0 .. size()-1.
class DiagonalBlock {
  public:
    // The whole matrix.
    explicit DiagonalBlock(CsrMatrix& matrix);
    // Rows and columns [first, last) of `matrix`.
    DiagonalBlock(CsrMatrix& matrix, std::int64_t first, std::int64_t last);

    [[nodiscard]] std::int64_t size() const noexcept { return last_ - first_; }
    [[nodiscard]] std::int64_t first() const noexcept { return first_; }
    [[nodiscard]] std::int64_t last() const noexcept { return last_; }
    [[nodiscard]] std::size_t stored() const noexcept { return stored_; }

    [[nodiscard]] std::int64_t row_begin(std::int64_t local_row) const {
        return row_begin_[static_cast<std::size_t>(local_row)];
    }
    [[nodiscard]] std::int64_t row_end(std::int64_t local_row) const {
        return row_end_[static_cast<std::size_t>(local_row)];
    }
    [[nodiscard]] const std::int32_t* indices() const noexcept { return indices_; }
    [[nodiscard]] double* data() const noexcept { return data_; }

    // Row partition boundaries, size partitions()+1, balanced by stored entry
    // count. Derived from the matrix alone.
    [[nodiscard]] const std::vector<std::int64_t>& partition() const noexcept {
        return partition_;
    }
    [[nodiscard]] std::size_t partitions() const noexcept {
        return partition_.size() - 1;
    }

  private:
    void build();

    std::int64_t first_ = 0;
    std::int64_t last_ = 0;
    const std::int32_t* indices_ = nullptr;
    double* data_ = nullptr;
    const std::int64_t* indptr_ = nullptr;
    std::size_t stored_ = 0;
    std::vector<std::int64_t> row_begin_;
    std::vector<std::int64_t> row_end_;
    std::vector<std::int64_t> partition_;
};

// A partition is only split off when the block holds at least this many stored
// entries. Below it a single partition is used, which keeps every small matrix
// bit identical to the Python reference and costs nothing, because a matrix
// this size is not where the time goes.
inline constexpr std::size_t kBlockSplitThreshold = 2'000'000;
inline constexpr std::size_t kMaxPartitions = 64;

// Test hook. The unit tests need a partitioned block without a two million
// entry matrix; nothing outside the tests calls this.
void set_block_split_threshold(std::size_t entries);
[[nodiscard]] std::size_t block_split_threshold();

// out[i] = sum over row i of the symmetric matrix the block represents, in
// scipy's coo_matvec order. `out` must hold block.size() elements.
void symmetric_marginals(const DiagonalBlock& block, double* out, int threads);

// One ICE pass, fused: every stored value is multiplied by scale[row] and then
// by scale[col], and the marginals of the *scaled* matrix are accumulated in
// the same traversal. Fusing the two halves the memory traffic of the
// iteration, which is what the loop is bound by.
//
// Returns the largest absolute value the scaling produced, so that the caller
// can reproduce the np.any(W.data > 1e100) guard without a third pass.
double scale_and_marginals(const DiagonalBlock& block, const double* scale, double* out,
                           int threads);

// The unfused reference for the kernel above, kept because the unit tests
// compare the two and because scale_and_marginals is easy to get subtly wrong.
double scale_rows_and_cols(const DiagonalBlock& block, const double* scale, int threads);

// y = (U + strict_upper(U)^T + diagonal_addend * I) * x, the matrix Knight-Ruiz
// iterates against. `y` must hold block.size() elements.
void symmetric_spmv(const DiagonalBlock& block, const double* x, double* y,
                    double diagonal_addend, int threads);

// value <- value * x[row] * x[col] for every stored entry of the block, which
// is what krbalancing's compute_normalised_matrix does, in the same pairing.
void scale_symmetric(const DiagonalBlock& block, const double* x, int threads);

}  // namespace hicx::kernels

#endif  // HICX_MATH_SPARSE_KERNELS_HPP
