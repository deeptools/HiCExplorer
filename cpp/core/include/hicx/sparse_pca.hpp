// hicPCA --eigenSolver lanczos: the leading eigenvectors of the covariance
// hicPCA takes, without the dense matrix (cpp/PLAN.md tier 11).
//
// hicPCA computes, per chromosome, C = np.cov(A) of the obs/exp block A and
// its eigenvectors. With the row means m and n bins,
//
//     X = A - m 1^T,        C = X X^T / (n - 1),
//
// so a product with C needs neither X nor C:
//
//     w = X^T v = A v - (m . v) 1        (A is symmetric)
//     C v = X w / (n - 1) = (A w - m (1 . w)) / (n - 1)
//
// that is two sparse products with A, one rank-one centering term per
// product and a diagonal scaling. The Krylov solver (Spectra, pinned in
// cpp/cmake/HicxSpectra.cmake) calls this product and keeps a Krylov basis of
// `ncv` vectors, so the memory is the sparse block plus O(ncv n) doubles.
//
// Rows the dense path cannot use. A stored NaN or infinity in row i of A
// makes the mean of that row, and with it row i and column i of np.cov,
// non-finite, and hicPCA.py:303-304 turns exactly those entries into zeros.
// The operator reproduces that by masking such rows and columns: the product
// is D C_good D with D the diagonal mask, which is the cleaned covariance.
// None of the test inputs contain one (measured on every chromosome of
// small_test_matrix.h5, small_test_matrix_50kb_res.h5 and
// mm9_reduced_chr1.cool under all three obs/exp methods).
//
// Order and sign. The eigenvectors come back ordered by descending
// eigenvalue, which is not always dgeev's column order (cpp/PLAN.md 5.4). The
// sign is made deterministic by making the largest magnitude component
// positive, the first one on a tie, which is what --compatMode v4 does;
// hicPCA then applies its gene or histone track rule on top, as it does for
// the dense solver.
//
// Determinism. The sparse products are split over rows and every row is
// accumulated sequentially, the two dot products of the centering term are
// taken sequentially, the start vector is a fixed function of n, and the
// solver's own dense arithmetic runs on one thread. The result is therefore
// byte-identical for any thread count.

#ifndef HICX_SPARSE_PCA_HPP
#define HICX_SPARSE_PCA_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "hicx/sparse_matrix.hpp"
#include "hicx/transform_ops.hpp"
#include "hicx/worker_pool.hpp"

namespace hicx {

struct LanczosDiagnostics {
    // Krylov subspace dimension used; 0 when the dense fallback ran.
    std::int64_t ncv = 0;
    std::int64_t iterations = 0;
    std::int64_t products = 0;
    // Chromosomes of at most kLanczosDenseLimit bins are solved densely.
    bool dense_fallback = false;
    // Rows masked because they carry a non-finite stored value.
    std::int64_t masked_rows = 0;
};

// At or below this many bins the covariance is formed and solved with dsyevr.
// The block is then at most 64 by 64 doubles, 32 KB, and a Krylov method has
// nothing to save: Spectra needs nev < ncv <= n.
inline constexpr std::int64_t kLanczosDenseLimit = 64;

// Relative residual tolerance handed to Spectra: a Ritz pair is accepted when
// its residual is below kLanczosTolerance * |eigenvalue|.
inline constexpr double kLanczosTolerance = 1e-13;

// The Krylov basis size for `nev` eigenvectors: max(2 nev + 1, 32), at most n.
[[nodiscard]] std::int64_t lanczos_subspace_size(std::int64_t nev, std::int64_t n);

// y = C x for the (masked) covariance of `obs_exp`, which must carry
// Symmetry::Full. Exposed for the unit test that checks it against the dense
// covariance.
class CovarianceOperator {
  public:
    CovarianceOperator(const CsrMatrix& obs_exp, int threads);
    [[nodiscard]] std::int64_t size() const noexcept { return n_; }
    [[nodiscard]] std::int64_t masked_rows() const noexcept { return masked_; }
    [[nodiscard]] const std::vector<unsigned char>& mask() const noexcept { return good_; }
    void apply(const double* x, double* y) const;

  private:
    void sparse_product(const double* x, double* y, bool masked_columns,
                        bool skip_masked_rows) const;

    const CsrMatrix* matrix_;
    std::int64_t n_;
    std::int64_t masked_ = 0;
    std::vector<unsigned char> good_;
    std::vector<double> means_;
    // Row ranges of roughly equal stored entries, one per piece of work. A
    // row's value does not depend on which range it is in.
    std::vector<std::int64_t> row_ranges_;
    std::unique_ptr<WorkerPool> pool_;
    mutable std::vector<double> u_;
    mutable std::vector<double> w_;
};

// The leading eigenvectors of the covariance of `obs_exp` for the 1-based
// `which`, ordered by descending eigenvalue. An index outside [1, n] gives an
// empty vector, as for the dense solvers. Throws when the solver does not
// converge.
[[nodiscard]] EigenResult covariance_eigenvectors_lanczos(const CsrMatrix& obs_exp,
                                                          const std::vector<int>& which,
                                                          int threads,
                                                          LanczosDiagnostics* diagnostics = nullptr);

}  // namespace hicx

#endif  // HICX_SPARSE_PCA_HPP
