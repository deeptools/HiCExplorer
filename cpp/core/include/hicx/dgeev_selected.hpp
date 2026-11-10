// scipy.linalg.eig's right eigenvectors for a few columns, bit for bit.
//
// hicPCA.py:305 calls the general eigensolver dgeev on the covariance of every
// chromosome and then keeps one or two columns of the n it computed
// (hicPCA.py:314-318). Because it keeps them by position, and the position of
// an eigenvalue in dgeev's output depends on the last bits of every
// intermediate (cpp/PLAN.md 5.4), the port has to reproduce dgeev's
// arithmetic exactly, not only its mathematics.
//
// **Why not a multithreaded LAPACK.** Measured on chrX of
// small_test_matrix.h5 (4,485 bins) with OpenBLAS 0.3.28 from the dependency
// prefix, by running each stage of dgeev at one and at eight OpenBLAS threads
// on the same bit-identical covariance: the Hessenberg reduction dgehrd, the
// orthogonal factor dorghr, the Schur iteration dhseqr and the back-transform
// dtrevc3 each change bits with the thread count, and so does dsyrk, which
// np.cov calls. The threaded dgemm driver (driver/level3/level3_thread.c) even
// splits the inner dimension differently from the serial one, (min_l + 1) / 2
// against a rounding to GEMM_UNROLL_M, so no thread count reproduces the
// serial result. A tool whose output changes with the thread count violates
// cpp/OPTIMIZATION.md section 3, so every LAPACK call here runs on one
// OpenBLAS thread, and the parallelism is limited to what is exact by
// construction.
//
// **What this file does.** It calls LAPACK's dgeev_ from the dependency
// prefix with the optimal workspace and returns the requested columns, so the
// eigensolver is the library's, unmodified. An earlier version replayed
// dtrevc3's blocked back-transform for the requested columns only, which saved
// about a fifth of the CPU time on large chromosomes; it was a transcription
// of LAPACK source and was replaced by the library call.
//
// The unit test test_dgeev_selected.cpp compares the columns with dgeev_ on a
// real chromosome block and, for the complex branch, on a non-symmetric matrix
// with complex eigenvalue pairs.

#ifndef HICX_DGEEV_SELECTED_HPP
#define HICX_DGEEV_SELECTED_HPP

#include <vector>

#include "hicx/transform_ops.hpp"

namespace hicx {

// dgeev_('N', 'V') on the n-by-n column major matrix `matrix`, returning for
// every requested 1-based index what hicPCA writes: the real part of dgeev's
// right eigenvector in that column (the first column of a complex conjugate
// pair for both of its members) and the real part of the eigenvalue. An index
// outside [1, n] gives an empty vector, as numpy's out of range slice does.
//
// dgeev_ with the optimal workspace, read at those columns. `matrix` is
// consumed and released.
//
// OpenBLAS must be pinned to one thread (pin_blas_to_one_thread). Several
// calls may run at once on different threads, which OpenBLAS supports.
[[nodiscard]] EigenResult dgeev_selected_eigenvectors(DenseSymmetric& matrix,
                                                      const std::vector<int>& which);

// The same for a general (not necessarily symmetric) column major matrix held
// in a plain buffer, exposed so that the complex pair branch can be tested
// against dgeev_ on a matrix that has complex eigenvalues. `values_imaginary`
// receives wi for every requested index.
[[nodiscard]] EigenResult dgeev_selected_general(int n, double* matrix,
                                                 const std::vector<int>& which,
                                                 std::vector<double>* values_imaginary);

}  // namespace hicx

#endif  // HICX_DGEEV_SELECTED_HPP
