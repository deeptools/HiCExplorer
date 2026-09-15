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
// **What is exact, and what this file does.**
//
//  1. Only the requested columns are back-transformed. dgeev spends 91 s of
//     its 405 s of CPU on mm9_reduced_chr1.cool (9,760 bins) in dtrevc3,
//     computing n eigenvectors of which hicPCA keeps two. dtrevc3's blocked
//     back-transform (NB > 1) processes the eigenvectors from the last to the
//     first in blocks and multiplies each block by the Schur vectors with one
//     dgemm. The block structure depends only on which diagonal blocks of the
//     Schur form are 2 by 2, so it is replayed first; then, for the blocks
//     that hold a wanted column only, the triangular solves are repeated with
//     the same LAPACK and BLAS calls and the block's dgemm is issued with
//     exactly dtrevc3's dimensions. dgebak and dgeev's normalisation act on
//     each column separately and are applied to the gathered columns.
// Tried and not kept: splitting the Hessenberg reduction's dominant product.
// dlahr2 multiplies the trailing block by a vector for every column. OpenBLAS's
// serial dgemv_n kernel (kernel/x86_64/dgemv_n_4.c) walks the rows in blocks of
// 2,048 and treats the last m % 4 rows with a scalar loop, so row ranges that
// start at multiples of 2,048 give a bit-identical product, and a transcription
// of dgehrd and dlahr2 with that one call split did reproduce dgeev_'s columns
// on mm9_reduced_chr1.cool. It did not pay: the product is limited by memory
// bandwidth, and the eigensolver took 254 s wall and 254 s CPU at one worker,
// 244 s and 296 s at two, 229 s and 312 s at three (load average 4 to 7). The
// dgemm calls cannot be split exactly at all: dgemm_kernel_16x2_skylakex.c
// accumulates a column differently depending on the size of the column group it
// falls in, and splitting them changed the result. So all of dgeev runs serial
// and the gain is the back-transform alone (cpp/OPTIMIZATION.md section 6:
// a change that does not measurably help is reverted).
//
// Verified bit for bit against dgeev_ itself: the unit test
// test_dgeev_selected.cpp compares the columns on a real chromosome block and,
// for the complex branch, on a non-symmetric matrix with complex eigenvalue
// pairs.

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
// Bit-identical to calling dgeev_ with the optimal workspace and reading the
// columns. `matrix` is consumed and released.
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
