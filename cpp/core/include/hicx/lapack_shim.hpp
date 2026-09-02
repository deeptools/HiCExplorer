// The handful of BLAS and LAPACK routines the transform kernels call.
//
// There is no cblas or lapacke header in the dependency prefix, only the
// library, so the Fortran symbols are declared here. That is deliberate as
// well as necessary: cpp/PLAN.md 3.4 requires the port to call the *same*
// OpenBLAS build that numpy and scipy are linked against
// (libopenblasp-r0.3.28.so in $HICX_DEPS/lib), because the column order and
// the eigenvector signs scipy.linalg.eig produces are properties of that
// build, and hicPCA compatibility mode reproduces them by calling into it
// rather than by imitating it.
//
// Fortran conventions that matter at every call site:
//   * arrays are column major, so an m-by-n matrix has leading dimension m and
//     element (i, j) lives at a[i + j * lda];
//   * every argument is passed by pointer;
//   * character arguments are passed as a pointer to the character, with a
//     hidden trailing length argument that every relevant compiler lets the
//     caller omit for length one strings.
//
// A symmetric matrix is its own transpose, so a symmetric matrix filled row by
// row in row major order is already the column major matrix LAPACK wants, and
// none of the callers transposes anything.

#ifndef HICX_LAPACK_SHIM_HPP
#define HICX_LAPACK_SHIM_HPP

extern "C" {

// dsyevr: eigenvalues and optionally eigenvectors of a real symmetric matrix,
// with range = 'I' selecting eigenvalues il..iu by index in ascending order.
// Computing only the requested pairs is what takes hicPCA's corrected mode
// from n eigenvectors to two.
void dsyevr_(const char* jobz, const char* range, const char* uplo, const int* n,
             double* a, const int* lda, const double* vl, const double* vu,
             const int* il, const int* iu, const double* abstol, int* m, double* w,
             double* z, const int* ldz, int* isuppz, double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);

// dgeev: eigenvalues and right eigenvectors of a general real matrix. This is
// what scipy.linalg.eig calls, and hicPCA calls scipy.linalg.eig on a
// symmetric covariance matrix (hicPCA.py:305). The eigenpairs come back in
// LAPACK's own order, which is not sorted; see cpp/PLAN.md 5.4.
//
// A complex conjugate pair (lambda, conj(lambda)) occupies two consecutive
// columns of vr: the real part in column j and the imaginary part in column
// j+1. scipy assembles them into two complex columns and hicPCA writes out
// only the real part, so a caller that wants what hicPCA writes takes column
// j for both members of the pair.
void dgeev_(const char* jobvl, const char* jobvr, const int* n, double* a,
            const int* lda, double* wr, double* wi, double* vl, const int* ldvl,
            double* vr, const int* ldvr, double* work, const int* lwork, int* info);

// dlamch: machine parameters. dsyevr's documented most-accurate abstol is
// 2 * dlamch('S'), the safe minimum.
double dlamch_(const char* cmach);

// OpenBLAS's own thread control, declared weak so that the library still links
// against a reference BLAS that does not provide it.
//
// This is not a performance knob, it is a correctness one. Measured on
// small_test_matrix_50kb_res.h5: the eigenvectors dgeev returns depend on the
// number of BLAS threads, and two chromosomes come back sign-inverted between
// OPENBLAS_NUM_THREADS=1 and 16, in the C++ port exactly as in the Python
// reference. A tool whose output changes with an environment variable is not
// reproducible, so the eigensolver is pinned to one BLAS thread and the
// parallelism is taken where it can be made deterministic, over the covariance
// rows. It is also cheaper: the same run costs 0.31 s of CPU at one BLAS
// thread and 75.7 s at sixteen, because OpenBLAS spins after every parallel
// region and the per chromosome blocks are small.
__attribute__((weak)) void openblas_set_num_threads(int num_threads);

}  // extern "C"

#endif  // HICX_LAPACK_SHIM_HPP
