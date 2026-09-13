// The obs/exp, Pearson, covariance and eigenvector machinery shared by
// hicPCA and hicTransform.
//
// hicPCA.py:284-305 and hicTransform.py:93-138 run the same pipeline on the
// same per chromosome submatrix, so it lives here once. hicPCA adds the
// eigendecomposition on top; hicTransform stops at the transformed matrix.
//
// Three properties of the Python that the port has to carry, because they are
// arithmetic and not style:
//
//  1. **The obs/exp kernels round through float32.** Every one of them does
//     `pSubmatrix.data = pSubmatrix.data.astype(np.float32)` before dividing
//     (utilities.py:503, :533, :583). obs_exp_matrix_lieberman and
//     obs_exp_matrix then divide the whole float32 array by a float64 array,
//     which numpy promotes back to float64, and finally cast to the dtype the
//     matrix had *before* the float32 step. obs_exp_matrix_non_zero divides
//     element by element and assigns each result back into the float32 array,
//     so its output stays float32 rounded. The two are not the same
//     computation and they are not interchangeable.
//  2. **The distance index is not the distance.** The expected value arrays
//     are accumulated at index abs(i - j) but read at index
//     ceil(abs(i - j) / 2) in obs_exp_matrix_lieberman and obs_exp_matrix
//     (utilities.py:499, :583). obs_exp_matrix_non_zero reads at abs(i - j).
//     Reproduced, not fixed.
//  3. **The covariance is computed from the sparse matrix, not from a dense
//     copy of it.** numpy centres a dense n-by-n block and calls dgemm on it,
//     which needs two dense blocks live and, on the designated hicPCA input,
//     costs 471 s and 4,070 MB. The port uses the algebraic identity
//
//         cov(A) = (A A^T - n m m^T) / (n - 1),   m = row means of A
//
//     which for a symmetric A is (A A - n m m^T) / (n - 1) and can be
//     evaluated straight out of the CSR into the single dense result block.
//     Measured against np.cov on the validation inputs: max relative
//     deviation 2.1e-12 with the covariance entries of
//     small_test_matrix.h5 chrX and 4.3e-14 relative to a floor of 1 on
//     mm9_reduced_chr1.cool, and the leading eigenvectors agree to 2.2e-16.
//     Critically, scipy.linalg.eig returns the same columns in the same order
//     with the same signs from either matrix, which is checked in
//     hicexplorer/test/general/test_hicPCA.py.
//
// Threading and SIMD follow cpp/OPTIMIZATION.md: the covariance is partitioned
// by output row, a fixed index range, each row is reduced sequentially, and
// rows are written into disjoint memory, so the result does not depend on the
// thread count. cpp/tests/test_transform_ops.cpp asserts that byte for byte at
// 1, 4 and 16 workers on a real chromosome block.
//
// The dense per-row finalisation has three implementations: an always compiled
// scalar reference, an AVX2+FMA variant and an AVX-512 variant, chosen once per
// process by __builtin_cpu_supports so that the choice cannot vary within a
// run. Measured in isolation on a Zen 4 core, over row lengths that occur in
// the corpus, with the arrays in cache:
//
//     row length   scalar    AVX2            AVX-512
//     512          0.155 s   0.041 s (3.8x)  0.028 s (5.5x, 1.45x over AVX2)
//     2,048        0.124 s   0.043 s (2.9x)  0.037 s (3.3x, 1.15x over AVX2)
//     9,760        0.119 s   0.038 s (3.1x)  0.036 s (3.3x, 1.05x over AVX2)
//     24,926       0.139 s   0.040 s (3.5x)  0.037 s (3.8x, 1.09x over AVX2)
//
// AVX-512 is faster than AVX2 at every size measured, which is the condition
// cpp/OPTIMIZATION.md section 2 sets for shipping a third path, but only
// barely at the large sizes: Zen 4 runs AVX-512 on a 256-bit datapath, so the
// gain comes from the wider registers rather than from doubled throughput, and
// it disappears into the memory traffic once the row no longer fits in L1.
// cpp/tests/test_transform_ops.cpp asserts the dispatched path agrees with the
// scalar reference; not bit for bit, because the fused multiply-add rounds
// once where the scalar expression rounds twice, but far inside the ED gate.

#ifndef HICX_TRANSFORM_OPS_HPP
#define HICX_TRANSFORM_OPS_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// --------------------------------------------------------------------------
// the two hiCMatrix mutators hicPCA and hicTransform share
//
// These belong next to reorder_bins and delete_bins in adjust_ops.hpp and
// should move there. They sit here so that the hicPCA port does not edit a
// header another workstream is changing at the same time.

// hiCMatrix.keepOnlyTheseChr (HiCMatrix.py:603). Keeps the bins of the named
// chromosomes in ascending bin order, carries the bin table, the correction
// factors and the NaN bins along, and clears distance_counts, which the Python
// does unconditionally at :677. Throws when a name is not in the matrix, which
// is the ValueError of :614.
void keep_only_chromosomes(MatrixData& data, const std::vector<std::string>& chromosomes);

// utilities.enlarge_bins lives in obsexp_ops.hpp, added by the hicFindTADs
// workstream; hicPCA calls it from there rather than shipping a second copy.

// --------------------------------------------------------------------------
// expected interactions per genomic distance

// utilities.expected_interactions_in_distance (utilities.py:293). Accumulates
// the stored values at index abs(row - col) and divides by
// length_chromosome - i * chromosome_count. No NaN or infinity cleanup, unlike
// its two siblings, so a distance whose divisor is zero yields an infinity
// that the caller has to deal with.
[[nodiscard]] std::vector<double> expected_interactions_in_distance(
    const CsrMatrix& matrix, std::int64_t length_chromosome,
    std::int64_t chromosome_count);

// utilities.expected_interactions_non_zero (utilities.py:317): the mean over
// the stored entries at each distance, NaN and infinity mapped to zero.
[[nodiscard]] std::vector<double> expected_interactions_non_zero(const CsrMatrix& matrix);

// utilities.expected_interactions (utilities.py:356) with pThreads=None: the
// sum at each distance divided by np.arange(n + 1, 1, -1), that is by
// n + 1 - i, which is one more than the number of cells on diagonal i.
// Reproduced as it stands.
[[nodiscard]] std::vector<double> expected_interactions(const CsrMatrix& matrix);

// --------------------------------------------------------------------------
// obs/exp

// utilities.obs_exp_matrix_lieberman (utilities.py:488).
void obs_exp_lieberman_in_place(CsrMatrix& matrix, std::int64_t length_chromosome,
                                std::int64_t chromosome_count);

// utilities.obs_exp_matrix_non_zero (utilities.py:510), including the Homer
// ligation factor and the trailing eliminate_zeros.
void obs_exp_non_zero_in_place(CsrMatrix& matrix, bool ligation_factor);

// utilities.obs_exp_matrix (utilities.py:554) with pDistance=None.
void obs_exp_in_place(CsrMatrix& matrix);

// utilities.convertNansToZeros followed by convertInfsToZeros over the stored
// values, which is what every caller applies to an obs/exp result.
void convert_nans_and_infs_to_zeros(CsrMatrix& matrix);

// --------------------------------------------------------------------------
// dense symmetric block

// A dense n-by-n symmetric matrix. Row major and column major agree for a
// symmetric matrix, so the same buffer is handed to LAPACK unchanged.
class DenseSymmetric {
  public:
    DenseSymmetric() = default;
    // Zero filled, which is the safe default for any caller that fills the
    // block partially.
    explicit DenseSymmetric(std::int64_t n) : DenseSymmetric(uninitialized(n)) {
        std::fill(values_.get(), values_.get() + element_count(), 0.0);
    }

    // Allocated and not touched. Measured on mm9_reduced_chr1.cool, a 9,760
    // bin block: zero filling the 762 MB the block occupies costs 0.44 s of
    // CPU and 0.56 s of wall time, which is the whole cost of the covariance
    // there, and every byte of it is overwritten immediately afterwards. Only
    // a caller that writes every element before anything reads it may use
    // this, which is covariance_of_symmetric and nothing else.
    [[nodiscard]] static DenseSymmetric uninitialized(std::int64_t n) {
        DenseSymmetric block;
        block.n_ = n;
        if (n > 0) {
            // new[] without an initialiser leaves the doubles indeterminate,
            // so the pages are first touched by whoever writes them, in
            // parallel rather than serially.
            block.values_.reset(new double[block.element_count()]);
        }
        return block;
    }

    [[nodiscard]] std::int64_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t element_count() const noexcept {
        return static_cast<std::size_t>(n_) * static_cast<std::size_t>(n_);
    }
    [[nodiscard]] double* data() noexcept { return values_.get(); }
    [[nodiscard]] const double* data() const noexcept { return values_.get(); }
    [[nodiscard]] double* row(std::int64_t i) noexcept {
        return values_.get() + static_cast<std::size_t>(i) * static_cast<std::size_t>(n_);
    }
    [[nodiscard]] const double* row(std::int64_t i) const noexcept {
        return values_.get() + static_cast<std::size_t>(i) * static_cast<std::size_t>(n_);
    }
    [[nodiscard]] double at(std::int64_t i, std::int64_t j) const noexcept {
        return values_[static_cast<std::size_t>(i) * static_cast<std::size_t>(n_) +
                       static_cast<std::size_t>(j)];
    }
    [[nodiscard]] double& at(std::int64_t i, std::int64_t j) noexcept {
        return values_[static_cast<std::size_t>(i) * static_cast<std::size_t>(n_) +
                       static_cast<std::size_t>(j)];
    }
    void release() {
        values_.reset();
        n_ = 0;
    }

  private:
    std::int64_t n_ = 0;
    std::unique_ptr<double[]> values_;
};

// np.cov(dense(matrix)) with the identity above. `matrix` must carry
// Symmetry::Full, which is what materialize_full() produces and what
// hiCMatrix's fillLowerTriangle leaves the Python with. `threads` partitions
// the output rows; the result is independent of it.
//
// The result is raw: no NaN or infinity cleanup, because hicPCA derives the
// Pearson matrix from the *uncleaned* covariance (hicPCA.py:296 computes
// np.corrcoef on the obs/exp block, not on the cleaned covariance) and only
// then cleans each of the two separately.
[[nodiscard]] DenseSymmetric covariance_of_symmetric(const CsrMatrix& matrix,
                                                     int threads);

// utilities.convertNansToZeros followed by convertInfsToZeros over a dense
// block, which is what hicPCA.py:303-304 applies to the covariance.
void zero_non_finite_in_place(DenseSymmetric& block, int threads);

// sqrt(diag(covariance)), the denominators numpy.corrcoef divides by.
[[nodiscard]] std::vector<double> pearson_scaling(const DenseSymmetric& covariance);

// One row of np.corrcoef: out[j] = clip(cov(i, j) / d[i] / d[j], -1, 1) with
// NaN mapped to zero. Reading a row without touching the block is what lets
// hicPCA write the Pearson matrix and then still run the eigensolver on the
// covariance, without a second dense block.
void pearson_row(const DenseSymmetric& covariance, const std::vector<double>& scaling,
                 std::int64_t i, double* out);

// np.corrcoef(dense(matrix)) over the whole block, in place, for a caller that
// does not need the covariance afterwards.
void covariance_to_pearson_in_place(DenseSymmetric& covariance, int threads);

// --------------------------------------------------------------------------
// dense correlation rows without the dense block

// The covariance or the Pearson correlation of a symmetric CSR matrix,
// produced one row at a time instead of as an n-by-n block.
//
// Why this exists beside covariance_of_symmetric. hicPCA works chromosome by
// chromosome and needs the whole covariance block resident anyway, because the
// eigensolver takes it. hicTransform does not: it writes the transform out and
// never looks at it again, and its --method pearson branch is the second worst
// memory ratio in the corpus, 5,150 MB against a 19.9 MB working set on
// Li_et_al_2015.h5 (cpp/PLAN.md 4.3). Removing the Python's five live copies
// still leaves two blocks that do not fit the budget of cpp/PLAN.md 4.5:
// the dense input block, 987 MB on that matrix, and the dense *output*, which
// as a CsrMatrix would be 123 M stored entries, another 1.5 GB. A row source
// removes both, so the resident set is the sparse matrix plus O(n) vectors
// whatever the density of the result.
//
// The arithmetic is covariance_of_symmetric's, element for element:
//     cov = (A A^T - n m m^T) / (n - 1),  m = row means of A,
// accumulated over the same index ranges in the same order and finalised by
// the same covariance_row_finalise, so the two agree bit for bit.
// cpp/tests/test_transform_ops.cpp asserts that on a real chromosome block.
//
// The Pearson scaling needs sqrt(diag(cov)) before any row can be scaled, so
// the covariance rows are computed twice for Kind::Pearson: once to collect
// the diagonal and once to emit. An O(nnz) shortcut for the diagonal was
// written and rejected. For a matrix in its explicit symmetric form the raw
// accumulator satisfies
//     acc(i, i) = sum_k A(i, k) * A(k, i) = sum_k A(i, k)^2
// summed in increasing k, which is the order the row accumulation visits it
// in, so the *accumulator* is reproduced exactly; but the finalisation is not.
// covariance_row_finalise dispatches to an AVX2 or AVX-512 body that uses one
// fused multiply-add where the scalar expression rounds twice, and whether
// element i of row i falls in the vector body or in the scalar tail depends on
// i modulo the lane count. A scalar shortcut therefore disagrees with the
// dispatched kernel in the last bit for most of the diagonal, which was
// measured on chrX of small_test_matrix_50kb_res.h5, and the second division
// by sqrt of it turns that into a last-bit difference across the whole Pearson
// matrix. Matching it would mean hard-coding the lane count of whichever path
// was dispatched, which is exactly the coupling cpp/OPTIMIZATION.md section 3
// exists to prevent. The extra pass was expected to be expensive and measured
// not to be: hicTransform --method pearson on Li_et_al_2015.h5 goes from 3.97
// to 4.06 s of CPU, 2.3 %, because the run is dominated by the blosc
// compression of a 373 MB result rather than by the accumulation, against a
// Python reference that spends 66.3 s. It buys an invariant a unit test can
// state: DenseCorrelationRows and covariance_of_symmetric are the same
// arithmetic, bit for bit.
//
// `blocks` are the bin ranges the correlation is taken within: one range over
// every bin is hicTransform's whole-matrix branch, one range per chromosome is
// --perChromosome. A row outside every range is all zeros.
class DenseCorrelationRows {
  public:
    enum class Kind {
        // np.cov(dense(block)), raw. No NaN or infinity cleanup, because
        // hicTransform.py:249 applies none to the covariance branch, unlike
        // the pearson one.
        Covariance,
        // np.corrcoef(dense(block)): cov divided by sqrt(diag) twice, clipped
        // into [-1, 1], with NaN mapped to zero as convertNansToZeros does.
        Pearson,
    };

    // `threads` splits the diagonal pass Kind::Pearson needs. It changes no
    // value: every row is produced by one worker over a fixed contiguous
    // range.
    DenseCorrelationRows(const CsrMatrix& matrix, std::vector<BinRange> blocks, Kind kind,
                         int threads = 1);

    [[nodiscard]] std::int64_t size() const noexcept { return n_; }

    // Fills out[0 .. size()-1] with row `i` of the result. Const and free of
    // shared mutable state, so several threads may fill different rows at
    // once; the values do not depend on how the rows were distributed.
    void fill_row(std::int64_t i, double* out) const;

    // The covariance row alone, before the Pearson rescaling. Exposed so that
    // the diagonal pass and fill_row cannot drift apart.
    void fill_covariance_row(std::int64_t i, double* out) const;

  private:
    const CsrMatrix* matrix_ = nullptr;
    std::int64_t n_ = 0;
    Kind kind_ = Kind::Covariance;
    std::vector<BinRange> blocks_;
    // Per bin: the block it belongs to, or -1.
    std::vector<std::int32_t> block_of_bin_;
    // Per bin: the mean of its row over its own block, the block's 1/(nb - 1)
    // and, for Pearson, sqrt of the covariance diagonal.
    std::vector<double> means_;
    std::vector<double> scaling_;
};

// --------------------------------------------------------------------------
// eigenvectors

// Pins OpenBLAS to one thread for the rest of the process, which is what makes
// the eigensolver reproducible: measured, dgeev returns sign-inverted
// eigenvectors for two chromosomes of small_test_matrix_50kb_res.h5 between
// OPENBLAS_NUM_THREADS=1 and 16. A tool that calls the eigensolver must call
// this first. It is a no operation when the linked BLAS is not OpenBLAS.
void pin_blas_to_one_thread();

enum class EigenSolver {
    // scipy.linalg.eig, the general non-symmetric solver, on the symmetric
    // covariance matrix, results taken in LAPACK's own column order. This is
    // what hicPCA.py:305 does and it is the only way to reproduce which
    // eigenvector the Python emits: dgeev's column order coincides with
    // descending eigenvalue order only when the spectrum is well separated,
    // and on small_test_matrix.h5 chrX the two leading eigenvalues differ by
    // 8.6e-04 relative and come back swapped. See cpp/PLAN.md 5.4.
    Dgeev,
    // dsyevr with range='I', only the requested eigenvectors, sorted by
    // descending eigenvalue, sign fixed by making the largest magnitude
    // component positive. Faster and smaller, and a different answer wherever
    // dgeev's order is not the sorted one.
    Dsyevr,
};

struct EigenResult {
    // One entry per requested index, in request order.
    std::vector<double> values;
    // vectors[k] is the eigenvector for values[k], length n.
    std::vector<std::vector<double>> vectors;
};

// The eigenvectors hicPCA writes out: for `which` = {1, 2} under Dgeev, the
// real parts of columns 0 and 1 of dgeev's right eigenvector matrix.
//
// `covariance` is consumed: both solvers overwrite their input, and the block
// is released before the result is returned so that the peak holds one dense
// block plus dgeev's eigenvector matrix rather than three.
//
// A requested index outside [1, n] yields an empty vector for that entry,
// which is what numpy's out-of-range column slice produces and what makes
// hicPCA skip the bin when writing (hicPCA.py:364, :398).
[[nodiscard]] EigenResult leading_eigenvectors(DenseSymmetric& covariance,
                                               const std::vector<int>& which,
                                               EigenSolver solver);

// --------------------------------------------------------------------------
// SIMD dispatch, exposed for the unit test of cpp/OPTIMIZATION.md section 2.

// Finalises one covariance row:
//
//     out[j] = (accumulator[j] - n * (mean_i * means[j])) * inverse
//
// which is the dense part of the identity. The bracketing is not cosmetic. The
// obvious form, (accumulator[j] - (n * mean_i) * means[j]), rounds n * mean_i
// before multiplying by means[j], and (n*m_i)*m_j is not (n*m_j)*m_i in
// floating point, so the result would not be exactly symmetric. It has to be:
// dgeev takes the whole matrix, not a triangle, and an input that is symmetric
// only to within an ulp can come back with complex eigenvalues. Multiplying
// the two means together first makes the expression symmetric by construction,
// because mean_i * means[j] is exactly mean_j * means[i].
//
// Both entry points compute the same expression; the scalar one is the
// reference implementation and is always compiled, the dispatched one picks
// the widest path the CPU supports.
void covariance_row_finalise_scalar(double* out, const double* accumulator,
                                    const double* means, double mean_i,
                                    double observation_count, double inverse,
                                    std::size_t n);
void covariance_row_finalise(double* out, const double* accumulator,
                             const double* means, double mean_i,
                             double observation_count, double inverse, std::size_t n);

// The name of the path covariance_row_finalise dispatches to on this CPU:
// "scalar", "avx2" or "avx512". Reported by the tools' --version so that a
// benchmark says which kernel it measured.
[[nodiscard]] std::string simd_path_name();

}  // namespace hicx

#endif  // HICX_TRANSFORM_OPS_HPP
