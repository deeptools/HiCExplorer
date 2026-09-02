#include "hicx/transform_ops.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/lapack_shim.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define HICX_X86 1
#endif

namespace hicx {

namespace {

// The stored entries of a matrix in CSR order, as (row, column, value). Every
// kernel below walks the matrix exactly once through this shape, which is
// cpp/OPTIMIZATION.md section 5's flattened iteration: the row index is
// recovered from the row pointer rather than from a nested loop with a
// variable trip count.
template <class F>
void for_each_entry(const CsrMatrix& matrix, F&& visit) {
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& data = matrix.data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            visit(row, static_cast<std::int64_t>(indices[k]), data[k], k);
        }
    }
}

void require_full_square(const CsrMatrix& matrix, const char* what) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error(std::string(what) + ": matrix is not square");
    }
    if (matrix.symmetry() != Symmetry::Full) {
        throw std::runtime_error(std::string(what) +
                                 ": needs the explicit symmetric matrix, call "
                                 "materialize_full() first");
    }
}

// float32 rounding, which every obs/exp kernel applies to the counts before
// dividing (utilities.py:503, :533, :583).
inline double to_float32(double value) {
    return static_cast<double>(static_cast<float>(value));
}

// utilities.convertInfsToZeros_ArrayFloat with pToEpsilon=False: NaN first,
// then infinity, both to zero.
inline double zero_non_finite(double value) {
    return std::isfinite(value) ? value : 0.0;
}

// The dtype the obs/exp kernels cast back to. numpy takes
// `type(pSubmatrix.data[0])` before the float32 step, so a float64 matrix
// stays float64 and an integer matrix is truncated back to an integer.
double cast_back(double value, DType kind) {
    switch (kind) {
        case DType::Integer:
            // numpy's float64 -> int32 cast truncates towards zero. A NaN or
            // an infinity is undefined there; both have already been mapped to
            // zero by the caller.
            return std::trunc(value);
        case DType::Float32:
            return to_float32(value);
        case DType::Float64:
        default:
            return value;
    }
}

}  // namespace

// --------------------------------------------------------------------------
// hiCMatrix mutators

void keep_only_chromosomes(MatrixData& data, const std::vector<std::string>& chromosomes) {
    const std::vector<std::pair<std::string, BinRange>> boundaries =
        chrom_bin_boundaries(data.cut_intervals);
    for (const std::string& name : chromosomes) {
        bool found = false;
        for (const std::pair<std::string, BinRange>& entry : boundaries) {
            if (entry.first == name) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("Chromosome name not in matrix. '" + name + "'");
        }
    }

    // The Python builds a boolean mask and takes np.flatnonzero of it, so the
    // selection is in ascending bin order whatever order the names were given
    // in. That is also what keeps the upper triangle representation valid.
    std::vector<bool> selected(data.cut_intervals.size(), false);
    for (const std::pair<std::string, BinRange>& entry : boundaries) {
        if (std::find(chromosomes.begin(), chromosomes.end(), entry.first) ==
            chromosomes.end()) {
            continue;
        }
        for (std::int64_t bin = entry.second.first; bin < entry.second.last; ++bin) {
            selected[static_cast<std::size_t>(bin)] = true;
        }
    }
    std::vector<std::int64_t> order;
    order.reserve(selected.size());
    for (std::size_t bin = 0; bin < selected.size(); ++bin) {
        if (selected[bin]) {
            order.push_back(static_cast<std::int64_t>(bin));
        }
    }

    data.matrix = select_bins(data.matrix, order);

    std::vector<CutInterval> kept;
    kept.reserve(order.size());
    for (std::int64_t bin : order) {
        kept.push_back(data.cut_intervals[static_cast<std::size_t>(bin)]);
    }
    data.cut_intervals = std::move(kept);

    if (data.correction_factors.has_value()) {
        std::vector<double> factors;
        factors.reserve(order.size());
        for (std::int64_t bin : order) {
            factors.push_back((*data.correction_factors)[static_cast<std::size_t>(bin)]);
        }
        data.correction_factors = std::move(factors);
    }

    if (!data.nan_bins.empty()) {
        std::vector<std::int64_t> mapped;
        std::vector<std::int64_t> position(selected.size(), -1);
        for (std::size_t k = 0; k < order.size(); ++k) {
            position[static_cast<std::size_t>(order[k])] = static_cast<std::int64_t>(k);
        }
        for (std::int64_t bin : data.nan_bins) {
            if (bin >= 0 && static_cast<std::size_t>(bin) < position.size() &&
                position[static_cast<std::size_t>(bin)] >= 0) {
                mapped.push_back(position[static_cast<std::size_t>(bin)]);
            }
        }
        std::sort(mapped.begin(), mapped.end());
        data.nan_bins = std::move(mapped);
    }

    // HiCMatrix.py:677 clears distance_counts unconditionally.
    data.distance_counts.reset();
}

// --------------------------------------------------------------------------
// expected interactions

std::vector<double> expected_interactions_in_distance(const CsrMatrix& matrix,
                                                      std::int64_t length_chromosome,
                                                      std::int64_t chromosome_count) {
    const std::size_t n = static_cast<std::size_t>(matrix.rows());
    std::vector<double> expected(n, 0.0);
    // The Python loops over the stored values in CSR order and accumulates
    // sequentially, so the reduction order is fixed by the storage and is
    // reproduced literally rather than replaced by a pairwise sum.
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t) {
        if (value == 0.0) {
            return;  // csr.nonzero() drops explicitly stored zeros
        }
        expected[static_cast<std::size_t>(std::llabs(row - col))] += value;
    });

    // count_times_i = -( arange(n) * int32(chromosome_count) - int32(length_chromosome) )
    // computed in float64 throughout, as np.arange(float(n)) is float64.
    const double count = static_cast<double>(static_cast<std::int32_t>(chromosome_count));
    const double length = static_cast<double>(static_cast<std::int32_t>(length_chromosome));
    for (std::size_t i = 0; i < n; ++i) {
        const double divisor = -((static_cast<double>(i) * count) - length);
        expected[i] = expected[i] / divisor;
    }
    // No NaN or infinity cleanup here: utilities.py:293-314 has none, unlike
    // its two siblings, and the infinities it can produce are what the obs/exp
    // caller turns into zeros afterwards.
    return expected;
}

std::vector<double> expected_interactions_non_zero(const CsrMatrix& matrix) {
    const std::size_t n = static_cast<std::size_t>(matrix.rows());
    std::vector<double> expected(n, 0.0);
    std::vector<double> occurrences(n, 0.0);
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t) {
        if (value == 0.0) {
            return;
        }
        const std::size_t distance = static_cast<std::size_t>(std::llabs(row - col));
        expected[distance] += value;
        occurrences[distance] += 1.0;
    });
    for (std::size_t i = 0; i < n; ++i) {
        expected[i] = zero_non_finite(expected[i] / occurrences[i]);
    }
    return expected;
}

std::vector<double> expected_interactions(const CsrMatrix& matrix) {
    const std::size_t n = static_cast<std::size_t>(matrix.rows());
    std::vector<double> expected(n, 0.0);
    if (n == 0) {
        return expected;
    }

    // The Python groups the stored values by distance and reduces each group
    // with np.sum, that is with numpy's pairwise summation over the values in
    // storage order. Collecting the groups first is what makes that order
    // reproducible.
    std::vector<std::size_t> counts(n, 0);
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t) {
        if (value == 0.0) {
            return;
        }
        counts[static_cast<std::size_t>(std::llabs(row - col))] += 1;
    });
    std::vector<std::size_t> offsets(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i) {
        offsets[i + 1] = offsets[i] + counts[i];
    }
    std::vector<double> grouped(offsets[n]);
    std::vector<std::size_t> cursor(offsets.begin(), offsets.end() - 1);
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t) {
        if (value == 0.0) {
            return;
        }
        const std::size_t distance = static_cast<std::size_t>(std::llabs(row - col));
        grouped[cursor[distance]++] = value;
    });

    // occurrences = np.arange(n + 1, 1, -1), that is n + 1 - i. Only the
    // distances between the smallest and the largest observed one are filled
    // in, exactly as the Python loop does; the rest stay zero and divide to
    // zero.
    for (std::size_t i = 0; i < n; ++i) {
        const double divisor = static_cast<double>(n + 1 - i);
        const double sum = counts[i] == 0
                               ? 0.0
                               : npy::pairwise_sum(grouped.data() + offsets[i], counts[i]);
        expected[i] = zero_non_finite(sum / divisor);
    }
    return expected;
}

// --------------------------------------------------------------------------
// obs/exp

void obs_exp_lieberman_in_place(CsrMatrix& matrix, std::int64_t length_chromosome,
                                std::int64_t chromosome_count) {
    if (matrix.stored_nnz() == 0) {
        return;
    }
    const std::vector<double> expected =
        expected_interactions_in_distance(matrix, length_chromosome, chromosome_count);
    const DType kind = matrix.dtype_kind();
    std::vector<double>& data = matrix.mutable_data();
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t k) {
        // distance = ceil(abs(row - col) / 2), which is not the index the
        // expected array was accumulated at. utilities.py:499, reproduced.
        const std::size_t distance =
            static_cast<std::size_t>((std::llabs(row - col) + 1) / 2);
        const double divided = to_float32(value) / expected[distance];
        data[k] = cast_back(zero_non_finite(divided), kind);
    });
    // The dtype does not change: astype(data_type) casts back to what it was.
}

void obs_exp_non_zero_in_place(CsrMatrix& matrix, bool ligation_factor) {
    if (matrix.stored_nnz() == 0) {
        return;
    }
    const std::vector<double> expected = expected_interactions_non_zero(matrix);

    // row_sums and total_interactions are taken before the float32 cast, on
    // the matrix as it stands.
    std::vector<double> row_sums(static_cast<std::size_t>(matrix.rows()), 0.0);
    {
        const std::vector<std::int64_t>& indptr = matrix.indptr();
        const std::vector<double>& values = matrix.data();
        for (std::int64_t row = 0; row < matrix.rows(); ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            // csr.sum(axis=1) is a matvec against a vector of ones, which
            // accumulates each row sequentially in storage order.
            double sum = 0.0;
            for (std::size_t k = begin; k < end; ++k) {
                sum += values[k];
            }
            row_sums[static_cast<std::size_t>(row)] = sum;
        }
    }
    const double total = matrix.sum().as_double();

    std::vector<double>& data = matrix.mutable_data();
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t k) {
        double divisor = expected[static_cast<std::size_t>(std::llabs(row - col))];
        if (ligation_factor) {
            divisor *= row_sums[static_cast<std::size_t>(row)] *
                       row_sums[static_cast<std::size_t>(col)] / total;
        }
        // The Python assigns each quotient back into a float32 array, so every
        // value is rounded to single precision one at a time.
        data[k] = to_float32(to_float32(value) / divisor);
    });
    for (double& value : data) {
        if (!std::isfinite(value)) {
            value = 0.0;
        }
    }
    matrix.set_dtype("float32");
    matrix.eliminate_zeros();
}

void obs_exp_in_place(CsrMatrix& matrix) {
    if (matrix.stored_nnz() == 0) {
        return;
    }
    const std::vector<double> expected = expected_interactions(matrix);
    const DType kind = matrix.dtype_kind();
    std::vector<double>& data = matrix.mutable_data();
    for_each_entry(matrix, [&](std::int64_t row, std::int64_t col, double value, std::size_t k) {
        const std::size_t distance =
            static_cast<std::size_t>((std::llabs(row - col) + 1) / 2);
        const double divided = to_float32(value) / expected[distance];
        data[k] = cast_back(zero_non_finite(divided), kind);
    });
}

void convert_nans_and_infs_to_zeros(CsrMatrix& matrix) {
    for (double& value : matrix.mutable_data()) {
        if (!std::isfinite(value)) {
            value = 0.0;
        }
    }
}

// --------------------------------------------------------------------------
// the dense per-row finalisation, scalar reference and SIMD variants

void covariance_row_finalise_scalar(double* out, const double* accumulator,
                                    const double* means, double mean_i,
                                    double observation_count, double inverse,
                                    std::size_t n) {
    for (std::size_t j = 0; j < n; ++j) {
        out[j] = (accumulator[j] - observation_count * (mean_i * means[j])) * inverse;
    }
}

#if defined(HICX_X86)

__attribute__((target("avx2,fma"))) static void covariance_row_finalise_avx2(
    double* out, const double* accumulator, const double* means, double mean_i,
    double observation_count, double inverse, std::size_t n) {
    const __m256d mean_i_vector = _mm256_set1_pd(mean_i);
    const __m256d count = _mm256_set1_pd(observation_count);
    const __m256d scale = _mm256_set1_pd(inverse);
    std::size_t j = 0;
    for (; j + 4 <= n; j += 4) {
        const __m256d a = _mm256_loadu_pd(accumulator + j);
        const __m256d m = _mm256_loadu_pd(means + j);
        // The product of the two means first, so that the expression stays
        // symmetric, then one fused multiply-add. fnmadd rounds once where the
        // scalar expression rounds twice, so the two paths are not bit
        // identical; cpp/OPTIMIZATION.md section 3 asks for agreement at the
        // ED tolerance, not for equality, for exactly this reason.
        const __m256d product = _mm256_mul_pd(mean_i_vector, m);
        const __m256d centred = _mm256_fnmadd_pd(count, product, a);
        _mm256_storeu_pd(out + j, _mm256_mul_pd(centred, scale));
    }
    for (; j < n; ++j) {
        out[j] = (accumulator[j] - observation_count * (mean_i * means[j])) * inverse;
    }
}

__attribute__((target("avx512f"))) static void covariance_row_finalise_avx512(
    double* out, const double* accumulator, const double* means, double mean_i,
    double observation_count, double inverse, std::size_t n) {
    const __m512d mean_i_vector = _mm512_set1_pd(mean_i);
    const __m512d count = _mm512_set1_pd(observation_count);
    const __m512d scale = _mm512_set1_pd(inverse);
    std::size_t j = 0;
    for (; j + 8 <= n; j += 8) {
        const __m512d a = _mm512_loadu_pd(accumulator + j);
        const __m512d m = _mm512_loadu_pd(means + j);
        const __m512d product = _mm512_mul_pd(mean_i_vector, m);
        const __m512d centred = _mm512_fnmadd_pd(count, product, a);
        _mm512_storeu_pd(out + j, _mm512_mul_pd(centred, scale));
    }
    for (; j < n; ++j) {
        out[j] = (accumulator[j] - observation_count * (mean_i * means[j])) * inverse;
    }
}

namespace {

enum class SimdPath { Scalar, Avx2, Avx512 };

SimdPath detect_simd_path() {
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f")) {
        return SimdPath::Avx512;
    }
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) {
        return SimdPath::Avx2;
    }
    return SimdPath::Scalar;
}

// Resolved once per process, so the dispatch choice cannot vary within a run,
// which cpp/OPTIMIZATION.md section 3 requires.
const SimdPath kSimdPath = detect_simd_path();

}  // namespace

void covariance_row_finalise(double* out, const double* accumulator, const double* means,
                             double mean_i, double observation_count, double inverse,
                             std::size_t n) {
    switch (kSimdPath) {
        case SimdPath::Avx512:
            covariance_row_finalise_avx512(out, accumulator, means, mean_i,
                                           observation_count, inverse, n);
            return;
        case SimdPath::Avx2:
            covariance_row_finalise_avx2(out, accumulator, means, mean_i, observation_count,
                                         inverse, n);
            return;
        case SimdPath::Scalar:
        default:
            covariance_row_finalise_scalar(out, accumulator, means, mean_i,
                                           observation_count, inverse, n);
            return;
    }
}

std::string simd_path_name() {
    switch (kSimdPath) {
        case SimdPath::Avx512:
            return "avx512";
        case SimdPath::Avx2:
            return "avx2";
        default:
            return "scalar";
    }
}

#else   // !HICX_X86

void covariance_row_finalise(double* out, const double* accumulator, const double* means,
                             double mean_i, double observation_count, double inverse,
                             std::size_t n) {
    covariance_row_finalise_scalar(out, accumulator, means, mean_i, observation_count,
                                   inverse, n);
}

std::string simd_path_name() { return "scalar"; }

#endif  // HICX_X86

// --------------------------------------------------------------------------
// covariance

DenseSymmetric covariance_of_symmetric(const CsrMatrix& matrix, int threads) {
    require_full_square(matrix, "covariance_of_symmetric");
    const std::int64_t n = matrix.rows();
    // Not zero filled: every element is written by compute_rows below, which
    // memsets the row it is about to accumulate into. Zero filling here would
    // touch the whole block a second time and serially; see
    // DenseSymmetric::uninitialized.
    DenseSymmetric result = DenseSymmetric::uninitialized(n);
    if (n == 0) {
        return result;
    }
    if (n < 2) {
        // np.cov with a single observation divides by zero and yields NaN,
        // which hicPCA maps to zero.
        result.at(0, 0) = 0.0;
        return result;
    }

    // m = A.mean(axis=1). The row sums go through numpy's pairwise reduction
    // so that they match the dense mean to the last bits the sparse layout
    // allows.
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    const std::vector<double>& data = matrix.data();
    std::vector<double> means(static_cast<std::size_t>(n));
    const double dn = static_cast<double>(n);
    for (std::int64_t i = 0; i < n; ++i) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i)]);
        const std::size_t end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(i) + 1]);
        means[static_cast<std::size_t>(i)] =
            npy::pairwise_sum(data.data() + begin, end - begin) / dn;
    }

    const double inverse = 1.0 / (dn - 1.0);

    // One task per contiguous row range. Rows write into disjoint memory and
    // each row's reduction runs sequentially in a fixed index order, so the
    // result does not depend on the number of threads
    // (cpp/OPTIMIZATION.md section 3).
    const int worker_count = std::max(1, threads);
    auto compute_rows = [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            double* out = result.row(i);
            std::memset(out, 0, static_cast<std::size_t>(n) * sizeof(double));
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(i)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(i) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                const double a = data[k];
                if (a == 0.0) {
                    continue;
                }
                const std::int64_t middle = indices[k];
                const std::size_t inner_begin =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(middle)]);
                const std::size_t inner_end =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(middle) + 1]);
                const std::int32_t* __restrict inner_cols = indices.data() + inner_begin;
                const double* __restrict inner_vals = data.data() + inner_begin;
                const std::size_t count = inner_end - inner_begin;
                for (std::size_t t = 0; t < count; ++t) {
                    out[inner_cols[t]] += a * inner_vals[t];
                }
            }
            const double mean_i = means[static_cast<std::size_t>(i)];
            covariance_row_finalise(out, out, means.data(), mean_i, dn, inverse,
                                    static_cast<std::size_t>(n));
        }
    };

    if (worker_count == 1 || n < 64) {
        compute_rows(0, n);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));
        const std::int64_t chunk = (n + worker_count - 1) / worker_count;
        for (int w = 0; w < worker_count; ++w) {
            const std::int64_t first = static_cast<std::int64_t>(w) * chunk;
            const std::int64_t last = std::min(n, first + chunk);
            if (first >= last) {
                break;
            }
            workers.emplace_back(compute_rows, first, last);
        }
        for (std::thread& worker : workers) {
            worker.join();
        }
    }
    return result;
}

void zero_non_finite_in_place(DenseSymmetric& block, int threads) {
    const std::int64_t n = block.size();
    if (n == 0) {
        return;
    }
    const int worker_count = std::max(1, threads);
    auto rows = [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            double* values = block.row(i);
            for (std::int64_t j = 0; j < n; ++j) {
                if (!std::isfinite(values[j])) {
                    values[j] = 0.0;
                }
            }
        }
    };
    if (worker_count == 1 || n < 64) {
        rows(0, n);
        return;
    }
    std::vector<std::thread> workers;
    const std::int64_t chunk = (n + worker_count - 1) / worker_count;
    for (int w = 0; w < worker_count; ++w) {
        const std::int64_t first = static_cast<std::int64_t>(w) * chunk;
        const std::int64_t last = std::min(n, first + chunk);
        if (first >= last) {
            break;
        }
        workers.emplace_back(rows, first, last);
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
}

std::vector<double> pearson_scaling(const DenseSymmetric& covariance) {
    const std::int64_t n = covariance.size();
    std::vector<double> scaling(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) {
        scaling[static_cast<std::size_t>(i)] = std::sqrt(covariance.at(i, i));
    }
    return scaling;
}

void pearson_row(const DenseSymmetric& covariance, const std::vector<double>& scaling,
                 std::int64_t i, double* out) {
    const std::int64_t n = covariance.size();
    const double* source = covariance.row(i);
    const double di = scaling[static_cast<std::size_t>(i)];
    for (std::int64_t j = 0; j < n; ++j) {
        // numpy.corrcoef divides twice, by d[:, None] and then by d[None, :],
        // and clips the result into [-1, 1]. Two divisions, not one multiply
        // by the reciprocal of the product.
        double value = source[j] / di;
        value = value / scaling[static_cast<std::size_t>(j)];
        if (value < -1.0) {
            value = -1.0;
        } else if (value > 1.0) {
            value = 1.0;
        }
        // np.clip leaves a NaN alone; hicPCA's convertNansToZeros then maps it
        // to zero, and no infinity can survive the clip.
        out[j] = std::isnan(value) ? 0.0 : value;
    }
}

void covariance_to_pearson_in_place(DenseSymmetric& covariance, int threads) {
    const std::int64_t n = covariance.size();
    if (n == 0) {
        return;
    }
    const std::vector<double> scaling = pearson_scaling(covariance);
    const int worker_count = std::max(1, threads);
    auto rows = [&](std::int64_t first, std::int64_t last) {
        for (std::int64_t i = first; i < last; ++i) {
            pearson_row(covariance, scaling, i, covariance.row(i));
        }
    };
    if (worker_count == 1 || n < 64) {
        rows(0, n);
        return;
    }
    std::vector<std::thread> workers;
    const std::int64_t chunk = (n + worker_count - 1) / worker_count;
    for (int w = 0; w < worker_count; ++w) {
        const std::int64_t first = static_cast<std::int64_t>(w) * chunk;
        const std::int64_t last = std::min(n, first + chunk);
        if (first >= last) {
            break;
        }
        workers.emplace_back(rows, first, last);
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
}

// --------------------------------------------------------------------------
// eigenvectors

namespace {

EigenResult eigenvectors_dgeev(DenseSymmetric& covariance, const std::vector<int>& which) {
    const int n = static_cast<int>(covariance.size());
    EigenResult result;
    result.values.assign(which.size(), 0.0);
    result.vectors.assign(which.size(), {});
    if (n == 0) {
        return result;
    }

    std::vector<double> wr(static_cast<std::size_t>(n));
    std::vector<double> wi(static_cast<std::size_t>(n));
    std::vector<double> vr(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));

    const char jobvl = 'N';
    const char jobvr = 'V';
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    double* vl = nullptr;
    const int ldvl = 1;
    dgeev_(&jobvl, &jobvr, &n, covariance.data(), &n, wr.data(), wi.data(), vl, &ldvl,
           vr.data(), &n, &work_query, &lwork, &info);
    if (info != 0) {
        throw std::runtime_error("dgeev workspace query failed");
    }
    lwork = static_cast<int>(work_query);
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    dgeev_(&jobvl, &jobvr, &n, covariance.data(), &n, wr.data(), wi.data(), vl, &ldvl,
           vr.data(), &n, work.data(), &lwork, &info);
    if (info != 0) {
        throw std::runtime_error("dgeev failed with info " + std::to_string(info));
    }
    // The input has been destroyed; release it before copying out, so that the
    // peak is one dense block plus vr rather than two plus vr.
    covariance.release();

    for (std::size_t entry = 0; entry < which.size(); ++entry) {
        const int index = which[entry] - 1;
        if (index < 0 || index >= n) {
            continue;  // numpy's out of range column slice is empty
        }
        // A complex conjugate pair occupies columns (j, j+1) of vr as (real,
        // imaginary). scipy assembles v_j = vr_j + i vr_{j+1} and
        // v_{j+1} = vr_j - i vr_{j+1}, and hicPCA writes only .real, so both
        // members of the pair yield column j.
        int source = index;
        if (wi[static_cast<std::size_t>(index)] < 0.0 && index > 0) {
            source = index - 1;
        }
        const double* column = vr.data() + static_cast<std::size_t>(source) *
                                               static_cast<std::size_t>(n);
        result.values[entry] = wr[static_cast<std::size_t>(index)];
        result.vectors[entry].assign(column, column + n);
    }
    return result;
}

EigenResult eigenvectors_dsyevr(DenseSymmetric& covariance, const std::vector<int>& which) {
    const int n = static_cast<int>(covariance.size());
    EigenResult result;
    result.values.assign(which.size(), 0.0);
    result.vectors.assign(which.size(), {});
    if (n == 0) {
        return result;
    }

    int highest = 0;
    for (int index : which) {
        highest = std::max(highest, index);
    }
    if (highest <= 0) {
        return result;
    }
    const int wanted = std::min(highest, n);

    const char jobz = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const double vl_unused = 0.0;
    const double vu_unused = 0.0;
    const int il = n - wanted + 1;
    const int iu = n;
    const char safe_minimum = 'S';
    const double abstol = 2.0 * dlamch_(&safe_minimum);
    int found = 0;
    std::vector<double> values(static_cast<std::size_t>(n));
    std::vector<double> vectors(static_cast<std::size_t>(n) * static_cast<std::size_t>(wanted));
    std::vector<int> isuppz(static_cast<std::size_t>(2 * std::max(wanted, 1)));

    int info = 0;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    dsyevr_(&jobz, &range, &uplo, &n, covariance.data(), &n, &vl_unused, &vu_unused, &il, &iu,
            &abstol, &found, values.data(), vectors.data(), &n, isuppz.data(), &work_query,
            &lwork, &iwork_query, &liwork, &info);
    if (info != 0) {
        throw std::runtime_error("dsyevr workspace query failed");
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
    std::vector<int> iwork(static_cast<std::size_t>(std::max(liwork, 1)));
    dsyevr_(&jobz, &range, &uplo, &n, covariance.data(), &n, &vl_unused, &vu_unused, &il, &iu,
            &abstol, &found, values.data(), vectors.data(), &n, isuppz.data(), work.data(),
            &lwork, iwork.data(), &liwork, &info);
    if (info != 0) {
        throw std::runtime_error("dsyevr failed with info " + std::to_string(info));
    }
    covariance.release();

    for (std::size_t entry = 0; entry < which.size(); ++entry) {
        const int index = which[entry] - 1;
        if (index < 0 || index >= n || index >= found) {
            continue;
        }
        // dsyevr returns the selected eigenvalues in ascending order; the
        // requested index counts from the largest.
        const int ascending = found - 1 - index;
        const double* column = vectors.data() + static_cast<std::size_t>(ascending) *
                                                    static_cast<std::size_t>(n);
        std::vector<double> vector(column, column + n);
        // A deterministic sign: the largest magnitude component is made
        // positive, with the first such component deciding a tie. LAPACK does
        // not fix the sign and this is the only stable convention available
        // without an external track.
        std::size_t extreme = 0;
        double best = -1.0;
        for (std::size_t j = 0; j < vector.size(); ++j) {
            const double magnitude = std::fabs(vector[j]);
            if (magnitude > best) {
                best = magnitude;
                extreme = j;
            }
        }
        if (vector[extreme] < 0.0) {
            for (double& value : vector) {
                value = -value;
            }
        }
        result.values[entry] = values[static_cast<std::size_t>(ascending)];
        result.vectors[entry] = std::move(vector);
    }
    return result;
}

}  // namespace

void pin_blas_to_one_thread() {
    if (openblas_set_num_threads != nullptr) {
        openblas_set_num_threads(1);
    }
}

EigenResult leading_eigenvectors(DenseSymmetric& covariance, const std::vector<int>& which,
                                 EigenSolver solver) {
    return solver == EigenSolver::Dgeev ? eigenvectors_dgeev(covariance, which)
                                        : eigenvectors_dsyevr(covariance, which);
}

}  // namespace hicx
