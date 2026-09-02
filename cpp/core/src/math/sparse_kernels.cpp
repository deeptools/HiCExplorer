#include "hicx/math/sparse_kernels.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace hicx::kernels {

namespace {

SimdPath detect_simd_path() {
    // HICX_SIMD pins the dispatch. It exists so that the scalar reference can
    // be measured and compared against on a machine that has AVX2, which is
    // what cpp/OPTIMIZATION.md section 6 asks for, and so that a result can be
    // reproduced on a machine with a different instruction set.
    if (const char* requested = std::getenv("HICX_SIMD"); requested != nullptr) {
        const std::string name(requested);
        if (name == "scalar") {
            return SimdPath::Scalar;
        }
    }
#if defined(__x86_64__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) {
        return SimdPath::Avx2;
    }
#endif
    return SimdPath::Scalar;
}

SimdPath& selected_path() {
    static SimdPath path = detect_simd_path();
    return path;
}

// Runs `body(partition_index)` over every partition. The partitions are a
// property of the block, so the work split does not change with the thread
// count; only how many of them run at once does.
template <class Body>
void for_each_partition(std::size_t partitions, int threads, Body&& body) {
    const std::size_t workers =
        std::min<std::size_t>(partitions, threads > 0 ? static_cast<std::size_t>(threads) : 1);
    if (workers <= 1) {
        for (std::size_t p = 0; p < partitions; ++p) {
            body(p);
        }
        return;
    }
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    // Static round robin, so which partition a worker takes is fixed too. It
    // does not affect the result, but it makes a profile reproducible.
    const auto run = [&](std::size_t worker) {
        for (std::size_t p = worker; p < partitions; p += workers) {
            body(p);
        }
    };
    for (std::size_t worker = 1; worker < workers; ++worker) {
        pool.emplace_back(run, worker);
    }
    run(0);
    for (std::thread& thread : pool) {
        thread.join();
    }
}

}  // namespace

SimdPath active_simd_path() { return selected_path(); }

void force_simd_path(SimdPath path) { selected_path() = path; }

namespace {
std::size_t& split_threshold() {
    static std::size_t threshold = kBlockSplitThreshold;
    return threshold;
}
}  // namespace

void set_block_split_threshold(std::size_t entries) { split_threshold() = entries; }

std::size_t block_split_threshold() { return split_threshold(); }

const char* simd_path_name(SimdPath path) {
    return path == SimdPath::Avx2 ? "avx2" : "scalar";
}

// --------------------------------------------------------------------------
// DiagonalBlock

DiagonalBlock::DiagonalBlock(CsrMatrix& matrix)
    : DiagonalBlock(matrix, 0, matrix.rows()) {}

DiagonalBlock::DiagonalBlock(CsrMatrix& matrix, std::int64_t first, std::int64_t last)
    : first_(first),
      last_(last),
      indices_(matrix.indices().data()),
      data_(matrix.mutable_data().data()),
      indptr_(matrix.indptr().data()) {
    build();
}

void DiagonalBlock::build() {
    const std::int64_t rows = size();
    row_begin_.resize(static_cast<std::size_t>(rows));
    row_end_.resize(static_cast<std::size_t>(rows));
    stored_ = 0;
    for (std::int64_t local = 0; local < rows; ++local) {
        const std::int64_t row = first_ + local;
        const std::int64_t begin = indptr_[row];
        const std::int64_t end = indptr_[row + 1];
        // Columns of an upper triangle row are sorted and at least the row
        // index, so only the upper bound has to be searched for.
        const std::int32_t* first_index = indices_ + begin;
        const std::int32_t* last_index = indices_ + end;
        const std::int32_t* stop =
            std::lower_bound(first_index, last_index, static_cast<std::int32_t>(last_));
        row_begin_[static_cast<std::size_t>(local)] = begin;
        row_end_[static_cast<std::size_t>(local)] = begin + (stop - first_index);
        stored_ += static_cast<std::size_t>(stop - first_index);
    }

    // Entry balanced partitions, computed from the stored counts and therefore
    // identical on every run and at every thread count.
    std::size_t partitions = 1;
    if (stored_ >= split_threshold()) {
        partitions = std::min<std::size_t>(kMaxPartitions, static_cast<std::size_t>(rows));
        partitions = std::max<std::size_t>(partitions, 1);
    }
    partition_.assign(partitions + 1, 0);
    partition_.front() = 0;
    partition_.back() = rows;
    if (partitions > 1) {
        std::size_t seen = 0;
        std::size_t next = 1;
        for (std::int64_t local = 0; local < rows && next < partitions; ++local) {
            seen += static_cast<std::size_t>(row_end_[static_cast<std::size_t>(local)] -
                                             row_begin_[static_cast<std::size_t>(local)]);
            while (next < partitions && seen * partitions >= stored_ * next) {
                partition_[next] = local + 1;
                ++next;
            }
        }
        for (; next < partitions; ++next) {
            partition_[next] = rows;
        }
        // Monotone, and no partition may start before the previous one.
        for (std::size_t p = 1; p < partition_.size(); ++p) {
            partition_[p] = std::max(partition_[p], partition_[p - 1]);
        }
    }
}

// --------------------------------------------------------------------------
// marginals

namespace {

// Scatter phase of one partition: every stored entry above the diagonal adds
// its value to the accumulator of its column. Within a partition the rows are
// visited in ascending order, which is the order scipy's COO has.
void marginal_scatter(const DiagonalBlock& block, std::int64_t from, std::int64_t to,
                      double* partial) {
    const std::int32_t* __restrict indices = block.indices();
    const double* __restrict data = block.data();
    const std::int64_t first = block.first();
    for (std::int64_t local = from; local < to; ++local) {
        const std::int64_t begin = block.row_begin(local);
        const std::int64_t end = block.row_end(local);
        for (std::int64_t k = begin; k < end; ++k) {
            const std::int64_t column = indices[k] - first;
            if (column > local) {
                partial[column] += data[k];
            }
        }
    }
}

}  // namespace

void symmetric_marginals(const DiagonalBlock& block, double* out, int threads) {
    const std::int64_t rows = block.size();
    const std::size_t partitions = block.partitions();
    const std::size_t width = static_cast<std::size_t>(rows);

    std::vector<double> scatter(width * partitions, 0.0);
    for_each_partition(partitions, threads, [&](std::size_t p) {
        marginal_scatter(block, block.partition()[p], block.partition()[p + 1],
                         scatter.data() + p * width);
    });

    // Combine in partition order, then continue the same accumulator with the
    // row's own entries. That is exactly coo_matvec's sequence.
    const std::int32_t* indices = block.indices();
    const double* data = block.data();
    (void)indices;
    for_each_partition(partitions, threads, [&](std::size_t p) {
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            double sum = scatter[static_cast<std::size_t>(local)];
            for (std::size_t q = 1; q < partitions; ++q) {
                sum += scatter[q * width + static_cast<std::size_t>(local)];
            }
            const std::int64_t begin = block.row_begin(local);
            const std::int64_t end = block.row_end(local);
            for (std::int64_t k = begin; k < end; ++k) {
                sum += data[k];
            }
            out[local] = sum;
        }
    });
}

// --------------------------------------------------------------------------
// scaling

namespace {

// data[k] *= scale[row]; data[k] *= scale[col]. Two separate multiplications,
// in that order, because iterativeCorrection.py:56-57 does exactly that and a
// single multiplication by the product rounds differently.
double scale_row_scalar(double* __restrict data, const std::int32_t* __restrict indices,
                        std::int64_t begin, std::int64_t end, const double* __restrict scale,
                        std::int64_t first, double row_scale) {
    double largest = 0.0;
    for (std::int64_t k = begin; k < end; ++k) {
        double value = data[k] * row_scale;
        value *= scale[indices[k] - first];
        data[k] = value;
        const double magnitude = std::fabs(value);
        if (magnitude > largest) {
            largest = magnitude;
        }
    }
    return largest;
}

#if defined(__x86_64__)
__attribute__((target("avx2"))) double scale_row_avx2(
    double* __restrict data, const std::int32_t* __restrict indices, std::int64_t begin,
    std::int64_t end, const double* __restrict scale, std::int64_t first, double row_scale) {
    const __m256d row = _mm256_set1_pd(row_scale);
    const __m128i offset = _mm_set1_epi32(static_cast<int>(first));
    __m256d largest = _mm256_setzero_pd();
    const __m256d absolute = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fffffffffffffffLL));
    std::int64_t k = begin;
    for (; k + 4 <= end; k += 4) {
        const __m256d values = _mm256_loadu_pd(data + k);
        __m128i columns = _mm_loadu_si128(reinterpret_cast<const __m128i*>(indices + k));
        columns = _mm_sub_epi32(columns, offset);
        const __m256d gathered = _mm256_i32gather_pd(scale, columns, 8);
        __m256d result = _mm256_mul_pd(values, row);
        result = _mm256_mul_pd(result, gathered);
        _mm256_storeu_pd(data + k, result);
        largest = _mm256_max_pd(largest, _mm256_and_pd(result, absolute));
    }
    alignas(32) double lanes[4];
    _mm256_store_pd(lanes, largest);
    double tail = std::max(std::max(lanes[0], lanes[1]), std::max(lanes[2], lanes[3]));
    const double rest = scale_row_scalar(data, indices, k, end, scale, first, row_scale);
    return std::max(tail, rest);
}
#endif

inline double scale_row(double* data, const std::int32_t* indices, std::int64_t begin,
                        std::int64_t end, const double* scale, std::int64_t first,
                        double row_scale) {
#if defined(__x86_64__)
    if (selected_path() == SimdPath::Avx2) {
        return scale_row_avx2(data, indices, begin, end, scale, first, row_scale);
    }
#endif
    return scale_row_scalar(data, indices, begin, end, scale, first, row_scale);
}

}  // namespace

double scale_rows_and_cols(const DiagonalBlock& block, const double* scale, int threads) {
    const std::size_t partitions = block.partitions();
    std::vector<double> largest(partitions, 0.0);
    for_each_partition(partitions, threads, [&](std::size_t p) {
        double local_largest = 0.0;
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            local_largest = std::max(
                local_largest,
                scale_row(block.data(), block.indices(), block.row_begin(local),
                          block.row_end(local), scale, block.first(),
                          scale[static_cast<std::size_t>(local)]));
        }
        largest[p] = local_largest;
    });
    return *std::max_element(largest.begin(), largest.end());
}

double scale_and_marginals(const DiagonalBlock& block, const double* scale, double* out,
                           int threads) {
    const std::int64_t rows = block.size();
    const std::size_t partitions = block.partitions();
    const std::size_t width = static_cast<std::size_t>(rows);

    // Phase one scales and scatters in the same traversal, so the 12 bytes per
    // stored entry are read and written once per ICE iteration instead of
    // twice.
    std::vector<double> scatter(width * partitions, 0.0);
    std::vector<double> largest(partitions, 0.0);
    for_each_partition(partitions, threads, [&](std::size_t p) {
        double* __restrict partial = scatter.data() + p * width;
        double* __restrict data = block.data();
        const std::int32_t* __restrict indices = block.indices();
        const std::int64_t first = block.first();
        double local_largest = 0.0;
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            const std::int64_t begin = block.row_begin(local);
            const std::int64_t end = block.row_end(local);
            local_largest =
                std::max(local_largest, scale_row(data, indices, begin, end, scale, first,
                                                  scale[static_cast<std::size_t>(local)]));
            for (std::int64_t k = begin; k < end; ++k) {
                const std::int64_t column = indices[k] - first;
                if (column > local) {
                    partial[column] += data[k];
                }
            }
        }
        largest[p] = local_largest;
    });

    const double* data = block.data();
    for_each_partition(partitions, threads, [&](std::size_t p) {
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            double sum = scatter[static_cast<std::size_t>(local)];
            for (std::size_t q = 1; q < partitions; ++q) {
                sum += scatter[q * width + static_cast<std::size_t>(local)];
            }
            for (std::int64_t k = block.row_begin(local); k < block.row_end(local); ++k) {
                sum += data[k];
            }
            out[local] = sum;
        }
    });
    return *std::max_element(largest.begin(), largest.end());
}

// --------------------------------------------------------------------------
// symmetric matrix vector product

void symmetric_spmv(const DiagonalBlock& block, const double* x, double* y,
                    double diagonal_addend, int threads) {
    const std::int64_t rows = block.size();
    const std::size_t partitions = block.partitions();
    const std::size_t width = static_cast<std::size_t>(rows);

    std::vector<double> scatter(width * partitions, 0.0);
    for_each_partition(partitions, threads, [&](std::size_t p) {
        double* __restrict partial = scatter.data() + p * width;
        const double* __restrict data = block.data();
        const std::int32_t* __restrict indices = block.indices();
        const std::int64_t first = block.first();
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            const double row_x = x[local];
            double row_sum = 0.0;
            for (std::int64_t k = block.row_begin(local); k < block.row_end(local); ++k) {
                const std::int64_t column = indices[k] - first;
                row_sum += data[k] * x[column];
                if (column > local) {
                    partial[column] += data[k] * row_x;
                }
            }
            // The row's own contribution goes into the private buffer as well,
            // so that the combine below is a single sequence per row.
            partial[local] += row_sum;
        }
    });

    for_each_partition(partitions, threads, [&](std::size_t p) {
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            double sum = scatter[static_cast<std::size_t>(local)];
            for (std::size_t q = 1; q < partitions; ++q) {
                sum += scatter[q * width + static_cast<std::size_t>(local)];
            }
            y[local] = sum + diagonal_addend * x[local];
        }
    });
}

// --------------------------------------------------------------------------
// symmetric elementwise scaling

namespace {

void scale_symmetric_row_scalar(double* __restrict data, const std::int32_t* __restrict indices,
                                std::int64_t begin, std::int64_t end,
                                const double* __restrict x, std::int64_t first,
                                double row_x) {
    for (std::int64_t k = begin; k < end; ++k) {
        data[k] = data[k] * row_x * x[indices[k] - first];
    }
}

#if defined(__x86_64__)
__attribute__((target("avx2"))) void scale_symmetric_row_avx2(
    double* __restrict data, const std::int32_t* __restrict indices, std::int64_t begin,
    std::int64_t end, const double* __restrict x, std::int64_t first, double row_x) {
    const __m256d row = _mm256_set1_pd(row_x);
    const __m128i offset = _mm_set1_epi32(static_cast<int>(first));
    std::int64_t k = begin;
    for (; k + 4 <= end; k += 4) {
        const __m256d values = _mm256_loadu_pd(data + k);
        __m128i columns = _mm_loadu_si128(reinterpret_cast<const __m128i*>(indices + k));
        columns = _mm_sub_epi32(columns, offset);
        const __m256d gathered = _mm256_i32gather_pd(x, columns, 8);
        _mm256_storeu_pd(data + k, _mm256_mul_pd(_mm256_mul_pd(values, row), gathered));
    }
    scale_symmetric_row_scalar(data, indices, k, end, x, first, row_x);
}
#endif

}  // namespace

void scale_symmetric(const DiagonalBlock& block, const double* x, int threads) {
    for_each_partition(block.partitions(), threads, [&](std::size_t p) {
        for (std::int64_t local = block.partition()[p]; local < block.partition()[p + 1];
             ++local) {
            const double row_x = x[local];
#if defined(__x86_64__)
            if (selected_path() == SimdPath::Avx2) {
                scale_symmetric_row_avx2(block.data(), block.indices(), block.row_begin(local),
                                         block.row_end(local), x, block.first(), row_x);
                continue;
            }
#endif
            scale_symmetric_row_scalar(block.data(), block.indices(), block.row_begin(local),
                                       block.row_end(local), x, block.first(), row_x);
        }
    });
}

}  // namespace hicx::kernels
