#include "hicx/simd_reduce.hpp"

#include <algorithm>
#include <array>
#include <cstdint>

#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace hicx::simd {

namespace {

// numpy's constants: PW_BLOCKSIZE from
// numpy/core/src/umath/loops_arithm_fp.dispatch.c.src, and the ufunc buffer
// size from np.getbufsize().
constexpr std::size_t kPairwiseBlockSize = 128;
constexpr std::size_t kReduceBufferSize = 8192;

// The tail every kernel shares: fewer than eight elements are summed left to
// right, exactly as numpy does.
double small_sum(const double* a, std::size_t n) {
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result += a[i];
    }
    return result;
}

// The fixed combination order of the eight accumulators. Every kernel ends
// here, so the grouping cannot drift between them.
double combine(const std::array<double, 8>& r) {
    return ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
}

double block_scalar(const double* a, std::size_t n) {
    if (n < 8) {
        return small_sum(a, n);
    }
    if (n <= kPairwiseBlockSize) {
        std::array<double, 8> r{a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]};
        std::size_t i = 8;
        const std::size_t limit = n - (n % 8);
        for (; i < limit; i += 8) {
            r[0] += a[i + 0];
            r[1] += a[i + 1];
            r[2] += a[i + 2];
            r[3] += a[i + 3];
            r[4] += a[i + 4];
            r[5] += a[i + 5];
            r[6] += a[i + 6];
            r[7] += a[i + 7];
        }
        double result = combine(r);
        for (; i < n; ++i) {
            result += a[i];
        }
        return result;
    }
    std::size_t half = n / 2;
    half -= half % 8;
    return block_scalar(a, half) + block_scalar(a + half, n - half);
}

template <class Block>
double buffered(const double* a, std::size_t n, Block&& block) {
    double accumulated = 0.0;
    std::size_t offset = 0;
    while (offset < n) {
        const std::size_t chunk = std::min(kReduceBufferSize, n - offset);
        accumulated += block(a + offset, chunk);
        offset += chunk;
    }
    return accumulated;
}

#if defined(__x86_64__)

// AVX2: the eight accumulators live in two 256 bit registers, r[0..3] in the
// first and r[4..7] in the second, so element i lands in accumulator i % 8
// exactly as in the scalar loop. The registers are stored back before the
// combination, so the final grouping is the scalar one, not a horizontal
// instruction's.
__attribute__((target("avx2"))) double block_avx2(const double* a, std::size_t n) {
    if (n < 8) {
        return small_sum(a, n);
    }
    if (n <= kPairwiseBlockSize) {
        __m256d lo = _mm256_loadu_pd(a);
        __m256d hi = _mm256_loadu_pd(a + 4);
        std::size_t i = 8;
        const std::size_t limit = n - (n % 8);
        for (; i < limit; i += 8) {
            lo = _mm256_add_pd(lo, _mm256_loadu_pd(a + i));
            hi = _mm256_add_pd(hi, _mm256_loadu_pd(a + i + 4));
        }
        std::array<double, 8> r{};
        _mm256_storeu_pd(r.data(), lo);
        _mm256_storeu_pd(r.data() + 4, hi);
        double result = combine(r);
        for (; i < n; ++i) {
            result += a[i];
        }
        return result;
    }
    std::size_t half = n / 2;
    half -= half % 8;
    return block_avx2(a, half) + block_avx2(a + half, n - half);
}

#endif  // __x86_64__

enum class Kernel { Scalar, Avx2 };

// AVX-512 is deliberately absent. A third kernel with all eight accumulators
// in one zmm register was written and measured against the AVX2 one on the
// development machine (Zen 4, which implements AVX-512 on a 256 bit datapath),
// over array lengths from 64 to 1,048,576 doubles:
//
//     n         scalar      AVX2       AVX-512
//     64        60.4 GB/s   72.7 GB/s  77.3 GB/s
//     256       88.1        133.9      122.9
//     1024      91.0        135.7      126.7
//     8192      84.4        103.2      102.8
//     65536     83.5        104.9      101.1
//     1048576   66.8        80.9       65.1
//
// It is faster than AVX2 only at n = 64, which no caller uses, and 24 percent
// slower at a megabyte. cpp/OPTIMIZATION.md 2 says not to ship a third path
// that is not faster, and 6 says to revert a change that does not measurably
// help, so it was removed rather than left in place.
Kernel detect() {
#if defined(__x86_64__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) {
        return Kernel::Avx2;
    }
#endif
    return Kernel::Scalar;
}

// Decided once, so the dispatch cannot vary within a run
// (cpp/OPTIMIZATION.md 3).
const Kernel kSelected = detect();

}  // namespace

double pairwise_sum_scalar(const double* data, std::size_t n) {
    return buffered(data, n, block_scalar);
}

bool avx2_available() {
#if defined(__x86_64__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0;
#else
    return false;
#endif
}

double pairwise_sum_avx2(const double* data, std::size_t n) {
#if defined(__x86_64__)
    return buffered(data, n, block_avx2);
#else
    return pairwise_sum_scalar(data, n);
#endif
}

double pairwise_sum(const double* data, std::size_t n) {
    switch (kSelected) {
#if defined(__x86_64__)
        case Kernel::Avx2:
            return buffered(data, n, block_avx2);
#endif
        default:
            return buffered(data, n, block_scalar);
    }
}

const char* active_kernel() {
    switch (kSelected) {
        case Kernel::Avx2:
            return "avx2";
        default:
            return "scalar";
    }
}

}  // namespace hicx::simd
