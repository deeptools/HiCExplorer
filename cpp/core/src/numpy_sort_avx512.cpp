// The AVX-512 half of np.argsort, see numpy_sort.hpp.
//
// This translation unit is compiled with the AVX512_SKX flags numpy uses for
// its simd_qsort dispatch target and includes the x86-simd-sort header at the
// commit numpy 1.26.4 vendors. Nothing here may run unless
// hicx::npy::detail::avx512_skx_available() said so; the rest of the library
// is built for the x86-64-v2 baseline.

#include <cstdint>

#include <stdexcept>

#include "hicx/numpy_sort.hpp"

#if defined(__x86_64__)
#include "avx512-64bit-argsort.hpp"
#endif

namespace hicx::npy::detail {

#if defined(__x86_64__)
void argsort_avx512_skx(const double* values, std::int64_t* index, std::int64_t n) {
    // numpy's ArgQSort<double> casts away the const; avx512_argsort only reads.
    avx512_argsort(const_cast<double*>(values), index, n);
}

void argsort_avx512_skx(const std::int64_t* values, std::int64_t* index, std::int64_t n) {
    avx512_argsort(const_cast<std::int64_t*>(values), index, n);
}
#else
// Not x86-64: avx512_skx_available() is false, nothing calls these.
void argsort_avx512_skx(const double*, std::int64_t*, std::int64_t) {
    throw std::logic_error("AVX-512 argsort is not available on this architecture");
}

void argsort_avx512_skx(const std::int64_t*, std::int64_t*, std::int64_t) {
    throw std::logic_error("AVX-512 argsort is not available on this architecture");
}
#endif

}  // namespace hicx::npy::detail
