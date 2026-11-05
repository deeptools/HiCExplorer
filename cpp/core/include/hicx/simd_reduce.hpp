// numpy's pairwise reduction, vectorised, with a scalar reference.
//
// hicx::npy::pairwise_sum reproduces np.add.reduce bit for bit: an 8192
// element ufunc buffer, and inside it eight accumulators over blocks of 128
// with a recursive split above that (cpp/PLAN.md 5.2). That structure is not
// an obstacle to vectorising, it is an invitation: eight independent float64
// accumulators are exactly two AVX2 registers or one AVX-512 register, and
// numpy's own SIMD loops use the same eight. Widening the accumulators
// therefore does not move a single value into a different accumulator, so
//
//     the vectorised result is bit-identical to the scalar one, on every input
//
// which is what the unit test asserts. That matters more here than the ED
// tolerance would suggest: hicFindTADs is byte-identical to the Python today,
// and a reduction that merely agreed to 1e-3 would throw that away.
//
// Dispatch follows cpp/OPTIMIZATION.md 2: the library is built at the plain
// x86-64 baseline, the AVX2 and AVX-512 kernels are compiled with function
// target attributes rather than a global -march, and the choice is made once
// per process with __builtin_cpu_supports. There is no -march=native anywhere.

#ifndef HICX_SIMD_REDUCE_HPP
#define HICX_SIMD_REDUCE_HPP

#include <cstddef>
#include <vector>

namespace hicx::simd {

// The scalar reference. Always compiled, always correct; the SIMD paths are
// optimisations of exactly this function.
[[nodiscard]] double pairwise_sum_scalar(const double* data, std::size_t n);

// The dispatched version. Identical results, chosen once per process.
[[nodiscard]] double pairwise_sum(const double* data, std::size_t n);

inline double pairwise_sum(const std::vector<double>& data) {
    return pairwise_sum(data.data(), data.size());
}

// Which kernel the dispatcher picked: "scalar" or "avx2".
//
// There is no AVX-512 kernel. One was written and measured against the AVX2
// one on the development machine, where AVX-512 runs on a 256 bit datapath: it
// was faster only at 64 elements and 24 percent slower at a megabyte, so
// cpp/OPTIMIZATION.md 2 and 6 say not to ship it. The measured table is in
// core/src/simd_reduce.cpp.
[[nodiscard]] const char* active_kernel();

// The individual kernels, for the unit test that compares them. Calling one
// the running CPU does not support is undefined, so the test guards with
// these predicates.
[[nodiscard]] bool avx2_available();
[[nodiscard]] double pairwise_sum_avx2(const double* data, std::size_t n);

}  // namespace hicx::simd

#endif  // HICX_SIMD_REDUCE_HPP
