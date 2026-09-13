// np.argsort and np.argpartition with numpy's own tie order.
//
// np.argsort's default kind, 'quicksort', is not stable, and which of several
// equal keys comes first is observable wherever a tool writes lines in argsort
// order. hicAggregateContacts writes its contact pairs that way, sorted by a
// centre value that is zero for a third of the lines, so the order of the tied
// lines is part of the output.
//
// numpy 1.26.4 does not use one algorithm for this. `aquicksort_double` and
// `aquicksort_longlong` first try `aquicksort_dispatch`, which on a CPU with
// the AVX512_SKX feature group (F, CD, BW, DQ, VL) calls
// `avx512_argsort` from the x86-simd-sort submodule pinned at commit
// 0631a88763a4a0a4c9e84d5eeb0ec5d36053730b; everywhere else it runs numpy's
// introsort, `aquicksort_` (median-of-3 quicksort, insertion sort below 16
// elements, heapsort when the recursion gets too deep). The two produce
// different tie orders, even for ten elements.
//
// This header does the same thing on the same criterion, so that the port
// reproduces numpy on whichever machine it runs on: the AVX-512 path compiles
// the pinned x86-simd-sort header itself (fetched at configure time, in its own
// translation unit built with the SKX flags), and the portable path is a line
// by line port of `aquicksort_` and `aheapsort_`. The dispatch is decided once
// per process with __builtin_cpu_supports, as numpy decides it once at import.
//
// A float64 array containing NaN goes, on the AVX-512 path, through
// x86-simd-sort's `std_argsort_withnan`, a std::sort with NaN ordered last;
// the portable path orders NaN last through numpy's DOUBLE_LT.
//
// np.argpartition has no SIMD dispatch in numpy 1.26 and is a port of
// `introselect_<Tag, arg=true>` from npysort/selection.cpp.

#ifndef HICX_NUMPY_SORT_HPP
#define HICX_NUMPY_SORT_HPP

#include <cstdint>
#include <vector>

namespace hicx::npy {

// np.argsort(values) for a contiguous float64 array.
[[nodiscard]] std::vector<std::int64_t> argsort(const std::vector<double>& values);
// np.argsort(values) for a contiguous int64 array.
[[nodiscard]] std::vector<std::int64_t> argsort(const std::vector<std::int64_t>& values);

// np.argpartition(values, kth) for a contiguous float64 array; a negative kth
// counts from the end, as numpy's does.
[[nodiscard]] std::vector<std::int64_t> argpartition(const std::vector<double>& values,
                                                     std::int64_t kth);

// Which argsort numpy would run here: "avx512_skx" or "introsort".
[[nodiscard]] const char* argsort_dispatch_name();

namespace detail {
// The two argsort implementations, exposed so that the unit tests can check
// both on every machine that can run them. The AVX-512 functions must only be
// called when avx512_skx_available() is true.
[[nodiscard]] bool avx512_skx_available();
void argsort_avx512_skx(const double* values, std::int64_t* index, std::int64_t n);
void argsort_avx512_skx(const std::int64_t* values, std::int64_t* index, std::int64_t n);
void argsort_introsort(const double* values, std::int64_t* index, std::int64_t n);
void argsort_introsort(const std::int64_t* values, std::int64_t* index, std::int64_t n);
}  // namespace detail

}  // namespace hicx::npy

#endif  // HICX_NUMPY_SORT_HPP
