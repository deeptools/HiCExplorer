// np.argsort and np.argpartition with numpy's own tie order.
//
// np.argsort's default kind, 'quicksort', is not stable, and which of several
// equal keys comes first is observable wherever a tool writes lines in argsort
// order. hicAggregateContacts writes its contact pairs that way, sorted by a
// centre value that is zero for a third of the lines, so the order of the tied
// lines is part of the output.
//
// argsort is an introspective quicksort (median-of-3 pivot, insertion sort
// for short ranges, heapsort when the recursion gets too deep) written from
// the published algorithm with the tie-deciding parameters numpy documents. It
// runs on every CPU, so the tie order does not depend on the machine, and it
// is the order numpy produces without its AVX-512 argsort. numpy 1.26 on a CPU with AVX512_SKX sorts with x86-simd-sort
// instead, which orders tied keys differently, so the tie order of the
// reference on such a CPU is not reproduced.
//
// A float64 array containing NaN is ordered with NaN last, through numpy's
// DOUBLE_LT.
//
// np.argpartition is an introspective selection (quickselect with a
// median-of-3 pivot, median of medians as the fall back), see numpy_sort.cpp.

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

namespace detail {
// The introsort argsort, exposed so that the unit tests can call it on an
// index array they prepare.
void argsort_introsort(const double* values, std::int64_t* index, std::int64_t n);
void argsort_introsort(const std::int64_t* values, std::int64_t* index, std::int64_t n);
// Selection of the kth smallest on a prepared index array (0 <= kth < n).
void argpartition_select(const double* values, std::int64_t* index, std::int64_t n,
                         std::int64_t kth);
}  // namespace detail

}  // namespace hicx::npy

#endif  // HICX_NUMPY_SORT_HPP
