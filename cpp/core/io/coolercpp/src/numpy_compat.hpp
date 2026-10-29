// Internal: bit exact reimplementations of the numpy and CPython behaviour
// that cooler's stored results and query semantics depend on.
//
// pairwise_sum and float_repr are derived from HiCExplorer v4
// cpp/core/src/numpy_compat.cpp.

#ifndef COOLERCPP_NUMPY_COMPAT_HPP
#define COOLERCPP_NUMPY_COMPAT_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace coolercpp::npy {

// np.add.reduce over a contiguous array: 8192 element buffers, each reduced
// with numpy's pairwise kernel (8 accumulators, blocks of 128), accumulated
// sequentially.
[[nodiscard]] double pairwise_sum(const double* data, std::size_t n);
[[nodiscard]] float pairwise_sum(const float* data, std::size_t n);

// repr() of a Python float.
[[nodiscard]] std::string float_repr(double value);
// str() of a numpy.float32 scalar: the shortest float32 round trip digits in
// the same layout.
[[nodiscard]] std::string float32_repr(float value);

// np.linspace(lo, hi, num, dtype=int) as numpy 1.26 computes it.
[[nodiscard]] std::vector<std::int64_t> linspace_int(std::int64_t lo, std::int64_t hi,
                                                     std::int64_t num);

// A resolved half open range [start, stop) with stop >= start.
struct SliceRange {
    std::int64_t start = 0;
    std::int64_t stop = 0;
    [[nodiscard]] std::int64_t size() const noexcept { return stop - start; }
};

// slice(start, stop).indices(length) for step 1, with stop raised to start
// when it lies before it. This is how numpy and h5py resolve a[start:stop].
[[nodiscard]] SliceRange slice_indices(std::optional<std::int64_t> start,
                                       std::optional<std::int64_t> stop,
                                       std::int64_t length) noexcept;

// numpy scalar indexing a[i]: negative indices count from the end, anything
// still outside [0, length) raises IndexError.
[[nodiscard]] std::int64_t wrap_index(std::int64_t i, std::int64_t length);

// datetime.now().isoformat(): local time, microseconds omitted when zero.
[[nodiscard]] std::string iso_now();

}  // namespace coolercpp::npy

#endif  // COOLERCPP_NUMPY_COMPAT_HPP
