// Order statistics, quantiles, linspace and searchsorted with the results of
// the numpy functions of the same name, written from their documented
// definitions:
//
//   * median: the middle value, for an even count the mean of the two middle
//     values; NaN when the input is empty or holds a NaN.
//   * quantile, method 'linear' (Hyndman and Fan, "Sample quantiles in
//     statistical packages", The American Statistician 50(4), 1996, definition
//     7): the sorted values x[0..n-1] are interpolated at the virtual index
//     (n - 1) q. The interpolation is evaluated from the nearer of the two
//     neighbours (a + (b - a) t for t < 1/2, b - (b - a)(1 - t) otherwise),
//     which keeps the result monotone in t and exact at both ends.
//   * linspace: start + k * step with step = (stop - start) / (num - 1) and the
//     last sample set to stop.
//   * searchsorted with side='right': the first position whose value is
//     greater than the key, NaN ordering after every number.
//
// Identity with numpy 1.26.4 was verified on randomised inputs including
// ties, NaN, infinities and unsorted search arrays, see
// cpp/tests/test_numpy_oracle.cpp.

#ifndef HICX_NUMPY_STATS_HPP
#define HICX_NUMPY_STATS_HPP

#include <cstdint>
#include <vector>

namespace hicx::npy {

// np.median(values) for float64 values.
[[nodiscard]] double median(std::vector<double> values);

// np.quantile(sorted, q) with method='linear' for sorted, NaN free values and
// 0 <= q <= 1. NaN for an empty input.
[[nodiscard]] double quantile_linear_sorted(const std::vector<double>& sorted, double q);

// np.nanquantile(values, quantiles) with method='linear'. Throws
// std::invalid_argument("Quantiles must be in the range [0, 1]") when a
// quantile is outside [0, 1] or NaN. All NaN input gives NaN quantiles.
[[nodiscard]] std::vector<double> nanquantile_linear(std::vector<double> values,
                                                     const std::vector<double>& quantiles);

// np.percentile(values, q) with method='linear': NaN as soon as one value is NaN.
[[nodiscard]] double percentile_linear(std::vector<double> values, double q);

// np.linspace(start, stop, num) for float64 scalars. Throws
// std::invalid_argument for a negative num.
[[nodiscard]] std::vector<double> linspace(double start, double stop, std::int64_t num);

// np.searchsorted(sorted, keys, side='right'). Each search starts from the
// window the previous result allows, so for an array that is not sorted the
// result depends on the order of the keys exactly as numpy's does.
[[nodiscard]] std::vector<std::int64_t> searchsorted_right(const std::vector<double>& sorted,
                                                           const std::vector<double>& keys);

}  // namespace hicx::npy

#endif  // HICX_NUMPY_STATS_HPP
