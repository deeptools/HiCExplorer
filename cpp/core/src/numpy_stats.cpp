#include "hicx/numpy_stats.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace hicx::npy {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Ordering of float64 for sorting and searching: NaN after every number.
bool nan_last_less(double a, double b) { return a < b || (b != b && a == a); }

bool has_nan(const std::vector<double>& values) {
    return std::any_of(values.begin(), values.end(), [](double v) { return std::isnan(v); });
}

// a + (b - a) t, evaluated from the nearer end.
double interpolate(double a, double b, double t) {
    const double span = b - a;
    return t >= 0.5 ? b - span * (1.0 - t) : a + span * t;
}

}  // namespace

double median(std::vector<double> values) {
    if (values.empty() || has_nan(values)) {
        return kNaN;
    }
    const std::size_t half = values.size() / 2;
    const auto middle = values.begin() + static_cast<std::ptrdiff_t>(half);
    std::nth_element(values.begin(), middle, values.end());
    const double upper = *middle;
    // The mean of the middle value(s), summed from +0 (so the median of
    // [-0.0] is +0.0, as in numpy).
    if (values.size() % 2 == 1) {
        return 0.0 + upper;
    }
    const double lower = *std::max_element(values.begin(), middle);
    return ((0.0 + lower) + upper) / 2.0;
}

double quantile_linear_sorted(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) {
        return kNaN;
    }
    const std::size_t last = sorted.size() - 1;
    const double virtual_index = static_cast<double>(last) * q;
    if (virtual_index >= static_cast<double>(last)) {
        // Top of the range: both neighbours are the last value. numpy gives the
        // weight one past the virtual index here (virtual index + 1, so a
        // weight above 1/2 and the upper-end formula). The value is the last
        // one either way; the weight decides only the sign of a zero and
        // whether an infinite value turns into NaN, and both are matched.
        return interpolate(sorted[last], sorted[last], virtual_index + 1.0);
    }
    const auto below = static_cast<std::size_t>(std::floor(virtual_index));
    return interpolate(sorted[below], sorted[below + 1], virtual_index - static_cast<double>(below));
}

std::vector<double> nanquantile_linear(std::vector<double> values,
                                       const std::vector<double>& quantiles) {
    for (const double q : quantiles) {
        if (!(q >= 0.0 && q <= 1.0)) {
            throw std::invalid_argument("Quantiles must be in the range [0, 1]");
        }
    }
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return std::isnan(v); }),
                 values.end());
    std::sort(values.begin(), values.end());
    std::vector<double> result;
    result.reserve(quantiles.size());
    for (const double q : quantiles) {
        result.push_back(quantile_linear_sorted(values, q));
    }
    return result;
}

double percentile_linear(std::vector<double> values, double q) {
    if (values.empty() || has_nan(values)) {
        return kNaN;
    }
    std::sort(values.begin(), values.end());
    return quantile_linear_sorted(values, q / 100.0);
}

std::vector<double> linspace(double start, double stop, std::int64_t num) {
    if (num < 0) {
        throw std::invalid_argument("Number of samples, " + std::to_string(num) +
                                    ", must be non-negative.");
    }
    std::vector<double> samples(static_cast<std::size_t>(num));
    const double divisions = static_cast<double>(num - 1);
    const double span = stop - start;
    const double step = num > 1 ? span / divisions : 0.0;
    for (std::int64_t k = 0; k < num; ++k) {
        double offset = 0.0;
        if (num > 1 && step == 0.0) {
            // A step that underflows to zero: divide the index first.
            offset = static_cast<double>(k) / divisions * span;
        } else if (num > 1) {
            offset = static_cast<double>(k) * step;
        } else {
            offset = static_cast<double>(k) * span;
        }
        samples[static_cast<std::size_t>(k)] = offset + start;
    }
    if (num > 1) {
        samples.back() = stop;
    }
    return samples;
}

std::vector<std::int64_t> searchsorted_right(const std::vector<double>& sorted,
                                             const std::vector<double>& keys) {
    std::vector<std::int64_t> result(keys.size());
    const auto size = static_cast<std::int64_t>(sorted.size());
    std::int64_t window_begin = 0;
    std::int64_t window_end = size;
    for (std::size_t k = 0; k < keys.size(); ++k) {
        const double key = keys[k];
        const double previous = k == 0 ? key : keys[k - 1];
        if (!nan_last_less(key, previous)) {
            // Not smaller than the previous key: the answer is at or behind it.
            window_end = size;
        } else {
            // Smaller: at or before the previous answer (one slot of slack).
            window_begin = 0;
            window_end = window_end < size ? window_end + 1 : size;
        }
        while (window_begin < window_end) {
            const std::int64_t middle = window_begin + ((window_end - window_begin) >> 1);
            if (!nan_last_less(key, sorted[static_cast<std::size_t>(middle)])) {
                window_begin = middle + 1;
            } else {
                window_end = middle;
            }
        }
        result[k] = window_begin;
    }
    return result;
}

}  // namespace hicx::npy
