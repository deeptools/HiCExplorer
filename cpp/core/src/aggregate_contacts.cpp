#include "hicx/aggregate_contacts.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"

namespace hicx::aggregate {

std::int64_t pick_bin(std::int64_t first, std::int64_t last, LargeRegions operation) {
    switch (operation) {
        case LargeRegions::First: return first;
        case LargeRegions::Last: return last;
        case LargeRegions::Center:
        default: {
            const double mean = (static_cast<double>(first) + static_cast<double>(last)) / 2.0;
            return static_cast<std::int64_t>(std::floor(mean));
        }
    }
}

Orientation orientation_of(const std::optional<std::string>& strand1, const std::optional<std::string>& strand2) {
    if (!strand1 && !strand2) {
        return Orientation::None;
    }
    const auto is = [](const std::optional<std::string>& s, const char* value) { return s && *s == value; };
    if (is(strand1, "+") && is(strand2, "+")) {
        return Orientation::None;
    }
    if (is(strand1, "+") && is(strand2, "-")) {
        return Orientation::FlipLR;
    }
    if (is(strand1, "-") && is(strand2, "+")) {
        return Orientation::FlipUD;
    }
    if (is(strand1, "-") && is(strand2, "-")) {
        return Orientation::Transpose;
    }
    return Orientation::Unchanged;
}

Oriented orient(const std::vector<double>& window, std::size_t w, Orientation orientation) {
    Oriented out;
    out.values.resize(w * w);
    for (std::size_t i = 0; i < w; ++i) {
        for (std::size_t j = 0; j < w; ++j) {
            std::size_t source = i * w + j;
            switch (orientation) {
                case Orientation::FlipLR: source = i * w + (w - 1 - j); break;
                case Orientation::FlipUD: source = (w - 1 - i) * w + j; break;
                case Orientation::Transpose: source = j * w + i; break;
                default: break;
            }
            out.values[i * w + j] = window[source];
        }
    }
    const bool logical = orientation == Orientation::FlipLR || orientation == Orientation::FlipUD;
    out.sum = logical ? npy::pairwise_sum(out.values.data(), out.values.size())
                      : npy::pairwise_sum(window.data(), window.size());
    return out;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    for (const double value : values) {
        if (std::isnan(value)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const std::size_t n = values.size();
    const std::size_t half = n / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half), values.end());
    const double upper = values[half];
    if (n % 2 == 1) {
        return upper;
    }
    const double lower = *std::max_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half));
    return (lower + upper) / 2.0;
}

std::vector<double> axis0(const std::vector<double>& values, std::size_t count, std::size_t cell, Average average) {
    std::vector<double> out(cell, 0.0);
    if (average == Average::Median) {
        std::vector<double> column(count);
        for (std::size_t c = 0; c < cell; ++c) {
            for (std::size_t k = 0; k < count; ++k) {
                column[k] = values[k * cell + c];
            }
            out[c] = median(column);
        }
        return out;
    }
    for (std::size_t k = 0; k < count; ++k) {
        for (std::size_t c = 0; c < cell; ++c) {
            out[c] = k == 0 ? values[c] : out[c] + values[k * cell + c];
        }
    }
    if (average == Average::Mean) {
        for (double& v : out) {
            v = count == 0 ? std::numeric_limits<double>::quiet_NaN() : v / static_cast<double>(count);
        }
    }
    return out;
}

std::optional<std::vector<bool>> outlier_rows(const std::vector<double>& values, std::size_t rows, std::size_t cols,
                                              double max_deviation) {
    const double data_median = median(values);
    std::vector<double> absolute(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        absolute[i] = std::fabs(values[i]);
    }
    const double mad = 1.4826 * median(absolute);
    if (mad == 0.0) {
        return std::nullopt;
    }
    std::vector<bool> outliers(rows, false);
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            const double deviation = std::fabs(values[r * cols + c] - data_median) / mad;
            if (deviation > max_deviation) {
                outliers[r] = true;
                break;
            }
        }
    }
    return outliers;
}

DiagnosticCluster diagnostic_cluster(const std::vector<double>& heatmap, std::size_t width,
                                     const std::vector<std::int64_t>& members, std::size_t half) {
    DiagnosticCluster out;
    std::vector<double> keys(members.size());
    for (std::size_t m = 0; m < members.size(); ++m) {
        keys[m] = heatmap[static_cast<std::size_t>(members[m]) * width + half];
    }
    const std::vector<std::int64_t> order = npy::argsort(keys);
    out.mean.assign(width, 0.0);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const auto row = static_cast<std::size_t>(members[static_cast<std::size_t>(*it)]);
        for (std::size_t c = 0; c < width; ++c) {
            const double value = heatmap[row * width + c];
            out.mean[c] = it == order.rbegin() ? value : out.mean[c] + value;
        }
        out.rows.push_back(row);
    }
    for (double& value : out.mean) {
        value = value / static_cast<double>(members.size());
    }
    return out;
}

}  // namespace hicx::aggregate
