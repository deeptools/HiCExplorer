#include "numpy_compat.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string_view>
#include <system_error>

#include "coolercpp/errors.hpp"

namespace coolercpp::npy {

namespace {

constexpr std::size_t kPairwiseBlockSize = 128;
constexpr std::size_t kReduceBufferSize = 8192;

template <typename T>
T pairwise_block(const T* a, std::size_t n) {
    if (n < 8) {
        T res = 0;
        for (std::size_t i = 0; i < n; ++i) {
            res += a[i];
        }
        return res;
    }
    if (n <= kPairwiseBlockSize) {
        // Eight accumulators, exactly as numpy does it; the unrolling is part
        // of the result.
        std::array<T, 8> r{a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]};
        std::size_t i = 8;
        for (; i < n - (n % 8); i += 8) {
            r[0] += a[i + 0];
            r[1] += a[i + 1];
            r[2] += a[i + 2];
            r[3] += a[i + 3];
            r[4] += a[i + 4];
            r[5] += a[i + 5];
            r[6] += a[i + 6];
            r[7] += a[i + 7];
        }
        T res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) {
            res += a[i];
        }
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return pairwise_block(a, n2) + pairwise_block(a + n2, n - n2);
}

template <typename T>
T buffered_pairwise_sum(const T* a, std::size_t n) {
    T accumulated = 0;
    std::size_t offset = 0;
    while (offset < n) {
        const std::size_t block = std::min(kReduceBufferSize, n - offset);
        accumulated += pairwise_block(a + offset, block);
        offset += block;
    }
    return accumulated;
}

}  // namespace

double pairwise_sum(const double* a, std::size_t n) { return buffered_pairwise_sum(a, n); }

float pairwise_sum(const float* a, std::size_t n) { return buffered_pairwise_sum(a, n); }

namespace {

template <typename T>
std::string repr_impl(T value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0 ? "-inf" : "inf";
    }
    const bool negative = std::signbit(value);
    const T magnitude = negative ? -value : value;

    std::array<char, 64> buffer{};
    auto [end, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), magnitude,
                                   std::chars_format::scientific);
    if (ec != std::errc()) {
        return std::to_string(value);
    }
    const std::string scientific(buffer.data(), end);
    const std::size_t e_pos = scientific.find('e');
    std::string digits = scientific.substr(0, e_pos);
    const int exponent = std::atoi(scientific.c_str() + e_pos + 1);
    const std::size_t dot = digits.find('.');
    if (dot != std::string::npos) {
        digits.erase(dot, 1);
    }
    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
    }
    const int decpt = (digits == "0") ? 1 : exponent + 1;

    std::string out;
    if (negative) {
        out.push_back('-');
    }
    if (decpt <= -4 || decpt > 16) {
        out.push_back(digits[0]);
        if (digits.size() > 1) {
            out.push_back('.');
            out.append(digits, 1, std::string::npos);
        }
        out.push_back('e');
        const int exp10 = (digits == "0") ? 0 : exponent;
        out.push_back(exp10 < 0 ? '-' : '+');
        std::string exp_digits = std::to_string(exp10 < 0 ? -exp10 : exp10);
        if (exp_digits.size() < 2) {
            exp_digits.insert(exp_digits.begin(), '0');
        }
        out.append(exp_digits);
        return out;
    }
    if (decpt <= 0) {
        out.append("0.");
        out.append(static_cast<std::size_t>(-decpt), '0');
        out.append(digits);
    } else if (static_cast<std::size_t>(decpt) >= digits.size()) {
        out.append(digits);
        out.append(static_cast<std::size_t>(decpt) - digits.size(), '0');
        out.append(".0");
    } else {
        out.append(digits, 0, static_cast<std::size_t>(decpt));
        out.push_back('.');
        out.append(digits, static_cast<std::size_t>(decpt), std::string::npos);
    }
    return out;
}

}  // namespace

std::string float_repr(double value) { return repr_impl(value); }

std::string float32_repr(float value) { return repr_impl(value); }

std::vector<std::int64_t> linspace_int(std::int64_t lo, std::int64_t hi, std::int64_t num) {
    if (num < 0) {
        throw ValueError("Number of samples, " + std::to_string(num) +
                         ", must be non-negative.");
    }
    // numpy/core/function_base.py linspace, endpoint=True, scalar bounds.
    const double start = static_cast<double>(lo) * 1.0;
    const double stop = static_cast<double>(hi) * 1.0;
    const std::int64_t div = num - 1;
    const double delta = stop - start;
    std::vector<double> y(static_cast<std::size_t>(num));
    for (std::int64_t k = 0; k < num; ++k) {
        y[static_cast<std::size_t>(k)] = static_cast<double>(k);
    }
    if (div > 0) {
        const double step = delta / static_cast<double>(div);
        if (step == 0.0) {
            for (double& v : y) {
                v /= static_cast<double>(div);
                v = v * delta;
            }
        } else {
            for (double& v : y) {
                v *= step;
            }
        }
    } else {
        for (double& v : y) {
            v = v * delta;
        }
    }
    for (double& v : y) {
        v += start;
    }
    if (num > 1) {
        y.back() = stop;
    }
    std::vector<std::int64_t> out(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) {
        out[i] = static_cast<std::int64_t>(std::floor(y[i]));
    }
    return out;
}

SliceRange slice_indices(std::optional<std::int64_t> start, std::optional<std::int64_t> stop,
                         std::int64_t length) noexcept {
    const auto clamp = [length](std::optional<std::int64_t> value,
                                std::int64_t fallback) -> std::int64_t {
        if (!value.has_value()) {
            return fallback;
        }
        std::int64_t v = *value;
        if (v < 0) {
            v += length;
            if (v < 0) {
                v = 0;
            }
        } else if (v > length) {
            v = length;
        }
        return v;
    };
    SliceRange out;
    out.start = clamp(start, 0);
    out.stop = clamp(stop, length);
    if (out.stop < out.start) {
        out.stop = out.start;
    }
    return out;
}

std::int64_t wrap_index(std::int64_t i, std::int64_t length) {
    const std::int64_t original = i;
    if (i < 0) {
        i += length;
    }
    if (i < 0 || i >= length) {
        throw IndexError("index " + std::to_string(original) +
                         " is out of bounds for axis 0 with size " + std::to_string(length));
    }
    return i;
}

std::string iso_now() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    const auto micros =
        std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()) %
        std::chrono::seconds(1);
    std::tm parts{};
    localtime_r(&seconds, &parts);
    char buffer[64];
    if (micros.count() == 0) {
        std::snprintf(buffer, sizeof(buffer), "%04d-%02d-%02dT%02d:%02d:%02d",
                      parts.tm_year + 1900, parts.tm_mon + 1, parts.tm_mday, parts.tm_hour,
                      parts.tm_min, parts.tm_sec);
    } else {
        std::snprintf(buffer, sizeof(buffer), "%04d-%02d-%02dT%02d:%02d:%02d.%06lld",
                      parts.tm_year + 1900, parts.tm_mon + 1, parts.tm_mday, parts.tm_hour,
                      parts.tm_min, parts.tm_sec, static_cast<long long>(micros.count()));
    }
    return buffer;
}

}  // namespace coolercpp::npy
