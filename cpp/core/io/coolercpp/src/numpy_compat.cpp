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

// Pairwise summation after N. J. Higham, "The accuracy of floating point
// summation", SIAM J. Sci. Comput. 14(4), 1993: a long sum is split in two
// halves that are summed separately and added, which keeps the rounding error
// growth logarithmic in the length. numpy's documentation (numpy.sum, Notes)
// states how its float summation is arranged, and the same arrangement is
// reproduced here because the order of the additions decides the last bits:
//   * below eight elements the values are added left to right;
//   * up to 128 elements, eight partial sums run side by side (element i goes
//     to partial sum i mod 8) and are combined as ((0+1)+(2+3))+((4+5)+(6+7)),
//     after which the leftover elements are added one by one;
//   * longer ranges are cut in two, the first half being a multiple of eight
//     long;
//   * a reduction over a whole array is not one pairwise pass but one pass per
//     buffer of 8192 elements (the documented ufunc buffer size), with the
//     buffer results added to a running total in order.
// Identity with np.add.reduce was verified bit for bit on random float32 and
// float64 arrays, see cpp/tests/test_numpy_oracle.cpp and
// cpp/scripts/gen_numpy_oracle_fixture.py.
constexpr std::size_t kUnroll = 8;
constexpr std::size_t kUnrolledLimit = 128;
constexpr std::size_t kReduceBufferSize = 8192;

template <typename T>
T pairwise_range(const T* a, std::size_t n) {
    if (n < kUnroll) {
        T total = 0;
        for (std::size_t i = 0; i < n; ++i) {
            total += a[i];
        }
        return total;
    }
    if (n <= kUnrolledLimit) {
        std::array<T, kUnroll> lane;
        for (std::size_t k = 0; k < kUnroll; ++k) {
            lane[k] = a[k];
        }
        const std::size_t whole = n - n % kUnroll;
        std::size_t i = kUnroll;
        for (; i < whole; i += kUnroll) {
            for (std::size_t k = 0; k < kUnroll; ++k) {
                lane[k] += a[i + k];
            }
        }
        T total = ((lane[0] + lane[1]) + (lane[2] + lane[3])) +
                  ((lane[4] + lane[5]) + (lane[6] + lane[7]));
        for (; i < n; ++i) {
            total += a[i];
        }
        return total;
    }
    std::size_t first = n / 2;
    first -= first % kUnroll;
    return pairwise_range(a, first) + pairwise_range(a + first, n - first);
}

// One pairwise pass per buffer, the buffer sums accumulated in order.
template <typename T>
T buffered_pairwise_sum(const T* a, std::size_t n) {
    T total = 0;
    for (std::size_t offset = 0; offset < n; offset += kReduceBufferSize) {
        total += pairwise_range(a + offset, std::min(kReduceBufferSize, n - offset));
    }
    return total;
}

}  // namespace

double pairwise_sum(const double* a, std::size_t n) { return buffered_pairwise_sum(a, n); }

float pairwise_sum(const float* a, std::size_t n) { return buffered_pairwise_sum(a, n); }

namespace {

// The shortest decimal string that converts back to the same value (the
// std::to_chars guarantee) laid out by the rules Python documents for
// repr(float): positional notation for 1e-4 <= |x| < 1e16 (compared as the stored value) with at least one
// digit after the point, otherwise d[.ddd]e[+-]XX with an exponent of at least
// two digits and no ".0" for a single digit. numpy prints float32 scalars with
// the same layout from the float32 shortest digits. Checked against Python's
// repr and numpy's str on random doubles and floats (all exponents including
// subnormals, zeros, infinities, NaN), see cpp/tests/test_numpy_oracle.cpp.
template <typename T>
std::string repr_impl(T value) {
    if (std::isnan(value)) {
        return "nan";
    }
    const bool negative = std::signbit(value);
    const std::string sign = negative ? "-" : "";
    if (std::isinf(value)) {
        return sign + "inf";
    }
    if (value == 0) {
        return sign + "0.0";
    }

    std::array<char, 64> buffer{};
    const auto converted = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                         negative ? -value : value, std::chars_format::scientific);
    if (converted.ec != std::errc()) {
        return std::to_string(static_cast<double>(value));
    }
    // "d[.ddd]e[+-]XX": collect the significant digits and the decimal exponent.
    const std::string text(buffer.data(), converted.ptr);
    const std::size_t e_pos = text.find('e');
    std::string digits;
    for (std::size_t i = 0; i < e_pos; ++i) {
        if (text[i] != '.') {
            digits.push_back(text[i]);
        }
    }
    const int exponent = std::atoi(text.c_str() + e_pos + 1);  // value = d.ddd * 10^exponent

    // The switch to exponent form is decided on the value itself, not on the
    // digits: float32(1e-4) is slightly below 1e-4 and prints as 1e-04.
    const T magnitude = negative ? -value : value;
    if (static_cast<double>(magnitude) < 1e-4 || static_cast<double>(magnitude) >= 1e16) {
        std::string out = sign;
        out.push_back(digits[0]);
        if (digits.size() > 1) {
            out.push_back('.');
            out.append(digits, 1, std::string::npos);
        }
        out.push_back('e');
        out.push_back(exponent < 0 ? '-' : '+');
        const std::string exp_text = std::to_string(exponent < 0 ? -exponent : exponent);
        if (exp_text.size() < 2) {
            out.push_back('0');
        }
        out.append(exp_text);
        return out;
    }

    const int point = exponent + 1;  // number of digits before the decimal point
    std::string out = sign;
    if (point <= 0) {
        out.append("0.");
        out.append(static_cast<std::size_t>(-point), '0');
        out.append(digits);
    } else if (static_cast<std::size_t>(point) >= digits.size()) {
        out.append(digits);
        out.append(static_cast<std::size_t>(point) - digits.size(), '0');
        out.append(".0");
    } else {
        out.append(digits, 0, static_cast<std::size_t>(point));
        out.push_back('.');
        out.append(digits, static_cast<std::size_t>(point), std::string::npos);
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
    // numpy.linspace as its documentation defines it, endpoint included: the
    // samples are start + k * step with step = (stop - start) / (num - 1), the
    // last sample being stop itself. (A step that underflows to zero is
    // handled by dividing k first, as the documentation notes.) The int
    // conversion floors. Verified against np.linspace(dtype=int).
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
