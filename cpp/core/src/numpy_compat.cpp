#include "hicx/numpy_compat.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <system_error>

namespace hicx::npy {

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

double pairwise_sum(const double* a, std::size_t n) {
    return buffered_pairwise_sum(a, n);
}

float pairwise_sum(const float* a, std::size_t n) {
    return buffered_pairwise_sum(a, n);
}

double pairwise_sum_unbuffered(const double* a, std::size_t n) {
    return pairwise_range(a, n);
}

float pairwise_sum_unbuffered(const float* a, std::size_t n) {
    return pairwise_range(a, n);
}

template <typename T>
PairwiseSumStream<T>::PairwiseSumStream() {
    buffer_.reserve(kReduceBufferSize);
}

template <typename T>
void PairwiseSumStream<T>::add(T value) {
    buffer_.push_back(value);
    ++count_;
    if (buffer_.size() == kReduceBufferSize) {
        accumulated_ += pairwise_range(buffer_.data(), buffer_.size());
        buffer_.clear();
    }
}

template <typename T>
T PairwiseSumStream<T>::result() const {
    // buffered_pairwise_sum adds a block only when there is one, so an empty
    // tail must not be added either: 0.0 + -0.0 would turn a -0.0 into 0.0.
    if (buffer_.empty()) {
        return accumulated_;
    }
    return accumulated_ + pairwise_range(buffer_.data(), buffer_.size());
}

template class PairwiseSumStream<double>;
template class PairwiseSumStream<float>;

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
std::string shortest_repr(T value) {
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

std::string float_repr(double value) {
    return shortest_repr(value);
}

std::string float32_repr(float value) {
    return shortest_repr(value);
}

std::string int_with_thousands_separator(std::int64_t value) {
    const bool negative = value < 0;
    // Unsigned magnitude, so that the most negative int64 is not negated.
    const std::uint64_t magnitude = negative ? 0 - static_cast<std::uint64_t>(value)
                                             : static_cast<std::uint64_t>(value);
    std::string digits = std::to_string(magnitude);
    std::string out;
    const std::size_t n = digits.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0 && (n - i) % 3 == 0) {
            out.push_back(',');
        }
        out.push_back(digits[i]);
    }
    return negative ? "-" + out : out;
}

std::string array_str(const std::vector<std::string>& items) {
    // str() of a numpy string array: quoted items separated by one space and
    // wrapped at numpy's default print width of 75 characters, continuation
    // lines indented by one space (the position after the opening bracket).
    // A word stays on the line when the line, the separating space, the word
    // and one more character (the separator after it, or the closing bracket)
    // fit in the width. Checked against str(np.array(items)) on 20000 random
    // string lists.
    constexpr std::size_t kPrintWidth = 75;
    if (items.empty()) {
        return "[]";
    }
    std::string out = "[";
    std::size_t line_length = 1;
    for (std::size_t i = 0; i < items.size(); ++i) {
        const std::string word = "'" + items[i] + "'";
        if (i > 0) {
            if (line_length + 1 + word.size() + 1 > kPrintWidth) {
                out.append("\n ");
                line_length = 1;
            } else {
                out.push_back(' ');
                ++line_length;
            }
        }
        out.append(word);
        line_length += word.size();
    }
    out.push_back(']');
    return out;
}

}  // namespace hicx::npy
