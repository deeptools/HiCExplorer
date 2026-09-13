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

constexpr std::size_t kPairwiseBlockSize = 128;

// numpy's ufunc reduction buffer, np.getbufsize(). np.add.reduce does not run
// the pairwise summation over the whole array: the inner loop is called once
// per buffer, and the block results are accumulated sequentially. Arrays
// longer than this therefore have a different rounding than a single pairwise
// pass would give, which is visible in the last digits of the sums hicInfo
// prints.
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
        // Eight accumulators, exactly as numpy does it. The unrolling is part
        // of the result, not an optimisation.
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
    // Divide in two, but keep both halves a multiple of the unroll factor.
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

double pairwise_sum(const double* a, std::size_t n) {
    return buffered_pairwise_sum(a, n);
}

float pairwise_sum(const float* a, std::size_t n) {
    return buffered_pairwise_sum(a, n);
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
        accumulated_ += pairwise_block(buffer_.data(), buffer_.size());
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
    return accumulated_ + pairwise_block(buffer_.data(), buffer_.size());
}

template class PairwiseSumStream<double>;
template class PairwiseSumStream<float>;

std::string float_repr(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value < 0 ? "-inf" : "inf";
    }

    const bool negative = std::signbit(value);
    const double magnitude = negative ? -value : value;

    // Shortest round tripping digits plus the decimal exponent.
    std::array<char, 64> buffer{};
    auto [end, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                   magnitude, std::chars_format::scientific);
    if (ec != std::errc()) {
        return std::to_string(value);
    }
    const std::string scientific(buffer.data(), end);

    // scientific looks like "d[.ddd]e[+-]dd"
    const std::size_t e_pos = scientific.find('e');
    std::string digits = scientific.substr(0, e_pos);
    const int exponent = std::atoi(scientific.c_str() + e_pos + 1);
    const std::size_t dot = digits.find('.');
    if (dot != std::string::npos) {
        digits.erase(dot, 1);
    }
    // Strip trailing zeros; to_chars shortest never emits them, but 0.0 gives
    // a single "0" which is handled below.
    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
    }

    // decpt is CPython's decimal point position: value == 0.<digits> * 10^decpt
    const int decpt = (digits == "0") ? 1 : exponent + 1;

    std::string out;
    if (negative) {
        out.push_back('-');
    }

    // CPython's 'r' format code: exponential when decpt <= -4 or decpt > 16.
    if (decpt <= -4 || decpt > 16) {
        out.push_back(digits[0]);
        // CPython does not pad the exponential form to a fractional digit:
        // repr(1e16) is '1e+16', not '1.0e+16'.
        if (digits.size() > 1) {
            out.push_back('.');
            out.append(digits, 1, std::string::npos);
        }
        out.push_back('e');
        const int exp10 = (digits == "0") ? 0 : exponent;
        out.push_back(exp10 < 0 ? '-' : '+');
        const int abs_exp = exp10 < 0 ? -exp10 : exp10;
        std::string exp_digits = std::to_string(abs_exp);
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
        out.append(".0");  // Py_DTSF_ADD_DOT_0
    } else {
        out.append(digits, 0, static_cast<std::size_t>(decpt));
        out.push_back('.');
        out.append(digits, static_cast<std::size_t>(decpt), std::string::npos);
    }
    return out;
}

std::string int_with_thousands_separator(std::int64_t value) {
    const bool negative = value < 0;
    std::string digits = std::to_string(negative ? -value : value);
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
    // numpy's default line width for str(array).
    constexpr std::size_t kLineWidth = 75;
    if (items.empty()) {
        return "[]";
    }
    std::string out = "[";
    std::size_t line_length = 1;
    for (std::size_t i = 0; i < items.size(); ++i) {
        const std::string word = "'" + items[i] + "'";
        const bool last = (i + 1 == items.size());
        const std::size_t needed = word.size() + (last ? 1 : 1);  // ']' or ' '
        if (i > 0 && line_length + needed > kLineWidth) {
            out.append("\n ");
            line_length = 1;
        }
        out.append(word);
        line_length += word.size();
        if (!last) {
            out.push_back(' ');
            line_length += 1;
        }
    }
    out.push_back(']');
    return out;
}

}  // namespace hicx::npy
