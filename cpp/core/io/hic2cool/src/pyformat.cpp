#include "pyformat.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <limits>

namespace hic2cool::detail {

std::string py_repr(const std::string& text) {
    const bool has_single = text.find('\'') != std::string::npos;
    const bool has_double = text.find('"') != std::string::npos;
    const char quote = has_single && !has_double ? '"' : '\'';
    std::string out(1, quote);
    for (unsigned char c : text) {
        if (c == '\\') {
            out += "\\\\";
        } else if (c == static_cast<unsigned char>(quote)) {
            out += '\\';
            out += static_cast<char>(c);
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\r') {
            out += "\\r";
        } else if (c == '\t') {
            out += "\\t";
        } else if (c < 0x20 || c == 0x7f) {
            char buffer[8];
            std::snprintf(buffer, sizeof(buffer), "\\x%02x", c);
            out += buffer;
        } else {
            out += static_cast<char>(c);
        }
    }
    return out + quote;
}

namespace {

template <class T, class F>
std::string join(const std::vector<T>& values, F format) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += format(values[i]);
    }
    return out + "]";
}

}  // namespace

std::string py_repr(const std::vector<std::string>& values) {
    return join(values, [](const std::string& v) { return py_repr(v); });
}

std::string py_repr(const std::vector<std::int32_t>& values) {
    return join(values, [](std::int32_t v) { return std::to_string(v); });
}

std::string py_repr(const std::vector<std::int64_t>& values) {
    return join(values, [](std::int64_t v) { return std::to_string(v); });
}

std::string py_repr(const std::vector<std::pair<std::int64_t, std::string>>& values) {
    std::string out = "{";
    for (std::size_t i = 0; i < values.size(); ++i) {
        out += (i > 0 ? ", " : "") + std::to_string(values[i].first) + ": " + py_repr(values[i].second);
    }
    return out + "}";
}

std::string py_repr(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    // float_repr_style 'short': the shortest digits that round trip, fixed
    // notation for decimal exponents in [-4, 16), scientific otherwise.
    char digits[64];
    auto result = std::to_chars(digits, digits + sizeof(digits), value, std::chars_format::scientific);
    std::string sci(digits, result.ptr);
    const auto e = sci.find('e');
    std::string mantissa = sci.substr(0, e);
    const int exponent = std::stoi(sci.substr(e + 1));
    const bool negative = mantissa[0] == '-';
    if (negative) {
        mantissa.erase(0, 1);
    }
    std::string significant;
    for (char c : mantissa) {
        if (c != '.') {
            significant += c;
        }
    }
    std::string out;
    if (exponent >= -4 && exponent < 16) {
        const int point = exponent + 1;
        if (point <= 0) {
            out = "0." + std::string(static_cast<std::size_t>(-point), '0') + significant;
        } else if (point >= static_cast<int>(significant.size())) {
            out = significant + std::string(static_cast<std::size_t>(point) - significant.size(), '0') + ".0";
        } else {
            out = significant.substr(0, static_cast<std::size_t>(point)) + "." +
                  significant.substr(static_cast<std::size_t>(point));
        }
    } else {
        out = significant.substr(0, 1);
        if (significant.size() > 1) {
            out += "." + significant.substr(1);
        }
        char buffer[16];
        std::snprintf(buffer, sizeof(buffer), "e%c%02d", exponent < 0 ? '-' : '+', std::abs(exponent));
        out += buffer;
    }
    return (negative ? "-" : "") + out;
}

std::string utcnow_isoformat() {
    const auto now = std::chrono::system_clock::now();
    const auto micros =
        std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count() % 1000000;
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
    gmtime_r(&seconds, &utc);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &utc);
    std::string text = buffer;
    if (micros != 0) {
        std::snprintf(buffer, sizeof(buffer), ".%06lld", static_cast<long long>(micros));
        text += buffer;
    }
    return text;
}

namespace {

std::string format_fixed(double value, int digits) {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer), "% .*f", digits, value);
    return buffer;
}

std::string format_sci(double value, int digits) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "% .*e", digits, value);
    return buffer;
}

std::string pad_left(const std::string& text, std::size_t width) {
    return text.size() >= width ? text : std::string(width - text.size(), ' ') + text;
}

}  // namespace

std::string pandas_bins_head(const std::string& norm, const std::vector<BinsPreviewRow>& rows) {
    // pandas FloatArrayFormatter with display.precision 6: "{: .6f}" with the
    // trailing zeros all numbers share trimmed, switching to "{: .6e}" when a
    // value would print as zero or a large value makes a string too long.
    constexpr int kDigits = 6;
    std::vector<std::string> values;
    bool has_large = false;
    bool has_small = false;
    for (const auto& row : rows) {
        if (std::isnan(row.value)) {
            values.emplace_back("NaN");
            continue;
        }
        const double magnitude = std::abs(row.value);
        has_large = has_large || magnitude > 1e6;
        has_small = has_small || (magnitude < 1e-6 && magnitude > 0);
        values.push_back(format_fixed(row.value, kDigits));
    }
    // _trim_zeros_float
    while (true) {
        bool trim = false;
        bool any = false;
        for (const auto& v : values) {
            if (v == "NaN") {
                continue;
            }
            any = true;
            const auto dot = v.find('.');
            if (dot == std::string::npos || v.back() != '0' || dot + 2 > v.size() - 1) {
                trim = false;
                break;
            }
            trim = true;
        }
        if (!any || !trim) {
            break;
        }
        for (auto& v : values) {
            if (v != "NaN") {
                v.pop_back();
            }
        }
    }
    std::size_t longest = 0;
    for (const auto& v : values) {
        longest = std::max(longest, v.size());
    }
    if (has_small || (longest > kDigits + 6 && has_large)) {
        for (std::size_t i = 0; i < rows.size(); ++i) {
            if (!std::isnan(rows[i].value)) {
                values[i] = format_sci(rows[i].value, kDigits);
            }
        }
    }

    // DataFrame.to_string: every cell carries a leading space (floats from
    // their "{: .6f}" format, NaN excepted), the names of numeric columns do
    // too, each column is right aligned to its widest entry, and columns
    // (the left aligned index first) are joined by one space.
    std::size_t index_width = 1;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        index_width = std::max(index_width, std::to_string(i).size());
    }
    std::size_t chrom_width = 5;              // "chrom"
    std::size_t start_width = 6;              // " start"
    std::size_t end_width = 4;                // " end"
    std::size_t norm_width = norm.size() + 1;  // " <norm>"
    for (std::size_t i = 0; i < rows.size(); ++i) {
        chrom_width = std::max(chrom_width, rows[i].chrom.size() + 1);
        start_width = std::max(start_width, std::to_string(rows[i].start).size() + 1);
        end_width = std::max(end_width, std::to_string(rows[i].end).size() + 1);
        norm_width = std::max(norm_width, values[i].size());
    }
    std::string out = std::string(index_width, ' ') + " " + pad_left("chrom", chrom_width) + " " +
                      pad_left("start", start_width) + " " + pad_left("end", end_width) + " " +
                      pad_left(norm, norm_width);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        std::string index = std::to_string(i);
        index += std::string(index_width - index.size(), ' ');
        out += "\n" + index + " " + pad_left(rows[i].chrom, chrom_width) + " " +
               pad_left(std::to_string(rows[i].start), start_width) + " " +
               pad_left(std::to_string(rows[i].end), end_width) + " " + pad_left(values[i], norm_width);
    }
    return out;
}

}  // namespace hic2cool::detail
