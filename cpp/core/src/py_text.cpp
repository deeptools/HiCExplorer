#include "hicx/py_text.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <cstdint>

namespace hicx::pytext {

namespace {

inline unsigned char byte_at(std::string_view text, std::size_t pos) noexcept {
    return static_cast<unsigned char>(text[pos]);
}

inline bool is_continuation(unsigned char byte) noexcept {
    return (byte & 0xC0U) == 0x80U;
}

struct CodePointRange {
    std::uint32_t first;
    std::uint32_t last;
};

// The code points str.isprintable() rejects (see gen_py_text_tables.py).
constexpr CodePointRange kNonPrintable[] = {
#include "py_text_nonprintable.inc"
};

// The zero of each decimal digit block (see gen_py_text_tables.py).
constexpr std::uint32_t kDigitZeros[] = {
#include "py_text_digit_zeros.inc"
};

// The value 0..9 of a Unicode decimal digit, or -1.
int decimal_value(std::uint32_t cp) noexcept {
    const auto* end = std::end(kDigitZeros);
    const auto* it = std::upper_bound(std::begin(kDigitZeros), end, cp);
    if (it == std::begin(kDigitZeros)) {
        return -1;
    }
    return cp - *(it - 1) < 10 ? static_cast<int>(cp - *(it - 1)) : -1;
}

bool is_printable(std::uint32_t cp) noexcept {
    const auto* end = std::end(kNonPrintable);
    const auto* it = std::lower_bound(std::begin(kNonPrintable), end, cp,
                                      [](const CodePointRange& range, std::uint32_t value) {
                                          return range.last < value;
                                      });
    return it == end || cp < it->first;
}

// Decodes the strictly valid UTF-8 sequence at text[pos]; returns its length,
// or 0 when the bytes there are not a valid sequence.
std::size_t decode_strict(std::string_view text, std::size_t pos, std::uint32_t& cp) noexcept {
    const unsigned char b0 = byte_at(text, pos);
    if (b0 < 0x80U) {
        cp = b0;
        return 1;
    }
    const std::string_view rest = text.substr(pos);
    const std::size_t length = b0 >= 0xF0U ? 4 : (b0 >= 0xE0U ? 3 : 2);
    if (b0 < 0xC2U || b0 > 0xF4U || rest.size() < length || first_invalid_utf8(rest.substr(0, length))) {
        return 0;
    }
    cp = length == 2 ? (b0 & 0x1FU) : (length == 3 ? (b0 & 0x0FU) : (b0 & 0x07U));
    for (std::size_t k = 1; k < length; ++k) {
        cp = (cp << 6) | (byte_at(text, pos + k) & 0x3FU);
    }
    return length;
}

void append_escape(std::string& out, std::uint32_t cp) {
    constexpr char kHex[] = "0123456789abcdef";
    const int digits = cp < 0x100 ? 2 : (cp < 0x10000 ? 4 : 8);
    out += cp < 0x100 ? "\\x" : (cp < 0x10000 ? "\\u" : "\\U");
    for (int shift = 4 * (digits - 1); shift >= 0; shift -= 4) {
        out += kHex[(cp >> shift) & 0xFU];
    }
}

}  // namespace

std::string numeric_text(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::uint32_t cp = 0;
        const std::size_t length = decode_strict(text, pos, cp);
        const int digit = length == 0 ? -1 : decimal_value(cp);
        if (digit >= 0) {
            out += static_cast<char>('0' + digit);
        } else if (length != 0 && cp >= 0x80 && whitespace_length(text, pos) == length) {
            out += ' ';
        } else {
            out.append(text.substr(pos, length == 0 ? 1 : length));
        }
        pos += length == 0 ? 1 : length;
    }
    return out;
}

std::string repr(std::string_view text) {
    const bool has_single = text.find('\'') != std::string_view::npos;
    const bool has_double = text.find('"') != std::string_view::npos;
    const char quote = has_single && !has_double ? '"' : '\'';
    std::string out(1, quote);
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::uint32_t cp = 0;
        const std::size_t length = decode_strict(text, pos, cp);
        if (length == 0) {
            append_escape(out, 0xDC00U + byte_at(text, pos));
            ++pos;
            continue;
        }
        if (cp == static_cast<unsigned char>(quote) || cp == '\\') {
            out += '\\';
            out += static_cast<char>(cp);
        } else if (cp == '\t') {
            out += "\\t";
        } else if (cp == '\n') {
            out += "\\n";
        } else if (cp == '\r') {
            out += "\\r";
        } else if (cp != ' ' && !is_printable(cp)) {
            append_escape(out, cp);
        } else {
            out.append(text.substr(pos, length));
        }
        pos += length;
    }
    out += quote;
    return out;
}

std::size_t whitespace_length(std::string_view text, std::size_t pos) noexcept {
    const std::size_t left = text.size() - pos;
    const unsigned char b0 = byte_at(text, pos);
    if (b0 < 0x80U) {
        return ((b0 >= 0x09U && b0 <= 0x0DU) || (b0 >= 0x1CU && b0 <= 0x1FU) || b0 == 0x20U)
                   ? 1
                   : 0;
    }
    if (b0 == 0xC2U && left >= 2) {
        const unsigned char b1 = byte_at(text, pos + 1);
        return (b1 == 0x85U || b1 == 0xA0U) ? 2 : 0;
    }
    if (left < 3) {
        return 0;
    }
    const unsigned char b1 = byte_at(text, pos + 1);
    const unsigned char b2 = byte_at(text, pos + 2);
    if (b0 == 0xE1U) {
        return (b1 == 0x9AU && b2 == 0x80U) ? 3 : 0;  // U+1680
    }
    if (b0 == 0xE2U) {
        if (b1 == 0x80U) {
            // U+2000-U+200A, U+2028, U+2029, U+202F
            return ((b2 >= 0x80U && b2 <= 0x8AU) || b2 == 0xA8U || b2 == 0xA9U || b2 == 0xAFU)
                       ? 3
                       : 0;
        }
        return (b1 == 0x81U && b2 == 0x9FU) ? 3 : 0;  // U+205F
    }
    if (b0 == 0xE3U) {
        return (b1 == 0x80U && b2 == 0x80U) ? 3 : 0;  // U+3000
    }
    return 0;
}

std::size_t previous_code_point_length(std::string_view text, std::size_t end) noexcept {
    std::size_t length = 1;
    while (length < 4 && length < end && is_continuation(byte_at(text, end - length))) {
        ++length;
    }
    return length;
}

std::size_t code_point_count(std::string_view text) noexcept {
    std::size_t count = 0;
    for (const char c : text) {
        if (!is_continuation(static_cast<unsigned char>(c))) {
            ++count;
        }
    }
    return count;
}

std::string_view strip(std::string_view text) noexcept {
    std::size_t begin = 0;
    while (begin < text.size()) {
        const std::size_t length = whitespace_length(text, begin);
        if (length == 0) {
            break;
        }
        begin += length;
    }
    std::size_t end = text.size();
    while (end > begin) {
        const std::size_t length = previous_code_point_length(text, end);
        if (whitespace_length(text, end - length) != length) {
            break;
        }
        end -= length;
    }
    return text.substr(begin, end - begin);
}

std::vector<std::string_view> split_whitespace(std::string_view text) {
    std::vector<std::string_view> fields;
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::size_t length = whitespace_length(text, pos);
        while (length > 0) {
            pos += length;
            length = pos < text.size() ? whitespace_length(text, pos) : 0;
        }
        if (pos >= text.size()) {
            break;
        }
        const std::size_t begin = pos;
        while (pos < text.size() && whitespace_length(text, pos) == 0) {
            pos += 1;
            while (pos < text.size() && is_continuation(byte_at(text, pos))) {
                ++pos;
            }
        }
        fields.push_back(text.substr(begin, pos - begin));
    }
    return fields;
}

std::vector<std::string_view> split(std::string_view text, std::string_view sep) {
    std::vector<std::string_view> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t found = text.find(sep, begin);
        if (found == std::string_view::npos) {
            fields.push_back(text.substr(begin));
            return fields;
        }
        fields.push_back(text.substr(begin, found - begin));
        begin = found + sep.size();
    }
}

std::string remove_whitespace(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    std::size_t pos = 0;
    while (pos < text.size()) {
        const std::size_t length = whitespace_length(text, pos);
        if (length > 0) {
            pos += length;
            continue;
        }
        out.push_back(text[pos++]);
        while (pos < text.size() && is_continuation(byte_at(text, pos))) {
            out.push_back(text[pos++]);
        }
    }
    return out;
}

std::optional<std::size_t> first_invalid_utf8(std::string_view bytes) noexcept {
    const std::size_t n = bytes.size();
    std::size_t i = 0;
    while (i < n) {
        const unsigned char b0 = byte_at(bytes, i);
        if (b0 < 0x80U) {
            ++i;
            continue;
        }
        std::size_t need = 0;
        unsigned char low = 0x80U;
        unsigned char high = 0xBFU;
        if (b0 >= 0xC2U && b0 <= 0xDFU) {
            need = 1;
        } else if (b0 == 0xE0U) {
            need = 2;
            low = 0xA0U;
        } else if ((b0 >= 0xE1U && b0 <= 0xECU) || b0 == 0xEEU || b0 == 0xEFU) {
            need = 2;
        } else if (b0 == 0xEDU) {
            need = 2;
            high = 0x9FU;  // no surrogates
        } else if (b0 == 0xF0U) {
            need = 3;
            low = 0x90U;
        } else if (b0 >= 0xF1U && b0 <= 0xF3U) {
            need = 3;
        } else if (b0 == 0xF4U) {
            need = 3;
            high = 0x8FU;
        } else {
            return i;
        }
        if (i + need >= n) {
            return i;  // cut off by the end of the text
        }
        const unsigned char b1 = byte_at(bytes, i + 1);
        if (b1 < low || b1 > high) {
            return i;
        }
        for (std::size_t k = 2; k <= need; ++k) {
            if (!is_continuation(byte_at(bytes, i + k))) {
                return i;
            }
        }
        i += need + 1;
    }
    return std::nullopt;
}

std::vector<std::string> universal_newline_lines(std::string_view text) {
    std::vector<std::string> lines;
    std::string current;
    const std::size_t n = text.size();
    for (std::size_t i = 0; i < n; ++i) {
        const char c = text[i];
        if (c == '\r' || c == '\n') {
            if (c == '\r' && i + 1 < n && text[i + 1] == '\n') {
                ++i;
            }
            current.push_back('\n');
            lines.push_back(std::move(current));
            current.clear();
            continue;
        }
        current.push_back(c);
    }
    if (!current.empty()) {
        lines.push_back(std::move(current));
    }
    return lines;
}

}  // namespace hicx::pytext
