#include "hicx/py_text.hpp"

#include <cstdint>

namespace hicx::pytext {

namespace {

inline unsigned char byte_at(std::string_view text, std::size_t pos) noexcept {
    return static_cast<unsigned char>(text[pos]);
}

inline bool is_continuation(unsigned char byte) noexcept {
    return (byte & 0xC0U) == 0x80U;
}

}  // namespace

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
