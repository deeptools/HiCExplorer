// Port of cooler/util.py parse_cooler_uri, parse_humanized,
// parse_region_string and parse_region (cooler 0.10.2, BSD-3-Clause).

#include "coolercpp/region.hpp"

#include <cerrno>
#include <cstdlib>

#include "coolercpp/errors.hpp"

namespace coolercpp {

ChromSizes::ChromSizes(std::vector<std::string> names, std::vector<std::int64_t> lengths,
                       DType dtype)
    : names_(std::move(names)), lengths_(std::move(lengths)), dtype_(dtype) {
    if (names_.size() != lengths_.size()) {
        throw ValueError("chromosome names and lengths differ in length");
    }
}

std::optional<std::int64_t> ChromSizes::find(std::string_view name) const {
    for (std::size_t i = 0; i < names_.size(); ++i) {
        if (names_[i] == name) {
            return lengths_[i];
        }
    }
    return std::nullopt;
}

std::int64_t ChromSizes::operator[](std::string_view name) const {
    const std::optional<std::int64_t> found = find(name);
    if (!found.has_value()) {
        throw KeyError(std::string(name));
    }
    return *found;
}

std::pair<std::string, std::string> parse_cooler_uri(std::string_view uri) {
    std::vector<std::string> parts;
    std::size_t pos = 0;
    while (true) {
        const std::size_t marker = uri.find("::", pos);
        if (marker == std::string_view::npos) {
            parts.emplace_back(uri.substr(pos));
            break;
        }
        parts.emplace_back(uri.substr(pos, marker - pos));
        pos = marker + 2;
    }
    if (parts.size() == 1) {
        return {parts[0], "/"};
    }
    if (parts.size() == 2) {
        std::string group = parts[1];
        if (group.empty() || group.front() != '/') {
            group.insert(group.begin(), '/');
        }
        return {parts[0], group};
    }
    throw ValueError("Invalid Cooler URI string");
}

namespace {

// Python's str.isspace over ASCII, which is what re's \s matches there.
bool py_space(char c) {
    const auto u = static_cast<unsigned char>(c);
    return u == ' ' || (u >= 0x09 && u <= 0x0D) || (u >= 0x1C && u <= 0x1F);
}

bool is_digit(char c) { return c >= '0' && c <= '9'; }
bool is_digit_or_comma(char c) { return is_digit(c) || c == ','; }
bool is_alpha(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }

std::string py_strip(std::string_view s) {
    std::size_t begin = 0;
    std::size_t end = s.size();
    while (begin < end && py_space(s[begin])) {
        ++begin;
    }
    while (end > begin && py_space(s[end - 1])) {
        --end;
    }
    return std::string(s.substr(begin, end - begin));
}

enum class TokenType { Hyphen, Coord, Other };

struct Token {
    TokenType type;
    std::string text;
};

// re.finditer over
//   \s*(?P<HYPHEN>-)|\s*(?P<COORD>[0-9,]+(\.[0-9]*)?(?:[a-z]+)?)|\s*(?P<OTHER>.+)
// with re.IGNORECASE, including the backtracking of the leading \s*.
std::vector<Token> tokenize(std::string_view s) {
    std::vector<Token> tokens;
    std::size_t p = 0;
    while (p < s.size()) {
        std::size_t q = p;
        while (q < s.size() && py_space(s[q])) {
            ++q;
        }
        if (q < s.size() && s[q] == '-') {
            tokens.push_back({TokenType::Hyphen, "-"});
            p = q + 1;
            continue;
        }
        if (q < s.size() && is_digit_or_comma(s[q])) {
            std::size_t e = q;
            while (e < s.size() && is_digit_or_comma(s[e])) {
                ++e;
            }
            if (e < s.size() && s[e] == '.') {
                ++e;
                while (e < s.size() && is_digit(s[e])) {
                    ++e;
                }
            }
            while (e < s.size() && is_alpha(s[e])) {
                ++e;
            }
            tokens.push_back({TokenType::Coord, std::string(s.substr(q, e - q))});
            p = e;
            continue;
        }
        // OTHER: '.' matches anything but a newline. When the greedy
        // whitespace run stops at a newline or at the end, the engine gives
        // whitespace back until '.' can match one character.
        std::size_t r = q;
        if (r >= s.size() || s[r] == '\n') {
            bool found = false;
            while (r > p) {
                --r;
                if (s[r] != '\n') {
                    found = true;
                    break;
                }
            }
            if (!found) {
                ++p;
                continue;
            }
        }
        std::size_t e = r;
        while (e < s.size() && s[e] != '\n') {
            ++e;
        }
        tokens.push_back({TokenType::Other, std::string(s.substr(r, e - r))});
        p = e;
    }
    return tokens;
}

const char* token_name(TokenType type) {
    switch (type) {
        case TokenType::Hyphen: return "HYPHEN";
        case TokenType::Coord: return "COORD";
        case TokenType::Other: return "OTHER";
    }
    return "";
}

// int(text) for the digit strings parse_humanized passes it.
std::int64_t py_int(const std::string& text) {
    bool valid = !text.empty();
    for (const char c : text) {
        valid = valid && is_digit(c);
    }
    if (!valid) {
        throw ValueError("invalid literal for int() with base 10: '" + text + "'");
    }
    errno = 0;
    const long long value = std::strtoll(text.c_str(), nullptr, 10);
    if (errno == ERANGE) {
        // Python has arbitrary precision integers; coolercpp does not.
        throw ValueError("integer out of the 64 bit range: '" + text + "'");
    }
    return value;
}

// float(text) for strings of digits and dots.
double py_float(const std::string& text) {
    std::size_t dots = 0;
    std::size_t digits = 0;
    for (const char c : text) {
        if (c == '.') {
            ++dots;
        } else if (is_digit(c)) {
            ++digits;
        }
    }
    if (dots > 1 || digits == 0) {
        throw ValueError("could not convert string to float: '" + text + "'");
    }
    return std::strtod(text.c_str(), nullptr);
}

std::string upper(std::string text) {
    for (char& c : text) {
        if (c >= 'a' && c <= 'z') {
            c = static_cast<char>(c - 'a' + 'A');
        }
    }
    return text;
}

}  // namespace

std::int64_t parse_humanized(std::string_view text) {
    std::string s;
    s.reserve(text.size());
    for (const char c : text) {
        if (c != ',') {
            s.push_back(c);
        }
    }
    // re.split("([0-9,.]+)", s)
    std::vector<std::string> parts;
    std::size_t pos = 0;
    while (true) {
        std::size_t i = pos;
        while (i < s.size() && !(is_digit(s[i]) || s[i] == '.')) {
            ++i;
        }
        if (i >= s.size()) {
            parts.push_back(s.substr(pos));
            break;
        }
        parts.push_back(s.substr(pos, i - pos));
        std::size_t j = i;
        while (j < s.size() && (is_digit(s[j]) || s[j] == '.')) {
            ++j;
        }
        parts.push_back(s.substr(i, j - i));
        pos = j;
    }
    if (parts.size() > 3) {
        throw ValueError("too many values to unpack (expected 3)");
    }
    if (parts.size() < 3) {
        throw ValueError("not enough values to unpack (expected 3, got " +
                         std::to_string(parts.size()) + ")");
    }
    const std::string& value = parts[1];
    const std::string& unit_text = parts[2];
    if (unit_text.empty()) {
        return py_int(value);
    }
    double number = py_float(value);
    const std::string unit = upper(py_strip(unit_text));
    if (unit == "K" || unit == "KB") {
        number *= 1000;
    } else if (unit == "M" || unit == "MB") {
        number *= 1000000;
    } else if (unit == "G" || unit == "GB") {
        number *= 1000000000;
    } else {
        throw ValueError("Unknown unit '" + unit + "'");
    }
    return static_cast<std::int64_t>(number);
}

GenomicRange parse_region_string(std::string_view text) {
    std::vector<std::string_view> parts;
    std::size_t pos = 0;
    while (true) {
        const std::size_t colon = text.find(':', pos);
        if (colon == std::string_view::npos) {
            parts.push_back(text.substr(pos));
            break;
        }
        parts.push_back(text.substr(pos, colon - pos));
        pos = colon + 1;
    }
    GenomicRange out;
    out.chrom = py_strip(parts[0]);
    if (out.chrom.empty()) {
        throw ValueError("Chromosome name cannot be empty");
    }
    if (parts.size() < 2) {
        return out;
    }
    const std::vector<Token> tokens = tokenize(parts[1]);
    std::size_t next = 0;
    const auto expect = [&](TokenType expected) -> const Token& {
        if (next >= tokens.size()) {
            throw ValueError(std::string("Expected ") + token_name(expected) +
                             " token missing");
        }
        const Token& token = tokens[next++];
        if (token.type != expected) {
            throw ValueError("Unexpected token \"" + token.text + "\"");
        }
        return token;
    };
    const std::int64_t start = parse_humanized(expect(TokenType::Coord).text);
    expect(TokenType::Hyphen);
    out.start = start;
    if (next >= tokens.size()) {
        return out;
    }
    const std::int64_t end = parse_humanized(expect(TokenType::Coord).text);
    if (end < start) {
        throw ValueError("End coordinate less than start");
    }
    out.end = end;
    return out;
}

RegionTuple parse_region(const Region& region, const ChromSizes* chromsizes) {
    GenomicRange range;
    if (region.is_string()) {
        range = parse_region_string(region.text());
    } else {
        range = GenomicRange{region.chrom(), region.start(), region.end()};
    }
    std::optional<std::int64_t> clen;
    if (chromsizes != nullptr) {
        clen = chromsizes->find(range.chrom);
        if (!clen.has_value()) {
            throw ValueError("Unknown sequence label: " + range.chrom);
        }
    }
    RegionTuple out;
    out.chrom = range.chrom;
    out.start = range.start.value_or(0);
    if (range.end.has_value()) {
        out.end = *range.end;
    } else {
        if (!clen.has_value()) {
            throw ValueError("Cannot determine end coordinate.");
        }
        out.end = *clen;
    }
    if (out.end < out.start) {
        throw ValueError("End cannot be less than start");
    }
    if (out.start < 0 || (clen.has_value() && out.end > *clen)) {
        throw ValueError("Genomic region out of bounds: [" + std::to_string(out.start) + ", " +
                         std::to_string(out.end) + ")");
    }
    return out;
}

}  // namespace coolercpp
