#include "coolercpp/json.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

#include "coolercpp/errors.hpp"
#include "numpy_compat.hpp"

namespace coolercpp::json {

Value::Value(bool value) : type_(Type::Bool), bool_(value), dtype_(DType::Bool) {}

Value::Value(double value) : type_(Type::Double), double_(value), dtype_(DType::Float64) {}

Value::Value(const char* value) : type_(Type::String), string_(value) {}

Value::Value(std::string value) : type_(Type::String), string_(std::move(value)) {}

Value::Value(std::string_view value) : type_(Type::String), string_(value) {}

Value Value::integer(std::int64_t value, DType dtype) {
    Value out(value);
    out.dtype_ = dtype;
    return out;
}

Value Value::unsigned_integer(std::uint64_t value, DType dtype) {
    Value out(value);
    out.dtype_ = dtype;
    return out;
}

Value Value::number(double value, DType dtype) {
    Value out(value);
    out.dtype_ = dtype;
    return out;
}

Value Value::bytes(std::string value) {
    Value out(std::move(value));
    out.type_ = Type::Bytes;
    return out;
}

Value Value::array(Array items) {
    Value out;
    out.type_ = Type::Array;
    out.array_ = std::make_shared<Array>(std::move(items));
    return out;
}

Value Value::object(Object members) {
    Value out;
    out.type_ = Type::Object;
    out.object_ = std::make_shared<Object>(std::move(members));
    return out;
}

bool Value::as_bool() const {
    if (type_ != Type::Bool) {
        throw TypeError("JSON value is not a bool");
    }
    return bool_;
}

std::int64_t Value::as_int() const {
    switch (type_) {
        case Type::Int: return int_;
        case Type::UInt: return static_cast<std::int64_t>(uint_);
        case Type::Bool: return bool_ ? 1 : 0;
        default: throw TypeError("JSON value is not an integer");
    }
}

std::uint64_t Value::as_uint() const {
    switch (type_) {
        case Type::Int: return static_cast<std::uint64_t>(int_);
        case Type::UInt: return uint_;
        case Type::Bool: return bool_ ? 1U : 0U;
        default: throw TypeError("JSON value is not an integer");
    }
}

double Value::as_double() const {
    switch (type_) {
        case Type::Int: return static_cast<double>(int_);
        case Type::UInt: return static_cast<double>(uint_);
        case Type::Double: return double_;
        case Type::Bool: return bool_ ? 1.0 : 0.0;
        default: throw TypeError("JSON value is not a number");
    }
}

const std::string& Value::as_string() const {
    if (type_ != Type::String && type_ != Type::Bytes) {
        throw TypeError("JSON value is not a string");
    }
    return string_;
}

const Value::Array& Value::as_array() const {
    if (type_ != Type::Array) {
        throw TypeError("JSON value is not an array");
    }
    return *array_;
}

const Value::Object& Value::as_object() const {
    if (type_ != Type::Object) {
        throw TypeError("JSON value is not an object");
    }
    return *object_;
}

const Value* Value::find(std::string_view key) const {
    if (type_ != Type::Object) {
        return nullptr;
    }
    for (const auto& [name, value] : *object_) {
        if (name == key) {
            return &value;
        }
    }
    return nullptr;
}

Value& Value::operator[](std::string_view key) {
    if (type_ == Type::Null) {
        *this = object();
    }
    if (type_ != Type::Object) {
        throw TypeError("JSON value is not an object");
    }
    if (object_.use_count() > 1) {
        object_ = std::make_shared<Object>(*object_);
    }
    for (auto& [name, value] : *object_) {
        if (name == key) {
            return value;
        }
    }
    object_->emplace_back(std::string(key), Value());
    return object_->back().second;
}

void Value::push_back(Value item) {
    if (type_ == Type::Null) {
        *this = array({});
    }
    if (type_ != Type::Array) {
        throw TypeError("JSON value is not an array");
    }
    if (array_.use_count() > 1) {
        array_ = std::make_shared<Array>(*array_);
    }
    array_->push_back(std::move(item));
}

bool operator==(const Value& a, const Value& b) {
    if (a.is_number() && b.is_number() && a.type_ != b.type_) {
        if (a.type_ == Type::Double || b.type_ == Type::Double) {
            return a.as_double() == b.as_double();
        }
        // Int against UInt.
        const std::int64_t i = a.type_ == Type::Int ? a.int_ : b.int_;
        const std::uint64_t u = a.type_ == Type::UInt ? a.uint_ : b.uint_;
        return i >= 0 && static_cast<std::uint64_t>(i) == u;
    }
    if (a.type_ != b.type_) {
        return false;
    }
    switch (a.type_) {
        case Type::Null: return true;
        case Type::Bool: return a.bool_ == b.bool_;
        case Type::Int: return a.int_ == b.int_;
        case Type::UInt: return a.uint_ == b.uint_;
        case Type::Double: return a.double_ == b.double_;
        case Type::String:
        case Type::Bytes: return a.string_ == b.string_;
        case Type::Array: return *a.array_ == *b.array_;
        case Type::Object: return *a.object_ == *b.object_;
    }
    return false;
}

namespace {

void append_utf8(std::string& out, std::uint32_t code) {
    if (code < 0x80) {
        out.push_back(static_cast<char>(code));
    } else if (code < 0x800) {
        out.push_back(static_cast<char>(0xC0 | (code >> 6)));
        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
    } else if (code < 0x10000) {
        out.push_back(static_cast<char>(0xE0 | (code >> 12)));
        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
    } else {
        out.push_back(static_cast<char>(0xF0 | (code >> 18)));
        out.push_back(static_cast<char>(0x80 | ((code >> 12) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
    }
}

// A parser for RFC 8259 JSON text with the choices of simplejson.loads that
// cooler's Cooler.info relies on (verified against simplejson 3.19.3 on
// Python 3.12.7, tests/test_cpython_oracles.cpp and
// tests/fixtures/cpython/json.txt):
//   * the NaN, Infinity and -Infinity literals are errors (allow_nan is off);
//   * control characters below U+0020 inside strings are errors;
//   * a repeated object key keeps its first position and takes the last value;
//   * integers outside the 64 bit range become doubles (Python keeps them
//     exact; nothing downstream needs more);
//   * a single U+FEFF in front of the text is skipped; only space, tab, line
//     feed and carriage return are whitespace;
//   * nesting deeper than kMaxDepth is an error (Python raises RecursionError
//     somewhere beyond a thousand levels).
class Reader {
  public:
    explicit Reader(std::string_view text) : in_(text) {}

    std::optional<Value> read_document() {
        Value result;
        if (accept_word("\xEF\xBB\xBF")) {
            // one byte order mark in front of the text is skipped
        }
        skip_space();
        if (!read_value(result)) {
            return std::nullopt;
        }
        skip_space();
        if (at_ != in_.size()) {
            return std::nullopt;
        }
        return result;
    }

  private:
    static constexpr int kMaxDepth = 5000;

    [[nodiscard]] bool done() const { return at_ >= in_.size(); }
    [[nodiscard]] char peek() const { return in_[at_]; }

    void skip_space() {
        while (!done() && (peek() == ' ' || peek() == '\t' || peek() == '\n' || peek() == '\r')) {
            ++at_;
        }
    }

    bool accept(char c) {
        if (!done() && peek() == c) {
            ++at_;
            return true;
        }
        return false;
    }

    bool accept_word(std::string_view word) {
        if (in_.substr(at_, word.size()) == word) {
            at_ += word.size();
            return true;
        }
        return false;
    }

    bool read_value(Value& out) {
        if (done()) {
            return false;
        }
        switch (peek()) {
            case '"': {
                ++at_;
                std::string text;
                if (!read_string_rest(text)) {
                    return false;
                }
                out = Value(std::move(text));
                return true;
            }
            case '{': return read_object(out);
            case '[': return read_array(out);
            case 'n':
                if (accept_word("null")) {
                    out = Value();
                    return true;
                }
                return false;
            case 't':
                if (accept_word("true")) {
                    out = Value(true);
                    return true;
                }
                return false;
            case 'f':
                if (accept_word("false")) {
                    out = Value(false);
                    return true;
                }
                return false;
            default: return read_number(out);
        }
    }

    // [ value { , value } ]
    bool read_array(Value& out) {
        if (++depth_ > kMaxDepth) {
            return false;
        }
        ++at_;  // '['
        Value::Array items;
        skip_space();
        if (!accept(']')) {
            do {
                skip_space();
                Value item;
                if (!read_value(item)) {
                    return false;
                }
                items.push_back(std::move(item));
                skip_space();
            } while (accept(','));
            if (!accept(']')) {
                return false;
            }
        }
        --depth_;
        out = Value::array(std::move(items));
        return true;
    }

    // { string : value { , string : value } }
    bool read_object(Value& out) {
        if (++depth_ > kMaxDepth) {
            return false;
        }
        ++at_;  // '{'
        Value::Object members;
        skip_space();
        if (!accept('}')) {
            do {
                skip_space();
                std::string key;
                if (!accept('"') || !read_string_rest(key)) {
                    return false;
                }
                skip_space();
                if (!accept(':')) {
                    return false;
                }
                skip_space();
                Value item;
                if (!read_value(item)) {
                    return false;
                }
                const auto same_key = [&key](const auto& member) { return member.first == key; };
                if (auto found = std::find_if(members.begin(), members.end(), same_key);
                    found != members.end()) {
                    found->second = std::move(item);
                } else {
                    members.emplace_back(std::move(key), std::move(item));
                }
                skip_space();
            } while (accept(','));
            if (!accept('}')) {
                return false;
            }
        }
        --depth_;
        out = Value::object(std::move(members));
        return true;
    }

    [[nodiscard]] bool digit_at(std::size_t i) const {
        return i < in_.size() && in_[i] >= '0' && in_[i] <= '9';
    }

    // -? ( 0 | [1-9] digits ) ( . digits )? ( [eE] [+-]? digits )?
    // A fraction or exponent that lacks its digits is not consumed: the number
    // ends before it and the caller meets the stray character.
    bool read_number(Value& out) {
        const std::size_t first = at_;
        std::size_t i = at_;
        if (i < in_.size() && in_[i] == '-') {
            ++i;
        }
        if (!digit_at(i)) {
            return false;
        }
        if (in_[i] == '0') {
            ++i;
        } else {
            while (digit_at(i)) {
                ++i;
            }
        }
        bool integral = true;
        if (i < in_.size() && in_[i] == '.' && digit_at(i + 1)) {
            integral = false;
            for (i += 1; digit_at(i); ++i) {
            }
        }
        if (i < in_.size() && (in_[i] == 'e' || in_[i] == 'E')) {
            std::size_t j = i + 1;
            if (j < in_.size() && (in_[j] == '+' || in_[j] == '-')) {
                ++j;
            }
            if (digit_at(j)) {
                integral = false;
                for (; digit_at(j); ++j) {
                }
                i = j;
            }
        }
        at_ = i;
        const char* begin = in_.data() + first;
        const char* end = in_.data() + i;
        if (integral) {
            std::int64_t signed_value = 0;
            if (std::from_chars(begin, end, signed_value).ec == std::errc()) {
                out = Value(signed_value);
                return true;
            }
            std::uint64_t unsigned_value = 0;
            if (*begin != '-' && std::from_chars(begin, end, unsigned_value).ec == std::errc()) {
                out = Value(unsigned_value);
                return true;
            }
        }
        double number = 0.0;
        if (std::from_chars(begin, end, number).ec == std::errc::result_out_of_range) {
            number = std::strtod(std::string(begin, end).c_str(), nullptr);  // inf or 0
        }
        out = Value(number);
        return true;
    }

    // Four hexadecimal digits of a \u escape.
    bool read_hex4(std::uint32_t& code) {
        if (in_.size() - at_ < 4) {
            return false;
        }
        code = 0;
        for (int k = 0; k < 4; ++k) {
            const char c = in_[at_ + static_cast<std::size_t>(k)];
            int digit = -1;
            if (c >= '0' && c <= '9') {
                digit = c - '0';
            } else if (c >= 'a' && c <= 'f') {
                digit = c - 'a' + 10;
            } else if (c >= 'A' && c <= 'F') {
                digit = c - 'A' + 10;
            }
            if (digit < 0) {
                return false;
            }
            code = code * 16 + static_cast<std::uint32_t>(digit);
        }
        at_ += 4;
        return true;
    }

    // The characters of a string after its opening quote, up to and including
    // the closing one. Unescaped bytes of UTF-8 text pass through unchanged.
    bool read_string_rest(std::string& out) {
        while (!done()) {
            const char c = in_[at_++];
            if (c == '"') {
                return true;
            }
            if (static_cast<unsigned char>(c) < 0x20) {
                return false;
            }
            if (c != '\\') {
                out.push_back(c);
                continue;
            }
            if (done()) {
                return false;
            }
            switch (in_[at_++]) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    std::uint32_t code = 0;
                    if (!read_hex4(code)) {
                        return false;
                    }
                    // A high surrogate followed by an escaped low surrogate is
                    // one code point; a lone surrogate is kept as such.
                    if (code >= 0xD800 && code <= 0xDBFF && in_.substr(at_, 2) == "\\u") {
                        const std::size_t resume = at_;
                        at_ += 2;
                        std::uint32_t low = 0;
                        if (read_hex4(low) && low >= 0xDC00 && low <= 0xDFFF) {
                            code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                        } else {
                            at_ = resume;
                        }
                    }
                    append_utf8(out, code);
                    break;
                }
                default: return false;
            }
        }
        return false;  // unterminated
    }

    std::string_view in_;
    std::size_t at_ = 0;
    int depth_ = 0;
};

void dump_string(std::string& out, const std::string& text) {
    out.push_back('"');
    std::size_t i = 0;
    const auto hex = [&out](std::uint32_t code) {
        static constexpr char kDigits[] = "0123456789abcdef";
        out += "\\u";
        out.push_back(kDigits[(code >> 12) & 0xF]);
        out.push_back(kDigits[(code >> 8) & 0xF]);
        out.push_back(kDigits[(code >> 4) & 0xF]);
        out.push_back(kDigits[code & 0xF]);
    };
    while (i < text.size()) {
        const auto byte = static_cast<unsigned char>(text[i]);
        std::uint32_t code = byte;
        std::size_t length = 1;
        if (byte >= 0x80) {
            // Decode one UTF-8 sequence; an invalid byte is emitted as is.
            if ((byte & 0xE0) == 0xC0 && i + 1 < text.size()) {
                code = ((byte & 0x1FU) << 6) | (static_cast<unsigned char>(text[i + 1]) & 0x3FU);
                length = 2;
            } else if ((byte & 0xF0) == 0xE0 && i + 2 < text.size()) {
                code = ((byte & 0x0FU) << 12) |
                       ((static_cast<unsigned char>(text[i + 1]) & 0x3FU) << 6) |
                       (static_cast<unsigned char>(text[i + 2]) & 0x3FU);
                length = 3;
            } else if ((byte & 0xF8) == 0xF0 && i + 3 < text.size()) {
                code = ((byte & 0x07U) << 18) |
                       ((static_cast<unsigned char>(text[i + 1]) & 0x3FU) << 12) |
                       ((static_cast<unsigned char>(text[i + 2]) & 0x3FU) << 6) |
                       (static_cast<unsigned char>(text[i + 3]) & 0x3FU);
                length = 4;
            }
        }
        i += length;
        switch (code) {
            case '"': out += "\\\""; continue;
            case '\\': out += "\\\\"; continue;
            case '\n': out += "\\n"; continue;
            case '\r': out += "\\r"; continue;
            case '\t': out += "\\t"; continue;
            case '\b': out += "\\b"; continue;
            case '\f': out += "\\f"; continue;
            default: break;
        }
        if (code >= 0x20 && code <= 0x7E) {
            out.push_back(static_cast<char>(code));
        } else if (code >= 0x10000) {
            const std::uint32_t n = code - 0x10000;
            hex(0xD800 | ((n >> 10) & 0x3FF));
            hex(0xDC00 | (n & 0x3FF));
        } else {
            hex(code);
        }
    }
    out.push_back('"');
}

void dump(std::string& out, const Value& value) {
    switch (value.type()) {
        case Type::Null: out += "null"; return;
        case Type::Bool: out += value.as_bool() ? "true" : "false"; return;
        case Type::Int: out += std::to_string(value.as_int()); return;
        case Type::UInt: out += std::to_string(value.as_uint()); return;
        case Type::Double: {
            const double d = value.as_double();
            if (std::isnan(d)) {
                out += "NaN";
            } else if (std::isinf(d)) {
                out += d > 0 ? "Infinity" : "-Infinity";
            } else {
                out += npy::float_repr(d);
            }
            return;
        }
        case Type::String:
        case Type::Bytes: dump_string(out, value.as_string()); return;
        case Type::Array: {
            out.push_back('[');
            bool first = true;
            for (const Value& item : value.as_array()) {
                if (!first) {
                    out += ", ";
                }
                first = false;
                dump(out, item);
            }
            out.push_back(']');
            return;
        }
        case Type::Object: {
            out.push_back('{');
            bool first = true;
            for (const auto& [key, item] : value.as_object()) {
                if (!first) {
                    out += ", ";
                }
                first = false;
                dump_string(out, key);
                out += ": ";
                dump(out, item);
            }
            out.push_back('}');
            return;
        }
    }
}

}  // namespace

std::optional<Value> parse(std::string_view text) {
    return Reader(text).read_document();
}

std::string dumps(const Value& value) {
    std::string out;
    dump(out, value);
    return out;
}

}  // namespace coolercpp::json
