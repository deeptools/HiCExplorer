#include "coolercpp/json.hpp"

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

bool is_ws(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }
bool is_digit(char c) { return c >= '0' && c <= '9'; }

// A port of simplejson's pure Python scanner (simplejson/scanner.py and
// decoder.py) in strict mode.
class Parser {
  public:
    explicit Parser(std::string_view text) : text_(text) {}

    std::optional<Value> document() {
        skip_ws();
        std::optional<Value> value = scan();
        if (!value.has_value()) {
            return std::nullopt;
        }
        skip_ws();
        if (pos_ != text_.size()) {
            return std::nullopt;  // "Extra data"
        }
        return value;
    }

  private:
    void skip_ws() {
        while (pos_ < text_.size() && is_ws(text_[pos_])) {
            ++pos_;
        }
    }

    bool literal(std::string_view word) {
        if (text_.substr(pos_, word.size()) == word) {
            pos_ += word.size();
            return true;
        }
        return false;
    }

    std::optional<Value> scan() {
        if (pos_ >= text_.size()) {
            return std::nullopt;
        }
        const char c = text_[pos_];
        if (c == '"') {
            ++pos_;
            std::optional<std::string> s = string_body();
            if (!s.has_value()) {
                return std::nullopt;
            }
            return Value(std::move(*s));
        }
        if (c == '{') {
            ++pos_;
            return object_body();
        }
        if (c == '[') {
            ++pos_;
            return array_body();
        }
        if (c == 'n' && literal("null")) {
            return Value();
        }
        if (c == 't' && literal("true")) {
            return Value(true);
        }
        if (c == 'f' && literal("false")) {
            return Value(false);
        }
        if (c == 'N' && literal("NaN")) {
            return Value(std::numeric_limits<double>::quiet_NaN());
        }
        if (c == 'I' && literal("Infinity")) {
            return Value(std::numeric_limits<double>::infinity());
        }
        if (c == '-' && text_.substr(pos_, 9) == "-Infinity") {
            pos_ += 9;
            return Value(-std::numeric_limits<double>::infinity());
        }
        return number();
    }

    std::optional<Value> number() {
        // NUMBER_RE = (-?(?:0|[1-9]\d*))(\.\d+)?([eE][-+]?\d+)?
        const std::size_t start = pos_;
        std::size_t p = pos_;
        if (p < text_.size() && text_[p] == '-') {
            ++p;
        }
        if (p >= text_.size() || !is_digit(text_[p])) {
            return std::nullopt;
        }
        if (text_[p] == '0') {
            ++p;
        } else {
            while (p < text_.size() && is_digit(text_[p])) {
                ++p;
            }
        }
        std::size_t integer_end = p;
        bool is_float = false;
        if (p + 1 < text_.size() && text_[p] == '.' && is_digit(text_[p + 1])) {
            p += 1;
            while (p < text_.size() && is_digit(text_[p])) {
                ++p;
            }
            is_float = true;
        }
        if (p < text_.size() && (text_[p] == 'e' || text_[p] == 'E')) {
            std::size_t q = p + 1;
            if (q < text_.size() && (text_[q] == '-' || text_[q] == '+')) {
                ++q;
            }
            if (q < text_.size() && is_digit(text_[q])) {
                while (q < text_.size() && is_digit(text_[q])) {
                    ++q;
                }
                p = q;
                is_float = true;
            }
        }
        pos_ = p;
        const std::string token(text_.substr(start, p - start));
        if (is_float) {
            return Value(std::strtod(token.c_str(), nullptr));
        }
        (void)integer_end;
        errno = 0;
        char* end = nullptr;
        const long long value = std::strtoll(token.c_str(), &end, 10);
        if (errno == 0) {
            return Value(static_cast<std::int64_t>(value));
        }
        errno = 0;
        if (token[0] != '-') {
            const unsigned long long u = std::strtoull(token.c_str(), &end, 10);
            if (errno == 0) {
                return Value(static_cast<std::uint64_t>(u));
            }
        }
        // Python keeps arbitrary precision; beyond 64 bits the nearest double
        // is kept.
        return Value(std::strtod(token.c_str(), nullptr));
    }

    std::optional<std::uint32_t> hex4() {
        if (pos_ + 4 > text_.size()) {
            return std::nullopt;
        }
        std::uint32_t code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = text_[pos_ + static_cast<std::size_t>(i)];
            code <<= 4;
            if (c >= '0' && c <= '9') {
                code |= static_cast<std::uint32_t>(c - '0');
            } else if (c >= 'a' && c <= 'f') {
                code |= static_cast<std::uint32_t>(c - 'a' + 10);
            } else if (c >= 'A' && c <= 'F') {
                code |= static_cast<std::uint32_t>(c - 'A' + 10);
            } else {
                return std::nullopt;
            }
        }
        pos_ += 4;
        return code;
    }

    std::optional<std::string> string_body() {
        std::string out;
        while (true) {
            if (pos_ >= text_.size()) {
                return std::nullopt;  // Unterminated string
            }
            const char c = text_[pos_++];
            if (c == '"') {
                return out;
            }
            if (c != '\\') {
                if (static_cast<unsigned char>(c) < 0x20) {
                    return std::nullopt;  // Invalid control character (strict)
                }
                out.push_back(c);
                continue;
            }
            if (pos_ >= text_.size()) {
                return std::nullopt;
            }
            const char e = text_[pos_++];
            switch (e) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    std::optional<std::uint32_t> code = hex4();
                    if (!code.has_value()) {
                        return std::nullopt;
                    }
                    std::uint32_t uni = *code;
                    if (uni >= 0xD800 && uni <= 0xDBFF && pos_ + 6 <= text_.size() &&
                        text_[pos_] == '\\' && text_[pos_ + 1] == 'u') {
                        const std::size_t save = pos_;
                        pos_ += 2;
                        std::optional<std::uint32_t> low = hex4();
                        if (low.has_value() && *low >= 0xDC00 && *low <= 0xDFFF) {
                            uni = 0x10000 + (((uni - 0xD800) << 10) | (*low - 0xDC00));
                        } else {
                            pos_ = save;
                        }
                    }
                    append_utf8(out, uni);
                    break;
                }
                default: return std::nullopt;  // Invalid \X escape sequence
            }
        }
    }

    std::optional<Value> object_body() {
        Value::Object members;
        skip_ws();
        if (pos_ < text_.size() && text_[pos_] == '}') {
            ++pos_;
            return Value::object(std::move(members));
        }
        while (true) {
            if (pos_ >= text_.size() || text_[pos_] != '"') {
                return std::nullopt;
            }
            ++pos_;
            std::optional<std::string> key = string_body();
            if (!key.has_value()) {
                return std::nullopt;
            }
            skip_ws();
            if (pos_ >= text_.size() || text_[pos_] != ':') {
                return std::nullopt;
            }
            ++pos_;
            skip_ws();
            std::optional<Value> value = scan();
            if (!value.has_value()) {
                return std::nullopt;
            }
            bool replaced = false;
            for (auto& [name, existing] : members) {
                if (name == *key) {
                    existing = std::move(*value);
                    replaced = true;
                    break;
                }
            }
            if (!replaced) {
                members.emplace_back(std::move(*key), std::move(*value));
            }
            skip_ws();
            if (pos_ >= text_.size()) {
                return std::nullopt;
            }
            const char c = text_[pos_++];
            if (c == '}') {
                return Value::object(std::move(members));
            }
            if (c != ',') {
                return std::nullopt;
            }
            skip_ws();
        }
    }

    std::optional<Value> array_body() {
        Value::Array items;
        skip_ws();
        if (pos_ < text_.size() && text_[pos_] == ']') {
            ++pos_;
            return Value::array(std::move(items));
        }
        while (true) {
            std::optional<Value> value = scan();
            if (!value.has_value()) {
                return std::nullopt;
            }
            items.push_back(std::move(*value));
            skip_ws();
            if (pos_ >= text_.size()) {
                return std::nullopt;
            }
            const char c = text_[pos_++];
            if (c == ']') {
                return Value::array(std::move(items));
            }
            if (c != ',') {
                return std::nullopt;
            }
            skip_ws();
        }
    }

    std::string_view text_;
    std::size_t pos_ = 0;
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
    Parser parser(text);
    return parser.document();
}

std::string dumps(const Value& value) {
    std::string out;
    dump(out, value);
    return out;
}

}  // namespace coolercpp::json
