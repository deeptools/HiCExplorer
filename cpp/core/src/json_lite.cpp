#include "hicx/json_lite.hpp"

#include <cctype>
#include <cstdlib>
#include <sstream>

#include "hicx/numpy_compat.hpp"

namespace hicx::json {

namespace {

class Parser {
  public:
    explicit Parser(const std::string& text) : text_(text) {}

    std::optional<Value> parse_document() {
        skip_whitespace();
        std::optional<Value> value = parse_value();
        if (!value.has_value()) {
            return std::nullopt;
        }
        skip_whitespace();
        if (position_ != text_.size()) {
            return std::nullopt;
        }
        return value;
    }

  private:
    void skip_whitespace() {
        while (position_ < text_.size() &&
               (text_[position_] == ' ' || text_[position_] == '\t' ||
                text_[position_] == '\n' || text_[position_] == '\r')) {
            ++position_;
        }
    }

    bool consume(char expected) {
        if (position_ < text_.size() && text_[position_] == expected) {
            ++position_;
            return true;
        }
        return false;
    }

    bool consume_literal(const char* literal) {
        const std::size_t length = std::char_traits<char>::length(literal);
        if (text_.compare(position_, length, literal) == 0) {
            position_ += length;
            return true;
        }
        return false;
    }

    std::optional<Value> parse_value() {
        skip_whitespace();
        if (position_ >= text_.size()) {
            return std::nullopt;
        }
        const char c = text_[position_];
        if (c == '{') {
            return parse_object();
        }
        if (c == '[') {
            return parse_array();
        }
        if (c == '"') {
            std::optional<std::string> text = parse_string();
            if (!text.has_value()) {
                return std::nullopt;
            }
            return Value::string(*text);
        }
        if (consume_literal("null")) {
            return Value::null();
        }
        if (consume_literal("true")) {
            return Value::boolean(true);
        }
        if (consume_literal("false")) {
            return Value::boolean(false);
        }
        return parse_number();
    }

    std::optional<Value> parse_object() {
        if (!consume('{')) {
            return std::nullopt;
        }
        Object object;
        skip_whitespace();
        if (consume('}')) {
            return Value::object(std::move(object));
        }
        while (true) {
            skip_whitespace();
            std::optional<std::string> key = parse_string();
            if (!key.has_value()) {
                return std::nullopt;
            }
            skip_whitespace();
            if (!consume(':')) {
                return std::nullopt;
            }
            std::optional<Value> value = parse_value();
            if (!value.has_value()) {
                return std::nullopt;
            }
            object.emplace(*key, *value);
            skip_whitespace();
            if (consume(',')) {
                continue;
            }
            if (consume('}')) {
                return Value::object(std::move(object));
            }
            return std::nullopt;
        }
    }

    std::optional<Value> parse_array() {
        if (!consume('[')) {
            return std::nullopt;
        }
        Array array;
        skip_whitespace();
        if (consume(']')) {
            return Value::array(std::move(array));
        }
        while (true) {
            std::optional<Value> value = parse_value();
            if (!value.has_value()) {
                return std::nullopt;
            }
            array.push_back(*value);
            skip_whitespace();
            if (consume(',')) {
                continue;
            }
            if (consume(']')) {
                return Value::array(std::move(array));
            }
            return std::nullopt;
        }
    }

    std::optional<std::string> parse_string() {
        if (!consume('"')) {
            return std::nullopt;
        }
        std::string out;
        while (position_ < text_.size()) {
            const char c = text_[position_++];
            if (c == '"') {
                return out;
            }
            if (c != '\\') {
                out.push_back(c);
                continue;
            }
            if (position_ >= text_.size()) {
                return std::nullopt;
            }
            const char escape = text_[position_++];
            switch (escape) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    if (position_ + 4 > text_.size()) {
                        return std::nullopt;
                    }
                    const std::string hex = text_.substr(position_, 4);
                    position_ += 4;
                    const long code = std::strtol(hex.c_str(), nullptr, 16);
                    // Only the BMP subset that fits into UTF-8 is needed here.
                    if (code < 0x80) {
                        out.push_back(static_cast<char>(code));
                    } else if (code < 0x800) {
                        out.push_back(static_cast<char>(0xC0 | (code >> 6)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    } else {
                        out.push_back(static_cast<char>(0xE0 | (code >> 12)));
                        out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    break;
                }
                default: return std::nullopt;
            }
        }
        return std::nullopt;
    }

    std::optional<Value> parse_number() {
        const std::size_t start = position_;
        if (position_ < text_.size() && (text_[position_] == '-' || text_[position_] == '+')) {
            ++position_;
        }
        bool has_digits = false;
        while (position_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[position_])) != 0) {
            ++position_;
            has_digits = true;
        }
        bool is_double = false;
        if (position_ < text_.size() && text_[position_] == '.') {
            is_double = true;
            ++position_;
            while (position_ < text_.size() &&
                   std::isdigit(static_cast<unsigned char>(text_[position_])) != 0) {
                ++position_;
                has_digits = true;
            }
        }
        if (position_ < text_.size() && (text_[position_] == 'e' || text_[position_] == 'E')) {
            is_double = true;
            ++position_;
            if (position_ < text_.size() &&
                (text_[position_] == '-' || text_[position_] == '+')) {
                ++position_;
            }
            while (position_ < text_.size() &&
                   std::isdigit(static_cast<unsigned char>(text_[position_])) != 0) {
                ++position_;
            }
        }
        if (!has_digits) {
            position_ = start;
            return std::nullopt;
        }
        const std::string token = text_.substr(start, position_ - start);
        if (is_double) {
            return Value::number(std::strtod(token.c_str(), nullptr));
        }
        return Value::integer(static_cast<std::int64_t>(std::strtoll(token.c_str(), nullptr, 10)));
    }

    const std::string& text_;
    std::size_t position_ = 0;
};

}  // namespace

Value Value::null() { return Value(); }

Value Value::boolean(bool value) {
    Value out;
    out.type_ = Type::Bool;
    out.bool_value_ = value;
    return out;
}

Value Value::integer(std::int64_t value) {
    Value out;
    out.type_ = Type::Int;
    out.int_value_ = value;
    return out;
}

Value Value::number(double value) {
    Value out;
    out.type_ = Type::Double;
    out.double_value_ = value;
    return out;
}

Value Value::string(std::string value) {
    Value out;
    out.type_ = Type::String;
    out.string_value_ = std::move(value);
    return out;
}

Value Value::array(Array value) {
    Value out;
    out.type_ = Type::Array;
    out.array_value_ = std::make_shared<Array>(std::move(value));
    return out;
}

Value Value::object(Object value) {
    Value out;
    out.type_ = Type::Object;
    out.object_value_ = std::make_shared<Object>(std::move(value));
    return out;
}

const Value* Value::find(const std::string& key) const {
    if (type_ != Type::Object || object_value_ == nullptr) {
        return nullptr;
    }
    const auto it = object_value_->find(key);
    return it == object_value_->end() ? nullptr : &it->second;
}

std::string Value::to_python_string() const {
    switch (type_) {
        case Type::Null: return "None";
        case Type::Bool: return bool_value_ ? "True" : "False";
        case Type::Int: return std::to_string(int_value_);
        case Type::Double: return npy::float_repr(double_value_);
        case Type::String: return string_value_;
        case Type::Array: {
            std::string out = "[";
            for (std::size_t i = 0; i < array_value_->size(); ++i) {
                if (i > 0) {
                    out += ", ";
                }
                const Value& item = (*array_value_)[i];
                out += item.is_string() ? "'" + item.as_string() + "'"
                                        : item.to_python_string();
            }
            return out + "]";
        }
        case Type::Object: {
            std::string out = "{";
            bool first = true;
            for (const auto& [key, item] : *object_value_) {
                if (!first) {
                    out += ", ";
                }
                first = false;
                out += "'" + key + "': ";
                out += item.is_string() ? "'" + item.as_string() + "'"
                                        : item.to_python_string();
            }
            return out + "}";
        }
    }
    return {};
}

std::optional<Value> parse(const std::string& text) {
    Parser parser(text);
    return parser.parse_document();
}

std::string dump_object(
    const std::vector<std::pair<std::string, std::string>>& fields) {
    // Only the escapes json.dumps emits for the characters that can occur in
    // a cooler metadata field are handled; anything else is passed through as
    // UTF-8, which is what json.dumps with ensure_ascii=False would do and
    // what the values in this dictionary always are.
    const auto quote = [](const std::string& text) {
        std::string out = "\"";
        for (const char character : text) {
            switch (character) {
                case '"': out += "\\\""; break;
                case '\\': out += "\\\\"; break;
                case '\n': out += "\\n"; break;
                case '\r': out += "\\r"; break;
                case '\t': out += "\\t"; break;
                default: out += character; break;
            }
        }
        return out + "\"";
    };

    std::string out = "{";
    bool first = true;
    for (const auto& [key, value] : fields) {
        if (!first) {
            out += ", ";
        }
        first = false;
        out += quote(key) + ": " + quote(value);
    }
    return out + "}";
}

}  // namespace hicx::json
