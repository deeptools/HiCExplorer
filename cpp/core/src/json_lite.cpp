#include "hicx/json_lite.hpp"

#include <utility>

#include "coolercpp/json.hpp"
#include "hicx/numpy_compat.hpp"

namespace hicx::json {

namespace {

// The same value in this header's types: objects as sorted maps, integers
// beyond 64 bits as doubles.
Value convert(const coolercpp::json::Value& value) {
    namespace cj = coolercpp::json;
    switch (value.type()) {
        case cj::Type::Null: return Value::null();
        case cj::Type::Bool: return Value::boolean(value.as_bool());
        case cj::Type::Int: return Value::integer(value.as_int());
        case cj::Type::UInt: return Value::number(static_cast<double>(value.as_uint()));
        case cj::Type::Double: return Value::number(value.as_double());
        case cj::Type::String:
        case cj::Type::Bytes: return Value::string(value.as_string());
        case cj::Type::Array: {
            Array items;
            for (const cj::Value& item : value.as_array()) {
                items.push_back(convert(item));
            }
            return Value::array(std::move(items));
        }
        case cj::Type::Object: {
            Object members;
            for (const auto& [key, item] : value.as_object()) {
                members.emplace(key, convert(item));
            }
            return Value::object(std::move(members));
        }
    }
    return Value::null();
}

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
    const std::optional<coolercpp::json::Value> parsed = coolercpp::json::parse(text);
    if (!parsed.has_value()) {
        return std::nullopt;
    }
    return convert(*parsed);
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
