// Minimal JSON reader.
//
// cooler stores every piece of file metadata as an HDF5 attribute and runs
// json.loads over string valued attributes when it builds Cooler.info. That is
// why 'bin-size' can come back as None from the literal string "null" and why
// 'metadata' turns into a dictionary. The reader below reproduces exactly that
// step; it is not meant as a general purpose JSON library.

#ifndef HICX_JSON_LITE_HPP
#define HICX_JSON_LITE_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace hicx::json {

class Value;
using Object = std::map<std::string, Value>;
using Array = std::vector<Value>;

enum class Type { Null, Bool, Int, Double, String, Array, Object };

class Value {
  public:
    Value() = default;
    static Value null();
    static Value boolean(bool value);
    static Value integer(std::int64_t value);
    static Value number(double value);
    static Value string(std::string value);
    static Value array(Array value);
    static Value object(Object value);

    [[nodiscard]] Type type() const noexcept { return type_; }
    [[nodiscard]] bool is_null() const noexcept { return type_ == Type::Null; }
    [[nodiscard]] bool is_string() const noexcept { return type_ == Type::String; }
    [[nodiscard]] bool is_object() const noexcept { return type_ == Type::Object; }
    [[nodiscard]] bool is_number() const noexcept {
        return type_ == Type::Int || type_ == Type::Double;
    }

    [[nodiscard]] bool as_bool() const { return bool_value_; }
    [[nodiscard]] std::int64_t as_int() const { return int_value_; }
    [[nodiscard]] double as_double() const {
        return type_ == Type::Int ? static_cast<double>(int_value_) : double_value_;
    }
    [[nodiscard]] const std::string& as_string() const { return string_value_; }
    [[nodiscard]] const Array& as_array() const { return *array_value_; }
    [[nodiscard]] const Object& as_object() const { return *object_value_; }

    [[nodiscard]] const Value* find(const std::string& key) const;

    // Python's str() of the parsed value, used when a metadata field is
    // interpolated into the tool output.
    [[nodiscard]] std::string to_python_string() const;

  private:
    Type type_ = Type::Null;
    bool bool_value_ = false;
    std::int64_t int_value_ = 0;
    double double_value_ = 0.0;
    std::string string_value_;
    std::shared_ptr<Array> array_value_;
    std::shared_ptr<Object> object_value_;
};

// Returns nothing when the text is not valid JSON, which is the case cooler
// swallows with `except ValueError: pass`.
[[nodiscard]] std::optional<Value> parse(const std::string& text);

}  // namespace hicx::json

#endif  // HICX_JSON_LITE_HPP
