// JSON values as cooler reads and writes them.
//
// cooler keeps every piece of metadata in HDF5 attributes. Cooler.info runs
// simplejson.loads over each attribute that h5py returns as a str and keeps
// the others as numpy values, and create_cooler stores the user metadata as
// simplejson.dumps(metadata). json::Value covers both sides: JSON types with
// ordered objects (Python dicts keep insertion order, and the key order is
// part of the dumped string), plus the numpy scalar type of numbers that came
// straight from a numeric HDF5 attribute, plus bytes for fixed length string
// attributes (h5py returns numpy.bytes_ for those, which cooler never decodes).
//
// Derived from HiCExplorer v4 cpp/core/src/json_lite.cpp; reworked for
// ordered objects, strict simplejson parsing and simplejson.dumps output.

#ifndef COOLERCPP_JSON_HPP
#define COOLERCPP_JSON_HPP

#include <concepts>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "coolercpp/dtype.hpp"

namespace coolercpp::json {

enum class Type { Null, Bool, Int, UInt, Double, String, Bytes, Array, Object };

class Value {
  public:
    using Array = std::vector<Value>;
    using Object = std::vector<std::pair<std::string, Value>>;

    Value() = default;
    Value(std::nullptr_t) {}  // NOLINT(google-explicit-constructor)
    Value(bool value);        // NOLINT
    template <std::signed_integral T>
    Value(T value)  // NOLINT
        : type_(Type::Int), int_(static_cast<std::int64_t>(value)), dtype_(DType::Int64) {}
    template <std::unsigned_integral T>
        requires(!std::same_as<T, bool>)
    Value(T value)  // NOLINT
        : type_(Type::UInt), uint_(static_cast<std::uint64_t>(value)), dtype_(DType::UInt64) {}
    Value(double value);             // NOLINT
    Value(const char* value);        // NOLINT
    Value(std::string value);        // NOLINT
    Value(std::string_view value);   // NOLINT

    [[nodiscard]] static Value integer(std::int64_t value, DType dtype = DType::Int64);
    [[nodiscard]] static Value unsigned_integer(std::uint64_t value,
                                                DType dtype = DType::UInt64);
    [[nodiscard]] static Value number(double value, DType dtype = DType::Float64);
    [[nodiscard]] static Value bytes(std::string value);
    [[nodiscard]] static Value array(Array items);
    [[nodiscard]] static Value object(Object members = {});

    [[nodiscard]] Type type() const noexcept { return type_; }
    [[nodiscard]] bool is_null() const noexcept { return type_ == Type::Null; }
    [[nodiscard]] bool is_bool() const noexcept { return type_ == Type::Bool; }
    [[nodiscard]] bool is_integer() const noexcept {
        return type_ == Type::Int || type_ == Type::UInt;
    }
    [[nodiscard]] bool is_number() const noexcept {
        return is_integer() || type_ == Type::Double;
    }
    [[nodiscard]] bool is_string() const noexcept { return type_ == Type::String; }
    [[nodiscard]] bool is_bytes() const noexcept { return type_ == Type::Bytes; }
    [[nodiscard]] bool is_array() const noexcept { return type_ == Type::Array; }
    [[nodiscard]] bool is_object() const noexcept { return type_ == Type::Object; }

    // The numpy scalar type: Int64/UInt64/Float64 for values that came from
    // JSON text or C++ literals, the attribute's own type for HDF5 numbers,
    // Bool for booleans.
    [[nodiscard]] DType dtype() const noexcept { return dtype_; }

    [[nodiscard]] bool as_bool() const;
    [[nodiscard]] std::int64_t as_int() const;
    [[nodiscard]] std::uint64_t as_uint() const;
    [[nodiscard]] double as_double() const;
    // String or Bytes content.
    [[nodiscard]] const std::string& as_string() const;
    [[nodiscard]] const Array& as_array() const;
    [[nodiscard]] const Object& as_object() const;

    // Object member lookup; nullptr when absent or not an object.
    [[nodiscard]] const Value* find(std::string_view key) const;
    // Object member access that inserts a null member when absent (dict[key]
    // assignment keeps the position of an existing key).
    Value& operator[](std::string_view key);
    // Array append.
    void push_back(Value item);

    friend bool operator==(const Value& a, const Value& b);

  private:
    Type type_ = Type::Null;
    bool bool_ = false;
    std::int64_t int_ = 0;
    std::uint64_t uint_ = 0;
    double double_ = 0.0;
    std::string string_;
    std::shared_ptr<Array> array_;
    std::shared_ptr<Object> object_;
    DType dtype_ = DType::Float64;
};

// simplejson.loads: returns nothing where simplejson raises JSONDecodeError
// (the case Cooler.info swallows). Accepts NaN, Infinity and -Infinity,
// rejects control characters inside strings, keeps object member order with
// the last value of a repeated key at the key's first position.
[[nodiscard]] std::optional<Value> parse(std::string_view text);

// simplejson.dumps with its defaults: separators ", " and ": ", ensure_ascii,
// NaN/Infinity literals, repr() for floats, insertion order.
[[nodiscard]] std::string dumps(const Value& value);

}  // namespace coolercpp::json

#endif  // COOLERCPP_JSON_HPP
