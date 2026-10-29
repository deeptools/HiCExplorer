// The element types of coolercpp's columnar tables, named after numpy.

#ifndef COOLERCPP_DTYPE_HPP
#define COOLERCPP_DTYPE_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <type_traits>

namespace coolercpp {

// numpy.bool_ occupies one byte. std::vector<bool> is a bitset without a
// contiguous buffer, so boolean columns store this type instead.
struct Bool8 {
    std::uint8_t value = 0;
    constexpr Bool8() = default;
    constexpr Bool8(bool v) : value(v ? 1 : 0) {}  // NOLINT(google-explicit-constructor)
    constexpr operator bool() const { return value != 0; }  // NOLINT
    friend constexpr bool operator==(Bool8 a, Bool8 b) { return a.value == b.value; }
};

enum class DType : std::uint8_t {
    Bool,
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float32,
    Float64,
    // A pandas object column holding str.
    String,
    // pandas.Categorical: integer codes into an ordered list of labels. This is
    // how cooler returns the chrom column of the bin table.
    Categorical,
};

// numpy's name for the type ("int32", "float64", "bool"); "object" for
// String and "category" for Categorical, as pandas prints them.
[[nodiscard]] std::string_view dtype_name(DType dtype) noexcept;

// Parses the spellings numpy.dtype() accepts for the supported types: "int32",
// "i4", "<i4", "int" (int64 on LP64), "float" (float64), "f8", "double",
// "uint16", "u2", "bool", "?", plus "object"/"str" and "category".
// Throws TypeError for anything else.
[[nodiscard]] DType dtype_from_name(std::string_view name);

[[nodiscard]] std::size_t itemsize(DType dtype) noexcept;  // 0 for String, Categorical
[[nodiscard]] bool is_integer(DType dtype) noexcept;
[[nodiscard]] bool is_signed_integer(DType dtype) noexcept;
[[nodiscard]] bool is_unsigned_integer(DType dtype) noexcept;
[[nodiscard]] bool is_float(DType dtype) noexcept;
// Bool, integers and floats.
[[nodiscard]] bool is_numeric(DType dtype) noexcept;

// numpy 1.26 numpy.result_type for two arrays (no value based casting):
// int32 + float32 -> float64, uint8 + int8 -> int16, uint64 + int64 -> float64.
[[nodiscard]] DType result_type(DType a, DType b);

// The result of a binary arithmetic operation between an array of `array`
// and a Python scalar of `scalar_kind` under numpy 1.26 value based casting,
// for the two cases cooler uses: `1 / weights` (a Python int, then true
// division) and `x * float`. True division of an integer array yields
// float64; a float32 array stays float32.
[[nodiscard]] DType true_divide_result(DType array) noexcept;

template <typename T>
struct dtype_of;
template <>
struct dtype_of<Bool8> {
    static constexpr DType value = DType::Bool;
};
template <>
struct dtype_of<std::int8_t> {
    static constexpr DType value = DType::Int8;
};
template <>
struct dtype_of<std::int16_t> {
    static constexpr DType value = DType::Int16;
};
template <>
struct dtype_of<std::int32_t> {
    static constexpr DType value = DType::Int32;
};
template <>
struct dtype_of<std::int64_t> {
    static constexpr DType value = DType::Int64;
};
template <>
struct dtype_of<std::uint8_t> {
    static constexpr DType value = DType::UInt8;
};
template <>
struct dtype_of<std::uint16_t> {
    static constexpr DType value = DType::UInt16;
};
template <>
struct dtype_of<std::uint32_t> {
    static constexpr DType value = DType::UInt32;
};
template <>
struct dtype_of<std::uint64_t> {
    static constexpr DType value = DType::UInt64;
};
template <>
struct dtype_of<float> {
    static constexpr DType value = DType::Float32;
};
template <>
struct dtype_of<double> {
    static constexpr DType value = DType::Float64;
};
template <>
struct dtype_of<std::string> {
    static constexpr DType value = DType::String;
};

template <typename T>
inline constexpr DType dtype_v = dtype_of<T>::value;

}  // namespace coolercpp

#endif  // COOLERCPP_DTYPE_HPP
