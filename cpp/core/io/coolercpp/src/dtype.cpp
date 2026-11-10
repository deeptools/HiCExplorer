#include "coolercpp/dtype.hpp"

#include <algorithm>
#include <array>

#include "coolercpp/errors.hpp"

namespace coolercpp {

std::string_view dtype_name(DType dtype) noexcept {
    switch (dtype) {
        case DType::Bool: return "bool";
        case DType::Int8: return "int8";
        case DType::Int16: return "int16";
        case DType::Int32: return "int32";
        case DType::Int64: return "int64";
        case DType::UInt8: return "uint8";
        case DType::UInt16: return "uint16";
        case DType::UInt32: return "uint32";
        case DType::UInt64: return "uint64";
        case DType::Float32: return "float32";
        case DType::Float64: return "float64";
        case DType::String: return "object";
        case DType::Categorical: return "category";
    }
    return "unknown";
}

DType dtype_from_name(std::string_view name) {
    std::string text(name);
    // A byte order prefix is accepted when it is the native (little endian) or
    // the "not applicable" one.
    if (!text.empty() && (text[0] == '<' || text[0] == '=' || text[0] == '|')) {
        text.erase(0, 1);
    }
    struct Entry {
        std::string_view name;
        DType dtype;
    };
    static constexpr std::array<Entry, 40> kNames{{
        {"bool", DType::Bool},       {"bool_", DType::Bool},     {"?", DType::Bool},
        {"b1", DType::Bool},         {"int8", DType::Int8},      {"i1", DType::Int8},
        {"int16", DType::Int16},     {"i2", DType::Int16},       {"int32", DType::Int32},
        {"i4", DType::Int32},        {"int64", DType::Int64},    {"i8", DType::Int64},
        {"int", DType::Int64},       {"uint8", DType::UInt8},    {"u1", DType::UInt8},
        {"uint16", DType::UInt16},   {"u2", DType::UInt16},      {"uint32", DType::UInt32},
        {"u4", DType::UInt32},       {"uint64", DType::UInt64},  {"u8", DType::UInt64},
        {"uint", DType::UInt64},     {"float32", DType::Float32}, {"f4", DType::Float32},
        {"single", DType::Float32},  {"float64", DType::Float64}, {"f8", DType::Float64},
        {"float", DType::Float64},   {"double", DType::Float64}, {"object", DType::String},
        {"O", DType::String},        {"str", DType::String},     {"category", DType::Categorical},
        {"int_", DType::Int64},      {"intc", DType::Int32},     {"uintc", DType::UInt32},
        {"short", DType::Int16},     {"ushort", DType::UInt16},  {"byte", DType::Int8},
        {"ubyte", DType::UInt8},
    }};
    for (const Entry& entry : kNames) {
        if (entry.name == text) {
            return entry.dtype;
        }
    }
    throw TypeError("data type '" + std::string(name) + "' not understood");
}

std::size_t itemsize(DType dtype) noexcept {
    switch (dtype) {
        case DType::Bool:
        case DType::Int8:
        case DType::UInt8: return 1;
        case DType::Int16:
        case DType::UInt16: return 2;
        case DType::Int32:
        case DType::UInt32:
        case DType::Float32: return 4;
        case DType::Int64:
        case DType::UInt64:
        case DType::Float64: return 8;
        case DType::String:
        case DType::Categorical: return 0;
    }
    return 0;
}

bool is_signed_integer(DType dtype) noexcept {
    return dtype == DType::Int8 || dtype == DType::Int16 || dtype == DType::Int32 ||
           dtype == DType::Int64;
}

bool is_unsigned_integer(DType dtype) noexcept {
    return dtype == DType::UInt8 || dtype == DType::UInt16 || dtype == DType::UInt32 ||
           dtype == DType::UInt64;
}

bool is_integer(DType dtype) noexcept {
    return is_signed_integer(dtype) || is_unsigned_integer(dtype);
}

bool is_float(DType dtype) noexcept {
    return dtype == DType::Float32 || dtype == DType::Float64;
}

bool is_numeric(DType dtype) noexcept {
    return dtype == DType::Bool || is_integer(dtype) || is_float(dtype);
}

namespace {

DType signed_of_size(std::size_t bytes) {
    switch (bytes) {
        case 1: return DType::Int8;
        case 2: return DType::Int16;
        case 4: return DType::Int32;
        default: return DType::Int64;
    }
}

}  // namespace

DType result_type(DType a, DType b) {
    if (!is_numeric(a) || !is_numeric(b)) {
        throw TypeError("no arithmetic promotion between " + std::string(dtype_name(a)) +
                        " and " + std::string(dtype_name(b)));
    }
    if (a == b) {
        return a;
    }
    if (a == DType::Bool) {
        return b;
    }
    if (b == DType::Bool) {
        return a;
    }
    if (is_float(a) || is_float(b)) {
        if (a == DType::Float64 || b == DType::Float64) {
            return DType::Float64;
        }
        // One float32 and one integer: float32 holds integers of up to 16 bits
        // exactly, anything wider promotes to float64.
        const DType other = is_float(a) ? b : a;
        return itemsize(other) <= 2 ? DType::Float32 : DType::Float64;
    }
    const bool a_signed = is_signed_integer(a);
    const bool b_signed = is_signed_integer(b);
    if (a_signed == b_signed) {
        return itemsize(a) >= itemsize(b) ? a : b;
    }
    const DType s = a_signed ? a : b;
    const DType u = a_signed ? b : a;
    if (itemsize(s) > itemsize(u)) {
        return s;
    }
    if (u == DType::UInt64) {
        return DType::Float64;
    }
    return signed_of_size(itemsize(u) * 2);
}

DType true_divide_result(DType array) noexcept {
    return array == DType::Float32 ? DType::Float32 : DType::Float64;
}

}  // namespace coolercpp
