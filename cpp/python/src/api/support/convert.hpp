// Conversions between Python objects and C++ for hicexplorer._hicx: str,
// numpy arrays and scalars, and scipy.sparse matrices.
//
// Every bound module includes this header (through api.hpp). Rules the
// helpers implement, so that a binding never has to repeat them:
//  * a str crosses as UTF-8 with the surrogatepass handler, so any str
//    survives the round trip (hicx/py_text.hpp works on that encoding);
//  * numpy arrays are read in C order whatever their strides, and written
//    back in place, as numpy's own in-place operations do; a read-only array
//    raises numpy's "assignment destination is read-only";
//  * the arithmetic kernels take the dtypes int32, int64, float32 and
//    float64 (dispatch_arithmetic); the element-wise NaN and infinity
//    replacements every real and complex dtype (dispatch_numeric). Any other
//    dtype raises TypeError naming it;
//  * a sparse matrix is read through its format attribute ('csr', 'csc',
//    'coo'), shape, data and index arrays, so scipy's *_matrix and *_array
//    classes are both accepted.

#ifndef HICX_PYTHON_API_CONVERT_HPP
#define HICX_PYTHON_API_CONVERT_HPP

#include <complex>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "errors.hpp"
#include "hicx/utilities.hpp"

namespace hicx::python {

namespace py = pybind11;

// --------------------------------------------------------------------------
// str

[[nodiscard]] bool is_str(py::handle object) noexcept;
// The UTF-8 (surrogatepass) bytes of a str.
[[nodiscard]] std::string utf8(py::handle str);
// A str from UTF-8 (surrogatepass) bytes.
[[nodiscard]] py::str to_str(std::string_view utf8_text);
// str(object), repr(object) as UTF-8.
[[nodiscard]] std::string str_of(py::handle object);
[[nodiscard]] std::string repr_of(py::handle object);
// type(object).__name__ as Python's messages use it (tp_name).
[[nodiscard]] std::string type_name(py::handle object);

// --------------------------------------------------------------------------
// numpy

enum class DType {
    Bool, Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64,
    Float16, Float32, Float64, LongDouble, Complex64, Complex128, Other
};

[[nodiscard]] bool is_ndarray(py::handle object) noexcept;
// Raises TypeError unless `object` is a numpy.ndarray.
[[nodiscard]] py::array as_ndarray(py::handle object);
// The element type; byte swapped arrays are Other.
[[nodiscard]] DType dtype_of(const py::array& array);
// str(array.dtype), e.g. "float32".
[[nodiscard]] std::string dtype_name(const py::array& array);

[[noreturn]] void raise_unsupported_dtype(const py::array& array, std::string_view function);

// Calls `f.template operator()<T>()` with T the C++ type of `dtype`, for the
// four dtypes the arithmetic kernels take. Raises TypeError for any other.
template <class F>
decltype(auto) dispatch_arithmetic(const py::array& array, std::string_view function, F&& f) {
    switch (dtype_of(array)) {
        case DType::Int32: return f.template operator()<std::int32_t>();
        case DType::Int64: return f.template operator()<std::int64_t>();
        case DType::Float32: return f.template operator()<float>();
        case DType::Float64: return f.template operator()<double>();
        default: raise_unsupported_dtype(array, function);
    }
}

// The same over every boolean, integer, real (but float16) and complex dtype;
// raises TypeError for any other.
template <class F>
decltype(auto) dispatch_numeric(const py::array& array, std::string_view function, F&& f) {
    switch (dtype_of(array)) {
        case DType::Bool: return f.template operator()<bool>();
        case DType::Int8: return f.template operator()<std::int8_t>();
        case DType::Int16: return f.template operator()<std::int16_t>();
        case DType::Int32: return f.template operator()<std::int32_t>();
        case DType::Int64: return f.template operator()<std::int64_t>();
        case DType::UInt8: return f.template operator()<std::uint8_t>();
        case DType::UInt16: return f.template operator()<std::uint16_t>();
        case DType::UInt32: return f.template operator()<std::uint32_t>();
        case DType::UInt64: return f.template operator()<std::uint64_t>();
        case DType::Float32: return f.template operator()<float>();
        case DType::Float64: return f.template operator()<double>();
        case DType::LongDouble: return f.template operator()<long double>();
        case DType::Complex64: return f.template operator()<std::complex<float>>();
        case DType::Complex128: return f.template operator()<std::complex<double>>();
        default: raise_unsupported_dtype(array, function);
    }
}

// Raises numpy's ValueError for an assignment into a read-only array.
void require_writeable(const py::array& array);

namespace detail {
// Visits every element address of `array` in C order.
template <class G>
void for_each_element(const py::array& array, G&& visit) {
    const py::ssize_t ndim = array.ndim();
    const char* base = static_cast<const char*>(array.data());
    if (ndim == 0) {
        visit(const_cast<char*>(base));
        return;
    }
    py::ssize_t total = 1;
    for (py::ssize_t d = 0; d < ndim; ++d) {
        total *= array.shape(d);
    }
    if (total == 0) {
        return;
    }
    std::vector<py::ssize_t> index(static_cast<std::size_t>(ndim), 0);
    const char* pointer = base;
    for (py::ssize_t n = 0; n < total; ++n) {
        visit(const_cast<char*>(pointer));
        for (py::ssize_t d = ndim - 1; d >= 0; --d) {
            if (++index[static_cast<std::size_t>(d)] < array.shape(d)) {
                pointer += array.strides(d);
                break;
            }
            pointer -= array.strides(d) * (array.shape(d) - 1);
            index[static_cast<std::size_t>(d)] = 0;
        }
    }
}
}  // namespace detail

// The elements in C order, as T (the array's own element type).
template <class T>
[[nodiscard]] std::vector<T> to_vector(const py::array& array) {
    std::vector<T> out;
    out.reserve(static_cast<std::size_t>(array.size()));
    detail::for_each_element(array, [&](char* p) {
        T value;
        std::memcpy(&value, p, sizeof(T));
        out.push_back(value);
    });
    return out;
}

// Writes `values` (C order, array.size() of them) into the array in place.
template <class T>
void write_back(py::array& array, const std::vector<T>& values) {
    require_writeable(array);
    std::size_t k = 0;
    detail::for_each_element(array, [&](char* p) {
        const T value = values[k++];
        std::memcpy(p, &value, sizeof(T));
    });
}

// A new contiguous one dimensional array.
template <class T>
[[nodiscard]] py::array from_vector(const std::vector<T>& values) {
    py::array_t<T> out(static_cast<py::ssize_t>(values.size()));
    if (!values.empty()) {
        std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(T));
    }
    return out;
}

// source.ravel()[keep]: a new one dimensional array of the source's dtype
// holding the elements at the C order positions `keep`.
[[nodiscard]] py::array take_elements(const py::array& source, const std::vector<std::size_t>& keep);

// numpy.T(value): np.float32, np.float64, np.int64 and the other scalar types.
template <class T>
[[nodiscard]] py::object numpy_scalar(T value) {
    return py::dtype::of<T>().attr("type")(value);
}

// --------------------------------------------------------------------------
// scipy.sparse

struct Sparse {
    py::object matrix;
    char format = 'r';  // 'r' csr, 'c' csc, 'o' coo
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    py::array data;
    // csr and csc
    py::array indptr;
    py::array indices;
    // coo
    py::array row;
    py::array col;

    [[nodiscard]] bool compressed() const noexcept { return format != 'o'; }
    // The stored entries in tocoo() order.
    [[nodiscard]] hicx::utilities::SparseEntries entries() const;
    [[nodiscard]] std::vector<std::int64_t> indptr_vector() const;
};

// Reads `object` as a scipy sparse matrix; raises TypeError for anything else.
// `formats` lists the accepted formats ("rco" all three).
[[nodiscard]] Sparse sparse(py::handle object, std::string_view function,
                            std::string_view formats = "rco");

// Replaces matrix.data with `data` (a new array, as `m.data = ...` does).
void set_data(Sparse& matrix, py::array data);

// eliminate_zeros(): `indptr` rewritten, indices and data reduced to the
// `keep` positions.
void apply_eliminate_zeros(Sparse& matrix, const std::vector<std::int64_t>& indptr,
                           const std::vector<std::size_t>& keep);

// --------------------------------------------------------------------------
// Python numbers

// operator.index(object) as int64 (Python int, numpy integers, bool).
[[nodiscard]] bool is_integer(py::handle object);
[[nodiscard]] std::int64_t to_int64(py::handle object);
// float(object)
[[nodiscard]] double to_double(py::handle object);
// int(text) with base 10, Python's own parser.
[[nodiscard]] py::object parse_int(py::handle text);

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_CONVERT_HPP
