// hicexplorer._hicx.utilities: the C++ bodies of hicexplorer.utilities.
//
// The computations are hicx::utilities (cpp/core/include/hicx/utilities.hpp);
// this file converts the Python arguments, calls them and builds the Python
// results of the 3.7.6 functions, of the same types, with the same errors.

#include "utilities.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <complex>
#include <limits>
#include <optional>
#include <set>

#include "api.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/py_text.hpp"
#include "hicx/utilities.hpp"

namespace hicx::python {

namespace util = hicx::utilities;

namespace {

constexpr const char* kLogger = "hicexplorer.utilities";

bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool equal(py::handle a, py::handle b) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_EQ);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

std::size_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

// ---------------------------------------------------------------------------
// element predicates and values for every numeric dtype

template <class T>
struct Element {
    static bool is_nan(T value) {
        if constexpr (std::is_floating_point_v<T>) {
            return std::isnan(value);
        } else if constexpr (std::is_same_v<T, std::complex<float>> ||
                             std::is_same_v<T, std::complex<double>>) {
            return std::isnan(value.real()) || std::isnan(value.imag());
        } else {
            return false;
        }
    }
    static bool is_inf(T value) {
        if constexpr (std::is_floating_point_v<T>) {
            return std::isinf(value);
        } else if constexpr (std::is_same_v<T, std::complex<float>> ||
                             std::is_same_v<T, std::complex<double>>) {
            return std::isinf(value.real()) || std::isinf(value.imag());
        } else {
            return false;
        }
    }
    static T from_double(double value) {
        if constexpr (std::is_same_v<T, std::complex<float>>) {
            return T(static_cast<float>(value), 0.0F);
        } else if constexpr (std::is_same_v<T, std::complex<double>>) {
            return T(value, 0.0);
        } else {
            return static_cast<T>(value);
        }
    }
};

// numpy's double to half conversion for the three values the functions
// assign (0.0, 1.0, 0.000001), round to nearest even.
std::uint16_t half_bits(double value) {
    if (value == 0.0) {
        return 0x0000U;
    }
    if (value == 1.0) {
        return 0x3C00U;
    }
    // 1e-6 is subnormal in half precision: units of 2**-24
    return static_cast<std::uint16_t>(std::nearbyint(value / std::ldexp(1.0, -24)));
}

enum class Which { Nan, Inf };

[[noreturn]] void raise_ufunc_type_error(Which which) {
    raise_error("TypeError", std::string("ufunc '") + (which == Which::Nan ? "isnan" : "isinf") +
                                 "' not supported for the input types, and the inputs could not "
                                 "be safely coerced to any supported types according to the "
                                 "casting rule ''safe''");
}

// The C order positions of the elements of `array` that are NaN (Inf).
std::vector<std::size_t> positions_of(const py::array& array, Which which) {
    std::vector<std::size_t> positions;
    if (dtype_of(array) == DType::Float16) {
        const std::vector<std::uint16_t> bits = to_vector<std::uint16_t>(array);
        for (std::size_t k = 0; k < bits.size(); ++k) {
            const bool exponent = (bits[k] & 0x7C00U) == 0x7C00U;
            const bool mantissa = (bits[k] & 0x03FFU) != 0;
            if (exponent && (which == Which::Nan ? mantissa : !mantissa)) {
                positions.push_back(k);
            }
        }
        return positions;
    }
    const char kind = array.dtype().kind();
    if (kind != 'b' && kind != 'i' && kind != 'u' && kind != 'f' && kind != 'c') {
        raise_ufunc_type_error(which);
    }
    dispatch_numeric(array, which == Which::Nan ? "isnan" : "isinf", [&]<class T>() {
        const std::vector<T> values = to_vector<T>(array);
        for (std::size_t k = 0; k < values.size(); ++k) {
            if (which == Which::Nan ? Element<T>::is_nan(values[k]) : Element<T>::is_inf(values[k])) {
                positions.push_back(k);
            }
        }
    });
    return positions;
}

// array[np.flatnonzero(pred(array))] = value: the flat positions index the
// first axis, as numpy's fancy indexing does with a one dimensional index.
void assign_rows(py::array& array, const std::vector<std::size_t>& positions, double value) {
    if (positions.empty()) {
        return;
    }
    const std::int64_t rows = array.shape(0);
    for (const std::size_t p : positions) {
        if (static_cast<std::int64_t>(p) >= rows) {
            raise_error("IndexError", util::index_out_of_bounds_message(static_cast<std::int64_t>(p), rows));
        }
    }
    require_writeable(array);
    const std::size_t row_size = rows == 0 ? 0 : static_cast<std::size_t>(array.size()) / static_cast<std::size_t>(rows);
    const std::set<std::size_t> selected(positions.begin(), positions.end());
    auto in_selected_row = [&](std::size_t flat) { return selected.count(flat / row_size) != 0; };
    if (dtype_of(array) == DType::Float16) {
        std::vector<std::uint16_t> bits = to_vector<std::uint16_t>(array);
        for (std::size_t k = 0; k < bits.size(); ++k) {
            if (in_selected_row(k)) {
                bits[k] = half_bits(value);
            }
        }
        write_back(array, bits);
        return;
    }
    dispatch_numeric(array, "assignment", [&]<class T>() {
        std::vector<T> values = to_vector<T>(array);
        for (std::size_t k = 0; k < values.size(); ++k) {
            if (in_selected_row(k)) {
                values[k] = Element<T>::from_double(value);
            }
        }
        write_back(array, values);
    });
}

// The same for an object that is not an ndarray: numpy finds the positions
// on np.asarray(object), and the object's own __setitem__ receives the index
// array, which for a list or a memoryview raises.
void assign_rows_foreign(py::handle object, const std::vector<std::size_t>& positions, double value) {
    if (positions.empty()) {
        return;
    }
    py::array_t<std::int64_t> index(static_cast<py::ssize_t>(positions.size()));
    std::copy(positions.begin(), positions.end(), index.mutable_data());
    const py::float_ number(value);
    if (PyObject_SetItem(object.ptr(), index.ptr(), number.ptr()) != 0) {
        throw_python_error();
    }
}

// The convert* functions: `if len(a) == 0: return a`, then replace NaN (and,
// for the second pass, infinities) by `value`.
py::object replace_non_finite(py::handle target, std::initializer_list<std::pair<Which, double>> passes) {
    if (length(target) == 0) {
        return py::reinterpret_borrow<py::object>(target);
    }
    if (is_ndarray(target)) {
        py::array array = py::reinterpret_borrow<py::array>(target);
        for (const auto& [which, value] : passes) {
            assign_rows(array, positions_of(array, which), value);
        }
    } else {
        for (const auto& [which, value] : passes) {
            const py::array view = py::array::ensure(target);
            if (!view) {
                throw_python_error();
            }
            assign_rows_foreign(target, positions_of(view, which), value);
        }
    }
    return py::reinterpret_borrow<py::object>(target);
}

// ---------------------------------------------------------------------------
// statistics over numeric arrays

template <class F>
decltype(auto) with_statistic_values(const py::array& array, std::string_view function, F&& f) {
    switch (dtype_of(array)) {
        case DType::Float32: return f(to_vector<float>(array), true);
        case DType::Float64: return f(to_vector<double>(array), false);
        case DType::UInt64: {
            const std::vector<std::uint64_t> v = to_vector<std::uint64_t>(array);
            return f(std::vector<double>(v.begin(), v.end()), false);
        }
        case DType::Bool: case DType::Int8: case DType::Int16: case DType::Int32: case DType::Int64:
        case DType::UInt8: case DType::UInt16: case DType::UInt32: {
            std::vector<std::int64_t> v;
            dispatch_numeric(array, function, [&]<class T>() {
                if constexpr (std::is_integral_v<T>) {
                    const std::vector<T> raw = to_vector<T>(array);
                    v.assign(raw.begin(), raw.end());
                }
            });
            return f(std::move(v), false);
        }
        default:
            raise_unsupported_dtype(array, function);
    }
}

py::object float_result(double value, bool single) {
    if (single) {
        return numpy_scalar<float>(static_cast<float>(value));
    }
    return numpy_scalar<double>(value);
}

// ---------------------------------------------------------------------------
// numpy array to str (ndarray.astype(str))

std::size_t numeric_str_width(DType dtype) {
    switch (dtype) {
        case DType::Bool: return 5;
        case DType::Int8: return 4;
        case DType::Int16: return 6;
        case DType::Int32: return 11;
        case DType::Int64: return 21;
        case DType::UInt8: return 3;
        case DType::UInt16: return 5;
        case DType::UInt32: return 10;
        case DType::UInt64: return 20;
        case DType::Float32: case DType::Float64: return 32;
        default: return 0;
    }
}

py::array unicode_array(const py::array& like, std::size_t width,
                        const std::vector<std::u32string>& values) {
    std::vector<py::ssize_t> shape(like.shape(), like.shape() + like.ndim());
    py::array out(py::dtype("<U" + std::to_string(std::max<std::size_t>(width, 1))), shape);
    auto* target = static_cast<char32_t*>(out.mutable_data());
    const std::size_t w = std::max<std::size_t>(width, 1);
    std::fill(target, target + values.size() * w, U'\0');
    for (std::size_t k = 0; k < values.size(); ++k) {
        std::copy(values[k].begin(), values[k].end(), target + k * w);
    }
    return out;
}

std::u32string ascii_u32(const std::string& text) {
    return std::u32string(text.begin(), text.end());
}

py::object array_as_str(const py::array& array) {
    const char kind = array.dtype().kind();
    const std::size_t itemsize = static_cast<std::size_t>(array.itemsize());
    std::vector<std::u32string> values;
    values.reserve(static_cast<std::size_t>(array.size()));
    if (kind == 'U') {
        const std::size_t width = itemsize / 4;
        detail::for_each_element(array, [&](char* p) {
            std::u32string value(width, U'\0');
            std::memcpy(value.data(), p, itemsize);
            values.push_back(std::move(value));
        });
        return unicode_array(array, width, values);
    }
    if (kind == 'S') {
        detail::for_each_element(array, [&](char* p) {
            std::size_t used = itemsize;
            while (used > 0 && p[used - 1] == '\0') {
                --used;
            }
            for (std::size_t k = 0; k < used; ++k) {
                if (static_cast<unsigned char>(p[k]) >= 0x80U) {
                    PyObject* decoded = PyUnicode_DecodeASCII(p, static_cast<Py_ssize_t>(used), "strict");
                    if (decoded == nullptr) {
                        throw_python_error();
                    }
                    Py_DECREF(decoded);
                }
            }
            values.push_back(std::u32string(reinterpret_cast<unsigned char*>(p),
                                            reinterpret_cast<unsigned char*>(p) + itemsize));
        });
        return unicode_array(array, itemsize, values);
    }
    if (kind == 'O') {
        std::size_t width = 1;
        detail::for_each_element(array, [&](char* p) {
            PyObject* item = *reinterpret_cast<PyObject**>(p);
            PyObject* text = PyObject_Str(item);
            if (text == nullptr) {
                throw_python_error();
            }
            Py_UCS4* buffer = PyUnicode_AsUCS4Copy(text);
            const Py_ssize_t n = PyUnicode_GetLength(text);
            Py_DECREF(text);
            if (buffer == nullptr) {
                throw_python_error();
            }
            values.emplace_back(reinterpret_cast<char32_t*>(buffer), reinterpret_cast<char32_t*>(buffer) + n);
            PyMem_Free(buffer);
            width = std::max(width, static_cast<std::size_t>(n));
        });
        return unicode_array(array, width, values);
    }
    const DType dtype = dtype_of(array);
    const std::size_t width = numeric_str_width(dtype);
    if (width == 0) {
        raise_unsupported_dtype(array, "toString");
    }
    dispatch_numeric(array, "toString", [&]<class T>() {
        if constexpr (std::is_same_v<T, bool>) {
            for (const bool v : to_vector<bool>(array)) {
                values.push_back(ascii_u32(v ? "True" : "False"));
            }
        } else if constexpr (std::is_integral_v<T>) {
            for (const T v : to_vector<T>(array)) {
                values.push_back(ascii_u32(std::to_string(v)));
            }
        } else if constexpr (std::is_same_v<T, float>) {
            for (const float v : to_vector<float>(array)) {
                values.push_back(ascii_u32(hicx::npy::float32_repr(v)));
            }
        } else if constexpr (std::is_same_v<T, double>) {
            for (const double v : to_vector<double>(array)) {
                values.push_back(ascii_u32(hicx::npy::float_repr(v)));
            }
        }
    });
    return unicode_array(array, width, values);
}

// ---------------------------------------------------------------------------
// sparse helpers

struct NonzeroData {
    util::SparseEntries entries;
    util::Nonzero nonzero;
    std::vector<std::int64_t> distance;
};

template <class T>
NonzeroData nonzero_of(const Sparse& matrix, const std::vector<T>& data) {
    NonzeroData out;
    out.entries = matrix.entries();
    out.nonzero = util::nonzero(out.entries, data.data());
    out.distance = out.nonzero.distances();
    return out;
}

// np.int32(value) of a Python number, as expected_interactions_in_distance
// applies it to the chromosome length and count.
std::int32_t to_int32(py::handle value) {
    if (is_integer(value)) {
        return static_cast<std::int32_t>(to_int64(value));
    }
    const double number = to_double(value);
    if (std::isnan(number)) {
        raise_error("ValueError", "cannot convert float NaN to integer");
    }
    if (std::isinf(number)) {
        raise_error("OverflowError", "cannot convert float infinity to integer");
    }
    return static_cast<std::int32_t>(static_cast<std::int64_t>(number));
}

// The IndexError the accumulation loop of expected_interactions_in_distance
// raises for a distance outside the matrix, before anything else happens.
void check_distances(const std::vector<std::int64_t>& distance, std::int64_t n) {
    for (const std::int64_t d : distance) {
        if (d >= n) {
            raise_error("IndexError", util::index_out_of_bounds_message(d, n));
        }
    }
}

// `a, b, c, d = item`
std::array<py::object, 4> unpack4(py::handle item) {
    PyObject* iterator = PyObject_GetIter(item.ptr());
    if (iterator == nullptr) {
        PyErr_Clear();
        raise_error("TypeError", "cannot unpack non-iterable " + type_name(item) + " object");
    }
    const py::object guard = py::reinterpret_steal<py::object>(iterator);
    std::array<py::object, 4> out;
    for (std::size_t k = 0; k < 5; ++k) {
        PyObject* value = PyIter_Next(iterator);
        if (value == nullptr) {
            if (PyErr_Occurred()) {
                throw_python_error();
            }
            if (k < 4) {
                raise_error("ValueError", "not enough values to unpack (expected 4, got " + std::to_string(k) + ")");
            }
            return out;
        }
        if (k == 4) {
            Py_DECREF(value);
            raise_error("ValueError", "too many values to unpack (expected 4)");
        }
        out[k] = py::reinterpret_steal<py::object>(value);
    }
    return out;
}

std::string thread_failure(const std::string& type, const std::string& message) {
    // traceback.format_exc() of the Python worker process; the C++ worker has
    // no Python frames to show.
    return "Fail: " + message + "Traceback (most recent call last):\n"
           "  File \"hicexplorer._hicx.utilities\", in expected_interactions_thread\n" +
           type + ": " + message + "\n";
}

// expected_interactions(pSubmatrix, pThreads): an array, None, or the
// failure message of a worker.
py::object expected_interactions_object(py::handle submatrix, py::handle threads) {
    Sparse matrix = sparse(submatrix, "expected_interactions");
    return dispatch_arithmetic(matrix.data, "expected_interactions", [&]<class T>() -> py::object {
        const std::vector<T> data = to_vector<T>(matrix.data);
        const NonzeroData nz = nonzero_of(matrix, data);
        const std::int64_t n = matrix.rows;
        if (nz.distance.empty()) {
            return py::none();
        }
        const bool parallel = !threads.is_none() && truthy(threads);
        if (!parallel) {
            return from_vector(*util::expected_interactions(n, nz.distance, data.data(), data.size()));
        }
        if (!is_integer(threads)) {
            raise_error("TypeError", "can't multiply sequence by non-int of type '" + type_name(threads) + "'");
        }
        const std::int64_t workers = to_int64(threads);
        std::vector<double> expected(static_cast<std::size_t>(n), 0.0);
        if (workers > 0) {
            const auto [lo, hi] = std::minmax_element(nz.distance.begin(), nz.distance.end());
            const std::int64_t per_worker = (*hi - *lo) / workers;
            std::optional<std::string> failure;
            for (std::int64_t i = 0; i < workers; ++i) {
                const std::int64_t low = *lo + i * per_worker;
                const std::int64_t high = i < workers - 1 ? *lo + (i + 1) * per_worker : *hi + 1;
                try {
                    const std::vector<double> part =
                        util::distance_sums(n, nz.distance, data.data(), data.size(), low, high);
                    for (std::size_t k = 0; k < part.size(); ++k) {
                        expected[k] += part[k];
                    }
                } catch (const hicx::PythonError& error) {
                    failure = thread_failure(error.type(), error.what());
                }
            }
            if (failure) {
                return py::str(to_str(*failure));
            }
        }
        for (std::int64_t i = 0; i < n; ++i) {
            const double v = expected[static_cast<std::size_t>(i)] / static_cast<double>(n + 1 - i);
            expected[static_cast<std::size_t>(i)] = std::isfinite(v) ? v : 0.0;
        }
        return from_vector(expected);
    });
}

py::object deepcopy(py::handle object) {
    return py::module_::import("copy").attr("deepcopy")(object);
}

}  // namespace

// ---------------------------------------------------------------------------
// toString and toBytes

py::object to_string(py::handle s) {
    if (PyUnicode_Check(s.ptr())) {
        return py::reinterpret_borrow<py::object>(s);
    }
    if (PyBytes_Check(s.ptr())) {
        PyObject* decoded = PyUnicode_DecodeASCII(PyBytes_AS_STRING(s.ptr()), PyBytes_GET_SIZE(s.ptr()), "strict");
        if (decoded == nullptr) {
            throw_python_error();
        }
        return py::reinterpret_steal<py::object>(decoded);
    }
    if (PyList_Check(s.ptr())) {
        py::list out;
        for (const py::handle item : py::reinterpret_borrow<py::list>(s)) {
            out.append(to_string(item));
        }
        return out;
    }
    if (is_ndarray(s)) {
        return array_as_str(py::reinterpret_borrow<py::array>(s));
    }
    return py::reinterpret_borrow<py::object>(s);
}

py::object to_bytes(py::handle s) {
    if (PyBytes_Check(s.ptr())) {
        return py::reinterpret_borrow<py::object>(s);
    }
    if (PyUnicode_Check(s.ptr())) {
        PyObject* encoded = PyUnicode_AsASCIIString(s.ptr());
        if (encoded == nullptr) {
            throw_python_error();
        }
        return py::reinterpret_steal<py::object>(encoded);
    }
    if (PyList_Check(s.ptr())) {
        py::list out;
        for (const py::handle item : py::reinterpret_borrow<py::list>(s)) {
            out.append(to_bytes(item));
        }
        return out;
    }
    return py::reinterpret_borrow<py::object>(s);
}

void bind_utilities(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.utilities.";

    m.def("readBed", [](py::handle pBedFile) {
        const std::vector<std::vector<std::string>> viewpoints = util::read_bed(read_text_lines(pBedFile));
        py::list out;
        for (const auto& v : viewpoints) {
            out.append(py::make_tuple(to_str(v[0]), to_str(v[1]), to_str(v[2])));
        }
        return out;
    }, py::arg("pBedFile"));

    m.def("readTargetBed", [](py::handle pBedFile) {
        util::TargetBed bed;
        try {
            bed = util::read_target_bed(read_text_lines(pBedFile));
        } catch (const hicx::PythonError& error) {
#if PY_VERSION_HEX < 0x030B0000
            if (error.type() == "UnboundLocalError") {
                raise_error("UnboundLocalError", "local variable 'gene' referenced before assignment");
            }
#endif
            throw;
        }
        const py::str outer("c_adj_norm");
        const py::str inner("t_adj_norm");
        py::list genes;
        py::dict target_dict;
        py::dict target_pos_dict;
        for (std::size_t k = 0; k < bed.genes.size(); ++k) {
            const py::str gene = to_str(bed.genes[k]);
            genes.append(gene);
            py::list entry;
            entry.append(outer);
            entry.append(inner);
            entry.append(py::str("genes"));
            entry.append(gene);
            target_dict[gene] = entry;
            py::dict positions;
            positions["chromosome"] = to_str(bed.positions[k].chromosome);
            py::list starts;
            for (const std::string& s : bed.positions[k].starts) {
                starts.append(to_str(s));
            }
            py::list ends;
            for (const std::string& e : bed.positions[k].ends) {
                ends.append(to_str(e));
            }
            positions["start_list"] = starts;
            positions["end_list"] = ends;
            target_pos_dict[gene] = positions;
        }
        py::dict inner_dict;
        inner_dict[inner] = genes;
        py::dict present_genes;
        present_genes[outer] = inner_dict;
        return py::make_tuple(present_genes, target_dict, target_pos_dict);
    }, py::arg("pBedFile"));

    m.def("writableFile", [](py::handle string) {
        int error = 0;
        if (PyLong_Check(string.ptr())) {
            // open(fd, 'w').close() closes the descriptor
            const long long fd = to_int64(string);
            if (fd < 0) {
                raise_error("ValueError", "negative file descriptor");
            }
            struct stat info {};
            if (::fstat(static_cast<int>(fd), &info) != 0) {
                error = errno;
            } else {
                ::close(static_cast<int>(fd));
            }
        } else {
            error = util::truncate_for_writing(fs_path(string));
        }
        if (error != 0) {
            const std::string message = str_of(string) + " file can be opened for writing";
            log(kLogger, "debug", message);
            raise_error("argparse.ArgumentTypeError", message);
        }
        return py::reinterpret_borrow<py::object>(string);
    }, py::arg("string"));

    m.def("remove_outliers", [](py::handle data, py::handle max_deviation) -> py::object {
        const double limit = to_double(max_deviation);
        py::array view = is_ndarray(data) ? py::reinterpret_borrow<py::array>(data) : py::array::ensure(data);
        if (!view) {
            throw_python_error();
        }
        const std::optional<std::vector<std::size_t>> keep =
            with_statistic_values(view, "remove_outliers", [&](auto values, bool) {
                return util::remove_outliers_keep(values, limit);
            });
        if (!keep) {
            return py::reinterpret_borrow<py::object>(data);
        }
        if (!is_ndarray(data)) {
            // data[deviation <= max_deviation] on a list
            raise_error("TypeError", "only integer scalar arrays can be converted to a scalar index");
        }
        return take_elements(view, *keep);
    }, py::arg("data"), py::arg("max_deviation") = 3.5);

    m.def("convertNansToZeros", [](py::handle ma) {
        replace_non_finite(ma.attr("data"), {{Which::Nan, 0.0}});
        return py::reinterpret_borrow<py::object>(ma);
    }, py::arg("ma"));

    m.def("convertInfsToZeros", [](py::handle ma) {
        replace_non_finite(ma.attr("data"), {{Which::Inf, 0.0}});
        return py::reinterpret_borrow<py::object>(ma);
    }, py::arg("ma"));

    m.def("convertInfsToZeros_ArrayFloat", [](py::handle pArray, py::handle pToEpsilon) {
        if (length(pArray) == 0) {
            return py::reinterpret_borrow<py::object>(pArray);
        }
        const double value = truthy(pToEpsilon) ? 0.000001 : 0.0;
        return replace_non_finite(pArray, {{Which::Nan, value}, {Which::Inf, value}});
    }, py::arg("pArray"), py::arg("pToEpsilon") = false);

    m.def("convertNansToOnes", [](py::handle pArray) {
        return replace_non_finite(pArray, {{Which::Nan, 1.0}});
    }, py::arg("pArray"));

    m.def("myAverage", [](py::handle valuesArray, py::handle avgType) -> py::object {
        py::array view = is_ndarray(valuesArray) ? py::reinterpret_borrow<py::array>(valuesArray)
                                                 : py::array::ensure(valuesArray);
        if (!view) {
            throw_python_error();
        }
        const std::vector<std::size_t> nan_positions = positions_of(view, Which::Nan);
        if (!is_ndarray(valuesArray)) {
            // valuesArray[np.logical_not(np.isnan(valuesArray))] on a list
            raise_error("TypeError", "only integer scalar arrays can be converted to a scalar index");
        }
        const bool use_mean = equal(avgType, py::str("mean"));
        const std::set<std::size_t> nans(nan_positions.begin(), nan_positions.end());
        return with_statistic_values(view, "myAverage", [&](auto values, bool single) {
            using V = typename decltype(values)::value_type;
            std::vector<V> kept;
            kept.reserve(values.size());
            for (std::size_t k = 0; k < values.size(); ++k) {
                if (nans.count(k) == 0) {
                    kept.push_back(values[k]);
                }
            }
            const double result = use_mean ? util::mean(kept) : util::median(kept);
            return float_result(result, single);
        });
    }, py::arg("valuesArray"), py::arg("avgType") = "mean");

    m.def("enlarge_bins", [](py::handle bin_intervals) {
        const std::size_t n = length(bin_intervals);
        if (n == 0) {
            // bin_intervals[-1] of an empty sequence
            if (PyObject_GetItem(bin_intervals.ptr(), py::int_(-1).ptr()) == nullptr) {
                throw_python_error();
            }
        }
        std::vector<std::array<py::object, 4>> items(n);
        std::vector<std::int64_t> starts(n);
        std::vector<std::int64_t> ends(n);
        for (std::size_t i = 0; i < n; ++i) {
            const py::object item = bin_intervals[py::int_(i)];
            items[i] = unpack4(item);
            if (!is_integer(items[i][1]) || !is_integer(items[i][2])) {
                raise_error("TypeError", "enlarge_bins: HiCExplorer 4 takes integer start and end positions");
            }
            starts[i] = to_int64(items[i][1]);
            ends[i] = to_int64(items[i][2]);
        }
        std::vector<bool> same_chrom(n > 0 ? n - 1 : 0);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            same_chrom[i] = equal(items[i][0], items[i + 1][0]);
        }
        const std::vector<util::EnlargedBin> bins = util::enlarge_bins(starts, ends, same_chrom);
        auto middle = [&](std::int64_t value, std::size_t type_from) -> py::object {
            if (PyLong_Check(items[type_from][1].ptr())) {
                return py::int_(value);
            }
            return numpy_scalar<std::int64_t>(value);  // numpy integer minus int is int64
        };
        auto coordinate = [&](const util::EnlargedBin& bin, std::size_t i, bool start) -> py::object {
            const util::Coordinate kind = start ? bin.start_kind : bin.end_kind;
            const std::int64_t value = start ? bin.start : bin.end;
            switch (kind) {
                case util::Coordinate::Zero: return py::int_(0);
                case util::Coordinate::Middle: return middle(value, start ? i : i + 1);
                default: return items[i][start ? 1 : 2];
            }
        };
        for (std::size_t i = 0; i < n; ++i) {
            if (!bins[i].rebuilt) {
                continue;
            }
            const py::tuple rebuilt = py::make_tuple(items[bins[i].chrom_from][0], coordinate(bins[i], i, true),
                                                     coordinate(bins[i], i, false), items[i][3]);
            if (PyObject_SetItem(bin_intervals.ptr(), py::int_(i).ptr(), rebuilt.ptr()) != 0) {
                throw_python_error();
            }
        }
        return py::reinterpret_borrow<py::object>(bin_intervals);
    }, py::arg("bin_intervals"));

    m.def("genomicRegion", [](py::handle string) -> py::object {
        std::string text;
        if (is_str(string)) {
            text = utf8(string);
        } else {
            // ''.join(string.split()) on something that is not a str
            const py::object joined = py::str("").attr("join")(string.attr("split")());
            text = utf8(joined);
        }
        try {
            const std::optional<std::string> region = util::genomic_region(text);
            if (!region) {
                return py::none();
            }
            return to_str(*region);
        } catch (const hicx::PythonError& error) {
            raise_error(error.type(), str_of(string) + " is not a valid region");
        }
    }, py::arg("string"));

    m.def("getUserRegion", [](py::handle chromSizes, py::handle regionString, py::handle max_chunk_size) {
        if (!is_str(regionString)) {
            regionString.attr("split")(":");
            raise_error("TypeError", "getUserRegion: HiCExplorer 4 takes the region as a str");
        }
        const std::string region_text = utf8(regionString);
        std::vector<py::str> region;
        for (const std::string_view part : hicx::pytext::split(region_text, ":")) {
            region.push_back(to_str(part));
        }
        const py::str chrom = region[0];
        PyObject* sizes_raw = PyObject_CallOneArg(reinterpret_cast<PyObject*>(&PyDict_Type), chromSizes.ptr());
        if (sizes_raw == nullptr) {
            throw_python_error();
        }
        const py::dict sizes = py::reinterpret_steal<py::dict>(sizes_raw);
        PyObject* size_raw = PyDict_GetItemWithError(sizes.ptr(), chrom.ptr());
        if (size_raw == nullptr) {
            if (PyErr_Occurred()) {
                throw_python_error();
            }
            raise_error("NameError", "Unknown chromosome: " + utf8(chrom) + "\nKnown chromosomes are: " +
                                         str_of(py::list(sizes)) + " ");
        }
        const py::object size = py::reinterpret_borrow<py::object>(size_raw);
        std::optional<std::int64_t> start;
        if (region.size() > 1) {
            start = to_int64(parse_int(region[1]));
        }
        std::optional<std::int64_t> end;
        py::object end_object;
        if (region.size() > 2) {
            end_object = parse_int(region[2]);
            end = to_int64(end_object);
        }
        if (!is_integer(size)) {
            raise_error("TypeError", "getUserRegion: HiCExplorer 4 takes integer chromosome sizes");
        }
        util::UserRegion bounds = util::user_region_bounds(to_int64(size), start, end, region_text);
        const py::object untiled_end = bounds.end_is_chrom_size ? size : end_object;
        std::optional<std::int64_t> tile;
        if (region.size() > 3) {
            tile = to_int64(parse_int(region[3]));
        }
        util::ChunkLimit limit;
        if (is_integer(max_chunk_size)) {
            limit.is_float = false;
            limit.integer = to_int64(max_chunk_size);
        } else {
            limit.floating = to_double(max_chunk_size);
        }
        util::user_region_finish(bounds, tile, limit);
        py::object region_end = untiled_end;
        if (bounds.end_tiled) {
            region_end = (bounds.end_is_chrom_size && !PyLong_Check(size.ptr()))
                             ? numpy_scalar<std::int64_t>(bounds.end)
                             : py::object(py::int_(bounds.end));
        }
        py::list chrom_sizes;
        chrom_sizes.append(py::make_tuple(chrom, untiled_end));
        return py::make_tuple(chrom_sizes, py::int_(bounds.start), region_end, py::int_(bounds.chunk_size));
    }, py::arg("chromSizes"), py::arg("regionString"), py::arg("max_chunk_size") = 1e6);

    m.def("expected_interactions_in_distance", [](py::handle pLength_chromosome, py::handle pChromosome_count,
                                                  py::handle pSubmatrix) {
        Sparse matrix = sparse(pSubmatrix, "expected_interactions_in_distance");
        return dispatch_arithmetic(matrix.data, "expected_interactions_in_distance", [&]<class T>() {
            const std::vector<T> data = to_vector<T>(matrix.data);
            const NonzeroData nz = nonzero_of(matrix, data);
            // the loop runs before np.int32() converts the two numbers
            check_distances(nz.distance, matrix.rows);
            const std::int32_t count = to_int32(pChromosome_count);
            const std::int32_t length = to_int32(pLength_chromosome);
            return from_vector(util::expected_interactions_in_distance(matrix.rows, nz.distance, data.data(),
                                                                       length, count));
        });
    }, py::arg("pLength_chromosome"), py::arg("pChromosome_count"), py::arg("pSubmatrix"));

    m.def("expected_interactions_non_zero", [](py::handle pSubmatrix) {
        Sparse matrix = sparse(pSubmatrix, "expected_interactions_non_zero");
        return dispatch_arithmetic(matrix.data, "expected_interactions_non_zero", [&]<class T>() {
            const std::vector<T> data = to_vector<T>(matrix.data);
            const NonzeroData nz = nonzero_of(matrix, data);
            return from_vector(util::expected_interactions_non_zero(matrix.rows, nz.distance, data.data()));
        });
    }, py::arg("pSubmatrix"));

    m.def("expected_interactions_thread", [](py::handle pData, py::handle pDistances, py::handle pMinDistance,
                                             py::handle pMaxDistance, py::handle pSize, py::handle pQueue) {
        py::object result;
        try {
            const py::array data = as_ndarray(pData);
            const py::array distances_array = as_ndarray(pDistances);
            std::vector<std::int64_t> distances;
            dispatch_numeric(distances_array, "expected_interactions_thread", [&]<class D>() {
                if constexpr (std::is_integral_v<D> && !std::is_same_v<D, bool>) {
                    const std::vector<D> raw = to_vector<D>(distances_array);
                    distances.assign(raw.begin(), raw.end());
                } else {
                    raise_unsupported_dtype(distances_array, "expected_interactions_thread");
                }
            });
            const std::int64_t size = to_int64(pSize);
            const std::int64_t low = to_int64(pMinDistance);
            const std::int64_t high = to_int64(pMaxDistance);
            result = dispatch_arithmetic(data, "expected_interactions_thread", [&]<class T>() {
                const std::vector<T> values = to_vector<T>(data);
                return from_vector(util::distance_sums(size, distances, values.data(), values.size(), low, high));
            });
        } catch (const hicx::PythonError& error) {
            result = to_str(thread_failure(error.type(), error.what()));
        } catch (py::error_already_set& error) {
            const std::string type = py::str(error.type().attr("__name__"));
            const std::string message = str_of(error.value());
            result = to_str(thread_failure(type, message));
        }
        pQueue.attr("put")(result);
        return py::none();
    }, py::arg("pData"), py::arg("pDistances"), py::arg("pMinDistance"), py::arg("pMaxDistance"),
       py::arg("pSize"), py::arg("pQueue"));

    m.def("expected_interactions", [](py::handle pSubmatrix, py::handle pThreads) {
        return expected_interactions_object(pSubmatrix, pThreads);
    }, py::arg("pSubmatrix"), py::arg("pThreads") = py::none());

    m.def("compute_zscore", [](py::handle pSubmatrix, py::handle pDepth, py::handle pThreads) {
        const py::object expected_object = expected_interactions_object(pSubmatrix, pThreads);
        Sparse matrix = sparse(pSubmatrix, "compute_zscore");
        dispatch_arithmetic(matrix.data, "compute_zscore", [&]<class T>() {
            std::vector<T> data = to_vector<T>(matrix.data);
            const NonzeroData nz = nonzero_of(matrix, data);
            std::optional<std::vector<double>> expected;
            if (is_str(expected_object)) {
                if (matrix.rows > 0) {
                    raise_error("TypeError", "ufunc 'subtract' did not contain a loop with signature matching "
                                             "types (dtype('float64'), dtype('<U1')) -> None");
                }
            } else if (!expected_object.is_none()) {
                expected = py::reinterpret_borrow<py::array>(expected_object).cast<std::vector<double>>();
            }
            std::optional<double> depth;
            bool depth_is_float = false;
            if (!pDepth.is_none()) {
                depth_is_float = !is_integer(pDepth);
                depth = depth_is_float ? to_double(pDepth) : static_cast<double>(to_int64(pDepth));
            }
            try {
                util::compute_zscore(matrix.rows, nz.distance, data.data(), data.size(), depth, depth_is_float,
                                     expected);
            } catch (...) {
                py::array target = matrix.data;
                write_back(target, data);  // the depth mask was applied before the error
                throw;
            }
            py::array target = matrix.data;
            write_back(target, data);
        });
        return py::reinterpret_borrow<py::object>(pSubmatrix);
    }, py::arg("pSubmatrix"), py::arg("pDepth"), py::arg("pThreads"));

    m.def("obs_exp_matrix_lieberman", [](py::handle pSubmatrix, py::handle pLength_chromosome,
                                         py::handle pChromosome_count) {
        Sparse matrix = sparse(pSubmatrix, "obs_exp_matrix_lieberman");
        dispatch_arithmetic(matrix.data, "obs_exp_matrix_lieberman", [&]<class T>() {
            const std::vector<T> data = to_vector<T>(matrix.data);
            const NonzeroData nz = nonzero_of(matrix, data);
            check_distances(nz.distance, matrix.rows);
            const std::int32_t count = to_int32(pChromosome_count);
            const std::int32_t length = to_int32(pLength_chromosome);
            const std::vector<double> expected =
                util::expected_interactions_in_distance(matrix.rows, nz.distance, data.data(), length, count);
            if (!data.empty()) {
                set_data(matrix, from_vector(util::obs_exp_divide(data.data(), data.size(), expected,
                                                                  nz.distance, 0.0)));
            }
        });
        return py::reinterpret_borrow<py::object>(pSubmatrix);
    }, py::arg("pSubmatrix"), py::arg("pLength_chromosome"), py::arg("pChromosome_count"));

    m.def("obs_exp_matrix_non_zero", [](py::handle pSubmatrix, py::handle ligation_factor, py::handle pInplace,
                                        py::handle pToEpsilon, py::handle /*pThreads*/) {
        const py::object target = truthy(pInplace) ? py::reinterpret_borrow<py::object>(pSubmatrix)
                                                   : deepcopy(pSubmatrix);
        Sparse matrix = sparse(target, "obs_exp_matrix_non_zero", "rc");
        const bool ligation = truthy(ligation_factor);
        const bool epsilon = truthy(pToEpsilon);
        dispatch_arithmetic(matrix.data, "obs_exp_matrix_non_zero", [&]<class T>() {
            const std::vector<T> data = to_vector<T>(matrix.data);
            const NonzeroData nz = nonzero_of(matrix, data);
            const std::vector<double> expected =
                util::expected_interactions_non_zero(matrix.rows, nz.distance, data.data());
            const std::vector<std::int64_t> indptr = matrix.indptr_vector();
            std::vector<float> values;
            try {
                values = util::obs_exp_non_zero_values(nz.entries, indptr, matrix.format == 'r', data.data(),
                                                       data.size(), nz.nonzero, expected, ligation, epsilon);
            } catch (...) {
                // submatrix.data was already the float32 copy when the loop raised
                std::vector<float> converted(data.begin(), data.end());
                set_data(matrix, from_vector(converted));
                throw;
            }
            set_data(matrix, from_vector(values));
            std::vector<std::int64_t> pruned = indptr;
            const std::vector<std::size_t> keep = util::eliminate_zeros(pruned, values.data());
            apply_eliminate_zeros(matrix, pruned, keep);
        });
        return target;
    }, py::arg("pSubmatrix"), py::arg("ligation_factor") = false, py::arg("pInplace") = true,
       py::arg("pToEpsilon") = false, py::arg("pThreads") = py::none());

    m.def("obs_exp_matrix", [](py::handle pSubmatrix, py::handle pInplace, py::handle pToEpsilon,
                               py::handle pThreads, py::handle pDistance) -> py::object {
        if (!pDistance.is_none()) {
            Sparse matrix = sparse(pSubmatrix, "obs_exp_matrix", "rc");
            const double limit = to_double(pDistance);
            dispatch_arithmetic(matrix.data, "obs_exp_matrix", [&]<class T>() {
                std::vector<T> data = to_vector<T>(matrix.data);
                const NonzeroData nz = nonzero_of(matrix, data);
                util::zero_from_distance(nz.distance, data.data(), data.size(), limit);
                write_back(matrix.data, data);
                std::vector<std::int64_t> indptr = matrix.indptr_vector();
                const std::vector<std::size_t> keep = util::eliminate_zeros(indptr, data.data());
                apply_eliminate_zeros(matrix, indptr, keep);
            });
        }
        const py::object expected_object = expected_interactions_object(pSubmatrix, pThreads);
        if (expected_object.is_none()) {
            return py::none();
        }
        if (is_str(expected_object)) {
            raise_error("TypeError", "only integer scalar arrays can be converted to a scalar index");
        }
        const std::vector<double> expected = expected_object.cast<std::vector<double>>();
        const bool inplace = truthy(pInplace);
        const double non_finite = truthy(pToEpsilon) ? 0.000001 : 0.0;
        Sparse source = sparse(pSubmatrix, "obs_exp_matrix");
        const py::object target = inplace ? py::reinterpret_borrow<py::object>(pSubmatrix) : deepcopy(pSubmatrix);
        Sparse matrix = sparse(target, "obs_exp_matrix");
        dispatch_arithmetic(source.data, "obs_exp_matrix", [&]<class T>() {
            const std::vector<T> data = to_vector<T>(source.data);
            if (data.empty()) {
                return;
            }
            const NonzeroData nz = nonzero_of(source, data);
            set_data(matrix, from_vector(util::obs_exp_divide(data.data(), data.size(), expected, nz.distance,
                                                              non_finite)));
        });
        return target;
    }, py::arg("pSubmatrix"), py::arg("pInplace") = true, py::arg("pToEpsilon") = false,
       py::arg("pThreads") = py::none(), py::arg("pDistance") = py::none());

    m.def("toString", [](py::handle s) { return to_string(s); }, py::arg("s"));
    m.def("toBytes", [](py::handle s) { return to_bytes(s); }, py::arg("s"));

    m.def("check_chrom_str_bytes", [](py::handle pIteratableObj, py::handle pObj) {
        py::object type_;
        if (PyList_Check(pObj.ptr()) && PyList_GET_SIZE(pObj.ptr()) > 0) {
            type_ = py::type::of(py::reinterpret_borrow<py::list>(pObj)[0]);
        } else {
            type_ = py::type::of(pObj);
        }
        auto first_type = [&]() -> py::object {
            py::object iterator = py::reinterpret_steal<py::object>(PyObject_GetIter(pIteratableObj.ptr()));
            if (!iterator) {
                throw_python_error();
            }
            PyObject* item = PyIter_Next(iterator.ptr());
            if (item == nullptr) {
                if (!PyErr_Occurred()) {
                    PyErr_SetNone(PyExc_StopIteration);
                }
                throw_python_error();
            }
            return py::type::of(py::reinterpret_steal<py::object>(item));
        };
        const int instance = PyObject_IsInstance(first_type().ptr(), type_.ptr());
        if (instance < 0) {
            throw_python_error();
        }
        if (instance == 0) {
            if (first_type().ptr() == reinterpret_cast<PyObject*>(&PyUnicode_Type)) {
                return to_string(pObj);
            }
            const py::object third = first_type();
            const py::object bytes_ = py::module_::import("numpy").attr("bytes_");
            if (third.ptr() == reinterpret_cast<PyObject*>(&PyBytes_Type) || equal(third, bytes_)) {
                return to_bytes(pObj);
            }
        }
        return py::reinterpret_borrow<py::object>(pObj);
    }, py::arg("pIteratableObj"), py::arg("pObj"));

    m.def("change_chrom_names", [](py::handle chrom) -> py::object {
        const py::object name = to_string(chrom);
        if (is_str(name)) {
            return to_str(util::change_chrom_name(utf8(name)));
        }
        // chrom.startswith('chr') on what toString left unchanged
        if (truthy(name.attr("startswith")("chr"))) {
            return name[py::slice(3, std::nullopt, std::nullopt)];
        }
        return py::str("chr") + name;
    }, py::arg("chrom"));

    m.def("opener", [](py::handle filename) -> py::object {
        const py::object handle = py::module_::import("builtins").attr("open")(filename, "rb");
        const py::bytes first = handle.attr("read")(2);
        handle.attr("seek")(0);
        if (util::gzip_magic(std::string(first))) {
            return py::module_::import("gzip").attr("GzipFile")(py::arg("fileobj") = handle);
        }
        return handle;
    }, py::arg("filename"));

    m.def("remove_non_ascii", [](py::handle pText) -> py::object {
        if (!is_str(pText)) {
            pText.attr("encode")("ASCII");
            raise_error("TypeError", "remove_non_ascii: HiCExplorer 4 takes a str");
        }
        const std::string text = utf8(pText);
        if (util::is_ascii(text)) {
            return py::reinterpret_borrow<py::object>(pText);
        }
        std::vector<char32_t> surrogates;
        const std::string out = util::unidecode(text, &surrogates);
        for (const char32_t cp : surrogates) {
            const py::str character = py::reinterpret_steal<py::str>(
                PyUnicode_FromOrdinal(static_cast<int>(cp)));
            warn(PyExc_RuntimeWarning, "Surrogate character " + repr_of(character) +
                                           " will be ignored. You might be using a narrow Python build.", 2);
        }
        return to_str(out);
    }, py::arg("pText"));

    m.def("check_cooler", [](py::handle pFileName) -> py::object {
        if (!is_str(pFileName)) {
            pFileName.attr("endswith")(".cool");
            raise_error("TypeError", "check_cooler: HiCExplorer 4 takes a str");
        }
        const std::string name = utf8(pFileName);
        if (name.size() >= 5 && name.compare(name.size() - 5, 5, ".cool") == 0) {
            return py::bool_(true);
        }
        if (hicx::pytext::split(name, "::").size() > 2) {
            raise_error("ValueError", "Invalid Cooler URI string");
        }
        return py::bool_(hicx::check_cooler(name));
    }, py::arg("pFileName"));

    m.def("in_units", [](py::handle pBasePosition) {
        return to_str(util::in_units(to_double(pBasePosition)));
    }, py::arg("pBasePosition"));
}

}  // namespace hicx::python
