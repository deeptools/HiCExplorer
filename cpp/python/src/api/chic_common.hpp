// Helpers shared by the bindings of hicexplorer.lib.viewpoint and the cHi-C
// modules (chicViewpoint, chicViewpointBackgroundModel, chicQualityControl,
// chicExportData, chicSignificantInteractions, chicDifferentialTest,
// chicAggregateStatistic, chicPlotViewpoint) and hicPlotViewpoint.
//
// Header only, so that the bindings need no new source in the extension's
// CMake list. What is here is Python's object protocol as the 3.7.6 code used
// it (unpacking, `in`, comparisons, scalar arithmetic that decides the type of
// a result, '{}'.format, exit(), Queue.put) and the h5py view of the files
// (hicx::chic::H5Reader turned into the objects h5py returns). Nothing here
// computes: the computations are hicx::chic and hicx::scipy in cpp/core.

#ifndef HICX_PYTHON_API_CHIC_COMMON_HPP
#define HICX_PYTHON_API_CHIC_COMMON_HPP

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <initializer_list>
#include <span>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

#include "api.hpp"
#include "hicx/chic_hdf5.hpp"

namespace hicx::python::chic_api {

namespace py = pybind11;

// ---------------------------------------------------------------------------
// modules

inline py::module_ numpy() {
    return py::module_::import("numpy");
}

inline py::object builtin(const char* name) {
    return py::module_::import("builtins").attr(name);
}

// hicexplorer._version.__version__, read from the package as the 3.7.6
// modules imported it.
inline py::object package_version() {
    return py::module_::import("hicexplorer._version").attr("__version__");
}

// ---------------------------------------------------------------------------
// object protocol

inline py::object steal(PyObject* object) {
    if (object == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(object);
}

inline bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool rich(py::handle a, py::handle b, int op) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), op);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}
inline bool eq(py::handle a, py::handle b) { return rich(a, b, Py_EQ); }
inline bool ne(py::handle a, py::handle b) { return rich(a, b, Py_NE); }
inline bool lt(py::handle a, py::handle b) { return rich(a, b, Py_LT); }
inline bool le(py::handle a, py::handle b) { return rich(a, b, Py_LE); }
inline bool gt(py::handle a, py::handle b) { return rich(a, b, Py_GT); }
inline bool ge(py::handle a, py::handle b) { return rich(a, b, Py_GE); }

// `item in container`
inline bool contains(py::handle container, py::handle item) {
    const int result = PySequence_Contains(container.ptr(), item.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline py::ssize_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return n;
}

inline py::object getitem(py::handle object, py::handle key) {
    return steal(PyObject_GetItem(object.ptr(), key.ptr()));
}
inline py::object getitem(py::handle object, py::ssize_t index) {
    return getitem(object, py::int_(index));
}
inline void setitem(py::handle object, py::handle key, py::handle value) {
    if (PyObject_SetItem(object.ptr(), key.ptr(), value.ptr()) < 0) {
        throw_python_error();
    }
}
inline void setitem(py::handle object, py::ssize_t index, py::handle value) {
    setitem(object, py::int_(index), value);
}

// Python's scalar operators, for the arithmetic whose result type the
// callers see (int, numpy.int64, numpy.float64 by the operands' types).
inline py::object add(py::handle a, py::handle b) { return steal(PyNumber_Add(a.ptr(), b.ptr())); }
inline py::object sub(py::handle a, py::handle b) { return steal(PyNumber_Subtract(a.ptr(), b.ptr())); }
inline py::object truediv(py::handle a, py::handle b) {
    return steal(PyNumber_TrueDivide(a.ptr(), b.ptr()));
}
inline py::object floordiv(py::handle a, py::handle b) {
    return steal(PyNumber_FloorDivide(a.ptr(), b.ptr()));
}
inline py::object neg(py::handle a) { return steal(PyNumber_Negative(a.ptr())); }

// int(object), float(object)
inline py::object int_of(py::handle object) { return steal(PyNumber_Long(object.ptr())); }
inline py::object float_of(py::handle object) { return steal(PyNumber_Float(object.ptr())); }

// str(object), format(object, spec)
inline py::object str_object(py::handle object) { return steal(PyObject_Str(object.ptr())); }
inline py::object format(py::handle object, const char* spec) {
    return steal(PyObject_Format(object.ptr(), py::str(spec).ptr()));
}
inline py::object format(py::handle object, py::handle spec) {
    return steal(PyObject_Format(object.ptr(), spec.ptr()));
}

// The items of `object` unpacked into `count` targets, with the errors of
// `a, b = object`.
inline std::vector<py::object> unpack(py::handle object, std::size_t count) {
    PyObject* iterator = PyObject_GetIter(object.ptr());
    if (iterator == nullptr) {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) {
            PyErr_Clear();
            raise_error("TypeError", "cannot unpack non-iterable " + type_name(object) + " object");
        }
        throw_python_error();
    }
    const py::object it = py::reinterpret_steal<py::object>(iterator);
    std::vector<py::object> out;
    for (std::size_t i = 0; i < count; ++i) {
        PyObject* item = PyIter_Next(it.ptr());
        if (item == nullptr) {
            if (PyErr_Occurred()) {
                throw_python_error();
            }
            raise_error("ValueError", "not enough values to unpack (expected " + std::to_string(count) +
                                          ", got " + std::to_string(i) + ")");
        }
        out.push_back(py::reinterpret_steal<py::object>(item));
    }
    PyObject* extra = PyIter_Next(it.ptr());
    if (extra != nullptr) {
        Py_DECREF(extra);
        raise_error("ValueError", "too many values to unpack (expected " + std::to_string(count) + ")");
    }
    if (PyErr_Occurred()) {
        throw_python_error();
    }
    return out;
}

// list(object)
inline py::list as_list(py::handle object) {
    return py::reinterpret_borrow<py::list>(steal(PySequence_List(object.ptr())));
}

// '/'.join(parts)
inline py::object join(const char* separator, py::handle parts) {
    return py::str(separator).attr("join")(parts);
}

// exit(code): the site builtin, which closes sys.stdin and raises SystemExit.
[[noreturn]] inline void python_exit(int code) {
    builtin("exit")(code);
    raise_error("SystemExit", std::to_string(code));
}

// str() of the exception a hicx::PythonError stands for, as `str(exp)` in an
// `except Exception as exp:` block sees it (a KeyError quotes its message).
inline std::string exception_text(const hicx::PythonError& error) {
    return str_of(exception_type(error.type())(to_str(error.what())));
}

// The str of the current Python exception from a py::error_already_set.
inline std::string exception_text(const py::error_already_set& error) {
    return str_of(error.value());
}

// `except Exception:` catches the error (not SystemExit, KeyboardInterrupt
// or GeneratorExit).
inline bool is_exception(py::error_already_set& error) {
    return error.matches(PyExc_Exception);
}

// 'Fail: ' + str(exp), the message the worker functions put into their queue
// when `except Exception as exp:` caught something.
inline std::string fail_text(const std::exception& error) {
    if (const auto* python = dynamic_cast<const hicx::PythonError*>(&error)) {
        return exception_text(*python);
    }
    if (const auto* set = dynamic_cast<const py::error_already_set*>(&error)) {
        return str_of(set->value());
    }
    return error.what();
}

// traceback.format_exc() in the `except` block of a worker function. For a
// Python error it is the traceback of the error as Python formats it (its
// frames are the Python ones the error passed, not 3.7.6's); for an error
// raised by the C++ it is the header and the exception line.
inline std::string traceback_text(const std::exception& error) {
    if (const auto* set = dynamic_cast<const py::error_already_set*>(&error)) {
        const py::handle trace = set->trace() ? py::handle(set->trace()) : py::handle(Py_None);
        const py::object lines = py::module_::import("traceback")
                                     .attr("format_exception")(set->type(), set->value(), trace);
        const std::string text = utf8(py::str("").attr("join")(lines));
        if (text.rfind("Traceback", 0) == 0) {
            return text;
        }
        return "Traceback (most recent call last):\n" + text;
    }
    std::string type = "Exception";
    std::string message = error.what();
    if (const auto* python = dynamic_cast<const hicx::PythonError*>(&error)) {
        type = python->type();
        message = exception_text(*python);
    }
    return "Traceback (most recent call last):\n" + type + ": " + message + "\n";
}

// numpy scalars
inline py::object np_int64(std::int64_t value) { return numpy_scalar<std::int64_t>(value); }
inline py::object np_float64(double value) { return numpy_scalar<double>(value); }

inline bool is_numpy_scalar(py::handle object) {
    return PyObject_IsInstance(object.ptr(), numpy().attr("generic").ptr()) == 1;
}

// The keys of a dict as `sorted(d)` orders them. Numbers are ordered by value
// in C++ (sorted is stable and dict keys are distinct, so the order is the
// same); anything else falls back to the builtin.
inline py::list sorted_keys(py::handle mapping) {
    py::list keys = as_list(mapping);
    const std::size_t n = keys.size();
    bool all_integers = true;
    bool all_numbers = true;
    std::vector<std::int64_t> integers;
    std::vector<double> reals;
    for (std::size_t i = 0; i < n && all_numbers; ++i) {
        const py::object key = keys[i];
        if (PyBool_Check(key.ptr()) || !(is_integer(key) || PyFloat_Check(key.ptr()) ||
                                         PyObject_IsInstance(key.ptr(), numpy().attr("floating").ptr()) == 1)) {
            all_numbers = false;
            break;
        }
        if (all_integers && is_integer(key)) {
            integers.push_back(to_int64(key));
        } else {
            all_integers = false;
        }
        reals.push_back(to_double(key));
    }
    if (!all_numbers) {
        return builtin("sorted")(mapping);
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (all_integers) {
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) { return integers[a] < integers[b]; });
    } else {
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) { return reals[a] < reals[b]; });
    }
    py::list out;
    for (const std::size_t i : order) {
        out.append(keys[i]);
    }
    return out;
}

// min(mapping) / max(mapping) of numeric keys
inline py::object extreme_key(py::handle mapping, bool maximum) {
    const py::list keys = as_list(mapping);
    if (keys.empty()) {
        raise_error("ValueError", std::string(maximum ? "max" : "min") + "() iterable argument is empty");
    }
    py::object best = keys[0];
    for (std::size_t i = 1; i < keys.size(); ++i) {
        const py::object key = keys[i];
        if (maximum ? gt(key, best) : lt(key, best)) {
            best = key;
        }
    }
    return best;
}

// file.write(text)
inline void write(py::handle file, const std::string& text) {
    file.attr("write")(to_str(text));
}

// open(path, mode), the Python file object a with block writes through.
inline py::object open_file(py::handle path, const char* mode) {
    return builtin("open")(path, mode);
}

// Runs `body(file)` inside `with open(path, mode) as file:`.
template <class F>
void with_open(py::handle path, const char* mode, F&& body) {
    py::object file = open_file(path, mode);
    try {
        body(file);
    } catch (...) {
        try {
            file.attr("close")();
        } catch (...) {
        }
        throw;
    }
    file.attr("close")();
}

// ---------------------------------------------------------------------------
// h5py

// h5py.File(path, 'r')
inline hicx::chic::H5Reader open_h5(py::handle path) {
    return hicx::chic::H5Reader(fs_path(path));
}

// The numpy dtype object of a value h5py reads.
inline py::object h5_dtype(const hicx::chic::H5Value& value) {
    if (value.kind == hicx::chic::H5Value::Kind::Bytes && value.variable_length) {
        return numpy().attr("dtype")("O");
    }
    return numpy().attr("dtype")(value.dtype_str());
}

// dataset[:] (as np.array(...) of it): an array of the stored dtype, bytes
// objects for variable length strings.
inline py::object h5_array(const hicx::chic::H5Value& value) {
    using Kind = hicx::chic::H5Value::Kind;
    if (value.kind == Kind::Bytes) {
        if (value.variable_length) {
            py::array out(py::dtype("O"), std::vector<py::ssize_t>{static_cast<py::ssize_t>(value.strings.size())});
            for (std::size_t i = 0; i < value.strings.size(); ++i) {
                py::bytes item(value.strings[i]);
                PyObject** slot = static_cast<PyObject**>(out.mutable_data()) + i;
                Py_XDECREF(*slot);
                *slot = item.release().ptr();
            }
            return out;
        }
        py::list items;
        for (const std::string& s : value.strings) {
            items.append(py::bytes(s));
        }
        return numpy().attr("array")(items, h5_dtype(value));
    }
    const py::dtype dtype = py::reinterpret_borrow<py::dtype>(h5_dtype(value));
    py::array out(dtype, std::vector<py::ssize_t>{static_cast<py::ssize_t>(value.length)});
    if (!value.raw.empty()) {
        std::memcpy(out.mutable_data(), value.raw.data(), value.raw.size());
    }
    return out;
}

inline py::object h5_element(const hicx::chic::H5Value& value, std::size_t i);

// dataset[()] of a scalar dataset: a numpy scalar, bytes for a variable length
// string, numpy.bytes_ for a fixed length one.
inline py::object h5_scalar(const hicx::chic::H5Value& value) {
    using Kind = hicx::chic::H5Value::Kind;
    if (!value.scalar) {
        return h5_array(value);
    }
    if (value.kind == Kind::Bytes) {
        const std::string text = value.strings.empty() ? std::string() : value.strings[0];
        if (value.variable_length) {
            return py::bytes(text);
        }
        return numpy().attr("bytes_")(py::bytes(text));
    }
    return h5_element(value, 0);
}

// The element i of an array h5py read, as indexing the numpy array gives it.
inline py::object h5_element(const hicx::chic::H5Value& value, std::size_t i) {
    using Kind = hicx::chic::H5Value::Kind;
    switch (value.kind) {
        case Kind::Bytes:
            if (value.variable_length) {
                return py::bytes(value.strings[i]);
            }
            return numpy().attr("bytes_")(py::bytes(value.strings[i]));
        case Kind::Float:
            if (value.item_size == 4) {
                return numpy_scalar<float>(static_cast<float>(value.as_double(i)));
            }
            return numpy_scalar<double>(value.as_double(i));
        case Kind::Integer:
            switch (value.item_size) {
                case 1: return numpy_scalar<std::int8_t>(static_cast<std::int8_t>(value.as_int64(i)));
                case 2: return numpy_scalar<std::int16_t>(static_cast<std::int16_t>(value.as_int64(i)));
                case 4: return numpy_scalar<std::int32_t>(static_cast<std::int32_t>(value.as_int64(i)));
                default: return numpy_scalar<std::int64_t>(value.as_int64(i));
            }
        default:
            return getitem(h5_array(value), static_cast<py::ssize_t>(i));
    }
}

// float(element), as the 3.7.6 code converts the values it reads
inline double h5_float(const hicx::chic::H5Value& value, std::size_t i) {
    if (value.kind == hicx::chic::H5Value::Kind::Bytes) {
        return to_double(float_of(py::bytes(value.strings[i])));
    }
    return value.as_double(i);
}

// x.decode("utf-8") of what dataset[()] returned
inline py::object decode_utf8(const py::object& value) {
    if (!PyBytes_Check(value.ptr())) {
        raise_error("AttributeError", "'" + type_name(value) + "' object has no attribute 'decode'");
    }
    return value.attr("decode")("utf-8");
}

// ---------------------------------------------------------------------------
// writing HDF5 files as the cHi-C functions wrote them with h5py

// h5py.File(path, 'w'), with h5py's error when the file cannot be created.
inline hicx::chic::Hdf5Writer create_h5(py::handle path_object) {
    const std::string path = fs_path(path_object);
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) {
        const int error = errno;
        const char* type = error == ENOENT ? "FileNotFoundError"
                           : (error == EACCES || error == EPERM) ? "PermissionError"
                           : error == EISDIR ? "IsADirectoryError" : "OSError";
        raise_error(type, "[Errno " + std::to_string(error) +
                              "] Unable to synchronously create file (unable to open file: name = '" + path +
                              "', errno = " + std::to_string(error) + ", error message = '" + std::strerror(error) +
                              "', flags = 13, o_flags = 242)");
    }
    ::close(fd);
    return hicx::chic::Hdf5Writer(path);
}

// "/" + '/'.join(parts)
inline std::string h5_path(std::initializer_list<std::string> parts) {
    std::string out;
    for (const std::string& part : parts) {
        out += "/" + part;
    }
    return out;
}

// group.create_group(name) of h5py: ValueError when it exists.
inline void create_group(hicx::chic::Hdf5Writer& writer, const std::string& path) {
    if (writer.exists(path)) {
        raise_error("ValueError", "Unable to synchronously create group (name already exists)");
    }
    writer.create_group(path);
}

// `if name not in group: group.create_group(name)`
inline void ensure_group(hicx::chic::Hdf5Writer& writer, const std::string& path) {
    if (!writer.exists(path)) {
        writer.create_group(path);
    }
}

// group[link] = file[target]: OSError when the name exists.
inline void hard_link(hicx::chic::Hdf5Writer& writer, const std::string& target, const std::string& link) {
    if (!writer.hard_link(target, link)) {
        raise_error("OSError", "Unable to create link (name already exists)");
    }
}

// The numpy dtype h5py would store `values` with (np.asarray(values)).
inline py::array as_array(py::handle values) {
    return py::reinterpret_borrow<py::array>(numpy().attr("asarray")(values));
}

[[noreturn]] inline void unsupported_dtype(const py::array& array) {
    const char kind = array.dtype().kind();
    if (kind == 'O') {
        raise_error("TypeError", "Object dtype dtype('O') has no native HDF5 equivalent");
    }
    if (kind == 'U') {
        raise_error("TypeError", "No conversion path for dtype: dtype('" + dtype_name(array) + "')");
    }
    raise_error("TypeError", "HiCExplorer 4 writes the cHi-C lists as int64 or float64, not " + dtype_name(array));
}

// create_dataset(path, data=values[, compression='gzip', compression_opts=level])
inline void write_values(hicx::chic::Hdf5Writer& writer, const std::string& path, py::handle values, int level) {
    if (writer.exists(path)) {
        raise_error("ValueError", "Unable to synchronously create dataset (name already exists)");
    }
    const py::array array = as_array(values);
    if (array.ndim() > 1) {
        raise_error("TypeError", "HiCExplorer 4 writes the cHi-C lists as one dimensional datasets");
    }
    const char kind = array.dtype().kind();
    if (kind == 'i' && array.dtype().itemsize() == 8) {
        const std::vector<std::int64_t> data = to_vector<std::int64_t>(array);
        writer.write_array(path, std::span<const std::int64_t>(data), level);
        return;
    }
    if (kind == 'f' && array.dtype().itemsize() == 8) {
        const std::vector<double> data = to_vector<double>(array);
        writer.write_array(path, std::span<const double>(data), level);
        return;
    }
    if (kind == 'S' || kind == 'O') {
        // a list of str becomes '<U', of bytes 'S': only variable length
        // strings (object arrays of str) are written as h5py would
        bool all_str = array.size() > 0;
        std::vector<std::string> strings;
        for (const py::handle item : values) {
            if (!is_str(item)) {
                all_str = false;
                break;
            }
            strings.push_back(utf8(item));
        }
        if (all_str && kind == 'O') {
            writer.write_strings(path, strings);
            return;
        }
    }
    unsupported_dtype(array);
}

// group[name] = value for a scalar
inline void write_scalar(hicx::chic::Hdf5Writer& writer, const std::string& path, py::handle value) {
    if (writer.exists(path)) {
        raise_error("ValueError", "Unable to synchronously create dataset (name already exists)");
    }
    if (is_str(value)) {
        writer.write_string(path, utf8(value));
        return;
    }
    if (PyFloat_Check(value.ptr()) ||
        PyObject_IsInstance(value.ptr(), numpy().attr("float64").ptr()) == 1) {
        writer.write_scalar(path, to_double(value));
        return;
    }
    if (!PyBool_Check(value.ptr()) && (PyLong_Check(value.ptr()) ||
                                       PyObject_IsInstance(value.ptr(), numpy().attr("int64").ptr()) == 1)) {
        writer.write_scalar(path, to_int64(value));
        return;
    }
    if (value.is_none()) {
        raise_error("TypeError", "Object dtype dtype('O') has no native HDF5 equivalent");
    }
    unsupported_dtype(as_array(value));
}

// attrs[name] = value
inline void write_attribute(hicx::chic::Hdf5Writer& writer, const std::string& object, const std::string& name,
                            py::handle value) {
    if (is_str(value)) {
        writer.set_attribute(object, name, utf8(value));
        return;
    }
    if (PyBool_Check(value.ptr()) || PyObject_IsInstance(value.ptr(), numpy().attr("bool_").ptr()) == 1) {
        writer.set_bool_attribute(object, name, truthy(value));
        return;
    }
    if (PyFloat_Check(value.ptr()) || PyObject_IsInstance(value.ptr(), numpy().attr("float64").ptr()) == 1) {
        writer.set_attribute(object, name, to_double(value));
        return;
    }
    if (PyLong_Check(value.ptr()) || PyObject_IsInstance(value.ptr(), numpy().attr("int64").ptr()) == 1) {
        writer.set_attribute(object, name, to_int64(value));
        return;
    }
    if (PyList_Check(value.ptr()) || PyTuple_Check(value.ptr()) || is_ndarray(value)) {
        const py::array array = as_array(value);
        if (array.dtype().kind() == 'i' && array.dtype().itemsize() == 8 && array.ndim() == 1) {
            const std::vector<std::int64_t> data = to_vector<std::int64_t>(array);
            writer.set_attribute(object, name, std::span<const std::int64_t>(data));
            return;
        }
        unsupported_dtype(array);
    }
    raise_error("TypeError", "Object dtype dtype('O') has no native HDF5 equivalent");
}

}  // namespace hicx::python::chic_api

#endif  // HICX_PYTHON_API_CHIC_COMMON_HPP
