// Helpers shared by the bindings of hicexplorer.iterativeCorrection,
// reduceMatrix, hicMergeMatrixBins, hicMergeLoops, hicCorrectMatrix,
// hicAdjustMatrix and hicDetectLoops (header only: every bound module is one
// translation unit of hicexplorer._hicx).

#ifndef HICX_PYTHON_API_LOOPS_CORRECT_COMMON_HPP
#define HICX_PYTHON_API_LOOPS_CORRECT_COMMON_HPP

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"

namespace hicx::python::lc {

inline bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool compare(py::handle a, py::handle b, int op) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), op);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool equal(py::handle a, py::handle b) {
    return compare(a, b, Py_EQ);
}

inline std::size_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

inline py::object item(py::handle object, py::handle key) {
    PyObject* value = PyObject_GetItem(object.ptr(), key.ptr());
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

inline py::object item(py::handle object, std::int64_t index) {
    return item(object, py::int_(index));
}

// exit(code) as the 3.7.6 code calls it: site's exit(), which closes
// sys.stdin and raises SystemExit(code).
[[noreturn]] inline void python_exit(int code) {
    py::module_::import("builtins").attr("exit")(code);
    throw std::logic_error("exit() returned");
}

// str() of a Python float, as str.format("{}") prints one.
inline std::string float_str(double value) {
    return str_of(py::float_(value));
}

// The values of a numpy array of any real or integer dtype as float64.
inline std::vector<double> doubles(const py::array& array, std::string_view function) {
    std::vector<double> out;
    dispatch_numeric(array, function, [&]<class T>() {
        if constexpr (std::is_arithmetic_v<T>) {
            const std::vector<T> values = to_vector<T>(array);
            out.assign(values.begin(), values.end());
        } else {
            raise_unsupported_dtype(array, function);
        }
    });
    return out;
}

// A monotonic clock in seconds, for the elapsed times the 3.7.6 code logs.
inline double now_seconds() {
    return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
}

}  // namespace hicx::python::lc

#endif  // HICX_PYTHON_API_LOOPS_CORRECT_COMMON_HPP
