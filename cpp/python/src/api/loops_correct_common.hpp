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
#include "hicx/coo_ops.hpp"
#include "hicx/utilities.hpp"

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

// matrix.eliminate_zeros() of a csr or csc matrix, in place.
inline void eliminate_zeros(py::handle matrix_object, std::string_view function) {
    Sparse matrix = sparse(matrix_object, function, "rc");
    dispatch_arithmetic(matrix.data, function, [&]<class T>() {
        const std::vector<T> data = to_vector<T>(matrix.data);
        std::vector<std::int64_t> indptr = matrix.indptr_vector();
        const std::vector<std::size_t> keep = hicx::utilities::eliminate_zeros(indptr, data.data());
        if (keep.size() != data.size()) {
            apply_eliminate_zeros(matrix, indptr, keep);
        }
    });
}

// A csr matrix in its element type, in storage order.
template <class T>
hicx::coo_ops::Csr<T> csr_of(const Sparse& matrix) {
    hicx::coo_ops::Csr<T> out;
    out.rows = matrix.rows;
    out.cols = matrix.cols;
    out.indptr = matrix.indptr_vector();
    const hicx::utilities::SparseEntries entries = matrix.entries();
    out.indices = entries.col;
    out.data = to_vector<T>(matrix.data);
    return out;
}

// scipy.sparse.csr_matrix((data, indices, indptr), shape) with int32 indices.
template <class T>
py::object csr_object(const hicx::coo_ops::Csr<T>& matrix) {
    std::vector<std::int32_t> indices(matrix.indices.begin(), matrix.indices.end());
    std::vector<std::int32_t> indptr(matrix.indptr.begin(), matrix.indptr.end());
    return py::module_::import("scipy.sparse").attr("csr_matrix")(
        py::make_tuple(from_vector(matrix.data), from_vector(indices), from_vector(indptr)),
        py::arg("shape") = py::make_tuple(matrix.rows, matrix.cols));
}

}  // namespace hicx::python::lc

#endif  // HICX_PYTHON_API_LOOPS_CORRECT_COMMON_HPP
