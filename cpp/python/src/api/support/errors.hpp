// Python exceptions from C++ for hicexplorer._hicx.
//
// The C++ code throws hicx::PythonError(type, message) (core or bindings);
// the translator registered by the module raises it as that Python type with
// that message. A pybind11::error_already_set (a Python error raised by a
// Python API call made from C++) passes through unchanged.

#ifndef HICX_PYTHON_API_ERRORS_HPP
#define HICX_PYTHON_API_ERRORS_HPP

#include <string>
#include <string_view>

#include <pybind11/pybind11.h>

#include "hicx/python_error.hpp"

namespace hicx::python {

namespace py = pybind11;

// Installs the hicx::PythonError translator. Called once by the module.
void register_exception_translator();

// The exception class of `type`: a builtin name ("ValueError") or a dotted
// path ("argparse.ArgumentTypeError"), imported on first use.
[[nodiscard]] py::object exception_type(const std::string& type);

// Throws hicx::PythonError, raised in Python as `type` with `message`.
[[noreturn]] void raise_error(const std::string& type, std::string_view message);

// Raises the OSError subclass for `error_number` as open() does, with
// `filename` (None for none) as its filename.
[[noreturn]] void raise_os_error(int error_number, py::handle filename);

// The current Python exception (after a failed C API call) as a C++ throw.
[[noreturn]] void throw_python_error();

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_ERRORS_HPP
