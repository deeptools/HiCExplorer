// Functions of the hicexplorer.utilities bindings that other bound modules
// call (hicexplorer.readBed calls toString).

#ifndef HICX_PYTHON_API_UTILITIES_HPP
#define HICX_PYTHON_API_UTILITIES_HPP

#include <pybind11/pybind11.h>

namespace hicx::python {

namespace py = pybind11;

// hicexplorer.utilities.toString(s)
[[nodiscard]] py::object to_string(py::handle s);

// hicexplorer.utilities.toBytes(s)
[[nodiscard]] py::object to_bytes(py::handle s);

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_UTILITIES_HPP
