// hicexplorer._hicx: the C++ implementation of HiCExplorer's public Python
// API. Each bound module of src/api/modules.cmake is a submodule,
// hicexplorer._hicx.<entry>, also importable under that name.

#include <pybind11/pybind11.h>

#include "support/convert.hpp"
#include "support/errors.hpp"

namespace py = pybind11;

namespace hicx::python {
#define HICX_API_MODULE(name) void bind_##name(py::module_& module);
#include "api_registry.inc"
#undef HICX_API_MODULE
}  // namespace hicx::python

PYBIND11_MODULE(_hicx, module) {
    module.doc() = "C++ implementation of the public Python API of HiCExplorer 4.";
    hicx::python::register_exception_translator();
    py::object modules = py::module_::import("sys").attr("modules");
    const std::string prefix = py::str(module.attr("__name__"));
    py::list names;
#define HICX_API_MODULE(name)                                       \
    {                                                               \
        py::module_ submodule = module.def_submodule(#name);        \
        hicx::python::bind_##name(submodule);                       \
        modules[py::str(prefix + "." #name)] = submodule;           \
        names.append(#name);                                        \
    }
#include "api_registry.inc"
#undef HICX_API_MODULE
    module.attr("MODULES") = py::tuple(names);

    // For the tests of the conversion helpers (hicexplorer/test/api/test_convert.py):
    // reads a hiCMatrix-like object into hicx::MatrixData and writes it back.
    module.def("_roundtrip_matrix_data", [](py::handle source, py::handle target) {
        hicx::python::assign_matrix_data(target, hicx::python::matrix_data(source));
    });
}
