// hicexplorer._hicx.hicBuildMatrixMicroC: the C++ body of
// hicexplorer.hicBuildMatrixMicroC.Tester.__init__ (hicx/build_matrix_api.hpp).
// The module's other names come from hicexplorer.lib.buildMatrixMethods
// through its star import and are bound in lib_buildMatrixMethods.cpp.

#include "api.hpp"
#include "build_matrix_support.hpp"

namespace hicx::python {

void bind_hicBuildMatrixMicroC(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.hicBuildMatrixMicroC.";

    // root is dirname(abspath(__file__)) + "/test/test_data/": one level up
    // from hicexplorer/hicBuildMatrixMicroC.py.
    m.def("Tester___init__", [](py::handle self, py::handle module_file) {
        tester_init(self, module_file, 1);
    }, py::arg("self"), py::arg("module_file"));
}

}  // namespace hicx::python
