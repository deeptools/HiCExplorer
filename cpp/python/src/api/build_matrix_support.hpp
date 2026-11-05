// Shared by the bindings of hicexplorer.lib.buildMatrixMethods and
// hicexplorer.hicBuildMatrixMicroC (lib_buildMatrixMethods.cpp,
// hicBuildMatrixMicroC.cpp).

#ifndef HICX_PYTHON_API_BUILD_MATRIX_SUPPORT_HPP
#define HICX_PYTHON_API_BUILD_MATRIX_SUPPORT_HPP

#include <unistd.h>

#include <cerrno>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/build_matrix_api.hpp"

namespace hicx::python {

// Tester.__init__ of both modules: `module_file` is the module's __file__,
// `levels` how many dirname()s lead from it to the hicexplorer package.
inline void tester_init(py::handle self, py::handle module_file, int levels) {
    const py::object value = py::module_::import("os").attr("environ").attr("get")("HIC_TEST_DATA_DIR", false);
    std::optional<std::string> test_data_dir;
    if (PyObject_IsTrue(value.ptr()) == 1) {
        if (!is_str(value)) {
            raise_error("TypeError", "Tester: HiCExplorer 4 takes HIC_TEST_DATA_DIR as a str");
        }
        test_data_dir = utf8(value);
    }
    std::vector<char> cwd(4096);
    while (::getcwd(cwd.data(), cwd.size()) == nullptr) {
        if (errno != ERANGE) {
            raise_os_error(errno, py::none());
        }
        cwd.resize(cwd.size() * 2);
    }
    const hicx::build_matrix_api::TesterPaths paths =
        hicx::build_matrix_api::tester_paths(test_data_dir, utf8(module_file), cwd.data(), levels);
    self.attr("root") = to_str(paths.root);
    self.attr("bam_file_1") = to_str(paths.bam_file_1);
}

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_BUILD_MATRIX_SUPPORT_HPP
