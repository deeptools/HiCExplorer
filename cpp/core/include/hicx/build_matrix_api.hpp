// The C++ behind the public Python names of hicexplorer.lib.buildMatrixMethods
// and hicexplorer.hicBuildMatrixMicroC of HiCExplorer 3.7.6 that the
// hicBuildMatrix machinery (hicx/build_matrix.hpp, hicx/build_matrix_run.hpp)
// does not already provide in the form the Python functions have
// (cpp/python/API_BINDINGS.md). No Python headers.
//
// Every function reproduces its Python counterpart, including its defects;
// those are marked QUIRK with the line of buildMatrixMethods.py.

#ifndef HICX_BUILD_MATRIX_API_HPP
#define HICX_BUILD_MATRIX_API_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hicx::build_matrix_api {

// --- paths -------------------------------------------------------------------

// posixpath.normpath, abspath (relative to `cwd`), dirname and join.
[[nodiscard]] std::string normpath(const std::string& path);
[[nodiscard]] std::string abspath(const std::string& path, const std::string& cwd);
[[nodiscard]] std::string dirname(const std::string& path);
[[nodiscard]] std::string join(const std::string& a, const std::string& b);

// Tester.__init__: root is $HIC_TEST_DATA_DIR when that is set and not empty,
// else the directory `levels` dirname()s above abspath(module_file) plus
// "/test/test_data/"; bam_file_1 is join(root, "hic.bam").
struct TesterPaths {
    std::string root;
    std::string bam_file_1;
};
[[nodiscard]] TesterPaths tester_paths(const std::optional<std::string>& test_data_dir,
                                       const std::string& module_file, const std::string& cwd,
                                       int levels);

}  // namespace hicx::build_matrix_api

#endif  // HICX_BUILD_MATRIX_API_HPP
