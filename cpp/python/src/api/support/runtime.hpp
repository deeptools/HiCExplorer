// What the 3.7.6 functions did through the Python runtime besides
// computing: logging, warnings and reading files with open()'s errors.

#ifndef HICX_PYTHON_API_RUNTIME_HPP
#define HICX_PYTHON_API_RUNTIME_HPP

#include <string>
#include <string_view>
#include <vector>

#include <pybind11/pybind11.h>

namespace hicx::python {

namespace py = pybind11;

// logging.getLogger(logger).<level>(message), level one of "debug", "info",
// "warning", "error". The logger is the module's (`log = getLogger(__name__)`
// of the Python module, e.g. "hicexplorer.utilities").
void log(const char* logger, const char* level, std::string_view message);

// warnings.warn(message, category, stacklevel).
void warn(PyObject* category, std::string_view message, int stacklevel = 1);

// The bytes of a file opened as open(path, 'r') or open(path, 'rb') opens
// it: str, bytes and os.PathLike paths and integer file descriptors (closed
// afterwards, as the with block of the Python closes them), with the same
// OSError subclasses, errno and filename, and ValueError for an embedded
// null byte.
[[nodiscard]] std::string read_file(py::handle path);

// file.readlines() of open(path, 'r'): decoded with the locale's encoding
// (strict) and with universal newlines, returned as UTF-8 (surrogatepass).
// A decoding error is the UnicodeDecodeError TextIOWrapper raises, with the
// positions it reports (it decodes 8192 byte chunks).
[[nodiscard]] std::vector<std::string> read_text_lines(py::handle path);

// The path as open() resolves it, for the C++ side: os.fspath() and, for
// bytes, os.fsdecode(). Raises TypeError like open() for anything else.
[[nodiscard]] std::string fs_path(py::handle path);

}  // namespace hicx::python

#endif  // HICX_PYTHON_API_RUNTIME_HPP
