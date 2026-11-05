// hicexplorer._hicx.parserCommon: the C++ body of
// hicexplorer.parserCommon.writableFile (hicx/parser_common.hpp).
//
// getParentArgParse and CustomFormatter stay Python: the first only builds an
// argparse.ArgumentParser with one option, the second is an argparse
// formatter class without a body. Neither computes anything; both are listed
// in cpp/scripts/public_api_allowlist.tsv with that reason.

#include <string>

#include "api.hpp"
#include "hicx/parser_common.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.parserCommon";

// format(value, ''), which "{}".format(value) calls.
std::string formatted(py::handle value) {
    PyObject* empty = PyUnicode_FromString("");
    if (empty == nullptr) {
        throw_python_error();
    }
    PyObject* text = PyObject_Format(value.ptr(), empty);
    Py_DECREF(empty);
    if (text == nullptr) {
        throw_python_error();
    }
    return utf8(py::reinterpret_steal<py::object>(text));
}

}  // namespace

void bind_parserCommon(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.parserCommon.";

    m.def("writableFile", [](py::object string) -> py::object {
        int error = 0;
        if (PyLong_Check(string.ptr())) {
            // open(fd, 'w'): a descriptor
            const long long fd = PyLong_AsLongLong(string.ptr());
            if (fd == -1 && PyErr_Occurred()) {
                throw_python_error();
            }
            if (fd < 0) {
                raise_error("ValueError", "negative file descriptor");
            }
            error = hicx::parser_common::close_for_writing(static_cast<int>(fd));
        } else {
            error = hicx::parser_common::truncate_for_writing(fs_path(string));
        }
        if (error != 0) {
            const std::string message = hicx::parser_common::writable_file_message(formatted(string));
            log(kLogger, "debug", message);
            raise_error("argparse.ArgumentTypeError", message);
        }
        return string;
    }, py::arg("string"));
}

}  // namespace hicx::python
