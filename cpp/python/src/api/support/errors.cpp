#include "errors.hpp"

#include <cerrno>

#include "convert.hpp"

namespace hicx::python {

py::object exception_type(const std::string& type) {
    const std::size_t dot = type.rfind('.');
    if (dot == std::string::npos) {
        return py::module_::import("builtins").attr(type.c_str());
    }
    return py::module_::import(type.substr(0, dot).c_str()).attr(type.substr(dot + 1).c_str());
}

void register_exception_translator() {
    py::register_exception_translator([](std::exception_ptr pointer) {
        if (!pointer) {
            return;
        }
        try {
            std::rethrow_exception(pointer);
        } catch (const hicx::PythonError& error) {
            try {
                const py::object type = exception_type(error.type());
                const py::str message = to_str(error.what());
                PyErr_SetObject(type.ptr(), message.ptr());
            } catch (py::error_already_set& lookup) {
                lookup.restore();
            }
        }
    });
}

void raise_error(const std::string& type, std::string_view message) {
    throw hicx::PythonError(type, std::string(message));
}

void raise_os_error(int error_number, py::handle filename) {
    errno = error_number;
    if (filename.is_none()) {
        PyErr_SetFromErrno(PyExc_OSError);
    } else {
        PyErr_SetFromErrnoWithFilenameObject(PyExc_OSError, filename.ptr());
    }
    throw py::error_already_set();
}

void throw_python_error() {
    throw py::error_already_set();
}

}  // namespace hicx::python
