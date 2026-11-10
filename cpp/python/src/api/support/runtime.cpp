#include "runtime.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>

#include "convert.hpp"
#include "errors.hpp"
#include "hicx/py_text.hpp"

namespace hicx::python {

void log(const char* logger, const char* level, std::string_view message) {
    py::module_::import("logging").attr("getLogger")(logger).attr(level)(to_str(message));
}

void warn(PyObject* category, std::string_view message, int stacklevel) {
    const std::string text(message);
    if (PyErr_WarnEx(category, text.c_str(), stacklevel) != 0) {
        throw_python_error();
    }
}

namespace {

// os.fspath(path), the object open() names in its errors.
py::object fspath_object(py::handle path) {
    PyObject* resolved = PyOS_FSPath(path.ptr());
    if (resolved == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(resolved);
}

std::string encoded_path(const py::object& resolved) {
    std::string bytes;
    if (PyBytes_Check(resolved.ptr())) {
        bytes.assign(PyBytes_AS_STRING(resolved.ptr()),
                     static_cast<std::size_t>(PyBytes_GET_SIZE(resolved.ptr())));
    } else {
        PyObject* encoded = PyUnicode_EncodeFSDefault(resolved.ptr());
        if (encoded == nullptr) {
            throw_python_error();
        }
        bytes.assign(PyBytes_AS_STRING(encoded), static_cast<std::size_t>(PyBytes_GET_SIZE(encoded)));
        Py_DECREF(encoded);
    }
    if (bytes.find('\0') != std::string::npos) {
        raise_error("ValueError", "embedded null byte");
    }
    return bytes;
}

std::string read_descriptor(int fd, py::handle filename) {
    std::string out;
    char buffer[1 << 16];
    while (true) {
        const ssize_t got = ::read(fd, buffer, sizeof(buffer));
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            const int error = errno;
            ::close(fd);
            raise_os_error(error, filename);
        }
        if (got == 0) {
            break;
        }
        out.append(buffer, static_cast<std::size_t>(got));
    }
    ::close(fd);
    return out;
}

bool locale_encoding_is_utf8(std::string* name) {
    static int cached = -1;
    static std::string cached_name;
    if (cached < 0) {
        const py::object encoding = py::module_::import("locale").attr("getpreferredencoding")(false);
        const py::object info = py::module_::import("codecs").attr("lookup")(encoding);
        cached_name = py::str(info.attr("name"));
        cached = cached_name == "utf-8" ? 1 : 0;
    }
    if (name != nullptr) {
        *name = cached_name;
    }
    return cached == 1;
}

// Raises the UnicodeDecodeError TextIOWrapper would: it decodes 8192 byte
// chunks with an incremental decoder that carries an incomplete sequence over
// to the next chunk and decodes the remainder with final=True at the end.
[[noreturn]] void raise_utf8_error(const std::string& bytes) {
    constexpr std::size_t kChunk = 8192;
    std::string pending;
    for (std::size_t offset = 0; offset < bytes.size(); offset += kChunk) {
        std::string data = pending + bytes.substr(offset, kChunk);
        Py_ssize_t consumed = 0;
        PyObject* decoded = PyUnicode_DecodeUTF8Stateful(
            data.data(), static_cast<Py_ssize_t>(data.size()), "strict", &consumed);
        if (decoded == nullptr) {
            throw_python_error();
        }
        Py_DECREF(decoded);
        pending = data.substr(static_cast<std::size_t>(consumed));
    }
    PyObject* decoded =
        PyUnicode_DecodeUTF8(pending.data(), static_cast<Py_ssize_t>(pending.size()), "strict");
    if (decoded == nullptr) {
        throw_python_error();
    }
    Py_DECREF(decoded);
    raise_error("UnicodeDecodeError", "invalid UTF-8 that Python accepts");
}

}  // namespace

std::string fs_path(py::handle path) {
    return encoded_path(fspath_object(path));
}

std::string read_file(py::handle path) {
    if (PyLong_Check(path.ptr())) {
        const long long fd = PyLong_AsLongLong(path.ptr());
        if (fd == -1 && PyErr_Occurred()) {
            throw_python_error();
        }
        if (fd < 0) {
            raise_error("ValueError", "negative file descriptor");
        }
        struct stat info {};
        if (::fstat(static_cast<int>(fd), &info) != 0) {
            raise_os_error(errno, py::none());
        }
        return read_descriptor(static_cast<int>(fd), py::none());
    }
    const py::object resolved = fspath_object(path);
    const std::string name = encoded_path(resolved);
    const int fd = ::open(name.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) {
        raise_os_error(errno, resolved);
    }
    struct stat info {};
    if (::fstat(fd, &info) == 0 && S_ISDIR(info.st_mode)) {
        ::close(fd);
        raise_os_error(EISDIR, resolved);
    }
    return read_descriptor(fd, resolved);
}

std::vector<std::string> read_text_lines(py::handle path) {
    const std::string bytes = read_file(path);
    std::string encoding;
    if (locale_encoding_is_utf8(&encoding)) {
        if (pytext::first_invalid_utf8(bytes)) {
            raise_utf8_error(bytes);
        }
        return pytext::universal_newline_lines(bytes);
    }
    PyObject* decoded = PyUnicode_Decode(bytes.data(), static_cast<Py_ssize_t>(bytes.size()),
                                         encoding.c_str(), "strict");
    if (decoded == nullptr) {
        throw_python_error();
    }
    const py::object text = py::reinterpret_steal<py::object>(decoded);
    return pytext::universal_newline_lines(utf8(text));
}

}  // namespace hicx::python
