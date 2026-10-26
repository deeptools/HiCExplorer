#include "convert.hpp"

namespace hicx::python {

// --------------------------------------------------------------------------
// str

bool is_str(py::handle object) noexcept {
    return PyUnicode_Check(object.ptr()) != 0;
}

std::string utf8(py::handle str) {
    PyObject* bytes = PyUnicode_AsEncodedString(str.ptr(), "utf-8", "surrogatepass");
    if (bytes == nullptr) {
        throw_python_error();
    }
    std::string out(PyBytes_AS_STRING(bytes), static_cast<std::size_t>(PyBytes_GET_SIZE(bytes)));
    Py_DECREF(bytes);
    return out;
}

py::str to_str(std::string_view utf8_text) {
    PyObject* str = PyUnicode_DecodeUTF8(utf8_text.data(), static_cast<Py_ssize_t>(utf8_text.size()),
                                         "surrogatepass");
    if (str == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::str>(str);
}

std::string str_of(py::handle object) {
    return utf8(py::str(object));
}

std::string repr_of(py::handle object) {
    return utf8(py::repr(object));
}

std::string type_name(py::handle object) {
    return Py_TYPE(object.ptr())->tp_name;
}

// --------------------------------------------------------------------------
// numpy

bool is_ndarray(py::handle object) noexcept {
    return py::isinstance<py::array>(object);
}

py::array as_ndarray(py::handle object) {
    if (!is_ndarray(object)) {
        raise_error("TypeError", "expected a numpy.ndarray, got " + type_name(object));
    }
    return py::reinterpret_borrow<py::array>(object);
}

DType dtype_of(const py::array& array) {
    const py::dtype dtype = array.dtype();
    const char kind = dtype.kind();
    const py::ssize_t size = dtype.itemsize();
    const char order = dtype.byteorder();
    const bool native = order == '=' || order == '|' || order == (PY_LITTLE_ENDIAN ? '<' : '>');
    if (!native) {
        return DType::Other;
    }
    switch (kind) {
        case 'b':
            return DType::Bool;
        case 'i':
            return size == 1 ? DType::Int8 : size == 2 ? DType::Int16 : size == 4 ? DType::Int32
                   : size == 8 ? DType::Int64 : DType::Other;
        case 'u':
            return size == 1 ? DType::UInt8 : size == 2 ? DType::UInt16 : size == 4 ? DType::UInt32
                   : size == 8 ? DType::UInt64 : DType::Other;
        case 'f':
            if (size == 2) return DType::Float16;
            if (size == 4) return DType::Float32;
            if (size == 8) return DType::Float64;
            if (static_cast<std::size_t>(size) == sizeof(long double)) return DType::LongDouble;
            return DType::Other;
        case 'c':
            if (size == 8) return DType::Complex64;
            if (size == 16) return DType::Complex128;
            return DType::Other;
        default:
            return DType::Other;
    }
}

std::string dtype_name(const py::array& array) {
    return str_of(array.dtype());
}

void raise_unsupported_dtype(const py::array& array, std::string_view function) {
    raise_error("TypeError", std::string(function) + ": arrays of dtype " + dtype_name(array) +
                                 " are not supported by HiCExplorer 4");
}

void require_writeable(const py::array& array) {
    if (!array.writeable()) {
        raise_error("ValueError", "assignment destination is read-only");
    }
}

// --------------------------------------------------------------------------
// scipy.sparse

namespace {

std::vector<std::int64_t> index_vector(const py::array& array) {
    switch (dtype_of(array)) {
        case DType::Int32: {
            const std::vector<std::int32_t> v = to_vector<std::int32_t>(array);
            return {v.begin(), v.end()};
        }
        case DType::Int64:
            return to_vector<std::int64_t>(array);
        case DType::UInt32: {
            const std::vector<std::uint32_t> v = to_vector<std::uint32_t>(array);
            return {v.begin(), v.end()};
        }
        default:
            raise_error("TypeError", "sparse index arrays of dtype " + dtype_name(array) +
                                         " are not supported by HiCExplorer 4");
    }
}

template <class I>
py::array take(const py::array& source, const std::vector<std::size_t>& keep) {
    const std::vector<I> values = to_vector<I>(source);
    py::array out(source.dtype(), std::vector<py::ssize_t>{static_cast<py::ssize_t>(keep.size())});
    I* target = static_cast<I*>(out.mutable_data());
    for (std::size_t i = 0; i < keep.size(); ++i) {
        target[i] = values[keep[i]];
    }
    return out;
}

}  // namespace

py::array take_elements(const py::array& source, const std::vector<std::size_t>& keep) {
    const std::size_t itemsize = static_cast<std::size_t>(source.itemsize());
    py::array contiguous = py::array::ensure(source, py::array::c_style);
    py::array out(source.dtype(), std::vector<py::ssize_t>{static_cast<py::ssize_t>(keep.size())});
    const char* from = static_cast<const char*>(contiguous.data());
    char* to = static_cast<char*>(out.mutable_data());
    for (std::size_t i = 0; i < keep.size(); ++i) {
        std::memcpy(to + i * itemsize, from + keep[i] * itemsize, itemsize);
    }
    return out;
}

Sparse sparse(py::handle object, std::string_view function, std::string_view formats) {
    const py::module_ scipy_sparse = py::module_::import("scipy.sparse");
    if (!scipy_sparse.attr("issparse")(object).cast<bool>()) {
        raise_error("TypeError", std::string(function) + ": expected a scipy.sparse matrix, got " +
                                     type_name(object));
    }
    Sparse m;
    m.matrix = py::reinterpret_borrow<py::object>(object);
    const std::string format = py::str(m.matrix.attr("format"));
    char code = 0;
    if (format == "csr") {
        code = 'r';
    } else if (format == "csc") {
        code = 'c';
    } else if (format == "coo") {
        code = 'o';
    }
    if (code == 0 || formats.find(code) == std::string_view::npos) {
        raise_error("TypeError", std::string(function) + ": scipy.sparse format '" + format +
                                     "' is not supported by HiCExplorer 4");
    }
    m.format = code;
    const py::tuple shape = m.matrix.attr("shape");
    m.rows = shape[0].cast<std::int64_t>();
    m.cols = shape[1].cast<std::int64_t>();
    m.data = m.matrix.attr("data");
    if (m.compressed()) {
        m.indptr = m.matrix.attr("indptr");
        m.indices = m.matrix.attr("indices");
    } else {
        m.row = m.matrix.attr("row");
        m.col = m.matrix.attr("col");
    }
    return m;
}

std::vector<std::int64_t> Sparse::indptr_vector() const {
    return index_vector(indptr);
}

hicx::utilities::SparseEntries Sparse::entries() const {
    hicx::utilities::SparseEntries e;
    e.rows = rows;
    e.cols = cols;
    const std::size_t stored = static_cast<std::size_t>(data.size());
    if (compressed()) {
        const std::vector<std::int64_t> pointer = index_vector(indptr);
        std::vector<std::int64_t> minor = index_vector(indices);
        std::vector<std::int64_t> major(minor.size(), 0);
        for (std::size_t i = 0; i + 1 < pointer.size(); ++i) {
            for (std::int64_t k = pointer[i]; k < pointer[i + 1] &&
                                              static_cast<std::size_t>(k) < major.size(); ++k) {
                major[static_cast<std::size_t>(k)] = static_cast<std::int64_t>(i);
            }
        }
        minor.resize(stored, 0);
        major.resize(stored, 0);
        if (format == 'r') {
            e.row = std::move(major);
            e.col = std::move(minor);
        } else {
            e.row = std::move(minor);
            e.col = std::move(major);
        }
    } else {
        e.row = index_vector(row);
        e.col = index_vector(col);
        e.row.resize(stored, 0);
        e.col.resize(stored, 0);
    }
    return e;
}

void set_data(Sparse& matrix, py::array data) {
    matrix.matrix.attr("data") = data;
    matrix.data = std::move(data);
}

void apply_eliminate_zeros(Sparse& matrix, const std::vector<std::int64_t>& indptr,
                           const std::vector<std::size_t>& keep) {
    // csr_eliminate_zeros rewrites indptr in place, prune() shortens the others
    if (dtype_of(matrix.indptr) == DType::Int32) {
        std::vector<std::int32_t> narrow(indptr.begin(), indptr.end());
        write_back(matrix.indptr, narrow);
    } else {
        write_back(matrix.indptr, indptr);
    }
    py::array indices = dtype_of(matrix.indices) == DType::Int32 ? take<std::int32_t>(matrix.indices, keep)
                                                                  : take<std::int64_t>(matrix.indices, keep);
    matrix.matrix.attr("indices") = indices;
    matrix.indices = indices;
    set_data(matrix, take_elements(matrix.data, keep));
}

// --------------------------------------------------------------------------
// Python numbers

bool is_integer(py::handle object) {
    return PyIndex_Check(object.ptr()) != 0;
}

std::int64_t to_int64(py::handle object) {
    PyObject* index = PyNumber_Index(object.ptr());
    if (index == nullptr) {
        throw_python_error();
    }
    const long long value = PyLong_AsLongLong(index);
    Py_DECREF(index);
    if (value == -1 && PyErr_Occurred()) {
        throw_python_error();
    }
    return value;
}

double to_double(py::handle object) {
    PyObject* number = PyNumber_Float(object.ptr());
    if (number == nullptr) {
        throw_python_error();
    }
    const double value = PyFloat_AS_DOUBLE(number);
    Py_DECREF(number);
    return value;
}

py::object parse_int(py::handle text) {
    PyObject* value = PyLong_FromUnicodeObject(text.ptr(), 10);
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

}  // namespace hicx::python
