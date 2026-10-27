// hicexplorer._hicx.lib_cnb: the C++ bodies of hicexplorer.lib.cnb.
//
// pdf and cdf are numpy expressions over scipy ufuncs in 3.7.6, so they take
// scalars, arrays and sequences, broadcast them against each other and
// return np.float64 for scalar (or 0-d) arguments and a float64 array
// otherwise. The element function is hicx::cnb (cpp/core/include/hicx/cnb.hpp).
//
// Deviation: the ufuncs evaluate bool, int8, int16, uint8, uint16, float16
// and float32 arguments in single precision; those dtypes (and complex,
// object and str arguments) raise TypeError naming HiCExplorer 4 here.

#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/cnb.hpp"

namespace hicx::python {

namespace {

struct Operand {
    std::vector<std::int64_t> shape;
    std::vector<double> values;  // C order
};

[[noreturn]] void unsupported(const char* function, const std::string& what) {
    raise_error("TypeError", std::string("hicexplorer.lib.cnb.") + function +
                                 ": HiCExplorer 4 evaluates float64 and int32/int64 arguments only, not " +
                                 what);
}

Operand operand(py::handle object, const char* function) {
    Operand out;
    if (PyBool_Check(object.ptr())) {
        unsupported(function, "bool");
    }
    if (PyLong_Check(object.ptr()) || PyFloat_Check(object.ptr())) {
        out.values.push_back(to_double(object));
        return out;
    }
    if (!is_ndarray(object) && !PyList_Check(object.ptr()) && !PyTuple_Check(object.ptr()) &&
        !py::isinstance(object, py::module_::import("numpy").attr("generic"))) {
        unsupported(function, "'" + type_name(object) + "'");
    }
    py::array array = py::array::ensure(object);
    if (!array) {
        throw_python_error();
    }
    switch (dtype_of(array)) {
        case DType::Int32: {
            for (const auto v : to_vector<std::int32_t>(array)) out.values.push_back(static_cast<double>(v));
            break;
        }
        case DType::Int64: {
            for (const auto v : to_vector<std::int64_t>(array)) out.values.push_back(static_cast<double>(v));
            break;
        }
        case DType::UInt32: {
            for (const auto v : to_vector<std::uint32_t>(array)) out.values.push_back(static_cast<double>(v));
            break;
        }
        case DType::UInt64: {
            for (const auto v : to_vector<std::uint64_t>(array)) out.values.push_back(static_cast<double>(v));
            break;
        }
        case DType::Float64:
            out.values = to_vector<double>(array);
            break;
        default:
            unsupported(function, "dtype " + dtype_name(array));
    }
    for (py::ssize_t d = 0; d < array.ndim(); ++d) {
        out.shape.push_back(array.shape(d));
    }
    return out;
}

std::string shape_text(const std::vector<std::int64_t>& shape) {
    std::string text = "(";
    for (std::size_t d = 0; d < shape.size(); ++d) {
        if (d > 0) {
            text += ",";
        }
        text += std::to_string(shape[d]);
    }
    if (shape.size() == 1) {
        text += ",";
    }
    return text + ")";
}

// numpy's broadcast of the shapes; std::nullopt-like empty flag on failure.
bool broadcast(const std::vector<const Operand*>& operands, std::vector<std::int64_t>& shape) {
    std::size_t ndim = 0;
    for (const Operand* o : operands) {
        ndim = std::max(ndim, o->shape.size());
    }
    shape.assign(ndim, 1);
    for (const Operand* o : operands) {
        const std::size_t offset = ndim - o->shape.size();
        for (std::size_t d = 0; d < o->shape.size(); ++d) {
            const std::int64_t n = o->shape[d];
            std::int64_t& target = shape[offset + d];
            if (n == target || n == 1) {
                continue;
            }
            if (target == 1) {
                target = n;
                continue;
            }
            return false;
        }
    }
    return true;
}

void require_broadcast(const std::vector<const Operand*>& operands) {
    std::vector<std::int64_t> shape;
    if (!broadcast(operands, shape)) {
        std::string message = "operands could not be broadcast together with shapes ";
        for (const Operand* o : operands) {
            message += shape_text(o->shape) + " ";
        }
        raise_error("ValueError", message);
    }
}

// The flat index into `o` of the output position `index` (C order of `shape`).
std::size_t source_index(const Operand& o, const std::vector<std::int64_t>& shape, std::int64_t index) {
    const std::size_t ndim = shape.size();
    const std::size_t offset = ndim - o.shape.size();
    std::size_t flat = 0;
    std::size_t stride = 1;
    for (std::size_t d = ndim; d-- > 0;) {
        const std::int64_t coordinate = index % shape[d];
        index /= shape[d];
        if (d >= offset) {
            const std::int64_t n = o.shape[d - offset];
            if (n != 1) {
                flat += static_cast<std::size_t>(coordinate) * stride;
            }
            stride *= static_cast<std::size_t>(n);
        }
    }
    return flat;
}

bool numeric_operand(py::handle object) {
    return PyLong_Check(object.ptr()) || PyFloat_Check(object.ptr()) || is_ndarray(object) ||
           py::isinstance(object, py::module_::import("numpy").attr("generic"));
}

// The Python `a + b` the 3.7.6 expression starts with, for operands that are
// not numbers or arrays: its TypeError (list + float, None + int), or, when
// it succeeds, nothing more is checked here.
void python_add(py::handle a, py::handle b, const char* function) {
    if (numeric_operand(a) && numeric_operand(b)) {
        return;
    }
    PyObject* sum = PyNumber_Add(a.ptr(), b.ptr());
    if (sum == nullptr) {
        throw_python_error();
    }
    const py::object result = py::reinterpret_steal<py::object>(sum);
    if (!numeric_operand(result)) {
        unsupported(function, "the concatenation '" + type_name(result) + "'");
    }
}

template <class F>
py::object evaluate(py::handle x_object, py::handle r_object, py::handle p_object, const char* function,
                    bool pdf, F&& element) {
    if (pdf) {
        python_add(r_object, x_object, function);
    }
    python_add(x_object, py::int_(1), function);
    const Operand x = operand(x_object, function);
    const Operand r = operand(r_object, function);
    const Operand p = operand(p_object, function);
    if (pdf) {
        // gammaln(pR + pX) first, then pR * np.log(pP), then the sum of all
        require_broadcast({&r, &x});
        require_broadcast({&r, &p});
        std::vector<std::int64_t> rx;
        broadcast({&r, &x}, rx);
        Operand rx_operand{rx, {}};
        std::vector<std::int64_t> rp;
        broadcast({&r, &p}, rp);
        Operand rp_operand{rp, {}};
        require_broadcast({&rx_operand, &rp_operand});
        if (!numeric_operand(p_object)) {
            // special.xlog1py(pX, -pP)
            PyObject* negated = PyNumber_Negative(p_object.ptr());
            if (negated == nullptr) {
                throw_python_error();
            }
            Py_DECREF(negated);
        }
    } else {
        require_broadcast({&r, &x, &p});
    }
    std::vector<std::int64_t> shape;
    broadcast({&x, &r, &p}, shape);
    std::int64_t total = 1;
    for (const std::int64_t n : shape) {
        total *= n;
    }
    std::vector<double> out(static_cast<std::size_t>(total));
    for (std::int64_t k = 0; k < total; ++k) {
        out[static_cast<std::size_t>(k)] =
            element(x.values[source_index(x, shape, k)], r.values[source_index(r, shape, k)],
                    p.values[source_index(p, shape, k)]);
    }
    if (shape.empty()) {
        return numpy_scalar<double>(out[0]);
    }
    py::array result = from_vector(out);
    std::vector<py::ssize_t> dims(shape.begin(), shape.end());
    return result.attr("reshape")(py::tuple(py::cast(dims)));
}

}  // namespace

void bind_lib_cnb(py::module_& m) {
    m.def("pdf", [](py::handle pX, py::handle pR, py::handle pP) {
        return evaluate(pX, pR, pP, "pdf", true, [](double x, double r, double p) { return hicx::cnb::pdf(x, r, p); });
    }, py::arg("pX"), py::arg("pR"), py::arg("pP"));

    m.def("cdf", [](py::handle pX, py::handle pR, py::handle pP) {
        return evaluate(pX, pR, pP, "cdf", false, [](double x, double r, double p) { return hicx::cnb::cdf(x, r, p); });
    }, py::arg("pX"), py::arg("pR"), py::arg("pP"));
}

}  // namespace hicx::python
