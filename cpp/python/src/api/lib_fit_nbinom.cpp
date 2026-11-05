// hicexplorer._hicx.lib_fit_nbinom: the C++ body of hicexplorer.lib.fit_nbinom.
//
// fit(X, initial_params=None) is the maximum likelihood fit of the negative
// binomial of the external fit_nbinom 1.2 package that 3.7.6 imported, which
// HiCExplorer 4 carries itself. The objective, the moment estimator start
// and the bounds are those of that package. A float64 or integer X is fitted
// by the direct maximum likelihood solution of hicx/nbinom_mle.hpp (the one
// chicViewpointBackgroundModel uses), where the package ran an L-BFGS-B
// optimiser; a float32 X takes the single precision objective and own
// projected optimiser of hicx/stats_ops.hpp (the one hicDetectLoops uses),
// because numpy 1.26 evaluates two of the five terms in float32 for such an
// array.
//
// Deviations from 3.7.6: int8, int16, uint8, uint16, uint32, uint64, float16,
// longdouble, complex and non-native byte order arrays raise TypeError
// naming HiCExplorer 4; a list or another object that is not an ndarray
// raises the AttributeError 3.7.6 raised for the missing `.size` without
// evaluating the earlier numpy call on it.

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/nbinom_mle.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::python {

namespace {

py::object fit(py::handle X, py::handle initial_params) {
    if (X.is_none()) {
        raise_error("TypeError", "unsupported operand type(s) for /: 'NoneType' and 'int'");
    }
    if (!is_ndarray(X)) {
        raise_error("AttributeError", "'" + type_name(X) + "' object has no attribute 'size'");
    }
    const py::array array = py::reinterpret_borrow<py::array>(X);

    std::optional<std::array<double, 2>> initial;
    if (!initial_params.is_none()) {
        py::object flat = py::module_::import("numpy").attr("asarray")(initial_params).attr("ravel")();
        if (py::len(flat) != 2) {
            raise_error("ValueError", "length of x0 != length of bounds");
        }
        const py::array values = py::reinterpret_steal<py::array>(
            flat.attr("astype")(py::module_::import("numpy").attr("float64")).release());
        const std::vector<double> pair = to_vector<double>(values);
        initial = std::array<double, 2>{pair[0], pair[1]};
    }

    std::vector<double> data;
    bool single = false;
    switch (dtype_of(array)) {
        case DType::Int32:
            for (const auto v : to_vector<std::int32_t>(array)) data.push_back(static_cast<double>(v));
            break;
        case DType::Int64:
            for (const auto v : to_vector<std::int64_t>(array)) data.push_back(static_cast<double>(v));
            break;
        case DType::Float32:
            single = true;
            for (const auto v : to_vector<float>(array)) data.push_back(static_cast<double>(v));
            break;
        case DType::Float64:
            data = to_vector<double>(array);
            break;
        case DType::Bool:
            raise_error("ValueError",
                        "Unsupported datatype for factorial: bool\n"
                        "Permitted data types are integers and floating point numbers");
        default:
            raise_error("TypeError", "hicexplorer.lib.fit_nbinom.fit: HiCExplorer 4 fits int32, int64, float32 and "
                                     "float64 arrays, not " + dtype_name(array));
    }

    hicx::stats::NBinomFit result;
    if (single && !data.empty()) {
        result = hicx::stats::fit_nbinom(data, hicx::stats::NBinomPrecision::Float32, initial);
    } else {
        result = hicx::stats::fit_nbinom_mle(data, initial);
    }
    py::dict out;
    out["size"] = numpy_scalar<double>(result.size);
    out["prob"] = numpy_scalar<double>(result.prob);
    return std::move(out);
}

}  // namespace

void bind_lib_fit_nbinom(py::module_& m) {
    m.def("fit", &fit, py::arg("X"), py::arg("initial_params") = py::none());
}

}  // namespace hicx::python
