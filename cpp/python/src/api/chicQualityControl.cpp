// hicexplorer._hicx.chicQualityControl: the C++ body of compute_sparsity
// (main() and parse_arguments() run the executable
// cpp/tools/chicQualityControl.cpp).
//
// The sparsity of a reference point is the share of non-zero elements of its
// viewpoint over the fixate range (np.count_nonzero(data) / len(data)), -1.0
// for a reference point on a chromosome the matrix does not have or whose
// viewpoint raises TypeError or IndexError. The viewpoint is the Viewpoint
// object's computeViewpoint, called as 3.7.6 called it.

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"

namespace hicx::python {

using namespace chic_api;

namespace {

// np.count_nonzero(data): NaN is not zero.
std::int64_t count_nonzero(py::handle data) {
    const py::array array = py::reinterpret_borrow<py::array>(
        numpy().attr("asarray")(data, py::arg("dtype") = numpy().attr("float64")));
    std::int64_t count = 0;
    for (const double value : to_vector<double>(array)) {
        count += value != 0.0 ? 1 : 0;
    }
    return count;
}

}  // namespace

void bind_chicQualityControl(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicQualityControl.";

    m.def("compute_sparsity", [](py::object pReferencePoints, py::object pViewpointObj, py::object pArgs,
                                 py::object pQueue) {
        py::list sparsity_list;
        try {
            const py::object chromosome_names = pViewpointObj.attr("hicMatrix").attr("getChrNames")();
            for (const py::handle referencePoint : pReferencePoints) {
                if (!referencePoint.is_none() && contains(chromosome_names, getitem(referencePoint, 0))) {
                    const py::object fixate = pArgs.attr("fixateRange");
                    const std::vector<py::object> region = unpack(
                        pViewpointObj.attr("calculateViewpointRange")(referencePoint, py::make_tuple(fixate, fixate)),
                        3);
                    py::object sparsity;
                    try {
                        const py::object data_list = unpack(pViewpointObj.attr("computeViewpoint")(
                                                                referencePoint, getitem(referencePoint, 0), region[0],
                                                                region[1]),
                                                            2)[0];
                        const std::int64_t nonzero = count_nonzero(data_list);
                        const py::ssize_t n = length(data_list);
                        if (n == 0) {
                            raise_error("ZeroDivisionError", "division by zero");
                        }
                        sparsity = py::float_(static_cast<double>(nonzero) / static_cast<double>(n));
                    } catch (py::error_already_set& error) {
                        if (!error.matches(PyExc_TypeError) && !error.matches(PyExc_IndexError)) {
                            throw;
                        }
                        sparsity = py::float_(-1.0);
                    } catch (const hicx::PythonError& error) {
                        if (error.type() != "TypeError" && error.type() != "IndexError") {
                            throw;
                        }
                        sparsity = py::float_(-1.0);
                    }
                    sparsity_list.append(sparsity);
                } else {
                    sparsity_list.append(py::float_(-1.0));
                }
            }
        } catch (py::error_already_set& error) {
            if (!is_exception(error)) {
                throw;
            }
            pQueue.attr("put")(to_str("Fail: " + str_of(error.value())));
            return;
        } catch (const std::exception& error) {
            pQueue.attr("put")(to_str("Fail: " + fail_text(error)));
            return;
        }
        pQueue.attr("put")(sparsity_list);
    }, py::arg("pReferencePoints"), py::arg("pViewpointObj"), py::arg("pArgs"), py::arg("pQueue"));
}

}  // namespace hicx::python
