// hicexplorer._hicx.chicViewpointBackgroundModel: the C++ body of
// compute_background (main() and parse_arguments() run the executable
// cpp/tools/chicViewpointBackgroundModel.cpp).
//
// compute_background collects, per relative position to the reference point,
// the (smoothed) interactions of every reference point; the viewpoints and
// the smoothing are the Viewpoint object's methods, called as 3.7.6 called
// them, and the values collected are the elements of the arrays they return.

#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"

namespace hicx::python {

using namespace chic_api;

void bind_chicViewpointBackgroundModel(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicViewpointBackgroundModel.";

    m.def("compute_background", [](py::object pReferencePoints, py::object pViewpointObj, py::object pArgs,
                                   py::object pQueue) {
        py::dict background_model_data;
        py::set relative_positions;
        try {
            for (const py::handle referencePoint : pReferencePoints) {
                const py::object fixate = pArgs.attr("fixateRange");
                const std::vector<py::object> region = unpack(
                    pViewpointObj.attr("calculateViewpointRange")(referencePoint, py::make_tuple(fixate, fixate)), 3);
                const py::object chromosome = getitem(referencePoint, 0);
                py::object data_list =
                    unpack(pViewpointObj.attr("computeViewpoint")(referencePoint, chromosome, region[0], region[1]), 2)[0];
                const py::object view_point_start =
                    unpack(pViewpointObj.attr("getReferencePointAsMatrixIndices")(referencePoint), 2)[0];
                const std::vector<py::object> range = unpack(
                    pViewpointObj.attr("getViewpointRangeAsMatrixIndices")(chromosome, region[0], region[1]), 2);
                const std::int64_t first = to_int64(range[0]);
                const std::int64_t last = to_int64(range[1]);
                // the first zip over the data only computes relative positions
                {
                    std::int64_t i = first;
                    for (const py::handle item : data_list) {
                        (void)item;
                        if (i >= last) {
                            break;
                        }
                        (void)sub(py::int_(i), view_point_start);
                        ++i;
                    }
                }
                const py::object average = pArgs.attr("averageContactBin");
                if (gt(average, py::int_(0))) {
                    data_list = pViewpointObj.attr("smoothInteractionValues")(data_list, average);
                }
                std::int64_t i = first;
                for (const py::handle data : data_list) {
                    if (i >= last) {
                        break;
                    }
                    const py::object relative_position = sub(py::int_(i), view_point_start);
                    if (contains(background_model_data, relative_position)) {
                        getitem(background_model_data, relative_position).attr("append")(data);
                    } else {
                        py::list values;
                        values.append(data);
                        background_model_data[relative_position] = values;
                        relative_positions.add(relative_position);
                    }
                    ++i;
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
        py::list result;
        result.append(background_model_data);
        result.append(relative_positions);
        pQueue.attr("put")(result);
    }, py::arg("pReferencePoints"), py::arg("pViewpointObj"), py::arg("pArgs"), py::arg("pQueue"));
}

}  // namespace hicx::python
