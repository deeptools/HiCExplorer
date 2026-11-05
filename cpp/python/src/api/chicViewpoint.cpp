// hicexplorer._hicx.chicViewpoint: the C++ bodies of the functions of
// hicexplorer.chicViewpoint besides main() and parse_arguments() (which run
// the executable cpp/tools/chicViewpoint.cpp).
//
// compute_viewpoint calls the methods of the Viewpoint object it is given, as
// 3.7.6 did, and its own module's adjustViewpointData and compute_x_fold as the
// C++ functions here. The sum of the fixate range is numpy's pairwise sum
// (hicx::npy::pairwise_sum), the x-fold a division with numpy's dtype rules,
// the units of the header hicx::utilities::in_units.

#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/py_text.hpp"
#include "hicx/utilities.hpp"

namespace hicx::python {

using namespace chic_api;

namespace {

constexpr const char* kLogger = "hicexplorer.chicViewpoint";

py::tuple adjust_viewpoint_data(py::handle pViewpointObj, py::handle pData, py::handle pBackground,
                                py::handle pReferencePoint, py::handle pRegionStart, py::handle pRegionEnd) {
    const std::vector<py::object> start =
        unpack(pViewpointObj.attr("getReferencePointAsMatrixIndices")(pReferencePoint), 2);
    const std::vector<py::object> range = unpack(
        pViewpointObj.attr("getViewpointRangeAsMatrixIndices")(getitem(pReferencePoint, 0), pRegionStart, pRegionEnd),
        2);
    const std::int64_t first = to_int64(range[0]);
    const std::int64_t last = to_int64(range[1]);
    const std::int64_t positions = last > first ? last - first : 0;
    // zip(range(...), pData): the first `positions` elements of each
    std::vector<py::object> data;
    for (const py::handle item : pData) {
        if (static_cast<std::int64_t>(data.size()) >= positions) {
            break;
        }
        data.push_back(py::reinterpret_borrow<py::object>(item));
    }
    py::list background;
    for (const py::handle item : pBackground) {
        if (static_cast<std::int64_t>(background.size()) >= positions) {
            break;
        }
        background.append(item);
    }
    // the relative positions are distinct, so the dict holds the data in
    // order followed by a zero for each background position beyond them
    std::vector<float> values;
    for (const py::object& value : data) {
        values.push_back(static_cast<float>(to_double(value)));
    }
    for (std::size_t k = data.size(); k < background.size(); ++k) {
        values.push_back(0.0F);
    }
    (void)start;
    return py::make_tuple(from_vector(values), background);
}

// pDataList / pBackgroundList for one dimensional numbers (numpy's true
// division and dtype rules).
py::object x_fold(py::handle pDataList, py::handle pBackgroundList) {
    if (!is_ndarray(pDataList) && !is_ndarray(pBackgroundList)) {
        raise_error("TypeError", "unsupported operand type(s) for /: '" + type_name(pDataList) + "' and '" +
                                     type_name(pBackgroundList) + "'");
    }
    const py::object np = numpy();
    const py::array a = py::reinterpret_borrow<py::array>(np.attr("asarray")(pDataList));
    const py::array b = py::reinterpret_borrow<py::array>(np.attr("asarray")(pBackgroundList));
    const py::dtype type = py::reinterpret_borrow<py::dtype>(np.attr("result_type")(a, b));
    const char kind = type.kind();
    if (kind != 'f' && kind != 'i' && kind != 'u' && kind != 'b') {
        raise_error("TypeError", "compute_x_fold: dtype " + str_of(type) + " is not supported by HiCExplorer 4");
    }
    if (a.ndim() > 1 || b.ndim() > 1) {
        raise_error("TypeError", "compute_x_fold: only one dimensional data is supported by HiCExplorer 4");
    }
    const bool single = kind == 'f' && type.itemsize() == 4;
    const py::array a64 = py::reinterpret_borrow<py::array>(np.attr("asarray")(a, py::arg("dtype") = np.attr("float64")));
    const py::array b64 = py::reinterpret_borrow<py::array>(np.attr("asarray")(b, py::arg("dtype") = np.attr("float64")));
    const std::vector<double> x = to_vector<double>(a64);
    const std::vector<double> y = to_vector<double>(b64);
    const std::size_t n = std::max(x.size(), y.size());
    const bool broadcast_x = x.size() == 1 && a.ndim() <= 1;
    const bool broadcast_y = y.size() == 1 && b.ndim() <= 1;
    if (x.size() != y.size() && !broadcast_x && !broadcast_y) {
        raise_error("ValueError", "operands could not be broadcast together with shapes (" + std::to_string(x.size()) +
                                      ",) (" + std::to_string(y.size()) + ",) ");
    }
    const std::size_t count = (x.empty() || y.empty()) ? 0 : n;
    std::vector<double> out(count);
    std::vector<float> narrow(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double p = x[broadcast_x ? 0 : i];
        const double q = y[broadcast_y ? 0 : i];
        if (single) {
            narrow[i] = static_cast<float>(p) / static_cast<float>(q);
        } else {
            out[i] = p / q;
        }
    }
    py::array result = single ? from_vector(narrow) : from_vector(out);
    if (a.ndim() == 0 && b.ndim() == 0) {
        return getitem(result, 0);
    }
    return result;
}

// np.sum of what computeViewpoint returned
py::object array_sum(py::handle data) {
    if (is_ndarray(data)) {
        const py::array array = py::reinterpret_borrow<py::array>(data);
        if (dtype_of(array) == DType::Float64) {
            const std::vector<double> values = to_vector<double>(array);
            return np_float64(hicx::npy::pairwise_sum(values));
        }
    }
    return numpy().attr("sum")(data);
}

// np.copy(a)
py::object array_copy(py::handle data) {
    if (is_ndarray(data)) {
        const py::array array = py::reinterpret_borrow<py::array>(data);
        if (array.ndim() == 1 && dtype_of(array) == DType::Float64) {
            return from_vector(to_vector<double>(array));
        }
        if (array.ndim() == 1 && dtype_of(array) == DType::Float32) {
            return from_vector(to_vector<float>(array));
        }
    }
    return numpy().attr("copy")(data);
}

// array.flatten() of interactionBackgroundData's result
py::object flatten(py::handle data) {
    if (is_ndarray(data)) {
        const py::array array = py::reinterpret_borrow<py::array>(data);
        if (dtype_of(array) == DType::Float64) {
            return from_vector(to_vector<double>(array));
        }
    }
    return data.attr("flatten")();
}

}  // namespace

void bind_chicViewpoint(py::module_& m) {
    m.doc() = "C++ bodies of hicexplorer.chicViewpoint.";

    m.def("adjustViewpointData", &adjust_viewpoint_data, py::arg("pViewpointObj"), py::arg("pData"),
          py::arg("pBackground"), py::arg("pReferencePoint"), py::arg("pRegionStart"), py::arg("pRegionEnd"));

    m.def("compute_x_fold", &x_fold, py::arg("pDataList"), py::arg("pBackgroundList"));

    m.def("compute_viewpoint", [](py::object pViewpointObj, py::object pArgs, py::object pQueue,
                                  py::object pReferencePoints, py::object pGeneList, py::object pMatrix,
                                  py::object pBackgroundModel, py::object pBackgroundModelRelativeInteractions) {
        py::list interaction_data_list;
        try {
            py::ssize_t i = 0;
            for (const py::handle referencePoint : pReferencePoints) {
                const py::object fixate = pArgs.attr("fixateRange");
                const std::vector<py::object> fixed = unpack(
                    pViewpointObj.attr("calculateViewpointRange")(referencePoint, py::make_tuple(fixate, fixate)), 3);
                const py::object chromosome = getitem(referencePoint, 0);
                const std::vector<py::object> intermediate =
                    unpack(pViewpointObj.attr("computeViewpoint")(referencePoint, chromosome, fixed[0], fixed[1]), 2);
                const py::object denominator = array_sum(intermediate[0]);

                const std::vector<py::object> region =
                    unpack(pViewpointObj.attr("calculateViewpointRange")(referencePoint, pArgs.attr("range")), 3);
                const std::vector<py::object> computed =
                    unpack(pViewpointObj.attr("computeViewpoint")(referencePoint, chromosome, region[0], region[1]), 2);
                py::object data_list = computed[0];
                const py::object index_reference_point = computed[1];

                py::object background = flatten(
                    pViewpointObj.attr("interactionBackgroundData")(pBackgroundModelRelativeInteractions, region[2]));
                if (length(data_list) != length(background)) {
                    const py::tuple adjusted = adjust_viewpoint_data(pViewpointObj, data_list, background, referencePoint,
                                                                     region[0], region[1]);
                    data_list = adjusted[0];
                    background = adjusted[1];
                }
                const py::object average = pArgs.attr("averageContactBin");
                if (gt(average, py::int_(0)) && ge(py::int_(length(data_list)), average)) {
                    data_list = pViewpointObj.attr("smoothInteractionValues")(data_list, average);
                }
                const py::object data_list_raw = array_copy(data_list);
                data_list = pViewpointObj.attr("computeRelativeValues")(data_list, denominator);
                const py::object x_fold_list = x_fold(data_list, background);
                const py::object p_value_list =
                    pViewpointObj.attr("pvalues")(pBackgroundModel, data_list_raw, index_reference_point);

                const py::object range = pArgs.attr("range");
                const std::vector<py::object> region_range = unpack(pViewpointObj.attr("calculateViewpointRange")(
                    referencePoint, py::make_tuple(getitem(range, 0), getitem(range, 1))), 3);
                const py::object gene = getitem(pGeneList, i);
                const py::object interaction_data = pViewpointObj.attr("createInteractionFileDataHDF5")(
                    referencePoint, chromosome, region_range[0], region_range[1], data_list, data_list_raw, gene,
                    denominator, p_value_list, x_fold_list);

                py::list parts;
                for (const py::handle part : referencePoint) {
                    parts.append(str_object(part));
                }
                const py::object reference_point_string = join("_", parts);
                const py::str region_start_in_units = to_str(hicx::utilities::in_units(to_double(region[0])));
                const py::str region_end_in_units = to_str(hicx::utilities::in_units(to_double(region[1])));
                const py::object denominator_text =
                    add(py::str("Sum of interactions in fixate range: "), str_object(denominator));
                py::object header = add(py::str("# Interaction file, created with HiCExplorer's chicViewpoint version "),
                                        package_version());
                header = add(header, py::str("\n# "));
                py::list fields;
                fields.append(pMatrix);
                fields.append(reference_point_string);
                fields.append(str_object(region_start_in_units));
                fields.append(str_object(region_end_in_units));
                fields.append(gene);
                fields.append(denominator_text);
                header = add(header, join("\t", fields));
                header = add(header, py::str("\n# Chromosome\tStart\tEnd\tGene\tSum of interactions\tRelative "
                                             "position\tRelative Interactions\tp-value\tx-fold\tRaw\n#"));
                // '.'.join(pMatrix.split('/')[-1].split('.')[:-1])
                const py::object base_name = getitem(pMatrix.attr("split")("/"), -1);
                const py::object stem =
                    join(".", getitem(base_name.attr("split")("."), py::slice(py::none(), py::int_(-1), py::none())));
                py::list name_parts;
                name_parts.append(stem);
                name_parts.append(reference_point_string);
                name_parts.append(gene);
                const py::object matrix_name = join("_", name_parts);
                py::list gene_parts;
                gene_parts.append(reference_point_string);
                gene_parts.append(gene);
                const py::object reference_point_gene = str_object(join("_", gene_parts));
                py::list entry;
                entry.append(matrix_name);
                entry.append(interaction_data);
                entry.append(header);
                entry.append(reference_point_gene);
                interaction_data_list.append(entry);
                ++i;
            }
        } catch (py::error_already_set& error) {
            if (!is_exception(error)) {
                throw;
            }
            const std::string message = str_of(error.value());
            log(kLogger, "debug", "Error! " + message);
            pQueue.attr("put")(to_str("Fail: " + message));
            return;
        } catch (const std::exception& error) {
            const std::string message = fail_text(error);
            log(kLogger, "debug", "Error! " + message);
            pQueue.attr("put")(to_str("Fail: " + message));
            return;
        }
        pQueue.attr("put")(interaction_data_list);
    }, py::arg("pViewpointObj"), py::arg("pArgs"), py::arg("pQueue"), py::arg("pReferencePoints"),
       py::arg("pGeneList"), py::arg("pMatrix"), py::arg("pBackgroundModel"),
       py::arg("pBackgroundModelRelativeInteractions"));
}

}  // namespace hicx::python
