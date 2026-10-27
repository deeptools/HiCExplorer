// hicexplorer._hicx.chicPlotViewpoint: the data plot_images computes besides
// drawing (main() and parse_arguments() run the executable
// cpp/tools/chicPlotViewpoint.cpp).
//
// plot_images stays matplotlib (cpp/scripts/public_api_allowlist.tsv); the
// values it draws come from the Viewpoint object's C++ methods and from the
// two functions here: the p-values of the heat strip (zero p-values set to
// one, the significant regions' p-values, clipped to the colour range, in
// float32 as numpy held them) and the ticks and labels of the x axis.

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"

namespace hicx::python {

using namespace chic_api;

namespace {

// np.array(values, dtype=np.float32)
std::vector<float> as_float32(py::handle values) {
    const py::array array = py::reinterpret_borrow<py::array>(
        numpy().attr("asarray")(values, py::arg("dtype") = numpy().attr("float32")));
    return to_vector<float>(array);
}

// np.maximum / np.minimum with NaN propagation, in float32
float maximum(float a, float b) { return (std::isnan(a) || std::isnan(b)) ? std::nanf("") : (a > b ? a : b); }
float minimum(float a, float b) { return (std::isnan(a) || std::isnan(b)) ? std::nanf("") : (a < b ? a : b); }

void require_index(py::handle value, const char* name) {
    if (!is_integer(value)) {
        raise_error("TypeError", "'" + type_name(value) + "' object cannot be interpreted as an integer");
    }
    (void)name;
}

}  // namespace

void bind_chicPlotViewpoint(py::module_& m) {
    m.doc() = "C++ data of hicexplorer.chicPlotViewpoint.plot_images.";

    // if truncateZeroPvalues: p_values = float32; p_values[p_values == 0.0] = 1.0
    // if pValue and (min or max set): float32, significant p-values, clip
    m.def("plot_images_p_values", [](py::object p_values, py::handle truncate_zero, py::handle p_value,
                                     py::handle min_p_value, py::handle max_p_value,
                                     py::handle significant_p_values) -> py::object {
        bool converted = false;
        std::vector<float> values;
        if (truthy(truncate_zero)) {
            values = as_float32(p_values);
            for (float& value : values) {
                if (value == 0.0F) {
                    value = 1.0F;
                }
            }
            converted = true;
        }
        if (truthy(p_value) && (!min_p_value.is_none() || !max_p_value.is_none())) {
            if (!converted) {
                values = as_float32(p_values);
            } else {
                values = as_float32(from_vector(values));
            }
            converted = true;
            if (!significant_p_values.is_none() && truthy(significant_p_values)) {
                const std::int64_t n = static_cast<std::int64_t>(values.size());
                for (const py::handle location : significant_p_values) {
                    const py::object first = getitem(location, 0);
                    const py::object last = getitem(location, 1);
                    require_index(first, "start");
                    require_index(last, "stop");
                    const float replacement = static_cast<float>(to_double(getitem(location, 2)));
                    for (std::int64_t x = to_int64(first); x < to_int64(last); ++x) {
                        if (x < n) {
                            const std::int64_t position = x < 0 ? x + n : x;
                            if (position < 0) {
                                raise_error("IndexError", "index " + std::to_string(x) +
                                                              " is out of bounds for axis 0 with size " +
                                                              std::to_string(n));
                            }
                            values[static_cast<std::size_t>(position)] = replacement;
                        }
                    }
                }
            }
            const bool low = !min_p_value.is_none();
            const bool high = !max_p_value.is_none();
            const float lower = low ? static_cast<float>(to_double(min_p_value)) : 0.0F;
            const float upper = high ? static_cast<float>(to_double(max_p_value)) : 0.0F;
            for (float& value : values) {
                if (low) {
                    value = maximum(value, lower);
                }
                if (high) {
                    value = minimum(value, upper);
                }
            }
        }
        if (!converted) {
            return p_values;
        }
        return from_vector(values);
    }, py::arg("p_values"), py::arg("truncateZeroPvalues"), py::arg("pValue"), py::arg("minPValue"),
       py::arg("maxPValue"), py::arg("significant_p_values"));

    // The ticks and labels of the x axis for pArgs.range and the resolution.
    m.def("plot_images_ticks", [](py::handle pRange, py::handle pResolution) {
        const py::object upstream = getitem(pRange, 0);
        const py::object downstream = getitem(pRange, 1);
        const double total = to_double(add(upstream, downstream));
        double divisor_legend = 0.0;
        double mod_legend = 0.0;
        std::string unit;
        if (total <= 2e6) {
            divisor_legend = 1e3;
            mod_legend = 2e5;
            if (total <= 1e4) {
                mod_legend = 5e3;
            } else if (total <= 5e4) {
                mod_legend = 1e4;
            } else if (total <= 1e5) {
                mod_legend = 5e4;
            } else if (total <= 5e5) {
                mod_legend = 1e5;
            }
            unit = "kb";
        } else if (total > 2e6) {
            divisor_legend = 1e6;
            mod_legend = 1e6;
            unit = "Mb";
        } else {
            raise_error("UnboundLocalError",
                        "cannot access local variable 'mod_legend' where it is not associated with a value");
        }
        require_index(upstream, "upstream");
        require_index(downstream, "downstream");
        const std::int64_t r0 = to_int64(upstream);
        const std::int64_t r1 = to_int64(downstream);
        const std::int64_t divisor = static_cast<std::int64_t>(divisor_legend);
        // Python's floor division of integers
        const auto floor_div = [](std::int64_t a, std::int64_t b) {
            std::int64_t q = a / b;
            if ((a % b != 0) && ((a < 0) != (b < 0))) {
                --q;
            }
            return q;
        };
        py::list ticks;
        py::list x_labels;
        // zip(range(r0), range(r0, 1, -1)): k = 0 .. r0 - 2, j = r0 - k
        for (std::int64_t k = 0; k + 1 < r0; ++k) {
            const std::int64_t j = r0 - k;
            if (std::fmod(static_cast<double>(j), mod_legend) == 0.0) {
                x_labels.append(py::str(std::to_string(floor_div(-j, divisor)) + unit));
                ticks.append(floordiv(py::int_(k), pResolution));
            }
        }
        x_labels.append(py::str("RP"));
        ticks.append(floordiv(upstream, pResolution));
        const py::object referencepoint_index = ticks[ticks.size() - 1];
        // zip(range(r1), range(1, r1 + 1)): k = 0 .. r1 - 1, j = k + 1
        for (std::int64_t k = 0; k < r1; ++k) {
            const std::int64_t j = k + 1;
            if (std::fmod(static_cast<double>(j), mod_legend) == 0.0) {
                x_labels.append(py::str(std::to_string(floor_div(j, divisor)) + unit));
                ticks.append(add(referencepoint_index, floordiv(py::int_(k), pResolution)));
            }
        }
        return py::make_tuple(ticks, x_labels);
    }, py::arg("pRange"), py::arg("pResolution"));
}

}  // namespace hicx::python
