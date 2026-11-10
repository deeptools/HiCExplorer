// hicexplorer._hicx.hicAverageRegions: calculateViewpointRange,
// getBinIndices and calculateViewpointRangeBins of
// hicexplorer/hicAverageRegions.py.
//
// The matrix queries are hicmatrix_queries.hpp, the windows
// hicx::average_regions (cpp/core/include/hicx/average_regions.hpp), which
// the tool hicAverageRegions runs too.
//
// Result types follow the Python arithmetic: a value computed from an
// integer of the range keeps the range's type (int for int, numpy.int64 for
// a numpy integer), a clamped value is the bound it was clamped to (0, the
// end of the chromosome's last bin minus one, or the chromosome's first or
// last bin). Ranges of other types (floats, unsigned numpy integers) raise
// TypeError naming HiCExplorer 4.

#include <string>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/average_regions.hpp"

namespace hicx::python {

namespace {

namespace ar = hicx::average_regions;

enum class IntKind { Python, Numpy };

IntKind int_kind(py::handle value, const char* function) {
    if (PyLong_Check(value.ptr())) {
        return IntKind::Python;
    }
    if (hicobj::py_isinstance(value, py::module_::import("numpy").attr("signedinteger"))) {
        return IntKind::Numpy;
    }
    raise_error("TypeError", std::string(function) + ": HiCExplorer 4 takes integer ranges, got " +
                                 type_name(value));
}

py::object make_int(std::int64_t value, IntKind kind) {
    if (kind == IntKind::Numpy) {
        return numpy_scalar<std::int64_t>(value);
    }
    return py::int_(value);
}

bool mapping_is(py::handle mapping, const char* name) {
    return hicobj::py_equal(mapping, py::str(name));
}

[[noreturn]] void unbound(const char* variable) {
    raise_error("UnboundLocalError", std::string("cannot access local variable '") + variable +
                                         "' where it is not associated with a value");
}

std::int64_t center(py::handle viewpoint) {
    return ar::center_of(to_double(hicobj::item(viewpoint, 1)), to_double(hicobj::item(viewpoint, 2)));
}

// getBinIndices(...)[position], with the TypeError of subscripting None.
py::object bin_index(py::handle region, int position) {
    if (region.is_none()) {
        raise_error("TypeError", "'NoneType' object is not subscriptable");
    }
    return hicobj::item(region, position);
}

py::object get_bin_indices(py::handle hic, py::handle viewpoint) {
    return hicobj::region_bin_range(hic, hicobj::item(viewpoint, 0), hicobj::item(viewpoint, 1),
                                    hicobj::item(viewpoint, 2));
}

}  // namespace

void bind_hicAverageRegions(py::module_& m) {
    m.def("calculateViewpointRange", [](py::handle pHiCMatrix, py::handle pViewpoint, py::handle pRange,
                                        py::handle pCoordinatesToBinMapping) {
        const py::object chrom_range = hicobj::chr_bin_range(pHiCMatrix, hicobj::item(pViewpoint, 0));
        const py::object last_index = py::int_(to_int64(hicobj::item(chrom_range, 1)) - 1);
        const py::object max_length_object = hicobj::item(hicobj::bin_pos(pHiCMatrix, last_index), 2);
        std::int64_t anchor = 0;
        if (mapping_is(pCoordinatesToBinMapping, "start")) {
            anchor = hicobj::int64_of(hicobj::item(pViewpoint, 1));
        } else if (mapping_is(pCoordinatesToBinMapping, "end")) {
            anchor = hicobj::int64_of(hicobj::item(pViewpoint, 2));
        } else if (mapping_is(pCoordinatesToBinMapping, "center")) {
            anchor = center(pViewpoint);
        } else {
            unbound("region_start");
        }
        const py::object up = hicobj::item(pRange, 0);
        const py::object down = hicobj::item(pRange, 1);
        const IntKind up_kind = int_kind(up, "calculateViewpointRange");
        const IntKind down_kind = int_kind(down, "calculateViewpointRange");
        const IntKind max_kind = int_kind(max_length_object, "calculateViewpointRange");
        const ar::Window window = ar::viewpoint_range(anchor, to_int64(up), to_int64(down), to_int64(max_length_object));
        const py::object start = window.start_out_of_range ? py::object(py::int_(0)) : make_int(window.start, up_kind);
        const py::object end = make_int(window.end, window.end_out_of_range ? max_kind : down_kind);
        return py::make_tuple(start, end, py::bool_(window.start_out_of_range), py::bool_(window.end_out_of_range));
    }, py::arg("pHiCMatrix"), py::arg("pViewpoint"), py::arg("pRange"), py::arg("pCoordinatesToBinMapping"));

    m.def("getBinIndices", [](py::handle pHiCMatrix, py::handle pViewpoint) {
        return get_bin_indices(pHiCMatrix, pViewpoint);
    }, py::arg("pHiCMatrix"), py::arg("pViewpoint"));

    m.def("calculateViewpointRangeBins", [](py::handle pHiCMatrix, py::handle pViewpoint, py::handle pRange,
                                            py::handle pCoordinatesToBinMapping) {
        py::object index;
        if (mapping_is(pCoordinatesToBinMapping, "start")) {
            index = bin_index(get_bin_indices(pHiCMatrix, pViewpoint), 0);
        } else if (mapping_is(pCoordinatesToBinMapping, "end")) {
            index = bin_index(get_bin_indices(pHiCMatrix, pViewpoint), 1);
        } else {
            const py::int_ middle(center(pViewpoint));
            const py::object region =
                hicobj::region_bin_range(pHiCMatrix, hicobj::item(pViewpoint, 0), middle, middle);
            index = bin_index(region, 1);
        }
        const py::object chrom_range = hicobj::chr_bin_range(pHiCMatrix, hicobj::item(pViewpoint, 0));
        const py::object first_bin = hicobj::item(chrom_range, 0);
        const py::object last_bin = hicobj::item(chrom_range, 1);
        if (hicobj::length(chrom_range) != 2) {
            raise_error("ValueError", "too many values to unpack (expected 2)");
        }
        const py::object up = hicobj::item(pRange, 0);
        const py::object down = hicobj::item(pRange, 1);
        const IntKind up_kind = int_kind(up, "calculateViewpointRangeBins");
        const IntKind down_kind = int_kind(down, "calculateViewpointRangeBins");
        const ar::Window window = ar::viewpoint_range_bins(to_int64(index), to_int64(up), to_int64(down),
                                                           to_int64(first_bin), to_int64(last_bin));
        const py::object start = window.start_out_of_range ? first_bin : make_int(window.start, up_kind);
        const py::object end = window.end_out_of_range ? last_bin : make_int(window.end, down_kind);
        return py::make_tuple(start, end, py::bool_(window.start_out_of_range), py::bool_(window.end_out_of_range));
    }, py::arg("pHiCMatrix"), py::arg("pViewpoint"), py::arg("pRange"), py::arg("pCoordinatesToBinMapping"));
}

}  // namespace hicx::python
