// The queries of a hicmatrix.HiCMatrix.hiCMatrix object that the functions
// of hicAverageRegions, hicCompartmentalization, hicPCA, hicPlotDistVsCounts,
// hicPlotSVL, hicPlotMatrix and hicAggregateContacts make: getChrBinRange,
// getBinPos, getRegionBinRange, getBinSize, getChrNames,
// get_chromosome_sizes and the interval tree overlap query.
//
// They are answered in C++ from the object's state (cut_intervals,
// chrBinBoundaries, the keys of interval_trees, bin_size), with the results,
// result types, messages and side effects of hicmatrix 17 (the version the
// 3.7.6 reference runs on):
//  * getChrBinRange returns the stored tuple or raises its ValueError;
//  * getBinPos returns cut_intervals[i] (negative i wraps, as the list does);
//  * getRegionBinRange converts the chromosome name to the type of the
//    interval tree keys, raises its ValueError for an unknown chromosome,
//    returns (start bin, end bin) of the first bin (by start, end, id) that
//    contains each position, and returns None after logging "Index error"
//    on the hicmatrix.HiCMatrix logger when a position lies outside every
//    bin. The bins of a chromosome are cut_intervals[first:last] of its
//    chrBinBoundaries entry, which is exactly what its interval tree holds:
//    intervalListToIntervalTree resets a chromosome's tree and overwrites its
//    boundary entry together when the chromosome reappears;
//  * getBinSize returns the cached bin_size, or computes it (hicx::BinTable)
//    and caches it on the object, logging the hicmatrix warning once when
//    the bins are not homogeneous.
//
// Deviation, documented in the cases files: a position that int() cannot
// convert raises int()'s ValueError; hicmatrix logs it with a '%d' format
// that itself fails inside logging and then indexes the tree with the
// unconverted value.

#ifndef HICX_PYTHON_API_HICMATRIX_QUERIES_HPP
#define HICX_PYTHON_API_HICMATRIX_QUERIES_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "api.hpp"
#include "hicx/bins.hpp"
#include "utilities.hpp"

namespace hicx::python::hicobj {

inline constexpr const char* kHicMatrixLogger = "hicmatrix.HiCMatrix";

inline bool py_isinstance(py::handle object, py::handle type) {
    const int result = PyObject_IsInstance(object.ptr(), type.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool py_contains(py::handle container, py::handle item) {
    const int result = PySequence_Contains(container.ptr(), item.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool py_truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool py_less(py::handle a, py::handle b) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_LT);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

inline bool py_equal(py::handle a, py::handle b) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_EQ);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

// int(object), Python's own conversion.
inline py::object py_int(py::handle object) {
    PyObject* value = PyNumber_Long(object.ptr());
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

inline std::int64_t int64_of(py::handle object) {
    return to_int64(py_int(object));
}

// obj[index] with Python's indexing.
inline py::object item(py::handle object, py::handle index) {
    PyObject* value = PyObject_GetItem(object.ptr(), index.ptr());
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

inline py::object item(py::handle object, std::int64_t index) {
    return item(object, py::int_(index));
}

inline std::size_t length(py::handle object) {
    const Py_ssize_t size = PyObject_Length(object.ptr());
    if (size < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(size);
}

// hic.getChrBinRange(chrName): the stored value.
inline py::object chr_bin_range(py::handle hic, py::handle chrName) {
    const py::object boundaries = hic.attr("chrBinBoundaries");
    if (py_contains(boundaries, chrName)) {
        return item(boundaries, chrName);
    }
    raise_error("ValueError", "chrName: " + str_of(chrName) + " not found in chrBinBoundaries" +
                                  "valid chromosomes are: " + str_of(boundaries.attr("keys")()));
}

// The two ends of a getChrBinRange tuple.
inline std::pair<std::int64_t, std::int64_t> range_values(py::handle range) {
    return {to_int64(item(range, 0)), to_int64(item(range, 1))};
}

// hic.getBinPos(binIndex)
inline py::object bin_pos(py::handle hic, py::handle binIndex) {
    const py::object intervals = hic.attr("cut_intervals");
    if (py_less(binIndex, py::int_(length(intervals)))) {
        return item(intervals, binIndex);
    }
    raise_error("ValueError", "binIndex: " + str_of(binIndex) + " not found");
}

// hic.getChrNames()
inline py::list chr_names(py::handle hic) {
    return py::list(hic.attr("chrBinBoundaries"));
}

// The chromosome name as getRegionBinRange converts it to the key type of
// interval_trees; raises its ValueError when the tree does not exist.
inline py::object region_chrom(py::handle hic, py::handle chrname) {
    const py::object trees = hic.attr("interval_trees");
    py::object iterator = py::reinterpret_steal<py::object>(PyObject_GetIter(trees.ptr()));
    if (!iterator) {
        throw_python_error();
    }
    PyObject* first_raw = PyIter_Next(iterator.ptr());
    if (first_raw == nullptr) {
        if (!PyErr_Occurred()) {
            PyErr_SetNone(PyExc_StopIteration);
        }
        throw_python_error();
    }
    const py::object first = py::reinterpret_steal<py::object>(first_raw);
    py::object name = py::reinterpret_borrow<py::object>(chrname);
    if (!py_isinstance(first, py::type::of(name))) {
        if (is_str(first)) {
            name = to_string(name);
        } else if (PyBytes_Check(first.ptr()) ||
                   py_isinstance(first, py::module_::import("numpy").attr("bytes_"))) {
            name = to_bytes(name);
        }
    }
    if (!py_contains(trees, name)) {
        raise_error("ValueError", "chromosome: " + str_of(name) + " name not found in matrix" +
                                      "valid names are: " + str_of(trees.attr("keys")()));
    }
    return name;
}

// The bins of the interval tree of `chrom` (a key of interval_trees): the
// cut intervals of its chrBinBoundaries range, as (start, end, id).
struct TreeBin {
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::int64_t id = 0;
};

inline std::vector<TreeBin> tree_bins(py::handle hic, py::handle chrom) {
    std::vector<TreeBin> bins;
    const py::object boundaries = hic.attr("chrBinBoundaries");
    if (!py_contains(boundaries, chrom)) {
        return bins;
    }
    const auto [first, last] = range_values(item(boundaries, chrom));
    const py::object intervals = hic.attr("cut_intervals");
    bins.reserve(static_cast<std::size_t>(last > first ? last - first : 0));
    for (std::int64_t i = first; i < last; ++i) {
        const py::object interval = item(intervals, i);
        bins.push_back(TreeBin{int64_of(item(interval, 1)), int64_of(item(interval, 2)), i});
    }
    return bins;
}

// sorted(tree[position:position + 1])[0].data, nullopt when nothing overlaps.
inline std::optional<std::int64_t> first_overlap(const std::vector<TreeBin>& bins, std::int64_t begin,
                                                 std::int64_t end) {
    const TreeBin* best = nullptr;
    for (const TreeBin& bin : bins) {
        if (bin.start < end && bin.end > begin && bin.start < bin.end) {
            if (best == nullptr || bin.start < best->start ||
                (bin.start == best->start && (bin.end < best->end || (bin.end == best->end && bin.id < best->id)))) {
                best = &bin;
            }
        }
    }
    if (best == nullptr) {
        return std::nullopt;
    }
    return best->id;
}

// hic.getRegionBinRange(chrname, startpos, endpos) as values; nullopt where
// the Python returns None (after its log record).
inline std::optional<std::pair<std::int64_t, std::int64_t>> region_bin_range_values(
    py::handle hic, py::handle chrname, py::handle startpos, py::handle endpos) {
    const py::object chrom = region_chrom(hic, chrname);
    const std::int64_t start = int64_of(startpos);
    const std::int64_t end = int64_of(endpos);
    const std::vector<TreeBin> bins = tree_bins(hic, chrom);
    const std::optional<std::int64_t> start_bin = first_overlap(bins, start, start + 1);
    const std::optional<std::int64_t> end_bin =
        start_bin ? first_overlap(bins, end, end + 1) : std::nullopt;
    if (!start_bin || !end_bin) {
        log(kHicMatrixLogger, "error", "Index error");
        return std::nullopt;
    }
    return std::make_pair(*start_bin, *end_bin);
}

inline py::object region_bin_range(py::handle hic, py::handle chrname, py::handle startpos,
                                   py::handle endpos) {
    const auto range = region_bin_range_values(hic, chrname, startpos, endpos);
    if (!range) {
        return py::none();
    }
    return py::make_tuple(range->first, range->second);
}

// hic.getBinSize(), with its caching on the object.
inline py::object bin_size(py::handle hic) {
    const py::object cached = hic.attr("bin_size");
    if (!cached.is_none()) {
        return cached;
    }
    const py::object intervals = hic.attr("cut_intervals");
    const std::vector<hicx::CutInterval> bins = cut_intervals(intervals);
    if (bins.empty()) {
        raise_error("ValueError", "not enough values to unpack (expected 4, got 0)");
    }
    py::object value;
    if (bins.size() == 1) {
        // diff[0] of np.array(end) - np.array(start)
        value = numpy_scalar<std::int64_t>(bins[0].end - bins[0].start);
    } else {
        hicx::BinTable table(bins);
        std::int64_t size = 0;
        try {
            size = table.bin_size();
        } catch (const std::runtime_error&) {
            raise_error("ValueError", "need at least one array to concatenate");
        }
        if (!table.bin_size_homogeneous()) {
            hic.attr("bin_size_homogeneous") = py::bool_(false);
            if (!py_truthy(hic.attr("non_homogeneous_warning_already_printed"))) {
                char text[64];
                std::snprintf(text, sizeof(text), "%f", static_cast<double>(size));
                log(kHicMatrixLogger, "warning",
                    std::string("Bin size is not homogeneous. ") + std::string(38, ' ') + "Median " +
                        text + "\n");
                hic.attr("non_homogeneous_warning_already_printed") = py::bool_(true);
            }
        }
        value = py::int_(size);
    }
    hic.attr("bin_size") = value;
    return value;
}

inline std::int64_t bin_size_value(py::handle hic) {
    return to_int64(bin_size(hic));
}

// hic.get_chromosome_sizes(): {chrom: end of its last bin}, an OrderedDict
// in chrBinBoundaries order.
inline py::object chromosome_sizes(py::handle hic) {
    const py::object boundaries = hic.attr("chrBinBoundaries");
    if (!py_truthy(boundaries)) {
        return py::none();
    }
    py::object out = py::module_::import("collections").attr("OrderedDict")();
    const py::object intervals = hic.attr("cut_intervals");
    for (const py::handle chrom : boundaries) {
        const auto [first, last] = range_values(item(boundaries, chrom));
        (void)first;
        const py::object interval = item(intervals, last - 1);
        out[item(interval, 0)] = item(interval, 2);
    }
    return out;
}

}  // namespace hicx::python::hicobj

#endif  // HICX_PYTHON_API_HICMATRIX_QUERIES_HPP
