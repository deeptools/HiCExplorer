// hicexplorer._hicx.hicPlotDistVsCounts: compute_distance_mean and
// from_bed_to_cut_interval of hicexplorer/hicPlotDistVsCounts.py.
//
// compute_distance_mean replaces hicmat.matrix with its upper triangle, cut
// at int(maxdepth * 1.5 / binsize) diagonals with maxdepth, as 3.7.6 does
// with triu(k=0) - triu(k=max_depth_in_bins) and eliminate_zeros (so an
// entry beyond the cut survives as NaN when it is NaN or infinite, since
// v - v is NaN for those), and reduces it per unit with hicx::dist_counts
// (cpp/core/include/hicx/dist_vs_counts.hpp), which the tool runs too. The
// result keeps the 3.7.6 types: {unit: OrderedDict((distance, (mean, nb)))}
// with numpy.float64 means and numpy.int64 bin counts, (nan, 0) for an empty
// diagonal, and the same log records.
//
// from_bed_to_cut_interval reads the lines of `fh` with str.strip().split()'s
// rules and finds the overlapped bins in the chromosome's bins as the
// interval tree does.
//
// Deviations: duplicate entries of hicmat.matrix are summed in storage
// order; when a BED region overlaps several bins that an earlier region
// overlapped already, the message names the first of them in bin order
// (3.7.6 the first in the interval tree's set order); `fh` must yield str.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/dist_vs_counts.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/py_text.hpp"
#include "hicx/utilities.hpp"

namespace hicx::python {

namespace {

namespace dc = hicx::dist_counts;

constexpr const char* kLogger = "hicexplorer.hicPlotDistVsCounts";

struct Triangle {
    hicx::CsrMatrix matrix;
    std::vector<dc::Entry> entries;  // row major, columns sorted
};

// triu(m, k=0) [- triu(m, k=cut)], zeros eliminated, as a canonical csr.
Triangle upper_triangle(py::handle matrix_object, std::optional<std::int64_t> cut) {
    const Sparse matrix = sparse(matrix_object, "compute_distance_mean");
    const hicx::utilities::SparseEntries coordinates = matrix.entries();
    std::vector<double> data;
    dispatch_numeric(matrix.data, "compute_distance_mean", [&]<class T>() {
        if constexpr (std::is_arithmetic_v<T>) {
            const std::vector<T> values = to_vector<T>(matrix.data);
            data.assign(values.begin(), values.end());
        } else {
            raise_unsupported_dtype(matrix.data, "compute_distance_mean");
        }
    });
    std::vector<std::map<std::int64_t, double>> rows(static_cast<std::size_t>(matrix.rows));
    for (std::size_t k = 0; k < data.size(); ++k) {
        const std::int64_t r = coordinates.row[k];
        const std::int64_t c = coordinates.col[k];
        if (c < r) {
            continue;
        }
        auto [it, inserted] = rows[static_cast<std::size_t>(r)].try_emplace(c, data[k]);
        if (!inserted) {
            it->second += data[k];
        }
    }
    Triangle out;
    std::vector<std::int64_t> indptr{0};
    std::vector<std::int32_t> indices;
    std::vector<double> values;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (const auto& [c, v] : rows[r]) {
            double value = v;
            if (cut && c - static_cast<std::int64_t>(r) >= *cut) {
                value = v - v;
            }
            if (value != 0.0) {
                indices.push_back(static_cast<std::int32_t>(c));
                values.push_back(value);
                out.entries.push_back({static_cast<std::int32_t>(r), static_cast<std::int32_t>(c), value});
            }
        }
        indptr.push_back(static_cast<std::int64_t>(indices.size()));
    }
    out.matrix = hicx::CsrMatrix(matrix.rows, matrix.cols, std::move(indptr), std::move(indices), std::move(values),
                                 dtype_name(matrix.data));
    return out;
}

std::string format_float(double value) {
    return hicx::npy::float_repr(value);
}

// The fields of a cut interval tuple as str.
std::vector<hicx::CutInterval> intervals_of(py::handle list) {
    return cut_intervals(list);
}

// type(like)(value) for the integer fields of a cut interval
py::object int_like(std::int64_t value, py::handle like) {
    if (PyLong_Check(like.ptr())) {
        return py::int_(value);
    }
    return py::type::of(like)(value);
}

}  // namespace

void bind_hicPlotDistVsCounts(py::module_& m) {
    m.def("compute_distance_mean", [](py::handle hicmat, py::handle maxdepth, py::handle perchr,
                                      py::handle custom_cut_intervals) {
        const py::object binsize_object = hicobj::bin_size(hicmat);
        const std::int64_t binsize = to_int64(binsize_object);
        const bool maxdepth_set = hicobj::py_truthy(maxdepth);
        std::optional<std::int64_t> cut;
        if (maxdepth_set) {
            if (hicobj::py_less(maxdepth, binsize_object)) {
                raise_error("SystemExit",
                            "Please specify a maxDepth larger than bin size (" + str_of(binsize_object) + ")");
            }
            cut = static_cast<std::int64_t>(to_double(maxdepth) * 1.5 / static_cast<double>(binsize));
        }
        Triangle triangle = upper_triangle(hicmat.attr("matrix"), cut);
        hicmat.attr("matrix") = csr_matrix_object(triangle.matrix);

        py::object genome_wide_object;
        if (custom_cut_intervals.is_none()) {
            genome_wide_object = hicmat.attr("cut_intervals");
        } else {
            if (hicobj::length(custom_cut_intervals) != hicobj::length(hicmat.attr("cut_intervals"))) {
                raise_error("AssertionError", "The custom_cut_intervals should have the same size as hicmat.cut_intervals");
            }
            genome_wide_object = py::reinterpret_borrow<py::object>(custom_cut_intervals);
        }
        const std::vector<hicx::CutInterval> genome_wide = intervals_of(genome_wide_object);

        struct Unit {
            py::object name;
            std::size_t begin = 0;
            std::size_t end = 0;
            std::int32_t offset = 0;
            std::int64_t bins = 0;
            std::vector<hicx::CutInterval> intervals;
        };
        std::vector<Unit> units;
        const std::vector<dc::Entry>& entries = triangle.entries;
        const std::int64_t n = triangle.matrix.rows();
        if (hicobj::py_truthy(perchr)) {
            for (const py::handle chrname : hicobj::chr_names(hicmat)) {
                const auto [first, last] = hicobj::range_values(hicobj::chr_bin_range(hicmat, chrname));
                Unit unit;
                unit.name = py::reinterpret_borrow<py::object>(chrname);
                const std::int64_t a = std::clamp<std::int64_t>(first, 0, n);
                const std::int64_t b = std::clamp<std::int64_t>(last, a, n);
                unit.offset = static_cast<std::int32_t>(a);
                unit.bins = b - a;
                for (std::int64_t x = first; x < last; ++x) {
                    if (x < 0 || x >= static_cast<std::int64_t>(genome_wide.size())) {
                        raise_error("IndexError", "list index out of range");
                    }
                    unit.intervals.push_back(genome_wide[static_cast<std::size_t>(x)]);
                }
                units.push_back(std::move(unit));
                units.back().begin = 0;
                units.back().end = 0;
            }
        } else {
            Unit unit;
            unit.name = py::str("all");
            unit.begin = 0;
            unit.end = entries.size();
            unit.bins = n;
            unit.intervals = genome_wide;
            units.push_back(std::move(unit));
        }

        py::dict mean_dict;
        const py::object ordered_dict = py::module_::import("collections").attr("OrderedDict");
        const double depth = maxdepth.is_none() ? std::numeric_limits<double>::infinity() : to_double(maxdepth);
        for (Unit& unit : units) {
            std::vector<dc::Entry> block;
            const std::vector<dc::Entry>* source = &entries;
            if (hicobj::py_truthy(perchr)) {
                const std::int32_t a = unit.offset;
                const std::int64_t b = unit.offset + unit.bins;
                for (const dc::Entry& e : entries) {
                    if (e.row >= a && e.row < b && e.col >= a && e.col < b) {
                        block.push_back(e);
                    }
                }
                source = &block;
                unit.begin = 0;
                unit.end = block.size();
            }
            log(kLogger, "info", "processing chromosome " + str_of(unit.name) + "\n");
            bool snapped = false;
            std::int64_t median = 0;
            std::vector<dc::Distance> distances;
            try {
                if (unit.intervals.size() > 1) {
                    (void)dc::fitted_starts(unit.intervals, &snapped, &median);
                }
                if (snapped) {
                    char text[64];
                    std::snprintf(text, sizeof(text), "%f", static_cast<double>(median));
                    log(hicobj::kHicMatrixLogger, "info",
                        std::string("[getCountsByDistance] Bin size is not homogeneous, setting \n"
                                    "the bin distance to the median: ") +
                            text + "\n");
                }
                distances = dc::scan(*source, unit.begin, unit.end, unit.offset, unit.bins, unit.intervals, binsize,
                                     maxdepth_set);
            } catch (const hicx::PythonError&) {
                throw;
            }
            py::list items;
            for (const dc::Distance& d : distances) {
                if (d.skipped) {
                    continue;
                }
                const auto k = static_cast<std::int64_t>(d.k);
                if (d.defined && d.zero_value) {
                    const std::string length = d.diagonal_is_float
                                                   ? format_float(d.diagonal_length)
                                                   : std::to_string(static_cast<std::int64_t>(d.diagonal_length));
                    log(kLogger, "info", "zero value for " + std::to_string(k) + ", diagonal len: " + length + "\n");
                }
                if (d.stop) {
                    log(kLogger, "info",
                        "skipping rest of chromosome " + str_of(unit.name) + ". Too many emtpy diagonals\n");
                } else if (d.nan_sum) {
                    log(kLogger, "info", "nan value found for distance " + std::to_string((k - 1) * binsize) + "\n");
                }
                if (k > 0 && static_cast<double>((k - 1) * binsize) <= depth) {
                    py::object value;
                    if (d.defined) {
                        value = py::make_tuple(numpy_scalar<double>(d.mean),
                                               numpy_scalar<std::int64_t>(static_cast<std::int64_t>(d.diagonal_length)));
                    } else {
                        value = py::make_tuple(py::float_(std::nan("")), py::int_(0));
                    }
                    items.append(py::make_tuple(py::int_((k - 1) * binsize), value));
                }
            }
            mean_dict[unit.name] = ordered_dict(items);
        }
        return mean_dict;
    }, py::arg("hicmat"), py::arg("maxdepth") = py::none(), py::arg("perchr") = false,
       py::arg("custom_cut_intervals") = py::none());

    m.def("from_bed_to_cut_interval", [](py::handle hicmat, py::handle fh) {
        const py::object original = hicmat.attr("cut_intervals");
        const std::size_t count = hicobj::length(original);
        std::vector<py::object> fresh(count);
        const py::list chrom_list = hicobj::chr_names(hicmat);
        std::int64_t id = 0;
        for (const py::handle line_object : fh) {
            if (!is_str(line_object)) {
                raise_error("TypeError", "from_bed_to_cut_interval: HiCExplorer 4 reads str lines");
            }
            const std::string line = utf8(line_object);
            if (line.empty()) {
                raise_error("IndexError", "string index out of range");
            }
            if (line[0] == '#') {
                continue;
            }
            std::vector<std::string> fields;
            for (const std::string_view field : hicx::pytext::split_whitespace(hicx::pytext::strip(line))) {
                fields.emplace_back(field);
            }
            if (fields.empty()) {
                raise_error("IndexError", "list index out of range");
            }
            if (!hicobj::py_contains(chrom_list, to_str(fields[0]))) {
                const std::string changed = hicx::utilities::change_chrom_name(fields[0]);
                if (hicobj::py_contains(chrom_list, to_str(changed))) {
                    fields[0] = changed;
                } else {
                    continue;
                }
            }
            if (fields.size() < 3) {
                raise_error("ValueError",
                            "not enough values to unpack (expected 3, got " + std::to_string(fields.size()) + ")");
            }
            const py::str chrom = to_str(fields[0]);
            const std::int64_t start = to_int64(parse_int(to_str(fields[1])));
            const std::int64_t end = to_int64(parse_int(to_str(fields[2])));
            const py::object trees = hicmat.attr("interval_trees");
            if (!hicobj::py_contains(trees, chrom)) {
                raise_error("KeyError", fields[0]);
            }
            std::vector<std::int64_t> overlapping;
            if (start < end) {
                for (const hicobj::TreeBin& bin : hicobj::tree_bins(hicmat, chrom)) {
                    if (bin.start < end && bin.end > start && bin.start < bin.end) {
                        overlapping.push_back(bin.id);
                    }
                }
            }
            if (overlapping.empty()) {
                raise_error("ValueError", "min() iterable argument is empty");
            }
            const std::int64_t first = *std::min_element(overlapping.begin(), overlapping.end());
            const std::int64_t start_pos = hicobj::int64_of(hicobj::item(hicobj::item(original, first), 1));
            for (const std::int64_t i : overlapping) {
                const py::object interval = hicobj::item(original, i);
                if (fresh[static_cast<std::size_t>(i)]) {
                    raise_error("Exception", "2 features must not overlap the same bin." + str_of(interval) +
                                                 " is overlapped twice");
                }
                const py::object s = hicobj::item(interval, 1);
                const py::object e = hicobj::item(interval, 2);
                fresh[static_cast<std::size_t>(i)] =
                    py::make_tuple(py::str(std::to_string(id)), int_like(hicobj::int64_of(s) - start_pos, s),
                                   int_like(hicobj::int64_of(e) - start_pos, e), hicobj::item(interval, 3));
            }
            ++id;
        }
        if (id <= 0) {
            raise_error("AssertionError", "No region overlapped with bins.");
        }
        py::list out;
        for (std::size_t i = 0; i < count; ++i) {
            if (fresh[i]) {
                out.append(fresh[i]);
                continue;
            }
            const py::object interval = hicobj::item(original, static_cast<std::int64_t>(i));
            const py::object s = hicobj::item(interval, 1);
            const py::object e = hicobj::item(interval, 2);
            out.append(py::make_tuple(py::str("_ignore_" + std::to_string(i)), py::int_(0),
                                      int_like(hicobj::int64_of(e) - hicobj::int64_of(s), e), hicobj::item(interval, 3)));
        }
        return out;
    }, py::arg("hicmat"), py::arg("fh"));
}

}  // namespace hicx::python
