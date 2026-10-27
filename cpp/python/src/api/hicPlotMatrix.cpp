// hicexplorer._hicx.hicPlotMatrix: the data computations of
// hicexplorer/hicPlotMatrix.py.
//
// Public functions answered here entirely: relabel_ticks, translate_region,
// getRegion, chromosome_start_end and make_start_pos_array. The drawing
// functions plotPerChr and plotBigwig stay matplotlib code in the Python
// module (cpp/scripts/public_api_allowlist.tsv) and get their data from the
// helpers below, whose names start with an underscore:
//  * _plotPerChr_chromosomes, _plotPerChr_unpack and _plotPerChr_matrix: the
//    chromosome list, the unpacking of the cut intervals, and the dense
//    chromosome block with the --log/--log1p clearing and its log
//    records (hicx::plot_matrix::clear_for_log);
//  * _bigwig_check and _bigwig_track: the file name check and, per file, the
//    scores read with libBigWig as pyBigWig reads them (hicx::bigwig::Track),
//    the x positions (np.linspace), the axis limits, the sign flip, scaling
//    and clipping, and every log record of plotBigwig in its order.
//
// chromosome names go through hicexplorer.utilities.check_chrom_str_bytes and
// change_chrom_names (both on C++) as in 3.7.6, the matrix queries through
// hicmatrix_queries.hpp.
//
// Deviations: translate_region and the bigwig file names take str (3.7.6
// raised its own TypeError or AttributeError for other types, which is kept
// for bytes and non-str objects); a negative number of bigwig bins raises
// ValueError (pyBigWig reads it as an unsigned long and does not return).

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/bigwig_track.hpp"
#include "hicx/compartmentalization.hpp"
#include "hicx/plot_matrix_data.hpp"
#include "utilities.hpp"

namespace hicx::python {

namespace {

namespace pm = hicx::plot_matrix;

constexpr const char* kLogger = "hicexplorer.hicPlotMatrix";

py::object utilities_function(const char* name) {
    return py::module_::import("hicexplorer.utilities").attr(name);
}

[[noreturn]] void system_exit(const std::string& message) {
    raise_error("SystemExit", message);
}

void require_str(py::handle value, const char* method) {
    if (is_str(value)) {
        return;
    }
    if (PyBytes_Check(value.ptr())) {
        raise_error("TypeError", "a bytes-like object is required, not 'str'");
    }
    raise_error("AttributeError", "'" + type_name(value) + "' object has no attribute '" + method + "'");
}

py::object translate_region(py::handle region_string, py::handle ma) {
    require_str(region_string, "replace");
    const std::vector<std::string> fields = pm::region_fields(utf8(region_string));
    const py::object check = utilities_function("check_chrom_str_bytes");
    const py::object change = utilities_function("change_chrom_names");
    const py::object trees = ma.attr("interval_trees");
    py::object chrom = check(trees, to_str(fields[0]));
    if (!hicobj::py_contains(py::list(trees), chrom)) {
        chrom = check(trees, change(chrom));
        if (!hicobj::py_contains(py::list(trees), chrom)) {
            system_exit("Chromosome name " + str_of(change(chrom)) + " in --region not in matrix");
        }
    }
    const py::object range = hicobj::chr_bin_range(ma, chrom);
    if (hicobj::length(range) != 2) {
        raise_error("ValueError", "too many values to unpack (expected 2)");
    }
    const py::object first_bin = hicobj::item(range, 0);
    const py::object last_bin = hicobj::item(range, 1);
    const py::object region_start = fields.size() > 1 ? parse_int(to_str(fields[1]))
                                                       : hicobj::item(hicobj::bin_pos(ma, first_bin), 1);
    const py::object region_end =
        fields.size() > 2 ? parse_int(to_str(fields[2]))
                          : hicobj::item(hicobj::bin_pos(ma, py::int_(to_int64(last_bin) - 1)), 2);
    return py::make_tuple(chrom, region_start, region_end);
}

// The bins of `chrom` getRegion keeps, as (idx, start) lists.
struct Selection {
    py::list idx;
    py::list start;
};

enum class Rule { CoolRegion, H5Region, CoolRegion2, H5Region2 };

Selection select_bins(py::handle ma, py::handle chrom, std::int64_t rs, std::int64_t re, Rule rule) {
    Selection out;
    const py::object intervals = ma.attr("cut_intervals");
    const bool text = is_str(chrom);
    const std::string name = text ? utf8(chrom) : std::string();
    const std::size_t n = hicobj::length(intervals);
    for (std::size_t idx = 0; idx < n; ++idx) {
        const py::object x = hicobj::item(intervals, static_cast<std::int64_t>(idx));
        const py::object x0 = hicobj::item(x, 0);
        const bool same = (text && is_str(x0)) ? utf8(x0) == name : hicobj::py_equal(x0, chrom);
        if (!same) {
            continue;
        }
        const py::object start_object = hicobj::item(x, 1);
        const std::int64_t s = to_int64(start_object);
        const std::int64_t e = to_int64(hicobj::item(x, 2));
        bool keep = false;
        switch (rule) {
            case Rule::CoolRegion:
                keep = (s >= rs && e <= re) || (s < re && e < re && e > rs) || (s > rs && s < re);
                break;
            case Rule::H5Region:
                keep = s >= rs && e <= re;
                break;
            case Rule::CoolRegion2:
                keep = (s >= rs && e < re) || (s < re && e < re && e > rs) || (s > rs && s < re);
                break;
            case Rule::H5Region2:
                keep = s >= rs && e < re;
                break;
        }
        if (keep) {
            out.idx.append(py::int_(idx));
            out.start.append(start_object);
        }
    }
    if (out.idx.empty()) {
        raise_error("ValueError", "not enough values to unpack (expected 2, got 0)");
    }
    return out;
}

py::object chromosome_start_end(py::handle hic) {
    py::dict out;
    const py::object intervals = hic.attr("cut_intervals");
    const py::object boundaries = hic.attr("chrBinBoundaries");
    for (const py::handle chromosome : boundaries) {
        const py::object range = hicobj::item(boundaries, chromosome);
        if (hicobj::length(range) != 2) {
            raise_error("ValueError", "too many values to unpack (expected 2)");
        }
        const std::int64_t start_bin = to_int64(hicobj::item(range, 0));
        const std::int64_t end_bin = to_int64(hicobj::item(range, 1));
        out[chromosome] = py::make_tuple(hicobj::item(hicobj::item(intervals, start_bin), 1),
                                         hicobj::item(hicobj::item(intervals, end_bin - 1), 2));
    }
    return out;
}

bool is_numpy_integer(py::handle value) {
    return !PyLong_Check(value.ptr()) && is_integer(value);
}

py::array dense_block(py::handle hic, std::int64_t first, std::int64_t last) {
    const Sparse matrix = sparse(hic.attr("matrix"), "plotPerChr");
    const std::int64_t a = std::clamp<std::int64_t>(first < 0 ? first + matrix.rows : first, 0, matrix.rows);
    const std::int64_t b = std::clamp<std::int64_t>(last < 0 ? last + matrix.rows : last, a, matrix.rows);
    const std::int64_t ca = std::clamp<std::int64_t>(first < 0 ? first + matrix.cols : first, 0, matrix.cols);
    const std::int64_t cb = std::clamp<std::int64_t>(last < 0 ? last + matrix.cols : last, ca, matrix.cols);
    const auto rows = static_cast<std::size_t>(b - a);
    const auto cols = static_cast<std::size_t>(cb - ca);
    std::vector<double> dense(rows * cols, 0.0);
    const hicx::utilities::SparseEntries entries = matrix.entries();
    std::vector<double> data;
    dispatch_numeric(matrix.data, "plotPerChr", [&]<class T>() {
        if constexpr (std::is_arithmetic_v<T>) {
            const std::vector<T> values = to_vector<T>(matrix.data);
            data.assign(values.begin(), values.end());
        } else {
            raise_unsupported_dtype(matrix.data, "plotPerChr");
        }
    });
    for (std::size_t k = 0; k < data.size(); ++k) {
        const std::int64_t r = entries.row[k];
        const std::int64_t c = entries.col[k];
        if (r >= a && r < b && c >= ca && c < cb) {
            dense[static_cast<std::size_t>(r - a) * cols + static_cast<std::size_t>(c - ca)] += data[k];
        }
    }
    py::array_t<double> out({static_cast<py::ssize_t>(rows), static_cast<py::ssize_t>(cols)});
    if (!dense.empty()) {
        std::memcpy(out.mutable_data(), dense.data(), dense.size() * sizeof(double));
    }
    return out;
}

const char* truth(bool value) {
    return value ? "True" : "False";
}

// Python's floor division of two int64 values
std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        --q;
    }
    return q;
}

// int(bigwig_end - start) // pResolution
std::int64_t bin_count(std::int64_t span, py::handle resolution) {
    if (resolution.is_none() || !is_integer(resolution)) {
        raise_error("TypeError", "unsupported operand type(s) for //: 'int' and '" + type_name(resolution) + "'");
    }
    const std::int64_t r = to_int64(resolution);
    if (r == 0) {
        raise_error("ZeroDivisionError", "integer division or modulo by zero");
    }
    return floor_div(span, r);
}

py::array float_array(const std::vector<double>& values) {
    return from_vector(values);
}

std::string first_ten(const std::vector<double>& values) {
    const std::vector<double> head(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(10, values.size())));
    return str_of(float_array(head));
}

// The bigwig chromosome name as check_chrom_str_bytes(bw.chroms(), chrom)
// leaves it: bw.chroms() has str keys.
py::object bigwig_name(py::handle chrom, const hicx::bigwig::Track& track) {
    if (track.chroms().empty()) {
        PyErr_SetNone(PyExc_StopIteration);
        throw_python_error();
    }
    if (is_str(chrom)) {
        return py::reinterpret_borrow<py::object>(chrom);
    }
    return to_string(chrom);
}

bool track_has(const hicx::bigwig::Track& track, py::handle name) {
    return is_str(name) && track.length(utf8(name)).has_value();
}

std::vector<double> track_mean(const hicx::bigwig::Track& track, py::handle chrom, std::int64_t start,
                               std::int64_t end, std::int64_t bins) {
    if (bins < 0) {
        raise_error("ValueError", "plotBigwig: HiCExplorer 4 refuses a negative number of bigwig bins (" +
                                      std::to_string(bins) + ")");
    }
    try {
        return track.mean(utf8(chrom), start, end, bins);
    } catch (const std::invalid_argument& error) {
        raise_error("RuntimeError", error.what());
    } catch (const std::runtime_error& error) {
        raise_error("RuntimeError", error.what());
    }
}

std::vector<double> linspace(double start, double stop, std::int64_t num) {
    if (num < 0) {
        raise_error("ValueError", "Number of samples, " + std::to_string(num) + ", must be non-negative.");
    }
    return hicx::compartments::linspace(start, stop, num);
}

}  // namespace

void bind_hicPlotMatrix(py::module_& m) {
    m.def("relabel_ticks", [](py::handle pXTicks) {
        (void)hicobj::item(pXTicks, -1);  // the IndexError of pXTicks[-1] for an empty sequence
        std::vector<double> ticks;
        for (const py::handle tick : pXTicks) {
            ticks.push_back(to_double(tick));
        }
        std::vector<std::string> labels;
        try {
            labels = pm::relabel_ticks(ticks);
        } catch (const std::out_of_range&) {
            raise_error("IndexError", "list index out of range");
        }
        py::list out;
        for (const std::string& label : labels) {
            out.append(to_str(label));
        }
        return out;
    }, py::arg("pXTicks"));

    m.def("translate_region", [](py::handle region_string, py::handle ma) {
        return translate_region(region_string, ma);
    }, py::arg("region_string"), py::arg("ma"));

    m.def("getRegion", [](py::handle args, py::handle ma) {
        const py::object region = translate_region(args.attr("region"), ma);
        const py::object chrom = hicobj::item(region, 0);
        const py::object region_start = hicobj::item(region, 1);
        const py::object region_end = hicobj::item(region, 2);
        py::list as_list;
        as_list.append(chrom);
        as_list.append(region_start);
        as_list.append(region_end);
        args.attr("region") = as_list;
        const bool is_cooler = hicobj::py_truthy(utilities_function("check_cooler")(args.attr("matrix")));
        const Selection first = select_bins(ma, chrom, hicobj::int64_of(region_start), hicobj::int64_of(region_end),
                                            is_cooler ? Rule::CoolRegion : Rule::H5Region);
        py::object chrom2 = chrom;
        py::object region_start2 = py::none();
        py::object region_end2 = py::none();
        py::object idx2 = py::tuple(first.idx);
        py::object start_pos2 = py::tuple(first.start);
        if (PyObject_HasAttrString(args.ptr(), "region2") && hicobj::py_truthy(args.attr("region2"))) {
            const py::object second = translate_region(args.attr("region2"), ma);
            const py::object check = utilities_function("check_chrom_str_bytes");
            const py::object change = utilities_function("change_chrom_names");
            const py::object trees = ma.attr("interval_trees");
            chrom2 = check(trees, hicobj::item(second, 0));
            region_start2 = hicobj::item(second, 1);
            region_end2 = hicobj::item(second, 2);
            if (!hicobj::py_contains(py::list(trees), chrom2)) {
                chrom2 = check(trees, change(chrom2));
                if (!hicobj::py_contains(py::list(trees), chrom2)) {
                    system_exit("Chromosome name " + str_of(change(chrom2)) + " in --region2 not in matrix");
                }
            }
            const Selection bins2 = select_bins(ma, chrom2, hicobj::int64_of(region_start2),
                                                hicobj::int64_of(region_end2),
                                                is_cooler ? Rule::CoolRegion2 : Rule::H5Region2);
            idx2 = py::tuple(bins2.idx);
            start_pos2 = py::tuple(bins2.start);
        }
        return py::make_tuple(chrom, region_start, region_end, py::tuple(first.idx), py::tuple(first.start), chrom2,
                              region_start2, region_end2, idx2, start_pos2);
    }, py::arg("args"), py::arg("ma"));

    m.def("chromosome_start_end", [](py::handle pHiCMatrix) { return chromosome_start_end(pHiCMatrix); },
          py::arg("pHiCMatrix"));

    m.def("make_start_pos_array", [](py::handle ma) {
        const py::object sizes_source = chromosome_start_end(ma);
        py::dict chrom_sizes;
        std::vector<std::pair<std::int64_t, bool>> unused;
        for (const py::handle chromosome : sizes_source) {
            const py::object v = hicobj::item(sizes_source, chromosome);
            const py::object start = hicobj::item(v, 0);
            const py::object end = hicobj::item(v, 1);
            const bool numpy = is_numpy_integer(start) || is_numpy_integer(end);
            const std::int64_t size = to_int64(end) - to_int64(start);
            chrom_sizes[chromosome] = numpy ? numpy_scalar<std::int64_t>(size) : py::object(py::int_(size));
        }
        const py::object intervals = ma.attr("cut_intervals");
        py::object prev_chrom = hicobj::item(hicobj::item(intervals, 0), 0);
        std::int64_t prev_sum = 0;
        bool sum_numpy = false;
        std::int64_t shift = 0;
        bool shift_numpy = false;
        py::list start_pos;
        const std::size_t n = hicobj::length(intervals);
        for (std::size_t index = 0; index < n; ++index) {
            const py::object interval = hicobj::item(intervals, static_cast<std::int64_t>(index));
            if (hicobj::length(interval) != 4) {
                raise_error("ValueError", hicobj::length(interval) < 4
                                              ? "not enough values to unpack (expected 4, got " +
                                                    std::to_string(hicobj::length(interval)) + ")"
                                              : std::string("too many values to unpack (expected 4)"));
            }
            const py::object chrom = hicobj::item(interval, 0);
            const py::object start_object = hicobj::item(interval, 1);
            const std::int64_t start = to_int64(start_object);
            const bool start_numpy = is_numpy_integer(start_object);
            if (index == 0 && start != 0) {
                shift = start;
                shift_numpy = start_numpy;
            }
            if (!hicobj::py_equal(chrom, prev_chrom)) {
                const py::object size = hicobj::item(chrom_sizes, prev_chrom);
                prev_sum += to_int64(size);
                sum_numpy = sum_numpy || is_numpy_integer(size);
                prev_chrom = chrom;
                if (start != 0) {
                    shift = start;
                    shift_numpy = start_numpy;
                } else {
                    shift = 0;
                    shift_numpy = false;
                }
            }
            const std::int64_t value = start - shift + prev_sum;
            start_pos.append((start_numpy || shift_numpy || sum_numpy) ? numpy_scalar<std::int64_t>(value)
                                                                       : py::object(py::int_(value)));
        }
        return start_pos;
    }, py::arg("ma"));

    m.def("_plotPerChr_chromosomes", [](py::handle hic_matrix) { return hicobj::chr_names(hic_matrix); },
          py::arg("hic_matrix"));

    m.def("_plotPerChr_unpack", [](py::handle hic_matrix) {
        // chrom, start, end, _ = zip(*hic_matrix.cut_intervals)
        const py::object intervals = hic_matrix.attr("cut_intervals");
        std::size_t shortest = std::numeric_limits<std::size_t>::max();
        for (const py::handle interval : intervals) {
            shortest = std::min(shortest, hicobj::length(interval));
        }
        if (shortest == std::numeric_limits<std::size_t>::max()) {
            shortest = 0;
        }
        if (shortest < 4) {
            raise_error("ValueError", "not enough values to unpack (expected 4, got " + std::to_string(shortest) + ")");
        }
        if (shortest > 4) {
            raise_error("ValueError", "too many values to unpack (expected 4)");
        }
        return py::none();
    }, py::arg("hic_matrix"));

    m.def("_plotPerChr_matrix", [](py::handle hic_matrix, py::handle chrname, py::handle args) {
        const auto [first, last] = hicobj::range_values(hicobj::chr_bin_range(hic_matrix, chrname));
        py::array matrix = dense_block(hic_matrix, first, last);
        const bool log_scale = hicobj::py_truthy(args.attr("log"));
        const bool log1p = hicobj::py_truthy(args.attr("log1p"));
        std::vector<double> values = to_vector<double>(matrix);
        if (log_scale || log1p) {
            const pm::LogClearing clearing = pm::clear_for_log(values);
            log(kLogger, "debug", std::string("any nan ") + truth(clearing.any_nan_before));
            log(kLogger, "debug", std::string("any inf ") + truth(clearing.any_inf_before));
            if (clearing.failed) {
                log(kLogger, "debug", "Clearing of matrix failed.");
            }
            log(kLogger, "debug", std::string("any nanafter remove of nan: ") + truth(clearing.any_nan_after));
            log(kLogger, "debug", std::string("any inf after remove of inf: ") + truth(clearing.any_inf_after));
        }
        if (log1p) {
            for (double& v : values) {
                v += 1.0;
            }
        }
        write_back(matrix, values);
        py::object chr_bin_boundary = py::module_::import("collections").attr("OrderedDict")();
        chr_bin_boundary[chrname] = hicobj::item(hicobj::chromosome_sizes(hic_matrix), chrname);
        return py::make_tuple(matrix, chr_bin_boundary);
    }, py::arg("hic_matrix"), py::arg("chrname"), py::arg("args"));

    m.def("_bigwig_check", [](py::handle pNameOfBigwigList) {
        bool any = false;
        for (const py::handle file : pNameOfBigwigList) {
            require_str(file, "split");
            any = true;
            const std::string name = utf8(file);
            const std::size_t dot = name.rfind('.');
            const std::string format = dot == std::string::npos ? name : name.substr(dot + 1);
            if (format != "bigwig" && format != "bw") {
                log(kLogger, "error", "Given files are not bigwig");
                PyErr_SetObject(PyExc_SystemExit, py::int_(1).ptr());
                throw_python_error();
            }
        }
        if (!any) {
            raise_error("UnboundLocalError",
                        "cannot access local variable 'file_format' where it is not associated with a value");
        }
        return py::none();
    }, py::arg("pNameOfBigwigList"));

    m.def("_bigwig_track", [](py::handle bigwigFile, py::handle pChromosomeSizes, py::handle pRegion,
                              py::handle pVertical, py::handle pFlipBigwigSign, py::handle pScaleFactorBigwig,
                              py::handle pValueMin, py::handle pValueMax, py::handle pResolution) -> py::object {
        const bool vertical = hicobj::py_truthy(pVertical);
        std::optional<hicx::bigwig::Track> track;
        try {
            track.emplace(fs_path(bigwigFile));
        } catch (const std::runtime_error& error) {
            raise_error("RuntimeError", error.what());
        }
        std::vector<double> x_values;
        std::vector<double> scores;
        py::object limits = py::none();
        if (hicobj::py_truthy(pRegion)) {
            if (hicobj::length(pRegion) != 3) {
                raise_error("ValueError", "too many values to unpack (expected 3)");
            }
            py::object chrom = hicobj::item(pRegion, 0);
            const py::object region_start = hicobj::item(pRegion, 1);
            py::object region_end = hicobj::item(pRegion, 2);
            const py::object chrom_end = hicobj::item(hicobj::item(pChromosomeSizes, chrom), 1);
            if (hicobj::py_less(chrom_end, region_end)) {
                region_end = chrom_end;
            }
            chrom = bigwig_name(chrom, *track);
            if (!track_has(*track, chrom)) {
                chrom = utilities_function("change_chrom_names")(chrom);
                if (!track_has(*track, chrom)) {
                    log(kLogger, "info", "bigwig file has no chromosome named: " + str_of(chrom) + ".");
                    return py::none();
                }
            }
            const py::object length = py::int_(*track->length(utf8(chrom)));
            const py::object bigwig_end = hicobj::py_less(region_end, length) ? region_end : length;
            const std::int64_t start = hicobj::int64_of(region_start);
            const std::int64_t num_bins = bin_count(hicobj::int64_of(bigwig_end) - start, pResolution);
            log(kLogger, "debug", "chrom " + str_of(chrom) + ", region_start " + str_of(region_start) + ", bigwig_end " +
                                      str_of(bigwig_end) + ", num_bins " + std::to_string(num_bins));
            scores = track_mean(*track, chrom, start, hicobj::int64_of(bigwig_end), num_bins);
            x_values = linspace(to_double(region_start), to_double(region_end), num_bins);
            log(kLogger, "debug", "_x_vals[:10] " + first_ten(x_values));
            if (x_values.size() != scores.size()) {
                raise_error("AssertionError", "");
            }
            if (!vertical) {
                log(kLogger, "debug", "limits region: " + str_of(region_start) + ", " + str_of(region_end));
            }
            limits = py::make_tuple(region_start, region_end);
        } else if (hicobj::py_truthy(pChromosomeSizes)) {
            log(kLogger, "debug", "pChromosomeSizes: " + str_of(pChromosomeSizes));
            std::int64_t chrom_length_sum = 0;
            py::object min_start = py::none();
            for (const py::handle chrom : pChromosomeSizes) {
                py::object name = bigwig_name(chrom, *track);
                if (!track_has(*track, name)) {
                    name = py::str("chr") + name;
                    if (!track_has(*track, name)) {
                        log(kLogger, "info", "bigwig file as no chromosome named: " + str_of(chrom) + ".");
                        return py::none();
                    }
                }
                const py::object bounds = hicobj::item(pChromosomeSizes, chrom);
                const py::object first = hicobj::item(bounds, 0);
                const py::object last = hicobj::item(bounds, 1);
                const std::int64_t s = to_int64(first);
                const std::int64_t e = to_int64(last);
                const std::int64_t length = *track->length(utf8(name));
                const std::int64_t bigwig_end = length < e ? length : e;
                const std::int64_t num_bins = bin_count(bigwig_end - s, pResolution);
                log(kLogger, "debug", "chrom " + str_of(name) + ", region_start " + str_of(first) + ", bigwig_end " +
                                          str_of(last) + ", num_bins " + std::to_string(num_bins));
                const std::vector<double> part = track_mean(*track, name, s, bigwig_end, num_bins);
                const std::vector<double> x = linspace(static_cast<double>(chrom_length_sum + s),
                                                       static_cast<double>(chrom_length_sum + e - s), num_bins);
                log(kLogger, "debug", "_x_vals[:10] " + first_ten(x));
                if (x.size() != part.size()) {
                    raise_error("AssertionError", "");
                }
                x_values.insert(x_values.end(), x.begin(), x.end());
                scores.insert(scores.end(), part.begin(), part.end());
                log(kLogger, "debug", "len bigwig score " + std::to_string(scores.size()) + ", x_values " +
                                          std::to_string(x_values.size()));
                chrom_length_sum += e - s;
                if (min_start.is_none() || hicobj::py_less(first, min_start)) {
                    min_start = first;
                }
            }
            if (!vertical) {
                log(kLogger, "debug", "limits: " + str_of(min_start) + ", " + std::to_string(chrom_length_sum));
            }
            limits = py::make_tuple(min_start, py::int_(chrom_length_sum));
        }
        log(kLogger, "debug", "Number of data points: " + std::to_string(scores.size()));
        if (hicobj::py_truthy(pFlipBigwigSign)) {
            log(kLogger, "info", "Flipping sign of bigwig values.");
            for (double& v : scores) {
                v *= -1.0;
            }
        }
        if (!pScaleFactorBigwig.is_none() && !hicobj::py_equal(pScaleFactorBigwig, py::float_(1.0))) {
            log(kLogger, "info", "Scaling bigwig values.");
            const double factor = to_double(pScaleFactorBigwig);
            for (double& v : scores) {
                v *= factor;
            }
        }
        if (!pValueMin.is_none() || !pValueMax.is_none()) {
            const double low = pValueMin.is_none() ? 0.0 : to_double(pValueMin);
            const double high = pValueMax.is_none() ? 0.0 : to_double(pValueMax);
            pm::clip(scores, pValueMin.is_none() ? nullptr : &low, pValueMax.is_none() ? nullptr : &high);
        }
        log(kLogger, "debug", "Set values");
        if (vertical) {
            log(kLogger, "debug", "Set values vertical");
            std::reverse(scores.begin(), scores.end());
        } else {
            log(kLogger, "debug", "Set values horizontal");
        }
        return py::make_tuple(float_array(x_values), float_array(scores), limits);
    }, py::arg("bigwigFile"), py::arg("pChromosomeSizes"), py::arg("pRegion"), py::arg("pVertical"),
       py::arg("pFlipBigwigSign"), py::arg("pScaleFactorBigwig"), py::arg("pValueMin"), py::arg("pValueMax"),
       py::arg("pResolution"));
}

}  // namespace hicx::python
