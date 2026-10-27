// hicexplorer._hicx.hicCorrectMatrix: the C++ bodies of correct_subparser,
// iterative_correction, fill_gaps, the methods of MAD, filter_by_zscore and
// the computations of plot_total_contact_dist.
//
// The MAD statistics, the row sums and the gap filling are
// hicx::correct_api (cpp/core/include/hicx/correct_api.hpp), in the dtype
// rules of numpy and scipy. The
// hiCMatrix is asked for its chromosome ranges through its own methods, MAD
// objects are built through the module's MAD class (whose methods are the
// functions here), and iterative_correction calls the module's
// iterativeCorrection (hicexplorer.iterativeCorrection, also C++).
// correct_subparser builds the argparse parser of 3.7.6 call by call.
//
// plot_total_contact_dist stays Python because it draws (allowlisted); the
// data it draws come from prepare_total_contact_dist,
// chromosome_contact_dist, genome_contact_dist, values_below and
// first_local_minimum below.
//
// Deviations: MAD takes int32, int64, uint32, uint64 and float64 points of one
// or two dimensions (float32 points, which numpy would keep in single
// precision, raise TypeError naming HiCExplorer 4), so filter_by_zscore and
// the contact distribution of a float32 matrix raise it too; fill_gaps reads
// the matrix with duplicate entries summed in storage order.

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/correct_api.hpp"
#include "hicx/numpy_compat.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace ca = hicx::correct_api;

namespace {

constexpr const char* kLogger = "hicexplorer.hicCorrectMatrix";

py::module_ this_module() {
    return py::module_::import("hicexplorer.hicCorrectMatrix");
}

[[noreturn]] void unsupported_points(const std::string& what) {
    raise_error("TypeError", "MAD: HiCExplorer 4 takes integer or float64 points, not " + what);
}

struct Points {
    std::vector<double> values;  // C order
    std::int64_t rows = 0;
    std::int64_t cols = 1;
};

Points points_of(py::handle object) {
    if (!is_ndarray(object)) {
        if (!py::hasattr(object, "shape")) {
            raise_error("AttributeError", "'" + type_name(object) + "' object has no attribute 'shape'");
        }
        unsupported_points("'" + type_name(object) + "'");
    }
    const py::array array = py::reinterpret_borrow<py::array>(object);
    Points p;
    switch (dtype_of(array)) {
        case DType::Int32:
        case DType::Int64:
        case DType::UInt32:
        case DType::UInt64:
        case DType::Float64:
            p.values = lc::doubles(array, "MAD");
            break;
        default:
            unsupported_points("dtype " + dtype_name(array));
    }
    if (array.ndim() == 1) {
        p.rows = array.shape(0);
    } else if (array.ndim() == 2) {
        p.rows = array.shape(0);
        p.cols = array.shape(1);
    } else {
        unsupported_points(std::to_string(array.ndim()) + " dimensional points");
    }
    return p;
}

// ---------------------------------------------------------------------------
// MAD

void mad_init(py::handle self, py::handle points_object) {
    self.attr("mad_b_value") = py::float_(ca::kMadB);
    const Points p = points_of(points_object);
    const ca::MadStats stats = ca::mad(p.values, p.rows, p.cols);
    self.attr("median") = numpy_scalar<double>(stats.median);
    self.attr("med_abs_deviation") = numpy_scalar<double>(stats.med_abs_deviation);
    self.attr("modified_z_score") = from_vector(stats.modified_z_score);
}

// A number or an array argument of value_to_mad / mad_to_value as float64,
// with its shape (empty for a scalar).
struct Values {
    std::vector<double> values;
    std::vector<py::ssize_t> shape;
    bool scalar = true;
};

Values values_of(py::handle object) {
    Values v;
    if (PyLong_Check(object.ptr()) || PyFloat_Check(object.ptr())) {
        v.values.push_back(to_double(object));
        return v;
    }
    py::array array = py::array::ensure(object);
    if (!array) {
        throw_python_error();
    }
    switch (dtype_of(array)) {
        case DType::Bool:
        case DType::Int8:
        case DType::Int16:
        case DType::Int32:
        case DType::Int64:
        case DType::UInt8:
        case DType::UInt16:
        case DType::UInt32:
        case DType::UInt64:
        case DType::Float64:
            v.values = lc::doubles(array, "MAD");
            break;
        default:
            unsupported_points("dtype " + dtype_name(array));
    }
    v.scalar = array.ndim() == 0;
    for (py::ssize_t d = 0; d < array.ndim(); ++d) {
        v.shape.push_back(array.shape(d));
    }
    return v;
}

py::object values_object(const Values& like, const std::vector<double>& values) {
    if (like.scalar) {
        return numpy_scalar<double>(values[0]);
    }
    return from_vector(values).attr("reshape")(py::tuple(py::cast(like.shape)));
}

py::object mad_value_to_mad(py::handle self, py::handle value) {
    const py::object median_object = self.attr("median");
    log(kLogger, "debug", "self.median: " + str_of(median_object));
    const Values v = values_of(value);
    const double med = to_double(median_object);
    std::vector<double> diff(v.values.size());
    for (std::size_t i = 0; i < diff.size(); ++i) {
        diff[i] = v.values[i] - med;
    }
    log(kLogger, "debug", "diff: " + str_of(values_object(v, diff)));
    const py::object mad_object = self.attr("med_abs_deviation");
    log(kLogger, "debug", "self.med_abs_deviation: " + str_of(mad_object));
    const py::object b_object = self.attr("mad_b_value");
    log(kLogger, "debug", "self.mad_b_value: " + str_of(b_object));
    const double mad = to_double(mad_object);
    const double b = to_double(b_object);
    std::vector<double> together(diff.size());
    std::vector<double> scaled(diff.size());
    for (std::size_t i = 0; i < diff.size(); ++i) {
        scaled[i] = b * diff[i];
        together[i] = scaled[i] / mad;
    }
    log(kLogger, "debug", "all together: " + str_of(values_object(v, together)));
    return values_object(v, mad == 0.0 ? scaled : together);
}

py::object mad_mad_to_value(py::handle self, py::handle mad_value) {
    if (PyList_Check(mad_value.ptr()) || PyTuple_Check(mad_value.ptr())) {
        // mad * self.med_abs_deviation on a sequence
        raise_error("TypeError", "can't multiply sequence by non-int of type '" +
                                     type_name(self.attr("med_abs_deviation")) + "'");
    }
    const Values v = values_of(mad_value);
    const double mad = to_double(self.attr("med_abs_deviation"));
    const double b = to_double(self.attr("mad_b_value"));
    const double med = to_double(self.attr("median"));
    std::vector<double> out(v.values.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = (v.values[i] * mad / b) + med;
    }
    return values_object(v, out);
}

py::object mad_is_outlier(py::handle self, py::handle lower, py::handle upper) {
    const Values z = values_of(self.attr("modified_z_score"));
    const double low = to_double(lower);
    const double high = to_double(upper);
    py::array_t<bool> out(static_cast<py::ssize_t>(z.values.size()));
    bool* flags = out.mutable_data();
    for (std::size_t i = 0; i < z.values.size(); ++i) {
        flags[i] = z.values[i] < low || z.values[i] > high;
    }
    if (z.scalar) {
        return py::module_::import("numpy").attr("bool_")(flags[0]);
    }
    return out.attr("reshape")(py::tuple(py::cast(z.shape)));
}

// ---------------------------------------------------------------------------
// row sums

// matrix[first:last, first:last].sum(axis=1) - .diagonal() (hicx::correct_api::coverage)
py::array coverage(const Sparse& matrix, std::int64_t first, std::int64_t last, bool zero_nan,
                   const char* function) {
    py::array out;
    dispatch_arithmetic(matrix.data, function, [&]<class T>() {
        if constexpr (std::is_same_v<T, float>) {
            unsupported_points("the float32 row sums of a float32 matrix");
        } else {
            out = from_vector(ca::coverage(matrix.indptr_vector(), matrix.entries().col, to_vector<T>(matrix.data),
                                           matrix.rows, first, last, zero_nan));
        }
    });
    return out;
}

py::object make_mad(py::handle row_sum) {
    return this_module().attr("MAD")(row_sum);
}

py::list flatnonzero(py::handle flags_object, std::int64_t offset) {
    const py::array flags = py::array::ensure(flags_object);
    const std::vector<bool> values = [&] {
        std::vector<bool> out;
        for (const auto v : to_vector<bool>(flags)) {
            out.push_back(v);
        }
        return out;
    }();
    py::list out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (values[i]) {
            out.append(numpy_scalar<std::int64_t>(static_cast<std::int64_t>(i) + offset));
        }
    }
    return out;
}

std::pair<std::int64_t, std::int64_t> chromosome_range(py::handle hic, py::handle chrname) {
    const py::object range = hic.attr("getChrBinRange")(chrname);
    return {to_int64(lc::item(range, 0)), to_int64(lc::item(range, 1))};
}

py::object filter_by_zscore(py::handle hic_ma, py::handle lower, py::handle upper, py::handle perchr) {
    log(kLogger, "info", "filtering by z-score");
    py::list to_remove;
    if (lc::truthy(perchr)) {
        const py::list chromosomes = py::list(hic_ma.attr("interval_trees"));
        for (const py::handle chrname : chromosomes) {
            const auto [first, last] = chromosome_range(hic_ma, chrname);
            const Sparse matrix = sparse(hic_ma.attr("matrix"), "filter_by_zscore", "r");
            const py::array row_sum = coverage(matrix, first, last, true, "filter_by_zscore");
            const py::object mad = make_mad(row_sum);
            const py::list problematic = flatnonzero(mad.attr("is_outlier")(lower, upper), first);
            if (problematic.empty()) {
                log(kLogger, "warning",
                    "Warning. No bins removed for chromosome " + str_of(chrname) + " using thresholds " +
                        str_of(lower) + " " + str_of(upper) + "\n");
            }
            for (const py::handle bin : problematic) {
                to_remove.append(bin);
            }
        }
    } else {
        const Sparse matrix = sparse(hic_ma.attr("matrix"), "filter_by_zscore", "r");
        const py::array row_sum = coverage(matrix, 0, matrix.rows, false, "filter_by_zscore");
        const py::object mad = make_mad(row_sum);
        to_remove = flatnonzero(mad.attr("is_outlier")(lower, upper), 0);
    }
    // sorted(to_remove)
    std::vector<std::int64_t> bins;
    for (const py::handle bin : to_remove) {
        bins.push_back(to_int64(bin));
    }
    std::stable_sort(bins.begin(), bins.end());
    py::list sorted_bins;
    for (const std::int64_t bin : bins) {
        sorted_bins.append(numpy_scalar<std::int64_t>(bin));
    }
    return sorted_bins;
}

// ---------------------------------------------------------------------------
// fill_gaps

py::object fill_gaps(py::handle hic_ma, py::handle failed_bins, py::handle fill_contiguous) {
    log(kLogger, "debug", "starting fill gaps");
    const Sparse matrix = sparse(hic_ma.attr("matrix"), "fill_gaps", "rco");
    // the failed bins and the positions of the consecutive ones
    std::vector<py::object> failed;
    for (const py::handle bin : failed_bins) {
        failed.push_back(py::reinterpret_borrow<py::object>(bin));
    }
    const bool contiguous = fill_contiguous.ptr() == Py_True;
    std::vector<std::int64_t> consecutive;
    std::vector<py::object> discontinuous;
    if (contiguous) {
        discontinuous = failed;
    } else {
        std::vector<std::int64_t> values;
        for (const py::object& bin : failed) {
            values.push_back(to_int64(bin));
        }
        std::vector<char> in_run(values.size(), 0);
        for (std::size_t k = 0; k + 1 < values.size(); ++k) {
            if (values[k + 1] - values[k] == 1) {
                in_run[k] = 1;
                in_run[k + 1] = 1;
            }
        }
        for (std::size_t k = 0; k < values.size(); ++k) {
            if (in_run[k]) {
                consecutive.push_back(static_cast<std::int64_t>(k));
            } else {
                discontinuous.push_back(failed[k]);
            }
        }
    }
    log(kLogger, "debug", "Filling " + std::to_string(discontinuous.size()) + " failed bins\n");
    std::vector<std::int64_t> bins;
    for (const py::object& bin : discontinuous) {
        bins.push_back(to_int64(bin));
    }
    py::object filled;
    dispatch_arithmetic(matrix.data, "fill_gaps", [&]<class T>() {
        const hicx::utilities::SparseEntries entries = matrix.entries();
        filled = lc::csr_object(
            ca::fill_gaps(matrix.rows, matrix.cols, entries.row, entries.col, to_vector<T>(matrix.data), bins));
    });
    // np.sort(failed_bins[consecutive_failed_idx])
    if (!is_ndarray(failed_bins)) {
        raise_error("TypeError", "only integer scalar arrays can be converted to a scalar index");
    }
    if (contiguous) {
        raise_error("IndexError", "arrays used as indices must be of integer (or boolean) type");
    }
    const py::array selected = py::reinterpret_borrow<py::array>(failed_bins);
    std::vector<std::size_t> positions(consecutive.begin(), consecutive.end());
    const py::object still_failed = take_elements(selected, positions);
    // the values of one dtype, sorted: integers stay exact
    py::array sorted_bins = py::array::ensure(still_failed).attr("copy")();
    dispatch_numeric(sorted_bins, "fill_gaps", [&]<class T>() {
        if constexpr (std::is_arithmetic_v<T>) {
            std::vector<T> values = to_vector<T>(sorted_bins);
            std::sort(values.begin(), values.end());
            write_back(sorted_bins, values);
        } else {
            raise_unsupported_dtype(sorted_bins, "fill_gaps");
        }
    });
    return py::make_tuple(filled, sorted_bins);
}

// ---------------------------------------------------------------------------
// correct_subparser

py::object correct_subparser() {
    const py::module_ argparse = py::module_::import("argparse");
    const py::object builtins = py::module_::import("builtins");
    const py::object str_type = builtins.attr("str");
    const py::object float_type = builtins.attr("float");
    const py::object int_type = builtins.attr("int");
    const py::object parser = argparse.attr("ArgumentParser")(py::arg("add_help") = false);
    const py::object required = parser.attr("add_argument_group")("Required arguments");
    required.attr("add_argument")("--matrix", "-m", py::arg("help") = "Name of the Hi-C matrix to correct in .h5 format.",
                                  py::arg("required") = true);
    required.attr("add_argument")("--outFileName", "-o",
                                  py::arg("help") = "File name to save the resulting matrix. The output is a .h5 file.",
                                  py::arg("required") = true);
    const py::object optional = parser.attr("add_argument_group")("Optional arguments");
    py::list choices;
    choices.append("KR");
    choices.append("ICE");
    optional.attr("add_argument")(
        "--correctionMethod",
        py::arg("help") = "Method to be used for matrix correction. It can be set to KR or ICE (Default: %(default)s).",
        py::arg("type") = str_type, py::arg("metavar") = "STR", py::arg("default") = "KR",
        py::arg("choices") = choices);
    optional.attr("add_argument")(
        "--filterThreshold", "-t",
        py::arg("help") =
            "Removes bins of low or large coverage. Usually these bins do not contain valid Hi-C data or represent "
            "regions that accumulate reads and thus must be discarded. Use hicCorrectMatrix diagnostic_plot to "
            "identify the modified z-value thresholds. A lower and upper threshold are required separated by "
            "space, e.g. --filterThreshold -1.5 5. Applied only for ICE!",
        py::arg("type") = float_type, py::arg("nargs") = 2, py::arg("default") = py::none());
    optional.attr("add_argument")("--iterNum", "-n",
                                  py::arg("help") = "Number of iterations to compute.only for ICE! (Default: %(default)s).",
                                  py::arg("type") = int_type, py::arg("metavar") = "INT", py::arg("default") = 500);
    optional.attr("add_argument")(
        "--inflationCutoff",
        py::arg("help") =
            "Value corresponding to the maximum number of times a bin can be scaled up during the iterative "
            "correction. For example, an inflation cutoff of 3 will filter out all bins that were expanded 3 times "
            "or more during the iterative correctionself.Only for ICE!",
        py::arg("type") = float_type);
    optional.attr("add_argument")(
        "--transCutoff", "-transcut",
        py::arg("help") = "Clip high counts in the top -transcut trans regions (i.e. between chromosomes). A usual "
                          "value is 0.05. Only for ICE! ",
        py::arg("type") = float_type);
    optional.attr("add_argument")(
        "--sequencedCountCutoff",
        py::arg("help") = "Each bin receives a value indicating the fraction that is covered by reads. A cutoff of "
                          "0.5 will discard all those bins that have less than half of the bin covered. Only for ICE!",
        py::arg("default") = py::none(), py::arg("type") = float_type);
    optional.attr("add_argument")(
        "--chromosomes",
        py::arg("help") = "List of chromosomes to be included in the iterative correction. The order of the given "
                          "chromosomes will be then kept for the resulting corrected matrix",
        py::arg("default") = py::none(), py::arg("nargs") = "+");
    optional.attr("add_argument")("--skipDiagonal", "-s",
                                  py::arg("help") = "If set, diagonal counts are not included. Only for ICE!",
                                  py::arg("action") = "store_true");
    optional.attr("add_argument")(
        "--perchr",
        py::arg("help") = "Normalize each chromosome separately. This is useful for samples from cells with uneven "
                          "number of chromosomes and/or translocations.",
        py::arg("action") = "store_true");
    optional.attr("add_argument")("--filteredBed",
                                  py::arg("help") = "Print bins filtered our by  --filterThreshold to this file");
    optional.attr("add_argument")("--verbose", py::arg("help") = "Print processing status.",
                                  py::arg("action") = "store_true");
    const std::string version = str_of(py::module_::import("hicexplorer._version").attr("__version__"));
    optional.attr("add_argument")("--version", py::arg("action") = "version",
                                  py::arg("version") = "%(prog)s " + version);
    return parser;
}

// ---------------------------------------------------------------------------
// plot_total_contact_dist's data

void prepare_total_contact_dist(py::handle hic_ma) {
    const py::module_ module = this_module();
    hic_ma.attr("matrix") = module.attr("convertNansToZeros")(hic_ma.attr("matrix"));
    hic_ma.attr("matrix") = module.attr("convertInfsToZeros")(hic_ma.attr("matrix"));
}

// row_sum[modified_z_score < 5] and the MAD of the unfiltered row sums
std::pair<py::object, py::object> contact_dist(const py::array& row_sum) {
    const py::object mad = make_mad(row_sum);
    const Values z = values_of(mad.attr("get_motified_zscores")());
    const std::vector<double> sums = lc::doubles(row_sum, "plot_total_contact_dist");
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < z.values.size(); ++i) {
        if (z.values[i] < 5) {
            keep.push_back(i);
        }
    }
    return {take_elements(row_sum, keep), mad};
}

py::object chromosome_contact_dist(py::handle hic_ma, py::handle chrname) {
    log(kLogger, "info", "Plotting chromosome " + str_of(chrname));
    const auto [first, last] = chromosome_range(hic_ma, chrname);
    const Sparse matrix = sparse(hic_ma.attr("matrix"), "plot_total_contact_dist", "r");
    const auto [row_sum, mad] = contact_dist(coverage(matrix, first, last, false, "plot_total_contact_dist"));
    if (lc::length(row_sum) == 0) {
        log(kLogger, "error", "No data for chromosome " + str_of(chrname) + ". Continue with next chromosome.");
    }
    return py::make_tuple(row_sum, mad);
}

py::object genome_contact_dist(py::handle hic_ma) {
    const Sparse matrix = sparse(hic_ma.attr("matrix"), "plot_total_contact_dist", "r");
    const auto [row_sum, mad] = contact_dist(coverage(matrix, 0, matrix.rows, false, "plot_total_contact_dist"));
    if (lc::length(row_sum) == 0) {
        log(kLogger, "error", "No data available, exit.");
        lc::python_exit(1);
    }
    return py::make_tuple(row_sum, mad);
}

// row_sum_values[row_sum_values < x_max]
py::object values_below(py::handle values_object, py::handle limit) {
    const py::array values = py::array::ensure(values_object);
    const std::vector<double> values_double = lc::doubles(values, "plot_total_contact_dist");
    const double bound = to_double(limit);
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < values_double.size(); ++i) {
        if (values_double[i] < bound) {
            keep.push_back(i);
        }
    }
    return take_elements(values, keep);
}

// bin_s[first x with 1 <= x < len(dist) - 1 and dist[x - 1] > dist[x] < dist[x + 1]], or None
py::object first_local_minimum(py::handle dist_object, py::handle bins_object) {
    const std::vector<double> dist = lc::doubles(py::array::ensure(dist_object), "plot_total_contact_dist");
    for (std::size_t x = 1; x + 1 < dist.size(); ++x) {
        if (dist[x - 1] > dist[x] && dist[x] < dist[x + 1]) {
            return lc::item(bins_object, static_cast<std::int64_t>(x));
        }
    }
    return py::none();
}

}  // namespace

void bind_hicCorrectMatrix(py::module_& m) {
    m.def("correct_subparser", &correct_subparser);
    m.def("iterative_correction", [](py::handle matrix, py::handle args) {
        const py::object result = this_module().attr("iterativeCorrection")(
            matrix, py::arg("M") = args.attr("iterNum"), py::arg("verbose") = args.attr("verbose"));
        const py::tuple pair = py::tuple(result);
        if (pair.size() != 2) {
            raise_error("ValueError", pair.size() < 2 ? "not enough values to unpack (expected 2, got " +
                                                           std::to_string(pair.size()) + ")"
                                                     : "too many values to unpack (expected 2)");
        }
        return py::make_tuple(pair[0], pair[1]);
    }, py::arg("matrix"), py::arg("args"));
    m.def("fill_gaps", &fill_gaps, py::arg("hic_ma"), py::arg("failed_bins"), py::arg("fill_contiguous") = false);
    m.def("MAD___init__", &mad_init, py::arg("self"), py::arg("points"));
    m.def("MAD_get_motified_zscores", [](py::handle self) { return self.attr("modified_z_score"); }, py::arg("self"));
    m.def("MAD_is_outlier", &mad_is_outlier, py::arg("self"), py::arg("lower_threshold"), py::arg("upper_threshold"));
    m.def("MAD_value_to_mad", &mad_value_to_mad, py::arg("self"), py::arg("value"));
    m.def("MAD_mad_to_value", &mad_mad_to_value, py::arg("self"), py::arg("mad"));
    m.def("filter_by_zscore", &filter_by_zscore, py::arg("hic_ma"), py::arg("lower_threshold"),
          py::arg("upper_threshold"), py::arg("perchr") = false);
    m.def("prepare_total_contact_dist", &prepare_total_contact_dist, py::arg("hic_ma"));
    m.def("chromosome_contact_dist", &chromosome_contact_dist, py::arg("hic_ma"), py::arg("chrname"));
    m.def("genome_contact_dist", &genome_contact_dist, py::arg("hic_ma"));
    m.def("values_below", &values_below, py::arg("values"), py::arg("limit"));
    m.def("first_local_minimum", &first_local_minimum, py::arg("dist"), py::arg("bins"));
}

}  // namespace hicx::python
