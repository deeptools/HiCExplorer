// hicexplorer._hicx.hicFindTADs: the C++ bodies of hicexplorer.hicFindTADs.
//
// The module functions and every method of HicFindTads keep their 3.7.6
// names, signatures, return types, side effects and exceptions. The
// computations are hicx::find_tads (cpp/core/include/hicx/find_tads.hpp), the
// z-score of hicx::convert_to_obs_exp_matrix and the matrix operations of
// libhicx4; this file converts the Python arguments, performs the hiCMatrix
// methods the Python called on the object's attributes (hicFindTADs.hpp) and
// builds the 3.7.6 results.
//
// The class keeps its attributes (hic_ma, num_processors, max_depth,
// min_depth, step, delta, min_boundary_distance, use_zscore, binsize,
// bedgraph_matrix, boundaries, correct_for_multiple_testing,
// threshold_comparisons, pvalueFDR); each method calls the function here
// with `self`, and calls between methods go through `self` or through the
// class HicFindTads, as in the Python, so a subclass still sees its
// overrides. The module global hic_ma is set as compute_spectra_matrix set
// it.
//
// Deviations from 3.7.6, all recorded in hicexplorer/test/api/cases_hicFindTADs.py:
//  * set_matrix with a file name loads the matrix through the C++ matrix
//    layer (hicx::ToolMatrix, the loader of the tools) into a hiCMatrix whose
//    matrixFileHandler is None, so a later hic_ma.save() picks its writer
//    from the file name instead of from the input.
//  * get_incremental_step_size refuses with ValueError the arguments on
//    which 3.7.6 loops forever (start_step_len <= 0 with min <= max, NaN or
//    an infinite maximum).
//  * Window lengths that are not whole numbers, floats in y_axis other than
//    float64, and matrices other than int and float raise TypeError naming
//    HiCExplorer 4; the Python accepted them only by accident.
//  * compute_spectra_matrix computes the scores with num_processors threads
//    instead of processes; the result is the same for every count, and the
//    count still decides, as in 3.7.6, whether a worker gets an empty share
//    of the bins and the call fails with ValueError.
//  * The CTCF boundary filter of the hicFindTADs tool (--ctcfFile, a
//    HiCExplorer 4 extension with the divergent orientation rule: '-' before,
//    '+' after a boundary) has no 3.7.6 API name; none of the names here
//    filters by CTCF.

#include "hicFindTADs.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "api.hpp"
#include "hicx/adjust_ops.hpp"
#include "hicx/find_tads.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/numpy_sort.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/parallel.hpp"
#include "hicx/tad_contacts.hpp"
#include "hicx/tool_matrix.hpp"
#include "utilities.hpp"

namespace hicx::python {

namespace ft = hicx::find_tads;

namespace {

constexpr const char* kLogger = "hicexplorer.hicFindTADs";
constexpr const char* kMatrixLogger = "hicmatrix.HiCMatrix";
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// --------------------------------------------------------------------------
// Python object protocol

bool truthy(py::handle object) {
    const int result = PyObject_IsTrue(object.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool compare(py::handle a, py::handle b, int op) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), op);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool contains(py::handle container, py::handle item) {
    const int result = PySequence_Contains(container.ptr(), item.ptr());
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

std::size_t length(py::handle object) {
    const Py_ssize_t n = PyObject_Length(object.ptr());
    if (n < 0) {
        throw_python_error();
    }
    return static_cast<std::size_t>(n);
}

py::object item(py::handle sequence, std::int64_t index) {
    PyObject* value = PySequence_GetItem(sequence.ptr(), static_cast<Py_ssize_t>(index));
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

py::object get_item(py::handle mapping, py::handle key) {
    PyObject* value = PyObject_GetItem(mapping.ptr(), key.ptr());
    if (value == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(value);
}

void set_item(py::handle mapping, py::handle key, py::handle value) {
    if (PyObject_SetItem(mapping.ptr(), key.ptr(), value.ptr()) != 0) {
        throw_python_error();
    }
}

py::list as_list(py::handle iterable) {
    PyObject* list = PySequence_List(iterable.ptr());
    if (list == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::list>(list);
}

// `a, b, c = value` with its errors.
std::vector<py::object> unpack(py::handle value, std::size_t expected) {
    PyObject* iterator = PyObject_GetIter(value.ptr());
    if (iterator == nullptr) {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) {
            PyErr_Clear();
            raise_error("TypeError", "cannot unpack non-iterable " + type_name(value) + " object");
        }
        throw_python_error();
    }
    py::object it = py::reinterpret_steal<py::object>(iterator);
    std::vector<py::object> out;
    while (true) {
        PyObject* next = PyIter_Next(it.ptr());
        if (next == nullptr) {
            if (PyErr_Occurred()) {
                throw_python_error();
            }
            break;
        }
        if (out.size() == expected) {
            Py_DECREF(next);
            raise_error("ValueError",
                        "too many values to unpack (expected " + std::to_string(expected) + ")");
        }
        out.push_back(py::reinterpret_steal<py::object>(next));
    }
    if (out.size() < expected) {
        raise_error("ValueError", "not enough values to unpack (expected " + std::to_string(expected) +
                                      ", got " + std::to_string(out.size()) + ")");
    }
    return out;
}

py::module_ numpy() {
    return py::module_::import("numpy");
}

py::module_ this_module() {
    return py::module_::import("hicexplorer.hicFindTADs");
}

// int(value)
py::object python_int(py::handle value) {
    PyObject* result = PyNumber_Long(value.ptr());
    if (result == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(result);
}

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

bool is_numpy_integer(py::handle value) {
    return py::isinstance(value, numpy().attr("integer"));
}

bool is_numpy_floating(py::handle value) {
    return py::isinstance(value, numpy().attr("floating"));
}

// '{}'.format(value)
std::string fmt(py::handle value) {
    return utf8(py::str(py::module_::import("builtins").attr("format")(value, "")));
}

// A whole number of base pairs (a window length or depth) as int64.
std::int64_t whole_number(py::handle value, const char* what) {
    if (is_integer(value)) {
        return to_int64(value);
    }
    const double real = to_double(value);
    if (std::isfinite(real) && std::trunc(real) == real && std::fabs(real) < 9.0e15) {
        return static_cast<std::int64_t>(real);
    }
    raise_error("TypeError", std::string("HiCExplorer 4 needs a whole number for ") + what +
                                 ", got " + repr_of(value));
}

// --------------------------------------------------------------------------
// numbers that keep their Python type (int, numpy.int64, float, numpy.float64)

struct Number {
    enum class Kind { PyInt, NpInt, PyFloat, NpFloat };
    Kind kind = Kind::PyInt;
    std::int64_t integer = 0;
    double real = 0.0;

    [[nodiscard]] bool integral() const { return kind == Kind::PyInt || kind == Kind::NpInt; }
    [[nodiscard]] double value() const { return integral() ? static_cast<double>(integer) : real; }

    static Number of(py::handle object) {
        Number n;
        if (is_numpy_integer(object)) {
            n.kind = Kind::NpInt;
            n.integer = to_int64(object);
        } else if (is_numpy_floating(object)) {
            n.kind = Kind::NpFloat;
            n.real = to_double(object);
        } else if (PyLong_Check(object.ptr())) {
            n.kind = Kind::PyInt;
            n.integer = to_int64(object);
        } else if (PyFloat_Check(object.ptr())) {
            n.kind = Kind::PyFloat;
            n.real = to_double(object);
        } else {
            raise_error("TypeError", "HiCExplorer 4 needs a number here, got " + repr_of(object));
        }
        return n;
    }

    [[nodiscard]] Number times(std::int64_t factor) const {
        Number n = *this;
        if (integral()) {
            n.integer = integer * factor;
        } else {
            n.real = real * static_cast<double>(factor);
        }
        return n;
    }

    [[nodiscard]] py::object object() const {
        switch (kind) {
            case Kind::PyInt: return py::int_(integer);
            case Kind::NpInt: return numpy_scalar<std::int64_t>(integer);
            case Kind::PyFloat: return py::float_(real);
            case Kind::NpFloat: return numpy_scalar<double>(real);
        }
        return py::none();
    }
};

// --------------------------------------------------------------------------
// numeric arrays

std::vector<double> doubles_of(py::handle object, const char* function) {
    if (is_ndarray(object)) {
        const py::array array = py::reinterpret_borrow<py::array>(object);
        std::vector<double> out;
        dispatch_numeric(array, function, [&]<class T>() {
            if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, long double>) {
                const std::vector<T> values = to_vector<T>(array);
                out.assign(values.begin(), values.end());
            } else {
                raise_unsupported_dtype(array, function);
            }
        });
        return out;
    }
    std::vector<double> out;
    for (const py::handle value : py::reinterpret_borrow<py::iterable>(object)) {
        out.push_back(to_double(value));
    }
    return out;
}

// A row major float64 2-D array, its rows and columns.
struct Matrix2D {
    std::vector<double> values;
    std::size_t rows = 0;
    std::size_t columns = 0;
};

Matrix2D matrix_2d(py::handle object, const char* function) {
    const py::array array = as_ndarray(object);
    if (array.ndim() != 2) {
        raise_error("TypeError", std::string("HiCExplorer 4: ") + function +
                                     " needs a two dimensional array, got " +
                                     std::to_string(array.ndim()) + " dimensions");
    }
    if (dtype_of(array) != DType::Float64) {
        raise_error("TypeError", std::string("HiCExplorer 4: ") + function +
                                     " needs a float64 array, got " + dtype_name(array));
    }
    Matrix2D m;
    m.values = to_vector<double>(array);
    m.rows = static_cast<std::size_t>(array.shape(0));
    m.columns = static_cast<std::size_t>(array.shape(1));
    return m;
}

py::array float64_array(const std::vector<double>& values) {
    return from_vector(values);
}

py::array int64_array(const std::vector<std::int64_t>& values) {
    return from_vector(values);
}

py::array matrix_object(const std::vector<double>& values, std::size_t rows, std::size_t columns,
                        bool integral) {
    py::array out = from_vector(values);
    out = out.attr("reshape")(py::make_tuple(rows, columns));
    if (integral) {
        out = out.attr("astype")("int64");
    }
    return out;
}

// np.array(list of str), a '<U' array.
py::array str_array(const py::list& values) {
    return numpy().attr("array")(values);
}

std::vector<std::string> strings_of(py::handle sequence) {
    std::vector<std::string> out;
    for (const py::handle value : py::reinterpret_borrow<py::iterable>(sequence)) {
        out.push_back(str_of(value));
    }
    return out;
}

}  // namespace

// ==========================================================================
// hiCMatrix methods on the object (hicFindTADs.hpp)

namespace hic {

std::pair<py::object, py::object> interval_trees(py::handle intervals) {
    py::dict trees;
    py::object boundaries = py::module_::import("collections").attr("OrderedDict")();
    if (length(intervals) == 0) {
        log(kMatrixLogger, "warning", "Interval list is empty");
        return {trees, boundaries};
    }
    const py::module_ intervaltree = py::module_::import("intervaltree");
    const py::object Interval = intervaltree.attr("Interval");
    std::int64_t id = 0;
    std::int64_t chrom_start = 0;
    py::object previous = py::none();
    py::object chrom = py::none();
    py::list current;
    for (const py::handle interval : py::reinterpret_borrow<py::iterable>(intervals)) {
        const py::object fields = py::reinterpret_borrow<py::object>(interval)[py::slice(0, 3, 1)];
        const std::vector<py::object> values = unpack(fields, 3);
        chrom = values[0];
        const py::object start = py::int_(values[1]);
        const py::object end = py::int_(values[2]);
        if (compare(previous, chrom, Py_NE)) {
            if (previous.is_none()) {
                previous = chrom;
            }
            set_item(boundaries, previous, py::make_tuple(chrom_start, id));
            chrom_start = id;
            current = py::list();
            set_item(trees, chrom, current);
            previous = chrom;
        }
        current.append(Interval(start, end, id));
        ++id;
    }
    set_item(boundaries, chrom, py::make_tuple(chrom_start, id));
    const py::object IntervalTree = intervaltree.attr("IntervalTree");
    for (const auto& [key, value] : py::dict(trees)) {
        set_item(trees, key, IntervalTree(value));
    }
    return {trees, boundaries};
}

void set_intervals(py::handle hic, py::handle intervals) {
    hic.attr("cut_intervals") = intervals;
    auto [trees, boundaries] = interval_trees(intervals);
    hic.attr("interval_trees") = trees;
    hic.attr("chrBinBoundaries") = boundaries;
}

hicx::BinTable bins(py::handle hic) {
    return hicx::BinTable(cut_intervals(hic.attr("cut_intervals")));
}

std::optional<std::pair<std::int64_t, std::int64_t>> region_bin_range(
    py::handle hic, const hicx::BinTable& table, py::handle chrname, py::handle startpos,
    py::handle endpos) {
    const py::object trees = hic.attr("interval_trees");
    py::object name = py::reinterpret_borrow<py::object>(chrname);
    PyObject* iterator = PyObject_GetIter(trees.ptr());
    if (iterator == nullptr) {
        throw_python_error();
    }
    py::object it = py::reinterpret_steal<py::object>(iterator);
    PyObject* first_raw = PyIter_Next(it.ptr());
    if (first_raw == nullptr) {
        if (!PyErr_Occurred()) {
            PyErr_SetNone(PyExc_StopIteration);
        }
        throw_python_error();
    }
    const py::object first = py::reinterpret_steal<py::object>(first_raw);
    if (!py::isinstance(first, py::type::of(name))) {
        if (py::isinstance<py::str>(first)) {
            name = to_string(name);
        } else if (PyBytes_Check(first.ptr()) || py::isinstance(first, numpy().attr("bytes_"))) {
            name = to_bytes(name);
        }
    }
    if (!contains(trees, name)) {
        raise_error("ValueError", "chromosome: " + str_of(name) +
                                      " name not found in matrixvalid names are: " +
                                      str_of(trees.attr("keys")()));
    }
    const std::int64_t start = to_int64(python_int(startpos));
    const std::int64_t end = to_int64(python_int(endpos));
    std::string key;
    if (PyBytes_Check(name.ptr())) {
        key = std::string(py::bytes(name));
    } else {
        key = str_of(name);
    }
    const auto range = table.region_bin_range(key, start, end);
    if (!range.has_value()) {
        log(kMatrixLogger, "error", "Index error");
    }
    return range;
}

py::object bin_size(py::handle hic) {
    py::object cached = hic.attr("bin_size");
    if (!cached.is_none()) {
        return cached;
    }
    const hicx::BinTable table = bins(hic);
    if (table.size() == 1) {
        cached = numpy_scalar<std::int64_t>(table.intervals()[0].end - table.intervals()[0].start);
    } else {
        cached = py::int_(table.bin_size());
        if (!table.bin_size_homogeneous()) {
            hic.attr("bin_size_homogeneous") = py::bool_(false);
            if (!truthy(hic.attr("non_homogeneous_warning_already_printed"))) {
                log(kMatrixLogger, "warning",
                    "Bin size is not homogeneous.                                       Median " +
                        hicx::find_tads::format_fixed(static_cast<double>(table.bin_size()), 6) +
                        "\n");
                hic.attr("non_homogeneous_warning_already_printed") = py::bool_(true);
            }
        }
    }
    hic.attr("bin_size") = cached;
    return cached;
}

py::object bin_pos(py::handle hic, py::handle index) {
    const py::object intervals = hic.attr("cut_intervals");
    if (compare(index, py::int_(length(intervals)), Py_LT)) {
        return get_item(intervals, index);
    }
    raise_error("ValueError", "binIndex: " + fmt(index) + " not found");
}

hicx::CsrMatrix matrix(py::handle hic) {
    hicx::CsrMatrix m = to_csr_matrix(sparse(hic.attr("matrix"), "hiCMatrix", "r"));
    hicx::tads::sort_row_indices(m);
    return m;
}

py::object csr_object(hicx::CsrMatrix m) {
    if (m.symmetry() == hicx::Symmetry::UpperTriangle) {
        m.materialize_full();
    }
    return csr_matrix_object(m);
}

void restore_masked_bins(py::handle hic) {
    const py::object orig_ids = hic.attr("orig_bin_ids");
    if (length(orig_ids) == 0) {
        return;
    }
    std::vector<std::int64_t> ids;
    for (const py::handle id : py::reinterpret_borrow<py::iterable>(orig_ids)) {
        ids.push_back(to_int64(id));
    }
    hicx::CsrMatrix m = matrix(hic);
    const std::int64_t kept = m.rows();
    const std::int64_t total = static_cast<std::int64_t>(ids.size());
    const std::int64_t added = total - kept;
    // vstack and hstack with empty float64 blocks, then rows = cols =
    // argsort(orig_bin_ids)
    std::vector<std::int64_t> indptr(m.indptr());
    indptr.resize(static_cast<std::size_t>(total) + 1, indptr.back());
    const std::string dtype = m.dtype() == "float64" || m.integral_dtype() || m.dtype() == "float32"
                                  ? std::string("float64")
                                  : m.dtype();
    hicx::CsrMatrix padded(total, total, std::move(indptr), std::vector<std::int32_t>(m.indices()),
                           std::vector<double>(m.data()), dtype);
    const std::vector<std::int64_t> order = hicx::npy::argsort(ids);
    hicx::CsrMatrix restored = hicx::select_bins(padded, order);
    hic.attr("matrix") = csr_object(restored);

    const py::object orig_intervals = hic.attr("orig_cut_intervals");
    py::list intervals;
    for (const std::int64_t k : order) {
        intervals.append(item(orig_intervals, k));
    }
    set_intervals(hic, intervals);
    hic.attr("nan_bins") = py::reinterpret_borrow<py::object>(orig_ids)[py::slice(kept, total, 1)];
    const py::object factors = hic.attr("correction_factors");
    if (!factors.is_none()) {
        std::vector<double> values = doubles_of(factors, "restoreMaskedBins");
        values.resize(values.size() + static_cast<std::size_t>(std::max<std::int64_t>(0, added)), kNaN);
        std::vector<double> reordered;
        for (const std::int64_t k : order) {
            reordered.push_back(values[static_cast<std::size_t>(k)]);
        }
        hic.attr("correction_factors") = float64_array(reordered);
    }
    hic.attr("orig_bin_ids") = py::list();
    hic.attr("orig_cut_intervals") = py::list();
    log(kMatrixLogger, "info", "masked bins were restored\n");
}

void write_text_file(py::handle path, const std::string& text) {
    const py::object handle = py::module_::import("builtins").attr("open")(path, "w");
    try {
        handle.attr("write")(to_str(text));
    } catch (...) {
        handle.attr("close")();
        throw;
    }
    handle.attr("close")();
}

void python_exit(py::object code) {
    py::module_::import("builtins").attr("exit")(code);
    raise_error("SystemExit", fmt(code));
}

}  // namespace hic

namespace {

// --------------------------------------------------------------------------
// The hiCMatrix operations compute_spectra_matrix performs

// hiCMatrix.maskBins(bin_ids) as compute_spectra_matrix calls it.
void mask_bins(py::handle hic, py::handle bin_ids_object) {
    if (bin_ids_object.is_none() || length(bin_ids_object) == 0) {
        return;
    }
    std::vector<std::int64_t> bin_ids;
    const auto read_ids = [](py::handle ids) {
        std::vector<std::int64_t> out;
        if (is_ndarray(ids)) {
            const py::array array = py::reinterpret_borrow<py::array>(ids);
            const DType kind = dtype_of(array);
            if (kind == DType::Float16 || kind == DType::Float32 || kind == DType::Float64 ||
                kind == DType::LongDouble) {
                raise_error("IndexError",
                            "arrays used as indices must be of integer (or boolean) type");
            }
        }
        for (const py::handle id : py::reinterpret_borrow<py::iterable>(ids)) {
            out.push_back(to_int64(id));
        }
        return out;
    };
    bin_ids = read_ids(bin_ids_object);
    // printchrtoremove(bin_ids, restore_masked_bins=False): the count per
    // chromosome on hicmatrix's logger, unless the same bins were reported last
    {
        bool same = false;
        if (py::hasattr(hic, "prev_to_remove")) {
            same = truthy(numpy().attr("array_equal")(hic.attr("prev_to_remove"), bin_ids_object));
        } else {
            log(kMatrixLogger, "debug", "No self.prev_to_remove defined, defining it now.");
        }
        if (!same) {
            const py::object intervals = hic.attr("cut_intervals");
            py::dict count;
            for (const std::int64_t id : bin_ids) {
                const py::object chrom = item(item(intervals, id), 0);
                const py::object previous = count.contains(chrom) ? py::object(count[chrom]) : py::int_(0);
                count[chrom] = py::int_(previous.cast<std::int64_t>() + 1);
            }
            log(kMatrixLogger, "info",
                "Number of poor regions to remove: " + std::to_string(bin_ids.size()) + " " +
                    repr_of(count));
            hic.attr("prev_to_remove") = bin_ids_object;
        }
    }
    hicx::CsrMatrix m = hic::matrix(hic);
    // an already masked matrix: restore first, mask the union
    const py::object orig_ids = hic.attr("orig_bin_ids");
    if (length(orig_ids) > 0) {
        try {
            std::vector<std::int64_t> ids = read_ids(orig_ids);
            std::set<std::int64_t> merged(ids.begin() + m.rows(), ids.end());
            for (const std::int64_t id : bin_ids) {
                const std::int64_t at = id < 0 ? id + static_cast<std::int64_t>(ids.size()) : id;
                if (at < 0 || at >= static_cast<std::int64_t>(ids.size())) {
                    throw std::out_of_range("index");
                }
                merged.insert(ids[static_cast<std::size_t>(at)]);
            }
            bin_ids.assign(merged.begin(), merged.end());
            hic::restore_masked_bins(hic);
            m = hic::matrix(hic);
        } catch (const std::out_of_range&) {
            // the Python's bare except: pass
        }
    }
    const py::object nan_bins = hic.attr("nan_bins");
    if (!nan_bins.is_none() && length(nan_bins) > 0) {
        log(kMatrixLogger, "info",
            "found existing " + std::to_string(length(nan_bins)) +
                " nan bins that will be included for masking ");
        std::set<std::int64_t> merged(bin_ids.begin(), bin_ids.end());
        for (const std::int64_t id : read_ids(nan_bins)) {
            merged.insert(id);
        }
        bin_ids.assign(merged.begin(), merged.end());
        hic.attr("nan_bins") = py::list();
    }
    const std::int64_t n = m.cols();
    std::vector<bool> drop(static_cast<std::size_t>(n), false);
    for (std::int64_t& id : bin_ids) {
        const std::int64_t at = id < 0 ? id + n : id;
        if (at < 0 || at >= n) {
            raise_error("IndexError", "index " + std::to_string(id) +
                                          " is out of bounds for axis 0 with size " + std::to_string(n));
        }
        drop[static_cast<std::size_t>(at)] = true;
    }
    std::vector<std::int64_t> rows;
    for (std::int64_t k = 0; k < n; ++k) {
        if (!drop[static_cast<std::size_t>(k)]) {
            rows.push_back(k);
        }
    }
    hic.attr("matrix") = hic::csr_object(hicx::select_bins(m, rows));
    std::vector<std::int64_t> orig = rows;
    orig.insert(orig.end(), bin_ids.begin(), bin_ids.end());
    hic.attr("orig_bin_ids") = int64_array(orig);
    const py::object intervals = hic.attr("cut_intervals");
    py::list kept;
    for (const std::int64_t r : rows) {
        kept.append(item(intervals, r));
    }
    py::list orig_intervals = as_list(kept);
    for (const std::int64_t id : bin_ids) {
        orig_intervals.append(item(intervals, id));
    }
    hic.attr("orig_cut_intervals") = orig_intervals;
    hic::set_intervals(hic, kept);
    const py::object factors = hic.attr("correction_factors");
    if (!factors.is_none()) {
        const std::vector<double> values = doubles_of(factors, "maskBins");
        std::vector<double> selected;
        for (const std::int64_t r : rows) {
            selected.push_back(values[static_cast<std::size_t>(r)]);
        }
        hic.attr("correction_factors") = float64_array(selected);
    }
}

// The context of the TAD-score kernels for the object's current matrix.
struct Scores {
    hicx::CsrMatrix matrix;
    std::vector<hicx::CutInterval> intervals;
    hicx::BinTable bins;
    ft::ScoreContext context;

    explicit Scores(py::handle hic)
        : matrix(hic::matrix(hic)), intervals(cut_intervals(hic.attr("cut_intervals"))),
          bins(intervals) {
        context.matrix = &matrix;
        context.cut_intervals = &intervals;
        context.bins = &bins;
        for (const auto& [chrom, size] : bins.chromosome_sizes()) {
            context.chromosome_sizes[chrom] = size;
        }
    }
    Scores(const Scores&) = delete;
    Scores& operator=(const Scores&) = delete;
};

// get_idx_of_bins_at_given_distance with the Python's errors.
std::pair<std::int64_t, std::int64_t> bins_at_given_distance(py::handle hic, const Scores& scores,
                                                             py::handle idx, std::int64_t window) {
    const std::vector<py::object> position = unpack(hic::bin_pos(hic, idx), 4);
    const py::object chrom = position[0];
    const std::int64_t cut_start = to_int64(position[1]);
    const std::int64_t cut_end = to_int64(position[2]);
    const std::int64_t left_start = std::max<std::int64_t>(0, cut_start - window);
    const auto left = hic::region_bin_range(hic, scores.bins, chrom, py::int_(left_start),
                                            py::int_(left_start + 1));
    if (!left.has_value()) {
        raise_error("TypeError", "'NoneType' object is not subscriptable");
    }
    const std::string name = str_of(chrom);
    const auto size = scores.context.chromosome_sizes.find(name);
    if (size == scores.context.chromosome_sizes.end()) {
        raise_error("KeyError", repr_of(chrom));
    }
    const std::int64_t right_end = std::min(size->second, cut_end + window) - 1;
    const auto right = hic::region_bin_range(hic, scores.bins, chrom, py::int_(right_end),
                                             py::int_(right_end));
    if (!right.has_value()) {
        raise_error("TypeError", "'NoneType' object is not subscriptable");
    }
    return {left->first, right->first};
}

// .todense().A1 of a block, in the matrix dtype.
py::array block_array(const std::vector<double>& values, const std::string& dtype) {
    py::array out = from_vector(values);
    return out.attr("astype")(dtype);
}

// numpy's mean of a dense block in the matrix dtype: float64 for float64 and
// integer matrices, float32 for float32.
py::object block_mean(const std::vector<double>& values, const std::string& dtype) {
    if (dtype == "float32") {
        std::vector<float> narrow(values.begin(), values.end());
        const float sum = hicx::npy::pairwise_sum(narrow.data(), narrow.size());
        return numpy_scalar<float>(sum / static_cast<float>(narrow.size()));
    }
    return numpy_scalar<double>(hicx::npy::pairwise_sum(values.data(), values.size()) /
                                static_cast<double>(values.size()));
}

void require_numeric_matrix(const hicx::CsrMatrix& m, const char* function) {
    const std::string& d = m.dtype();
    if (d != "float64" && d != "float32" && !m.integral_dtype()) {
        raise_error("TypeError", std::string("HiCExplorer 4: ") + function +
                                     " needs an integer, float32 or float64 matrix, got " + d);
    }
}

// get_cut_weight(hic_matrix, cut, window_len, return_mean): None, 0, the mean
// or the block. `scores` is the object's matrix.
py::object cut_weight(py::handle hic, const Scores& scores, py::handle cut, py::handle window_len,
                      bool return_mean) {
    const std::int64_t n = scores.matrix.rows();
    if (compare(cut, py::int_(0), Py_LT) || compare(cut, py::int_(n), Py_GT)) {
        return py::none();
    }
    const std::int64_t window = whole_number(window_len, "window_len");
    std::pair<std::int64_t, std::int64_t> range;
    try {
        range = bins_at_given_distance(hic, scores, cut, window);
    } catch (const hicx::PythonError& error) {
        if (error.type() == "TypeError") {
            return py::none();
        }
        throw;
    }
    const std::int64_t at = to_int64(cut);
    const auto [r0, r1] = hicx::tads::normalise_slice(range.first, at, n);
    const auto [c0, c1] = hicx::tads::normalise_slice(at, range.second, scores.matrix.cols());
    ft::DenseBlock block;
    ft::extract_block(scores.matrix, r0, r1, c0, c1, block);
    if (return_mean) {
        if (block.stored == 0) {
            return py::int_(0);
        }
        return block_mean(block.values, scores.matrix.dtype());
    }
    return block_array(block.values, scores.matrix.dtype());
}

// The rows compute_matrix computes for `bins_list`: for every bin whose cut
// weights are all defined and not NaN, its position objects and its scores.
struct ScoreRows {
    py::list chrom;
    py::list start;
    py::list end;
    std::vector<double> values;
    std::size_t rows = 0;
    bool any_float = false;
    std::vector<char> kept;  // per entry of bins_list
};

ScoreRows score_rows(py::handle hic, const Scores& scores, const std::vector<std::int64_t>& bins_list,
                     const std::vector<std::int64_t>& windows, unsigned int threads) {
    const py::object intervals = hic.attr("cut_intervals");
    const std::int64_t count = static_cast<std::int64_t>(length(intervals));
    for (const std::int64_t cut : bins_list) {
        // cut_intervals[cut] first, then getBinPos(cut) inside get_cut_weight
        const std::int64_t at = cut < 0 ? cut + count : cut;
        if (at < 0 || at >= count) {
            raise_error("IndexError", "list index out of range");
        }
    }
    std::vector<double> values(bins_list.size() * windows.size(), 0.0);
    std::vector<char> keep(bins_list.size(), 0);
    std::vector<char> floats(bins_list.size(), 0);
    require_numeric_matrix(scores.matrix, "compute_matrix");
    hicx::parallel_for(bins_list.size(), threads, [&](std::size_t index) {
        ft::DenseBlock scratch;
        const std::int64_t cut = bins_list[index];
        bool usable = true;
        bool any_float = false;
        for (std::size_t w = 0; w < windows.size(); ++w) {
            bool defined = true;
            double value = 0.0;
            if (cut < 0) {
                defined = false;
            } else {
                value = ft::cut_weight_mean(scores.context, cut, windows[w], scratch, defined);
                if (defined && scratch.stored != 0) {
                    any_float = true;
                }
            }
            if (!defined || std::isnan(value)) {
                usable = false;
                break;
            }
            values[index * windows.size() + w] = value;
        }
        keep[index] = usable ? 1 : 0;
        floats[index] = any_float ? 1 : 0;
    });
    ScoreRows out;
    out.kept = keep;
    for (std::size_t index = 0; index < bins_list.size(); ++index) {
        if (keep[index] == 0) {
            continue;
        }
        const std::vector<py::object> position = unpack(item(intervals, bins_list[index]), 4);
        out.chrom.append(position[0]);
        out.start.append(position[1]);
        out.end.append(position[2]);
        out.values.insert(out.values.end(),
                          values.begin() + static_cast<std::ptrdiff_t>(index * windows.size()),
                          values.begin() + static_cast<std::ptrdiff_t>((index + 1) * windows.size()));
        out.any_float = out.any_float || floats[index] != 0;
        ++out.rows;
    }
    return out;
}

// get_incremental_step_size with the result type of min_win_size + int(...).
py::list step_sizes(py::handle min_win_size, py::handle max_win_size, py::handle start_step_len) {
    const Number minimum = Number::of(min_win_size);
    const Number maximum = Number::of(max_win_size);
    const Number step = Number::of(start_step_len);
    const std::vector<double> values =
        ft::incremental_step_size_real(minimum.value(), maximum.value(), step.value());
    py::list out;
    for (const double v : values) {
        Number n = minimum;
        if (n.integral()) {
            n.integer = static_cast<std::int64_t>(v);
        } else {
            n.real = v;
        }
        out.append(n.object());
    }
    return out;
}

std::vector<std::int64_t> windows_of(const py::list& steps) {
    std::vector<std::int64_t> out;
    for (const py::handle s : steps) {
        out.push_back(whole_number(s, "the window sizes"));
    }
    return out;
}

// --------------------------------------------------------------------------
// module functions

py::object compute_matrix_wrapper(py::handle args) {
    PyObject* tuple = PySequence_Tuple(args.ptr());
    if (tuple == nullptr) {
        if (PyErr_ExceptionMatches(PyExc_TypeError)) {
            PyErr_Clear();
            raise_error("TypeError", "hicexplorer.hicFindTADs.compute_matrix() argument after * must be an "
                                     "iterable, not " +
                                         type_name(args));
        }
        throw_python_error();
    }
    const py::object arguments = py::reinterpret_steal<py::object>(tuple);
    PyObject* result = PyObject_Call(this_module().attr("compute_matrix").ptr(), arguments.ptr(), nullptr);
    if (result == nullptr) {
        throw_python_error();
    }
    return py::reinterpret_steal<py::object>(result);
}

py::object get_cut_weight_by_bin_id(py::handle matrix_object, py::handle cut, py::handle depth,
                                    py::handle return_mean) {
    Sparse s = sparse(matrix_object, "get_cut_weight_by_bin_id", "r");
    if (compare(cut, py::int_(0), Py_LT) || compare(cut, py::int_(s.rows), Py_GT)) {
        return py::none();
    }
    hicx::CsrMatrix m = to_csr_matrix(s);
    hicx::tads::sort_row_indices(m);
    require_numeric_matrix(m, "get_cut_weight_by_bin_id");
    const std::int64_t at = to_int64(cut);
    const std::int64_t d = whole_number(depth, "depth");
    const std::int64_t start = std::max<std::int64_t>(0, at - d);
    const std::int64_t end = std::min<std::int64_t>(s.rows, at + d);
    const auto [r0, r1] = hicx::tads::normalise_slice(start, at, s.rows);
    const auto [c0, c1] = hicx::tads::normalise_slice(at, end, s.cols);
    if (return_mean.ptr() == Py_True) {
        const std::int64_t size = (r1 - r0) * (c1 - c0);
        if (size == 0) {
            raise_error("ZeroDivisionError", "float division by zero");
        }
        // (m.astype(float64) * (1. / size)) @ ones, row by row, then numpy's sum
        const double scale = 1.0 / static_cast<double>(size);
        std::vector<double> row_sums;
        for (std::int64_t row = r0; row < r1; ++row) {
            double sum = 0.0;
            const std::size_t b = static_cast<std::size_t>(m.indptr()[static_cast<std::size_t>(row)]);
            const std::size_t e = static_cast<std::size_t>(m.indptr()[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = b; k < e; ++k) {
                const std::int64_t column = m.indices()[k];
                if (column >= c0 && column < c1) {
                    sum += m.data()[k] * scale;
                }
            }
            row_sums.push_back(sum);
        }
        if (m.dtype() == "float32") {
            std::vector<float> narrow(row_sums.begin(), row_sums.end());
            return numpy_scalar<float>(hicx::npy::pairwise_sum(narrow.data(), narrow.size()));
        }
        return numpy_scalar<double>(hicx::npy::pairwise_sum(row_sums.data(), row_sums.size()));
    }
    ft::DenseBlock block;
    ft::extract_block(m, r0, r1, c0, c1, block);
    return block_array(block.values, m.dtype());
}

py::object get_idx_of_bins_at_given_distance(py::handle hic_matrix, py::handle idx, py::handle window_len) {
    const Scores scores(hic_matrix);
    const auto [left, right] =
        bins_at_given_distance(hic_matrix, scores, idx, whole_number(window_len, "window_len"));
    return py::make_tuple(left, right);
}

py::object get_cut_weight(py::handle hic_matrix, py::handle cut, py::handle window_len, py::handle return_mean) {
    const Scores scores(hic_matrix);
    require_numeric_matrix(scores.matrix, "get_cut_weight");
    return cut_weight(hic_matrix, scores, cut, window_len, return_mean.ptr() == Py_True);
}

py::object get_triangle(py::handle hic_matrix, py::handle cut, py::handle window_len, py::handle /*return_mean*/) {
    const Scores scores(hic_matrix);
    const std::int64_t n = scores.matrix.rows();
    if (compare(cut, py::int_(0), Py_LT) || compare(cut, py::int_(n), Py_GT)) {
        return py::none();
    }
    require_numeric_matrix(scores.matrix, "get_triangle");
    const auto [left, right] =
        bins_at_given_distance(hic_matrix, scores, cut, whole_number(window_len, "window_len"));
    const std::int64_t at = to_int64(cut);
    std::vector<double> values;
    const auto upper = [&](std::int64_t from, std::int64_t to) {
        const auto [b, e] = hicx::tads::normalise_slice(from, to, n);
        ft::DenseBlock block;
        ft::extract_block(scores.matrix, b, e, b, e, block);
        for (std::int64_t r = 0; r < block.rows; ++r) {
            for (std::int64_t c = r; c < block.columns; ++c) {
                values.push_back(block.values[static_cast<std::size_t>(r * block.columns + c)]);
            }
        }
    };
    upper(left, at);
    upper(at, right);
    return block_array(values, scores.matrix.dtype());
}

py::object get_incremental_step_size(py::handle min_win_size, py::handle max_win_size, py::handle start_step_len) {
    return step_sizes(min_win_size, max_win_size, start_step_len);
}

std::vector<std::int64_t> index_list(py::handle values) {
    std::vector<std::int64_t> out;
    for (const py::handle v : py::reinterpret_borrow<py::iterable>(values)) {
        out.push_back(to_int64(v));
    }
    return out;
}

py::object compute_matrix(py::handle bins_list, py::handle min_win_size, py::handle max_win_size,
                          py::handle step_len) {
    const py::object hic = this_module().attr("hic_ma");
    const py::list steps = step_sizes(min_win_size, max_win_size, step_len);
    const std::vector<std::int64_t> bins = index_list(bins_list);
    if (bins.empty()) {
        raise_error("ValueError", "not enough values to unpack (expected 3, got 0)");
    }
    const Scores scores(hic);
    const ScoreRows rows = score_rows(hic, scores, bins, windows_of(steps), 1);
    if (rows.rows == 0) {
        raise_error("ValueError", "not enough values to unpack (expected 3, got 0)");
    }
    const std::size_t columns = rows.rows == 0 ? 0 : rows.values.size() / rows.rows;
    return py::make_tuple(py::tuple(rows.chrom), py::tuple(rows.start), py::tuple(rows.end),
                          matrix_object(rows.values, rows.rows, columns, !rows.any_float));
}

void print_args(py::handle args) {
    for (const py::handle pair : py::reinterpret_borrow<py::iterable>(args.attr("_get_kwargs")())) {
        const std::vector<py::object> kv = unpack(pair, 2);
        log(kLogger, "info", fmt(kv[0]) + ":\t" + fmt(kv[1]) + "\n");
    }
}

// --------------------------------------------------------------------------
// HicFindTads

py::object the_class() {
    return this_module().attr("HicFindTads");
}

// hm.hiCMatrix(path) through the C++ matrix layer.
py::object load_matrix(py::handle path) {
    const std::string file = fs_path(path);
    hicx::ToolMatrix loaded = hicx::ToolMatrix::load(file);
    const hicx::MatrixData& data = loaded.data();
    const py::object hic = py::module_::import("hicmatrix.HiCMatrix").attr("hiCMatrix")();
    hic.attr("matrix") = hic::csr_object(data.matrix);
    const py::module_ np = numpy();
    py::object like = py::none();
    if (loaded.input_is_h5()) {
        like = py::make_tuple(np.attr("str_")(""), numpy_scalar<std::int64_t>(0), numpy_scalar<std::int64_t>(0),
                              numpy_scalar<double>(0.0));
    }
    const py::list intervals = cut_intervals_object(data.cut_intervals, like);
    hic::set_intervals(hic, intervals);
    if (data.nan_bins.empty() && loaded.input_is_h5()) {
        hic.attr("nan_bins") = float64_array({});
    } else {
        hic.attr("nan_bins") = int64_array(data.nan_bins);
    }
    hic.attr("correction_factors") =
        data.correction_factors ? py::object(float64_array(*data.correction_factors)) : py::object(py::none());
    hic.attr("distance_counts") =
        data.distance_counts ? py::object(float64_array(*data.distance_counts)) : py::object(py::none());
    return hic;
}

void HicFindTads_set_matrix(py::handle self, py::handle pMatrix, py::handle pChromosomes) {
    if (py::isinstance<py::str>(pMatrix)) {
        self.attr("hic_ma") = load_matrix(pMatrix);
    } else {
        self.attr("hic_ma") = pMatrix;
    }
    if (pChromosomes.is_none()) {
        return;
    }
    const py::object hic = self.attr("hic_ma");
    log(kLogger, "debug", "args.chromosomeOrder: " + fmt(pChromosomes));
    log(kLogger, "debug", "ma.chrBinBoundaries " + fmt(hic.attr("chrBinBoundaries")));
    const py::module_ module = this_module();
    const py::object chromosomes = module.attr("toBytes")(pChromosomes);
    py::list valid;
    py::list invalid;
    const py::object boundaries = hic.attr("chrBinBoundaries");
    for (const py::handle chrom : py::reinterpret_borrow<py::iterable>(module.attr("toString")(chromosomes))) {
        if (contains(boundaries, chrom)) {
            valid.append(chrom);
        } else {
            invalid.append(chrom);
        }
    }
    if (length(invalid) > 0) {
        log(kLogger, "warning",
            "WARNING: The following chromosome/scaffold names were not found. Please checkthe correct "
            "spelling of the chromosome names. \n");
        log(kLogger, "warning", utf8(py::str("\n").attr("join")(invalid)));
    }
    // hiCMatrix.reorderChromosomes(valid)
    const py::object names =
        py::module_::import("hicmatrix.utilities").attr("check_chrom_str_bytes")(boundaries, valid);
    std::vector<std::int64_t> order;
    for (const py::handle name : py::reinterpret_borrow<py::iterable>(names)) {
        if (!contains(boundaries, name)) {
            raise_error("ValueError", "Chromosome name '" + str_of(name) +
                                          "' not found. Please check the correct spelling of the "
                                          "chromosomes and try again");
        }
        const std::vector<py::object> range = unpack(get_item(boundaries, name), 2);
        for (std::int64_t b = to_int64(range[0]); b < to_int64(range[1]); ++b) {
            order.push_back(b);
        }
    }
    // reorderBins(order)
    const hicx::CsrMatrix m = hic::matrix(hic);
    const std::int64_t original_rows = m.rows();
    hic.attr("matrix") = hic::csr_object(hicx::select_bins(m, order));
    const py::object intervals = hic.attr("cut_intervals");
    py::list reordered;
    for (const std::int64_t b : order) {
        reordered.append(item(intervals, b));
    }
    const py::object nan_bins = hic.attr("nan_bins");
    if (length(nan_bins) > 0) {
        std::vector<char> flagged(static_cast<std::size_t>(original_rows), 0);
        for (const py::handle b : py::reinterpret_borrow<py::iterable>(nan_bins)) {
            std::int64_t at = to_int64(b);
            if (at < 0) {
                at += original_rows;
            }
            flagged[static_cast<std::size_t>(at)] = 1;
        }
        std::vector<std::int64_t> kept;
        for (std::size_t k = 0; k < order.size(); ++k) {
            if (flagged[static_cast<std::size_t>(order[k])] != 0) {
                kept.push_back(static_cast<std::int64_t>(k));
            }
        }
        hic.attr("nan_bins") = int64_array(kept);
    } else {
        hic.attr("nan_bins") = py::list();
    }
    hic::set_intervals(hic, reordered);
}

void HicFindTads_set_variables(py::handle self) {
    const Number binsize = Number::of(self.attr("binsize"));
    const auto bp = [&](const char* name, std::int64_t small, std::int64_t mid, std::int64_t large) {
        const Number n = binsize.value() < 1000 ? binsize.times(small)
                         : binsize.value() < 20000 ? binsize.times(mid)
                                                   : binsize.times(large);
        self.attr(name) = n.object();
    };
    if (self.attr("max_depth").is_none()) {
        bp("max_depth", 60, 40, 10);
    } else if (Number::of(self.attr("max_depth")).value() < binsize.times(5).value()) {
        log(kLogger, "error", "Please specify a --maxDepth that is at least 5 times larger than the matrix bin size");
        hic::python_exit(py::int_(1));
    }
    if (self.attr("min_depth").is_none()) {
        bp("min_depth", 30, 10, 5);
    } else if (Number::of(self.attr("min_depth")).value() < binsize.times(3).value()) {
        log(kLogger, "error", "Please specify a --minDepth that is at least 3 times larger than the matrix bin size");
        hic::python_exit(py::int_(1));
    }
    if (self.attr("step").is_none()) {
        self.attr("step") = binsize.value() < 1000 ? binsize.times(4).object() : binsize.times(2).object();
    } else if (Number::of(self.attr("step")).value() < binsize.value()) {
        log(kLogger, "error", "Please specify a --step that is at least the size of the matrix bin size");
        hic::python_exit(py::int_(1));
    }
    log(kLogger, "debug", "max depth:\t" + fmt(self.attr("max_depth")) + "\n");
    log(kLogger, "debug", "min depth:\t" + fmt(self.attr("min_depth")) + "\n");
    log(kLogger, "debug", "step:\t" + fmt(self.attr("step")) + "\n");
    log(kLogger, "debug", "bin size:\t" + fmt(self.attr("binsize")) + "\n");
}

void HicFindTads___init__(py::handle self, py::handle matrix, py::handle num_processors, py::handle max_depth,
                          py::handle min_depth, py::handle step, py::handle delta,
                          py::handle min_boundary_distance, py::handle use_zscore,
                          py::handle p_correct_for_multiple_testing, py::handle p_threshold_comparisons,
                          py::handle pChromosomes) {
    self.attr("set_matrix")(matrix, pChromosomes);
    if (!max_depth.is_none() && !min_depth.is_none() && compare(max_depth, min_depth, Py_LE)) {
        log(kLogger, "error", "Please check that maxDepth is larger than minDepth.");
        hic::python_exit(py::none());
    }
    self.attr("num_processors") = num_processors;
    self.attr("max_depth") = max_depth;
    self.attr("min_depth") = min_depth;
    self.attr("step") = step;
    self.attr("delta") = delta;
    self.attr("min_boundary_distance") = min_boundary_distance;
    self.attr("use_zscore") = use_zscore;
    self.attr("binsize") = hic::bin_size(self.attr("hic_ma"));
    self.attr("bedgraph_matrix") = py::none();
    self.attr("boundaries") = py::none();
    self.attr("set_variables")();
    self.attr("correct_for_multiple_testing") = p_correct_for_multiple_testing;
    self.attr("threshold_comparisons") = p_threshold_comparisons;
}

py::object HicFindTads_peakdetect(py::handle y_axis, py::handle x_axis_in, py::handle lookahead_in,
                                  py::handle delta, py::handle chrom) {
    const std::size_t n = length(y_axis);
    py::object x_axis = py::reinterpret_borrow<py::object>(x_axis_in);
    if (x_axis.is_none()) {
        std::vector<std::int64_t> range(n);
        for (std::size_t i = 0; i < n; ++i) {
            range[i] = static_cast<std::int64_t>(i);
        }
        x_axis = int64_array(range);
    }
    if (n != length(x_axis)) {
        raise_error("ValueError", "Input vectors y_axis and x_axis must have same length");
    }
    const bool scalar = PyLong_Check(delta.ptr()) || PyFloat_Check(delta.ptr()) ||
                        is_numpy_integer(delta) || is_numpy_floating(delta) ||
                        py::isinstance(delta, numpy().attr("bool_"));
    if (!(scalar && compare(delta, py::int_(0), Py_GE))) {
        raise_error("ValueError", "delta must be a positive number");
    }
    const bool y_array = is_ndarray(y_axis);
    if (y_array) {
        const py::array array = py::reinterpret_borrow<py::array>(y_axis);
        const DType kind = dtype_of(array);
        if (kind == DType::Float32 || kind == DType::Float16 || kind == DType::LongDouble) {
            raise_error("TypeError", "HiCExplorer 4: peakdetect needs a float64 or integer y_axis, got " +
                                         dtype_name(array));
        }
    }
    const std::vector<double> y = doubles_of(y_axis, "peakdetect");
    std::optional<std::vector<std::string>> names;
    if (!chrom.is_none()) {
        names = strings_of(chrom);
    }
    const ft::PeakResult peaks =
        ft::peakdetect_general(y, names ? &*names : nullptr, to_int64(lookahead_in), to_double(delta),
                               y_array, assertions_enabled(), type_name(y_axis));
    const auto build = [&](const std::vector<ft::PeakAt>& list) {
        py::list out;
        for (const ft::PeakAt& peak : list) {
            py::list entry;
            entry.append(peak.position < 0 ? py::object(py::none()) : item(x_axis, peak.position));
            entry.append(item(y_axis, peak.value_index));
            out.append(entry);
        }
        return out;
    };
    py::list result;
    result.append(build(peaks.maxima));
    result.append(build(peaks.minima));
    return result;
}

py::object HicFindTads_delta_wrt_window(py::handle min_idx_list, py::handle matrix_avg, py::handle chrom,
                                        py::handle window_len) {
    const std::vector<double> scores = doubles_of(matrix_avg, "delta_wrt_window");
    const std::vector<std::string> names = strings_of(chrom);
    std::vector<std::int64_t> minima;
    py::list keys;
    for (const py::handle m : py::reinterpret_borrow<py::iterable>(min_idx_list)) {
        minima.push_back(to_int64(m));
        keys.append(m);
    }
    const std::vector<double> values =
        ft::delta_wrt_window_values(minima, scores, names, to_int64(window_len));
    py::dict out;
    for (std::size_t i = 0; i < minima.size(); ++i) {
        // np.nan near a border, a numpy.float64 difference otherwise
        out[keys[i]] = std::isnan(values[i]) ? py::object(py::float_(kNaN))
                                             : py::object(numpy_scalar<double>(values[i]));
    }
    return out;
}

py::object HicFindTads_find_consensus_minima(py::handle tad_score_matrix, py::handle lookahead, py::handle chrom) {
    const Matrix2D m = matrix_2d(tad_score_matrix, "find_consensus_minima");
    const py::array average = float64_array(ft::row_means(m.values, m.rows, m.columns));
    const py::object cls = the_class();
    const py::object peaks =
        cls.attr("peakdetect")(average, py::arg("lookahead") = lookahead, py::arg("chrom") = chrom);
    const std::vector<py::object> both = unpack(peaks, 2);
    py::object min_idx;
    if (truthy(both[1])) {
        py::list positions;
        for (const py::handle peak : py::reinterpret_borrow<py::iterable>(both[1])) {
            positions.append(item(peak, 0));
        }
        min_idx = py::tuple(positions);
    } else {
        min_idx = py::list();
    }
    const py::object delta = cls.attr("delta_wrt_window")(min_idx, average, chrom);
    return py::make_tuple(min_idx, delta);
}

// One value as numpy holds it in the array the Python builds: np.int64 for
// an integer column, np.float64 otherwise.
py::object column_scalar(double value, bool integer) {
    if (integer) {
        return numpy_scalar<std::int64_t>(static_cast<std::int64_t>(value));
    }
    return numpy_scalar<double>(value);
}

bool all_integers(const py::list& values) {
    for (const py::handle v : values) {
        if (!is_integer(v) || PyBool_Check(v.ptr())) {
            return false;
        }
    }
    return true;
}

py::object HicFindTads_hierarchical_clustering(py::handle /*self*/, py::handle boundary_list,
                                               py::handle clusters_cutoff_in) {
    // clusters_cutoff = np.sort(np.unique(clusters_cutoff))[::-1]
    std::vector<double> cutoffs;
    std::vector<std::string> cutoff_names;
    const bool have_cutoffs = truthy(clusters_cutoff_in);
    if (have_cutoffs) {
        const py::list given = as_list(clusters_cutoff_in);
        const bool integer = all_integers(given);
        std::set<double> unique;
        for (const py::handle v : given) {
            unique.insert(to_double(v));
        }
        for (auto it = unique.rbegin(); it != unique.rend(); ++it) {
            cutoffs.push_back(*it);
            cutoff_names.push_back(integer ? std::to_string(static_cast<std::int64_t>(*it))
                                           : hicx::npy::float_repr(*it));
        }
    } else if (clusters_cutoff_in.is_none()) {
        raise_error("TypeError", "'NoneType' object is not subscriptable");
    }
    // chrom, start, value = zip(*boundary_list)
    const py::list entries = as_list(boundary_list);
    std::size_t width = std::numeric_limits<std::size_t>::max();
    std::vector<py::list> rows;
    for (const py::handle entry : entries) {
        rows.push_back(as_list(entry));
        width = std::min(width, length(rows.back()));
    }
    if (rows.empty() || width < 3) {
        raise_error("ValueError", "not enough values to unpack (expected 3, got " +
                                      std::to_string(rows.empty() ? 0 : width) + ")");
    }
    if (width > 3) {
        raise_error("ValueError", "too many values to unpack (expected 3)");
    }
    py::list chrom;
    py::list start;
    py::list value;
    for (const py::list& row : rows) {
        chrom.append(row[0]);
        start.append(row[1]);
        value.append(row[2]);
    }
    const bool start_integer = all_integers(start);
    const bool value_integer = all_integers(value);
    std::vector<double> starts;
    std::vector<double> values;
    for (const py::handle s : start) {
        starts.push_back(to_double(s));
    }
    for (const py::handle v : value) {
        values.push_back(to_double(v));
    }
    // unique_chr, indices = np.unique(chrom, return_index=True); indices[1:]
    const std::vector<std::string> names = strings_of(chrom);
    std::map<std::string, std::int64_t> first;
    for (std::size_t i = 0; i < names.size(); ++i) {
        first.emplace(names[i], static_cast<std::int64_t>(i));
    }
    std::vector<std::int64_t> split;
    for (const auto& [name, index] : first) {
        split.push_back(index);
    }
    split.erase(split.begin());
    const auto pieces = ft::numpy_split_bounds(starts.size(), split);
    const py::object str_type = numpy().attr("str_");
    py::dict z_value;
    std::vector<py::list> levels_per_chrom;
    std::size_t chrom_index = 0;
    for (const auto& [name, index] : first) {
        (void)index;
        const py::object chrom_name = str_type(to_str(name));
        ft::ClusterInput input;
        input.start.assign(starts.begin() + static_cast<std::ptrdiff_t>(pieces[chrom_index].first),
                           starts.begin() + static_cast<std::ptrdiff_t>(pieces[chrom_index].second));
        input.value.assign(values.begin() + static_cast<std::ptrdiff_t>(pieces[chrom_index].first),
                           values.begin() + static_cast<std::ptrdiff_t>(pieces[chrom_index].second));
        input.start_integer = start_integer;
        input.value_integer = value_integer;
        z_value[chrom_name] = py::list();
        const ft::ClusterResult result = ft::cluster_chromosome(input, cutoffs);
        py::list merges;
        for (const ft::ClusterMerge& merge : result.merges) {
            merges.append(py::make_tuple(merge.left, merge.right,
                                         column_scalar(input.value[merge.value_index], value_integer),
                                         merge.size, merge.x_left, merge.x_right));
        }
        z_value[chrom_name] = merges;
        // cluster_to_regions: (chrom_name, start, end) per cluster, by start
        py::list levels;
        for (const auto& level : result.levels) {
            py::list regions;
            for (const ft::ClusterRegion& region : level) {
                regions.append(py::make_tuple(
                    chrom_name, column_scalar(input.start[region.start_domain], start_integer),
                    column_scalar(input.start[region.end_domain + 1], start_integer)));
            }
            levels.append(regions);
        }
        levels_per_chrom.push_back(levels);
        ++chrom_index;
    }
    py::dict clusters;
    for (std::size_t c = 0; c < cutoff_names.size(); ++c) {
        py::list merged;
        for (const py::list& levels : levels_per_chrom) {
            if (c < length(levels)) {
                for (const py::handle region : levels[c]) {
                    merged.append(region);
                }
            }
        }
        clusters[to_str(cutoff_names[c])] = merged;
    }
    return py::make_tuple(z_value, clusters);
}

void HicFindTads_save_linkage(py::handle Z, py::handle file_name) {
    py::object handle;
    try {
        handle = py::module_::import("builtins").attr("open")(file_name, "w");
    } catch (py::error_already_set& error) {
        if (!error.matches(PyExc_OSError)) {
            throw;
        }
        log(kLogger, "error", "Can't save linkage file:\n" + fmt(file_name));
        return;
    }
    std::string text;
    std::int64_t count = 0;
    const auto flush = [&]() {
        handle.attr("write")(to_str(text));
        handle.attr("close")();
    };
    try {
        for (const py::handle pair : py::reinterpret_borrow<py::iterable>(Z.attr("items")())) {
            const std::vector<py::object> kv = unpack(pair, 2);
            for (const py::handle values : py::reinterpret_borrow<py::iterable>(kv[1])) {
                const std::vector<py::object> v = unpack(values, 6);
                ++count;
                text += fmt(kv[0]) + "\t" + fmt(python_int(v[4])) + "\t" + fmt(python_int(v[5])) + "\tclust_" +
                        std::to_string(count) + "\t" + fmt(v[2]) + "\t.\t" + fmt(v[0]) + "\t" + fmt(v[1]) +
                        "\t" + fmt(v[3]) + "\n";
            }
        }
    } catch (...) {
        flush();
        throw;
    }
    flush();
}

py::object HicFindTads_get_domains(py::handle boundary_list) {
    py::object prev_start = py::none();
    py::object prev_chrom = item(item(boundary_list, 0), 0);
    py::list domains;
    for (const py::handle entry : py::reinterpret_borrow<py::iterable>(boundary_list)) {
        const std::vector<py::object> v = unpack(entry, 3);
        if (v[1].is_none()) {
            prev_start = v[1];
            prev_chrom = v[0];
            continue;
        }
        if (compare(prev_chrom, v[0], Py_NE)) {
            prev_chrom = v[0];
            prev_start = py::none();
            continue;
        }
        domains.append(py::make_tuple(v[0], prev_start, v[1], v[2]));
        prev_start = v[1];
        prev_chrom = v[0];
    }
    return domains;
}

// json.dumps of one value of the bedgraph matrix header.
std::string json_value(py::handle value) {
    if (value.is_none()) {
        return "null";
    }
    if (PyBool_Check(value.ptr())) {
        return value.ptr() == Py_True ? "true" : "false";
    }
    if (PyLong_Check(value.ptr())) {
        return str_of(value);
    }
    if (PyFloat_Check(value.ptr())) {
        const double v = to_double(value);
        if (std::isnan(v)) {
            return "NaN";
        }
        if (std::isinf(v)) {
            return v > 0 ? "Infinity" : "-Infinity";
        }
        return hicx::npy::float_repr(v);
    }
    if (py::isinstance<py::str>(value)) {
        return utf8(py::module_::import("json").attr("dumps")(value));
    }
    raise_error("TypeError", "Object of type " + str_of(py::type::of(value).attr("__name__")) +
                                 " is not JSON serializable");
}

std::string toString_text(py::handle value) {
    return str_of(to_string(value));
}

void HicFindTads_save_bedgraph_matrix(py::handle self, py::handle outfile) {
    const std::string header = "#{\"step\":" + json_value(self.attr("step")) +
                               ",\"minDepth\":" + json_value(self.attr("min_depth")) +
                               ",\"maxDepth\":" + json_value(self.attr("max_depth")) +
                               ",\"binsize\":" + json_value(self.attr("binsize")) + "}\n";
    const py::object bedgraph = self.attr("bedgraph_matrix");
    const py::object chrom = get_item(bedgraph, py::str("chrom"));
    const py::object start = get_item(bedgraph, py::str("chr_start"));
    const py::object end = get_item(bedgraph, py::str("chr_end"));
    const py::object matrix = get_item(bedgraph, py::str("matrix"));
    std::string text = header;
    const std::size_t rows = length(chrom);
    for (std::size_t r = 0; r < rows; ++r) {
        const std::vector<double> row = doubles_of(item(matrix, static_cast<std::int64_t>(r)), "save_bedgraph_matrix");
        std::string values;
        for (std::size_t c = 0; c < row.size(); ++c) {
            if (c > 0) {
                values += '\t';
            }
            values += ft::format_fixed(row[c], 6);
        }
        text += toString_text(item(chrom, static_cast<std::int64_t>(r))) + "\t" +
                toString_text(item(start, static_cast<std::int64_t>(r))) + "\t" +
                toString_text(item(end, static_cast<std::int64_t>(r))) + "\t" + values + "\n";
    }
    hic::write_text_file(outfile, text);
}

void HicFindTads_save_clusters(py::handle clusters, py::handle file_prefix) {
    for (const py::handle pair : py::reinterpret_borrow<py::iterable>(clusters.attr("items")())) {
        const std::vector<py::object> kv = unpack(pair, 2);
        const py::object path = py::str("{}_{}.bed").attr("format")(file_prefix, kv[0]);
        const py::object handle = py::module_::import("builtins").attr("open")(path, "w");
        std::string text;
        try {
            for (const py::handle interval : py::reinterpret_borrow<py::iterable>(kv[1])) {
                const std::vector<py::object> v = unpack(interval, 3);
                text += fmt(v[0]) + "\t" + fmt(v[1]) + "\t" + fmt(v[2]) + "\t.\t0\t.\n";
            }
        } catch (...) {
            handle.attr("write")(to_str(text));
            handle.attr("close")();
            throw;
        }
        handle.attr("write")(to_str(text));
        handle.attr("close")();
    }
}

std::string py_format(py::handle value, const char* spec) {
    return utf8(py::str(py::module_::import("builtins").attr("format")(value, spec)));
}

void HicFindTads_save_domains_and_boundaries(py::handle self, py::handle prefix_in) {
    const py::object bedgraph = self.attr("bedgraph_matrix");
    const py::object chrom = get_item(bedgraph, py::str("chrom"));
    const py::object chr_start = get_item(bedgraph, py::str("chr_start"));
    const py::object chr_end = get_item(bedgraph, py::str("chr_end"));
    const Matrix2D matrix = matrix_2d(get_item(bedgraph, py::str("matrix")), "save_domains_and_boundaries");
    const py::object boundaries = self.attr("boundaries");
    const py::object min_idx_in = get_item(boundaries, py::str("min_idx"));
    const py::object delta_of_min = get_item(boundaries, py::str("delta"));
    const py::object pvalue_of_min = get_item(boundaries, py::str("pvalues"));

    const std::vector<std::string> names = strings_of(chrom);
    const std::size_t n = names.size();
    std::map<std::string, std::int64_t> first;
    for (std::size_t i = 0; i < n; ++i) {
        first.emplace(names[i], static_cast<std::int64_t>(i));
    }
    std::vector<std::int64_t> chr_start_idx;
    std::vector<std::int64_t> chr_end_idx;
    if (first.size() == 1) {
        chr_start_idx = {0};
        chr_end_idx = {static_cast<std::int64_t>(n) - 1};
    } else {
        for (const auto& [name, index] : first) {
            chr_end_idx.push_back((index == 0 ? static_cast<std::int64_t>(n) : index) - 1);
        }
        chr_start_idx = chr_end_idx;  // the same array object in the Python
    }
    std::vector<double> candidates;
    bool integral = true;
    for (const std::int64_t v : chr_start_idx) {
        candidates.push_back(static_cast<double>(v));
    }
    for (const std::int64_t v : chr_end_idx) {
        candidates.push_back(static_cast<double>(v));
    }
    std::size_t given = 0;
    for (const py::handle v : py::reinterpret_borrow<py::iterable>(min_idx_in)) {
        ++given;
        if (!is_integer(v)) {
            integral = false;
        }
        candidates.push_back(to_double(v));
    }
    if (given == 0) {
        integral = false;  // np.concatenate with an empty list gives float64
    }
    std::sort(candidates.begin(), candidates.end());
    const std::vector<double> means = ft::row_means(matrix.values, matrix.rows, matrix.columns);

    const std::string mode = py::isinstance<py::str>(self.attr("correct_for_multiple_testing"))
                                 ? utf8(self.attr("correct_for_multiple_testing"))
                                 : std::string("\x01");
    const double delta = to_double(self.attr("delta"));
    std::optional<double> threshold;
    const auto threshold_of = [&]() {
        if (!threshold.has_value()) {
            threshold = to_double(mode == "fdr" ? self.attr("pvalueFDR") : self.attr("threshold_comparisons"));
        }
        return *threshold;
    };
    std::vector<py::object> filtered;
    std::vector<std::int64_t> filtered_index;
    for (const double c : candidates) {
        const py::object key = column_scalar(c, integral);
        if (!contains(delta_of_min, key)) {
            set_item(delta_of_min, key, py::float_(kNaN));
        }
        if (mode != "fdr" && mode != "bonferroni" && mode != "None") {
            continue;
        }
        if (to_double(get_item(delta_of_min, key)) >= delta && contains(pvalue_of_min, key) &&
            to_double(get_item(pvalue_of_min, key)) <= threshold_of()) {
            filtered.push_back(key);
            filtered_index.push_back(static_cast<std::int64_t>(c));
        }
    }
    const std::string count_text = std::to_string(filtered.size());
    if (mode == "fdr") {
        log(kLogger, "info", "FDR correction. Number of boundaries for delta " + fmt(self.attr("delta")) +
                                 ", qval " + fmt(self.attr("threshold_comparisons")) + ": " + count_text);
    } else if (mode == "bonferroni") {
        log(kLogger, "info", "Bonferroni correction. Number of boundaries for delta " + fmt(self.attr("delta")) +
                                 " and pval " + fmt(self.attr("threshold_comparisons")) + ": " + count_text);
    } else {
        log(kLogger, "info", "No multiple testing correction. Number of boundaries for delta " +
                                 fmt(self.attr("delta")) + ": " + count_text + ", used threshold: " +
                                 fmt(self.attr("threshold_comparisons")));
    }

    const std::vector<std::int64_t> starts = index_list(chr_start);
    const std::vector<std::int64_t> ends = index_list(chr_end);
    const std::int64_t rows = static_cast<std::int64_t>(n);
    const auto wrap = [rows](std::int64_t i) { return i < 0 ? i + rows : i; };
    const std::string prefix = utf8(py::str(prefix_in));
    const py::module_ builtins = py::module_::import("builtins");
    const py::object boundary_file = builtins.attr("open")(to_str(prefix + "_boundaries.bed"), "w");
    py::object domains_file;
    py::object gff_file;
    try {
        domains_file = builtins.attr("open")(to_str(prefix + "_domains.bed"), "w");
        gff_file = builtins.attr("open")(to_str(prefix + "_boundaries.gff"), "w");
    } catch (...) {
        if (domains_file) {
            domains_file.attr("close")();
        }
        boundary_file.attr("close")();
        throw;
    }
    std::string boundary_text;
    std::string domains_text;
    std::string gff_text;
    const std::string delta_text = fmt(self.attr("delta"));
    const std::set<std::int64_t> chromosome_ends(chr_end_idx.begin(), chr_end_idx.end());
    std::int64_t count = 1;
    for (std::size_t position = 0; position < filtered.size(); ++position) {
        const std::int64_t index = filtered_index[position];
        if (chromosome_ends.count(index) != 0) {
            continue;
        }
        const py::object key = filtered[position];
        if (!contains(delta_of_min, key)) {
            set_item(delta_of_min, key, py::float_(kNaN));
        }
        const std::size_t at = static_cast<std::size_t>(wrap(index));
        const std::size_t before = static_cast<std::size_t>(wrap(index - 1));
        const std::int64_t right_center =
            starts[at] + static_cast<std::int64_t>(static_cast<double>(ends[at] - starts[at]) / 2.0);
        const std::int64_t left_center =
            starts[before] + static_cast<std::int64_t>(static_cast<double>(ends[before] - starts[before]) / 2.0);
        if (names[at] != names[before]) {
            continue;
        }
        const std::string name = toString_text(item(chrom, wrap(index)));
        const std::string score = ft::format_fixed(means[at], 12);
        const std::string id = py_format(key, "05d");
        boundary_text += name + "\t" + std::to_string(left_center) + "\t" + std::to_string(right_center) +
                         "\tB" + id + "\t" + score + "\t.\n";
        gff_text += name + "\tHiCExplorer\tboundary\t" + std::to_string(left_center) + "\t" +
                    std::to_string(right_center) + "\t" + score + "\t.\t.\tID=B" + id +
                    ";delta=" + ft::format_fixed(to_double(get_item(delta_of_min, key)), 12) +
                    ";pvalue=" + ft::format_fixed(to_double(get_item(pvalue_of_min, key)), 12) +
                    ";tad_sep=" + score + "\n";
        if (position + 1 == filtered.size() ||
            names[at] != names[static_cast<std::size_t>(wrap(filtered_index[position + 1]))]) {
            continue;
        }
        const std::int64_t start = starts[at];
        const std::int64_t end = starts[static_cast<std::size_t>(wrap(filtered_index[position + 1]))];
        const char* rgb = count % 2 == 0 ? "51,160,44" : "31,120,180";
        domains_text += name + "\t" + std::to_string(start) + "\t" + std::to_string(end) + "\tID_" + delta_text +
                        "_" + std::to_string(count) + "\t" + score + "\t.\t" + std::to_string(start) + "\t" +
                        std::to_string(end) + "\t" + rgb + "\n";
        ++count;
    }
    boundary_file.attr("write")(to_str(boundary_text));
    domains_file.attr("write")(to_str(domains_text));
    gff_file.attr("write")(to_str(gff_text));
    boundary_file.attr("close")();
    domains_file.attr("close")();
    gff_file.attr("close")();

    std::string score_text;
    for (std::int64_t index = 1; index < rows; ++index) {
        const std::size_t at = static_cast<std::size_t>(index);
        const std::int64_t right_center =
            starts[at] + static_cast<std::int64_t>(static_cast<double>(ends[at] - starts[at]) / 2.0);
        const std::int64_t left_center =
            starts[at - 1] + static_cast<std::int64_t>(static_cast<double>(ends[at - 1] - starts[at - 1]) / 2.0);
        if (right_center <= left_center) {
            continue;
        }
        score_text += toString_text(item(chrom, index)) + "\t" + std::to_string(left_center) + "\t" +
                      std::to_string(right_center) + "\t" + ft::format_fixed(means[at], 12) + "\n";
    }
    hic::write_text_file(to_str(prefix + "_score.bedgraph"), score_text);
}

void HicFindTads_compute_spectra_matrix(py::handle self, py::handle perchr) {
    const py::object hic = self.attr("hic_ma");
    log(kLogger, "info", "removing diagonal values\n");
    hic.attr("matrix") = hic::csr_object(ft::diagflat_zero(hic::matrix(hic)));
    mask_bins(hic, hic.attr("nan_bins"));
    const py::list orig_intervals = as_list(hic.attr("cut_intervals"));
    if (truthy(self.attr("use_zscore"))) {
        log(kLogger, "info", "Computing z-score matrix...\n");
        hicx::MatrixData data;
        data.matrix = hic::matrix(hic);
        data.cut_intervals = cut_intervals(hic.attr("cut_intervals"));
        const py::object binsize = hic::bin_size(hic);
        hicx::ObsExpOptions options;
        options.max_depth_bp = Number::of(self.attr("max_depth")).value() * 2.5;
        options.zscore = true;
        options.perchr = truthy(perchr);
        hicx::convert_to_obs_exp_matrix(data, to_int64(binsize), options);
        hic.attr("matrix") = hic::csr_object(data.matrix);
    }
    const py::object new_intervals = this_module().attr("enlarge_bins")(hic.attr("cut_intervals"));
    if (compare(new_intervals, orig_intervals, Py_NE)) {
        hic::set_intervals(hic, new_intervals);
        hic.attr("orig_cut_intervals") = new_intervals;
        py::list ids;
        for (std::size_t k = 0; k < length(new_intervals); ++k) {
            ids.append(py::int_(k));
        }
        hic.attr("orig_bin_ids") = ids;
        hic.attr("nan_bins") = py::list();
    }
    const py::object binsize_object = hic::bin_size(hic);
    const Number binsize = Number::of(binsize_object);
    const auto modulo_nonzero = [&](py::handle value) {
        const py::object rest = py::reinterpret_steal<py::object>(PyNumber_Remainder(value.ptr(), binsize_object.ptr()));
        if (!rest) {
            throw_python_error();
        }
        return compare(rest, py::int_(0), Py_NE);
    };
    if (modulo_nonzero(self.attr("min_depth"))) {
        log(kLogger, "warning", "Warning. specified *depth* is not multiple of the Hi-C matrix bin size (" +
                                    fmt(binsize_object) + ")\n");
    }
    if (modulo_nonzero(self.attr("step"))) {
        log(kLogger, "warning", "Warning. specified *step* is not multiple of the Hi-C matrix bin size (" +
                                    fmt(binsize_object) + ")\n");
    }
    self.attr("binsize") = binsize_object;
    log(kLogger, "info", "Computing TAD-separation scores...\n");
    const auto in_bins = [&](const char* name) {
        return static_cast<std::int64_t>(Number::of(self.attr(name)).value() / binsize.value());
    };
    const std::int64_t min_depth_in_bins = in_bins("min_depth");
    const std::int64_t max_depth_in_bins = in_bins("max_depth");
    const std::int64_t step_in_bins = in_bins("step");
    if (step_in_bins == 0) {
        log(kLogger, "error", "Please select a step size larger than " + fmt(binsize_object));
        hic::python_exit(py::int_(1));
    }
    const py::list steps = step_sizes(self.attr("min_depth"), self.attr("max_depth"), self.attr("step"));
    log(kLogger, "info",
        "computing spectrum for window sizes between " + std::to_string(min_depth_in_bins) + " (" +
            fmt(binsize.times(min_depth_in_bins).object()) + " bp)and " + std::to_string(max_depth_in_bins) +
            " (" + fmt(binsize.times(max_depth_in_bins).object()) + " bp) at the following window sizes " +
            std::to_string(step_in_bins) + " " + repr_of(steps) + "\n");
    if (min_depth_in_bins <= 1) {
        log(kLogger, "error",
            "ERROR\nminDepth length too small. Use a value that is at least twice as large as the bin size "
            "which is: " +
                fmt(binsize_object) + "\n");
        hic::python_exit(py::int_(0));
    }
    if (max_depth_in_bins <= 1) {
        log(kLogger, "error",
            "ERROR\nmaxDepth length too small. Use a value that is larger than the bin size which is: " +
                fmt(binsize_object) + "\n");
        hic::python_exit(py::int_(0));
    }
    {
        hicx::CsrMatrix truncated = hic::matrix(hic);
        ft::truncate_to_distance(truncated, 2 * max_depth_in_bins);
        hic.attr("matrix") = hic::csr_object(truncated);
    }
    this_module().attr("hic_ma") = hic;

    std::vector<std::int64_t> bins_to_consider;
    for (const py::handle chrom : as_list(hic.attr("chrBinBoundaries"))) {
        const std::vector<py::object> range = unpack(get_item(hic.attr("chrBinBoundaries"), chrom), 2);
        for (std::int64_t b = to_int64(range[0]); b < to_int64(range[1]); ++b) {
            bins_to_consider.push_back(b);
        }
    }
    const std::int64_t processors = to_int64(self.attr("num_processors"));
    if (processors <= 0) {
        raise_error("ValueError", "number sections must be larger than 0.");
    }
    if (processors > 1) {
        log(kLogger, "info", "Using " + std::to_string(processors) + " processors\n");
    }
    // np.array_split(bins_to_consider, num_processors)
    const std::size_t total = bins_to_consider.size();
    const std::size_t sections = static_cast<std::size_t>(processors);
    const std::size_t each = total / sections;
    const std::size_t extras = total % sections;
    const Scores scores(hic);
    const std::vector<std::int64_t> windows = windows_of(steps);
    ScoreRows all = score_rows(hic, scores, bins_to_consider, windows,
                               static_cast<unsigned int>(std::min<std::int64_t>(processors, 64)));
    // a share whose bins all fail ends the call, as zip(*[]) does in a worker
    {
        std::size_t begin = 0;
        for (std::size_t share = 0; share < sections; ++share) {
            const std::size_t size = each + (share < extras ? 1 : 0);
            bool produced = false;
            for (std::size_t k = begin; k < begin + size; ++k) {
                produced = produced || all.kept[k] != 0;
            }
            if (!produced) {
                raise_error("ValueError", "not enough values to unpack (expected 3, got 0)");
            }
            begin += size;
        }
    }
    py::list chrom_list;
    for (const py::handle c : all.chrom) {
        chrom_list.append(to_string(c));
    }
    std::vector<std::int64_t> start_values;
    std::vector<std::int64_t> end_values;
    for (const py::handle s : all.start) {
        start_values.push_back(to_int64(s));
    }
    for (const py::handle e : all.end) {
        end_values.push_back(to_int64(e));
    }
    py::dict bedgraph;
    bedgraph["chrom"] = str_array(chrom_list);
    bedgraph["chr_start"] = int64_array(start_values);
    bedgraph["chr_end"] = int64_array(end_values);
    bedgraph["matrix"] = matrix_object(all.values, all.rows, windows.size(), !all.any_float);
    self.attr("bedgraph_matrix") = bedgraph;
}

void HicFindTads_load_bedgraph_matrix(py::handle self, py::handle filename, py::handle pChromosomes) {
    const std::vector<std::string> lines = read_text_lines(filename);
    py::object parameters;
    py::list chrom_list;
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    std::vector<std::vector<double>> rows;
    const py::object float_type = py::module_::import("builtins").attr("float");
    for (const std::string& raw : lines) {
        const py::str line = to_str(raw);
        if (raw.rfind("#", 0) == 0) {
            parameters = py::module_::import("json").attr("loads")(line[py::slice(1, std::nullopt, 1)].attr("strip")());
            continue;
        }
        const py::list fields = line.attr("strip")().attr("split")("\t");
        const std::vector<py::object> first = unpack(fields[py::slice(0, 3, 1)], 3);
        if (!pChromosomes.is_none() && !contains(pChromosomes, first[0])) {
            continue;
        }
        chrom_list.append(first[0]);
        starts.push_back(to_int64(py::int_(float_type(first[1]))));
        ends.push_back(to_int64(py::int_(float_type(first[2]))));
        std::vector<double> values;
        for (std::size_t k = 3; k < length(fields); ++k) {
            values.push_back(to_double(float_type(fields[k])));
        }
        rows.push_back(std::move(values));
    }
    if (!parameters) {
        raise_error("UnboundLocalError",
                    "cannot access local variable 'parameters' where it is not associated with a value");
    }
    self.attr("min_depth") = get_item(parameters, py::str("minDepth"));
    self.attr("max_depth") = get_item(parameters, py::str("maxDepth"));
    self.attr("step") = get_item(parameters, py::str("step"));
    self.attr("binsize") = get_item(parameters, py::str("binsize"));
    if (rows.empty()) {
        raise_error("ValueError", "need at least one array to concatenate");
    }
    const std::size_t columns = rows[0].size();
    std::vector<double> flat;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        if (rows[r].size() != columns) {
            raise_error("ValueError",
                        "all the input array dimensions except for the concatenation axis must match exactly, "
                        "but along dimension 1, the array at index 0 has size " +
                            std::to_string(columns) + " and the array at index " + std::to_string(r) +
                            " has size " + std::to_string(rows[r].size()));
        }
        flat.insert(flat.end(), rows[r].begin(), rows[r].end());
    }
    py::dict bedgraph;
    bedgraph["chrom"] = str_array(chrom_list);
    bedgraph["chr_start"] = int64_array(starts);
    bedgraph["chr_end"] = int64_array(ends);
    bedgraph["matrix"] = matrix_object(flat, rows.size(), columns, false);
    self.attr("bedgraph_matrix") = bedgraph;
}

py::object HicFindTads_min_pvalue(py::handle self, py::handle min_idx) {
    log(kLogger, "info", "Computing p-values for window length: " + fmt(self.attr("min_depth")) + "\n");
    const py::object bedgraph = self.attr("bedgraph_matrix");
    const py::object chrom = get_item(bedgraph, py::str("chrom"));
    const py::object chr_start = get_item(bedgraph, py::str("chr_start"));
    const py::object chr_end = get_item(bedgraph, py::str("chr_end"));
    const py::object hic = self.attr("hic_ma");
    const py::object window_object = self.attr("min_depth");
    const std::int64_t window = whole_number(window_object, "min_depth");
    const Scores scores(hic);
    require_numeric_matrix(scores.matrix, "min_pvalue");
    const py::module_ module = this_module();
    const bool checks = assertions_enabled();

    py::list new_min_idx;
    std::vector<double> pvalues;
    std::vector<char> tested;  // the p-value came from ranksums (a numpy.float64)
    const auto values_of = [&](py::handle cut) -> std::optional<std::vector<double>> {
        const py::object block = cut_weight(hic, scores, cut, window_object, false);
        if (block.is_none()) {
            return std::nullopt;
        }
        return doubles_of(block, "min_pvalue");
    };
    for (const py::handle idx : py::reinterpret_borrow<py::iterable>(min_idx)) {
        const auto range = hic::region_bin_range(hic, scores.bins, get_item(chrom, idx),
                                                 get_item(chr_start, idx), get_item(chr_end, idx));
        if (!range.has_value()) {
            continue;
        }
        const std::int64_t matrix_idx = range->first;
        new_min_idx.append(idx);
        const std::vector<py::object> position = unpack(hic::bin_pos(hic, py::int_(matrix_idx)), 4);
        if (checks && !(str_of(to_string(get_item(chrom, idx))) == str_of(to_string(position[0])) &&
                        compare(get_item(chr_start, idx), position[1], Py_EQ) &&
                        compare(get_item(chr_end, idx), position[2], Py_EQ))) {
            raise_error("AssertionError", "");
        }
        const auto [left_idx, right_idx] = bins_at_given_distance(hic, scores, py::int_(matrix_idx), window);
        std::vector<double> left = values_of(py::int_(left_idx)).value_or(std::vector<double>{});
        std::vector<double> right = values_of(py::int_(right_idx)).value_or(std::vector<double>{});
        const std::optional<std::vector<double>> boundary = values_of(py::int_(matrix_idx));
        double pval = kNaN;
        bool from_test = false;
        if ((left.empty() && right.empty()) || !boundary.has_value() || boundary->empty() || left.empty() ||
            right.empty()) {
            pval = kNaN;
        } else {
            from_test = true;
            const double p1 = hicx::tads::rank_sum_test(*boundary, left).pvalue;
            const double p2 = hicx::tads::rank_sum_test(*boundary, right).pvalue;
            pval = p2 < p1 ? p2 : p1;  // min(pval1, pval2)
        }
        pvalues.push_back(pval);
        tested.push_back(from_test ? 1 : 0);
    }
    (void)module;
    const py::object mode = self.attr("correct_for_multiple_testing");
    py::list values;
    const bool fdr = py::isinstance<py::str>(mode) && utf8(mode) == "fdr";
    const bool bonferroni = py::isinstance<py::str>(mode) && utf8(mode) == "bonferroni";
    if (fdr) {
        for (double& p : pvalues) {
            if (std::isnan(p)) {
                p = 1.0;
            }
        }
        const std::optional<double> cutoff = ft::fdr_cutoff(pvalues, to_double(self.attr("threshold_comparisons")));
        self.attr("pvalueFDR") = cutoff.has_value() ? py::object(numpy_scalar<double>(*cutoff)) : py::object(py::int_(0));
        for (const double p : pvalues) {
            values.append(numpy_scalar<double>(p));
        }
    } else if (bonferroni) {
        const double count = static_cast<double>(pvalues.size());
        for (double& p : pvalues) {
            p *= count;
            if (!std::isnan(p) && p > 1) {
                p = 1.0;
            }
            values.append(numpy_scalar<double>(p));
        }
    } else {
        for (std::size_t k = 0; k < pvalues.size(); ++k) {
            values.append(tested[k] != 0 ? py::object(numpy_scalar<double>(pvalues[k]))
                                         : py::object(py::float_(kNaN)));
        }
    }
    py::object ordered = py::module_::import("collections").attr("OrderedDict")();
    for (std::size_t k = 0; k < length(new_min_idx); ++k) {
        set_item(ordered, new_min_idx[k], values[k]);
    }
    return ordered;
}

void HicFindTads_find_boundaries(py::handle self) {
    const py::object bedgraph = self.attr("bedgraph_matrix");
    const std::vector<std::int64_t> starts = index_list(get_item(bedgraph, py::str("chr_start")));
    const std::vector<std::int64_t> ends = index_list(get_item(bedgraph, py::str("chr_end")));
    std::vector<double> widths;
    for (std::size_t k = 0; k < starts.size(); ++k) {
        widths.push_back(static_cast<double>(ends[k] - starts[k]));
    }
    const double average_bin_size = ft::numpy_median(widths);
    if (self.attr("min_boundary_distance").is_none()) {
        self.attr("min_boundary_distance") = numpy_scalar<double>(average_bin_size * 4);
    }
    const double distance = to_double(self.attr("min_boundary_distance"));
    const double ratio = distance / average_bin_size;
    if (!std::isfinite(ratio)) {
        raise_error(std::isnan(ratio) ? "ValueError" : "OverflowError",
                    std::isnan(ratio) ? "cannot convert float NaN to integer"
                                      : "cannot convert float infinity to integer");
    }
    const std::int64_t lookahead = static_cast<std::int64_t>(ratio);
    if (lookahead < 1) {
        raise_error("ValueError", "minBoundaryDistance must be '1' or above in value");
    }
    const py::object matrix = get_item(bedgraph, py::str("matrix"));
    const py::object found = the_class().attr("find_consensus_minima")(
        matrix, py::arg("lookahead") = lookahead, py::arg("chrom") = get_item(bedgraph, py::str("chrom")));
    const std::vector<py::object> parts = unpack(found, 2);
    const py::object min_idx = parts[0];
    const py::object delta = parts[1];
    const py::object pvalues = self.attr("min_pvalue")(min_idx);
    const std::size_t found_count = length(min_idx);
    if (found_count <= 10) {
        const Matrix2D m = matrix_2d(matrix, "find_boundaries");
        const std::vector<double> means = ft::row_means(m.values, m.rows, m.columns);
        const double mean = hicx::npy::pairwise_sum(means.data(), means.size()) / static_cast<double>(means.size());
        const std::string msg =
            "Please check the parameters:\n delta: " + fmt(self.attr("delta")) + "\n minBoundaryDistance: " +
            fmt(self.attr("min_boundary_distance")) + "\n\nTAD Score values:\n mean: " + ft::format_fixed(mean, 3) +
            "\n median: " + ft::format_fixed(ft::numpy_median(means), 3) + "\n 1st quartile: " +
            ft::format_fixed(ft::numpy_percentile(means, 25), 3) + "\n 3rd quartile: " +
            ft::format_fixed(ft::numpy_percentile(means, 75), 3) + "\n";
        if (found_count == 0) {
            log(kLogger, "error", "\n*ERROR*\nNo boundaries were found. " + msg);
            hic::python_exit(py::int_(1));
        }
        log(kLogger, "info", "Only " + std::to_string(found_count) + " boundaries found. " + msg);
    }
    py::dict boundaries;
    boundaries["min_idx"] = min_idx;
    boundaries["delta"] = delta;
    boundaries["pvalues"] = pvalues;
    self.attr("boundaries") = boundaries;
}

}  // namespace

void bind_hicFindTADs(py::module_& m) {
    m.def("compute_matrix_wrapper", &compute_matrix_wrapper, py::arg("args"));
    m.def("get_cut_weight_by_bin_id", &get_cut_weight_by_bin_id, py::arg("matrix"), py::arg("cut"),
          py::arg("depth"), py::arg("return_mean") = false);
    m.def("get_idx_of_bins_at_given_distance", &get_idx_of_bins_at_given_distance, py::arg("hic_matrix"),
          py::arg("idx"), py::arg("window_len"));
    m.def("get_cut_weight", &get_cut_weight, py::arg("hic_matrix"), py::arg("cut"), py::arg("window_len"),
          py::arg("return_mean") = false);
    m.def("get_triangle", &get_triangle, py::arg("hic_matrix"), py::arg("cut"), py::arg("window_len"),
          py::arg("return_mean") = false);
    m.def("get_incremental_step_size", &get_incremental_step_size, py::arg("min_win_size"),
          py::arg("max_win_size"), py::arg("start_step_len"));
    m.def("compute_matrix", &compute_matrix, py::arg("bins_list"), py::arg("min_win_size") = 8,
          py::arg("max_win_size") = 50, py::arg("step_len") = 2);
    m.def("print_args", &print_args, py::arg("args"));

    m.def("HicFindTads___init__", &HicFindTads___init__, py::arg("self"), py::arg("matrix"),
          py::arg("num_processors"), py::arg("max_depth"), py::arg("min_depth"), py::arg("step"), py::arg("delta"),
          py::arg("min_boundary_distance"), py::arg("use_zscore"), py::arg("p_correct_for_multiple_testing"),
          py::arg("p_threshold_comparisons"), py::arg("pChromosomes"));
    m.def("HicFindTads_set_matrix", &HicFindTads_set_matrix, py::arg("self"), py::arg("pMatrix"),
          py::arg("pChromosomes"));
    m.def("HicFindTads_set_variables", &HicFindTads_set_variables, py::arg("self"));
    m.def("HicFindTads_peakdetect", &HicFindTads_peakdetect, py::arg("y_axis"), py::arg("x_axis"),
          py::arg("lookahead"), py::arg("delta"), py::arg("chrom"));
    m.def("HicFindTads_delta_wrt_window", &HicFindTads_delta_wrt_window, py::arg("min_idx_list"),
          py::arg("matrix_avg"), py::arg("chrom"), py::arg("window_len"));
    m.def("HicFindTads_find_consensus_minima", &HicFindTads_find_consensus_minima, py::arg("tad_score_matrix"),
          py::arg("lookahead"), py::arg("chrom"));
    m.def("HicFindTads_hierarchical_clustering", &HicFindTads_hierarchical_clustering, py::arg("self"),
          py::arg("boundary_list"), py::arg("clusters_cutoff"));
    m.def("HicFindTads_save_linkage", &HicFindTads_save_linkage, py::arg("Z"), py::arg("file_name"));
    m.def("HicFindTads_get_domains", &HicFindTads_get_domains, py::arg("boundary_list"));
    m.def("HicFindTads_save_bedgraph_matrix", &HicFindTads_save_bedgraph_matrix, py::arg("self"),
          py::arg("outfile"));
    m.def("HicFindTads_save_clusters", &HicFindTads_save_clusters, py::arg("clusters"), py::arg("file_prefix"));
    m.def("HicFindTads_save_domains_and_boundaries", &HicFindTads_save_domains_and_boundaries, py::arg("self"),
          py::arg("prefix"));
    m.def("HicFindTads_compute_spectra_matrix", &HicFindTads_compute_spectra_matrix, py::arg("self"),
          py::arg("perchr"));
    m.def("HicFindTads_load_bedgraph_matrix", &HicFindTads_load_bedgraph_matrix, py::arg("self"),
          py::arg("filename"), py::arg("pChromosomes"));
    m.def("HicFindTads_min_pvalue", &HicFindTads_min_pvalue, py::arg("self"), py::arg("min_idx"));
    m.def("HicFindTads_find_boundaries", &HicFindTads_find_boundaries, py::arg("self"));
}

}  // namespace hicx::python
