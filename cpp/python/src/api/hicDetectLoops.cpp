// hicexplorer._hicx.hicDetectLoops: the C++ bodies of the public functions
// of hicexplorer.hicDetectLoops.
//
// The candidate selection is hicx::loops (cpp/core/include/hicx/
// detect_loops.hpp, shared with the tool): is_neighborhood_maximum and
// test_candidate per candidate, the negative binomial fit of
// hicx::stats::fit_nbinom and the tail 1 - betainc of scipy 1.14
// (hicx::stats::betainc). What 3.7.6 did with processes and queues is done
// here in this process: the functions that put their result into a queue
// still put it there, with the 'Fail: <message><traceback>' strings the
// workers sent on an error, and the functions that started workers compute
// the same partitions' results without them. Log records a 3.7.6 worker
// process emitted never reached the parent's handlers and are not emitted
// by the functions that started them here either. The functions the Python
// calls through module globals (obs_exp_matrix, read_threshold_file,
// compute_long_range_contacts, neighborhood_merge, candidate_region_test,
// cluster_to_genome_position_mapping) are called through this module, and
// the hiCMatrix through its own methods (keepOnlyTheseChr, getBinSize,
// getBinPos). compute_loops with a queue loads the chromosome with the C++
// matrix layer (hic_object.hpp).
//
// Deviations from 3.7.6:
//  * the traceback text of a 'Fail: ' string names this module instead of
//    the Python frames of the worker (the message before it is the same);
//  * the one second sleeps of the queue polling loops are gone;
//  * candidate_region_test with a peak width larger than the window size
//    raises ValueError; 3.7.6 logged a warning in the worker, which then
//    returned without a result, and waited for it forever;
//  * with float32 obs/exp values the tail probability is computed in
//    double precision and rounded to float32, where scipy evaluates its
//    single precision betainc (the fit is the float32 one of the tool);
//  * create_distance_distribution and compute_p_values_mask take numpy
//    arrays for the obs/exp values and distances; other sequences give a
//    'Fail: ' TypeError naming HiCExplorer 4.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "hic_object.hpp"
#include "hicx/detect_loops.hpp"
#include "hicx/py_text.hpp"
#include "hicx/stats_ops.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicDetectLoops";
namespace lp = hicx::loops;

py::module_ this_module() {
    return py::module_::import("hicexplorer.hicDetectLoops");
}

// ---------------------------------------------------------------------------
// 'Fail: ' strings

struct Failure {
    std::string type;
    std::string message;
};

// The exception currently being handled, as type name and str(); SystemExit
// and KeyboardInterrupt are not Exceptions and are rethrown.
Failure current_failure() {
    try {
        throw;
    } catch (py::error_already_set& error) {
        if (!error.matches(PyExc_Exception)) {
            throw;
        }
        return {py::str(error.type().attr("__qualname__")), str_of(error.value())};
    } catch (const hicx::PythonError& error) {
        return {error.type(), error.what()};
    }
}

std::string fail_text(const Failure& failure, const char* function) {
    return "Fail: " + failure.message + "Traceback (most recent call last):\n  File \"hicexplorer._hicx.hicDetectLoops\", in " +
           function + "\n" + failure.type + ": " + failure.message + "\n";
}

bool is_fail(py::handle value) {
    if (!PyUnicode_Check(value.ptr())) {
        return false;
    }
    return utf8(value).find("Fail: ") != std::string::npos;
}

// a, b = value, with Python's unpacking errors
std::pair<py::object, py::object> unpack_two(py::handle value) {
    const py::object iterable = py::reinterpret_borrow<py::object>(value);
    std::vector<py::object> items;
    for (const py::handle item : iterable) {
        items.push_back(py::reinterpret_borrow<py::object>(item));
        if (items.size() > 2) {
            raise_error("ValueError", "too many values to unpack (expected 2)");
        }
    }
    if (items.size() < 2) {
        raise_error("ValueError", "not enough values to unpack (expected 2, got " + std::to_string(items.size()) + ")");
    }
    return {items[0], items[1]};
}

// ---------------------------------------------------------------------------
// matrices

// A scipy matrix as the canonical CSR the kernels read (duplicates summed,
// columns sorted), which is what scipy's slicing and element access see.
hicx::CsrMatrix canonical(py::handle object, const char* function) {
    const Sparse matrix = sparse(object, function);
    const hicx::utilities::SparseEntries entries = matrix.entries();
    std::vector<std::int32_t> row(entries.row.begin(), entries.row.end());
    std::vector<std::int32_t> col(entries.col.begin(), entries.col.end());
    return hicx::CsrMatrix::from_coo(matrix.rows, matrix.cols, row, col, lc::doubles(matrix.data, function),
                                     dtype_name(matrix.data));
}

lp::Candidate candidate_of(py::handle item) {
    const lp::Candidate c{to_int64(lc::item(item, 0)), to_int64(lc::item(item, 1))};
    if (c.row < 0 || c.col < 0) {
        raise_error("ValueError", "HiCExplorer 4 takes non-negative candidate indices, not (" + std::to_string(c.row) +
                                      ", " + std::to_string(c.col) + ")");
    }
    return c;
}

// matrix[row, col] with scipy's index check, for the element the Python reads
void check_element(const hicx::CsrMatrix& matrix, const lp::Candidate& c) {
    if (c.row < -matrix.rows() || c.row >= matrix.rows()) {
        raise_error("IndexError", "row index (" + std::to_string(c.row) + ") out of range");
    }
    if (c.col < -matrix.cols() || c.col >= matrix.cols()) {
        raise_error("IndexError", "column index (" + std::to_string(c.col) + ") out of range");
    }
}

// nonzero(): the coordinates of the stored values != 0 in storage order
struct NonzeroEntries {
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
};

NonzeroEntries nonzero(const Sparse& matrix) {
    const hicx::utilities::SparseEntries entries = matrix.entries();
    const std::vector<double> data = lc::doubles(matrix.data, "nonzero");
    NonzeroEntries out;
    for (std::size_t k = 0; k < data.size(); ++k) {
        if (data[k] != 0.0) {
            out.row.push_back(entries.row[k]);
            out.col.push_back(entries.col[k]);
        }
    }
    return out;
}

[[noreturn]] void boolean_index_mismatch(std::size_t dimension, std::size_t mask) {
    raise_error("IndexError", "boolean index did not match indexed array along dimension 0; dimension is " +
                                  std::to_string(dimension) + " but corresponding boolean dimension is " +
                                  std::to_string(mask));
}

// matrix.data[mask] = 0 with mask over the nonzero entries, in place
void zero_where(py::handle matrix_object, const std::vector<char>& mask) {
    Sparse matrix = sparse(matrix_object, "compute_loops");
    if (static_cast<std::size_t>(matrix.data.size()) != mask.size()) {
        boolean_index_mismatch(static_cast<std::size_t>(matrix.data.size()), mask.size());
    }
    dispatch_arithmetic(matrix.data, "compute_loops", [&]<class T>() {
        std::vector<T> data = to_vector<T>(matrix.data);
        bool changed = false;
        for (std::size_t k = 0; k < data.size(); ++k) {
            if (mask[k]) {
                data[k] = T(0);
                changed = true;
            }
        }
        if (changed) {
            write_back(matrix.data, data);
        }
    });
}

// triu(matrix, format='csr')
py::object triu_csr(py::handle matrix_object) {
    const Sparse matrix = sparse(matrix_object, "compute_loops");
    const hicx::utilities::SparseEntries entries = matrix.entries();
    const std::vector<double> data = lc::doubles(matrix.data, "compute_loops");
    std::vector<std::int32_t> row;
    std::vector<std::int32_t> col;
    std::vector<double> values;
    for (std::size_t k = 0; k < data.size(); ++k) {
        if (entries.row[k] <= entries.col[k]) {
            row.push_back(static_cast<std::int32_t>(entries.row[k]));
            col.push_back(static_cast<std::int32_t>(entries.col[k]));
            values.push_back(data[k]);
        }
    }
    return csr_matrix_object(
        hicx::CsrMatrix::from_coo(matrix.rows, matrix.cols, row, col, std::move(values), dtype_name(matrix.data)));
}

// ---------------------------------------------------------------------------
// get_linenumber, read_threshold_file, write_bedgraph

py::object get_linenumber() {
    PyFrameObject* frame = PyEval_GetFrame();
    if (frame == nullptr) {
        raise_error("AttributeError", "'NoneType' object has no attribute 'f_back'");
    }
    PyFrameObject* back = PyFrame_GetBack(frame);
    if (back == nullptr) {
        raise_error("AttributeError", "'NoneType' object has no attribute 'f_lineno'");
    }
    const int line = PyFrame_GetLineNumber(back);
    Py_DECREF(back);
    return py::int_(line);
}

py::object read_threshold_file(py::handle pFile) {
    py::dict distance_value;
    for (const std::string& raw : read_text_lines(pFile)) {
        const std::string line(hicx::pytext::strip(raw));
        if (line.rfind('#', 0) == 0) {
            continue;
        }
        if (line.empty()) {
            break;
        }
        const py::list fields = to_str(line).attr("split")("\t");
        if (fields.size() != 2) {
            raise_error("ValueError", fields.size() < 2 ? "not enough values to unpack (expected 2, got " +
                                                              std::to_string(fields.size()) + ")"
                                                        : "too many values to unpack (expected 2)");
        }
        const py::object key = parse_int(fields[0]);
        PyObject* value = PyFloat_FromString(fields[1].ptr());
        if (value == nullptr) {
            throw_python_error();
        }
        distance_value[key] = py::reinterpret_steal<py::object>(value);
    }
    return distance_value;
}

std::string percent_s(py::handle loop_item) {
    // "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % loop_item
    std::vector<py::object> values;
    if (PyTuple_Check(loop_item.ptr())) {
        for (const py::handle v : loop_item) {
            values.push_back(py::reinterpret_borrow<py::object>(v));
        }
    } else {
        values.push_back(py::reinterpret_borrow<py::object>(loop_item));
    }
    if (values.size() < 7) {
        raise_error("TypeError", "not enough arguments for format string");
    }
    if (values.size() > 7) {
        raise_error("TypeError", "not all arguments converted during string formatting");
    }
    std::string line;
    for (std::size_t k = 0; k < 7; ++k) {
        line += str_of(values[k]);
        line += k < 6 ? "\t" : "\n";
    }
    return line;
}

void write_bedgraph(py::handle pLoops, py::handle pOutFileName, py::handle pStartRegion, py::handle pEndRegion) {
    // open(pOutFileName, 'w') through Python, which gives open()'s errors and
    // the locale's encoding; the lines are formatted here
    const py::object handle = py::module_::import("builtins").attr("open")(pOutFileName, "w");
    try {
        const bool region = lc::truthy(pStartRegion) && lc::truthy(pEndRegion);
        for (const py::handle loop_item : pLoops) {
            if (region) {
                if (lc::compare(lc::item(loop_item, 1), pStartRegion, Py_GE) &&
                    lc::compare(lc::item(loop_item, 2), pEndRegion, Py_LE) &&
                    lc::compare(lc::item(loop_item, 4), pStartRegion, Py_GE) &&
                    lc::compare(lc::item(loop_item, 5), pEndRegion, Py_LE)) {
                    handle.attr("write")(to_str(percent_s(loop_item)));
                }
            } else {
                handle.attr("write")(to_str(percent_s(loop_item)));
            }
        }
    } catch (...) {
        handle.attr("close")();
        throw;
    }
    handle.attr("close")();
}

// ---------------------------------------------------------------------------
// distance distributions and p-values

struct Distributions {
    py::dict positions;
    py::dict obs_exp;
};

// create_distance_distribution without the queue
Distributions distance_distribution(py::handle pDataObsExp, py::handle pDistances, py::handle pMinDistance,
                                    py::handle pMaxDistance) {
    Distributions out;
    const std::int64_t first = to_int64(pMinDistance);
    const std::int64_t last = to_int64(pMaxDistance);
    if (first >= last) {
        return out;
    }
    if (!is_ndarray(pDataObsExp) || !is_ndarray(pDistances)) {
        raise_error("TypeError", "create_distance_distribution: HiCExplorer 4 takes numpy arrays, not '" +
                                     type_name(is_ndarray(pDataObsExp) ? pDistances : pDataObsExp) + "'");
    }
    const py::array data = py::reinterpret_borrow<py::array>(pDataObsExp);
    const std::vector<double> distances = lc::doubles(py::reinterpret_borrow<py::array>(pDistances),
                                                      "create_distance_distribution");
    for (std::int64_t distance = first; distance < last; ++distance) {
        std::vector<std::size_t> keep;
        for (std::size_t k = 0; k < distances.size(); ++k) {
            if (distances[k] == static_cast<double>(distance)) {
                keep.push_back(k);
            }
        }
        if (static_cast<std::size_t>(data.size()) != distances.size()) {
            boolean_index_mismatch(static_cast<std::size_t>(data.size()), distances.size());
        }
        out.obs_exp[py::int_(distance)] = take_elements(data, keep);
        std::vector<std::int64_t> positions(keep.begin(), keep.end());
        out.positions[py::int_(distance)] = from_vector(positions);
    }
    return out;
}

// compute_p_values_mask without the queue: the positions that pass
py::list p_values_mask(py::handle obs_exp, py::handle keys, py::handle pPValuePreselection, py::handle positions,
                       py::handle pResolution, py::handle pObsExpThreshold) {
    py::list true_values;
    if (lc::length(keys) == 0) {
        return true_values;
    }
    const bool float_threshold = PyFloat_Check(pPValuePreselection.ptr());
    const double obs_exp_threshold = to_double(pObsExpThreshold);
    for (const py::handle key : keys) {
        const py::object values_object = lc::item(obs_exp, key);
        const py::array values = py::reinterpret_steal<py::array>(
            py::module_::import("numpy").attr("array")(values_object).release());
        const DType dtype = dtype_of(values);
        const bool single = dtype == DType::Float32;
        std::vector<double> data;
        switch (dtype) {
            case DType::Int32:
            case DType::Int64:
            case DType::Float32:
            case DType::Float64:
                data = lc::doubles(values, "compute_p_values_mask");
                break;
            default:
                raise_error("TypeError", "compute_p_values_mask: HiCExplorer 4 takes int32, int64, float32 and float64 "
                                         "obs/exp values, not " + dtype_name(values));
        }
        // fit_nbinom.fit of an empty array gives size 10 and prob NaN, and then
        // nothing is tested
        hicx::stats::NBinomFit fit;
        if (!data.empty()) {
            fit = hicx::stats::fit_nbinom(
                data, single ? hicx::stats::NBinomPrecision::Float32 : hicx::stats::NBinomPrecision::Float64);
        }
        double limit = 0.0;
        if (float_threshold) {
            limit = to_double(pPValuePreselection);
        } else {
            // pPValuePreselection[int(key * pResolution)]
            PyObject* product = PyNumber_Multiply(key.ptr(), pResolution.ptr());
            if (product == nullptr) {
                throw_python_error();
            }
            const py::object genomic = py::reinterpret_steal<py::object>(product);
            PyObject* as_int = PyNumber_Long(genomic.ptr());
            if (as_int == nullptr) {
                throw_python_error();
            }
            limit = to_double(lc::item(pPValuePreselection, py::reinterpret_steal<py::object>(as_int)));
        }
        py::object key_positions;
        for (std::size_t k = 0; k < data.size(); ++k) {
            bool selected = false;
            if (single) {
                const float value = static_cast<float>(data[k]);
                if (!(value >= static_cast<float>(obs_exp_threshold))) {
                    continue;
                }
                const float cdf = static_cast<float>(hicx::stats::betainc(
                    static_cast<double>(static_cast<float>(fit.size)), static_cast<double>(value + 1.0F),
                    static_cast<double>(static_cast<float>(fit.prob))));
                selected = 1.0F - cdf <= static_cast<float>(limit);
            } else {
                if (!(data[k] >= obs_exp_threshold)) {
                    continue;
                }
                const double p = 1.0 - hicx::stats::betainc(fit.size, data[k] + 1.0, fit.prob);
                selected = p <= limit;
            }
            if (selected) {
                if (!key_positions) {
                    key_positions = lc::item(positions, key);
                }
                true_values.append(lc::item(key_positions, static_cast<std::int64_t>(k)));
            }
        }
    }
    return true_values;
}

// ---------------------------------------------------------------------------
// neighborhood_merge and candidate_region_test

// neighborhood_merge_thread without the queue
py::list neighborhood_keep(py::handle candidates, py::handle pWindowSize, const hicx::CsrMatrix& matrix) {
    py::list kept;
    if (lc::length(candidates) == 0) {
        return kept;
    }
    const std::int64_t window = to_int64(pWindowSize);
    for (const py::handle item : candidates) {
        const lp::Candidate c = candidate_of(item);
        // np.max(neighborhood) == pInteractionCountMatrix[candidate[0], candidate[1]]
        const std::int64_t start_x = (c.row - window) > 0 ? c.row - window : 0;
        const std::int64_t start_y = (c.col - window) > 0 ? c.col - window : 0;
        const std::int64_t end_x = c.row + window + 1 < matrix.rows() ? c.row + window + 1 : matrix.rows();
        const std::int64_t end_y = c.col + window + 1 < matrix.cols() ? c.col + window + 1 : matrix.cols();
        if (end_x > start_x && end_y > start_y) {
            check_element(matrix, c);
        }
        if (lp::is_neighborhood_maximum(matrix, c, window)) {
            kept.append(item);
        }
    }
    return kept;
}

struct RegionResults {
    py::list mask;
    py::list pvalues;
    bool peak_wider = false;
};

// candidate_region_test_thread without the queue
RegionResults region_test(const hicx::CsrMatrix& matrix, py::handle candidates, py::handle pWindowSize,
                          py::handle pPValue, py::handle pPeakWindowSize) {
    RegionResults out;
    if (lc::length(candidates) == 0) {
        return out;
    }
    const std::int64_t window = to_int64(pWindowSize);
    const std::int64_t peak = to_int64(pPeakWindowSize);
    const double p_value = to_double(pPValue);
    for (const py::handle item : candidates) {
        const lp::Candidate c = candidate_of(item);
        const std::int64_t start_x = (c.row - window) > 0 ? c.row - window : 0;
        const std::int64_t end_x = c.row + window + 1 < matrix.rows() ? c.row + window + 1 : matrix.rows();
        if (end_x > start_x) {
            check_element(matrix, c);
        }
        const lp::CandidateTest test = lp::test_candidate(matrix, c, window, p_value, peak);
        switch (test.outcome) {
            case lp::CandidateTest::Outcome::Rejected:
                out.mask.append(py::bool_(false));
                break;
            case lp::CandidateTest::Outcome::Accepted:
                out.mask.append(py::bool_(true));
                out.pvalues.append(numpy_scalar<double>(test.pvalue));
                break;
            case lp::CandidateTest::Outcome::NoPeak:
                raise_error("IndexError", "index 0 is out of bounds for axis 0 with size 0");
            case lp::CandidateTest::Outcome::EmptyMaximum:
                raise_error("ValueError", "zero-size array to reduction operation maximum which has no identity");
            case lp::CandidateTest::Outcome::PeakWiderThanWindow:
                log(kLogger, "warning",
                    "Neighborhood window size (" + std::to_string(window) + ") needs to be larger than peak width(" +
                        std::to_string(peak) + ").");
                out.peak_wider = true;
                return out;
        }
    }
    return out;
}

// pCandidates[i * n:(i + 1) * n] for thread i of pThreads, the last one to the end
std::vector<py::object> thread_slices(py::handle candidates, py::handle pThreads) {
    const std::int64_t threads = to_int64(pThreads);
    const std::int64_t total = static_cast<std::int64_t>(lc::length(candidates));
    if (threads == 0) {
        raise_error("ZeroDivisionError", "integer division or modulo by zero");
    }
    std::int64_t per_thread = total / threads;
    if ((total % threads != 0) && ((total < 0) != (threads < 0))) {
        --per_thread;
    }
    std::vector<py::object> slices;
    for (std::int64_t i = 0; i < threads; ++i) {
        const py::object stop = i < threads - 1 ? py::object(py::int_((i + 1) * per_thread)) : py::object(py::none());
        slices.push_back(lc::item(candidates, py::slice(py::int_(i * per_thread), stop, py::none())));
    }
    return slices;
}

py::object neighborhood_merge(py::handle pCandidates, py::handle pWindowSize, py::handle pInteractionCountMatrix,
                              py::handle pThreads) {
    log(kLogger, "debug",
        "pCandidates " + str_of(lc::item(pCandidates, py::slice(py::none(), py::int_(10), py::none()))));
    const std::vector<py::object> slices = thread_slices(pCandidates, pThreads);
    std::optional<std::string> fail;
    py::list merged;
    std::optional<hicx::CsrMatrix> matrix;
    for (const py::object& slice : slices) {
        if (lc::length(slice) == 0) {
            continue;
        }
        try {
            if (!matrix) {
                matrix = canonical(pInteractionCountMatrix, "neighborhood_merge");
            }
            for (const py::handle item : neighborhood_keep(slice, pWindowSize, *matrix)) {
                merged.append(item);
            }
        } catch (...) {
            fail = fail_text(current_failure(), "neighborhood_merge_thread");
        }
    }
    if (fail) {
        return to_str(*fail);
    }
    return py::make_tuple(merged, true);
}

py::object candidate_region_test(py::handle pHiCMatrix, py::handle pCandidates, py::handle pWindowSize,
                                 py::handle pPValue, py::handle pPeakWindowSize, py::handle pThreads) {
    log(kLogger, "debug", "candidate_region_test initial: " + std::to_string(lc::length(pCandidates)));
    if (lc::length(pCandidates) == 0) {
        return py::make_tuple(py::none(), py::none());
    }
    const py::object candidates = py::module_::import("numpy").attr("array")(pCandidates);
    const std::vector<py::object> slices = thread_slices(candidates, pThreads);
    std::optional<std::string> fail;
    std::vector<bool> mask;
    py::list pvalues;
    std::optional<hicx::CsrMatrix> matrix;
    for (const py::object& slice : slices) {
        if (lc::length(slice) == 0) {
            continue;
        }
        try {
            if (!matrix) {
                matrix = canonical(pHiCMatrix, "candidate_region_test");
            }
            const RegionResults results = region_test(*matrix, slice, pWindowSize, pPValue, pPeakWindowSize);
            if (results.peak_wider) {
                raise_error("ValueError",
                            "candidate_region_test: the peak width is larger than the window size; HiCExplorer "
                            "3.7.6 waits forever here");
            }
            for (const py::handle flag : results.mask) {
                mask.push_back(lc::truthy(flag));
            }
            for (const py::handle p : results.pvalues) {
                pvalues.append(p);
            }
        } catch (py::error_already_set& error) {
            if (error.matches(PyExc_ValueError) && utf8(error.value()).find("waits forever") != std::string::npos) {
                throw;
            }
            fail = fail_text(current_failure(), "candidate_region_test_thread");
        } catch (const hicx::PythonError& error) {
            if (std::string(error.what()).find("waits forever") != std::string::npos) {
                throw;
            }
            fail = fail_text(current_failure(), "candidate_region_test_thread");
        }
    }
    if (fail) {
        return py::make_tuple(to_str(*fail), py::none());
    }
    if (mask.empty()) {
        return py::make_tuple(py::none(), py::none());
    }
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
            keep.push_back(i);
        }
    }
    // pCandidates[mask] of the 2 dimensional candidate array: whole rows
    const py::array all = py::reinterpret_borrow<py::array>(candidates);
    py::array_t<bool> flags(static_cast<py::ssize_t>(mask.size()));
    for (std::size_t i = 0; i < mask.size(); ++i) {
        flags.mutable_data()[i] = mask[i];
    }
    const py::object selected = lc::item(all, flags);
    if (keep.empty()) {
        return py::make_tuple(py::none(), py::none());
    }
    return py::make_tuple(selected, py::module_::import("numpy").attr("array")(pvalues));
}

// ---------------------------------------------------------------------------
// compute_long_range_contacts

py::object compute_long_range_contacts(py::handle pHiCMatrix, py::handle pObsExpMatrix, py::handle pWindowSize,
                                       py::handle pPValue, py::handle pPeakWindowSize,
                                       py::handle pPValuePreselection, py::handle pMinimumInteractionsThreshold,
                                       py::handle pObsExpThreshold, py::handle pThreads) {
    const Sparse obs_exp = sparse(pObsExpMatrix, "compute_long_range_contacts");
    NonzeroEntries entries = nonzero(obs_exp);
    std::vector<std::int64_t> distance(entries.row.size());
    for (std::size_t k = 0; k < distance.size(); ++k) {
        distance[k] = std::llabs(entries.row[k] - entries.col[k]);
    }
    // pHiCMatrix.matrix.data >= pMinimumInteractionsThreshold
    const Sparse hic_matrix = sparse(pHiCMatrix.attr("matrix"), "compute_long_range_contacts");
    const std::vector<double> hic_data = lc::doubles(hic_matrix.data, "compute_long_range_contacts");
    const double minimum_interactions = to_double(pMinimumInteractionsThreshold);
    std::vector<char> hard(hic_data.size(), 0);
    for (std::size_t k = 0; k < hic_data.size(); ++k) {
        hard[k] = hic_data[k] >= minimum_interactions;
    }
    // queue = [None] * pThreads
    if (!is_integer(pThreads)) {
        raise_error("TypeError", "can't multiply sequence by non-int of type '" + type_name(pThreads) + "'");
    }
    const std::int64_t threads = to_int64(pThreads);
    if (distance.empty()) {
        raise_error("ValueError", "zero-size array to reduction operation minimum which has no identity");
    }
    const std::int64_t min_distance = *std::min_element(distance.begin(), distance.end());
    const std::int64_t max_distance = *std::max_element(distance.begin(), distance.end());
    if (threads == 0) {
        raise_error("ZeroDivisionError", "integer division or modulo by zero");
    }
    if (threads < 0) {
        // no worker is started, so no distance has a distribution and nothing
        // is selected
        if (PyObject_DelAttrString(pHiCMatrix.ptr(), "matrix") != 0) {
            throw_python_error();
        }
        if (distance.size() != hard.size()) {
            raise_error("ValueError", "operands could not be broadcast together with shapes (" +
                                          std::to_string(distance.size()) + ",) (" + std::to_string(hard.size()) + ",) ");
        }
        return py::make_tuple(py::none(), py::none());
    }
    // the distance distributions, one range per thread in 3.7.6
    std::optional<std::string> fail;
    Distributions distributions;
    {
        const py::array distances = from_vector(distance);
        try {
            distributions = distance_distribution(obs_exp.data, distances, py::int_(min_distance),
                                                  py::int_(max_distance + 1));
        } catch (...) {
            fail = fail_text(current_failure(), "create_distance_distribution");
        }
    }
    if (PyObject_DelAttrString(pHiCMatrix.ptr(), "matrix") != 0) {
        throw_python_error();
    }
    if (fail) {
        return py::make_tuple(to_str(*fail), py::none());
    }
    // the p-values, the sorted keys in chunks of len // pThreads + 1
    py::list keys = py::list(distributions.obs_exp);
    keys.attr("sort")();
    const std::size_t chunk = keys.size() / static_cast<std::size_t>(threads) + 1;
    std::vector<char> mask(distance.size(), 0);
    py::object resolution = py::none();
    for (std::int64_t i = 0; i < threads; ++i) {
        const std::size_t begin = static_cast<std::size_t>(i) * chunk;
        const std::size_t end = i < threads - 1 ? std::min(keys.size(), begin + chunk) : keys.size();
        if (begin >= end) {
            continue;
        }
        if (resolution.is_none()) {
            resolution = pHiCMatrix.attr("getBinSize")();
        }
        py::list part;
        for (std::size_t k = begin; k < end; ++k) {
            part.append(keys[k]);
        }
        try {
            for (const py::handle position : p_values_mask(distributions.obs_exp, part, pPValuePreselection,
                                                           distributions.positions, resolution, pObsExpThreshold)) {
                mask[static_cast<std::size_t>(to_int64(position))] = 1;
            }
        } catch (...) {
            fail = fail_text(current_failure(), "compute_p_values_mask");
        }
    }
    if (fail) {
        return py::make_tuple(to_str(*fail), py::none());
    }
    if (mask.size() != hard.size()) {
        raise_error("ValueError", "operands could not be broadcast together with shapes (" + std::to_string(mask.size()) +
                                      ",) (" + std::to_string(hard.size()) + ",) ");
    }
    py::list candidate_rows;
    std::vector<std::int32_t> flat;
    for (std::size_t k = 0; k < mask.size(); ++k) {
        if (mask[k] && hard[k]) {
            flat.push_back(static_cast<std::int32_t>(entries.row[k]));
            flat.push_back(static_cast<std::int32_t>(entries.col[k]));
        }
    }
    if (flat.empty()) {
        return py::make_tuple(py::none(), py::none());
    }
    const py::object candidates =
        from_vector(flat).attr("reshape")(py::make_tuple(static_cast<py::ssize_t>(flat.size() / 2), 2));
    const py::module_ module = this_module();
    const auto [merged, flag] =
        unpack_two(module.attr("neighborhood_merge")(candidates, pWindowSize, pObsExpMatrix, pThreads));
    (void)flag;
    if (is_fail(merged)) {
        return py::make_tuple(merged, py::none());
    }
    if (lc::length(merged) == 0) {
        return py::make_tuple(py::none(), py::none());
    }
    const auto [tested, pvalues] = unpack_two(module.attr("candidate_region_test")(
        pObsExpMatrix, merged, pWindowSize, pPValue, pPeakWindowSize, pThreads));
    return py::make_tuple(tested, pvalues);
}

// ---------------------------------------------------------------------------
// cluster_to_genome_position_mapping

py::object cluster_to_genome_position_mapping(py::handle pHicMatrix, py::handle pCandidates, py::handle pPValueList,
                                              py::handle pMaxLoopDistance) {
    py::list mapped;
    std::int64_t i = 0;
    for (const py::handle candidate : pCandidates) {
        const py::object x = pHicMatrix.attr("getBinPos")(lc::item(candidate, 0));
        const py::object y = pHicMatrix.attr("getBinPos")(lc::item(candidate, 1));
        auto four = [](py::handle position) {
            std::vector<py::object> fields;
            for (const py::handle v : position) {
                fields.push_back(py::reinterpret_borrow<py::object>(v));
                if (fields.size() > 4) {
                    raise_error("ValueError", "too many values to unpack (expected 4)");
                }
            }
            if (fields.size() < 4) {
                raise_error("ValueError",
                            "not enough values to unpack (expected 4, got " + std::to_string(fields.size()) + ")");
            }
            return fields;
        };
        const std::vector<py::object> fx = four(x);
        const std::vector<py::object> fy = four(y);
        const std::int64_t distance = std::llabs(to_int64(py::int_(fx[1])) - to_int64(py::int_(fy[1])));
        if (!pMaxLoopDistance.is_none() && lc::compare(py::int_(distance), pMaxLoopDistance, Py_GT)) {
            ++i;
            continue;
        }
        mapped.append(py::make_tuple(fx[0], fx[1], fx[2], fy[0], fy[1], fy[2], lc::item(pPValueList, i)));
        ++i;
    }
    return mapped;
}

// ---------------------------------------------------------------------------
// compute_loops

py::object compute_loops(py::handle hic_object, py::handle pRegion, py::handle pArgs, py::handle pIsCooler,
                         py::handle pQueue) {
    const py::module_ module = this_module();
    const bool queue = !pQueue.is_none();
    py::object pHiCMatrix = py::reinterpret_borrow<py::object>(hic_object);
    py::object mapped_loops;
    try {
        const bool cooler = lc::truthy(pIsCooler);
        if (queue) {
            py::list region_list;
            region_list.append(pRegion);
            const py::object distance_object = pArgs.attr("maxLoopDistance");
            std::optional<std::int64_t> distance;
            if (!distance_object.is_none()) {
                distance = to_int64(distance_object);
            }
            const py::object matrix_path = pArgs.attr("matrix");
            const hic::Loaded loaded = hic::load(fs_path(matrix_path), utf8(pRegion), cooler ? distance : std::nullopt);
            pHiCMatrix = hic::object(loaded, matrix_path, region_list, !cooler);
        }
        if (!cooler) {
            py::list region_list;
            region_list.append(pRegion);
            pHiCMatrix.attr("keepOnlyTheseChr")(region_list);
            PyObject* quotient =
                PyNumber_TrueDivide(pArgs.attr("maxLoopDistance").ptr(), pHiCMatrix.attr("getBinSize")().ptr());
            if (quotient == nullptr) {
                throw_python_error();
            }
            const double max_loop_distance = to_double(py::reinterpret_steal<py::object>(quotient));
            const NonzeroEntries entries = nonzero(sparse(pHiCMatrix.attr("matrix"), "compute_loops"));
            std::vector<char> mask(entries.row.size(), 0);
            for (std::size_t k = 0; k < mask.size(); ++k) {
                mask[k] = static_cast<double>(std::llabs(entries.row[k] - entries.col[k])) > max_loop_distance;
            }
            zero_where(pHiCMatrix.attr("matrix"), mask);
            lc::eliminate_zeros(pHiCMatrix.attr("matrix"), "compute_loops");
        }
        if (lc::length(pHiCMatrix.attr("matrix").attr("data")) == 0) {
            py::list none;
            none.append(py::none());
            pQueue.attr("put")(none);
            return py::none();
        }
        const py::tuple shape = pHiCMatrix.attr("matrix").attr("shape");
        auto no_loops = [&]() -> py::object {
            if (!queue) {
                return py::none();
            }
            py::list none;
            none.append(py::none());
            pQueue.attr("put")(none);
            return py::none();
        };
        if (to_int64(shape[0]) < 5 || to_int64(shape[1]) < 5) {
            log(kLogger, "debug", "Computed loops for " + str_of(pRegion) + ": 0");
            return no_loops();
        }
        if (pArgs.attr("windowSize").is_none()) {
            const py::object bin_size = pHiCMatrix.attr("getBinSize")();
            const double b = to_double(bin_size);
            std::int64_t window = 5;
            if (0 < b && b <= 5000) {
                window = 10;
            }
            pArgs.attr("windowSize") = py::int_(window);
            log(kLogger, "debug", "Setting window size to: " + str_of(pArgs.attr("windowSize")));
        }
        if (pArgs.attr("peakWidth").is_none()) {
            PyObject* width = PyNumber_Subtract(pArgs.attr("windowSize").ptr(), py::int_(3).ptr());
            if (width == nullptr) {
                throw_python_error();
            }
            pArgs.attr("peakWidth") = py::reinterpret_steal<py::object>(width);
        }
        log(kLogger, "debug", "Setting peak width to: " + str_of(pArgs.attr("peakWidth")));
        pHiCMatrix.attr("matrix") = triu_csr(pHiCMatrix.attr("matrix"));
        lc::eliminate_zeros(pHiCMatrix.attr("matrix"), "compute_loops");
        {
            const NonzeroEntries entries = nonzero(sparse(pHiCMatrix.attr("matrix"), "compute_loops"));
            std::vector<char> mask(entries.row.size(), 0);
            for (std::size_t k = 0; k < mask.size(); ++k) {
                mask[k] = entries.row[k] == entries.col[k];
            }
            zero_where(pHiCMatrix.attr("matrix"), mask);
            lc::eliminate_zeros(pHiCMatrix.attr("matrix"), "compute_loops");
        }
        const py::object expected = pArgs.attr("expected");
        const py::object threads = pArgs.attr("threadsPerChromosome");
        py::object obs_exp;
        bool bound = false;
        if (lc::equal(expected, py::str("mean"))) {
            obs_exp = module.attr("obs_exp_matrix")(pHiCMatrix.attr("matrix"), py::arg("pInplace") = false,
                                                    py::arg("pToEpsilon") = true, py::arg("pThreads") = threads);
            bound = true;
        } else if (lc::equal(expected, py::str("mean_nonzero"))) {
            obs_exp = module.attr("obs_exp_matrix_non_zero")(pHiCMatrix.attr("matrix"), py::arg("ligation_factor") = false,
                                                             py::arg("pInplace") = false, py::arg("pToEpsilon") = true,
                                                             py::arg("pThreads") = threads);
            bound = true;
        } else if (lc::equal(expected, py::str("mean_nonzero_ligation"))) {
            obs_exp = module.attr("obs_exp_matrix_non_zero")(pHiCMatrix.attr("matrix"), py::arg("ligation_factor") = true,
                                                             py::arg("pInplace") = false, py::arg("pToEpsilon") = true,
                                                             py::arg("pThreads") = threads);
            bound = true;
        }
        if (!bound) {
            raise_error("UnboundLocalError",
                        "cannot access local variable 'obs_exp_csr_matrix' where it is not associated with a value");
        }
        if (!py::isinstance(obs_exp, py::module_::import("scipy.sparse").attr("csr_matrix"))) {
            return no_loops();
        }
        lc::eliminate_zeros(pHiCMatrix.attr("matrix"), "compute_loops");
        lc::eliminate_zeros(obs_exp, "compute_loops");
        if (lc::length(pHiCMatrix.attr("matrix").attr("data")) != lc::length(obs_exp.attr("data"))) {
            return no_loops();
        }
        {
            PyObject* as_float = PyNumber_Float(pArgs.attr("pValuePreselection").ptr());
            if (as_float != nullptr) {
                pArgs.attr("pValuePreselection") = py::reinterpret_steal<py::object>(as_float);
            } else {
                PyErr_Clear();
                pArgs.attr("pValuePreselection") = module.attr("read_threshold_file")(pArgs.attr("pValuePreselection"));
            }
        }
        const auto [candidates, pvalues] = unpack_two(module.attr("compute_long_range_contacts")(
            pHiCMatrix, obs_exp, pArgs.attr("windowSize"), pArgs.attr("pValue"), pArgs.attr("peakWidth"),
            pArgs.attr("pValuePreselection"), pArgs.attr("peakInteractionsThreshold"),
            pArgs.attr("obsExpThreshold"), threads));
        if (candidates.is_none()) {
            log(kLogger, "info", "Computed loops for " + str_of(pRegion) + ": 0");
            return no_loops();
        }
        if (PyUnicode_Check(candidates.ptr()) && is_fail(candidates)) {
            if (queue) {
                pQueue.attr("put")(candidates);
                return py::none();
            }
            return candidates;
        }
        mapped_loops = module.attr("cluster_to_genome_position_mapping")(pHiCMatrix, candidates, pvalues,
                                                                        pArgs.attr("maxLoopDistance"));
        log(kLogger, "debug", "Computed loops for " + str_of(pRegion) + ": " + std::to_string(lc::length(mapped_loops)));
    } catch (...) {
        const std::string text = fail_text(current_failure(), "compute_loops");
        if (queue) {
            pQueue.attr("put")(to_str(text));
            return py::none();
        }
        return to_str(text);
    }
    if (!queue) {
        return mapped_loops;
    }
    py::list result;
    result.append(mapped_loops);
    pQueue.attr("put")(result);
    return py::none();
}

}  // namespace

void bind_hicDetectLoops(py::module_& m) {
    m.def("get_linenumber", &get_linenumber);
    m.def("create_distance_distribution", [](py::handle pDataObsExp, py::handle pDistances, py::handle pWindowSize,
                                             py::handle pMinDistance, py::handle pMaxDistance, py::handle pQueue) {
        (void)pWindowSize;
        py::object result;
        try {
            const Distributions d = distance_distribution(pDataObsExp, pDistances, pMinDistance, pMaxDistance);
            py::list pair;
            pair.append(d.positions);
            pair.append(d.obs_exp);
            result = pair;
        } catch (...) {
            result = to_str(fail_text(current_failure(), "create_distance_distribution"));
        }
        pQueue.attr("put")(result);
    }, py::arg("pDataObsExp"), py::arg("pDistances"), py::arg("pWindowSize"), py::arg("pMinDistance"),
       py::arg("pMaxDistance"), py::arg("pQueue"));
    m.def("compute_p_values_mask", [](py::handle obs_exp, py::handle keys, py::handle threshold, py::handle positions,
                                      py::handle resolution, py::handle minimum, py::handle obs_exp_threshold,
                                      py::handle pQueue) {
        (void)minimum;
        py::object result;
        try {
            result = p_values_mask(obs_exp, keys, threshold, positions, resolution, obs_exp_threshold);
        } catch (...) {
            result = to_str(fail_text(current_failure(), "compute_p_values_mask"));
        }
        pQueue.attr("put")(result);
    }, py::arg("pGenomicDistanceDistributionsObsExp"), py::arg("pGenomicDistanceDistributionsKeyList"),
       py::arg("pPValuePreselection"), py::arg("pGenomicDistanceDistributionPosition"), py::arg("pResolution"),
       py::arg("pMinimumInteractionsThreshold"), py::arg("pObsExpThreshold"), py::arg("pQueue"));
    m.def("compute_long_range_contacts", &compute_long_range_contacts, py::arg("pHiCMatrix"),
          py::arg("pObsExpMatrix"), py::arg("pWindowSize"), py::arg("pPValue"), py::arg("pPeakWindowSize"),
          py::arg("pPValuePreselection"), py::arg("pMinimumInteractionsThreshold"), py::arg("pObsExpThreshold"),
          py::arg("pThreads"));
    m.def("neighborhood_merge_thread", [](py::handle pCandidateList, py::handle pWindowSize,
                                          py::handle pInteractionCountMatrix, py::handle pQueue) {
        py::object result;
        try {
            if (lc::length(pCandidateList) == 0) {
                result = py::list();
            } else {
                result = neighborhood_keep(pCandidateList, pWindowSize,
                                           canonical(pInteractionCountMatrix, "neighborhood_merge_thread"));
            }
        } catch (...) {
            result = to_str(fail_text(current_failure(), "neighborhood_merge_thread"));
        }
        pQueue.attr("put")(result);
    }, py::arg("pCandidateList"), py::arg("pWindowSize"), py::arg("pInteractionCountMatrix"), py::arg("pQueue"));
    m.def("neighborhood_merge", &neighborhood_merge, py::arg("pCandidates"), py::arg("pWindowSize"),
          py::arg("pInteractionCountMatrix"), py::arg("pThreads"));
    m.def("candidate_region_test_thread", [](py::handle pHiCMatrix, py::handle pCandidates, py::handle pWindowSize,
                                             py::handle pPValue, py::handle pPeakWindowSize,
                                             py::handle pQueue) -> py::object {
        py::object result;
        try {
            if (lc::length(pCandidates) == 0) {
                py::list pair;
                pair.append(py::list());
                pair.append(py::list());
                result = pair;
            } else {
                const RegionResults r = region_test(canonical(pHiCMatrix, "candidate_region_test_thread"), pCandidates,
                                                    pWindowSize, pPValue, pPeakWindowSize);
                if (r.peak_wider) {
                    return py::make_tuple(py::none(), py::none());
                }
                py::list pair;
                pair.append(r.mask);
                pair.append(r.pvalues);
                result = pair;
            }
        } catch (...) {
            result = to_str(fail_text(current_failure(), "candidate_region_test_thread"));
        }
        pQueue.attr("put")(result);
        return py::none();
    }, py::arg("pHiCMatrix"), py::arg("pCandidates"), py::arg("pWindowSize"), py::arg("pPValue"),
       py::arg("pPeakWindowSize"), py::arg("pQueue"));
    m.def("candidate_region_test", &candidate_region_test, py::arg("pHiCMatrix"), py::arg("pCandidates"),
          py::arg("pWindowSize"), py::arg("pPValue"), py::arg("pPeakWindowSize"), py::arg("pThreads"));
    m.def("cluster_to_genome_position_mapping", &cluster_to_genome_position_mapping, py::arg("pHicMatrix"),
          py::arg("pCandidates"), py::arg("pPValueList"), py::arg("pMaxLoopDistance"));
    m.def("write_bedgraph", &write_bedgraph, py::arg("pLoops"), py::arg("pOutFileName"),
          py::arg("pStartRegion") = py::none(), py::arg("pEndRegion") = py::none());
    m.def("compute_loops", &compute_loops, py::arg("pHiCMatrix"), py::arg("pRegion"), py::arg("pArgs"),
          py::arg("pIsCooler"), py::arg("pQueue") = py::none());
    m.def("read_threshold_file", &read_threshold_file, py::arg("pFile"));
}

}  // namespace hicx::python
