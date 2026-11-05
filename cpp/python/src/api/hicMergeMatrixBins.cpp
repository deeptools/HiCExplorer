// hicexplorer._hicx.hicMergeMatrixBins: the C++ bodies of
// remove_nans_if_needed, running_window_merge and merge_bins.
//
// The functions take a hiCMatrix and change it in place, as in 3.7.6. What
// the Python computes is computed here: the running window
// (hicx::coo_ops::running_window, with scipy's summation order), the bin
// groups and coverages of merge_bins, the nan_bins recomputed from the
// column sums, eliminate_zeros. What the Python asks the hiCMatrix to do is
// asked of it through its own methods, so a subclass keeps working:
// maskBins for the NaN bins and setCutIntervals for the merged bin table.
// merge_bins reduces the matrix through this module's reduce_matrix
// (hicexplorer.reduceMatrix.reduce_matrix, also C++), as the Python does.
//
// Deviations: matrix dtypes other than int32, int64, float32 and float64,
// and a num_bins that is not an integer in merge_bins, raise TypeError naming
// HiCExplorer 4.

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/coo_ops.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicMergeMatrixBins";

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

py::object remove_nans_if_needed(py::object hic) {
    const py::object nan_bins = hic.attr("nan_bins");
    if (!nan_bins.is_none() && lc::length(nan_bins) > 0) {
        hic.attr("maskBins")(hic.attr("nan_bins"));
        hic.attr("orig_bin_ids") = py::list();
        hic.attr("orig_cut_intervals") = py::list();
        hic.attr("correction_factors") = py::none();
        log(kLogger, "warning",
            "*WARNING*: The matrix is probably a corrected matrix that contains NaN bins. This bins can not be "
            "merged and are removed. It is preferable to first merge bins in a uncorrected  matrix and then "
            "correct the matrix. Correction factors, if present, are removed as well.");
    }
    return hic;
}

// hic.nan_bins = np.flatnonzero(hic.matrix.sum(0).A == 0)
void recompute_nan_bins(py::handle hic) {
    const Sparse matrix = sparse(hic.attr("matrix"), "merge_bins", "r");
    py::object nan_bins;
    dispatch_arithmetic(matrix.data, "merge_bins", [&]<class T>() {
        nan_bins = from_vector(hicx::coo_ops::zero_columns(lc::csr_of<T>(matrix)));
    });
    hic.attr("nan_bins") = nan_bins;
}

py::object running_window_merge(py::object hic_matrix, py::handle num_bins) {
    hic_matrix = remove_nans_if_needed(hic_matrix);
    if (lc::equal(num_bins, py::int_(1))) {
        return hic_matrix;
    }
    if (assertions_enabled()) {
        PyObject* remainder = PyNumber_Remainder(num_bins.ptr(), py::int_(2).ptr());
        if (remainder == nullptr) {
            throw_python_error();
        }
        const py::object rest = py::reinterpret_steal<py::object>(remainder);
        if (!lc::equal(rest, py::int_(1))) {
            raise_error("AssertionError", "num_bins has to be an odd number");
        }
    }
    const Sparse matrix = sparse(hic_matrix.attr("matrix"), "running_window_merge", "r");
    const std::int64_t bins = to_int64(num_bins);
    py::object merged;
    dispatch_arithmetic(matrix.data, "running_window_merge", [&]<class T>() {
        merged = lc::csr_object(hicx::coo_ops::running_window(lc::csr_of<T>(matrix), bins));
    });
    hic_matrix.attr("matrix") = merged;
    recompute_nan_bins(hic_matrix);
    lc::eliminate_zeros(hic_matrix.attr("matrix"), "running_window_merge");
    return hic_matrix;
}

py::object merge_bins(py::object hic, py::handle num_bins_object) {
    hic = remove_nans_if_needed(hic);
    const py::list intervals = py::list(hic.attr("cut_intervals"));
    std::vector<py::object> names;
    std::vector<py::object> starts;
    std::vector<py::object> ends;
    std::vector<py::object> coverages;
    std::size_t fields = intervals.empty() ? 0 : std::numeric_limits<std::size_t>::max();
    for (const py::handle interval : intervals) {
        fields = std::min(fields, lc::length(interval));
    }
    if (fields < 4) {
        raise_error("ValueError", "not enough values to unpack (expected 4, got " + std::to_string(fields) + ")");
    }
    if (fields > 4) {
        raise_error("ValueError", "too many values to unpack (expected 4)");
    }
    for (const py::handle interval : intervals) {
        names.push_back(lc::item(interval, 0));
        starts.push_back(lc::item(interval, 1));
        ends.push_back(lc::item(interval, 2));
        coverages.push_back(lc::item(interval, 3));
    }
    if (!is_integer(num_bins_object) || PyBool_Check(num_bins_object.ptr())) {
        raise_error("TypeError", "merge_bins: HiCExplorer 4 takes an integer num_bins, not '" +
                                     type_name(num_bins_object) + "'");
    }
    const std::int64_t num_bins = to_int64(num_bins_object);
    auto mean = [&](std::size_t first, std::size_t last) {
        std::vector<double> values;
        for (std::size_t i = first; i < last; ++i) {
            values.push_back(to_double(coverages[i]));
        }
        const double value = values.empty() ? std::nan("")
                                            : hicx::npy::pairwise_sum(values) / static_cast<double>(values.size());
        return numpy_scalar<double>(value);
    };
    py::list new_bins;
    py::list bins_to_merge;
    py::object prev_ref = names[0];
    std::size_t idx_start = 0;
    py::object new_start = starts[0];
    std::int64_t count = 0;
    std::size_t idx = 0;
    py::object ref;
    for (idx = 0; idx < names.size(); ++idx) {
        ref = names[idx];
        // (count > 0 and count % num_bins == 0) or ref != prev_ref
        bool cut = false;
        if (count > 0) {
            if (num_bins == 0) {
                raise_error("ZeroDivisionError", "integer modulo by zero");
            }
            std::int64_t rest = count % num_bins;
            if (rest != 0 && ((rest < 0) != (num_bins < 0))) {
                rest += num_bins;  // Python's modulo takes the sign of the divisor
            }
            cut = rest == 0;
        }
        if (!cut) {
            cut = lc::compare(ref, prev_ref, Py_NE);
        }
        if (cut) {
            if (static_cast<double>(count) < static_cast<double>(num_bins) / 2.0) {
                log(kLogger, "debug",
                    str_of(prev_ref) + " has few bins (" + std::to_string(count) + "). Skipping it\n");
            } else {
                new_bins.append(py::make_tuple(names[idx_start], new_start, ends[idx - 1], mean(idx_start, idx)));
                py::list members;
                for (std::size_t i = idx_start; i < idx; ++i) {
                    members.append(py::int_(i));
                }
                bins_to_merge.append(members);
            }
            idx_start = idx;
            new_start = starts[idx];
            count = 0;
        }
        prev_ref = ref;
        count += 1;
    }
    idx = names.size() - 1;
    new_bins.append(py::make_tuple(ref, new_start, ends[idx], mean(idx_start, names.size())));
    py::list members;
    for (std::size_t i = idx_start; i <= idx; ++i) {
        members.append(py::int_(i));
    }
    bins_to_merge.append(members);
    const py::object reduce = py::module_::import("hicexplorer.hicMergeMatrixBins").attr("reduce_matrix");
    hic.attr("matrix") = reduce(hic.attr("matrix"), bins_to_merge, py::arg("diagonal") = true);
    lc::eliminate_zeros(hic.attr("matrix"), "merge_bins");
    hic.attr("setCutIntervals")(new_bins);
    recompute_nan_bins(hic);
    return hic;
}

}  // namespace

void bind_hicMergeMatrixBins(py::module_& m) {
    m.def("remove_nans_if_needed", &remove_nans_if_needed, py::arg("hic"));
    m.def("running_window_merge", &running_window_merge, py::arg("hic_matrix"), py::arg("num_bins"));
    m.def("merge_bins", &merge_bins, py::arg("hic"), py::arg("num_bins"));
}

}  // namespace hicx::python
