// hicexplorer._hicx.iterativeCorrection: the C++ body of
// hicexplorer.iterativeCorrection.iterativeCorrection.
//
// The balancing is hicx::iterative_correction (cpp/core/include/hicx/
// iterative_correction.hpp) on the COO of the matrix; this file reproduces
// the side effects of the Python: the NaN values of the argument's data set
// to zero in place, the logger level set to INFO with verbose, the log
// records, ValueError for an asymmetric matrix and exit(1) (SystemExit) for
// the two large value guards.
//
// Deviations: the symmetry measure is summed in a different order than
// scipy's (it only matters within an ulp of the 1e-10 threshold); scipy
// sparse formats other than csr, csc and coo raise TypeError.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicx/iterative_correction.hpp"
#include "hicx/sparse_matrix.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace ic = hicx::iterative_correction;

namespace {

constexpr const char* kLogger = "hicexplorer.iterativeCorrection";

bool has_nan_sum(const py::array& data) {
    const char kind = data.dtype().kind();
    if (kind != 'f' && kind != 'c') {
        return false;
    }
    const std::vector<double> values = lc::doubles(data, "iterativeCorrection");
    bool positive = false;
    bool negative = false;
    for (const double v : values) {
        if (std::isnan(v)) {
            return true;
        }
        positive = positive || v == INFINITY;
        negative = negative || v == -INFINITY;
    }
    return positive && negative;
}

void zero_nans_in_place(Sparse& matrix) {
    py::array data = matrix.data;
    dispatch_numeric(data, "iterativeCorrection", [&]<class T>() {
        if constexpr (std::is_floating_point_v<T>) {
            std::vector<T> values = to_vector<T>(data);
            bool changed = false;
            for (T& v : values) {
                if (std::isnan(v)) {
                    v = T(0);
                    changed = true;
                }
            }
            if (changed) {
                write_back(data, values);
            }
        }
    });
}

std::string format_estimate(std::int64_t iternum, std::int64_t M, double elapsed) {
    const double estimated = (static_cast<double>(M - iternum) * elapsed) / static_cast<double>(iternum);
    const double minutes_total = std::floor(estimated / 60.0);
    const double sec = std::fmod(estimated, 60.0);
    const double h = std::floor(minutes_total / 60.0);
    const double m = std::fmod(minutes_total, 60.0);
    char buffer[128];
    std::snprintf(buffer, sizeof buffer, "pass %lld Estimated time %.0f:%.0f:%.0f",
                  static_cast<long long>(iternum), h, m, sec);
    return buffer;
}

}  // namespace

void bind_iterativeCorrection(py::module_& m) {
    m.def("iterativeCorrection", [](py::handle matrix, py::handle /*v*/, py::handle M_object,
                                    py::handle tolerance_object, py::handle verbose_object) -> py::object {
        const bool verbose = lc::truthy(verbose_object);
        if (verbose) {
            py::module_::import("logging").attr("getLogger")(kLogger).attr("setLevel")(20);
        }
        Sparse input = sparse(matrix, "iterativeCorrection");
        if (has_nan_sum(input.data)) {
            log(kLogger, "warning",
                "[iterative correction] the matrix contains nans, they will be replaced by zeros.");
            zero_nans_in_place(input);
        }
        const hicx::utilities::SparseEntries entries = input.entries();
        ic::Coo coo;
        coo.rows = input.rows;
        coo.cols = input.cols;
        coo.row = entries.row;
        coo.col = entries.col;
        coo.data = lc::doubles(input.data, "iterativeCorrection");
        if (ic::asymmetry(coo) > 1e-10) {
            raise_error("ValueError", "Please provide symmetric matrix!");
        }
        const double start_time = lc::now_seconds();
        log(kLogger, "info", "starting iterative correction");
        const std::int64_t M = to_int64(M_object);
        const double tolerance = to_double(tolerance_object);
        if (coo.rows == 0 && M > 0) {
            raise_error("ValueError", "zero-size array to reduction operation maximum which has no identity");
        }
        const ic::Result result = ic::correct(coo, M, tolerance, [&](std::int64_t iternum, double deviation) {
            if (verbose && iternum % 5 == 0) {
                log(kLogger, "info", format_estimate(iternum, M, lc::now_seconds() - start_time));
                log(kLogger, "info", "max delta - 1 = " + lc::float_str(deviation) + " ");
            }
            if (deviation < tolerance) {
                log(kLogger, "info",
                    "[iterative correction] " + std::to_string(iternum + 1) + " iterations used\n");
            }
        });
        if (result.outcome == ic::Outcome::LargeValuesInLoop) {
            log(kLogger, "error",
                "*Error* matrix correction is producing extremely large values. This is often caused by bins "
                "of low counts. Use a more stringent filtering of bins.");
            lc::python_exit(1);
        }
        if (result.outcome == ic::Outcome::LargeValuesAtEnd) {
            log(kLogger, "error",
                "*Error* matrix correction produced extremely large values. This is often caused by bins of "
                "low counts. Use a more stringent filtering of bins.");
            lc::python_exit(1);
        }
        std::vector<std::int32_t> row(coo.row.begin(), coo.row.end());
        std::vector<std::int32_t> col(coo.col.begin(), coo.col.end());
        const hicx::CsrMatrix corrected =
            hicx::CsrMatrix::from_coo(coo.rows, coo.cols, row, col, std::move(coo.data), "float64");
        return py::make_tuple(csr_matrix_object(corrected), from_vector(result.total_bias));
    }, py::arg("matrix"), py::arg("v") = py::none(), py::arg("M") = 50, py::arg("tolerance") = 1e-5,
       py::arg("verbose") = false);
}

}  // namespace hicx::python
