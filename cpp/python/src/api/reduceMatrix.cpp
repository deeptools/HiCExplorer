// hicexplorer._hicx.reduceMatrix: the C++ body of
// hicexplorer.reduceMatrix.reduce_matrix.
//
// The reduction is hicx::reduce_entries (cpp/core/include/hicx/
// reduce_matrix_entries.hpp) over the stored entries in the order the
// Python's triu(..., format='coo') or tocoo() lists them; this file converts
// the arguments, keeps the 3.7.6 early return of the argument itself, the
// NaN warning and the two timing debug records, and builds the result: a
// csr_matrix with use_triu, a coo_matrix without (also for a sparse array
// argument, as the Python builds a coo_matrix).
//
// Deviations: data dtypes other than int32, int64, float32 and float64, and
// bool indices in bins_to_merge, raise TypeError naming HiCExplorer 4.

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "api.hpp"
#include "hicx/reduce_matrix_entries.hpp"
#include "loops_correct_common.hpp"

namespace hicx::python {

namespace re = hicx::reduce_entries;

namespace {

constexpr const char* kLogger = "hicexplorer.reduceMatrix";

// numpy's float64 -> T cast (astype), a C conversion.
template <class T>
double cast_value(double value) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<double>(static_cast<T>(value));
    } else {
        // x86 cvttsd2si: NaN and out of range values become the minimum
        const double low = static_cast<double>(std::numeric_limits<T>::min());
        const double high = -low;
        if (std::isnan(value) || value < low || value >= high) {
            return static_cast<double>(std::numeric_limits<T>::min());
        }
        return static_cast<double>(static_cast<T>(value));
    }
}

template <class T>
T wrap(double value) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(value);
    } else {
        return static_cast<T>(static_cast<std::int64_t>(value));
    }
}

template <class T>
re::Arithmetic arithmetic() {
    re::Arithmetic a;
    a.cast = [](double v) { return cast_value<T>(v); };
    a.add = [](double x, double y) {
        if constexpr (std::is_floating_point_v<T>) {
            return static_cast<double>(static_cast<T>(static_cast<T>(x) + static_cast<T>(y)));
        } else {
            using U = std::make_unsigned_t<T>;
            return static_cast<double>(static_cast<T>(static_cast<U>(wrap<T>(x)) + static_cast<U>(wrap<T>(y))));
        }
    };
    a.subtract = [](double x, double y) {
        if constexpr (std::is_floating_point_v<T>) {
            return static_cast<double>(static_cast<T>(static_cast<T>(x) - static_cast<T>(y)));
        } else {
            using U = std::make_unsigned_t<T>;
            return static_cast<double>(static_cast<T>(static_cast<U>(wrap<T>(x)) - static_cast<U>(wrap<T>(y))));
        }
    };
    a.twice = [](double x) {
        if constexpr (std::is_floating_point_v<T>) {
            return static_cast<double>(static_cast<T>(static_cast<T>(x) * static_cast<T>(2)));
        } else {
            using U = std::make_unsigned_t<T>;
            return static_cast<double>(static_cast<T>(static_cast<U>(wrap<T>(x)) * static_cast<U>(2)));
        }
    };
    return a;
}

std::string seconds(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.5f", value);
    return buffer;
}

template <class T>
py::array typed(const std::vector<double>& values) {
    std::vector<T> out;
    out.reserve(values.size());
    for (const double v : values) {
        out.push_back(wrap<T>(v));
    }
    return from_vector(out);
}

py::object reduce_matrix(py::handle matrix, py::handle bins_to_merge, py::handle use_triu_object,
                         py::handle diagonal_object) {
    const bool use_triu = lc::truthy(use_triu_object);
    Sparse input = sparse(matrix, "reduce_matrix");
    const hicx::utilities::SparseEntries stored = input.entries();
    re::Entries entries;
    entries.rows = input.rows;
    entries.cols = input.cols;
    entries.row = stored.row;
    entries.col = stored.col;
    dispatch_arithmetic(input.data, "reduce_matrix", [&]<class T>() {
        for (const T v : to_vector<T>(input.data)) {
            entries.data.push_back(static_cast<double>(v));
        }
    });
    if (use_triu) {
        entries = re::upper_triangle(entries);
    }
    const std::int64_t M = static_cast<std::int64_t>(lc::length(bins_to_merge));
    double start_time = lc::now_seconds();
    if (M == entries.rows) {
        return py::reinterpret_borrow<py::object>(matrix);
    }
    const std::size_t num_nan = re::count_nan(entries);
    if (num_nan > 0) {
        log(kLogger, "warning", "*Warning*\nmatrix contains " + std::to_string(num_nan) + " NaN values.");
    }
    std::vector<std::int64_t> group_of(static_cast<std::size_t>(entries.rows), -1);
    std::int64_t k = 0;
    for (const py::handle group : bins_to_merge) {
        for (const py::handle x : group) {
            if (PyBool_Check(x.ptr())) {
                raise_error("TypeError", "reduce_matrix: HiCExplorer 4 takes integer bin indices, not bool");
            }
            if (!is_integer(x)) {
                raise_error("IndexError",
                            "only integers, slices (`:`), ellipsis (`...`), numpy.newaxis (`None`) and integer or "
                            "boolean arrays are valid indices");
            }
            std::int64_t index = to_int64(x);
            if (index < -entries.rows || index >= entries.rows) {
                raise_error("IndexError", "index " + std::to_string(index) + " is out of bounds for axis 0 with size " +
                                              std::to_string(entries.rows));
            }
            if (index < 0) {
                index += entries.rows;
            }
            group_of[static_cast<std::size_t>(index)] = k;
        }
        ++k;
    }
    log(kLogger, "debug", "time complex unique: " + seconds(lc::now_seconds() - start_time));
    start_time = lc::now_seconds();
    const bool diagonal = diagonal_object.ptr() != Py_False;
    re::Cells cells;
    py::array data;
    dispatch_arithmetic(input.data, "reduce_matrix", [&]<class T>() {
        cells = re::reduce(entries, group_of, M, use_triu, diagonal, arithmetic<T>());
        data = typed<T>(cells.data);
    });
    const py::module_ sp = py::module_::import("scipy.sparse");
    const py::tuple shape = py::make_tuple(M, M);
    if (!use_triu) {
        std::vector<std::int32_t> row(cells.row.begin(), cells.row.end());
        std::vector<std::int32_t> col(cells.col.begin(), cells.col.end());
        py::object result = sp.attr("coo_matrix")(
            py::make_tuple(data, py::make_tuple(from_vector(row), from_vector(col))), py::arg("shape") = shape);
        log(kLogger, "debug", "time create sparse matrix: " + seconds(lc::now_seconds() - start_time));
        return result;
    }
    log(kLogger, "debug", "time create sparse matrix: " + seconds(lc::now_seconds() - start_time));
    std::vector<std::int32_t> indptr(static_cast<std::size_t>(M) + 1, 0);
    std::vector<std::int32_t> indices(cells.col.begin(), cells.col.end());
    for (const std::int64_t r : cells.row) {
        ++indptr[static_cast<std::size_t>(r) + 1];
    }
    for (std::size_t i = 1; i < indptr.size(); ++i) {
        indptr[i] += indptr[i - 1];
    }
    return sp.attr("csr_matrix")(
        py::make_tuple(data, from_vector(indices), from_vector(indptr)), py::arg("shape") = shape);
}

}  // namespace

void bind_reduceMatrix(py::module_& m) {
    m.def("reduce_matrix", &reduce_matrix, py::arg("matrix"), py::arg("bins_to_merge"),
          py::arg("use_triu") = true, py::arg("diagonal") = false);
}

}  // namespace hicx::python
