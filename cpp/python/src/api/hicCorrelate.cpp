// hicexplorer._hicx.hicCorrelate: the C++ bodies of the public functions of
// hicexplorer.hicCorrelate.
//
// get_vectors is hicx::correlate::get_vectors (scipy.sparse's arithmetic on
// the stored entries). plot_correlation draws with matplotlib and scipy's
// dendrogram; the numbers it draws from, the complete linkage of the
// correlation matrix and the colour range, are plot_correlation_data here.
//
// Differences from 3.7.6 (cpp/python/API_BINDINGS.md):
//  * get_vectors takes csr, csc and coo matrices of dtype int32, int64,
//    float32 or float64; any other sparse format or dtype raises TypeError
//    naming HiCExplorer 4, and a dense array raises TypeError where 3.7.6
//    failed later on `_mat.data += 1`.
//  * plot_correlation takes a two dimensional ndarray; a one dimensional one
//    (scipy reads it as a condensed distance matrix) raises TypeError.

#include <algorithm>
#include <cmath>
#include <string>

#include "api.hpp"
#include "hicx/correlate_api.hpp"

namespace hicx::python {

namespace corr = hicx::correlate;

namespace {

bool equal(py::handle a, py::handle b) {
    const int result = PyObject_RichCompareBool(a.ptr(), b.ptr(), Py_EQ);
    if (result < 0) {
        throw_python_error();
    }
    return result == 1;
}

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

// scipy.sparse's upcast of the four supported dtypes.
DType upcast(DType a, DType b) {
    const auto is_float = [](DType d) { return d == DType::Float32 || d == DType::Float64; };
    if (a == DType::Float64 || b == DType::Float64) {
        return DType::Float64;
    }
    if (a == DType::Float32 && b == DType::Float32) {
        return DType::Float32;
    }
    if (is_float(a) || is_float(b)) {
        return DType::Float64;
    }
    return (a == DType::Int64 || b == DType::Int64) ? DType::Int64 : DType::Int32;
}

DType arithmetic_dtype(const py::array& data) {
    return dispatch_arithmetic(data, "get_vectors", [&]<class T>() { return dtype_of(data); });
}

template <class T>
std::vector<T> converted_data(const py::array& data) {
    return dispatch_arithmetic(data, "get_vectors", [&]<class S>() {
        const std::vector<S> source = to_vector<S>(data);
        return std::vector<T>(source.begin(), source.end());
    });
}

template <class T>
corr::SparseInput<T> sparse_input(const Sparse& matrix) {
    corr::SparseInput<T> out;
    out.rows = matrix.rows;
    out.cols = matrix.cols;
    out.data = converted_data<T>(matrix.data);
    const hicx::utilities::SparseEntries entries = matrix.entries();
    switch (matrix.format) {
        case 'r':
            out.format = corr::SparseFormat::Csr;
            out.indptr = matrix.indptr_vector();
            out.minor = entries.col;
            break;
        case 'c':
            out.format = corr::SparseFormat::Csc;
            out.indptr = matrix.indptr_vector();
            out.minor = entries.row;
            break;
        default:
            out.format = corr::SparseFormat::Coo;
            out.major = entries.row;
            out.minor = entries.col;
            break;
    }
    return out;
}

template <class T>
py::object vectors_of(const Sparse& first, const Sparse& second) {
    const auto result = corr::get_vectors(sparse_input<T>(first), sparse_input<T>(second));
    return py::make_tuple(from_vector(result.first), from_vector(result.second));
}

py::object get_vectors(py::handle mat1, py::handle mat2) {
    const py::object shape1 = mat1.attr("shape");
    const py::object shape2 = mat2.attr("shape");
    if (!equal(shape1, shape2)) {
        if (assertions_enabled()) {
            raise_error("AssertionError",
                        "Matrices have different shapes. Computation of correlation is not "
                        "possible.");
        }
        raise_error("ValueError", "inconsistent shapes");
    }
    const Sparse first = sparse(mat1, "get_vectors");
    const Sparse second = sparse(mat2, "get_vectors");
    const DType type = upcast(arithmetic_dtype(first.data), arithmetic_dtype(second.data));
    switch (type) {
        case DType::Int32: return vectors_of<std::int32_t>(first, second);
        case DType::Int64: return vectors_of<std::int64_t>(first, second);
        case DType::Float32: return vectors_of<float>(first, second);
        default: return vectors_of<double>(first, second);
    }
}

py::object plot_correlation_data(py::handle corr_matrix, py::handle vmax, py::handle vmin) {
    if (!is_ndarray(corr_matrix)) {
        corr_matrix.attr("shape");
        raise_error("TypeError", "plot_correlation: expected a numpy.ndarray, got " +
                                     type_name(corr_matrix));
    }
    const py::array array = py::reinterpret_borrow<py::array>(corr_matrix);
    if (array.ndim() == 0) {
        raise_error("IndexError", "tuple index out of range");
    }
    const std::vector<double> values = dispatch_numeric(array, "plot_correlation", [&]<class T>() {
        if constexpr (std::is_same_v<T, std::complex<float>> ||
                      std::is_same_v<T, std::complex<double>> || std::is_same_v<T, long double>) {
            raise_error("TypeError", "plot_correlation: HiCExplorer 4 takes a real matrix");
            return std::vector<double>();
        } else {
            const std::vector<T> source = to_vector<T>(array);
            return std::vector<double>(source.begin(), source.end());
        }
    });
    py::object max_value = vmax.is_none() ? py::object(py::int_(1)) : py::reinterpret_borrow<py::object>(vmax);
    py::object min_value = py::reinterpret_borrow<py::object>(vmin);
    if (vmin.is_none()) {
        min_value = py::int_(corr::matrix_minimum(values) >= 0 ? 0 : -1);
    }
    if (array.ndim() == 1) {
        raise_error("TypeError",
                    "plot_correlation: HiCExplorer 4 takes a two dimensional matrix");
    }
    if (array.ndim() > 2) {
        raise_error("ValueError", "`y` must be 1 or 2 dimensional.");
    }
    const std::int64_t rows = array.shape(0);
    const std::int64_t cols = array.shape(1);
    if (corr::looks_like_distance_matrix(values, rows, cols)) {
        warn(py::module_::import("scipy.cluster.hierarchy").attr("ClusterWarning").ptr(),
             "The symmetric non-negative hollow observation matrix looks suspiciously like an "
             "uncondensed distance matrix",
             1);
    }
    const std::vector<double> linkage = corr::complete_linkage_of_rows(values, rows, cols);
    py::array out = from_vector(linkage);
    out = out.reshape(std::vector<py::ssize_t>{rows - 1, 4});
    return py::make_tuple(out, max_value, min_value);
}

}  // namespace

void bind_hicCorrelate(py::module_& m) {
    m.def("get_vectors", &get_vectors, py::arg("mat1"), py::arg("mat2"));
    m.def("plot_correlation_data", &plot_correlation_data, py::arg("corr_matrix"),
          py::arg("vmax"), py::arg("vmin"));
}

}  // namespace hicx::python
