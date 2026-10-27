// hicexplorer._hicx.hicCompartmentalization: get_indices,
// count_interactions and within_vs_between_compartments of
// hicexplorer/hicCompartmentalization.py.
//
// The arithmetic is hicx::compartments (cpp/core/include/hicx/
// compartmentalization.hpp), which the tool hicCompartmentalization runs too.
// pc1 is the pandas DataFrame of the Python, read through its columns
// ("chr", "quantile", "bin_id"); obs_exp a hiCMatrix whose matrix is a scipy
// csr matrix.
//
// count_interactions with an offset writes NaN into the named diagonals of
// obs_exp.matrix in place, as the scipy assignment of 3.7.6 does: existing
// entries are overwritten in the data array, missing ones are inserted
// (scipy's _insert_many, with its SparseEfficiencyWarning), the index arrays
// keep int32 unless they were int64, and a matrix with sorted indices is
// sorted again. An integer matrix raises numpy's ValueError for a NaN.
//
// Deviations: count_interactions reads float64 and integer matrices (a
// float32 matrix raises TypeError naming HiCExplorer 4, 3.7.6 summed its
// blocks in float32), and the matrix must hold no duplicate entries.
// within_vs_between_compartments takes a two dimensional float64 array.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/compartmentalization.hpp"

namespace hicx::python {

namespace {

namespace cm = hicx::compartments;

bool assertions_enabled() {
    return py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0;
}

// row['end'] - 1
py::object minus_one(py::handle value) {
    if (PyFloat_Check(value.ptr())) {
        return py::float_(to_double(value) - 1.0);
    }
    return py::int_(hicobj::int64_of(value) - 1);
}

py::object region(py::handle obs_exp, py::handle row) {
    py::object range = hicobj::region_bin_range(obs_exp, hicobj::item(row, py::str("chr")),
                                                hicobj::item(row, py::str("start")),
                                                minus_one(hicobj::item(row, py::str("end"))));
    if (range.is_none()) {
        raise_error("TypeError", "'NoneType' object is not subscriptable");
    }
    return range;
}

std::vector<std::int64_t> index_values(const py::array& array) {
    std::vector<std::int64_t> out;
    if (dtype_of(array) == DType::Int32) {
        for (const std::int32_t v : to_vector<std::int32_t>(array)) {
            out.push_back(v);
        }
    } else if (dtype_of(array) == DType::Int64) {
        out = to_vector<std::int64_t>(array);
    } else {
        raise_error("TypeError", "count_interactions: HiCExplorer 4 takes int32 or int64 sparse indices");
    }
    return out;
}

template <class T>
py::array index_array(const std::vector<std::int64_t>& values) {
    std::vector<T> narrow(values.begin(), values.end());
    return from_vector(narrow);
}

// obs_exp.matrix[rows, cols] = np.nan for the positions (rows[k], cols[k]),
// as scipy's _set_many does it on a csr matrix.
void assign_nan(py::handle matrix_object, const std::vector<std::int64_t>& rows,
                const std::vector<std::int64_t>& cols) {
    if (rows.empty()) {
        return;
    }
    Sparse matrix = sparse(matrix_object, "count_interactions", "r");
    const DType kind = dtype_of(matrix.data);
    if (kind != DType::Float64 && kind != DType::Float32) {
        raise_error("ValueError", "cannot convert float NaN to integer");
    }
    std::vector<std::int64_t> indptr = index_values(matrix.indptr);
    std::vector<std::int64_t> indices = index_values(matrix.indices);
    const bool wide = dtype_of(matrix.indptr) == DType::Int64 || dtype_of(matrix.indices) == DType::Int64;
    // positions of the existing entries; -1 for the ones to insert
    std::vector<std::int64_t> offsets(rows.size(), -1);
    for (std::size_t k = 0; k < rows.size(); ++k) {
        for (std::int64_t p = indptr[static_cast<std::size_t>(rows[k])];
             p < indptr[static_cast<std::size_t>(rows[k]) + 1]; ++p) {
            if (indices[static_cast<std::size_t>(p)] == cols[k]) {
                offsets[k] = p;
                break;
            }
        }
    }
    // self.data[offsets[mask]] = x[mask], in place
    py::array data = matrix.data;
    const auto set_existing = [&]<class T>() {
        std::vector<T> values = to_vector<T>(data);
        for (const std::int64_t p : offsets) {
            if (p >= 0) {
                values[static_cast<std::size_t>(p)] = std::numeric_limits<T>::quiet_NaN();
            }
        }
        write_back(data, values);
    };
    if (kind == DType::Float64) {
        set_existing.template operator()<double>();
    } else {
        set_existing.template operator()<float>();
    }
    std::vector<std::size_t> insert;
    for (std::size_t k = 0; k < rows.size(); ++k) {
        if (offsets[k] < 0) {
            insert.push_back(k);
        }
    }
    if (insert.empty()) {
        return;
    }
    warn(py::module_::import("scipy.sparse").attr("SparseEfficiencyWarning").ptr(),
         "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.", 3);
    const bool sorted = hicobj::py_truthy(matrix_object.attr("has_sorted_indices"));
    const std::size_t n_rows = indptr.size() - 1;
    // the new columns per row, in the stable order of the insertion
    std::vector<std::vector<std::int64_t>> added(n_rows);
    for (const std::size_t k : insert) {
        added[static_cast<std::size_t>(rows[k])].push_back(cols[k]);
    }
    const auto rebuild = [&]<class T>() {
        const std::vector<T> old = to_vector<T>(data);
        std::vector<T> values;
        std::vector<std::int64_t> new_indices;
        std::vector<std::int64_t> new_indptr(n_rows + 1, 0);
        values.reserve(old.size() + insert.size());
        new_indices.reserve(old.size() + insert.size());
        for (std::size_t r = 0; r < n_rows; ++r) {
            std::vector<std::pair<std::int64_t, T>> entries;
            for (const std::int64_t c : added[r]) {
                entries.emplace_back(c, std::numeric_limits<T>::quiet_NaN());
            }
            for (std::int64_t p = indptr[r]; p < indptr[r + 1]; ++p) {
                entries.emplace_back(indices[static_cast<std::size_t>(p)], old[static_cast<std::size_t>(p)]);
            }
            if (sorted) {
                std::stable_sort(entries.begin(), entries.end(),
                                 [](const auto& a, const auto& b) { return a.first < b.first; });
            }
            for (const auto& [c, v] : entries) {
                new_indices.push_back(c);
                values.push_back(v);
            }
            new_indptr[r + 1] = static_cast<std::int64_t>(new_indices.size());
        }
        matrix_object.attr("indices") = wide ? index_array<std::int64_t>(new_indices)
                                              : index_array<std::int32_t>(new_indices);
        matrix_object.attr("data") = from_vector(values);
        matrix_object.attr("indptr") = wide ? index_array<std::int64_t>(new_indptr)
                                             : index_array<std::int32_t>(new_indptr);
    };
    if (kind == DType::Float64) {
        rebuild.template operator()<double>();
    } else {
        rebuild.template operator()<float>();
    }
    if (sorted) {
        matrix_object.attr("has_sorted_indices") = py::bool_(true);
    }
}

py::array square_array(const std::vector<double>& values, std::int64_t size) {
    py::array_t<double> out({static_cast<py::ssize_t>(size), static_cast<py::ssize_t>(size)});
    if (!values.empty()) {
        std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(double));
    }
    return out;
}

}  // namespace

void bind_hicCompartmentalization(py::module_& m) {
    m.def("get_indices", [](py::handle obs_exp, py::handle row) {
        const std::int64_t first = to_int64(hicobj::item(region(obs_exp, row), 0));
        const std::int64_t last = to_int64(hicobj::item(region(obs_exp, row), 1)) + 1;
        std::vector<std::int64_t> indices;
        for (std::int64_t i = first; i < last; ++i) {
            indices.push_back(i);
        }
        return from_vector(indices);
    }, py::arg("obs_exp"), py::arg("row"));

    m.def("count_interactions", [](py::handle obs_exp, py::handle pc1, py::handle quantiles_number,
                                   py::handle offset) {
        std::set<std::string> chromosomes;
        for (const py::handle name : hicobj::item(pc1, py::str("chr"))) {
            chromosomes.insert(str_of(name));
        }
        const std::int64_t quantiles = to_int64(quantiles_number);
        if (quantiles < 0) {
            raise_error("ValueError", "negative dimensions are not allowed");
        }
        if (hicobj::py_truthy(offset)) {
            const py::object matrix = obs_exp.attr("matrix");
            const std::int64_t n = hicobj::int64_of(hicobj::item(matrix.attr("shape"), 0));
            for (const py::handle dist_object : offset) {
                const std::int64_t dist = hicobj::int64_of(dist_object);
                if (assertions_enabled() && dist < 0) {
                    raise_error("AssertionError", "");
                }
                std::vector<std::int64_t> rows;
                std::vector<std::int64_t> cols;
                for (std::int64_t i = 0; i < n - dist; ++i) {
                    rows.push_back(i);
                    cols.push_back(i + dist);
                }
                assign_nan(obs_exp.attr("matrix"), rows, cols);
                assign_nan(obs_exp.attr("matrix"), cols, rows);
            }
        }
        std::vector<std::vector<std::int64_t>> bin_ids;
        std::vector<std::int64_t> quantile_of_row;
        if (!chromosomes.empty() && quantiles > 0) {
            for (const py::handle q : hicobj::item(pc1, py::str("quantile"))) {
                quantile_of_row.push_back(hicobj::int64_of(q));
            }
            for (const py::handle ids : hicobj::item(pc1, py::str("bin_id"))) {
                std::vector<std::int64_t> row;
                for (const py::handle id : ids) {
                    row.push_back(to_int64(id));
                }
                bin_ids.push_back(std::move(row));
            }
        }
        const Sparse matrix = sparse(obs_exp.attr("matrix"), "count_interactions", "r");
        const DType kind = dtype_of(matrix.data);
        if (kind == DType::Float32) {
            raise_error("TypeError", "count_interactions: HiCExplorer 4 takes float64 or integer matrices");
        }
        const hicx::CsrMatrix csr = to_csr_matrix(matrix);
        for (const auto& row : bin_ids) {
            for (const std::int64_t id : row) {
                if (id < 0 || id >= csr.rows() || id >= csr.cols()) {
                    raise_error("IndexError", "index (" + std::to_string(id) + ") out of range");
                }
            }
        }
        const std::vector<double> result = cm::count_interactions(
            csr, bin_ids, quantile_of_row, quantiles, {}, static_cast<std::int64_t>(chromosomes.size()));
        return square_array(result, quantiles);
    }, py::arg("obs_exp"), py::arg("pc1"), py::arg("quantiles_number"), py::arg("offset"));

    m.def("within_vs_between_compartments", [](py::handle normalised_sum_per_quantile, py::handle quantiles_number) {
        const py::array array = as_ndarray(normalised_sum_per_quantile);
        if (array.ndim() != 2) {
            raise_error("IndexError", "too many indices for array: array is " + std::to_string(array.ndim()) +
                                          "-dimensional, but 2 were indexed");
        }
        if (dtype_of(array) != DType::Float64) {
            raise_unsupported_dtype(array, "within_vs_between_compartments");
        }
        const std::vector<double> ratios =
            cm::within_vs_between(to_vector<double>(array), array.shape(0), array.shape(1),
                                  to_int64(quantiles_number));
        py::list out;
        for (const double ratio : ratios) {
            out.append(numpy_scalar<double>(ratio));
        }
        return out;
    }, py::arg("normalised_sum_per_quantile"), py::arg("quantiles_number"));
}

}  // namespace hicx::python
