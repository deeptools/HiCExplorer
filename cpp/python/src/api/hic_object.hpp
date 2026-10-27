// hicmatrix.HiCMatrix.hiCMatrix(path, pChrnameList) for the bindings of
// hicexplorer.hicAdjustMatrix and hicexplorer.hicDetectLoops: the matrix is
// read by the C++ matrix layer (hicx::read_cool, hicx::read_hicexplorer_h5,
// with the weights applied there) and handed over as a hiCMatrix whose
// attributes are what the 3.7.6 constructor leaves: the symmetric matrix
// after fillLowerTriangle, cut_intervals with the field types of the format
// (numpy.str_, numpy.int64, numpy.int64, numpy.float64 for h5; str, int, int,
// float for cool), nan_bins (float64 when empty, as np.array([]) is),
// correction_factors and distance_counts with the load time swap, the
// interval trees and chrBinBoundaries (built by the object's own
// intervalListToIntervalTree), and a MatrixFileHandler of the input's format
// carrying the Cool state a load leaves (correction operator, versions,
// metadata), so that save() writes what it wrote after a 3.7.6 load.
//
// Deviation: the cool metadata is the info dictionary with every value as
// its str(); hiCMatrix.save reads only three string entries of it.

#ifndef HICX_PYTHON_API_HIC_OBJECT_HPP
#define HICX_PYTHON_API_HIC_OBJECT_HPP

#include <optional>
#include <string>

#include "api.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/matrix_data.hpp"

namespace hicx::python::hic {

inline bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

struct Loaded {
    hicx::MatrixData data;  // hiCMatrix meaning of correction_factors / distance_counts
    bool h5 = false;
    hicx::CoolLoadResult cool_state;
};

// What the hiCMatrix constructor computes, in C++.
inline Loaded load(const std::string& path, const std::optional<std::string>& chromosome,
                   const std::optional<std::int64_t>& distance) {
    Loaded out;
    if (ends_with(path, ".h5")) {
        out.data = hicx::read_hicexplorer_h5(path);
        out.h5 = true;
    } else {
        hicx::CoolLoadOptions options;
        options.chrom_name = chromosome;
        options.distance = distance;
        out.cool_state = hicx::read_cool(path, options);
        out.data = std::move(out.cool_state.data);
    }
    std::swap(out.data.correction_factors, out.data.distance_counts);
    // fillLowerTriangle: when nothing is stored below the diagonal, M + triu(M, 1).T
    if (out.data.matrix.lower_triangle_is_zero()) {
        out.data.matrix.materialize_full();
        out.data.matrix.eliminate_zeros();
    }
    return out;
}

inline py::object optional_array(const std::optional<std::vector<double>>& values, bool column) {
    if (!values) {
        return py::none();
    }
    py::array array = from_vector(*values);
    if (column) {
        return array.attr("reshape")(py::make_tuple(static_cast<py::ssize_t>(values->size()), 1));
    }
    return array;
}

// The hiCMatrix object.
inline py::object object(const Loaded& loaded, py::handle path, py::handle chrname_list, bool interval_tree) {
    const py::module_ hm = py::module_::import("hicmatrix.HiCMatrix");
    py::object hic = hm.attr("hiCMatrix")();
    const py::object handler_class = py::module_::import("hicmatrix.lib").attr("MatrixFileHandler");
    py::object handler = handler_class(py::arg("pFileType") = loaded.h5 ? "h5" : "cool", py::arg("pMatrixFile") = path,
                                       py::arg("pChrnameList") = chrname_list);
    if (!loaded.h5) {
        py::object file = handler.attr("matrixFile");
        file.attr("minValue") = py::none();
        file.attr("maxValue") = py::none();
        py::dict metadata;
        for (const auto& [key, value] : loaded.cool_state.metadata) {
            metadata[to_str(key)] = to_str(value);
        }
        file.attr("hic_metadata") = metadata;
        if (loaded.cool_state.correction_operator) {
            file.attr("correctionOperator") = py::str(std::string(1, *loaded.cool_state.correction_operator));
        }
        if (loaded.cool_state.hic2cool_version) {
            file.attr("hic2cool_version") = to_str(*loaded.cool_state.hic2cool_version);
        }
        if (loaded.cool_state.hicmatrix_version) {
            file.attr("hicmatrix_version") = to_str(*loaded.cool_state.hicmatrix_version);
        }
    }
    hic.attr("matrixFileHandler") = handler;
    hic.attr("matrix") = csr_matrix_object(loaded.data.matrix);
    py::object like = py::none();
    if (loaded.h5) {
        const py::module_ np = py::module_::import("numpy");
        like = py::make_tuple(np.attr("str_")(""), np.attr("int64")(0), np.attr("int64")(0), np.attr("float64")(0.0));
    }
    hic.attr("cut_intervals") = cut_intervals_object(loaded.data.cut_intervals, like);
    if (loaded.data.nan_bins.empty()) {
        hic.attr("nan_bins") = py::module_::import("numpy").attr("array")(py::list());
    } else {
        hic.attr("nan_bins") = from_vector(loaded.data.nan_bins);
    }
    hic.attr("correction_factors") =
        optional_array(loaded.data.correction_factors, loaded.data.correction_factors_are_column);
    hic.attr("distance_counts") = optional_array(loaded.data.distance_counts, false);
    if (interval_tree) {
        const py::tuple trees = hic.attr("intervalListToIntervalTree")(hic.attr("cut_intervals"));
        hic.attr("interval_trees") = trees[0];
        hic.attr("chrBinBoundaries") = trees[1];
    }
    return hic;
}

}  // namespace hicx::python::hic

#endif  // HICX_PYTHON_API_HIC_OBJECT_HPP
