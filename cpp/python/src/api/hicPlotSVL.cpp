// hicexplorer._hicx.hicPlotSVL: compute_relation_short_long_range of
// hicexplorer/hicPlotSVL.py.
//
// The sums are hicx::svl (cpp/core/include/hicx/svl.hpp), which the tool
// hicPlotSVL runs too. With pIsCooler the chromosome is read from the cool
// file by the C++ matrix layer (hicx::ToolMatrix, as hiCMatrix(pMatrixFile,
// pChrnameList=[chromosome]) reads it); otherwise pHiCMatrix is a hiCMatrix
// object whose matrix and bins are read from its attributes, and its
// getBinSize cache is filled as the Python fills it. The result list
// [ratios, short range sums, long range sums] goes to pQueue.put and the
// function returns None, as in 3.7.6. The values are numpy scalars of the
// types numpy gives: int64 sums and float64 ratios for integer matrices,
// float32 for float32 and float64 for float64 matrices.
//
// Deviations: a cool file the C++ reader cannot open raises OSError with the
// reader's message (hicmatrix logs two warnings and h5py names the file "not
// an HDF5 file"); the matrix of a hiCMatrix object must be csr (3.7.6 sliced
// whatever scipy can slice).

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/svl.hpp"
#include "hicx/tool_matrix.hpp"

namespace hicx::python {

namespace {

constexpr const char* kFunction = "compute_relation_short_long_range";

// pDistance / bin size with Python's true division.
double max_distance(py::handle distance, py::handle bin_size) {
    if (PyLong_Check(bin_size.ptr()) && to_int64(bin_size) == 0) {
        raise_error("ZeroDivisionError", "division by zero");
    }
    return to_double(distance) / to_double(bin_size);
}

py::object value_of(const hicx::svl::BlockSums& sums, double value) {
    if (sums.kind == hicx::DType::Float32) {
        return numpy_scalar<float>(static_cast<float>(value));
    }
    return numpy_scalar<double>(value);
}

}  // namespace

void bind_hicPlotSVL(py::module_& m) {
    m.def(kFunction, [](py::handle pHiCMatrix, py::handle pChromosomes, py::handle pDistance, py::handle pIsCooler,
                        py::handle pQueue) {
        py::list svl_relations;
        py::list sum_smaller;
        py::list sum_greater;
        const bool is_cooler = hicobj::py_truthy(pIsCooler);
        std::optional<hicx::CsrMatrix> whole;
        for (const py::handle chromosome : pChromosomes) {
            hicx::svl::BlockSums sums;
            try {
                if (is_cooler) {
                    const std::string path = fs_path(pHiCMatrix);
                    const std::string name = utf8(py::str(chromosome));
                    std::optional<hicx::ToolMatrix> block;
                    try {
                        if (name.find(':') == std::string::npos) {
                            const std::vector<std::string> names = hicx::CoolFile(path).chrom_names();
                            if (std::find(names.begin(), names.end(), name) == names.end()) {
                                raise_error("ValueError", "Unknown sequence label: " + name);
                            }
                        }
                        block = hicx::ToolMatrix::load(path, name);
                    } catch (const hicx::PythonError&) {
                        throw;
                    } catch (const std::exception& error) {
                        raise_error("OSError", error.what());
                    }
                    const hicx::BinTable bins(block->cut_intervals());
                    const std::int64_t size = bins.bin_size();
                    if (!bins.bin_size_homogeneous()) {
                        char text[64];
                        std::snprintf(text, sizeof(text), "%f", static_cast<double>(size));
                        log(hicobj::kHicMatrixLogger, "warning",
                            std::string("Bin size is not homogeneous. ") + std::string(38, ' ') + "Median " + text +
                                "\n");
                    }
                    const py::object bin_size = bins.size() == 1 ? numpy_scalar<std::int64_t>(size)
                                                                 : py::object(py::int_(size));
                    sums = hicx::svl::block_sums(block->matrix(), 0, block->matrix().rows(),
                                                 max_distance(pDistance, bin_size));
                } else {
                    const auto [first, last] =
                        hicobj::range_values(hicobj::chr_bin_range(pHiCMatrix, chromosome));
                    if (!whole) {
                        whole = to_csr_matrix(sparse(pHiCMatrix.attr("matrix"), kFunction, "r"));
                    }
                    const double distance = max_distance(pDistance, hicobj::bin_size(pHiCMatrix));
                    const std::int64_t n = whole->rows();
                    const std::int64_t begin = std::clamp<std::int64_t>(first < 0 ? first + n : first, 0, n);
                    const std::int64_t end = std::clamp<std::int64_t>(last < 0 ? last + n : last, begin, n);
                    sums = hicx::svl::block_sums(*whole, begin, end, distance);
                }
            } catch (const hicx::svl::ExplicitZeros& zeros) {
                raise_error("IndexError",
                            "boolean index did not match indexed array along dimension 0; dimension is " +
                                std::to_string(zeros.stored()) + " but corresponding boolean dimension is " +
                                std::to_string(zeros.nonzero()));
            }
            if (!sums.kept) {
                continue;
            }
            if (sums.kind == hicx::DType::Integer) {
                svl_relations.append(numpy_scalar<double>(sums.ratio));
                sum_smaller.append(numpy_scalar<std::int64_t>(sums.smaller_int));
                sum_greater.append(numpy_scalar<std::int64_t>(sums.greater_int));
            } else {
                svl_relations.append(value_of(sums, sums.ratio));
                sum_smaller.append(value_of(sums, sums.smaller));
                sum_greater.append(value_of(sums, sums.greater));
            }
        }
        py::list result;
        result.append(svl_relations);
        result.append(sum_smaller);
        result.append(sum_greater);
        pQueue.attr("put")(result);
        return py::none();
    }, py::arg("pHiCMatrix"), py::arg("pChromosomes"), py::arg("pDistance"), py::arg("pIsCooler"),
       py::arg("pQueue"));
}

}  // namespace hicx::python
