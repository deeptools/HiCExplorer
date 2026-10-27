// hicexplorer._hicx.hicPCA: correlateEigenvectorWithGeneTrack and
// correlateEigenvectorWithHistonMarkTrack of hicexplorer/hicPCA.py.
//
// correlateEigenvectorWithGeneTrack reads the gene track through
// hicexplorer.utilities.opener and hicexplorer.readBed.ReadBed (both on C++),
// counts the genes per bin with the hiCMatrix queries of
// hicmatrix_queries.hpp and decides each chromosome's sign with
// hicx::pca::pearsonr (scipy 1.14's arithmetic, cpp/core/include/hicx/
// pca_tracks.hpp). The eigenvectors are an (n, k) array or nested list as
// main passes them; the result is a new array, the transpose of a C ordered
// (k, n) array, as in 3.7.6.
//
// correlateEigenvectorWithHistonMarkTrack queries the histone mark track it
// is given, a pyBigWig object, through its chroms() and stats() methods,
// exactly as often as 3.7.6 does, and flips the vectors of pEigenvector in
// place (float or complex arrays: complex vectors are compared
// lexicographically as numpy 1.26 does, and -1 * z is the complex product).
//
// Deviation: eigenvectors of other dtypes than float32, float64, complex64
// and complex128 raise TypeError naming HiCExplorer 4.

#include <algorithm>
#include <complex>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "api.hpp"
#include "hicmatrix_queries.hpp"
#include "hicx/pca_tracks.hpp"

namespace hicx::python {

namespace {

constexpr const char* kLogger = "hicexplorer.hicPCA";

py::array new_array(const std::vector<double>& values, py::ssize_t rows, py::ssize_t cols, DType dtype) {
    py::array_t<double> out({rows, cols});
    if (!values.empty()) {
        std::memcpy(out.mutable_data(), values.data(), values.size() * sizeof(double));
    }
    if (dtype == DType::Float32) {
        return out.attr("astype")("float32");
    }
    return out;
}

void size_warnings(py::handle size, py::handle interval) {
    log(kLogger, "warning", "Your chromosome sizes do not match the chromosome sizes of the extraTrack data!");
    log(kLogger, "warning",
        "Your chromosome " + str_of(interval.attr("chromosome")) + "; Size " + str_of(size) + ". ExtraTrack data " +
            str_of(interval.attr("chromosome")) + " " + str_of(interval.attr("start")) + " " +
            str_of(interval.attr("end")));
    log(kLogger, "warning",
        "Please create your interaction matrix with a chromosome size file! However, if the sizes are intended and "
        "it is accepted that certain regions are not part of the correlation, you can ignore this message.");
}

void warn_pearson(const hicx::pca::Pearson& result) {
    const py::object stats = py::module_::import("scipy.stats");
    if (result.constant) {
        if (PyErr_WarnEx(stats.attr("ConstantInputWarning").ptr(),
                         "An input array is constant; the correlation coefficient is not defined.", 2) != 0) {
            throw_python_error();
        }
    }
    if (result.nearly_constant) {
        if (PyErr_WarnEx(stats.attr("NearConstantInputWarning").ptr(),
                         "An input array is nearly constant; the computed correlation coefficient may be "
                         "inaccurate.",
                         2) != 0) {
            throw_python_error();
        }
    }
}

// vector > 0 and vector < 0 of a real or complex vector (numpy 1.26 orders
// complex numbers lexicographically).
template <class T>
bool positive(const T& value) {
    if constexpr (std::is_floating_point_v<T>) {
        return value > 0;
    } else {
        return value.real() > 0 || (value.real() == 0 && value.imag() > 0);
    }
}

template <class T>
bool negative(const T& value) {
    if constexpr (std::is_floating_point_v<T>) {
        return value < 0;
    } else {
        return value.real() < 0 || (value.real() == 0 && value.imag() < 0);
    }
}

// -1 * value: an exact negation for reals, numpy's complex product with
// (-1 + 0j) for complex numbers.
template <class T>
T times_minus_one(const T& value) {
    if constexpr (std::is_floating_point_v<T>) {
        return -value;
    } else {
        using R = typename T::value_type;
        const R a = value.real();
        const R b = value.imag();
        return T(a * R(-1) - b * R(0), a * R(0) + b * R(-1));
    }
}

}  // namespace

void bind_hicPCA(py::module_& m) {
    m.def("correlateEigenvectorWithGeneTrack", [](py::handle pMatrix, py::handle pEigenvector, py::handle pGeneTrack) {
        log(kLogger, "debug", "correlate eigenvector!");
        const py::object file_h = py::module_::import("hicexplorer.utilities").attr("opener")(pGeneTrack);
        const py::object bed = py::module_::import("hicexplorer.readBed").attr("ReadBed")(file_h);

        std::vector<double> gene_occurrence(hicobj::length(pMatrix.attr("cut_intervals")), 0.0);
        const py::list chromosome_list = hicobj::chr_names(pMatrix);

        for (const py::handle interval : bed) {
            const py::object chromosome_name = interval.attr("chromosome");
            if (!hicobj::py_contains(chromosome_list, chromosome_name)) {
                continue;
            }
            const py::object start = interval.attr("start");
            if (hicobj::py_less(hicobj::item(hicobj::chromosome_sizes(pMatrix), chromosome_name), start)) {
                size_warnings(hicobj::item(hicobj::chromosome_sizes(pMatrix), chromosome_name), interval);
                continue;
            }
            const auto bin_id = hicobj::region_bin_range_values(pMatrix, chromosome_name, start, interval.attr("end"));
            if (!bin_id) {
                size_warnings(hicobj::item(hicobj::chromosome_sizes(pMatrix), chromosome_name), interval);
                continue;
            }
            const std::int64_t n = static_cast<std::int64_t>(gene_occurrence.size());
            std::int64_t index = bin_id->second;
            if (index < 0) {
                index += n;
            }
            if (index < 0 || index >= n) {
                raise_error("IndexError", "index " + std::to_string(bin_id->second) + " is out of bounds for axis 0 with size " +
                                              std::to_string(n));
            }
            gene_occurrence[static_cast<std::size_t>(index)] += 1.0;
        }

        std::vector<std::pair<std::int64_t, std::int64_t>> ranges;
        for (const py::handle chromosome : chromosome_list) {
            ranges.push_back(hicobj::range_values(hicobj::chr_bin_range(pMatrix, chromosome)));
        }

        // np.array(pEigenvector).real.transpose()
        const py::array eigen = py::array::ensure(pEigenvector);
        if (!eigen) {
            throw_python_error();
        }
        const DType kind = dtype_of(eigen);
        if (kind != DType::Float64 && kind != DType::Float32 && kind != DType::Complex128 &&
            kind != DType::Complex64) {
            raise_unsupported_dtype(eigen, "correlateEigenvectorWithGeneTrack");
        }
        if (eigen.ndim() != 2) {
            if (eigen.ndim() == 1 && !chromosome_list.empty() && eigen.size() > 0) {
                raise_error("IndexError", "invalid index to scalar variable.");
            }
            raise_error("TypeError", "correlateEigenvectorWithGeneTrack: HiCExplorer 4 takes (bins, vectors) "
                                     "eigenvectors");
        }
        const auto rows = static_cast<std::size_t>(eigen.shape(0));
        const auto cols = static_cast<std::size_t>(eigen.shape(1));
        std::vector<double> real(rows * cols);
        if (kind == DType::Float64) {
            real = to_vector<double>(eigen);
        } else if (kind == DType::Float32) {
            const std::vector<float> values = to_vector<float>(eigen);
            real.assign(values.begin(), values.end());
        } else if (kind == DType::Complex128) {
            const std::vector<std::complex<double>> values = to_vector<std::complex<double>>(eigen);
            for (std::size_t k = 0; k < values.size(); ++k) {
                real[k] = values[k].real();
            }
        } else {
            const std::vector<std::complex<float>> values = to_vector<std::complex<float>>(eigen);
            for (std::size_t k = 0; k < values.size(); ++k) {
                real[k] = static_cast<double>(values[k].real());
            }
        }
        // transposed: vector i is column i
        std::vector<double> vectors(cols * rows);
        for (std::size_t r = 0; r < rows; ++r) {
            for (std::size_t c = 0; c < cols; ++c) {
                vectors[c * rows + r] = real[r * cols + c];
            }
        }
        for (const auto& [first, last] : ranges) {
            const auto clamp = [rows](std::int64_t index) {
                const auto n = static_cast<std::int64_t>(rows);
                if (index < 0) {
                    index += n;
                }
                return static_cast<std::size_t>(std::clamp<std::int64_t>(index, 0, n));
            };
            const std::size_t a = clamp(first);
            const std::size_t b = std::max(a, clamp(last));
            const std::size_t ga = std::min(static_cast<std::size_t>(std::max<std::int64_t>(first, 0)), gene_occurrence.size());
            const std::size_t gb = std::max(ga, std::min(static_cast<std::size_t>(std::max<std::int64_t>(last, 0)),
                                                         gene_occurrence.size()));
            const std::vector<double> genes(gene_occurrence.begin() + static_cast<std::ptrdiff_t>(ga),
                                            gene_occurrence.begin() + static_cast<std::ptrdiff_t>(gb));
            for (std::size_t i = 0; i < cols; ++i) {
                double* vector = vectors.data() + i * rows;
                std::vector<double> part(vector + a, vector + b);
                if (kind == DType::Float32 || kind == DType::Complex64) {
                    // the float32 values, widened by pearsonr's float64 dtype
                    for (double& value : part) {
                        value = static_cast<double>(static_cast<float>(value));
                    }
                }
                hicx::pca::Pearson correlation;
                try {
                    correlation = hicx::pca::pearsonr(part, genes);
                } catch (const std::invalid_argument& error) {
                    raise_error("ValueError", error.what());
                }
                warn_pearson(correlation);
                if (correlation.statistic < 0) {
                    for (std::size_t k = a; k < b; ++k) {
                        vector[k] = -vector[k];
                    }
                }
            }
        }
        // np.array(pEigenvector).transpose(): a C ordered (k, n) copy, transposed
        return new_array(vectors, static_cast<py::ssize_t>(cols), static_cast<py::ssize_t>(rows), (kind == DType::Float32 || kind == DType::Complex64) ? DType::Float32 : DType::Float64)
            .attr("transpose")();
    }, py::arg("pMatrix"), py::arg("pEigenvector"), py::arg("pGeneTrack"));

    m.def("correlateEigenvectorWithHistonMarkTrack", [](py::handle pEigenvector, py::handle bwTrack, py::handle chromosome,
                                                        py::handle start, py::handle end, py::handle pHistonMarkTrack,
                                                        py::handle pHistonMarkType) {
        (void)pHistonMarkTrack;
        const std::size_t count = hicobj::length(pEigenvector);
        for (std::size_t index = 0; index < count; ++index) {
            py::object vector_object = hicobj::item(pEigenvector, static_cast<std::int64_t>(index));
            py::array vector = as_ndarray(vector_object);
            const auto process = [&]<class T>() {
                std::vector<T> values = to_vector<T>(vector);
                std::vector<std::size_t> pos_indices;
                std::vector<std::size_t> neg_indices;
                for (std::size_t k = 0; k < values.size(); ++k) {
                    if (positive(values[k])) {
                        pos_indices.push_back(k);
                    }
                    if (negative(values[k])) {
                        neg_indices.push_back(k);
                    }
                }
                if (!hicobj::py_contains(bwTrack.attr("chroms")().attr("keys")(), chromosome)) {
                    return;
                }
                // the Python sums a Python float from the int 0
                const auto mean_of = [&](const std::vector<std::size_t>& indices) {
                    double sum = 0.0;
                    bool any = false;
                    for (const std::size_t ind : indices) {
                        const py::int_ position(ind);
                        const py::object first = hicobj::item(
                            bwTrack.attr("stats")(chromosome, hicobj::item(start, position), hicobj::item(end, position)), 0);
                        if (hicobj::py_truthy(first)) {
                            const py::object again = hicobj::item(
                                bwTrack.attr("stats")(chromosome, hicobj::item(start, position), hicobj::item(end, position)),
                                0);
                            sum += to_double(again);
                            any = true;
                        }
                    }
                    (void)any;
                    return sum != 0.0 ? sum / static_cast<double>(indices.size()) : 0.0;
                };
                const double pos_mean = mean_of(pos_indices);
                const double neg_mean = mean_of(neg_indices);
                bool active = hicobj::py_equal(pHistonMarkType, py::str("active"));
                if (!active) {
                    if (py::module_::import("sys").attr("flags").attr("optimize").cast<int>() == 0 &&
                        !hicobj::py_equal(pHistonMarkType, py::str("inactive"))) {
                        raise_error("AssertionError", "");
                    }
                }
                if (hicx::pca::histone_flip(pos_mean, neg_mean, active)) {
                    for (const std::size_t k : pos_indices) {
                        values[k] = active ? -values[k] : times_minus_one(values[k]);
                    }
                    for (const std::size_t k : neg_indices) {
                        values[k] = active ? -values[k] : times_minus_one(values[k]);
                    }
                    write_back(vector, values);
                }
            };
            switch (dtype_of(vector)) {
                case DType::Float64: process.template operator()<double>(); break;
                case DType::Float32: process.template operator()<float>(); break;
                case DType::Complex128: process.template operator()<std::complex<double>>(); break;
                case DType::Complex64: process.template operator()<std::complex<float>>(); break;
                default: raise_unsupported_dtype(vector, "correlateEigenvectorWithHistonMarkTrack");
            }
            // pEigenvector[index] = vector
            if (PyObject_SetItem(pEigenvector.ptr(), py::int_(index).ptr(), vector_object.ptr()) != 0) {
                throw_python_error();
            }
        }
        return py::none();
    }, py::arg("pEigenvector"), py::arg("bwTrack"), py::arg("chromosome"), py::arg("start"), py::arg("end"),
       py::arg("pHistonMarkTrack"), py::arg("pHistonMarkType"));
}

}  // namespace hicx::python
