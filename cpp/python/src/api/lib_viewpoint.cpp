// hicexplorer._hicx.lib_viewpoint: the C++ bodies of the methods of
// hicexplorer.lib.viewpoint.Viewpoint.
//
// The Python class keeps its attribute (hicMatrix) and its methods; each
// method calls the function here with `self`. Calls from one method to
// another go through `self`, as in the Python, so a subclass that overrides a
// method still sees it called.
//
// What is C++ and what is the object protocol. The computations are C++: the
// row sums and the collapsing of the reference point's bins
// (hicx::chic::collapse_viewpoint), the smoothing, the relative values, the
// negative binomial p-values (hicx::scipy::betainc, which is what
// hicexplorer.lib.cnb.cdf evaluates), the merging of neighbours, the rbz-score
// and the parsing and formatting of the text files; the HDF5 files the
// methods open by path are read with hicx::chic::H5Reader. The arguments are
// Python objects, and what they are asked for goes through their own
// interface as in 3.7.6: the bin lookups of the hiCMatrix object
// (getRegionBinRange, getBinPos, getChrBinRange, getBinSize, with the state
// and log records hicmatrix keeps), the csr arrays of its matrix, the dicts
// and lists the methods read and change in place, and the h5py group that
// writeInteractionFileHDF5 and createUniqueHDFGroup write into. Scalar
// arithmetic whose result type the caller sees (int, numpy.int64,
// numpy.float64 depending on the operands) is done with Python's number
// protocol.
//
// Behaviour kept from 3.7.6 that looks wrong, each where it happens:
//  * createUniqueHDFGroup never retries another name (file_name is not
//    updated), so a name that exists loops forever;
//  * readBackgroundDataFile extends the model with numpy.int64 keys and
//    shares one value list between the extended keys;
//  * getDataForPlotting drops the reference point from the data when its
//    width does not exceed the resolution, and removes the keys outside the
//    range from the caller's background model;
//  * merge_neighbors never keeps an unmerged last element and changes the
//    caller's lists in place;
//  * readInteractionFile's `list(keys) == 0` test is never true.

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "api.hpp"
#include "chic_common.hpp"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/py_text.hpp"
#include "hicx/scipy_special.hpp"

namespace hicx::python {

using namespace chic_api;
namespace hc = hicx::chic;

namespace {

constexpr const char* kLogger = "hicexplorer.lib.viewpoint";

void debug(const std::string& message) { log(kLogger, "debug", message); }

// Runs a core function, raising its hicx::chic exceptions as the Python ones.
template <class F>
decltype(auto) core(F&& f) {
    try {
        return f();
    } catch (const hc::TypeError& error) {
        raise_error("TypeError", error.what());
    } catch (const hc::IndexError& error) {
        raise_error("IndexError", error.what());
    } catch (const hc::ValueError& error) {
        raise_error("ValueError", error.what());
    }
}

// '{}'.format(x)
std::string fmt(py::handle x) { return utf8(format(x, "")); }

// ---------------------------------------------------------------------------
// numbers

// An ndarray or list of numbers as float64, and whether numpy holds it as
// float32 (np.mean and betainc then work in single precision).
struct Numbers {
    std::vector<double> values;
    bool float32 = false;
};

Numbers numbers(py::handle data) {
    Numbers out;
    py::array array;
    if (is_ndarray(data)) {
        array = py::reinterpret_borrow<py::array>(data);
    } else {
        array = py::reinterpret_borrow<py::array>(numpy().attr("asarray")(data));
    }
    if (dtype_of(array) == DType::Float32) {
        out.float32 = true;
        for (const float v : to_vector<float>(array)) {
            out.values.push_back(static_cast<double>(v));
        }
        return out;
    }
    const py::array as_double = py::reinterpret_borrow<py::array>(
        numpy().attr("asarray")(array, py::arg("dtype") = numpy().attr("float64")));
    out.values = to_vector<double>(as_double);
    return out;
}

// self.hicMatrix.matrix[row, column] for a csr (or csc) matrix: the sum of the
// stored entries at that position, as a double.
struct MatrixRows {
    Sparse matrix;
    std::vector<std::int64_t> indptr;
    std::vector<std::int64_t> indices;
    std::vector<double> data;

    explicit MatrixRows(py::handle object) : matrix(sparse(object, "computeViewpoint", "rc")) {
        indptr = matrix.indptr_vector();
        const py::array idx = py::reinterpret_borrow<py::array>(
            numpy().attr("asarray")(matrix.indices, py::arg("dtype") = numpy().attr("int64")));
        indices = to_vector<std::int64_t>(idx);
        const py::array values = py::reinterpret_borrow<py::array>(
            numpy().attr("asarray")(matrix.data, py::arg("dtype") = numpy().attr("float64")));
        data = to_vector<double>(values);
    }

    // scipy's checks of an integer index
    static std::int64_t checked(std::int64_t index, std::int64_t n, const char* what) {
        if (index < -n || index >= n) {
            raise_error("IndexError", std::string(what) + " index (" + std::to_string(index) +
                                          ") out of range");
        }
        return index < 0 ? index + n : index;
    }

    // Adds row `row` over the columns [first, first + count) to out.
    void add_row(std::int64_t row, std::int64_t first, std::int64_t count, std::vector<double>& out) const {
        const std::int64_t rows = matrix.rows;
        const std::int64_t cols = matrix.cols;
        if (count <= 0) {
            return;
        }
        const std::int64_t r = checked(row, rows, "row");
        std::vector<double> dense(static_cast<std::size_t>(count), 0.0);
        if (matrix.format == 'r') {
            for (std::int64_t k = indptr[r]; k < indptr[r + 1]; ++k) {
                const std::int64_t c = indices[k];
                if (c >= first && c < first + count) {
                    dense[static_cast<std::size_t>(c - first)] += data[k];
                }
            }
            for (std::int64_t j = 0; j < count; ++j) {
                (void)checked(first + j, cols, "column");
            }
        } else {
            for (std::int64_t j = 0; j < count; ++j) {
                const std::int64_t c = checked(first + j, cols, "column");
                for (std::int64_t k = indptr[c]; k < indptr[c + 1]; ++k) {
                    if (indices[k] == r) {
                        dense[static_cast<std::size_t>(j)] += data[k];
                    }
                }
            }
        }
        for (std::int64_t j = 0; j < count; ++j) {
            out[static_cast<std::size_t>(j)] += dense[static_cast<std::size_t>(j)];
        }
    }
};

// ---------------------------------------------------------------------------
// HDF5

// interactionFileHDF5Object[path][:] wrapped in the try block of
// readInteractionFile: a scalar dataset raises h5py's error, which is logged.
std::optional<hc::H5Value> read_array_logged(const hc::H5Reader& reader, const std::string& path) {
    const hc::H5Value value = reader.read(path);
    if (value.scalar) {
        debug(path);
        debug("Illegal slicing argument for scalar dataspace");
        return std::nullopt;
    }
    return value;
}

py::tuple read_interaction_file(py::handle self, py::handle pFilePath, py::handle pInternalIdentifierTriplet) {
    (void)self;
    hc::H5Reader reader = open_h5(pFilePath);
    const std::string internal_path = utf8(join("/", pInternalIdentifierTriplet));
    if (!reader.contains(internal_path)) {
        debug("internal path not found: " + internal_path);
        return py::make_tuple(py::dict(), py::dict(), py::list());
    }
    if (reader.kind(internal_path) == hc::H5Reader::Kind::Dataset) {
        raise_error("AttributeError", "'Dataset' object has no attribute 'keys'");
    }
    std::vector<hc::H5Value> data;
    for (const char* name : {"relative_position_list", "interaction_data_list", "pvalue", "raw", "xfold",
                             "start_list", "end_list"}) {
        const std::string path = internal_path + "/" + name;
        if (reader.contains(path)) {
            if (auto value = read_array_logged(reader, path)) {
                data.push_back(std::move(*value));
            }
        }
    }
    std::optional<py::object> chromosome;
    std::optional<py::object> gene;
    std::optional<py::object> sum_of_interactions;
    if (reader.contains(internal_path + "/chromosome")) {
        chromosome = decode_utf8(h5_scalar(reader.read(internal_path + "/chromosome")));
    } else {
        debug("internal path not found: " + internal_path + "/chromosome");
    }
    if (reader.contains(internal_path + "/gene")) {
        gene = decode_utf8(h5_scalar(reader.read(internal_path + "/gene")));
    } else {
        debug("internal path not found: " + internal_path + "/gene");
    }
    if (reader.contains(internal_path + "/sum_of_interactions")) {
        sum_of_interactions = h5_scalar(reader.read(internal_path + "/sum_of_interactions"));
    } else {
        debug("internal path not found: " + internal_path + "/sum_of_interactions");
    }
    py::object reference_point_start = py::none();
    py::object reference_point_end = py::none();
    if (reader.contains(internal_path + "/reference_point_start")) {
        reference_point_start = h5_scalar(reader.read(internal_path + "/reference_point_start"));
    }
    if (reader.contains(internal_path + "/reference_point_end")) {
        reference_point_end = h5_scalar(reader.read(internal_path + "/reference_point_end"));
    }
    reader.close();

    py::dict interaction_data;
    py::dict interaction_file_data;
    const auto empty = [] { return py::make_tuple(py::dict(), py::dict(), py::list()); };
    if (data.empty()) {
        return empty();
    }
    try {
        for (std::size_t i = 0; i < data[0].length; ++i) {
            for (std::size_t column = 1; column <= 4; ++column) {
                if (column >= data.size() || i >= data[column].length) {
                    return empty();
                }
            }
            const py::object key = h5_element(data[0], i);
            py::list values;
            for (std::size_t column = 1; column <= 4; ++column) {
                values.append(py::float_(h5_float(data[column], i)));
            }
            interaction_data[key] = values;
            if (!chromosome || data.size() < 7 || i >= data[5].length || i >= data[6].length) {
                return empty();
            }
            const auto truncated = [&](std::size_t column) {
                return int_of(py::float_(h5_float(data[column], i)));
            };
            py::list line;
            line.append(str_object(*chromosome));
            line.append(truncated(5));
            line.append(truncated(6));
            if (!gene) {
                return empty();
            }
            line.append(str_object(*gene));
            if (!sum_of_interactions) {
                return empty();
            }
            line.append(float_of(*sum_of_interactions));
            line.append(truncated(0));
            line.append(py::float_(h5_float(data[1], i)));
            line.append(py::float_(h5_float(data[2], i)));
            line.append(py::float_(h5_float(data[4], i)));
            line.append(py::float_(h5_float(data[3], i)));
            interaction_file_data[key] = line;
        }
    } catch (py::error_already_set& error) {
        if (!is_exception(error)) {
            throw;
        }
        return empty();
    } catch (const hicx::PythonError&) {
        return empty();
    }
    py::list reference_point;
    reference_point.append(reference_point_start);
    reference_point.append(reference_point_end);
    return py::make_tuple(interaction_data, interaction_file_data, reference_point);
}

// ---------------------------------------------------------------------------
// shared pieces of the methods

// self.hicMatrix.getRegionBinRange(...)
py::object region_bin_range(py::handle self, py::handle chromosome, py::handle start, py::handle end) {
    return self.attr("hicMatrix").attr("getRegionBinRange")(chromosome, start, end);
}

// The width of the reference point's genomic coordinates and whether it
// exceeds the resolution: `np.abs(int(start) - int(end)) > pResolution`.
struct Width {
    bool wide = false;
    py::object peak_width;  // np.abs(...) // pResolution
};

Width reference_width(py::handle coordinates, py::handle resolution) {
    const std::vector<py::object> parts = unpack(coordinates, 3);
    const std::int64_t start = to_int64(int_of(parts[1]));
    const std::int64_t end = to_int64(int_of(parts[2]));
    const py::object width = np_int64(std::llabs(start - end));
    Width out;
    out.wide = gt(width, resolution);
    if (out.wide) {
        out.peak_width = floordiv(width, resolution);
    }
    return out;
}

// `while i < peak_width: ... i += 1`
std::int64_t repetitions(const py::object& peak_width) {
    std::int64_t count = 0;
    py::object i = py::int_(0);
    while (lt(i, peak_width)) {
        ++count;
        i = add(i, py::int_(1));
    }
    return count;
}

}  // namespace

void bind_lib_viewpoint(py::module_& m) {
    m.doc() = "C++ bodies of the methods of hicexplorer.lib.viewpoint.Viewpoint.";

    m.def("Viewpoint___init__", [](py::object self, py::object pHiCMatrix) {
        self.attr("hicMatrix") = pHiCMatrix;
    }, py::arg("self"), py::arg("pHiCMatrix") = py::none());

    m.def("Viewpoint_readReferencePointFile", [](py::object self, py::handle pBedFile, py::handle pGene) -> py::object {
        (void)self;
        py::list viewpoints;
        py::list gene_list;
        const bool with_gene = truthy(pGene);
        for (const std::string& line : read_text_lines(pBedFile)) {
            const std::vector<std::string_view> fields = hicx::pytext::split(hicx::pytext::strip(line), "\t");
            if (line.empty() || fields.empty()) {
                continue;
            }
            py::object chrom;
            py::object start;
            py::object end;
            if (fields.size() == 3) {
                chrom = to_str(fields[0]);
                start = to_str(fields[1]);
                end = to_str(fields[1]);
                if (with_gene) {
                    gene_list.append(to_str(fields[2]));
                }
            } else if (fields.size() > 3) {
                chrom = to_str(fields[0]);
                start = to_str(fields[1]);
                end = to_str(fields[2]);
                if (with_gene) {
                    gene_list.append(to_str(fields[3]));
                }
            } else {
                continue;
            }
            viewpoints.append(py::make_tuple(chrom, start, end));
        }
        if (with_gene) {
            return py::make_tuple(viewpoints, gene_list);
        }
        return viewpoints;
    }, py::arg("self"), py::arg("pBedFile"), py::arg("pGene") = true);

    m.def("Viewpoint_readInteractionFile", &read_interaction_file,
          py::arg("self"), py::arg("pFilePath"), py::arg("pInternalIdentifierTriplet"));

    m.def("Viewpoint_interactionDataForPlot", [](py::object self, py::handle pFilePath, py::handle pInternalIdentifierTriplet) {
        const std::vector<py::object> read = unpack(self.attr("readInteractionFile")(pFilePath, pInternalIdentifierTriplet), 3);
        const py::object interaction_file_data = read[1];
        const py::object reference_point = read[2];
        py::dict interaction_data;
        py::dict p_score;
        py::dict genomic_coordinates;
        for (const py::handle item : interaction_file_data.attr("items")()) {
            const std::vector<py::object> kv = unpack(item, 2);
            interaction_data[kv[0]] = getitem(kv[1], 6);
            p_score[kv[0]] = getitem(kv[1], 7);
            genomic_coordinates[kv[0]] = getitem(kv[1], py::slice(py::none(), py::int_(3), py::none()));
        }
        return py::make_tuple(interaction_data, p_score, genomic_coordinates, reference_point);
    }, py::arg("self"), py::arg("pFilePath"), py::arg("pInternalIdentifierTriplet"));

    m.def("Viewpoint_readBackgroundDataFile", [](py::object self, py::handle pBedFile, py::handle pRange,
                                                 py::handle pFixateRange, py::handle pMean) {
        (void)self;
        py::dict distance;
        const std::vector<std::string> lines = read_text_lines(pBedFile);
        const bool mean = truthy(pMean);
        for (std::size_t n = 1; n < lines.size(); ++n) {
            const std::vector<std::string_view> fields = hicx::pytext::split(lines[n], "\t");
            const auto field = [&](py::ssize_t index) -> py::object {
                const py::ssize_t size = static_cast<py::ssize_t>(fields.size());
                const py::ssize_t position = index < 0 ? index + size : index;
                if (position < 0 || position >= size) {
                    raise_error("IndexError", "list index out of range");
                }
                return to_str(fields[static_cast<std::size_t>(position)]);
            };
            if (!mean) {
                const py::object key = parse_int(field(0));
                py::list values;
                values.append(float_of(field(1)));
                values.append(float_of(field(2)));
                values.append(float_of(field(3)));
                distance[key] = values;
            } else {
                const py::object key = parse_int(field(0));
                py::list values;
                values.append(float_of(field(-1)));
                distance[key] = values;
            }
        }
        py::object max_key = extreme_key(distance, true);
        py::object min_key = extreme_key(distance, false);
        if (gt(max_key, pFixateRange)) {
            max_key = py::reinterpret_borrow<py::object>(pFixateRange);
        }
        const py::object negative_fixate = neg(pFixateRange);
        if (lt(min_key, negative_fixate)) {
            min_key = negative_fixate;
        }
        const py::list keys = as_list(distance);
        if (keys.size() < 2) {
            raise_error("IndexError", "list index out of range");
        }
        // np.absolute(np.absolute(keys[0]) - np.absolute(keys[1])): numpy.int64
        const std::int64_t increment = std::llabs(std::llabs(to_int64(keys[0])) - std::llabs(to_int64(keys[1])));
        const py::object inc = np_int64(increment);
        const py::object upper = getitem(pRange, 1);
        const py::object lower = neg(getitem(pRange, 0));
        if (lt(max_key, upper)) {
            py::object i = max_key;
            while (lt(i, upper)) {
                i = add(i, inc);
                distance[i] = getitem(distance, max_key);
            }
        }
        if (gt(min_key, lower)) {
            py::object i = min_key;
            while (gt(i, lower)) {
                i = sub(i, inc);
                distance[i] = getitem(distance, min_key);
            }
        }
        return distance;
    }, py::arg("self"), py::arg("pBedFile"), py::arg("pRange"), py::arg("pFixateRange"), py::arg("pMean") = false);

    m.def("Viewpoint_writeInteractionFile", [](py::object self, py::handle pBedFile, py::handle pData, py::handle pHeader,
                                               py::handle pPValueData, py::handle pXfold, py::handle pDecimalPlaces) {
        (void)self;
        const py::object path = add(pBedFile, py::str(".txt")).attr("strip")();
        const py::object spec = py::str("." + fmt(pDecimalPlaces) + "f");
        with_open(path, "w", [&](py::object& fh) {
            write(fh, fmt(pHeader) + "\n");
            py::ssize_t j = 0;
            for (const py::handle interaction : pData) {
                std::string line;
                for (py::ssize_t k = 0; k < 6; ++k) {
                    line += fmt(getitem(interaction, k)) + "\t";
                }
                line += utf8(format(getitem(interaction, 6), spec)) + "\t";
                line += utf8(format(getitem(pPValueData, j), spec)) + "\t";
                line += utf8(format(getitem(pXfold, j), spec)) + "\t";
                line += utf8(format(getitem(interaction, 7), spec)) + "\n";
                write(fh, line);
                ++j;
            }
        });
    }, py::arg("self"), py::arg("pBedFile"), py::arg("pData"), py::arg("pHeader"), py::arg("pPValueData"),
       py::arg("pXfold"), py::arg("pDecimalPlaces") = 12);

    m.def("Viewpoint_writeInteractionFileHDF5", [](py::object self, py::handle pInteractionFileGroupH5Object,
                                                   py::handle pFileName, py::handle pData, py::handle pReferencePoint) {
        const std::vector<py::object> created =
            unpack(self.attr("createUniqueHDFGroup")(pInteractionFileGroupH5Object, pFileName), 2);
        const py::object group = created[0];
        const py::object file_name = created[1];
        setitem(group, py::str("chromosome"), str_object(getitem(pData, 0)));
        const auto compressed = [&](const char* name, py::ssize_t index) {
            group.attr("create_dataset")(name, py::arg("data") = getitem(pData, index),
                                         py::arg("compression") = "gzip", py::arg("compression_opts") = 9);
        };
        compressed("start_list", 1);
        compressed("end_list", 2);
        group.attr("create_dataset")("gene", py::arg("data") = getitem(pData, 3));
        group.attr("create_dataset")("sum_of_interactions", py::arg("data") = getitem(pData, 4));
        compressed("relative_position_list", 5);
        compressed("interaction_data_list", 6);
        compressed("pvalue", 7);
        compressed("xfold", 8);
        compressed("raw", 9);
        group.attr("create_dataset")("reference_point_start", py::arg("data") = int_of(getitem(pReferencePoint, 0)));
        group.attr("create_dataset")("reference_point_end", py::arg("data") = int_of(getitem(pReferencePoint, 1)));
        return file_name;
    }, py::arg("self"), py::arg("pInteractionFileGroupH5Object"), py::arg("pFileName"), py::arg("pData"),
       py::arg("pReferencePoint"));

    m.def("Viewpoint_createUniqueHDFGroup", [](py::object self, py::handle pGroupObject, py::handle pAdditionalGroupName) {
        (void)self;
        py::object name = py::reinterpret_borrow<py::object>(pAdditionalGroupName);
        py::object file_name;
        py::object group;
        std::int64_t counter = 0;
        // 3.7.6 never retries another name: file_name keeps the first one,
        // so an existing group makes this loop forever, as it did.
        while (true) {
            if (counter != 0) {
                name = add(add(name, py::str("_")), py::str(std::to_string(counter)));
            } else {
                file_name = name;
            }
            try {
                group = pGroupObject.attr("create_group")(file_name);
                break;
            } catch (py::error_already_set& error) {
                if (!error.matches(PyExc_ValueError)) {
                    throw;
                }
                ++counter;
            }
        }
        return py::make_tuple(group, file_name);
    }, py::arg("self"), py::arg("pGroupObject"), py::arg("pAdditionalGroupName"));

    m.def("Viewpoint_computeViewpoint", [](py::object self, py::handle pReferencePoint, py::handle pChromViewpoint,
                                           py::handle pRegion_start, py::handle pRegion_end) {
        const std::vector<py::object> indices = unpack(self.attr("getReferencePointAsMatrixIndices")(pReferencePoint), 2);
        const py::object range_object =
            self.attr("getViewpointRangeAsMatrixIndices")(pChromViewpoint, pRegion_start, pRegion_end);
        py::list view_point_range = as_list(range_object);
        setitem(view_point_range, 1, add(getitem(view_point_range, 1), py::int_(1)));
        const std::int64_t view_point_start = to_int64(indices[0]);
        const std::int64_t view_point_end = to_int64(indices[1]);
        const std::int64_t range_start = to_int64(getitem(view_point_range, 0));
        const std::int64_t range_end = to_int64(getitem(view_point_range, 1));
        const std::int64_t elements = range_end - range_start;
        if (elements < 0) {
            raise_error("ValueError", "negative dimensions are not allowed");
        }
        std::vector<double> data(static_cast<std::size_t>(elements), 0.0);
        const py::object hic = self.attr("hicMatrix");
        if (view_point_start <= view_point_end) {
            const MatrixRows rows(hic.attr("matrix"));
            for (std::int64_t row = view_point_start; row <= view_point_end; ++row) {
                (void)unpack(hic.attr("getBinPos")(row), 4);
                rows.add_row(row, range_start, elements, data);
            }
        }
        const hc::ComputedViewpoint result = core([&] {
            return hc::collapse_viewpoint(data, view_point_start, view_point_end, range_start);
        });
        // view_point_start - view_point_range[0], in the operands' type
        const py::object index_before_viewpoint = sub(indices[0], getitem(view_point_range, 0));
        return py::make_tuple(from_vector(result.data), index_before_viewpoint);
    }, py::arg("self"), py::arg("pReferencePoint"), py::arg("pChromViewpoint"), py::arg("pRegion_start"),
       py::arg("pRegion_end"));

    // createInteractionFileData and createInteractionFileDataHDF5 walk the same
    // positions; `hdf5` selects the result layout.
    const auto interaction_positions = [](py::object self, py::handle pReferencePoint, py::handle pChromViewpoint,
                                          py::handle pRegion_start, py::handle pRegion_end,
                                          py::handle pInteractionData, py::handle pInteractionDataRaw,
                                          const std::function<void(py::handle, py::handle, py::handle, py::handle,
                                                                   py::ssize_t)>& emit) {
        const std::vector<py::object> indices = unpack(self.attr("getReferencePointAsMatrixIndices")(pReferencePoint), 2);
        py::list view_point_range =
            as_list(self.attr("getViewpointRangeAsMatrixIndices")(pChromViewpoint, pRegion_start, pRegion_end));
        setitem(view_point_range, 1, add(getitem(view_point_range, 1), py::int_(1)));
        const py::object hic = self.attr("hicMatrix");
        const std::vector<py::object> first = unpack(hic.attr("getBinPos")(indices[0]), 4);
        const std::vector<py::object> last = unpack(hic.attr("getBinPos")(indices[1]), 4);
        const py::object start = first[1];
        const py::object end = last[2];
        const std::int64_t view_point_start = to_int64(indices[0]);
        const std::int64_t view_point_end = to_int64(indices[1]);
        std::vector<std::int64_t> positions;
        for (std::int64_t i = to_int64(getitem(view_point_range, 0)); i < view_point_start; ++i) {
            positions.push_back(i);
        }
        positions.push_back(view_point_start);
        for (std::int64_t i = view_point_end + 1; i < to_int64(getitem(view_point_range, 1)); ++i) {
            positions.push_back(i);
        }
        const py::ssize_t n = length(pInteractionData);
        std::int64_t relative_position = -1;
        for (py::ssize_t j = 0; j < n && static_cast<std::size_t>(j) < positions.size(); ++j) {
            const std::int64_t idx = positions[static_cast<std::size_t>(j)];
            try {
                const std::vector<py::object> bin = unpack(hic.attr("getBinPos")(idx), 4);
                if (relative_position < 0) {
                    relative_position = to_int64(int_of(bin[1])) - to_int64(int_of(start));
                } else {
                    relative_position = to_int64(int_of(bin[2])) - to_int64(int_of(end));
                }
                (void)pInteractionDataRaw;
                emit(bin[0], bin[1], bin[2], py::int_(relative_position), j);
            } catch (py::error_already_set& error) {
                if (!is_exception(error)) {
                    throw;
                }
                log(kLogger, "error", "Failed to get bin position of index " + std::to_string(idx));
                python_exit(1);
            } catch (const hicx::PythonError&) {
                log(kLogger, "error", "Failed to get bin position of index " + std::to_string(idx));
                python_exit(1);
            }
        }
    };

    m.def("Viewpoint_createInteractionFileData", [interaction_positions](
              py::object self, py::handle pReferencePoint, py::handle pChromViewpoint, py::handle pRegion_start,
              py::handle pRegion_end, py::handle pInteractionData, py::handle pInteractionDataRaw, py::handle pGene,
              py::handle pSumOfInteractions) {
        py::list interactions_list;
        interaction_positions(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData,
                              pInteractionDataRaw,
                              [&](py::handle chrom, py::handle start, py::handle end, py::handle relative,
                                  py::ssize_t j) {
                                  interactions_list.append(py::make_tuple(
                                      chrom, start, end, pGene, str_object(pSumOfInteractions), relative,
                                      float_of(getitem(pInteractionData, j)), float_of(getitem(pInteractionDataRaw, j))));
                              });
        return interactions_list;
    }, py::arg("self"), py::arg("pReferencePoint"), py::arg("pChromViewpoint"), py::arg("pRegion_start"),
       py::arg("pRegion_end"), py::arg("pInteractionData"), py::arg("pInteractionDataRaw"), py::arg("pGene"),
       py::arg("pSumOfInteractions"));

    m.def("Viewpoint_createInteractionFileDataHDF5", [interaction_positions](
              py::object self, py::handle pReferencePoint, py::handle pChromViewpoint, py::handle pRegion_start,
              py::handle pRegion_end, py::handle pInteractionData, py::handle pInteractionDataRaw, py::handle pGene,
              py::handle pSumOfInteractions, py::handle pPValueList, py::handle xFoldList) {
        py::object chrom = py::str("");
        py::list start_list;
        py::list end_list;
        py::list relative_position_list;
        py::list interaction_data_list;
        py::list raw_data_list;
        interaction_positions(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData,
                              pInteractionDataRaw,
                              [&](py::handle chrom_second, py::handle start, py::handle end, py::handle relative,
                                  py::ssize_t j) {
                                  chrom = py::reinterpret_borrow<py::object>(chrom_second);
                                  start_list.append(start);
                                  end_list.append(end);
                                  relative_position_list.append(relative);
                                  interaction_data_list.append(float_of(getitem(pInteractionData, j)));
                                  raw_data_list.append(float_of(getitem(pInteractionDataRaw, j)));
                              });
        py::list out;
        out.append(chrom);
        out.append(start_list);
        out.append(end_list);
        out.append(pGene);
        out.append(pSumOfInteractions);
        out.append(relative_position_list);
        out.append(interaction_data_list);
        out.append(pPValueList);
        out.append(xFoldList);
        out.append(raw_data_list);
        return out;
    }, py::arg("self"), py::arg("pReferencePoint"), py::arg("pChromViewpoint"), py::arg("pRegion_start"),
       py::arg("pRegion_end"), py::arg("pInteractionData"), py::arg("pInteractionDataRaw"), py::arg("pGene"),
       py::arg("pSumOfInteractions"), py::arg("pPValueList"), py::arg("xFoldList"));

    m.def("Viewpoint_getViewpointRangeAsMatrixIndices", [](py::object self, py::handle pChromViewpoint,
                                                           py::handle pRegion_start, py::handle pRegion_end) {
        return region_bin_range(self, pChromViewpoint, pRegion_start, pRegion_end);
    }, py::arg("self"), py::arg("pChromViewpoint"), py::arg("pRegion_start"), py::arg("pRegion_end"));

    m.def("Viewpoint_getReferencePointAsMatrixIndices", [](py::object self, py::handle pReferencePoint) {
        const py::ssize_t n = length(pReferencePoint);
        std::vector<py::object> range;
        if (n == 2) {
            const py::object position = int_of(getitem(pReferencePoint, 1));
            range = unpack(region_bin_range(self, getitem(pReferencePoint, 0), position, int_of(getitem(pReferencePoint, 1))), 2);
        } else if (n == 3) {
            const py::object start = int_of(getitem(pReferencePoint, 1));
            range = unpack(region_bin_range(self, getitem(pReferencePoint, 0), start, int_of(getitem(pReferencePoint, 2))), 2);
        } else {
            log(kLogger, "error", "No valid reference point given. " + fmt(pReferencePoint));
            python_exit(1);
        }
        return py::make_tuple(range[0], range[1]);
    }, py::arg("self"), py::arg("pReferencePoint"));

    m.def("Viewpoint_smoothInteractionValues", [](py::object self, py::handle pData, py::handle pWindowSize) {
        (void)self;
        // window_size = np.int32(np.floor(pWindowSize / 2)); even windows hold
        // one element fewer upstream. The core takes the window that has the
        // same half and parity.
        const double window = to_double(pWindowSize);
        const std::int64_t half = static_cast<std::int32_t>(std::floor(window / 2.0));
        const bool even = eq(steal(PyNumber_Remainder(pWindowSize.ptr(), py::int_(2).ptr())), py::int_(0));
        const std::int64_t equivalent = 2 * half + (even ? 0 : 1);
        const py::ssize_t n = length(pData);
        (void)n;
        const Numbers data = numbers(pData);
        std::vector<double> result;
        if (data.float32) {
            std::vector<float> single(data.values.begin(), data.values.end());
            result = core([&] { return hc::smooth_interaction_values(std::span<const float>(single), equivalent); });
        } else {
            result = core([&] { return hc::smooth_interaction_values(std::span<const double>(data.values), equivalent); });
        }
        return from_vector(result);
    }, py::arg("self"), py::arg("pData"), py::arg("pWindowSize"));

    m.def("Viewpoint_computeRelativeValues", [](py::object self, py::handle pData, py::handle pDenominator) {
        (void)self;
        py::array output = py::reinterpret_borrow<py::array>(
            numpy().attr("array")(pData, py::arg("dtype") = builtin("float")));
        const std::vector<double> values = to_vector<double>(output);
        const double denominator = truthy(pDenominator) ? to_double(pDenominator) : 0.0;
        const std::vector<double> result = hc::compute_relative_values(std::span<const double>(values), denominator);
        write_back(output, result);
        return output;
    }, py::arg("self"), py::arg("pData"), py::arg("pDenominator") = py::none());

    m.def("Viewpoint_calculateViewpointRange", [](py::object self, py::handle pViewpoint, py::handle pRange) {
        const py::object hic = self.attr("hicMatrix");
        const py::object chr_bin_range = hic.attr("getChrBinRange")(getitem(pViewpoint, 0));
        if (chr_bin_range.is_none()) {
            return py::make_tuple(py::none(), py::none(), py::none());
        }
        const py::object max_length = getitem(hic.attr("getBinPos")(sub(getitem(chr_bin_range, 1), py::int_(1))), 2);
        const py::object bin_size = hic.attr("getBinSize")();
        py::list range;
        range.append(getitem(pRange, 0));
        range.append(getitem(pRange, 1));
        py::object region_start = sub(int_of(getitem(pViewpoint, 1)), getitem(pRange, 0));
        if (lt(region_start, py::int_(0))) {
            region_start = py::int_(0);
            setitem(range, 0, int_of(getitem(pViewpoint, 1)));
        }
        py::object region_end = add(int_of(getitem(pViewpoint, 2)), getitem(pRange, 1));
        if (gt(region_end, max_length)) {
            // -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
            region_end = sub(max_length, py::int_(1));
            setitem(range, 1, add(sub(max_length, int_of(getitem(pViewpoint, 2))), bin_size));
        }
        return py::make_tuple(region_start, region_end, range);
    }, py::arg("self"), py::arg("pViewpoint"), py::arg("pRange"));

    m.def("Viewpoint_getDataForPlotting", [](py::object self, py::handle pFilePath, py::handle pInteractionFile,
                                             py::handle pRange, py::handle pBackgroundModel, py::handle pResolution) {
        const std::vector<py::object> plot = unpack(self.attr("interactionDataForPlot")(pFilePath, pInteractionFile), 4);
        const py::object interaction_data = plot[0];
        const py::object p_value_data = plot[1];
        const py::object genomic_coordinates = plot[2];
        const py::object viewpoint = plot[3];

        py::list data;
        py::list p_value;
        py::object data_background = py::none();
        py::object viewpoint_index_end = py::none();
        py::object viewpoint_index_start;
        std::optional<py::list> background_data_keys_sorted;
        const auto in_range = [&](py::handle key) {
            return ge(key, neg(getitem(pRange, 0))) && le(key, getitem(pRange, 1));
        };
        if (truthy(pRange)) {
            for (const py::handle key : as_list(interaction_data)) {
                if (in_range(key)) {
                    continue;
                }
                interaction_data.attr("pop")(key, py::none());
            }
            if (truthy(pBackgroundModel)) {
                for (const py::handle key : sorted_keys(pBackgroundModel)) {
                    if (in_range(key)) {
                        continue;
                    }
                    pBackgroundModel.attr("pop")(key, py::none());
                }
                background_data_keys_sorted = sorted_keys(pBackgroundModel);
            }
        }
        const py::int_ zero(0);
        const auto add_end = [&]() {
            if (viewpoint_index_end.is_none()) {
                viewpoint_index_end = viewpoint_index_start;
            } else {
                viewpoint_index_end = add(viewpoint_index_end, viewpoint_index_start);
            }
        };
        if (truthy(pBackgroundModel)) {
            if (!background_data_keys_sorted) {
                raise_error("UnboundLocalError",
                            "cannot access local variable 'background_data_keys_sorted' where it is not "
                            "associated with a value");
            }
            viewpoint_index_start = background_data_keys_sorted->attr("index")(zero);
            py::list background;
            data_background = background;
            for (const py::handle key : *background_data_keys_sorted) {
                if (!contains(interaction_data, key)) {
                    continue;
                }
                if (eq(key, zero)) {
                    const Width width = reference_width(getitem(genomic_coordinates, key), pResolution);
                    if (width.wide) {
                        viewpoint_index_end = width.peak_width;
                        for (std::int64_t i = repetitions(width.peak_width); i > 0; --i) {
                            data.append(getitem(interaction_data, key));
                        }
                    }
                } else {
                    data.append(getitem(interaction_data, key));
                }
                if (contains(p_value_data, key)) {
                    if (eq(key, zero)) {
                        const Width width = reference_width(getitem(genomic_coordinates, key), pResolution);
                        if (width.wide) {
                            for (std::int64_t i = repetitions(width.peak_width); i > 0; --i) {
                                p_value.append(getitem(p_value_data, key));
                            }
                        }
                    } else {
                        p_value.append(getitem(p_value_data, key));
                    }
                }
                if (eq(key, zero)) {
                    const Width width = reference_width(getitem(genomic_coordinates, key), pResolution);
                    if (width.wide) {
                        for (std::int64_t i = repetitions(width.peak_width); i > 0; --i) {
                            background.append(getitem(getitem(pBackgroundModel, key), 0));
                            debug("pBackgroundModel[key][0] " + fmt(getitem(getitem(pBackgroundModel, key), 0)));
                            debug("peak_width " + fmt(width.peak_width));
                            debug("data_background[-1] " + fmt(background[background.size() - 1]));
                        }
                    }
                } else {
                    background.append(getitem(getitem(pBackgroundModel, key), 0));
                }
            }
            add_end();
        } else {
            data = py::list();
            const py::list interaction_key = sorted_keys(interaction_data);
            for (const py::handle key : interaction_key) {
                if (eq(key, zero)) {
                    const Width width = reference_width(getitem(genomic_coordinates, key), pResolution);
                    if (width.wide) {
                        viewpoint_index_end = width.peak_width;
                        for (std::int64_t i = repetitions(width.peak_width); i > 0; --i) {
                            data.append(getitem(interaction_data, key));
                        }
                    }
                } else {
                    data.append(getitem(interaction_data, key));
                }
                if (contains(p_value_data, key)) {
                    if (eq(key, zero)) {
                        const Width width = reference_width(getitem(genomic_coordinates, key), pResolution);
                        if (width.wide) {
                            for (std::int64_t i = repetitions(width.peak_width); i > 0; --i) {
                                p_value.append(getitem(p_value_data, key));
                            }
                        }
                    } else {
                        p_value.append(getitem(p_value_data, key));
                    }
                }
            }
            viewpoint_index_start = interaction_key.attr("index")(zero);
            add_end();
        }
        return py::make_tuple(data, data_background, p_value, viewpoint_index_start, viewpoint_index_end, viewpoint);
    }, py::arg("self"), py::arg("pFilePath"), py::arg("pInteractionFile"), py::arg("pRange"),
       py::arg("pBackgroundModel"), py::arg("pResolution"));

    // The data the drawing methods plotBackgroundModel and plotPValue compute.
    m.def("Viewpoint_backgroundModelBand", [](py::handle pBackgroundData, py::handle pXFold) {
        // pBackgroundData = np.array(pBackgroundData) and, when pXFold is set,
        // the upper border of the band, pBackgroundData * pXFold
        const py::array array = py::reinterpret_borrow<py::array>(numpy().attr("array")(pBackgroundData));
        if (!truthy(pXFold)) {
            return py::make_tuple(array, py::none());
        }
        const DType dtype = dtype_of(array);
        const bool integer_factor = is_integer(pXFold) && !PyFloat_Check(pXFold.ptr());
        py::array upper;
        if (dtype == DType::Float64 || ((dtype == DType::Int64 || dtype == DType::Int32) && !integer_factor)) {
            const Numbers values = numbers(array);
            std::vector<double> product = values.values;
            const double factor = to_double(pXFold);
            for (double& value : product) {
                value *= factor;
            }
            upper = from_vector(product);
        } else if (dtype == DType::Float32) {
            const Numbers values = numbers(array);
            std::vector<float> product;
            const float factor = static_cast<float>(to_double(pXFold));
            for (const double value : values.values) {
                product.push_back(static_cast<float>(value) * factor);
            }
            upper = from_vector(product);
        } else if (dtype == DType::Int64 && integer_factor) {
            std::vector<std::int64_t> product = to_vector<std::int64_t>(array);
            const std::int64_t factor = to_int64(pXFold);
            for (std::int64_t& value : product) {
                value *= factor;
            }
            upper = from_vector(product);
        } else {
            raise_error("TypeError", "plotBackgroundModel: background data of dtype " + dtype_name(array) +
                                         " is not supported by HiCExplorer 4");
        }
        upper = py::reinterpret_borrow<py::array>(upper.attr("reshape")(array.attr("shape")));
        return py::make_tuple(array, py::object(upper));
    }, py::arg("pBackgroundData"), py::arg("pXFold"));

    m.def("Viewpoint_pValueRows", [](py::handle pPValueData) {
        // _z_score = np.empty([2, len(pPValueData)]); _z_score[:, :] = pPValueData
        const py::ssize_t n = length(pPValueData);
        const Numbers values = numbers(pPValueData);
        if (static_cast<py::ssize_t>(values.values.size()) != n) {
            raise_error("ValueError", "could not broadcast input array");
        }
        py::array_t<double> out({py::ssize_t{2}, n});
        double* p = out.mutable_data();
        for (py::ssize_t row = 0; row < 2; ++row) {
            for (py::ssize_t i = 0; i < n; ++i) {
                p[row * n + i] = values.values[static_cast<std::size_t>(i)];
            }
        }
        return out;
    }, py::arg("pPValueData"));

    m.def("Viewpoint_pValueTicks", [](py::handle pPValueData) {
        // [min(pPValueData), max(pPValueData)]
        py::list items = as_list(pPValueData);
        if (items.empty()) {
            raise_error("ValueError", "min() iterable argument is empty");
        }
        py::object lowest = items[0];
        py::object highest = items[0];
        for (std::size_t i = 1; i < items.size(); ++i) {
            if (lt(items[i], lowest)) {
                lowest = items[i];
            }
            if (gt(items[i], highest)) {
                highest = items[i];
            }
        }
        py::list ticks;
        ticks.append(lowest);
        ticks.append(highest);
        return ticks;
    }, py::arg("pPValueData"));

    m.def("Viewpoint_pValueLevels", [](py::handle pValueSignificanceLevels) {
        // pValueSignificanceLevels.insert(0, -1); .append(1), in place
        pValueSignificanceLevels.attr("insert")(0, -1);
        pValueSignificanceLevels.attr("append")(1);
        py::list ticks;
        ticks.append(getitem(pValueSignificanceLevels, 1));
        ticks.append(getitem(pValueSignificanceLevels, 2));
        ticks.append(getitem(pValueSignificanceLevels, 3));
        return ticks;
    }, py::arg("pValueSignificanceLevels"));

    m.def("Viewpoint_writePlotData", [](py::object self, py::handle pInteractionFileDataRaw, py::handle pFileName,
                                        py::handle pBackgroundModel) {
        (void)self;
        const py::list keys = sorted_keys(pInteractionFileDataRaw);
        with_open(add(pFileName, py::str(".txt")), "w", [&](py::object& output_file) {
            write(output_file,
                  "#ChrViewpoint\tStart\tEnd\tGene\tChrInteraction\tStart\tEnd\tRelative position\tRelative "
                  "Interactions\trbz-score\tRaw");
            if (truthy(pBackgroundModel)) {
                write(output_file, "\tbackground model\tbackground model SEM\n");
            } else {
                write(output_file, "\n");
            }
            for (const py::handle key : keys) {
                const py::object array = getitem(pInteractionFileDataRaw, key);
                std::string line;
                for (py::ssize_t k = 0; k < 11; ++k) {
                    line += (k > 0 ? "\t" : "") + fmt(getitem(array, k));
                }
                write(output_file, line);
                if (truthy(pBackgroundModel)) {
                    write(output_file, "\t" + fmt(getitem(array, 11)) + "\t" + fmt(getitem(array, 12)) + "\n");
                } else {
                    write(output_file, "\n");
                }
            }
        });
    }, py::arg("self"), py::arg("pInteractionFileDataRaw"), py::arg("pFileName"), py::arg("pBackgroundModel"));

    m.def("Viewpoint_interactionBackgroundData", [](py::object self, py::handle pBackground, py::handle pRange) {
        (void)self;
        py::list background_model;
        for (const py::handle key : sorted_keys(pBackground)) {
            if (ge(key, neg(getitem(pRange, 0))) && le(key, getitem(pRange, 1))) {
                background_model.append(getitem(pBackground, key));
            }
        }
        return numpy().attr("array")(background_model);
    }, py::arg("self"), py::arg("pBackground"), py::arg("pRange"));

    m.def("Viewpoint_rbz_score", [](py::object self, py::handle pRelativeInteractions, py::handle pBackgroundModel,
                                    py::handle pBackgroundModelSEM) -> py::object {
        (void)self;
        const py::ssize_t n = length(pRelativeInteractions);
        const py::ssize_t nb = length(pBackgroundModel);
        const py::ssize_t ns = length(pBackgroundModelSEM);
        if (n != nb || n != ns) {
            log(kLogger, "info",
                "Computing of rbz-score failed, data is having different size. \nrelative interactions " +
                    std::to_string(n) + " background model " + std::to_string(nb) + " background model SEM " +
                    std::to_string(ns));
            return py::none();
        }
        if (!is_ndarray(pRelativeInteractions) && !is_ndarray(pBackgroundModel)) {
            raise_error("TypeError", "unsupported operand type(s) for -: '" + type_name(pRelativeInteractions) +
                                         "' and '" + type_name(pBackgroundModel) + "'");
        }
        const py::object np = numpy();
        const py::array a = py::reinterpret_borrow<py::array>(np.attr("asarray")(pRelativeInteractions));
        const py::array b = py::reinterpret_borrow<py::array>(np.attr("asarray")(pBackgroundModel));
        const py::array s = py::reinterpret_borrow<py::array>(np.attr("asarray")(pBackgroundModelSEM));
        // the dtype of pRelativeInteractions - pBackgroundModel, then of the
        // in-place division by the SEM
        const py::dtype difference = py::reinterpret_borrow<py::dtype>(np.attr("result_type")(a, b));
        const py::dtype quotient = py::reinterpret_borrow<py::dtype>(np.attr("result_type")(difference, s));
        const char kind = difference.kind();
        if (kind != 'f' || (difference.itemsize() != 4 && difference.itemsize() != 8)) {
            if (kind == 'i' || kind == 'u' || kind == 'b') {
                PyErr_SetObject(np.attr("core").attr("_exceptions").attr("_UFuncOutputCastingError").ptr(),
                                np.attr("core").attr("_exceptions").attr("_UFuncOutputCastingError")(
                                    np.attr("divide"), "same_kind",
                                    quotient.kind() == 'f' ? py::object(quotient) : np.attr("dtype")("float64"),
                                    difference, 0).ptr());
                throw_python_error();
            }
            raise_error("TypeError", "rbz_score: dtype " + str_of(difference) + " is not supported by HiCExplorer 4");
        }
        const bool single = difference.itemsize() == 4;
        const bool single_division = single && quotient.itemsize() == 4;
        const Numbers x = numbers(a);
        const Numbers y = numbers(b);
        const Numbers z = numbers(s);
        if (static_cast<py::ssize_t>(x.values.size()) != n || static_cast<py::ssize_t>(y.values.size()) != n ||
            static_cast<py::ssize_t>(z.values.size()) != n) {
            raise_error("TypeError", "rbz_score: only one dimensional data is supported by HiCExplorer 4");
        }
        std::vector<double> out(static_cast<std::size_t>(n));
        for (py::ssize_t i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i);
            double value;
            if (single) {
                const float d = static_cast<float>(x.values[k]) - static_cast<float>(y.values[k]);
                value = single_division ? static_cast<double>(d / static_cast<float>(z.values[k]))
                                        : static_cast<double>(static_cast<float>(static_cast<double>(d) / z.values[k]));
            } else {
                value = (x.values[k] - y.values[k]) / z.values[k];
            }
            if (std::isnan(value) || std::isinf(value)) {
                value = -1.0;
            }
            out[k] = value;
        }
        if (single) {
            std::vector<float> narrowed(out.begin(), out.end());
            return from_vector(narrowed);
        }
        return from_vector(out);
    }, py::arg("self"), py::arg("pRelativeInteractions"), py::arg("pBackgroundModel"), py::arg("pBackgroundModelSEM"));

    m.def("Viewpoint_readRejectedFile", [](py::object self, py::handle pFilePath, py::handle pDifferentialHighlightTriplet,
                                           py::handle pViewpointIndexStart, py::handle pViewpointIndexEnd,
                                           py::handle pResolution, py::handle pRange, py::handle pViewpoint) -> py::object {
        (void)self;
        if (pDifferentialHighlightTriplet.is_none() || length(pDifferentialHighlightTriplet) == 0) {
            return py::none();
        }
        hc::H5Reader reader = open_h5(pFilePath);
        const std::string internal_path = utf8(join("/", pDifferentialHighlightTriplet));
        if (!reader.contains(internal_path)) {
            debug("internal path not found: " + internal_path);
            return py::none();
        }
        if (reader.kind(internal_path) == hc::H5Reader::Kind::Dataset) {
            raise_error("AttributeError", "'Dataset' object has no attribute 'keys'");
        }
        hc::H5Value start_list;
        hc::H5Value end_list;
        hc::H5Value relative_distance_list;
        try {
            start_list = reader.read(internal_path + "/start_list");
            end_list = reader.read(internal_path + "/end_list");
            relative_distance_list = reader.read(internal_path + "/relative_distance_list");
            if (start_list.scalar || end_list.scalar || relative_distance_list.scalar) {
                return py::none();
            }
        } catch (const hicx::PythonError&) {
            return py::none();
        }
        if (start_list.length == 0) {
            return py::none();
        }
        py::list highlight_areas_list;
        const std::vector<py::object> reference = unpack(pViewpoint, 2);
        const std::size_t n = std::min({start_list.length, end_list.length, relative_distance_list.length});
        for (std::size_t i = 0; i < n; ++i) {
            const py::object start = h5_element(start_list, i);
            const py::object end = h5_element(end_list, i);
            const py::object relative_distance = int_of(h5_element(relative_distance_list, i));
            if (ge(relative_distance, neg(getitem(pRange, 0))) && le(relative_distance, getitem(pRange, 1))) {
                const py::object width = truediv(sub(end, start), pResolution);
                py::object relative_genomic;
                py::object viewpoint_index;
                if (lt(relative_distance, py::int_(0))) {
                    relative_genomic = sub(start, int_of(reference[0]));
                    viewpoint_index = py::reinterpret_borrow<py::object>(pViewpointIndexStart);
                } else {
                    relative_genomic = sub(start, int_of(reference[1]));
                    viewpoint_index = py::reinterpret_borrow<py::object>(pViewpointIndexEnd);
                }
                const py::object relative_position = add(viewpoint_index, truediv(relative_genomic, pResolution));
                py::list area;
                area.append(relative_position);
                area.append(add(relative_position, width));
                highlight_areas_list.append(area);
            }
        }
        if (highlight_areas_list.empty()) {
            return py::none();
        }
        debug("highlight_areas_list " + fmt(highlight_areas_list));
        return highlight_areas_list;
    }, py::arg("self"), py::arg("pFilePath"), py::arg("pDifferentialHighlightTriplet"), py::arg("pViewpointIndexStart"),
       py::arg("pViewpointIndexEnd"), py::arg("pResolution"), py::arg("pRange"), py::arg("pViewpoint"));

    m.def("Viewpoint_pvalues", [](py::object self, py::handle pBackgroundModel, py::handle pDataList,
                                  py::handle pIndexReferencePoint) {
        (void)self;
        std::vector<double> p_value_list;
        py::ssize_t i = 0;
        const bool single = is_ndarray(pDataList) &&
                            dtype_of(py::reinterpret_borrow<py::array>(pDataList)) == DType::Float32;
        for (const py::handle data_element : pDataList) {
            py::object relative_distance = sub(py::int_(i), pIndexReferencePoint);
            if (!contains(pBackgroundModel, relative_distance)) {
                relative_distance = extreme_key(pBackgroundModel, !lt(relative_distance, py::int_(0)));
            }
            if (eq(data_element, py::float_(0.0))) {
                p_value_list.push_back(0.0);
            } else {
                const py::object model = getitem(pBackgroundModel, relative_distance);
                const double size = to_double(getitem(model, 0));
                const double probability = to_double(getitem(model, 1));
                const double x = to_double(data_element);
                const double shifted = single ? static_cast<double>(static_cast<float>(x) + 1.0F) : x + 1.0;
                p_value_list.push_back(hicx::scipy::betainc(size, shifted, probability));
            }
            ++i;
        }
        for (double& value : p_value_list) {
            value = 1.0 - value;
            if (std::isnan(value) || std::isinf(value)) {
                value = 1.0;
            }
        }
        return from_vector(p_value_list);
    }, py::arg("self"), py::arg("pBackgroundModel"), py::arg("pDataList"), py::arg("pIndexReferencePoint"));

    m.def("Viewpoint_merge_neighbors", [](py::object self, py::handle pScoresDictionary, py::handle pMergeThreshold) -> py::object {
        (void)self;
        if (pScoresDictionary.is_none() || length(pScoresDictionary) == 0) {
            debug("dict is None or empty");
            return py::none();
        }
        const py::list key_list = as_list(pScoresDictionary);
        const std::size_t n = key_list.size();
        std::vector<std::vector<py::object>> merge_ids;
        std::vector<py::object> non_merge;
        for (std::size_t i = 0; i + 1 < n; ++i) {
            const py::object key_pre = key_list[i];
            const py::object key_suc = key_list[i + 1];
            const std::int64_t pre = to_int64(int_of(getitem(getitem(pScoresDictionary, key_pre), 5)));
            const std::int64_t suc = to_int64(int_of(getitem(getitem(pScoresDictionary, key_suc), 5)));
            if (le(np_int64(std::llabs(pre - suc)), pMergeThreshold)) {
                if (!merge_ids.empty() && eq(merge_ids.back().back(), key_pre)) {
                    merge_ids.back().push_back(key_suc);
                } else {
                    merge_ids.push_back({key_pre, key_suc});
                }
            } else {
                // `if i == len(key_list) - 1` never holds
                if (!merge_ids.empty() && ne(merge_ids.back().back(), key_pre)) {
                    non_merge.push_back(key_pre);
                } else if (merge_ids.empty()) {
                    non_merge.push_back(key_pre);
                }
            }
        }
        py::dict scores_dict;
        py::dict merged_lines_dict;
        // np.array(list(map(float, line[-4:])))
        const auto last_four = [](py::handle line) {
            std::vector<double> values;
            for (const py::handle item : getitem(line, py::slice(py::int_(-4), py::none(), py::none()))) {
                values.push_back(to_double(float_of(item)));
            }
            return values;
        };
        const auto shape_error = [](std::size_t a, std::size_t b) {
            raise_error("ValueError", "operands could not be broadcast together with shapes (" + std::to_string(a) +
                                          ",) (" + std::to_string(b) + ",) (" + std::to_string(a) + ",)");
        };
        for (const auto& element : merge_ids) {
            py::list lines;
            lines.append(getitem(pScoresDictionary, element[0]));
            std::size_t index_maximum_element = 0;
            py::object base_element = getitem(pScoresDictionary, element[0]);
            std::vector<double> values = last_four(base_element);
            double max_value = to_double(float_of(getitem(base_element, -1)));
            for (std::size_t i = 1; i < element.size(); ++i) {
                const py::object line = getitem(pScoresDictionary, element[i]);
                setitem(base_element, -6, getitem(line, -6));
                const std::vector<double> more = last_four(line);
                if (more.size() != values.size() && more.size() != 1) {
                    shape_error(values.size(), more.size());
                }
                for (std::size_t k = 0; k < values.size(); ++k) {
                    values[k] += more.size() == 1 ? more[0] : more[k];
                }
                lines.append(line);
                const double raw = to_double(float_of(getitem(line, -1)));
                if (max_value < raw) {
                    max_value = raw;
                    index_maximum_element = i;
                }
            }
            base_element = getitem(pScoresDictionary, element[index_maximum_element]);
            for (std::size_t k = 0; k < 4; ++k) {
                if (k >= values.size()) {
                    raise_error("IndexError", "index " + std::to_string(k) + " is out of bounds for axis 0 with size " +
                                                  std::to_string(values.size()));
                }
                setitem(base_element, static_cast<py::ssize_t>(k) - 4, np_float64(values[k]));
            }
            setitem(base_element, 2, getitem(getitem(pScoresDictionary, element.back()), 2));
            setitem(base_element, 1, getitem(getitem(pScoresDictionary, element.front()), 1));
            scores_dict[element[index_maximum_element]] = base_element;
            merged_lines_dict[element[index_maximum_element]] = lines;
        }
        for (const py::object& key : non_merge) {
            scores_dict[key] = getitem(pScoresDictionary, key);
            py::list single;
            single.append(getitem(pScoresDictionary, key));
            merged_lines_dict[key] = single;
        }
        return py::make_tuple(scores_dict, merged_lines_dict);
    }, py::arg("self"), py::arg("pScoresDictionary"), py::arg("pMergeThreshold") = 1000);

    m.def("Viewpoint_readSignificantRegionsFile", [](py::object self, py::handle pFilePath,
                                                     py::handle pInternalIdentifierTriplet, py::handle pViewpointIndexStart,
                                                     py::handle pViewpointIndexEnd, py::handle pResolution,
                                                     py::handle pRange, py::handle pViewpoint) {
        const std::vector<py::object> read = unpack(self.attr("readInteractionFile")(pFilePath, pInternalIdentifierTriplet), 3);
        py::list highlight_areas_list;
        py::list p_values;
        std::vector<py::object> reference;
        if (length(pViewpoint) == 2) {
            reference = unpack(pViewpoint, 2);
        } else {
            debug("viewpoint_split " + fmt(pViewpoint) + ", file: " + fmt(pFilePath) + " " +
                  fmt(pInternalIdentifierTriplet));
            return py::make_tuple(py::none(), py::none());
        }
        for (const py::handle line : read[1].attr("values")()) {
            const py::object start = int_of(getitem(line, 1));
            const py::object end = int_of(getitem(line, 2));
            const py::object relative = int_of(getitem(line, 5));
            if (ge(relative, neg(getitem(pRange, 0))) && le(relative, getitem(pRange, 1))) {
                const py::object width = truediv(sub(end, start), pResolution);
                py::object relative_genomic;
                py::object viewpoint_index;
                if (lt(relative, py::int_(0))) {
                    relative_genomic = sub(start, int_of(reference[0]));
                    viewpoint_index = py::reinterpret_borrow<py::object>(pViewpointIndexStart);
                } else {
                    relative_genomic = sub(start, int_of(reference[1]));
                    viewpoint_index = py::reinterpret_borrow<py::object>(pViewpointIndexEnd);
                }
                const py::object relative_position = add(viewpoint_index, truediv(relative_genomic, pResolution));
                const py::object stop = add(relative_position, width);
                py::list area;
                area.append(relative_position);
                area.append(stop);
                highlight_areas_list.append(area);
                py::list p;
                p.append(int_of(relative_position));
                p.append(int_of(stop));
                p.append(float_of(getitem(line, -3)));
                p_values.append(p);
            }
        }
        if (highlight_areas_list.empty()) {
            return py::make_tuple(py::none(), py::none());
        }
        return py::make_tuple(highlight_areas_list, p_values);
    }, py::arg("self"), py::arg("pFilePath"), py::arg("pInternalIdentifierTriplet"), py::arg("pViewpointIndexStart"),
       py::arg("pViewpointIndexEnd"), py::arg("pResolution"), py::arg("pRange"), py::arg("pViewpoint"));

    m.def("Viewpoint_readTargetHDFFile", [](py::object self, py::handle pFile) {
        (void)self;
        py::dict present_genes;
        py::dict target_dict;
        hc::H5Reader reader = open_h5(pFile);
        const std::vector<std::string> keys_target_file = reader.keys("/");
        py::list keys_list;
        for (const std::string& key : keys_target_file) {
            keys_list.append(to_str(key));
        }
        debug("keys_interactionFile " + fmt(keys_list));
        for (const std::string& outer : keys_target_file) {
            const py::str outer_name = to_str(outer);
            if (!contains(present_genes, outer_name)) {
                present_genes[outer_name] = py::dict();
            }
            if (reader.kind(outer) == hc::H5Reader::Kind::Dataset) {
                raise_error("AttributeError", "'Dataset' object has no attribute 'keys'");
            }
            for (const std::string& inner : reader.keys(outer)) {
                const py::str inner_name = to_str(inner);
                const py::object genes_of_outer = getitem(present_genes, outer_name);
                if (!contains(genes_of_outer, inner_name)) {
                    setitem(genes_of_outer, inner_name, py::list());
                }
                const std::string inner_path = outer + "/" + inner;
                if (reader.kind(inner_path) == hc::H5Reader::Kind::Dataset) {
                    raise_error("TypeError", "Accessing a group is done with bytes or str, not <class 'str'>");
                }
                const std::string genes_path = inner_path + "/genes";
                if (!reader.contains(genes_path)) {
                    raise_error("KeyError", reader.missing_message(genes_path));
                }
                for (const std::string& gene : reader.keys(genes_path)) {
                    const py::str gene_name = to_str(gene);
                    py::list path;
                    path.append(outer_name);
                    path.append(inner_name);
                    path.append(py::str("genes"));
                    path.append(gene_name);
                    target_dict[gene_name] = path;
                    getitem(genes_of_outer, inner_name).attr("append")(gene_name);
                }
            }
        }
        reader.close();
        return py::make_tuple(target_dict, present_genes);
    }, py::arg("self"), py::arg("pFile"));

    m.def("Viewpoint_readAggregatedFileHDF", [](py::object self, py::handle pAggregatedFileName, py::handle pInternalPath) {
        (void)self;
        hc::H5Reader reader = open_h5(pAggregatedFileName);
        const std::string internal_path = utf8(join("/", pInternalPath));
        // h5py's File.get: None for a missing object, and None[()] raises
        const auto get_scalar = [&](const std::string& name) {
            const std::string path = internal_path + "/" + name;
            if (!reader.contains(path)) {
                raise_error("TypeError", "'NoneType' object is not subscriptable");
            }
            return h5_scalar(reader.read(path));
        };
        const py::object chromosome = decode_utf8(get_scalar("chromosome"));
        const py::object gene_name = decode_utf8(get_scalar("gene_name"));
        const hc::H5Value start_list = reader.read(internal_path + "/start_list");
        const hc::H5Value end_list = reader.read(internal_path + "/end_list");
        const hc::H5Value relative_distance_list = reader.read(internal_path + "/relative_distance_list");
        const hc::H5Value raw_target_list = reader.read(internal_path + "/raw_target_list");
        const py::object sum_of_interactions = float_of(get_scalar("sum_of_interactions"));
        reader.close();
        py::list line_content;
        py::list data;
        for (std::size_t i = 0; i < start_list.length; ++i) {
            const auto element = [&](const hc::H5Value& value) {
                if (i >= value.length) {
                    raise_error("IndexError", "index " + std::to_string(i) + " is out of bounds for axis 0 with size " +
                                                  std::to_string(value.length));
                }
                return h5_element(value, i);
            };
            py::list line;
            line.append(chromosome);
            line.append(element(start_list));
            line.append(element(end_list));
            line.append(gene_name);
            line.append(sum_of_interactions);
            line.append(element(relative_distance_list));
            line.append(element(raw_target_list));
            line_content.append(line);
            py::list datum;
            datum.append(sum_of_interactions);
            datum.append(element(raw_target_list));
            data.append(datum);
        }
        return py::make_tuple(line_content, data);
    }, py::arg("self"), py::arg("pAggregatedFileName"), py::arg("pInternalPath"));

    m.def("Viewpoint_readDifferentialFile", [](py::object self, py::handle pFilePath, py::handle pQuadruple) {
        (void)self;
        hc::H5Reader reader = open_h5(pFilePath);
        const std::string internal_path = utf8(join("/", pQuadruple));
        py::list output_list;
        for (const char* item : {"accepted", "all", "rejected"}) {
            try {
                const std::string base = internal_path + "/" + item;
                if (!reader.contains(base)) {
                    raise_error("KeyError", reader.missing_message(base));
                }
                const auto get_scalar = [&](const std::string& name) {
                    const std::string path = base + "/" + name;
                    if (!reader.contains(path)) {
                        raise_error("TypeError", "'NoneType' object is not subscriptable");
                    }
                    return h5_scalar(reader.read(path));
                };
                const py::object chromosome = decode_utf8(get_scalar("chromosome"));
                const py::object start_list = h5_array(reader.read(base + "/start_list"));
                const py::object end_list = h5_array(reader.read(base + "/end_list"));
                const py::object relative_distance_list = h5_array(reader.read(base + "/relative_distance_list"));
                const py::object gene = decode_utf8(get_scalar("gene"));
                const py::object pvalue_list = h5_array(reader.read(base + "/pvalue_list"));
                const py::object raw_target_list_1 = h5_array(reader.read(base + "/raw_target_list_1"));
                const py::object raw_target_list_2 = h5_array(reader.read(base + "/raw_target_list_2"));
                const py::object sum_1 = float_of(get_scalar("sum_of_interactions_1"));
                const py::object sum_2 = float_of(get_scalar("sum_of_interactions_2"));
                const py::ssize_t n = length(start_list);
                const auto repeated = [n](const py::object& value) {
                    py::list out;
                    for (py::ssize_t k = 0; k < n; ++k) {
                        out.append(value);
                    }
                    return out;
                };
                output_list.append(builtin("zip")(repeated(chromosome), start_list, end_list, repeated(gene),
                                                  relative_distance_list, repeated(sum_1), raw_target_list_1,
                                                  repeated(sum_2), raw_target_list_2, pvalue_list));
            } catch (const hicx::PythonError& error) {
                debug("exception " + exception_text(error));
                output_list.append(py::list());
            } catch (py::error_already_set& error) {
                if (!is_exception(error)) {
                    throw;
                }
                debug("exception " + exception_text(error));
                output_list.append(py::list());
            }
        }
        return output_list;
    }, py::arg("self"), py::arg("pFilePath"), py::arg("pQuadruple"));
}

}  // namespace hicx::python
