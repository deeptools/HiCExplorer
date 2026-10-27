// Port of hicexplorer/hicPlotAverageRegions.py (cpp/PLAN.md tier 7, option (a)).
//
// The data of this figure is the npz that hicAverageRegions writes, which is
// ported in C++ (cpp/tools/hicAverageRegions.cpp) and compared with the
// Python output by the npz comparator. This binary parses the command line,
// reads the archive as scipy.sparse.load_npz does (csr, csc, coo, dia or bsr,
// in its stored dtype), and computes toarray() and np.triu(); the dense upper
// triangle goes to the drawing module plot/hicexplorer_plot/
// hicPlotAverageRegions.py as an .npy file in the same dtype. The module does
// what the reference's main does from there: scipy.ndimage.rotate by 45
// degrees with cval NaN, the upper half, +1 for --log1p, and matshow with the
// chosen norm and a colorbar. The rotation is a spline interpolation of the
// image for display, so it stays with the scipy of the drawing layer, and the
// +1 follows it there because it is added to the rotated values.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/npz_file.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kDescription =
    "\n"
    "        hicPlotAverage regions plots the data computed by hicAverageRegions. It shows the "
    "summed up and averaged regions around\n"
    "        all given reference points. This tool is useful to plot differences at certain "
    "reference points as for example TAD boundaries between samples.\n";

const char* const kUsage =
    "usage: hicPlotAverageRegions --matrix MATRIX --outputFile OUTPUTFILE [--log1p]\n"
    "                             [--log] [--colorMap COLORMAP] [--vMin VMIN]\n"
    "                             [--vMax VMAX] [--dpi DPI] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "        hicPlotAverage regions plots the data computed by hicAverageRegions. It shows the "
    "summed up and averaged regions around\n"
    "        all given reference points. This tool is useful to plot differences at certain "
    "reference points as for example TAD boundaries between samples.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The averaged regions file computed by\n"
    "                        hicAverageRegions (npz file).\n"
    "  --outputFile OUTPUTFILE, -o OUTPUTFILE\n"
    "                        The averaged regions plot.\n"
    "\n"
    "Optional arguments:\n"
    "  --log1p               Plot log1p of the matrix values.\n"
    "  --log                 Plot log of the matrix values.\n"
    "  --colorMap COLORMAP   Color map to use for the heatmap. Available values can\n"
    "                        be seen here: http://matplotlib.org/examples/color/col\n"
    "                        ormaps_reference.html (Default: hot_r).\n"
    "  --vMin VMIN           Minimum score value.\n"
    "  --vMax VMAX           Maximum score value.\n"
    "  --dpi DPI             Resolution of image ifouput is a raster graphics image\n"
    "                        (e.g png, jpg) (Default: 300).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the figure is drawn by the hicexplorer_plot drawing layer with the\n"
    "calls of the Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only\n"
    "option --plotData FILE writes the data of the figure as JSON to FILE instead\n"
    "of drawing it.\n";


// ---------------------------------------------------------------------------
// scipy.sparse.load_npz(path).toarray(), then np.triu.

const hicx::npz::Array* find_array(const std::vector<hicx::npz::Array>& arrays,
                                   const char* name) {
    for (const hicx::npz::Array& array : arrays) {
        if (array.name == name) {
            return &array;
        }
    }
    return nullptr;
}

const hicx::npz::Array& require_array(const std::vector<hicx::npz::Array>& arrays,
                                      const char* name, const std::string& path) {
    const hicx::npz::Array* array = find_array(arrays, name);
    if (array == nullptr) {
        throw std::runtime_error(std::string("KeyError: '") + name +
                                 ".npy is not a file in the archive' (" + path + ")");
    }
    return *array;
}

std::int64_t element_count(const hicx::npz::Array& array) {
    std::int64_t count = 1;
    for (const std::int64_t extent : array.shape) {
        count *= extent;
    }
    return count;
}

// An integer array (indices, offsets, shape) as int64 values.
std::vector<std::int64_t> integers(const hicx::npz::Array& array) {
    const std::string& dtype = array.dtype;
    const std::int64_t count = element_count(array);
    std::vector<std::int64_t> out(static_cast<std::size_t>(count));
    auto load = [&](auto tag) {
        using T = decltype(tag);
        if (array.data.size() < static_cast<std::size_t>(count) * sizeof(T)) {
            throw std::runtime_error("ValueError: truncated array " + array.name);
        }
        for (std::int64_t i = 0; i < count; ++i) {
            T value;
            std::memcpy(&value, array.data.data() + i * static_cast<std::int64_t>(sizeof(T)),
                        sizeof(T));
            out[static_cast<std::size_t>(i)] = static_cast<std::int64_t>(value);
        }
    };
    if (dtype == "<i8") {
        load(std::int64_t{});
    } else if (dtype == "<i4") {
        load(std::int32_t{});
    } else if (dtype == "<i2") {
        load(std::int16_t{});
    } else if (dtype == "|i1" || dtype == "<i1") {
        load(std::int8_t{});
    } else if (dtype == "<u8") {
        load(std::uint64_t{});
    } else if (dtype == "<u4") {
        load(std::uint32_t{});
    } else {
        throw std::runtime_error("TypeError: the index array " + array.name +
                                 " has the unsupported dtype " + dtype);
    }
    return out;
}

// The dense matrix in the dtype of `data`: every stored element added into
// its cell in storage order (bool cells, T = Bool, are or-ed, as scipy's bool
// wrapper adds), then everything below the diagonal set to zero.
struct Bool {
    std::uint8_t value = 0;
};

template <typename T>
std::string dense_upper(const hicx::npz::Array& data, std::int64_t rows, std::int64_t cols,
                        const std::vector<std::int64_t>& row_of,
                        const std::vector<std::int64_t>& col_of,
                        const std::vector<std::int64_t>& element_of) {
    std::vector<T> dense(static_cast<std::size_t>(rows * cols), T{});
    const std::int64_t count = element_count(data);
    if (data.data.size() < static_cast<std::size_t>(count) * sizeof(T)) {
        throw std::runtime_error("ValueError: truncated array data");
    }
    for (std::size_t k = 0; k < row_of.size(); ++k) {
        const std::int64_t r = row_of[k];
        const std::int64_t c = col_of[k];
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::runtime_error("ValueError: a stored index lies outside the shape");
        }
        T value;
        std::memcpy(&value, data.data.data() + element_of[k] * static_cast<std::int64_t>(sizeof(T)),
                    sizeof(T));
        T& cell = dense[static_cast<std::size_t>(r * cols + c)];
        if constexpr (std::is_same_v<T, Bool>) {
            cell.value = (cell.value != 0 || value.value != 0) ? 1 : 0;
        } else {
            cell = static_cast<T>(cell + value);
        }
    }
    for (std::int64_t r = 1; r < rows; ++r) {
        for (std::int64_t c = 0; c < std::min(r, cols); ++c) {
            dense[static_cast<std::size_t>(r * cols + c)] = T{};
        }
    }
    std::string bytes(dense.size() * sizeof(T), '\0');
    if (!dense.empty()) {
        std::memcpy(bytes.data(), dense.data(), bytes.size());
    }
    return bytes;
}

struct DenseMatrix {
    std::string dtype;
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::string bytes;
};

DenseMatrix load_npz_upper_triangle(const std::string& path) {
    const std::vector<hicx::npz::Array> arrays = hicx::npz::read_npz(path);
    const hicx::npz::Array* format_array = find_array(arrays, "format");
    if (format_array == nullptr) {
        throw std::runtime_error("ValueError: The file " + path +
                                 " does not contain a sparse array or matrix.");
    }
    std::string format = format_array->data;
    while (!format.empty() && format.back() == '\0') {
        format.pop_back();
    }
    if (format != "csr" && format != "csc" && format != "bsr" && format != "dia" &&
        format != "coo") {
        throw std::runtime_error("NotImplementedError: Load is not implemented for sparse "
                                 "matrix of format " + format + ".");
    }
    const std::vector<std::int64_t> shape = integers(require_array(arrays, "shape", path));
    if (shape.size() != 2) {
        throw std::runtime_error("ValueError: the sparse matrix shape is not two dimensional");
    }
    DenseMatrix dense;
    dense.rows = shape[0];
    dense.cols = shape[1];
    const hicx::npz::Array& data = require_array(arrays, "data", path);
    dense.dtype = data.dtype;

    // Every stored element as (row, column, position in data), in the order
    // toarray() adds them.
    std::vector<std::int64_t> row_of;
    std::vector<std::int64_t> col_of;
    std::vector<std::int64_t> element_of;
    auto add = [&](std::int64_t r, std::int64_t c, std::int64_t element) {
        row_of.push_back(r);
        col_of.push_back(c);
        element_of.push_back(element);
    };
    if (format == "csr" || format == "csc") {
        const std::vector<std::int64_t> indptr = integers(require_array(arrays, "indptr", path));
        const std::vector<std::int64_t> indices = integers(require_array(arrays, "indices", path));
        for (std::size_t major = 0; major + 1 < indptr.size(); ++major) {
            for (std::int64_t k = indptr[major]; k < indptr[major + 1]; ++k) {
                const std::int64_t minor = indices.at(static_cast<std::size_t>(k));
                if (format == "csr") {
                    add(static_cast<std::int64_t>(major), minor, k);
                } else {
                    add(minor, static_cast<std::int64_t>(major), k);
                }
            }
        }
    } else if (format == "coo") {
        const std::vector<std::int64_t> row = integers(require_array(arrays, "row", path));
        const std::vector<std::int64_t> col = integers(require_array(arrays, "col", path));
        for (std::size_t k = 0; k < row.size(); ++k) {
            add(row[k], col.at(k), static_cast<std::int64_t>(k));
        }
    } else if (format == "dia") {
        // dia's tocoo: diagonal d holds A[n - offsets[d], n] at data[d, n].
        const std::vector<std::int64_t> offsets = integers(require_array(arrays, "offsets", path));
        const std::int64_t length = data.shape.size() == 2 ? data.shape[1] : 0;
        for (std::size_t d = 0; d < offsets.size(); ++d) {
            for (std::int64_t n = 0; n < length; ++n) {
                const std::int64_t r = n - offsets[d];
                if (r >= 0 && r < dense.rows && n < dense.cols) {
                    add(r, n, static_cast<std::int64_t>(d) * length + n);
                }
            }
        }
    } else {  // bsr
        const std::vector<std::int64_t> indptr = integers(require_array(arrays, "indptr", path));
        const std::vector<std::int64_t> indices = integers(require_array(arrays, "indices", path));
        if (data.shape.size() != 3) {
            throw std::runtime_error("ValueError: bsr data is not three dimensional");
        }
        const std::int64_t block_rows = data.shape[1];
        const std::int64_t block_cols = data.shape[2];
        for (std::size_t block_row = 0; block_row + 1 < indptr.size(); ++block_row) {
            for (std::int64_t k = indptr[block_row]; k < indptr[block_row + 1]; ++k) {
                const std::int64_t block_col = indices.at(static_cast<std::size_t>(k));
                for (std::int64_t r = 0; r < block_rows; ++r) {
                    for (std::int64_t c = 0; c < block_cols; ++c) {
                        add(static_cast<std::int64_t>(block_row) * block_rows + r,
                            block_col * block_cols + c, (k * block_rows + r) * block_cols + c);
                    }
                }
            }
        }
    }

    const std::string& t = dense.dtype;
    if (t == "<f8") {
        dense.bytes = dense_upper<double>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<f4") {
        dense.bytes = dense_upper<float>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<i8") {
        dense.bytes = dense_upper<std::int64_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<i4") {
        dense.bytes = dense_upper<std::int32_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<i2") {
        dense.bytes = dense_upper<std::int16_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "|i1") {
        dense.bytes = dense_upper<std::int8_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<u8") {
        dense.bytes = dense_upper<std::uint64_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<u4") {
        dense.bytes = dense_upper<std::uint32_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "<u2") {
        dense.bytes = dense_upper<std::uint16_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "|u1") {
        dense.bytes = dense_upper<std::uint8_t>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else if (t == "|b1") {
        dense.bytes = dense_upper<Bool>(data, dense.rows, dense.cols, row_of, col_of, element_of);
    } else {
        throw std::runtime_error("TypeError: the matrix values have the dtype " + t +
                                 ", which this port does not read (real or integer values only)");
    }
    return dense;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("hicPlotAverageRegions", kDescription);
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"npz"})
        .help("The averaged regions file computed by hicAverageRegions (npz file).");
    required.add({"--outputFile", "-o"})
        .required()
        .output({"png", "pdf", "svg"})
        .help("The averaged regions plot.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--log1p"}).action(cli::Action::StoreTrue).help("Plot log1p of the matrix values.");
    optional.add({"--log"}).action(cli::Action::StoreTrue).help("Plot log of the matrix values.");
    optional.add({"--colorMap"})
        .default_value("hot_r")
        .help("Color map to use for the heatmap. Available values can be seen here: "
              "http://matplotlib.org/examples/color/colormaps_reference.html (Default: "
              "%(default)s).");
    optional.add({"--vMin"}).type("float").help("Minimum score value.");
    optional.add({"--vMax"}).type("float").help("Maximum score value.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{300})
        .help("Resolution of image ifouput is a raster graphics image (e.g png, jpg) (Default: "
              "%(default)s).");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the figure is drawn from as JSON to this file and do not draw "
              "the figure.");

    const cli::Namespace ns = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotAverageRegions", !ns.given("plotData")); refused != 0) {
        return refused;
    }
    const std::string matrix = ns.str("matrix");

    std::vector<std::string> temporary_files;
    std::string dense_path;
    try {
        const DenseMatrix dense = load_npz_upper_triangle(matrix);
        dense_path = plot::temporary_file();
        temporary_files.push_back(dense_path);
        plot::write_npy(dense_path, dense.dtype, {dense.rows, dense.cols}, dense.bytes);
    } catch (const std::exception& error) {
        for (const std::string& path : temporary_files) {
            std::remove(path.c_str());
        }
        std::fprintf(stderr, "hicPlotAverageRegions: %s\n", error.what());
        return 1;
    }

    auto optional_number = [&](const std::string& dest) {
        const std::optional<double> value = ns.opt_real(dest);
        return value.has_value() ? plot::json_number(*value) : std::string("null");
    };
    plot::JsonObject data;
    data.add("matrix", plot::json_string(dense_path));
    data.add("outputFile", plot::json_string(ns.str("outputFile")));
    data.add("log1p", plot::json_bool(ns.flag("log1p")));
    data.add("log", plot::json_bool(ns.flag("log")));
    data.add("colorMap", plot::json_string(ns.str("colorMap")));
    data.add("vMin", optional_number("vMin"));
    data.add("vMax", optional_number("vMax"));
    data.add("dpi", plot::json_int(ns.integer("dpi")));
    data.add("temporary_files", plot::json_strings(temporary_files));

    hicx::report_resource_usage("hicPlotAverageRegions");
    return plot::draw("hicPlotAverageRegions", data.str(), ns.opt_str("plotData"));
}
