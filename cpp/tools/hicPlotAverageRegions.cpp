// Port of hicexplorer/hicPlotAverageRegions.py (cpp/PLAN.md tier 7, option (a)).
//
// The data of this figure is the npz that hicAverageRegions writes, which is
// ported in C++ (cpp/tools/hicAverageRegions.cpp) and compared with the
// Python output by the npz comparator. This binary parses the command line
// and checks that the file is a scipy sparse npz archive, as load_npz does,
// and the drawing module plot/hicexplorer_plot/hicPlotAverageRegions.py does
// what the reference's main does after parsing: load_npz, toarray, np.triu,
// scipy.ndimage.rotate by 45 degrees with cval NaN, the upper half, +1 for
// --log1p, and matshow with the chosen norm and a colorbar. The rotation is a
// spline interpolation for display only, so it stays with the scipy of the
// reference environment in the drawing layer.

#include <cstdint>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
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
    const std::string matrix = ns.str("matrix");

    // scipy.sparse.load_npz: the archive must name a sparse format.
    try {
        const std::vector<hicx::npz::Array> arrays = hicx::npz::read_npz(matrix);
        std::optional<std::string> format;
        for (const hicx::npz::Array& array : arrays) {
            if (array.name == "format") {
                format = array.data;
                while (!format->empty() && format->back() == '\0') {
                    format->pop_back();
                }
            }
        }
        if (!format.has_value()) {
            throw std::runtime_error("ValueError: The file " + matrix +
                                     " does not contain a sparse array or matrix.");
        }
        if (*format != "csr" && *format != "csc" && *format != "bsr" && *format != "dia" &&
            *format != "coo") {
            throw std::runtime_error("NotImplementedError: Load is not implemented for sparse "
                                     "matrix of format " + *format + ".");
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotAverageRegions: %s\n", error.what());
        return 1;
    }

    auto optional_number = [&](const std::string& dest) {
        const std::optional<double> value = ns.opt_real(dest);
        return value.has_value() ? plot::json_number(*value) : std::string("null");
    };
    plot::JsonObject data;
    data.add("matrix", plot::json_string(matrix));
    data.add("outputFile", plot::json_string(ns.str("outputFile")));
    data.add("log1p", plot::json_bool(ns.flag("log1p")));
    data.add("log", plot::json_bool(ns.flag("log")));
    data.add("colorMap", plot::json_string(ns.str("colorMap")));
    data.add("vMin", optional_number("vMin"));
    data.add("vMax", optional_number("vMax"));
    data.add("dpi", plot::json_int(ns.integer("dpi")));

    hicx::report_resource_usage("hicPlotAverageRegions");
    return plot::draw("hicPlotAverageRegions", data.str(), ns.opt_str("plotData"));
}
