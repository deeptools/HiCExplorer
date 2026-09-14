// Port of hicexplorer/hicPlotViewpoint.py (cpp/PLAN.md tier 7, option (a)).
//
// The interactions of a reference point or region with every bin of a region
// are summed here, per matrix, and the --interactionOutFileName bedgraphs are
// written. The figure is drawn by plot/hicexplorer_plot/hicPlotViewpoint.py
// with the matplotlib calls of the Python tool (hicx::plot::draw).
//
// What is reproduced, pinned by the harness cases:
//
//  1. The region and the reference point lose every ',', ';' and '!', every
//     '-' becomes ':', and the text is split at ':'. The region needs three
//     fields; the reference point two (a point) or three (a region).
//  2. getRegionBinRange: the bins containing the start and the end position,
//     so the region's end bin is excluded from the sums and a reference
//     point is the single bin containing it. A position outside every bin
//     ends the reference with a TypeError; a chromosome that is not in the
//     matrix with a ValueError.
//  3. --chromosome reaches keepOnlyTheseChr as a string, and a string is
//     iterated by character: every character must be a chromosome name, and
//     the chromosomes kept are those whose name occurs in the string. So `-C X`
//     keeps X and `-C chr2` fails.
//  4. The sums add the full symmetric matrix entries in float64, viewpoint
//     bin by viewpoint bin; the bedgraph values print with '{:.12f}'. One
//     bedgraph per matrix, <-i>_<basename>.bedgraph, later matrices with the
//     same basename overwriting earlier ones.
//  5. The ticks use the bin ranges of the last matrix.
//
// Deliberate deviation: a region that does not split into three fields makes
// the reference log an error and exit 0 without writing anything. The port
// exits 1, because a tool never exits 0 without writing the files the user
// asked for (cpp/AGENTS_CONTRACT.md rule 7).

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kUsage =
    "usage: hicPlotViewpoint --matrix MATRIX [MATRIX ...] --region REGION\n"
    "                        --outFileName OUTFILENAME --referencePoint\n"
    "                        REFERENCEPOINT [--chromosome CHROMOSOME]\n"
    "                        [--interactionOutFileName INTERACTIONOUTFILENAME]\n"
    "                        [--dpi DPI] [--version] [--help]\n";

const char* const kHelp =
    "\n"
    "Plots the number of interactions around a given reference point in a region.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX [MATRIX ...], -m MATRIX [MATRIX ...]\n"
    "                        Hi-C matrix to plot.\n"
    "  --region REGION       The format is chr:start-end.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name of the image to save.\n"
    "  --referencePoint REFERENCEPOINT, -rp REFERENCEPOINT\n"
    "                        Reference point. Needs to be in the format: 'chr:100'\n"
    "                        for a single reference point or 'chr:100-200' for a\n"
    "                        reference region.\n"
    "\n"
    "Optional arguments:\n"
    "  --chromosome CHROMOSOME, -C CHROMOSOME\n"
    "                        Optional parameter: Only show results for this\n"
    "                        chromosome.\n"
    "  --interactionOutFileName INTERACTIONOUTFILENAME, -i INTERACTIONOUTFILENAME\n"
    "                        Optional parameter: If set, a bedgraph file with all\n"
    "                        interaction will be created.\n"
    "  --dpi DPI             Optional parameter: Resolution for the image in case\n"
    "                        theouput is a raster graphics image (e.g png, jpg)\n"
    "                        (Default: 300).\n"
    "  --version             show program's version number and exit\n"
    "  --help, -h            show this help message and exit\n"
    "\n"
    "C++ port: the interactions are summed in C++ and the figure is drawn by the\n"
    "hicexplorer_plot drawing layer with the matplotlib calls of the Python tool\n"
    "(HICX_PLOT_PYTHON names the interpreter). The C++-only option --plotData FILE\n"
    "writes the data of the figure as JSON to FILE instead of drawing it.\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string region;
    std::string out_file_name;
    std::string reference_point;
    std::optional<std::string> chromosome;
    std::optional<std::string> interaction_out_file_name;
    std::int64_t dpi = 300;
    std::optional<std::string> plot_data;
};

Arguments parse_arguments(int argc, char** argv) {
    cli::Parser parser("hicPlotViewpoint",
                       "Plots the number of interactions around a given reference point in a "
                       "region.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .nargs("+")
        .input({"h5", "cool"})
        .help("Hi-C matrix to plot.");
    required.add({"--region"}).required().help("The format is chr:start-end. ");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name of the image to save.");
    required.add({"--referencePoint", "-rp"})
        .required()
        .help("Reference point. Needs to be in the format: 'chr:100' for a single reference point "
              "or 'chr:100-200' for a reference region.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--chromosome", "-C"})
        .help("Optional parameter: Only show results for this chromosome.");
    optional.add({"--interactionOutFileName", "-i"})
        .required(false)
        .output({"bedgraph"}, "prefix")
        .help("Optional parameter:  If set, a bedgraph file with all interaction will be "
              "created.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(std::int64_t{300})
        .help("Optional parameter: Resolution for the image in case theouput is a raster "
              "graphics image (e.g png, jpg) (Default: %(default)s).");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the figure is drawn from as JSON to this file and do not draw "
              "the figure.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("matrix");
    args.region = ns.str("region");
    args.out_file_name = ns.str("outFileName");
    args.reference_point = ns.str("referencePoint");
    args.chromosome = ns.opt_str("chromosome");
    args.interaction_out_file_name = ns.opt_str("interactionOutFileName");
    args.dpi = ns.integer("dpi");
    args.plot_data = ns.opt_str("plotData");
    return args;
}

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

std::string normalised(std::string text) {
    std::string out;
    for (const char c : text) {
        if (c == ',' || c == ';' || c == '!') {
            continue;
        }
        out += c == '-' ? ':' : c;
    }
    return out;
}

std::vector<std::string> split_colon(const std::string& text) {
    std::vector<std::string> parts;
    std::string current;
    for (const char c : text) {
        if (c == ':') {
            parts.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    parts.push_back(current);
    return parts;
}

std::int64_t to_int(const std::string& text) {
    std::int64_t value = 0;
    if (!cli::python_int(text, &value)) {
        throw PythonError("ValueError: invalid literal for int() with base 10: '" + text + "'");
    }
    return value;
}

// The bins of hiCMatrix after keepOnlyTheseChr, by their original index.
struct View {
    std::vector<std::int64_t> original;  // kept bin -> matrix bin
    std::vector<hicx::CutInterval> intervals;
};

View keep_chromosomes(const hicx::ToolMatrix& matrix, const std::optional<std::string>& chromosome) {
    View view;
    const std::vector<hicx::CutInterval>& all = matrix.cut_intervals();
    if (!chromosome.has_value()) {
        view.intervals = all;
        view.original.resize(all.size());
        for (std::size_t i = 0; i < all.size(); ++i) {
            view.original[i] = static_cast<std::int64_t>(i);
        }
        return view;
    }
    const auto& boundaries = matrix.boundaries();
    auto known = [&](const std::string& name) {
        return std::any_of(boundaries.begin(), boundaries.end(),
                           [&](const auto& entry) { return entry.first == name; });
    };
    for (const char c : *chromosome) {
        if (!known(std::string(1, c))) {
            throw PythonError("ValueError: Chromosome name not in matrix. '" + std::string(1, c) +
                              "' (--chromosome is iterated by character, hicPlotViewpoint.py:77)");
        }
    }
    std::vector<char> selected(all.size(), 0);
    for (const auto& entry : boundaries) {
        if (chromosome->find(entry.first) == std::string::npos) {
            continue;
        }
        for (std::int64_t i = entry.second.first; i < entry.second.last; ++i) {
            selected[static_cast<std::size_t>(i)] = 1;
        }
    }
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (selected[i]) {
            view.original.push_back(static_cast<std::int64_t>(i));
            view.intervals.push_back(all[i]);
        }
    }
    if (view.intervals.empty()) {
        throw PythonError("StopIteration: no chromosome is left after --chromosome");
    }
    return view;
}

// sorted(interval_trees[chrom][position:position + 1])[0].data
std::optional<std::int64_t> bin_containing(const View& view, const std::string& chrom,
                                           std::int64_t position) {
    std::optional<std::tuple<std::int64_t, std::int64_t, std::int64_t>> best;
    for (std::size_t i = 0; i < view.intervals.size(); ++i) {
        const hicx::CutInterval& bin = view.intervals[i];
        if (bin.chrom == chrom && bin.start <= position && position < bin.end) {
            const auto key = std::make_tuple(bin.start, bin.end, static_cast<std::int64_t>(i));
            if (!best.has_value() || key < *best) {
                best = key;
            }
        }
    }
    if (!best.has_value()) {
        return std::nullopt;
    }
    return std::get<2>(*best);
}

std::pair<std::int64_t, std::int64_t> region_bin_range(const View& view, const std::string& chrom,
                                                       std::int64_t start, std::int64_t end) {
    if (std::none_of(view.intervals.begin(), view.intervals.end(),
                     [&](const hicx::CutInterval& bin) { return bin.chrom == chrom; })) {
        throw PythonError("ValueError: chromosome: " + chrom + " name not found in matrix");
    }
    const std::optional<std::int64_t> first = bin_containing(view, chrom, start);
    const std::optional<std::int64_t> last = bin_containing(view, chrom, end);
    if (!first.has_value() || !last.has_value()) {
        throw PythonError("TypeError: cannot unpack non-iterable NoneType object (no bin of " +
                          chrom + " contains " + std::to_string(first ? end : start) + ")");
    }
    return {*first, *last};
}

struct Viewpoint {
    std::int64_t view_point_start = 0;
    std::int64_t view_point_end = 0;
    std::pair<std::int64_t, std::int64_t> view_point_range;
    std::vector<double> data;
    std::string bedgraph;
};

Viewpoint viewpoint_values(const std::string& path, const std::vector<std::string>& reference,
                           const std::string& chrom, std::int64_t region_start,
                           std::int64_t region_end, const Arguments& args) {
    const hicx::ToolMatrix matrix = hicx::ToolMatrix::load(path);
    const View view = keep_chromosomes(matrix, args.chromosome);

    Viewpoint out;
    if (reference.size() == 2) {
        std::tie(out.view_point_start, out.view_point_end) =
            region_bin_range(view, reference[0], to_int(reference[1]), to_int(reference[1]));
    } else if (reference.size() == 3) {
        std::tie(out.view_point_start, out.view_point_end) =
            region_bin_range(view, reference[0], to_int(reference[1]), to_int(reference[2]));
    } else {
        std::string text = "[";
        for (std::size_t i = 0; i < reference.size(); ++i) {
            text += (i > 0 ? ", '" : "'") + reference[i] + "'";
        }
        std::fprintf(stderr, "No valid reference point given. %s]\n", text.c_str());
        std::exit(1);
    }
    out.view_point_range = region_bin_range(view, chrom, region_start, region_end);
    const std::int64_t elements = out.view_point_range.second - out.view_point_range.first;
    if (elements < 0) {
        throw PythonError("ValueError: negative dimensions are not allowed (the region ends before "
                          "it starts)");
    }
    out.data.assign(static_cast<std::size_t>(elements), 0.0);
    const hicx::CsrMatrix& csr = matrix.matrix();
    const bool interactions = args.interaction_out_file_name.has_value();
    for (std::int64_t vp = out.view_point_start; vp <= out.view_point_end; ++vp) {
        if (vp >= static_cast<std::int64_t>(view.intervals.size())) {
            throw PythonError("ValueError: binIndex: " + std::to_string(vp) + " not found");
        }
        const hicx::CutInterval& first = view.intervals[static_cast<std::size_t>(vp)];
        const std::int64_t row = view.original[static_cast<std::size_t>(vp)];
        for (std::int64_t j = 0; j < elements; ++j) {
            const std::int64_t idx = out.view_point_range.first + j;
            const double value =
                csr.at(row, view.original[static_cast<std::size_t>(idx)]);
            out.data[static_cast<std::size_t>(j)] += value;
            if (interactions) {
                const hicx::CutInterval& second = view.intervals[static_cast<std::size_t>(idx)];
                char number[64];
                std::snprintf(number, sizeof number, "%.12f", value);
                out.bedgraph += first.chrom + "\t" + std::to_string(first.start) + "\t" +
                                std::to_string(first.end) + "\t" + second.chrom + "\t" +
                                std::to_string(second.start) + "\t" +
                                std::to_string(second.end) + "\t" + number + "\n";
            }
        }
    }
    return out;
}

std::string basename_of(const std::string& path) {
    const std::string::size_type slash = path.rfind('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotViewpoint", !args.plot_data.has_value()); refused != 0) {
        return refused;
    }
    std::fputs("This tool is deprecated. Please use chicViewpoint, chicViewpointBackgroundModel and "
               "chicPlotViewpoint.\n",
               stderr);

    std::vector<Viewpoint> viewpoints;
    std::vector<std::string> reference;
    std::int64_t region_start = 0;
    std::int64_t region_end = 0;
    try {
        if (args.region.empty()) {
            throw PythonError("NameError: name 'chrom' is not defined (an empty --region)");
        }
        const std::string region_text = normalised(args.region);
        const std::vector<std::string> region = split_colon(region_text);
        if (region.size() != 3) {
            std::fprintf(stderr,
                         "Region format is invalid %s. The reference exits 0 here without "
                         "writing anything; the C++ port exits 1 (cpp/AGENTS_CONTRACT.md rule "
                         "7).\n",
                         region_text.c_str());
            return 1;
        }
        const std::string chrom = region[0];
        region_start = to_int(region[1]);
        region_end = to_int(region[2]);
        reference = split_colon(normalised(args.reference_point));
        for (const std::string& path : args.matrices) {
            viewpoints.push_back(
                viewpoint_values(path, reference, chrom, region_start, region_end, args));
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotViewpoint: %s\n", error.what());
        return 1;
    }

    std::vector<std::string> data_json;
    std::vector<std::string> legend;
    for (std::size_t i = 0; i < viewpoints.size(); ++i) {
        data_json.push_back(plot::json_numbers(viewpoints[i].data));
        legend.push_back(basename_of(args.matrices[i]));
    }
    if (args.interaction_out_file_name.has_value()) {
        for (std::size_t i = 0; i < viewpoints.size(); ++i) {
            const std::string name =
                *args.interaction_out_file_name + "_" + legend[i] + ".bedgraph";
            std::ofstream out(name, std::ios::binary | std::ios::trunc);
            out << viewpoints[i].bedgraph;
            out.close();
            if (!out) {
                std::fprintf(stderr, "hicPlotViewpoint: cannot write %s\n", name.c_str());
                return 1;
            }
        }
    }

    const Viewpoint& last = viewpoints.back();
    plot::JsonObject data;
    data.add("outFileName", plot::json_string(args.out_file_name));
    data.add("dpi", plot::json_int(args.dpi));
    data.add("referencePoint", plot::json_strings(reference));
    data.add("region_start", plot::json_int(region_start));
    data.add("region_end", plot::json_int(region_end));
    data.add("view_point_start", plot::json_int(last.view_point_start));
    data.add("view_point_end", plot::json_int(last.view_point_end));
    data.add("view_point_range", plot::json_ints({last.view_point_range.first,
                                                  last.view_point_range.second}));
    data.add("data", plot::json_list(data_json));
    data.add("legend", plot::json_strings(legend));

    hicx::report_resource_usage("hicPlotViewpoint");
    return plot::draw("hicPlotViewpoint", data.str(), args.plot_data);
}
