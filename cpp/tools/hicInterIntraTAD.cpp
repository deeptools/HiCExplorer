// Port of hicexplorer/hicInterIntraTAD.py.
//
// For every TAD of a hicFindTADs domain file, writes the contact sum, the
// number of cells, the number of stored cells and the density of the TAD
// itself and of its left and right inter-TAD blocks, and the ratios of the
// inter sums to the intra sum.
//
// The geometry, and every quirk of it, is shared with hicDifferentialTAD and
// documented in tad_contacts_impl.hpp. What is specific to this tool is the
// text, which is Python's str() of whatever object the Python happens to hold:
//
//  * A column that was not computed because the neighbour does not exist is
//    the literal int 0 (hicInterIntraTAD.py:179-186) and prints as "0", while
//    a computed density is a Python float and prints as "1.0". The same
//    column therefore switches format between rows.
//  * A sum is the numpy scalar scipy's `.sum()` returns, so its type follows
//    the matrix dtype: int64 for an integer matrix, float32 or float64
//    otherwise, and a ratio follows numpy's promotion of the two operands.
//    tad_contacts_impl's PyNumber carries that.
//  * A block with no cells makes the density a ZeroDivisionError. The worker
//    reports it and the tool exits with status 1 before writing anything. On
//    an h5 matrix this happens for the last TAD of every chromosome with at
//    least three TADs (point 3 of tad_contacts_impl.hpp), so the Python tool
//    cannot process an h5 matrix at all; pinned by
//    test_hicInterIntraTAD.py::test_h5_input_raises_zero_division and
//    reproduced here.
//
// The scatter plot (cpp/PLAN.md tier 7, option (a)): the Python writes it on
// every run, to ratio.png when --outFileNameRatioPlot is not given
// (hicInterIntraTAD.py:39-42, :509-514). After the table the two ratio
// columns go to plot/hicexplorer_plot/hicInterIntraTAD.py, which draws them
// with the reference's matplotlib calls. The C++-only option --plotData
// writes that data as JSON instead of drawing.
//
// Threading: one independent problem per TAD, hicx::parallel_for into
// preallocated slots, written in file order, so the output does not depend on
// --threads. The Python's does in the degenerate cases listed in
// tad_contacts_impl.hpp; this port reproduces --threads 1.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/argparse.hpp"
#include "hicx/parallel.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/version.hpp"
#include "tad_contacts_impl.hpp"

namespace {

const char* const kUsage =
    "usage: hicInterIntraTAD [--matrix MATRIX] [--tadDomains TADDOMAINS]\n"
    "                        [--outFileName OUTFILENAME]\n"
    "                        [--outFileNameRatioPlot OUTFILENAMERATIOPLOT]\n"
    "                        [--fontsize FONTSIZE] [--dpi DPI] [--threads THREADS]\n"
    "                        [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Extracts and computes different inter and intra TAD values and ratios.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix which was used to compute the TADs\n"
    "  --tadDomains TADDOMAINS, -td TADDOMAINS\n"
    "                        The TADs domain file computed by hicFindTADs.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileNameRatioPlot OUTFILENAMERATIOPLOT, -op OUTFILENAMERATIOPLOT\n"
    "                        Outfile name for the inter-left/intra vs inter-\n"
    "                        right/intra ratio plot\n"
    "  --fontsize FONTSIZE   Fontsize in the plot for x and y axis.\n"
    "  --dpi DPI             The dpi of the scatter plot.\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads to use, the parallelization is\n"
    "                        implemented per chromosome (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the table and the ratios are computed in C++, and the ratio plot is\n"
    "drawn by the hicexplorer_plot drawing layer with the matplotlib calls of the\n"
    "Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option\n"
    "--plotData FILE writes the data of the plot as JSON to FILE instead of drawing\n"
    "it.\n";

struct Arguments {
    std::optional<std::string> matrix;
    std::optional<std::string> domains;
    std::string out_file = "output_interintra_tad.tzt";
    std::string plot_file = "ratio.png";
    double fontsize = 15;
    long long dpi = 300;
    std::optional<std::string> plot_data;
    long long threads = 4;
};

// hicInterIntraTAD.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicInterIntraTAD",
                       "Extracts and computes different inter and intra TAD values and ratios.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .input({"cool", "h5"})
        .help("The matrix which was used to compute the TADs");
    required.add({"--tadDomains", "-td"})
        .input({"bed"})
        .help("The TADs domain file computed by hicFindTADs.");
    required.add({"--outFileName", "-o"})
        .default_value("output_interintra_tad.tzt")
        .output({"txt"})
        .help("Outfile name");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileNameRatioPlot", "-op"})
        .default_value("ratio.png")
        .output({"png", "pdf", "svg"})
        .help("Outfile name for the inter-left/intra vs inter-right/intra ratio plot");
    optional.add({"--fontsize"})
        .type("float")
        .default_value(15)
        .help("Fontsize in the plot for x and y axis.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(300)
        .help("The dpi of the scatter plot.");
    optional.add({"--threads", "-t"})
        .type("int")
        .default_value(4)
        .help("Number of threads to use, the parallelization is implemented per chromosome.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the ratio plot is drawn from as JSON to this file and do not "
              "draw the plot.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.opt_str("matrix");
    args.domains = ns.opt_str("tadDomains");
    args.out_file = ns.str("outFileName");
    args.plot_file = ns.str("outFileNameRatioPlot");
    args.fontsize = ns.real("fontsize");
    args.dpi = ns.integer("dpi");
    args.plot_data = ns.opt_str("plotData");
    args.threads = ns.integer("threads");
    return args;
}

struct Side {
    hicx::tads::PyNumber sum = hicx::tads::PyNumber::py_int(0);
    hicx::tads::PyNumber density = hicx::tads::PyNumber::py_int(0);
    std::int64_t contacts = 0;
    std::int64_t nnz = 0;
};

struct TadResult {
    Side left;
    Side right;
    Side intra;
    hicx::tads::PyNumber left_ratio;
    hicx::tads::PyNumber right_ratio;
    hicx::tads::PyNumber both_ratio;
    std::optional<std::string> error;
};

// sum, shape[0] * shape[1], nnz and nnz / (shape[0] * shape[1]) of one block.
Side measure(const hicx::tads::ContactMatrix& matrix, const hicx::tads::Block& block) {
    namespace tads = hicx::tads;
    const tads::DenseBlock values = tads::extract_block(matrix, block);
    Side side;
    side.sum = tads::block_sum(values, tads::block_dtype(matrix, block));
    side.contacts = values.rows * values.cols;
    side.nnz = values.nnz;
    if (side.contacts == 0) {
        throw std::runtime_error("division by zero");
    }
    side.density = tads::PyNumber::float64(static_cast<double>(side.nnz) /
                                           static_cast<double>(side.contacts));
    return side;
}

bool write_file(const std::string& path, const std::string& content) {
    std::FILE* handle = std::fopen(path.c_str(), "wb");
    if (handle == nullptr) {
        return false;
    }
    const bool written = std::fwrite(content.data(), 1, content.size(), handle) ==
                         content.size();
    return std::fclose(handle) == 0 && written;
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    namespace tads = hicx::tads;

    std::string plot_json;
    try {
        if (!args.domains.has_value() || !args.matrix.has_value()) {
            throw std::runtime_error("--matrix and --tadDomains must both be given");
        }
        const std::vector<tads::Domain> domains = tads::read_domains(*args.domains);
        const std::vector<std::vector<tads::Domain>> chromosomes =
            tads::group_by_chromosome(domains);
        const bool is_cooler = hicx::check_cooler(*args.matrix);
        if (args.threads == 0) {
            throw std::runtime_error("integer division or modulo by zero");
        }
        const tads::ContactMatrix matrix = tads::load_contact_matrix(*args.matrix, is_cooler);

        struct Item {
            std::size_t chromosome = 0;
            std::size_t index = 0;
        };
        std::vector<Item> items;
        if (args.threads > 0) {
            for (std::size_t c = 0; c < chromosomes.size(); ++c) {
                for (std::size_t i = 0; i < chromosomes[c].size(); ++i) {
                    items.push_back(Item{c, i});
                }
            }
        }

        std::vector<TadResult> results(items.size());
        const unsigned int workers = static_cast<unsigned int>(
            std::clamp<long long>(args.threads, 1, hicx::hardware_threads()));
        hicx::parallel_for(items.size(), workers, [&](std::size_t k) {
            TadResult& result = results[k];
            try {
                const tads::TadGeometry geometry =
                    tads::tad_geometry(matrix, chromosomes[items[k].chromosome], items[k].index);
                result.intra = measure(matrix, geometry.intra);
                if (geometry.left.has_value()) {
                    result.left = measure(matrix, *geometry.left);
                }
                if (geometry.right.has_value()) {
                    result.right = measure(matrix, *geometry.right);
                }
                result.left_ratio = tads::py_divide(result.left.sum, result.intra.sum);
                result.right_ratio = tads::py_divide(result.right.sum, result.intra.sum);
                result.both_ratio = tads::py_divide(
                    tads::py_add(result.left.sum, result.right.sum), result.intra.sum);
            } catch (const std::exception& error) {
                result.error = error.what();
            }
        });
        for (const TadResult& result : results) {
            if (result.error.has_value()) {
                std::fprintf(stderr, "ERROR:hicexplorer.hicInterIntraTAD:%s\n",
                             result.error->c_str());
                return 1;
            }
        }

        std::string text = "# Created with HiCExplorer's hicInterIntraTAD version ";
        text += hicx::kVersion;
        text += "\n";
        text += "# Chromosome\tstart\tend\tname\tscore\tstrand\tinter_left_sum\t"
                "inter_right_sum\tinter_left_density\tinter_right_density\t"
                "inter_left_number_of_contacts\tinter_right_number_of_contacts\t"
                "inter_left_number_of_contacts_nnz\tinter_right_number_of_contacts_nnz\t"
                "intra_sum\tintra_number_of_contacts\tintra_number_of_contacts_nnz\t"
                "intra_density\tinter_left_intra_ratio\tinter_right_intra_ratio\t"
                "inter_left_inter_right_intra_ratio\n";
        for (std::size_t k = 0; k < items.size(); ++k) {
            const tads::Domain& domain = chromosomes[items[k].chromosome][items[k].index];
            const TadResult& r = results[k];
            std::string line;
            for (std::size_t column = 0; column < domain.text.size(); ++column) {
                if (column > 0) {
                    line += '\t';
                }
                line += domain.text[column];
            }
            const auto field = [&line](const std::string& value) {
                line += '\t';
                line += value;
            };
            field(r.left.sum.str());
            field(r.right.sum.str());
            field(r.left.density.str());
            field(r.right.density.str());
            field(std::to_string(r.left.contacts));
            field(std::to_string(r.right.contacts));
            field(std::to_string(r.left.nnz));
            field(std::to_string(r.right.nnz));
            field(r.intra.sum.str());
            field(std::to_string(r.intra.contacts));
            field(std::to_string(r.intra.nnz));
            field(r.intra.density.str());
            field(r.left_ratio.str());
            field(r.right_ratio.str());
            field(r.both_ratio.str());
            line += '\n';
            text += line;
        }
        if (!write_file(args.out_file, text)) {
            throw std::runtime_error("cannot write '" + args.out_file + "'");
        }
        // plt.scatter(inter_left_intra_ratio_list, inter_right_intra_ratio_list)
        std::vector<double> left_ratios;
        std::vector<double> right_ratios;
        left_ratios.reserve(results.size());
        right_ratios.reserve(results.size());
        for (const TadResult& r : results) {
            left_ratios.push_back(r.left_ratio.as_double());
            right_ratios.push_back(r.right_ratio.as_double());
        }
        hicx::plot::JsonObject data;
        data.add("plotFile", hicx::plot::json_string(args.plot_file));
        data.add("fontsize", hicx::plot::json_number(args.fontsize));
        data.add("dpi", hicx::plot::json_int(args.dpi));
        data.add("x", hicx::plot::json_numbers(left_ratios));
        data.add("y", hicx::plot::json_numbers(right_ratios));
        plot_json = data.str();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicInterIntraTAD: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicInterIntraTAD");
    return hicx::plot::draw("hicInterIntraTAD", plot_json, args.plot_data);
}
