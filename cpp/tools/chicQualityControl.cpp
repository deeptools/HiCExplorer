// Port of hicexplorer/chicQualityControl.py.
//
// For every reference point and every matrix: the fraction of non zero
// elements in the viewpoint over --fixateRange on either side. A reference
// point is written to --outFileName when it passes the threshold in at least
// one matrix, to *_rejected_filter when it passes in none, and to
// *_failed_reference_points when it is faulty (sparsity -1.0) in any.
//
// Behaviour reproduced as the Python has it, not as its help text says:
//
//   * --sparsity is documented as removing a viewpoint as soon as it is of bad
//     quality in at least one matrix, but chicQualityControl.py:228-240 accepts
//     it when any matrix is above the threshold. Reported as a defect.
//   * The output loop re-reads the reference point file line by line and
//     indexes the sparsity table with the line number. A line the reader
//     skipped (a blank line, or one with fewer than three fields) therefore
//     shifts every later line onto the wrong sparsity, and once the line count
//     exceeds the table the Python raises IndexError after having written the
//     earlier lines. Reproduced: the same lines are written, then exit 1.
//   * A reference point on a chromosome the matrix has, but outside its bins,
//     is a TypeError inside computeViewpoint and becomes -1.0. A reference
//     point whose start lies after its end is a ValueError there, which is not
//     caught, and aborts the run with exit 1.
//
// Figures (cpp/PLAN.md tier 7, option (a)). The Python also draws
// sparsity.png and histogram.png with matplotlib (chicQualityControl.py:262-307).
// After the report the sparsity of every reference point that is faulty in
// no matrix goes, per matrix, to plot/hicexplorer_plot/chicQualityControl.py,
// which draws both figures with the reference's calls. The C++-only option
// --plotData writes that data as JSON instead.
//
// Threading: measured, and not used. Spreading the 45 reference points of
// the test data over worker threads does not shorten the run (0.090 s of wall
// time at --threads 1, 0.093 s at 16) but raises the CPU time from 0.080 s to
// 0.377 s and the peak RSS from 45 MB to 120 MB, above this tool's memory
// budget: the work per reference point is a few thousand sparse lookups, and
// every worker thread brings its own malloc arena. Per cpp/OPTIMIZATION.md 6
// the loop is sequential. --threads is accepted, and 0 still fails as the
// Python's integer division by it does. The output is written in file order;
// the Python collects its worker chunks in index order, which is the same
// order, so neither output depends on the thread count.

#include <cstdio>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/version.hpp"

namespace {

using hicx::chic::ReferencePoint;

const char* const kUsage =
    "usage: chicQualityControl --matrices MATRICES [MATRICES ...] --referencePoints\n"
    "                          REFERENCEPOINTS --sparsity SPARSITY\n"
    "                          [--outFileName OUTFILENAME]\n"
    "                          [--outFileNameHistogram OUTFILENAMEHISTOGRAM]\n"
    "                          [--outFileNameSparsity OUTFILENAMESPARSITY]\n"
    "                          [--threads THREADS] [--fixateRange FIXATERANGE]\n"
    "                          [--dpi DPI] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Computes the sparsity of each viewpoint to determine the quality. A viewpoint is "
    "considered to be of bad quality if it is too sparse i.e. if there are too many "
    "locations with no interactions recorded.\n"
    "\n"
    "This script creates three output files: a plot with the sparsity distribution per "
    "matrix, a plot with the sparsity distribution as histograms and a filtered reference "
    "points file.\n"
    "\n"
    "An example usage is:\n"
    "\n"
    "$ chicQualityControl -m matrix1.cool matrix2.cool -rp referencePointsFile.bed --range "
    "20000 40000 --sparsity 0.01 -o referencePointFile_QC_passed.bed\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        The input matrices to apply the QC on.\n"
    "  --referencePoints REFERENCEPOINTS, -rp REFERENCEPOINTS\n"
    "                        Bed file contains all reference points which are\n"
    "                        checked for a sufficient number of interactions.\n"
    "  --sparsity SPARSITY, -s SPARSITY\n"
    "                        Viewpoints with a sparsity less than the value given\n"
    "                        are considered of bad quality. If multiple matrices\n"
    "                        are given, the viewpoint is removed as soon as it is\n"
    "                        of bad quality in at least one matrix.\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The output file name of the passed reference points.\n"
    "                        Used as prefix for the plots as well (Default:\n"
    "                        new_referencepoints.bed).\n"
    "  --outFileNameHistogram OUTFILENAMEHISTOGRAM, -oh OUTFILENAMEHISTOGRAM\n"
    "                        The output file for the histogram plot (Default:\n"
    "                        histogram.png).\n"
    "  --outFileNameSparsity OUTFILENAMESPARSITY, -os OUTFILENAMESPARSITY\n"
    "                        The output file for the sparsity distribution plot\n"
    "                        (Default: sparsity.png).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads (Default: 4).\n"
    "  --fixateRange FIXATERANGE, -fs FIXATERANGE\n"
    "                        Fixate score of background model starting at distance\n"
    "                        x. E.g. all values greater than 500kb are set to the\n"
    "                        value of the 500kb bin (Default: 500000).\n"
    "  --dpi DPI             Optional parameter: Resolution for the image if\n"
    "                        theoutput is a raster graphics image (e.g png, jpg)\n"
    "                        (Default: 300).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the sparsity and the reference point files are computed in C++, and\n"
    "the two figures are drawn by the hicexplorer_plot drawing layer with the\n"
    "matplotlib calls of the Python tool (HICX_PLOT_PYTHON names the interpreter).\n"
    "The C++-only option --plotData FILE writes the data of the figures as JSON to\n"
    "FILE instead of drawing them.\n";

std::string basename(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

struct Slot {
    double sparsity = -1.0;
    bool failed = false;
    std::string message;
};

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("chicQualityControl",
                       "Computes the sparsity of each viewpoint to determine the quality.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .nargs("+")
        .required()
        .input({"cool", "h5"})
        .help("The input matrices to apply the QC on.");
    required.add({"--referencePoints", "-rp"})
        .type("str")
        .required()
        .input({"bed"})
        .help("Bed file contains all reference points which are checked for a sufficient number "
              "of interactions.");
    required.add({"--sparsity", "-s"})
        .type("float")
        .required()
        .help("Viewpoints with a sparsity less than the value given are considered of bad "
              "quality.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("new_referencepoints.bed")
        .output({"bed"})
        .help("The output file name of the passed reference points.");
    optional.add({"--outFileNameHistogram", "-oh"})
        .default_value("histogram.png")
        .output({"png", "pdf", "svg"})
        .help("The output file for the histogram plot");
    optional.add({"--outFileNameSparsity", "-os"})
        .default_value("sparsity.png")
        .output({"png", "pdf", "svg"})
        .help("The output file for the sparsity distribution plot");
    optional.add({"--threads", "-t"}).type("int").default_value(4).help("Number of threads");
    optional.add({"--fixateRange", "-fs"})
        .type("int")
        .default_value(500000)
        .help("Range on either side of a reference point the sparsity is computed over.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(300)
        .help("Resolution for the image if the output is a raster graphics image.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figures as JSON, without drawing them (cpp/PLAN.md tier 7).")
        .help("Write the data the two figures are drawn from as JSON to this file and do not "
              "draw them.");
    const cli::Namespace args = parser.parse(argc, argv);
    if (const int refused = hicx::plot::preflight("chicQualityControl", !args.given("plotData")); refused != 0) {
        return refused;
    }
    std::string plot_json;

    const std::vector<std::string>& matrices = args.strs("matrices");
    const std::string reference_point_file = args.str("referencePoints");
    const double sparsity_threshold = args.real("sparsity");
    const std::string out_file_name = args.str("outFileName");
    const std::int64_t threads = args.integer("threads");
    const std::int64_t fixate_range = args.integer("fixateRange");

    try {
        const hicx::chic::ReferencePoints reference_points =
            hicx::chic::read_reference_points(reference_point_file);
        const std::size_t count = reference_points.points.size();
        if (threads == 0) {
            // len(referencePoints) // args.threads
            throw std::runtime_error("ZeroDivisionError: integer division or modulo by zero");
        }

        // sparsity[matrix][reference point]
        std::vector<std::vector<double>> sparsity;
        for (const std::string& matrix_path : matrices) {
            const hicx::chic::ViewpointMatrix matrix = hicx::chic::ViewpointMatrix::load(matrix_path);
            const std::vector<std::string> names = matrix.chromosome_names();
            std::vector<Slot> slots(count);
            const auto compute_slot = [&](std::size_t index) {
                Slot& slot = slots[index];
                const ReferencePoint& point = reference_points.points[index];
                try {
                    if (!matrix.has_chromosome(point.chromosome)) {
                        slot.sparsity = -1.0;
                        return;
                    }
                    const hicx::chic::ViewpointRange range = hicx::chic::calculate_viewpoint_range(
                        matrix, point, fixate_range, fixate_range);
                    try {
                        const hicx::chic::ComputedViewpoint viewpoint = hicx::chic::compute_viewpoint(
                            matrix, point, point.chromosome, range.region_start, range.region_end);
                        std::size_t nonzero = 0;
                        for (const double value : viewpoint.data) {
                            nonzero += value != 0.0 ? 1 : 0;
                        }
                        if (viewpoint.data.empty()) {
                            throw std::runtime_error("ZeroDivisionError: division by zero");
                        }
                        slot.sparsity = static_cast<double>(nonzero) /
                                        static_cast<double>(viewpoint.data.size());
                    } catch (const hicx::chic::TypeError&) {
                        slot.sparsity = -1.0;
                    } catch (const hicx::chic::IndexError&) {
                        slot.sparsity = -1.0;
                    }
                } catch (const std::exception& error) {
                    slot.failed = true;
                    slot.message = error.what();
                }
            };
            for (std::size_t index = 0; index < count; ++index) {
                compute_slot(index);
            }
            for (const Slot& slot : slots) {
                if (slot.failed) {
                    std::fprintf(stderr, "chicQualityControl: %s\n", slot.message.c_str());
                    return 1;
                }
            }
            std::vector<double> column(count);
            for (std::size_t i = 0; i < count; ++i) {
                column[i] = slots[i].sparsity;
            }
            sparsity.push_back(std::move(column));
        }

        const std::vector<std::string> lines = hicx::chic::read_lines(reference_point_file);
        std::size_t accepted = 0;
        std::size_t rejected = 0;
        std::size_t failures = 0;
        {
            std::ofstream raw(out_file_name + "_raw_filter", std::ios::binary);
            raw << "# Created with chicQualityControl version " << hicx::kVersion << "\n";
            raw << "# A sparsity of -1.0 indicates a faulty reference point e.g. no data for "
                   "this reference point was in the matrix.\n";
            raw << "# Used Matrices ";
            for (const std::string& matrix_path : matrices) {
                raw << matrix_path << "\t";
            }
            raw << "\n# Chromosome\tStart\tEnd";
            for (const std::string& matrix_path : matrices) {
                raw << "\tSparsity " << basename(matrix_path);
            }
            raw << "\n";

            std::ofstream failed(out_file_name + "_failed_reference_points", std::ios::binary);
            std::ofstream rejected_file(out_file_name + "_rejected_filter", std::ios::binary);
            std::ofstream accepted_file(out_file_name, std::ios::binary);
            for (std::size_t i = 0; i < lines.size(); ++i) {
                if (i >= count) {
                    std::fprintf(stderr,
                                 "chicQualityControl: IndexError: index %zu is out of bounds for "
                                 "axis 0 with size %zu\n",
                                 i, count);
                    return 1;
                }
                std::string sparsity_text;
                std::size_t passing = 0;
                std::size_t negative = 0;
                for (std::size_t j = 0; j < sparsity.size(); ++j) {
                    const double value = sparsity[j][i];
                    sparsity_text += (j == 0 ? "" : "\t") + hicx::npy::float_repr(value);
                    if (value == -1.0) {
                        ++negative;
                    } else if (value > sparsity_threshold) {
                        ++passing;
                    }
                }
                raw << hicx::chic::strip(lines[i]) << "\t" << sparsity_text << "\n";
                if (negative > 0) {
                    failed << lines[i];
                    ++failures;
                } else if (passing > 0) {
                    accepted_file << lines[i];
                    ++accepted;
                } else {
                    rejected_file << lines[i];
                    ++rejected;
                }
            }
        }

        {
            std::ofstream report(out_file_name + "_report", std::ios::binary);
            report << "# Created with chicQualityControl version " << hicx::kVersion << "\n";
            report << "# QC report for matrices: ";
            for (const std::string& matrix_path : matrices) {
                report << matrix_path << " ";
            }
            report << "\n";
            report << "#Sparsity threshold for rejection: sparsity <= "
                   << hicx::npy::float_repr(sparsity_threshold) << " are rejected.\n";
            report << "\nNumber of reference points: " << (accepted + rejected + failures) << "\n";
            report << "Number of accepted reference points: " << accepted << "\n";
            report << "Number of rejected reference points: " << rejected << "\n";
            report << "Number of faulty reference points: " << failures << "\n";
            report << "\n\nA faulty reference point is caused by the non-presence of the "
                      "chromosome in one of the given matrices.\n";
            report << "It can also be caused by the non-presence of valid Hi-C reads in a "
                      "region, especially at the chromosome ends.\n";
            report << "Please check the results of hicInfo to validate this for your data.\n";
        }

        // chicQualityControl.py:262-281: the reference points faulty in no
        // matrix, then x[i] = the sparsity column of matrix i.
        std::vector<std::string> labels;
        std::vector<std::string> columns;
        for (std::size_t j = 0; j < sparsity.size(); ++j) {
            std::vector<double> kept;
            for (std::size_t i = 0; i < count; ++i) {
                bool faulty = false;
                for (const auto& other : sparsity) {
                    faulty = faulty || other[i] == -1.0;
                }
                if (!faulty) {
                    kept.push_back(sparsity[j][i]);
                }
            }
            columns.push_back(hicx::plot::json_numbers(kept));
            labels.push_back(basename(matrices[j]));
        }
        hicx::plot::JsonObject data;
        data.add("outFileNameSparsity", hicx::plot::json_string(args.str("outFileNameSparsity")));
        data.add("outFileNameHistogram",
                 hicx::plot::json_string(args.str("outFileNameHistogram")));
        data.add("dpi", hicx::plot::json_int(args.integer("dpi")));
        data.add("sparsity", hicx::plot::json_number(sparsity_threshold));
        data.add("labels", hicx::plot::json_strings(labels));
        data.add("x", hicx::plot::json_list(columns));
        plot_json = data.str();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicQualityControl: %s\n", error.what());
        return 1;
    }

    return hicx::plot::draw("chicQualityControl", plot_json, args.opt_str("plotData"));
}
