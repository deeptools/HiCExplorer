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
// Figures. The Python also draws sparsity.png and histogram.png with
// matplotlib. The C++ port draws no figures (tiers 7 and 8 await the project
// owner's decision). Under the project rule, a figure file named explicitly
// with --outFileNameHistogram or --outFileNameSparsity makes the tool refuse
// before it writes anything; a figure the user did not name is skipped with a
// note on stderr.
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

#include "chic_arguments.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/version.hpp"

namespace {

using hicx::chic::ReferencePoint;

const char* const kUsage =
    "usage: chicQualityControl --matrices MATRICES [MATRICES ...]\n"
    "                          --referencePoints REFERENCEPOINTS --sparsity SPARSITY\n"
    "                          [--outFileName OUTFILENAME]\n"
    "                          [--outFileNameHistogram OUTFILENAMEHISTOGRAM]\n"
    "                          [--outFileNameSparsity OUTFILENAMESPARSITY]\n"
    "                          [--threads THREADS] [--fixateRange FIXATERANGE]\n"
    "                          [--dpi DPI] [--help] [--version]\n";

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
    using hicx::chic_cli::Arity;
    using hicx::chic_cli::Option;
    using hicx::chic_cli::Type;

    hicx::chic_cli::Parser parser(
        "chicQualityControl", kUsage,
        "Computes the sparsity of each viewpoint to determine the quality.",
        hicx::kVersion);
    parser.add(Option{{"--matrices", "-m"}, "matrices", Arity::OneOrMore, Type::String, true, {},
                      "The input matrices to apply the QC on."});
    parser.add(Option{{"--referencePoints", "-rp"}, "referencePoints", Arity::One, Type::String,
                      true, {}, "Bed file with the reference points."});
    parser.add(Option{{"--sparsity", "-s"}, "sparsity", Arity::One, Type::Float, true, {},
                      "Sparsity threshold."});
    parser.add(Option{{"--outFileName", "-o"}, "outFileName", Arity::One, Type::String, false,
                      {"new_referencepoints.bed"}, "Output file of the passed reference points."});
    parser.add(Option{{"--outFileNameHistogram", "-oh"}, "outFileNameHistogram", Arity::One,
                      Type::String, false, {"histogram.png"}, "Histogram plot."});
    parser.add(Option{{"--outFileNameSparsity", "-os"}, "outFileNameSparsity", Arity::One,
                      Type::String, false, {"sparsity.png"}, "Sparsity distribution plot."});
    parser.add(Option{{"--threads", "-t"}, "threads", Arity::One, Type::Int, false, {"4"},
                      "Number of threads."});
    parser.add(Option{{"--fixateRange", "-fs"}, "fixateRange", Arity::One, Type::Int, false,
                      {"500000"}, "Range on either side of a reference point."});
    parser.add(Option{{"--dpi"}, "dpi", Arity::One, Type::Int, false, {"300"},
                      "Resolution of the plots."});
    const hicx::chic_cli::Parsed args = parser.parse(argc, argv);

    // The project rule for figures: refuse before writing anything when one
    // was requested by name.
    for (const char* dest : {"outFileNameHistogram", "outFileNameSparsity"}) {
        if (args.explicitly_given(dest)) {
            std::fprintf(stderr,
                         "chicQualityControl: error: --%s %s was requested, but plotting is not "
                         "yet available in the C++ port of HiCExplorer. No output was written. "
                         "Use the Python chicQualityControl for the figures.\n",
                         dest, args.str(dest).c_str());
            return 1;
        }
    }

    const std::vector<std::string>& matrices = args.list("matrices");
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
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicQualityControl: %s\n", error.what());
        return 1;
    }

    std::fprintf(stderr,
                 "chicQualityControl: note: the sparsity and histogram figures (%s, %s) were not "
                 "drawn; plotting is not yet available in the C++ port of HiCExplorer.\n",
                 args.str("outFileNameSparsity").c_str(), args.str("outFileNameHistogram").c_str());
    return 0;
}
