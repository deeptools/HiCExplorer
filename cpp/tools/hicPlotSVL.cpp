// Port of hicexplorer/hicPlotSVL.py.
//
// For every chromosome of every matrix: the sum of the contacts at a distance
// of at most --distance, the sum of those further apart, and their ratio. With
// two or more matrices, the Wilcoxon rank sum p-value between every pair of
// per chromosome ratio distributions.
//
// What is reproduced, all pinned by hicexplorer/test/general/test_hicPlotSVL.py
// and by the harness cases:
//
//  1. **The data file misaligns chromosomes after a skipped one (reference
//     defect, pinned).** A chromosome whose ratio is inf or nan is dropped from
//     the value list at hicPlotSVL.py:121-122, but the writer at :253-259 still
//     indexes that list by the chromosome's position. Every value after a
//     skipped chromosome is written against the name of an earlier chromosome,
//     and the last rows stay blank. test_data/hicPlotSVL/data.txt itself shows 7
//     values for 12 chromosome names.
//  2. **The row names come from the last matrix.** chromosomes_list is the
//     loop variable of the matrix loop, so with an h5 and a cool matrix of
//     different chromosome sets the rows are named after the last one.
//  3. **The sums keep the matrix dtype.** np.sum over an int32 or int64 matrix
//     is an integer and prints as one; over a float32 matrix it is a float32
//     accumulated in float32, and the ratio of two float32 sums is a float32
//     division. The file still shows the float64 repr of those float32
//     values, because '{}'.format of a numpy float32 goes through
//     float.__format__. Float sums are numpy's pairwise sums
//     over the entries in scipy's row major order of the symmetric matrix,
//     reproduced with npy::PairwiseSumStream, which never materialises the
//     masked arrays.
//  4. **--threads does not change the output**: the chromosome list is cut
//     into contiguous pieces and the results are concatenated in order.
//     --threads 0 raises ZeroDivisionError (:159). A negative value makes
//     range(threads) empty, so no chromosome is evaluated and every data row is
//     blank; that is reproduced too.
//
// The box plot (cpp/PLAN.md tier 7, option (a)): after the p-value and data
// files, the kept ratios of every matrix go to
// plot/hicexplorer_plot/hicPlotSVL.py, which draws them with the reference's
// boxplot calls (hicPlotSVL.py:207-222), under --plotFileName or its default
// plot.png. The C++-only option --plotData writes that data as JSON instead.
//
// Deliberate deviations:
//
//  * **A failing chromosome exits instead of hanging.** Every chromosome is
//    evaluated in a multiprocessing.Process, and an exception there (a
//    chromosome that is not in the matrix, a bin size of zero, explicit zeros
//    stored in the matrix) ends the worker without putting anything on its
//    queue, so the parent polls forever (:183-196). The port exits 1 with the
//    cause and writes nothing. Approved by the orchestrating session
//    2026-09-13.
//
// Threading: none, deliberately. The work per chromosome is one pass over its
// stored entries; the time goes into reading the matrix, which the HDF5
// library serialises anyway. Measured on gm12878_chr1.cool, the pass is a
// small fraction of the load. cpp/OPTIMIZATION.md 6 asks for a measured gain
// before threading, and there is none to measure. --threads is accepted and
// validated for compatibility.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicPlotSVL --matrices MATRICES [MATRICES ...]\n"
    "                  [--plotFileName PLOTFILENAME] [--outFileName OUTFILENAME]\n"
    "                  [--outFileNameData OUTFILENAMEDATA] [--distance DISTANCE]\n"
    "                  [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                  [--threads THREADS] [--dpi DPI]\n"
    "                  [--colorList COLORLIST [COLORLIST ...]] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Plots the relation between short and long range interactions as boxplots and if "
    "more than one matrix is given, p-values of the distributions are computed.\n"
    "An example usage is:\n"
    "$ hicPlotSVL -m hmec_10kb.cool nhek_10kb.cool\n"
    "\n"
    "The datapoints per sample are the ratios per chromosome.\n"
    "\n"
    "Required arguments:\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        The matrix (or multiple matrices) to use for the\n"
    "                        comparison\n"
    "\n"
    "Optional arguments:\n"
    "  --plotFileName PLOTFILENAME, -pfn PLOTFILENAME\n"
    "                        Plot name (Default: plot.png).\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File the p-values are written to, p-values are only\n"
    "                        computed if at least two matrices are given (Default:\n"
    "                        p_values.txt).\n"
    "  --outFileNameData OUTFILENAMEDATA, -od OUTFILENAMEDATA\n"
    "                        File the computed ratios are written to (Default:\n"
    "                        data.txt).\n"
    "  --distance DISTANCE, -d DISTANCE\n"
    "                        Distance (in bp) which should be considered as short\n"
    "                        range. Default 2MB (2000000).\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes to include in the analysis. If not set,\n"
    "                        all chromosomes are included.\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Number of threads. Using the python multiprocessing\n"
    "                        module (Default: 4).\n"
    "  --dpi DPI             Optional parameter: Resolution for the image in case\n"
    "                        theoutput is a raster graphics image (e.g png, jpg)\n"
    "                        (Default: 300).\n"
    "  --colorList COLORLIST [COLORLIST ...], -cl COLORLIST [COLORLIST ...]\n"
    "                        Colorlist for the boxplots (Default: g b c m y k).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ port: the ratios and p-values are computed in C++, and the box plot is\n"
    "drawn by the hicexplorer_plot drawing layer with the matplotlib calls of the\n"
    "Python tool (HICX_PLOT_PYTHON names the interpreter). The C++-only option\n"
    "--plotData FILE writes the data of the plot as JSON to FILE instead of drawing\n"
    "it.\n";

struct Arguments {
    std::vector<std::string> matrices;
    std::string plot_file_name = "plot.png";
    std::optional<std::string> plot_data;
    std::string out_file_name = "p_values.txt";
    std::string out_file_name_data = "data.txt";
    std::int64_t distance = 2000000;
    std::optional<std::vector<std::string>> chromosomes;
    std::int64_t threads = 4;
    std::int64_t dpi = 300;
    std::vector<std::string> color_list = {"g", "b", "c", "m", "y", "k"};
};

// hicPlotSVL.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    namespace json = hicx::json;
    cli::Parser parser("hicPlotSVL",
                       "Plots the relation between short and long range interactions as boxplots "
                       "and if more than one matrix is given, p-values of the distributions are "
                       "computed.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrices", "-m"})
        .nargs("+")
        .required()
        .input({"h5", "cool"})
        .help("The matrix (or multiple matrices) to use for the comparison");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--plotFileName", "-pfn"})
        .default_value("plot.png")
        .output({"png", "pdf", "svg"})
        .help("Plot name.");
    optional.add({"--outFileName", "-o"})
        .default_value("p_values.txt")
        .output({"txt"})
        .help("File the p-values are written to, p-values are only computed if at least two "
              "matrices are given.");
    optional.add({"--outFileNameData", "-od"})
        .default_value("data.txt")
        .output({"txt"})
        .help("File the computed ratios are written to.");
    optional.add({"--distance", "-d"})
        .default_value(2000000)
        .type("int")
        .help("Distance (in bp) which should be considered as short range.");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("Chromosomes to include in the analysis. If not set, all chromosomes are included.");
    optional.add({"--threads", "-t"}).default_value(4).type("int").help("Number of threads.");
    optional.add({"--dpi"})
        .type("int")
        .default_value(300)
        .help("Resolution for the image in case the output is a raster graphics image.");
    optional.add({"--colorList", "-cl"})
        .default_value(json::Value::array({json::Value::string("g"), json::Value::string("b"),
                                           json::Value::string("c"), json::Value::string("m"),
                                           json::Value::string("y"), json::Value::string("k")}))
        .type("str")
        .nargs("+")
        .help("Colorlist for the boxplots.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the box plot as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the box plot is drawn from as JSON to this file and do not draw "
              "the plot.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrices = ns.strs("matrices");
    args.plot_file_name = ns.str("plotFileName");
    args.plot_data = ns.opt_str("plotData");
    args.out_file_name = ns.str("outFileName");
    args.out_file_name_data = ns.str("outFileNameData");
    args.distance = ns.integer("distance");
    if (ns.given("chromosomes")) {
        args.chromosomes = ns.strs("chromosomes");
    }
    args.threads = ns.integer("threads");
    args.dpi = ns.integer("dpi");
    args.color_list = ns.strs("colorList");
    return args;
}

// A failure inside a Python worker process, after which the reference never
// terminates.
class WorkerFailure : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

struct ChromosomeSvl {
    bool kept = false;
    double ratio_value = 0.0;  // the ratio as scipy's ranksums sees it
    bool float32 = false;      // a numpy float32 in the reference
    std::string ratio;
    std::string smaller;
    std::string greater;
};

// The entries of matrix[first:last, first:last] with the lower triangle
// filled in, in scipy's row major order, exact zeros skipped, passed to
// visit(row, column, value).
//
// For Symmetry::UpperTriangle a row's entries below the diagonal are the
// column's entries above it, so they are collected per column first. That
// transposed copy of the block is the one allocation this tool makes beyond
// the matrix itself, and it is only needed when the order of summation
// matters, which is for float matrices.
template <class Visit>
void visit_block_row_major(const hicx::CsrMatrix& matrix, std::int64_t first,
                           std::int64_t last, Visit&& visit) {
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    const auto n = static_cast<std::size_t>(last - first);
    if (matrix.symmetry() != hicx::Symmetry::UpperTriangle) {
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                const std::int64_t column = indices[k];
                if (column >= first && column < last && data[k] != 0.0) {
                    visit(row, column, data[k]);
                }
            }
        }
        return;
    }
    std::vector<std::size_t> column_start(n + 1, 0);
    for (std::int64_t row = first; row < last; ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column > row && column < last && data[k] != 0.0) {
                ++column_start[static_cast<std::size_t>(column - first) + 1];
            }
        }
    }
    for (std::size_t c = 0; c < n; ++c) {
        column_start[c + 1] += column_start[c];
    }
    std::vector<std::int32_t> lower_column(column_start[n]);
    std::vector<double> lower_value(column_start[n]);
    std::vector<std::size_t> cursor(column_start.begin(), column_start.end() - 1);
    for (std::int64_t row = first; row < last; ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column > row && column < last && data[k] != 0.0) {
                const std::size_t slot = cursor[static_cast<std::size_t>(column - first)]++;
                lower_column[slot] = static_cast<std::int32_t>(row);
                lower_value[slot] = data[k];
            }
        }
    }
    for (std::int64_t row = first; row < last; ++row) {
        const auto local = static_cast<std::size_t>(row - first);
        for (std::size_t p = column_start[local]; p < column_start[local + 1]; ++p) {
            visit(row, static_cast<std::int64_t>(lower_column[p]), lower_value[p]);
        }
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int64_t column = indices[k];
            if (column >= row && column < last && data[k] != 0.0) {
                visit(row, column, data[k]);
            }
        }
    }
}

// compute_relation_short_long_range (:95-127) for one chromosome block.
ChromosomeSvl svl_for_block(const hicx::CsrMatrix& matrix, std::int64_t first,
                            std::int64_t last, double max_distance,
                            const std::string& chromosome, const std::string& path) {
    const auto& indptr = matrix.indptr();
    const auto& indices = matrix.indices();
    const auto& data = matrix.data();
    const bool upper = matrix.symmetry() == hicx::Symmetry::UpperTriangle;

    if (!upper) {
        // hic_matrix.nonzero() skips stored zeros but hic_matrix.data does
        // not, so the boolean mask is shorter than the array it indexes and
        // numpy raises IndexError inside the worker.
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                if (indices[k] >= first && indices[k] < last && data[k] == 0.0) {
                    throw WorkerFailure(
                        "chromosome '" + chromosome + "' of " + path +
                        " stores explicit zeros, which make hicPlotSVL.py:118 index the "
                        "matrix data with a mask of a different length");
                }
            }
        }
    }
    const auto is_short = [max_distance](std::int64_t row, std::int64_t column) {
        const std::int64_t distance = row > column ? row - column : column - row;
        return static_cast<double>(distance) <= max_distance;
    };

    ChromosomeSvl result;
    const hicx::DType kind = matrix.dtype_kind();
    if (kind == hicx::DType::Integer) {
        // Integer sums do not depend on the order of the additions, so the
        // stored triangle is enough. numpy's int64 accumulator wraps, hence
        // the unsigned arithmetic.
        std::uint64_t smaller = 0;
        std::uint64_t greater = 0;
        for (std::int64_t row = first; row < last; ++row) {
            for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                 k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
                const std::int64_t column = indices[k];
                if (column < first || column >= last || data[k] == 0.0) {
                    continue;
                }
                const auto value = static_cast<std::uint64_t>(static_cast<std::int64_t>(data[k]));
                const std::uint64_t times = (upper && column != row) ? 2U : 1U;
                (is_short(row, column) ? smaller : greater) += times * value;
            }
        }
        const auto small_sum = static_cast<std::int64_t>(smaller);
        const auto large_sum = static_cast<std::int64_t>(greater);
        const double ratio = static_cast<double>(small_sum) / static_cast<double>(large_sum);
        if (std::isinf(ratio) || std::isnan(ratio)) {
            return result;
        }
        result.kept = true;
        result.ratio_value = ratio;
        result.ratio = hicx::npy::float_repr(ratio);
        result.smaller = std::to_string(small_sum);
        result.greater = std::to_string(large_sum);
        return result;
    }

    if (kind == hicx::DType::Float32) {
        hicx::npy::PairwiseSumStream<float> smaller;
        hicx::npy::PairwiseSumStream<float> greater;
        visit_block_row_major(matrix, first, last,
                              [&](std::int64_t row, std::int64_t column, double value) {
                                  (is_short(row, column) ? smaller : greater)
                                      .add(static_cast<float>(value));
                              });
        const float small_sum = smaller.result();
        const float large_sum = greater.result();
        const float ratio = small_sum / large_sum;
        if (std::isinf(ratio) || std::isnan(ratio)) {
            return result;
        }
        result.kept = true;
        result.ratio_value = static_cast<double>(ratio);
        result.float32 = true;
        // '{}'.format(np.float32(x)) goes through float.__format__, so the
        // file shows the float64 repr of the widened value
        // (48.16666793823242), not numpy's float32 str (48.166668).
        result.ratio = hicx::npy::float_repr(static_cast<double>(ratio));
        result.smaller = hicx::npy::float_repr(static_cast<double>(small_sum));
        result.greater = hicx::npy::float_repr(static_cast<double>(large_sum));
        return result;
    }

    hicx::npy::PairwiseSumStream<double> smaller;
    hicx::npy::PairwiseSumStream<double> greater;
    visit_block_row_major(matrix, first, last,
                          [&](std::int64_t row, std::int64_t column, double value) {
                              (is_short(row, column) ? smaller : greater).add(value);
                          });
    const double small_sum = smaller.result();
    const double large_sum = greater.result();
    const double ratio = small_sum / large_sum;
    if (std::isinf(ratio) || std::isnan(ratio)) {
        return result;
    }
    result.kept = true;
    result.ratio_value = ratio;
    result.ratio = hicx::npy::float_repr(ratio);
    result.smaller = hicx::npy::float_repr(small_sum);
    result.greater = hicx::npy::float_repr(large_sum);
    return result;
}

double max_distance_for(std::int64_t distance, std::int64_t bin_size,
                        const std::string& chromosome, const std::string& path) {
    if (bin_size == 0) {
        throw WorkerFailure("the bin size of " + path + " at chromosome '" + chromosome +
                            "' is 0, so hicPlotSVL.py:105/112 divides by zero");
    }
    return static_cast<double>(distance) / static_cast<double>(bin_size);
}

std::string missing_chromosome(const std::string& chromosome, const std::string& path) {
    return "chromosome '" + chromosome + "' is not in " + path;
}

void write_or_throw(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("cannot open " + path + " for writing");
    }
    out << content;
    if (!out.flush()) {
        throw std::runtime_error("cannot write " + path);
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    if (const int refused = hicx::plot::preflight("hicPlotSVL", !args.plot_data.has_value()); refused != 0) {
        return refused;
    }

    std::vector<std::vector<ChromosomeSvl>> per_matrix;
    std::vector<std::string> chromosomes_list;
    try {
        for (const std::string& path : args.matrices) {
            const bool is_cooler = hicx::check_cooler(path);
            std::optional<hicx::ToolMatrix> whole;
            std::vector<std::string> chromosomes;
            if (!is_cooler) {
                whole = hicx::ToolMatrix::load(path);
                for (const auto& entry : whole->boundaries()) {
                    chromosomes.push_back(entry.first);
                }
            } else {
                chromosomes = hicx::CoolFile(path).chrom_names();
            }
            if (args.chromosomes.has_value()) {
                chromosomes = *args.chromosomes;
            }
            if (args.threads == 0) {
                throw std::runtime_error(
                    "ZeroDivisionError: integer division or modulo by zero "
                    "(hicPlotSVL.py:159 divides the number of chromosomes by --threads)");
            }

            std::vector<ChromosomeSvl> kept;
            // range(args.threads) is empty for a negative --threads, so the
            // reference starts no worker and evaluates no chromosome.
            if (args.threads > 0) {
                std::optional<hicx::BinTable> whole_bins;
                std::optional<std::vector<std::string>> cool_chromosomes;
                if (whole.has_value()) {
                    whole_bins.emplace(whole->cut_intervals());
                } else {
                    cool_chromosomes = hicx::CoolFile(path).chrom_names();
                }
                for (const std::string& chromosome : chromosomes) {
                    ChromosomeSvl svl;
                    if (whole.has_value()) {
                        const std::optional<hicx::BinRange> range =
                            whole_bins->chrom_bin_range(chromosome);
                        if (!range.has_value()) {
                            throw WorkerFailure(missing_chromosome(chromosome, path));
                        }
                        const double max_distance = max_distance_for(
                            args.distance, whole_bins->bin_size(), chromosome, path);
                        svl = svl_for_block(whole->matrix(), range->first, range->last,
                                            max_distance, chromosome, path);
                    } else {
                        if (std::find(cool_chromosomes->begin(), cool_chromosomes->end(),
                                      chromosome) == cool_chromosomes->end()) {
                            throw WorkerFailure(missing_chromosome(chromosome, path));
                        }
                        const hicx::ToolMatrix block = hicx::ToolMatrix::load(path, chromosome);
                        const hicx::BinTable bins(block.cut_intervals());
                        const double max_distance =
                            max_distance_for(args.distance, bins.bin_size(), chromosome, path);
                        svl = svl_for_block(block.matrix(), 0, block.matrix().rows(),
                                            max_distance, chromosome, path);
                    }
                    if (svl.kept) {
                        kept.push_back(std::move(svl));
                    }
                }
            }
            per_matrix.push_back(std::move(kept));
            chromosomes_list = std::move(chromosomes);
        }
    } catch (const WorkerFailure& failure) {
        std::fprintf(stderr,
                     "hicPlotSVL: %s. The Python reference raises inside a worker process "
                     "at this point and then waits for that worker forever "
                     "(hicPlotSVL.py:183-196); the C++ port exits instead. Nothing was "
                     "written.\n",
                     failure.what());
        return 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotSVL: %s\n", error.what());
        return 1;
    }

    const std::string distance = std::to_string(args.distance);
    try {
        if (args.matrices.size() > 1) {
            std::string text = "# Created with HiCExplorer's hicPlotSVL " +
                               std::string(hicx::kVersion) + "\n";
            text += "# Short range vs long range contacts per chromosome, p-values of each "
                    "distribution against each other distribution with Wilcoxon rank-sum\n";
            text += "# Short range contacts: <= " + distance + "\n";
            for (std::size_t i = 0; i < per_matrix.size(); ++i) {
                std::vector<double> sample;
                for (const ChromosomeSvl& svl : per_matrix[i]) {
                    sample.push_back(svl.ratio_value);
                }
                for (std::size_t j = i + 1; j < per_matrix.size(); ++j) {
                    std::vector<double> sample2;
                    for (const ChromosomeSvl& svl : per_matrix[j]) {
                        sample2.push_back(svl.ratio_value);
                    }
                    const hicx::stats::RanksumsResult test = hicx::stats::ranksums(
                        std::span<const double>(sample), std::span<const double>(sample2));
                    text += args.matrices[i] + "\t" + args.matrices[j] + "\t" +
                            hicx::npy::float_repr(test.pvalue) + "\n";
                }
            }
            write_or_throw(args.out_file_name, text);
        }

        std::string text =
            "# Created with HiCExplorer's hicPlotSVL " + std::string(hicx::kVersion) + "\n";
        text += "# Short range vs long range contacts per chromosome: raw data\n";
        text += "# Short range contacts: <= " + distance + "\n";
        text += "#\t";
        for (std::size_t i = 0; i < args.matrices.size(); ++i) {
            text += (i > 0 ? "\t\t\t" : "") + args.matrices[i];
        }
        text += "\n# Chromosome\t";
        for (std::size_t i = 0; i < args.matrices.size(); ++i) {
            text += (i > 0 ? "\t" : "") + std::string("Ratio\tSum <= ") + distance +
                    "\tSum > " + distance;
        }
        text += "\n";
        // Pinned reference defect 1 of the file comment: row i takes the i-th
        // *kept* value of every matrix, whatever chromosome it came from.
        for (std::size_t i = 0; i < chromosomes_list.size(); ++i) {
            text += chromosomes_list[i] + "\t";
            for (const auto& values : per_matrix) {
                if (i < values.size()) {
                    text += values[i].ratio + "\t" + values[i].smaller + "\t" +
                            values[i].greater + "\t";
                } else {
                    text += "\t";
                }
            }
            text += "\n";
        }
        write_or_throw(args.out_file_name_data, text);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicPlotSVL: %s\n", error.what());
        return 1;
    }

    // plt.boxplot(short_v_long_range) at :216, one sample per matrix.
    std::vector<std::string> samples;
    for (const auto& values : per_matrix) {
        std::vector<double> ratios;
        bool float32 = false;
        for (const ChromosomeSvl& svl : values) {
            ratios.push_back(svl.ratio_value);
            float32 = float32 || svl.float32;
        }
        hicx::plot::JsonObject sample;
        sample.add("values", hicx::plot::json_numbers(ratios));
        sample.add("float32", hicx::plot::json_bool(float32));
        samples.push_back(sample.str());
    }
    hicx::plot::JsonObject data;
    data.add("plotFileName", hicx::plot::json_string(args.plot_file_name));
    data.add("dpi", hicx::plot::json_int(args.dpi));
    data.add("colorList", hicx::plot::json_strings(args.color_list));
    data.add("matrices", hicx::plot::json_strings(args.matrices));
    data.add("samples", hicx::plot::json_list(samples));
    hicx::report_resource_usage("hicPlotSVL");
    return hicx::plot::draw("hicPlotSVL", data.str(), args.plot_data);
}
