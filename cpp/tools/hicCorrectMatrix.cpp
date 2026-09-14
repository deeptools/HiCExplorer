// Port of hicexplorer/hicCorrectMatrix.py, the `correct` subcommand.
//
// Two balancing methods, Knight-Ruiz and iterative correction, over four
// combinations of --perchr and output format. The mathematics lives in
// core/src/math/{ice,kr_balancing,sparse_kernels}.cpp; this file is the
// argument handling, the bin filtering and, above all, the call ordering,
// because the call ordering is observable.
//
// The four combinations, and why they differ (finding F4 of cpp/STATUS.md,
// confirmed as intended behaviour by the project owner on 2026-09-01):
//
//   The output format decides what a corrected matrix means. HiCExplorer's h5
//   format historically could not carry correction factors, so for an output
//   name ending in '.h5' the correction is applied to the values and the file
//   holds an already corrected matrix. cool carries the factors as a weight
//   column, so for any other name the raw matrix is written and the factors are
//   stored beside it to be applied on read.
//
//   With --perchr that split has a second, mechanical consequence.
//   hicCorrectMatrix.py:726 calls get_normalised_matrix(True) only for an '.h5'
//   name, and that call is also what triggers krbalancing's rescale_norm_vector.
//   The following get_normalisation_vector(False) at :731 therefore returns
//   *rescaled* factors for '.h5' and *unrescaled* ones for anything else. The
//   whole matrix path at :754 calls get_normalisation_vector(True)
//   unconditionally, so there the factors are rescaled for both formats and only
//   the matrix values differ. This port reproduces the call ordering, not only
//   the arithmetic, because the rescaling is an in place side effect on the
//   normalisation vector.
//
//   ICE has its own asymmetry: the whole matrix path calls setMatrixValues
//   unconditionally (:740), so a whole matrix ICE run writes corrected values
//   whatever the output name says, while --perchr assembles its result into a
//   lil_matrix that is only read back for an '.h5' name.
//
// Two further consequences of --perchr that the port has to reproduce: the
// lil_matrix is only ever written inside the diagonal chromosome blocks
// (:711,728), so every inter chromosomal contact is absent from a --perchr
// output; and Knight-Ruiz balances A = M + 1e-5 * I and returns triu(A), so its
// output carries an entry on every position of the main diagonal even where the
// input had none.
//
// Deliberate deviations, all recorded for cpp/STATUS.md:
//
//  * krbalancing calls exit(0) after 300 outer iterations, having printed the
//    whole normalisation vector to stdout. That leaves a caller with a success
//    status and no output file. Reproduced in neither mode: this port reports
//    the failure and exits 1.
//  * --compatMode {v3,v4} does not exist in the Python. v4 is the default and
//    is float64 throughout; v3 reproduces krbalancing's float32 input rounding
//    and float32 rescaling accumulators, so that the cost of those two defects
//    can be measured on real data. v3 sums in a fixed order, so it is
//    deterministic where the reference is not.
//  * --threads does not exist in the Python either. Any thread count produces
//    the same bytes, see cpp/OPTIMIZATION.md section 3.
//  * `diagnostic_plot` is a matplotlib subcommand and belongs to tier 7. It
//    fails with an explicit message rather than half working.
//
// Reproduced faithfully, including the failures:
//
//  * --inflationCutoff without --transCutoff raises UnboundLocalError in the
//    Python, because pre_row_sum is only assigned inside the --transCutoff
//    branch (hicCorrectMatrix.py:699 against :770). Exit 1 with the reason.
//  * --inflationCutoff with --transCutoff then reaches the final
//    printchrtoremove with bin ids from before the masking, which raises
//    IndexError as soon as the union it is given differs from the previous
//    call's argument. Exit 1 with the reason.
//  * --sequencedCountCutoff asserts that the per bin coverage is a np.float64
//    (:678). A cool file's cut intervals carry the Python float 1.0, so the
//    assertion fails and the tool exits non zero on every cool input.
//  * --transCutoff calls truncTrans, which is a no operation (finding F11).

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/correct_ops.hpp"
#include "hicx/math/ice.hpp"
#include "hicx/math/kr_balancing.hpp"
#include "hicx/math/sparse_kernels.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: hicCorrectMatrix [-h] [--version]  ...\n";

const char* const kCorrectUsage =
    "usage: hicCorrectMatrix correct --matrix MATRIX --outFileName OUTFILENAME\n"
    "                                [--correctionMethod STR]\n"
    "                                [--filterThreshold FILTERTHRESHOLD FILTERTHRESHOLD]\n"
    "                                [--iterNum INT] [--inflationCutoff INFLATIONCUTOFF]\n"
    "                                [--transCutoff TRANSCUTOFF]\n"
    "                                [--sequencedCountCutoff SEQUENCEDCOUNTCUTOFF]\n"
    "                                [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                                [--skipDiagonal] [--perchr] [--filteredBed FILTEREDBED]\n"
    "                                [--verbose] [--compatMode {v3,v4}] [--threads INT]\n"
    "                                [--help]\n";

const char* const kHelp =
    "\n"
    "This function provides 2 balancing methods which can be applied on a raw\n"
    "matrix.\n"
    "\n"
    "I. KR: It balances a matrix using a fast balancing algorithm introduced by\n"
    "Knight and Ruiz (2012).\n"
    "\n"
    "II. ICE: Iterative correction of a Hi-C matrix (see Imakaev et al. 2012\n"
    "Nature Methods for details).\n"
    "\n"
    "Options:\n"
    "  correct               Run Knight-Ruiz matrix balancing algorithm (KR) or the\n"
    "                        iterative matrix correction (ICE).\n"
    "  diagnostic_plot       Plots a histogram of the coverage per bin together with\n"
    "                        the modified z-score. Not implemented in the C++ port,\n"
    "                        see cpp/PLAN.md tier 7.\n"
    "\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

const char* const kCorrectHelp =
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        Name of the Hi-C matrix to correct in .h5 format.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        File name to save the resulting matrix. The output is a\n"
    "                        .h5 file.\n"
    "\n"
    "Optional arguments:\n"
    "  --correctionMethod STR\n"
    "                        Method to be used for matrix correction. It can be set to\n"
    "                        KR or ICE (Default: KR).\n"
    "  --filterThreshold FILTERTHRESHOLD FILTERTHRESHOLD, -t ...\n"
    "                        Removes bins of low or large coverage. Applied only for\n"
    "                        ICE!\n"
    "  --iterNum INT, -n INT\n"
    "                        Number of iterations to compute. Only for ICE!\n"
    "                        (Default: 500).\n"
    "  --inflationCutoff INFLATIONCUTOFF\n"
    "                        Maximum number of times a bin can be scaled up during the\n"
    "                        iterative correction. Only for ICE!\n"
    "  --transCutoff TRANSCUTOFF, -transcut TRANSCUTOFF\n"
    "                        Clip high counts in the top -transcut trans regions.\n"
    "                        Only for ICE!\n"
    "  --sequencedCountCutoff SEQUENCEDCOUNTCUTOFF\n"
    "                        Discard bins covered by fewer reads than this fraction.\n"
    "                        Only for ICE!\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to be included in the correction.\n"
    "  --skipDiagonal, -s    If set, diagonal counts are not included. Only for ICE!\n"
    "  --perchr              Normalize each chromosome separately.\n"
    "  --filteredBed FILTEREDBED\n"
    "                        Print bins filtered out by --filterThreshold to this file\n"
    "  --verbose             Print processing status.\n"
    "  --compatMode {v3,v4}  v4 (default) balances in float64. v3 reproduces\n"
    "                        krbalancing's float32 input rounding and float32\n"
    "                        rescaling accumulators. Only for KR. Not a Python option.\n"
    "  --threads INT         Worker threads. The result does not depend on this.\n"
    "                        Not a Python option (Default: 4).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

enum class Method { Kr, Ice };

struct Arguments {
    std::string command;
    std::string matrix;
    std::string out_file_name;
    std::string plot_name;
    Method method = Method::Kr;
    std::optional<std::pair<double, double>> filter_threshold;
    std::int64_t iter_num = 500;
    std::optional<double> inflation_cutoff;
    std::optional<double> trans_cutoff;
    std::optional<double> sequenced_count_cutoff;
    std::optional<double> x_max;
    std::vector<std::string> chromosomes;
    bool has_chromosomes = false;
    bool skip_diagonal = false;
    bool perchr = false;
    std::optional<std::string> filtered_bed;
    bool verbose = false;
    bool compat_v3 = false;
    int threads = 4;
};

const char* const kDiagnosticUsage =
    "usage: hicCorrectMatrix diagnostic_plot --matrix hic_matrix.h5 -o file.png\n";

const char* const kDiagnosticHelp =
    "\n"
    "options:\n"
    "  -h, --help            show this help message and exit\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        Name of the Hi-C matrix to correct in .h5 format.\n"
    "                        (default: None)\n"
    "  --plotName PLOTNAME, -o PLOTNAME\n"
    "                        File name to save the diagnostic plot. (default: None)\n"
    "\n"
    "Optional arguments:\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        List of chromosomes to be included in the iterative\n"
    "                        correction. The order of the given chromosomes will be\n"
    "                        then kept for the resulting corrected matrix.\n"
    "                        (default: None)\n"
    "  --xMax XMAX           Max value for the x-axis in counts per bin. (default:\n"
    "                        None)\n"
    "  --perchr              Compute histogram per chromosome. For samples from\n"
    "                        cells with uneven number of chromosomes and/or\n"
    "                        translocations it is advisable to check the histograms\n"
    "                        per chromosome to find the most conservative\n"
    "                        `filterThreshold`. (default: False)\n"
    "  --verbose             Print processing status. (default: False)\n";

[[noreturn]] void fail(const std::string& message) {
    std::fputs(kCorrectUsage, stderr);
    std::fprintf(stderr, "hicCorrectMatrix correct: error: %s\n", message.c_str());
    std::exit(2);
}

// hicCorrectMatrix.py parse_arguments and correct_subparser.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicCorrectMatrix",
                       "This function provides 2 balancing methods which can be applied on a raw "
                       "matrix: the Knight-Ruiz balancing (KR) and the iterative correction of "
                       "Imakaev et al. (ICE).");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& options = parser.group("options");
    options.add({"-h", "--help"})
        .action(cli::Action::Help)
        .help("show this help message and exit");
    options.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    parser.subcommands("command", false, std::string());

    cli::Parser& plot = parser.add_subcommand(
        "diagnostic_plot",
        "Plots a histogram of the coverage per bin together with the modified z-score based on "
        "the median absolute deviation method.");
    plot.set_usage(kDiagnosticUsage).set_help(kDiagnosticHelp);
    plot.group("options")
        .add({"-h", "--help"})
        .action(cli::Action::Help)
        .help("show this help message and exit");
    cli::ArgumentGroup& plot_required = plot.group("Required arguments");
    plot_required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Name of the Hi-C matrix to correct in .h5 format.");
    plot_required.add({"--plotName", "-o"})
        .required()
        .output({"png", "pdf", "svg"})
        .help("File name to save the diagnostic plot.");
    cli::ArgumentGroup& plot_optional = plot.group("Optional arguments");
    plot_optional.add({"--chromosomes"})
        .nargs("+")
        .help("List of chromosomes to be included in the iterative correction.");
    plot_optional.add({"--xMax"}).type("float").help("Max value for the x-axis in counts per bin.");
    plot_optional.add({"--perchr"})
        .action(cli::Action::StoreTrue)
        .help("Compute histogram per chromosome.");
    plot_optional.add({"--verbose"})
        .action(cli::Action::StoreTrue)
        .help("Print processing status.");

    cli::Parser& correct = parser.add_subcommand(
        "correct",
        "Run Knight-Ruiz matrix balancing algorithm (KR) or the iterative matrix correction (ICE).");
    correct.set_usage(kCorrectUsage).set_help(kCorrectHelp);
    correct.group("options")
        .add({"-h", "--help"})
        .action(cli::Action::Help)
        .help("show this help message and exit");
    cli::ArgumentGroup& required = correct.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Name of the Hi-C matrix to correct in .h5 format.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"h5", "cool"})
        .help("File name to save the resulting matrix.");
    cli::ArgumentGroup& optional = correct.group("Optional arguments");
    optional.add({"--correctionMethod"})
        .type("str")
        .metavar("STR")
        .default_value("KR")
        .choices({"KR", "ICE"})
        .help("Method to be used for matrix correction.");
    optional.add({"--filterThreshold", "-t"})
        .type("float")
        .nargs(2)
        .help("Removes bins of low or large coverage. Applied only for ICE!");
    optional.add({"--iterNum", "-n"})
        .type("int")
        .metavar("INT")
        .default_value(500)
        .help("Number of iterations to compute. Only for ICE!");
    optional.add({"--inflationCutoff"})
        .type("float")
        .help("Maximum number of times a bin can be scaled up during the iterative correction.");
    optional.add({"--transCutoff", "-transcut"})
        .type("float")
        .help("Clip high counts in the top -transcut trans regions. Only for ICE!");
    optional.add({"--sequencedCountCutoff"})
        .type("float")
        .help("Discard bins covered by reads on less than this fraction. Only for ICE!");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("List of chromosomes to be included in the iterative correction.");
    optional.add({"--skipDiagonal", "-s"})
        .action(cli::Action::StoreTrue)
        .help("If set, diagonal counts are not included. Only for ICE!");
    optional.add({"--perchr"})
        .action(cli::Action::StoreTrue)
        .help("Normalize each chromosome separately.");
    optional.add({"--filteredBed"})
        .output({"bed"})
        .help("Print bins filtered out by --filterThreshold to this file.");
    optional.add({"--verbose"})
        .action(cli::Action::StoreTrue)
        .help("Print processing status.");
    // argparse expands %(prog)s to the sub-parser's prog here.
    optional.add({"--version"}).version(std::string("hicCorrectMatrix correct ") + hicx::kVersion);
    optional.add({"--compatMode"})
        .choices({"v3", "v4"})
        .default_value("v4")
        .cpp_only("v3 reproduces krbalancing's float32 input rounding and float32 rescaling "
                  "accumulators, to measure the cost of those defects.")
        .help("v4 balances in float64; v3 reproduces krbalancing's float32 arithmetic. Only for KR.");
    optional.add({"--threads"})
        .type("int")
        .metavar("INT")
        .default_value(4)
        .cpp_only("Worker threads of the C++ port; the result does not depend on the number.")
        .help("Worker threads. The result does not depend on this.");

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.command = ns.command();
    if (args.command.empty()) {
        // Without a subcommand argparse sets no option at all, and main reads
        // args.verbose first.
        std::fputs("hicCorrectMatrix: AttributeError: 'Namespace' object has no attribute "
                   "'verbose'. Without a subcommand no option is set, and "
                   "hicCorrectMatrix.py:599 reads args.verbose.\n",
                   stderr);
        std::exit(1);
    }
    args.matrix = ns.str("matrix");
    args.chromosomes = ns.strs("chromosomes");
    args.has_chromosomes = ns.given("chromosomes");
    args.perchr = ns.flag("perchr");
    args.verbose = ns.flag("verbose");
    if (args.command == "diagnostic_plot") {
        args.plot_name = ns.str("plotName");
        args.x_max = ns.opt_real("xMax");
        return args;
    }
    args.out_file_name = ns.str("outFileName");
    args.method = ns.str("correctionMethod") == "ICE" ? Method::Ice : Method::Kr;
    if (ns.given("filterThreshold")) {
        const std::vector<double> threshold = ns.reals("filterThreshold");
        args.filter_threshold = std::make_pair(threshold[0], threshold[1]);
    }
    args.iter_num = ns.integer("iterNum");
    args.inflation_cutoff = ns.opt_real("inflationCutoff");
    args.trans_cutoff = ns.opt_real("transCutoff");
    args.sequenced_count_cutoff = ns.opt_real("sequencedCountCutoff");
    args.skip_diagonal = ns.flag("skipDiagonal");
    args.filtered_bed = ns.opt_str("filteredBed");
    args.compat_v3 = ns.str("compatMode") == "v3";
    const std::int64_t threads = ns.integer("threads");
    if (threads < 1) {
        fail("argument --threads: must be at least 1");
    }
    args.threads = static_cast<int>(std::min<std::int64_t>(threads, 1 << 16));
    return args;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void log_info(bool verbose, const std::string& message) {
    if (verbose) {
        std::fprintf(stderr, "INFO:hicexplorer.hicCorrectMatrix:%s\n", message.c_str());
    }
}

// convertNansToZeros followed by convertInfsToZeros. The Python leaves the
// entries stored as explicit zeros and both writers drop them again with
// eliminate_zeros, so dropping them here produces the same file and the same
// marginals.
void drop_non_finite(hicx::CsrMatrix& matrix) {
    bool any = false;
    for (const double value : matrix.data()) {
        if (!std::isfinite(value)) {
            any = true;
            break;
        }
    }
    if (!any) {
        return;
    }
    for (double& value : matrix.mutable_data()) {
        if (!std::isfinite(value)) {
            value = 0.0;
        }
    }
    matrix.eliminate_zeros();
}

// The mask accumulated over the successive maskBins calls, expressed once in
// the coordinates of the original matrix. hiCMatrix reaches the same state by
// restoring the previous mask inside every further maskBins call.
struct MaskChain {
    bool active = false;
    std::vector<std::int64_t> kept;
    std::vector<hicx::CutInterval> original_intervals;
};

void apply_mask(hicx::MatrixData& data, MaskChain& chain, const std::vector<char>& masked) {
    const hicx::correct::MaskState state = hicx::correct::mask_bins_in_place(data, masked);
    if (!chain.active) {
        chain.active = true;
        chain.original_intervals = state.original_intervals;
        chain.kept = state.kept;
        return;
    }
    std::vector<std::int64_t> kept;
    kept.reserve(state.kept.size());
    for (const std::int64_t bin : state.kept) {
        kept.push_back(chain.kept[static_cast<std::size_t>(bin)]);
    }
    chain.kept = std::move(kept);
}

void restore_mask(hicx::MatrixData& data, const MaskChain& chain) {
    if (!chain.active) {
        return;
    }
    hicx::correct::MaskState state;
    state.kept = chain.kept;
    state.original_intervals = chain.original_intervals;
    std::vector<char> kept_flag(chain.original_intervals.size(), 0);
    for (const std::int64_t bin : chain.kept) {
        kept_flag[static_cast<std::size_t>(bin)] = 1;
    }
    for (std::size_t bin = 0; bin < kept_flag.size(); ++bin) {
        if (kept_flag[bin] == 0) {
            state.removed.push_back(static_cast<std::int64_t>(bin));
        }
    }
    hicx::correct::restore_masked_bins(data, state);
}

// A = M + addend * I, materialised. Only needed when the balanced matrix is
// what gets written, because that matrix carries an entry on every position of
// the diagonal.
void materialise_kr_matrix(hicx::CsrMatrix& matrix, double addend, bool round_float32) {
    if (round_float32) {
        for (double& value : matrix.mutable_data()) {
            value = static_cast<double>(static_cast<float>(value));
        }
    }
    hicx::correct::add_missing_diagonal(matrix, 0.0);
    const std::vector<std::int64_t>& indptr = matrix.indptr();
    const std::vector<std::int32_t>& indices = matrix.indices();
    std::vector<double>& data = matrix.mutable_data();
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        if (begin < end && indices[begin] == static_cast<std::int32_t>(row)) {
            data[begin] += addend;
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        if (args.command == "diagnostic_plot") {
            std::fputs(
                "hicCorrectMatrix: diagnostic_plot is not implemented in the C++ port. It "
                "draws a matplotlib histogram of the per bin coverage with a modified "
                "z-score axis, which belongs to tier 7 of cpp/PLAN.md and stays a Python "
                "plotting shell over the C++ core. Use the Python hicCorrectMatrix for "
                "this subcommand.\n",
                stderr);
            return 1;
        }
        if (args.method == Method::Kr && args.compat_v3 && args.perchr) {
            // Nothing forbids it; the note is here so that the mode is not
            // silently taken to be a whole matrix only feature.
            log_info(args.verbose, "compatMode v3 applies per chromosome as well");
        }

        const bool output_is_h5 = ends_with(args.out_file_name, ".h5");
        const int threads = args.threads;

        // hicCorrectMatrix.py:603-609.
        std::optional<hicx::ToolMatrix> loaded;
        if (hicx::check_cooler(args.matrix) && args.has_chromosomes &&
            args.chromosomes.size() == 1) {
            loaded = hicx::ToolMatrix::load(args.matrix, args.chromosomes[0]);
        } else {
            loaded = hicx::ToolMatrix::load(args.matrix);
            if (args.has_chromosomes) {
                std::vector<std::int64_t> order;
                for (const std::string& chromosome : args.chromosomes) {
                    bool found = false;
                    for (const auto& entry : loaded->boundaries()) {
                        if (entry.first != chromosome) {
                            continue;
                        }
                        for (std::int64_t bin = entry.second.first; bin < entry.second.last;
                             ++bin) {
                            order.push_back(bin);
                        }
                        found = true;
                        break;
                    }
                    if (!found) {
                        std::fprintf(stderr,
                                     "hicCorrectMatrix: ValueError: Chromosome name '%s' not "
                                     "found. Please check the correct spelling of the "
                                     "chromosomes and try again\n",
                                     chromosome.c_str());
                        return 1;
                    }
                }
                hicx::reorder_bins(loaded->data(), order);
                loaded->refresh_boundaries();
            }
        }
        hicx::ToolMatrix& hic = *loaded;
        hicx::MatrixData& data = hic.data();

        MaskChain chain;

        // ICE masks the zero coverage bins before anything else (:613-617).
        if (args.method == Method::Ice) {
            const std::vector<std::int64_t> zero_bins =
                hicx::correct::zero_coverage_bins(data.matrix, threads);
            log_info(args.verbose,
                     "Removing " + std::to_string(zero_bins.size()) + " zero value bins");
            // hiCMatrix.maskBins returns before touching anything when the
            // list is empty, so an empty zero-coverage set leaves even the
            // existing NaN bins alone. Only a non-empty call folds them in.
            if (!zero_bins.empty()) {
                std::vector<char> masked(static_cast<std::size_t>(data.matrix.rows()), 0);
                for (const std::int64_t bin : zero_bins) {
                    masked[static_cast<std::size_t>(bin)] = 1;
                }
                for (const std::int64_t bin : data.nan_bins) {
                    masked[static_cast<std::size_t>(bin)] = 1;
                }
                apply_mask(data, chain, masked);
                hic.refresh_boundaries();
            }
        }

        drop_non_finite(data.matrix);
        data.matrix.set_dtype("float64");

        if (args.skip_diagonal) {
            hicx::correct::remove_diagonal(data.matrix);
        }

        if (args.method == Method::Ice) {
            if (!args.filter_threshold.has_value()) {
                std::fputs("ERROR:hicexplorer.hicCorrectMatrix:min and max filtering "
                           "thresholds should be set\n",
                           stderr);
                return 1;
            }
            std::vector<std::string> warnings;
            const std::vector<std::int64_t> outliers = hicx::correct::filter_by_zscore(
                data.matrix, hic.boundaries(), args.filter_threshold->first,
                args.filter_threshold->second, args.perchr, threads, &warnings);
            for (const std::string& chromosome : warnings) {
                std::fprintf(stderr,
                             "WARNING:hicexplorer.hicCorrectMatrix:Warning. No bins removed "
                             "for chromosome %s using thresholds %g %g\n",
                             chromosome.c_str(), args.filter_threshold->first,
                             args.filter_threshold->second);
            }

            if (args.filtered_bed.has_value()) {
                std::ofstream bed(*args.filtered_bed);
                if (!bed) {
                    std::fprintf(stderr, "hicCorrectMatrix: cannot open %s for writing\n",
                                 args.filtered_bed->c_str());
                    return 1;
                }
                // hicCorrectMatrix.py:668 iterates set(outlier_regions), so
                // the lines come out in CPython set slot order and not in bin
                // order. The checked-in reference filtered.bed has that order
                // and the existing Python test compares it line by line.
                for (const std::int64_t bin : hicx::correct::cpython_set_order(outliers)) {
                    const hicx::CutInterval& interval =
                        data.cut_intervals[static_cast<std::size_t>(bin)];
                    bed << interval.chrom << '\t' << interval.start << '\t' << interval.end
                        << "\t.\t" << hicx::npy::float_repr(interval.extra) << "\t.\n";
                }
            }

            if (!outliers.empty()) {
                std::vector<char> masked(static_cast<std::size_t>(data.matrix.rows()), 0);
                for (const std::int64_t bin : outliers) {
                    masked[static_cast<std::size_t>(bin)] = 1;
                }
                for (const std::int64_t bin : data.nan_bins) {
                    masked[static_cast<std::size_t>(bin)] = 1;
                }
                apply_mask(data, chain, masked);
                hic.refresh_boundaries();
            }

            if (args.sequenced_count_cutoff.has_value() && *args.sequenced_count_cutoff != 0.0 &&
                *args.sequenced_count_cutoff > 0.0 && *args.sequenced_count_cutoff < 1.0) {
                if (!ends_with(args.matrix, ".h5")) {
                    // hicCorrectMatrix.py:678. A cool file's cut intervals
                    // carry the Python float 1.0, not a np.float64, so the
                    // assertion fails before anything is filtered.
                    std::fputs("hicCorrectMatrix: AssertionError: "
                               "hicCorrectMatrix.py:678 asserts that the per bin coverage is "
                               "a numpy.float64. A cool file's cut intervals carry the "
                               "Python float 1.0, so --sequencedCountCutoff fails on every "
                               "cool input.\n",
                               stderr);
                    return 1;
                }
                std::vector<char> masked(static_cast<std::size_t>(data.matrix.rows()), 0);
                std::size_t failed = 0;
                for (std::size_t bin = 0; bin < data.cut_intervals.size(); ++bin) {
                    if (data.cut_intervals[bin].extra < *args.sequenced_count_cutoff) {
                        masked[bin] = 1;
                        ++failed;
                    }
                }
                log_info(args.verbose,
                         "Bins with low coverage: " + std::to_string(failed));
                if (failed > 0) {
                    for (const std::int64_t bin : data.nan_bins) {
                        masked[static_cast<std::size_t>(bin)] = 1;
                    }
                    apply_mask(data, chain, masked);
                    hic.refresh_boundaries();
                }
            }

            if (args.trans_cutoff.has_value() && *args.trans_cutoff > 0.0 &&
                *args.trans_cutoff < 100.0) {
                // ma.truncTrans is a no operation, finding F11: it unpacks a
                // 3-tuple into 2 names and then compares where it meant to
                // assign. Reproduced as the no operation it is.
                log_info(args.verbose, "truncTrans is a no operation in the reference");
            }
        }

        if (args.inflation_cutoff.has_value() && *args.inflation_cutoff > 0.0 &&
            args.method == Method::Ice) {
            if (!args.trans_cutoff.has_value() || !(*args.trans_cutoff > 0.0) ||
                !(*args.trans_cutoff < 100.0)) {
                std::fputs("hicCorrectMatrix: UnboundLocalError: cannot access local "
                           "variable 'pre_row_sum' where it is not associated with a value. "
                           "hicCorrectMatrix.py:699 assigns pre_row_sum only inside the "
                           "--transCutoff branch, and :770 reads it for --inflationCutoff, "
                           "so --inflationCutoff without --transCutoff cannot work.\n",
                           stderr);
                return 1;
            }
            std::fputs("hicCorrectMatrix: IndexError: list index out of range. "
                       "hicCorrectMatrix.py:776 passes the union of the MAD outliers and the "
                       "inflated bins to printchrtoremove, whose ids index the matrix as it "
                       "was before the masking, so the lookup runs past the end of the "
                       "shortened cut_intervals. --inflationCutoff cannot complete.\n",
                       stderr);
            return 1;
        }

        std::vector<double> correction_factors;
        bool factors_are_column = false;
        std::string failure;

        if (args.perchr) {
            const std::vector<std::pair<std::string, hicx::BinRange>> boundaries =
                hic.boundaries();
            std::vector<double> raw_values;
            const bool needs_raw_matrix = !output_is_h5;
            if (args.method == Method::Kr && output_is_h5) {
                materialise_kr_matrix(data.matrix, 0.00001, args.compat_v3);
            } else if (args.method == Method::Ice && needs_raw_matrix) {
                // --perchr ICE assembles its result into a lil_matrix that is
                // never read back for a non '.h5' name, so the raw values have
                // to survive the in place correction.
                raw_values = data.matrix.data();
            } else if (args.method == Method::Kr && args.compat_v3 && needs_raw_matrix) {
                raw_values = data.matrix.data();
            }

            correction_factors.reserve(static_cast<std::size_t>(data.matrix.rows()));
            for (const auto& entry : boundaries) {
                hicx::kernels::DiagonalBlock block(data.matrix, entry.second.first,
                                                   entry.second.last);
                if (args.method == Method::Ice) {
                    hicx::ice::Options options;
                    options.max_iterations = args.iter_num;
                    options.threads = threads;
                    const hicx::ice::Result result = hicx::ice::correct(block, options);
                    if (result.failed) {
                        std::fprintf(stderr, "ERROR:hicexplorer.iterativeCorrection:%s\n",
                                     result.message.c_str());
                        return 1;
                    }
                    if (result.converged) {
                        log_info(args.verbose, "[iterative correction] " +
                                                   std::to_string(result.iterations + 1) +
                                                   " iterations used");
                    }
                    log_info(args.verbose, entry.first + ": [iterative correction] passes " +
                                               std::to_string(result.iterations));
                    correction_factors.insert(correction_factors.end(),
                                              result.correction_factors.begin(),
                                              result.correction_factors.end());
                } else {
                    hicx::kr::Options options;
                    options.threads = threads;
                    options.float32_input = args.compat_v3 && !output_is_h5;
                    options.float32_rescale = args.compat_v3;
                    options.diagonal_addend = output_is_h5 ? 0.0 : 0.00001;
                    hicx::kr::Balancer balancer(block, options);
                    if (!balancer.compute()) {
                        std::fprintf(stderr, "hicCorrectMatrix: %s\n",
                                     balancer.message().c_str());
                        return 1;
                    }
                    // hicCorrectMatrix.py:726-732. The matrix is only asked for
                    // when the output is '.h5', and that request is what
                    // rescales the vector; the vector is then read with
                    // rescale=False either way.
                    if (output_is_h5) {
                        balancer.normalise_matrix(true);
                        log_info(args.verbose,
                                 entry.first + ": normalisation factor is " +
                                     hicx::npy::float_repr(balancer.normalisation_factor()));
                    }
                    const std::vector<double>& vector = balancer.normalisation_vector(false);
                    correction_factors.insert(correction_factors.end(), vector.begin(),
                                              vector.end());
                    factors_are_column = true;
                }
            }
            if (output_is_h5) {
                hicx::correct::keep_only_diagonal_blocks(data.matrix, boundaries);
            } else if (!raw_values.empty()) {
                data.matrix.mutable_data() = std::move(raw_values);
            }
        } else {
            if (args.method == Method::Ice) {
                hicx::kernels::DiagonalBlock block(data.matrix);
                hicx::ice::Options options;
                options.max_iterations = args.iter_num;
                options.threads = threads;
                const hicx::ice::Result result = hicx::ice::correct(block, options);
                if (result.failed) {
                    std::fprintf(stderr, "ERROR:hicexplorer.iterativeCorrection:%s\n",
                                 result.message.c_str());
                    return 1;
                }
                if (result.converged) {
                    // iterativeCorrection.py:73 logs iternum + 1, an off by one
                    // in the message, and logs nothing when the loop runs to
                    // its iteration limit.
                    log_info(args.verbose, "[iterative correction] " +
                                               std::to_string(result.iterations + 1) +
                                               " iterations used");
                }
                // Not a Python line. The harness and cpp/PLAN.md 5.3 want the
                // pass count compared between the two implementations, and the
                // reference only prints it on convergence.
                log_info(args.verbose,
                         "[iterative correction] passes " + std::to_string(result.iterations));
                correction_factors = result.correction_factors;
            } else {
                std::vector<double> raw_values;
                if (args.compat_v3 && !output_is_h5) {
                    raw_values = data.matrix.data();
                }
                if (output_is_h5) {
                    materialise_kr_matrix(data.matrix, 0.00001, args.compat_v3);
                }
                hicx::kernels::DiagonalBlock kr_block(data.matrix);
                hicx::kr::Options options;
                options.threads = threads;
                options.float32_input = args.compat_v3 && !output_is_h5;
                options.float32_rescale = args.compat_v3;
                options.diagonal_addend = output_is_h5 ? 0.0 : 0.00001;
                hicx::kr::Balancer balancer(kr_block, options);
                if (!balancer.compute()) {
                    std::fprintf(stderr, "hicCorrectMatrix: %s\n", balancer.message().c_str());
                    return 1;
                }
                // hicCorrectMatrix.py:754 asks for the vector with rescale=True
                // first, so the factors are rescaled for both output formats
                // here, unlike the --perchr path.
                correction_factors = balancer.normalisation_vector(true);
                log_info(args.verbose, "normalisation factor is " +
                                           hicx::npy::float_repr(
                                               balancer.normalisation_factor()));
                if (output_is_h5) {
                    balancer.normalise_matrix(true);
                } else if (!raw_values.empty()) {
                    data.matrix.mutable_data() = std::move(raw_values);
                }
                factors_are_column = true;
            }
        }

        data.correction_factors = std::move(correction_factors);
        data.correction_factors_are_column = factors_are_column;
        restore_mask(data, chain);

        if (!hic.save(args.out_file_name)) {
            std::fprintf(stderr,
                         "hicCorrectMatrix: --outFileName '%s' ends in neither 'cool' nor "
                         "'h5'; hiCMatrix.save writes nothing in that case\n",
                         args.out_file_name.c_str());
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicCorrectMatrix: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicCorrectMatrix");
    return 0;
}
