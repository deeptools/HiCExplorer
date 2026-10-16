// hicArrowhead: contact domain calling with Juicer/Juicebox's Arrowhead
// algorithm (Rao et al. 2014, Cell, doi: 10.1016/j.cell.2014.11.021).
// HiCExplorer has no Arrowhead caller of its own (hicFindTADs uses an
// independent insulation-score method, PLAN.md tier 5), so this is a
// faithful reimplementation of the real tool, not an invented stand-in --
// see tools/arrowhead_impl.hpp for the full algorithm notes, references to
// the real Java source read for this port, and its documented deviations
// (normalization is whatever the input matrix already carries rather than a
// runtime NONE/VC/VC_SQRT/KR selector; connected components use a plain
// flood fill; the optional feature/control list scoring is not ported; the
// output format is trimmed to the fields that carry real information; one
// upstream edge-case window-sizing inconsistency at a chromosome's last
// window is not reproduced).
//
// Output format (new, no Python precedent; a simplified BEDPE, not
// byte-identical to Juicer's own "*_blocks.bedpe"): tab separated,
//
//   chrom  x1  x2  score  uVarScore  lVarScore  upSign  loSign
//
// Domains are square blocks on the diagonal (x1, x2 apply to both axes), so
// there is one coordinate pair, not two.

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "arrowhead_impl.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;

const char* const kUsage =
    "usage: hicArrowhead --matrix MATRIX --outFileName OUTFILENAME\n"
    "                    [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                    [--matrixWidth MATRIXWIDTH] [--varThreshold VARTHRESHOLD]\n"
    "                    [--highSignThreshold HIGHSIGNTHRESHOLD]\n"
    "                    [--maxLowSignThreshold MAXLOWSIGNTHRESHOLD]\n"
    "                    [--minLowSignThreshold MINLOWSIGNTHRESHOLD]\n"
    "                    [--decrementLowSignThreshold DECREMENTLOWSIGNTHRESHOLD]\n"
    "                    [--minBlockSize MINBLOCKSIZE] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "Calls contact domains on the given contact matrix with Juicer/Juicebox's\n"
    "Arrowhead algorithm (Rao et al. 2014). New in HiCExplorer v4 (PLAN.md tier\n"
    "9): a faithful reimplementation of the real tool's own method, not an\n"
    "invented stand-in -- see the source file header for the full algorithm\n"
    "notes and its documented deviations (normalization, connected components,\n"
    "the feature/control list scoring option, the output format).\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        The matrix to call contact domains on.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        Outfile name to store the called domains (tab\n"
    "                        separated): chrom, x1, x2, score, uVarScore,\n"
    "                        lVarScore, upSign, loSign.\n"
    "\n"
    "Optional arguments:\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes to include in the analysis. If not set,\n"
    "                        all chromosomes are included.\n"
    "  --matrixWidth MATRIXWIDTH\n"
    "                        Size of the sliding window along the diagonal, in\n"
    "                        bins (Arrowhead's -m). (Default: 2000).\n"
    "  --varThreshold VARTHRESHOLD\n"
    "                        Variance-sum threshold for the high-confidence\n"
    "                        pass. (Default: 0.2).\n"
    "  --highSignThreshold HIGHSIGNTHRESHOLD\n"
    "                        Sign threshold for the high-confidence pass.\n"
    "                        (Default: 0.5).\n"
    "  --maxLowSignThreshold MAXLOWSIGNTHRESHOLD\n"
    "                        Starting sign threshold for the low-confidence\n"
    "                        sweep. (Default: 0.4).\n"
    "  --minLowSignThreshold MINLOWSIGNTHRESHOLD\n"
    "                        Lowest sign threshold the low-confidence sweep\n"
    "                        will try. (Default: 0).\n"
    "  --decrementLowSignThreshold DECREMENTLOWSIGNTHRESHOLD\n"
    "                        Step size of the low-confidence sweep.\n"
    "                        (Default: 0.1).\n"
    "  --minBlockSize MINBLOCKSIZE\n"
    "                        Minimum width, in bins, for a low-confidence-only\n"
    "                        block to be kept. (Default: 60).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n"
    "\n"
    "C++ only, no Python counterpart (PLAN.md tier 9).\n";

class PythonError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

bool contains(const std::vector<std::string>& list, const std::string& value) {
    return std::find(list.begin(), list.end(), value) != list.end();
}

// A dense, symmetric window [lim_start, lim_end) of one chromosome's
// stored matrix, local bin indices offset by the chromosome's own first
// bin: Juicer's extractLocalBoundedRegion(..., fillUnderDiagonal=false)
// followed by MatrixTools.fillLowerLeftTriangle(observed).
hicx::arrowhead::Matrix dense_window(const hicx::CsrMatrix& csr, std::int64_t chrom_first,
                                     std::int64_t lim_start, std::int64_t lim_end) {
    const int n = static_cast<int>(lim_end - lim_start);
    hicx::arrowhead::Matrix out(n);
    const auto& indptr = csr.indptr();
    const auto& indices = csr.indices();
    const auto& data = csr.data();
    for (std::int64_t row = lim_start; row < lim_end; ++row) {
        const std::int64_t grow = chrom_first + row;
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(grow)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(grow) + 1]); ++k) {
            const std::int64_t col = indices[k] - chrom_first;
            if (col < lim_start || col >= lim_end || data[k] == 0.0) {
                continue;
            }
            const int i = static_cast<int>(row - lim_start);
            const int j = static_cast<int>(col - lim_start);
            out.at(i, j) = data[k];
            out.at(j, i) = data[k];
        }
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("hicArrowhead",
                       "Calls contact domains with Juicer/Juicebox's Arrowhead algorithm.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"}).required().input({"h5", "cool"}).help(
        "The matrix to call contact domains on.");
    required.add({"--outFileName", "-o"})
        .required()
        .output({"bedpe", "txt", "tsv"})
        .help("Outfile name to store the called domains.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("Chromosomes to include in the analysis.");
    optional.add({"--matrixWidth"})
        .type("int")
        .default_value(std::int64_t{2000})
        .help("Size of the sliding window along the diagonal, in bins (Default: %(default)s).");
    optional.add({"--varThreshold"})
        .type("float")
        .default_value(0.2)
        .help("Variance-sum threshold for the high-confidence pass (Default: %(default)s).");
    optional.add({"--highSignThreshold"})
        .type("float")
        .default_value(0.5)
        .help("Sign threshold for the high-confidence pass (Default: %(default)s).");
    optional.add({"--maxLowSignThreshold"})
        .type("float")
        .default_value(0.4)
        .help("Starting sign threshold for the low-confidence sweep (Default: %(default)s).");
    optional.add({"--minLowSignThreshold"})
        .type("float")
        .default_value(0.0)
        .help("Lowest sign threshold the low-confidence sweep will try (Default: %(default)s).");
    optional.add({"--decrementLowSignThreshold"})
        .type("float")
        .default_value(0.1)
        .help("Step size of the low-confidence sweep (Default: %(default)s).");
    optional.add({"--minBlockSize"})
        .type("int")
        .default_value(std::int64_t{60})
        .help("Minimum width, in bins, for a low-confidence-only block to be kept (Default: "
              "%(default)s).");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);

    try {
        const std::string matrix_path = ns.str("matrix");
        const std::optional<std::vector<std::string>> chromosomes =
            ns.given("chromosomes") ? std::optional(ns.strs("chromosomes")) : std::nullopt;

        hicx::arrowhead::Options options;
        options.matrix_width = static_cast<int>(ns.integer("matrixWidth"));
        options.var_threshold = ns.real("varThreshold");
        options.high_sign_threshold = ns.real("highSignThreshold");
        options.max_low_sign_threshold = ns.real("maxLowSignThreshold");
        options.min_low_sign_threshold = ns.real("minLowSignThreshold");
        options.decrement_low_sign_threshold = ns.real("decrementLowSignThreshold");
        options.min_block_size = static_cast<int>(ns.integer("minBlockSize"));
        if (options.matrix_width < 2) {
            throw PythonError("ValueError: --matrixWidth must be at least 2");
        }

        const hicx::ToolMatrix matrix = hicx::ToolMatrix::load(matrix_path);
        const hicx::CsrMatrix& csr = matrix.matrix();

        std::ofstream out(ns.str("outFileName"));
        if (!out) {
            throw PythonError("IOError: cannot write to " + ns.str("outFileName"));
        }
        out << "#chrom\tx1\tx2\tscore\tuVarScore\tlVarScore\tupSign\tloSign\n";

        std::size_t total_domains = 0;
        for (const auto& [chrom, range] : matrix.boundaries()) {
            if (chromosomes.has_value() && !contains(*chromosomes, chrom)) {
                continue;
            }
            const std::int64_t chrom_first = range.first;
            const std::int64_t chrom_length = range.last - range.first;
            if (chrom_length < 2) {
                continue;
            }
            const std::vector<hicx::arrowhead::HighScore> domains = hicx::arrowhead::call_domains(
                static_cast<int>(chrom_length), options,
                [&](int lim_start, int lim_end) {
                    return dense_window(csr, chrom_first, lim_start, lim_end);
                });
            for (const hicx::arrowhead::HighScore& h : domains) {
                const std::int64_t x1 = matrix.cut_intervals()[static_cast<std::size_t>(chrom_first + h.i)].start;
                const std::int64_t x2 = matrix.cut_intervals()[static_cast<std::size_t>(chrom_first + h.j)].start;
                out << chrom << '\t' << x1 << '\t' << x2 << '\t' << h.score << '\t' << h.u_var_score
                    << '\t' << h.l_var_score << '\t' << h.up_sign << '\t' << h.lo_sign << '\n';
                ++total_domains;
            }
        }
        out.close();

        hicx::report_resource_usage("hicArrowhead");
        std::fprintf(stderr, "hicArrowhead: %zu domains written to %s\n", total_domains,
                    ns.str("outFileName").c_str());
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicArrowhead: %s\n", error.what());
        return 1;
    }
}
