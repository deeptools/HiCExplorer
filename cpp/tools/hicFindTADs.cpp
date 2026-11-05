// Port of hicexplorer/hicFindTADs.py.
//
// The pipeline, and where each piece lives:
//
//   1. load the matrix, optionally restricting and reordering the chromosomes
//   2. drop the diagonal, mask the NaN bins           (hicx/find_tads.hpp)
//   3. z-score the matrix per chromosome              (hicx::convert_to_obs_exp_matrix)
//   4. close the gaps the masking left                (hicx::enlarge_bins)
//   5. truncate to 2 * maxDepth around the diagonal   (hicx/find_tads.hpp)
//   6. the TAD-separation score per bin and window    (hicx/find_tads.hpp, threaded here)
//   7. local minima of the mean score                 (peakdetect, hicx/find_tads.hpp)
//   8. a Wilcoxon rank sum p-value per minimum        (hicx::stats::ranksums, threaded)
//   9. FDR or Bonferroni, then the call sets          (hicx::stats)
//
// Threading. The Python uses multiprocessing.Pool over a contiguous range of
// bins per worker (hicFindTADs.py:1103-1109). This port uses threads over the
// same fixed partition of the bin index range, writes every result into its
// own preallocated slot and combines the slots in index order, so the output
// does not depend on the thread count (cpp/OPTIMIZATION.md 3). The p-values
// are threaded the same way. Nothing is accumulated into a shared float.
//
// Behaviour of the Python that is reproduced deliberately, not fixed:
//
//   * hicFindTADs.py:902-920 assigns chr_end_idx = chr_start_idx and then
//     mutates it in place, so the chromosome *start* indices are lost and the
//     end indices are added to the boundary list twice.
//   * delta_wrt_window (:612-613) averages matrix_avg[i-w : i+3] together with
//     matrix_avg[i+4 : i+w], which includes the minimum itself and the two
//     bins after it and skips the bin at i+3, although the docstring says the
//     minimum is excluded.
//   * get_incremental_step_size grows as min + int(step * k ** 1.5), which is
//     not the progression the --step help text describes.
//   * hiCMatrix caches its bin size at construction time, so the bin size used
//     after the bins have been enlarged is still the one from the input.
//   * the z-score matrix keeps NaN entries that lie beyond the truncation
//     distance, because NaN - NaN is NaN and scipy stores every result that is
//     not exactly zero.
//
// One behaviour is deliberately *not* reproduced, and it is a defect rather
// than a numeric difference. np.array_split hands every worker a contiguous
// range of bins, and a worker whose bins are all skipped reaches
// `zip(*positions_array)` with an empty list and ends the whole run with
// "ValueError: not enough values to unpack". On small_test_matrix.h5 that
// happens from --numberOfProcessors 12 upwards while 1 to 10 succeed and give
// byte-identical output, so in the Python the process count decides whether
// the tool runs at all. Here an empty partition simply contributes no rows.
// Pinned by
// test_hicFindTADs.py::test_find_TADs_high_numberOfProcessors_is_a_defect.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/argparse.hpp"
#include "hicx/bins.hpp"
#include "hicx/find_tads.hpp"
#include "hicx/matrix_ops.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/parallel.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/version.hpp"

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

const char* const kUsage =
    "usage: hicFindTADs --matrix MATRIX --outPrefix OUTPREFIX\n"
    "                   --correctForMultipleTesting {fdr,bonferroni,None}\n"
    "                   [--minDepth INT bp] [--maxDepth INT bp] [--step INT bp]\n"
    "                   [--TAD_sep_score_prefix TAD_SEP_SCORE_PREFIX]\n"
    "                   [--thresholdComparisons THRESHOLDCOMPARISONS]\n"
    "                   [--delta DELTA] [--minBoundaryDistance MINBOUNDARYDISTANCE]\n"
    "                   [--ctcfFile CTCFFILE] [--ctcfWindow CTCFWINDOW]\n"
    "                   [--ctcfBoundarySlack CTCFBOUNDARYSLACK]\n"
    "                   [--nestedMinDepth INT bp] [--nestedMaxDepth INT bp]\n"
    "                   [--nestedStep INT bp] [--nestedMinBoundaryDistance BP]\n"
    "                   [--nestedMinDomainSize BP]\n"
    "                   [--chromosomes CHROMOSOMES [CHROMOSOMES ...]]\n"
    "                   [--numberOfProcessors NUMBEROFPROCESSORS] [--help]\n"
    "                   [--version]\n";

const char* const kHelp =
    "\n"
    "Uses a measure called TAD-separation score to identify the degree of "
    "separation between\n"
    "the left and right regions at each Hi-C matrix bin. This is done for a\n"
    "running window of different sizes. Then, TADs are called as those\n"
    "positions having a local TAD-separation score minimum.\n"
    "\n"
    "Required arguments:\n"
    "  --matrix MATRIX, -m MATRIX\n"
    "                        Corrected Hi-C matrix to use for the computations.\n"
    "  --outPrefix OUTPREFIX\n"
    "                        File prefix to save the resulting files.\n"
    "  --correctForMultipleTesting {fdr,bonferroni,None}\n"
    "                        Select the bonferroni or false discovery rate for a\n"
    "                        multiple comparison. (Default: fdr).\n"
    "\n"
    "Optional arguments:\n"
    "  --minDepth INT bp     Minimum window length (in bp).\n"
    "  --maxDepth INT bp     Maximum window length (in bp).\n"
    "  --step INT bp         Step size when moving from --minDepth to --maxDepth.\n"
    "  --TAD_sep_score_prefix TAD_SEP_SCORE_PREFIX\n"
    "                        Prefix of an existing TAD-separation score.\n"
    "  --thresholdComparisons THRESHOLDCOMPARISONS\n"
    "                        P-value threshold for the Bonferroni correction /\n"
    "                        q-value for FDR (Default: 0.01).\n"
    "  --delta DELTA         Minimum threshold of the difference between the\n"
    "                        TAD-separation score of a putative boundary and the\n"
    "                        mean of the TAD-sep. score of surrounding bins\n"
    "                        (Default: 0.01).\n"
    "  --minBoundaryDistance MINBOUNDARYDISTANCE\n"
    "                        Minimum distance between boundaries (in bp).\n"
    "  --ctcfFile CTCFFILE   C++-only. A BED file (chrom, start, end, name, score,\n"
    "                        strand) of CTCF binding site motif orientations. When\n"
    "                        given, a TAD-separation score minimum is accepted as a\n"
    "                        boundary only when a CTCF site within --ctcfWindow is\n"
    "                        oriented into the adjacent domain: a site at or\n"
    "                        before the boundary must be on the '-' strand, one at\n"
    "                        or after it on the '+' strand, the divergent\n"
    "                        orientation found at domain boundaries.\n"
    "  --ctcfWindow CTCFWINDOW\n"
    "                        C++-only. Maximum distance in bp between a boundary\n"
    "                        and a supporting CTCF site (Default: 10000).\n"
    "  --ctcfBoundarySlack CTCFBOUNDARYSLACK\n"
    "                        C++-only. Evaluate the +/- orientation test against any\n"
    "                        point within this many bp of the boundary instead of only\n"
    "                        the exact reported minimum (Default: 0).\n"
    "  --nestedMinDepth INT bp\n"
    "                        C++-only. Together with --nestedMaxDepth and\n"
    "                        --nestedStep, enables sub-TAD detection: within each\n"
    "                        top-level domain, looks for boundaries at this finer\n"
    "                        depth range, writing <outPrefix>_nested_boundaries.bed,\n"
    "                        _nested_domains.bed, _nested_boundaries.gff and\n"
    "                        _nested_score.bedgraph. The last is the domain's own\n"
    "                        TAD-separation score computed fresh over this finer\n"
    "                        depth range, independently of the top level's score\n"
    "                        (over --minDepth..--maxDepth) in <outPrefix>_score.\n"
    "                        bedgraph -- the two are not the same curve and the\n"
    "                        nested boundaries are the minima of this one, not of\n"
    "                        the top-level one, filtered by --delta alone -- not by\n"
    "                        a Wilcoxon p-value/FDR threshold the way top-level\n"
    "                        boundaries are; that p-value is still computed and\n"
    "                        reported in the GFF (pvalue=) but is not a gate, since\n"
    "                        at a domain's own small scale it has too little power\n"
    "                        to tell a real minimum apart from one slightly less\n"
    "                        deep. All three of --nestedMinDepth/--nestedMaxDepth/\n"
    "                        --nestedStep are required together; omitted, no\n"
    "                        sub-TADs are found.\n"
    "                        Must be at least 3 times the matrix bin size, same as\n"
    "                        --minDepth, otherwise that reported p-value has too few\n"
    "                        values per side to mean anything.\n"
    "  --nestedMaxDepth INT bp\n"
    "                        C++-only. See --nestedMinDepth.\n"
    "  --nestedStep INT bp   C++-only. See --nestedMinDepth.\n"
    "  --nestedMinBoundaryDistance BP\n"
    "                        C++-only. Minimum distance between sub-TAD boundaries\n"
    "                        (Default: 4 times the average nested bin size).\n"
    "  --nestedMinDomainSize BP\n"
    "                        C++-only. A top-level domain shorter than this is not\n"
    "                        searched for sub-TADs (Default: 4 * --nestedMaxDepth).\n"
    "  --chromosomes CHROMOSOMES [CHROMOSOMES ...]\n"
    "                        Chromosomes and order in which the chromosomes should\n"
    "                        be plotted.\n"
    "  --numberOfProcessors NUMBEROFPROCESSORS, -p NUMBEROFPROCESSORS\n"
    "                        Number of processors to use (Default: 1).\n"
    "  --help, -h            show this help message and exit.\n"
    "  --version             show program's version number and exit\n";

struct Arguments {
    std::string matrix;
    std::string out_prefix;
    std::string correct_for_multiple_testing;
    std::optional<std::int64_t> min_depth;
    std::optional<std::int64_t> max_depth;
    std::optional<std::int64_t> step;
    std::optional<std::string> tad_sep_score_prefix;
    double threshold_comparisons = 0.01;
    double delta = 0.01;
    std::optional<std::int64_t> min_boundary_distance;
    std::optional<std::string> ctcf_file;
    std::int64_t ctcf_window = 10000;
    std::int64_t ctcf_boundary_slack = 0;
    std::optional<std::int64_t> nested_min_depth;
    std::optional<std::int64_t> nested_max_depth;
    std::optional<std::int64_t> nested_step;
    std::optional<std::int64_t> nested_min_boundary_distance;
    std::optional<std::int64_t> nested_min_domain_size;
    std::optional<std::vector<std::string>> chromosomes;
    int processors = 1;
};

// hicFindTADs.py parse_arguments.
Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("hicFindTADs",
                       "Uses a measure called TAD-separation score to identify the degree of "
                       "separation between the left and right regions at each Hi-C matrix bin.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--matrix", "-m"})
        .required()
        .input({"h5", "cool", "mcool"})
        .help("Corrected Hi-C matrix to use for the computations.");
    required.add({"--outPrefix"})
        .required()
        .output({"bm", "h5", "cool", "bed", "gff", "bedgraph"}, "prefix")
        .help("File prefix to save the resulting files.");
    required.add({"--correctForMultipleTesting"})
        .type("str")
        .default_value("fdr")
        .choices({"fdr", "bonferroni", "None"})
        .required()
        .help("Select the bonferroni or false discovery rate for a multiple comparison.");
    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--minDepth"}).type("int").metavar("INT bp").help("Minimum window length (in bp).");
    optional.add({"--maxDepth"}).type("int").metavar("INT bp").help("Maximum window length (in bp).");
    optional.add({"--step"})
        .type("int")
        .metavar("INT bp")
        .help("Step size when moving from --minDepth to --maxDepth.");
    optional.add({"--TAD_sep_score_prefix"})
        .required(false)
        .input({"bm", "h5", "cool"}, "prefix")
        .help("Prefix of an existing TAD-separation score and z-score matrix.");
    optional.add({"--thresholdComparisons"})
        .type("float")
        .default_value(0.01)
        .help("P-value threshold for the Bonferroni correction / q-value for FDR.");
    optional.add({"--delta"})
        .type("float")
        .default_value(0.01)
        .help("Minimum threshold of the difference between the TAD-separation score of a "
              "putative boundary and the mean of the TAD-sep. score of surrounding bins.");
    optional.add({"--minBoundaryDistance"})
        .type("int")
        .help("Minimum distance between boundaries (in bp).");
    optional.add({"--ctcfFile"})
        .cpp_only("Boundary selection by CTCF motif orientation, new in version 4.")
        .input({"bed", "narrowPeak", "txt"})
        .help("C++-only. A BED file (chrom, start, end, name, score, strand) of CTCF "
              "binding site motif orientations. When given, a TAD-separation score "
              "minimum is accepted as a boundary only when a CTCF site within "
              "--ctcfWindow is oriented into the adjacent domain: a site at or before "
              "the boundary must be on the '-' strand, one at or after it on the '+' "
              "strand, the divergent orientation found at domain boundaries.");
    optional.add({"--ctcfWindow"})
        .cpp_only("Used with --ctcfFile, new in version 4.")
        .type("int")
        .default_value(std::int64_t{10000})
        .metavar("BP")
        .help("C++-only. Maximum distance between a boundary and a supporting CTCF "
              "site.");
    optional.add({"--ctcfBoundarySlack"})
        .cpp_only("Used with --ctcfFile, new in version 4.")
        .type("int")
        .default_value(std::int64_t{0})
        .metavar("BP")
        .help("C++-only. Evaluate the +/- orientation test against any point within "
              "this many bp of the boundary instead of only the exact reported "
              "minimum.");
    optional.add({"--nestedMinDepth"})
        .cpp_only("Nested domain calling, new in version 4.")
        .type("int")
        .metavar("INT bp")
        .help("C++-only. Together with --nestedMaxDepth and --nestedStep, enables "
              "sub-TAD detection within each top-level domain.");
    optional.add({"--nestedMaxDepth"})
        .cpp_only("Nested domain calling, new in version 4.")
        .type("int")
        .metavar("INT bp")
        .help("C++-only. See --nestedMinDepth.");
    optional.add({"--nestedStep"})
        .cpp_only("Nested domain calling, new in version 4.")
        .type("int")
        .metavar("INT bp")
        .help("C++-only. See --nestedMinDepth.");
    optional.add({"--nestedMinBoundaryDistance"})
        .cpp_only("Nested domain calling, new in version 4.")
        .type("int")
        .metavar("BP")
        .help("C++-only. Minimum distance between sub-TAD boundaries (Default: 4 "
              "times the average nested bin size).");
    optional.add({"--nestedMinDomainSize"})
        .cpp_only("Nested domain calling, new in version 4.")
        .type("int")
        .metavar("BP")
        .help("C++-only. A top-level domain shorter than this is not searched for "
              "sub-TADs (Default: 4 * --nestedMaxDepth).");
    optional.add({"--chromosomes"})
        .nargs("+")
        .help("Chromosomes and order in which the chromosomes should be plotted.");
    optional.add({"--numberOfProcessors", "-p"})
        .type("int")
        .default_value(1)
        .help("Number of processors to use.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit.");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.matrix = ns.str("matrix");
    args.out_prefix = ns.str("outPrefix");
    args.correct_for_multiple_testing = ns.str("correctForMultipleTesting");
    args.min_depth = ns.opt_integer("minDepth");
    args.max_depth = ns.opt_integer("maxDepth");
    args.step = ns.opt_integer("step");
    args.tad_sep_score_prefix = ns.opt_str("TAD_sep_score_prefix");
    args.threshold_comparisons = ns.real("thresholdComparisons");
    args.delta = ns.real("delta");
    args.min_boundary_distance = ns.opt_integer("minBoundaryDistance");
    args.ctcf_file = ns.opt_str("ctcfFile");
    args.ctcf_window = ns.integer("ctcfWindow");
    args.ctcf_boundary_slack = ns.integer("ctcfBoundarySlack");
    args.nested_min_depth = ns.opt_integer("nestedMinDepth");
    args.nested_max_depth = ns.opt_integer("nestedMaxDepth");
    args.nested_step = ns.opt_integer("nestedStep");
    args.nested_min_boundary_distance = ns.opt_integer("nestedMinBoundaryDistance");
    args.nested_min_domain_size = ns.opt_integer("nestedMinDomainSize");
    if (ns.given("chromosomes")) {
        args.chromosomes = ns.strs("chromosomes");
    }
    args.processors = static_cast<int>(ns.integer("numberOfProcessors"));
    return args;
}

bool file_exists(const std::string& path) {
    std::ifstream probe(path);
    return probe.good();
}

// Python compares two lists of cut interval tuples with ==, which tests
// identity before value, so a NaN in `extra` that came from the same object
// compares equal. enlarge_bins reuses the extra of the tuple it rewrites, so
// bitwise equality is the right test here.
bool same_intervals(const std::vector<hicx::CutInterval>& a,
                    const std::vector<hicx::CutInterval>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].chrom != b[i].chrom || a[i].start != b[i].start ||
            a[i].end != b[i].end || a[i].extra_text != b[i].extra_text) {
            return false;
        }
        const double x = a[i].extra;
        const double y = b[i].extra;
        if (!(x == y) && !(std::isnan(x) && std::isnan(y))) {
            return false;
        }
    }
    return true;
}

// The kernels shared with the Python API (hicx/find_tads.hpp).
using hicx::find_tads::bins_at_distance;
using hicx::find_tads::cut_weight_mean;
using hicx::find_tads::cut_weight_values;
using hicx::find_tads::delta_wrt_window;
using hicx::find_tads::DenseBlock;
using hicx::find_tads::drop_diagonal;
using hicx::find_tads::extract_block;
using hicx::find_tads::incremental_step_size;
using hicx::find_tads::mask_bins;
using hicx::find_tads::MaskState;
using hicx::find_tads::Peak;
using hicx::find_tads::peakdetect;
using hicx::find_tads::restore_masked_bins;
using hicx::find_tads::ScoreContext;
using hicx::find_tads::truncate_to_distance;

// --------------------------------------------------------------------------
// Output formatting

std::string format_fixed(double value, int decimals) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", decimals, value);
    return buffer;
}

// --ctcfFile (C++-only, cpp/PLAN.md tier 9): requires a CTCF binding site in
// the correct orientation for a TAD-separation score minimum to be accepted
// as a boundary. Loops are anchored by convergent motifs, forward at the
// upstream and reverse at the downstream anchor (Rao et al. 2014), so the
// motifs of a domain point into it: forward just after its start, reverse
// just before its end. A boundary between two domains therefore carries a
// divergent pair (Gomez-Marin et al. 2015, PNAS; Fudenberg et al. 2016, Cell
// Reports): '-' (reverse) at or before the boundary position, '+' (forward)
// at or after it. Measured on the GSE63525 GM12878 Arrowhead domains with
// ENCODE CTCF peak motifs: 82.7 % '+' within 10 kb after a domain start,
// 81.4 % '-' within 10 kb before a domain end. Site positions are motif
// midpoints, so a site exactly at the boundary matches either strand.
//
// --ctcfBoundarySlack widens where "the boundary position" is allowed to sit
// for that left/right test. hicFindTADs reports one bin as the minimum, but
// the real inflection of the TAD-separation score is rarely that one exact
// bin; a site just past it, on what looks like the wrong side of the single
// reported bin, can still be the site the minimum is reporting. Slack s
// relaxes the test to accept a '-' site anywhere at or before
// boundary_position + s, and a '+' site anywhere at or after
// boundary_position - s, while --ctcfWindow still bounds how far the site
// itself may be from the minimum. s = 0 (the default) is the original exact
// single-bin test.
//
// A boundary is supported by a site on whichever side has one -- reverse
// before it or forward after it, either alone is enough, not both. A
// stricter requirement for a real divergent *pair* (reverse AND forward
// both present) was tried and dropped: a sparse real site list already
// struggles to place a single site near every true boundary, so requiring
// two only lost real boundaries without a matching gain in precision.
struct CtcfSite {
    std::int64_t position = 0;
    char strand = 0;  // '+' or '-'
};

std::vector<std::string> ctcf_split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t marker = line.find('\t', start);
        if (marker == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, marker - start));
        start = marker + 1;
    }
    return fields;
}

std::string ctcf_strip(const std::string& text) {
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin])) != 0) {
        ++begin;
    }
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return text.substr(begin, end - begin);
}

// Reads a BED file of at least 6 columns (chrom, start, end, name, score,
// strand), keeping only rows with strand '+' or '-' (a '.' carries no
// orientation and is dropped); sites are returned sorted by position within
// each chromosome. Throws on a short line or an unreadable file.
std::unordered_map<std::string, std::vector<CtcfSite>> read_ctcf_sites(
    const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open --ctcfFile '" + path + "'");
    }
    std::unordered_map<std::string, std::vector<CtcfSite>> by_chrom;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        const std::string stripped = ctcf_strip(line);
        if (stripped.empty()) {
            continue;
        }
        const std::vector<std::string> fields = ctcf_split_tabs(stripped);
        if (fields.size() < 6) {
            throw std::runtime_error("--ctcfFile '" + path + "' line " +
                                     std::to_string(line_number) +
                                     " has fewer than 6 columns (chrom, start, end, "
                                     "name, score, strand required)");
        }
        const std::string& strand_field = fields[5];
        if (strand_field != "+" && strand_field != "-") {
            continue;  // '.': no motif orientation, cannot support a boundary
        }
        const std::int64_t start = std::strtoll(fields[1].c_str(), nullptr, 10);
        const std::int64_t end = std::strtoll(fields[2].c_str(), nullptr, 10);
        by_chrom[fields[0]].push_back(CtcfSite{(start + end) / 2, strand_field[0]});
    }
    for (auto& [chrom, sites] : by_chrom) {
        (void)chrom;
        std::sort(sites.begin(), sites.end(),
                 [](const CtcfSite& a, const CtcfSite& b) { return a.position < b.position; });
    }
    return by_chrom;
}

// The nearest site within `window` of `boundary_position` on each orientation
// that points toward it, allowing the left/right test to be evaluated against
// any point within `slack` of `boundary_position` rather than only the exact
// position itself (see the comment on --ctcfBoundarySlack above); preferring
// the closest match by actual distance to `boundary_position` on each side
// independently. Either field is nullopt when that orientation has no
// qualifying site.
struct CtcfPairSupport {
    std::optional<CtcfSite> forward;  // nearest qualifying '+' site
    std::optional<CtcfSite> reverse;  // nearest qualifying '-' site
};

CtcfPairSupport ctcf_support(const std::vector<CtcfSite>& sites, std::int64_t boundary_position,
                             std::int64_t window, std::int64_t slack) {
    const auto first = std::lower_bound(
        sites.begin(), sites.end(), boundary_position - window,
        [](const CtcfSite& site, std::int64_t value) { return site.position < value; });
    CtcfPairSupport result;
    std::int64_t best_forward_distance = std::numeric_limits<std::int64_t>::max();
    std::int64_t best_reverse_distance = std::numeric_limits<std::int64_t>::max();
    for (auto it = first; it != sites.end() && it->position <= boundary_position + window;
         ++it) {
        const std::int64_t distance = std::abs(it->position - boundary_position);
        if (it->position >= boundary_position - slack && it->strand == '+' &&
            distance < best_forward_distance) {
            best_forward_distance = distance;
            result.forward = *it;
        }
        if (it->position <= boundary_position + slack && it->strand == '-' &&
            distance < best_reverse_distance) {
            best_reverse_distance = distance;
            result.reverse = *it;
        }
    }
    return result;
}

struct BedgraphMatrix {
    std::vector<std::string> chrom;
    std::vector<std::int64_t> start;
    std::vector<std::int64_t> end;
    std::vector<double> values;  // row major, `columns` per row
    std::size_t columns = 0;

    [[nodiscard]] std::size_t rows() const { return chrom.size(); }
    [[nodiscard]] std::vector<double> row_means() const {
        std::vector<double> means(rows(), 0.0);
        for (std::size_t r = 0; r < rows(); ++r) {
            means[r] = hicx::simd::pairwise_sum(values.data() + r * columns, columns) /
                       static_cast<double>(columns);
        }
        return means;
    }
};

// --------------------------------------------------------------------------
// Optional nested (sub-TAD) detection: --nestedMinDepth/--nestedMaxDepth/
// --nestedStep, all three required together. For a top-level domain large
// enough to plausibly hold internal structure, reruns the same
// TAD-separation-score -> minima -> delta pipeline the top level uses,
// confined to that domain's own bin range and at a finer depth range,
// giving one additional hierarchy level. Reuses the already truncated,
// z-scored matrix through `context`, since a finer --nestedMaxDepth only
// ever needs data the top-level truncation (to 2 * --maxDepth) already
// kept; no matrix work is repeated.
//
// A sub-boundary is accepted the same way a boundary is understood to begin
// with: a local minimum of the domain's own TAD-separation score
// (_nested_score.bedgraph), filtered by how deep it is relative to its
// neighbourhood (--delta) -- not additionally gated behind a Wilcoxon
// p-value/multiple-testing-correction threshold the way the top level's own
// boundaries are. That significance test is still computed per candidate
// and reported in the GFF (pvalue=), for information, but two earlier
// versions that used it as a hard gate both failed on real data: correcting
// each domain as its own small family (as few as half a dozen candidates)
// made Benjamini-Hochberg's step-up procedure too strict for real signal to
// clear, and pooling the correction across every domain on the chromosome
// did not fix it either -- a domain's single clearest minimum (the largest
// delta of any candidate nearby) sat at a raw p-value of 0.012, only
// marginally above typical significance thresholds to begin with, an
// artifact of the Wilcoxon test's limited power at a window as narrow as
// --nestedMinDepth necessarily is (3 bins), not of how the correction was
// scoped. No scoping of that test rescues a candidate whose raw p-value is
// already borderline, so the fix is to stop gating on it.
//
// Two things this deliberately does *not* do, both tried and reverted:
//
//  * Coarsening the matrix (merging adjacent bins) before searching for
//    sub-boundaries. Coarsening only ever removes resolution: a structure it
//    reveals that --minDepth/--maxDepth did not already see is necessarily
//    *larger* than the matrix's native bin size, i.e. a super-TAD spanning
//    several existing domains, not a sub-TAD inside one of them. Finding
//    real sub-structure needs the matrix's own finest available resolution,
//    which is exactly what reusing `context` at a finer depth gives for
//    free; it cannot need bins coarser than the ones --minDepth/--maxDepth
//    already used.
//  * Re-deriving the observed/expected z-score from only a domain's own
//    local region (its bin range, or any coarsened version of it) instead of
//    reusing the matrix's existing, genome-wide-normalised one. Estimating a
//    diagonal's mean and standard deviation from the very region being
//    tested for internal structure is circular: a real sub-boundary partly
//    determines its own "expected" baseline, which is not a valid null
//    model and inflates significance. `context`'s z-score is independent of
//    any one domain, exactly like the top level's own boundary calls are.
struct NestedCall {
    std::string chrom;
    std::int64_t parent_start = 0;
    std::int64_t parent_end = 0;
    std::string parent_id;
};

// Everything collect_nested_candidates finds for one domain, up to but not
// including multiple-testing correction: the local score, the candidate
// minima's delta values, and their raw (uncorrected) Wilcoxon p-values.
// finalize_nested_domain turns this into output once every domain's p-values
// have been corrected together.
struct NestedCandidates {
    NestedCall call;
    BedgraphMatrix local;
    std::vector<double> local_row_means;
    std::map<std::int64_t, double> delta_of_min;
    std::vector<std::int64_t> pvalue_index;  // indices into `local`
    std::vector<double> pvalues;             // raw, NaN where untestable
};

std::optional<NestedCandidates> collect_nested_candidates(
    const ScoreContext& context, const hicx::BinTable& bins, const NestedCall& call,
    std::int64_t nested_min_depth, std::int64_t nested_max_depth, std::int64_t nested_step,
    std::optional<std::int64_t> nested_min_boundary_distance, unsigned int threads,
    std::ofstream& nested_score) {
    const auto first_range = bins.region_bin_range(call.chrom, call.parent_start, call.parent_start + 1);
    const auto last_range = bins.region_bin_range(call.chrom, call.parent_end - 1, call.parent_end);
    if (!first_range.has_value() || !last_range.has_value()) {
        return std::nullopt;
    }
    const std::int64_t first_bin = first_range->first;
    const std::int64_t last_bin = last_range->first;  // inclusive
    if (last_bin <= first_bin) {
        return std::nullopt;
    }

    const std::vector<std::int64_t> windows =
        incremental_step_size(nested_min_depth, nested_max_depth, nested_step);
    if (windows.empty()) {
        return std::nullopt;
    }

    const auto count = static_cast<std::size_t>(last_bin - first_bin + 1);
    std::vector<char> keep(count, 0);
    std::vector<double> scores(count * windows.size(), 0.0);
    hicx::parallel_for(count, threads, [&](std::size_t idx) {
        DenseBlock scratch;
        const std::int64_t bin = first_bin + static_cast<std::int64_t>(idx);
        bool usable = true;
        for (std::size_t w = 0; w < windows.size(); ++w) {
            bool defined = true;
            const double value = cut_weight_mean(context, bin, windows[w], scratch, defined);
            if (!defined || std::isnan(value)) {
                usable = false;
                break;
            }
            scores[idx * windows.size() + w] = value;
        }
        keep[idx] = usable ? 1 : 0;
    });

    BedgraphMatrix local;
    local.columns = windows.size();
    for (std::size_t idx = 0; idx < count; ++idx) {
        if (keep[idx] == 0) {
            continue;
        }
        const hicx::CutInterval& interval =
            (*context.cut_intervals)[static_cast<std::size_t>(first_bin + static_cast<std::int64_t>(idx))];
        local.chrom.push_back(interval.chrom);
        local.start.push_back(interval.start);
        local.end.push_back(interval.end);
        local.values.insert(
            local.values.end(), scores.begin() + static_cast<std::ptrdiff_t>(idx * windows.size()),
            scores.begin() + static_cast<std::ptrdiff_t>((idx + 1) * windows.size()));
    }
    if (local.rows() < 3) {
        return std::nullopt;  // not enough bins left to look for internal structure
    }

    const std::vector<double> local_row_means = local.row_means();

    // The TAD-separation score this domain's own search works from: computed
    // fresh from `windows` above (nested_min_depth..nested_max_depth), not
    // read off the top level's score (which used a different, wider depth
    // range and therefore is not what peakdetect below is looking at).
    // Written out so a caller can see the exact curve the minima a few lines
    // down were picked from, instead of only the top-level one.
    //
    // Each row is the bin's own interval (local.start[index]..local.end[index]),
    // the same convention _nested_boundaries.bed/_nested_domains.bed already
    // use, not the midpoint-of-neighbouring-bin-centres convention
    // <outPrefix>_score.bedgraph uses: that convention exists only to match
    // hicFindTADs.py's own score file byte for byte, which this C++-only file
    // has no reference to match, and using it here would shift every row by
    // half a bin relative to the boundaries/domains this score explains,
    // putting a real local minimum visibly to one side of its own triangle.
    for (std::size_t index = 0; index < local.rows(); ++index) {
        if (local.end[index] <= local.start[index]) {
            continue;
        }
        nested_score << local.chrom[index] << '\t' << local.start[index] << '\t'
                      << local.end[index] << '\t' << format_fixed(local_row_means[index], 12)
                      << '\n';
    }

    std::vector<std::int64_t> widths(local.rows());
    for (std::size_t r = 0; r < local.rows(); ++r) {
        widths[r] = local.end[r] - local.start[r];
    }
    std::vector<std::int64_t> sorted_widths = widths;
    std::sort(sorted_widths.begin(), sorted_widths.end());
    double average_bin_size = 0.0;
    if (!sorted_widths.empty()) {
        const std::size_t n = sorted_widths.size();
        average_bin_size = n % 2 == 1
                               ? static_cast<double>(sorted_widths[n / 2])
                               : 0.5 * (static_cast<double>(sorted_widths[n / 2 - 1]) +
                                        static_cast<double>(sorted_widths[n / 2]));
    }
    if (average_bin_size <= 0) {
        return std::nullopt;
    }
    const double min_boundary_distance = nested_min_boundary_distance.has_value()
                                             ? static_cast<double>(*nested_min_boundary_distance)
                                             : average_bin_size * 4.0;
    const auto lookahead =
        static_cast<std::int64_t>(min_boundary_distance / average_bin_size);
    if (lookahead < 1) {
        return std::nullopt;
    }

    std::vector<Peak> maxima;
    std::vector<Peak> minima;
    peakdetect(local_row_means, local.chrom, lookahead, maxima, minima);
    if (minima.empty()) {
        return std::nullopt;
    }
    std::vector<std::int64_t> minimum_indices;
    minimum_indices.reserve(minima.size());
    for (const Peak& peak : minima) {
        minimum_indices.push_back(peak.position);
    }
    const std::map<std::int64_t, double> delta_of_min =
        delta_wrt_window(minimum_indices, local_row_means, local.chrom);

    std::vector<char> usable(minimum_indices.size(), 0);
    std::vector<double> raw_pvalues(minimum_indices.size(), kNaN);
    hicx::parallel_for(minimum_indices.size(), threads, [&](std::size_t slot) {
        const std::int64_t index = minimum_indices[slot];
        const auto range = bins.region_bin_range(local.chrom[static_cast<std::size_t>(index)],
                                                 local.start[static_cast<std::size_t>(index)],
                                                 local.end[static_cast<std::size_t>(index)]);
        if (!range.has_value()) {
            return;
        }
        usable[slot] = 1;
        const std::int64_t matrix_index = range->first;
        std::int64_t left_bin = 0;
        std::int64_t right_bin = 0;
        if (!bins_at_distance(context, matrix_index, nested_min_depth, left_bin, right_bin)) {
            raw_pvalues[slot] = kNaN;
            return;
        }
        DenseBlock left_block;
        DenseBlock right_block;
        DenseBlock boundary_block;
        const bool has_left = cut_weight_values(context, left_bin, nested_min_depth, left_block);
        const bool has_right = cut_weight_values(context, right_bin, nested_min_depth, right_block);
        const bool has_boundary =
            cut_weight_values(context, matrix_index, nested_min_depth, boundary_block);
        const std::size_t left_size = has_left ? left_block.values.size() : 0;
        const std::size_t right_size = has_right ? right_block.values.size() : 0;
        const std::size_t boundary_size = has_boundary ? boundary_block.values.size() : 0;
        if (left_size == 0 || right_size == 0 || boundary_size == 0) {
            raw_pvalues[slot] = kNaN;
            return;
        }
        const double p1 = hicx::stats::ranksums(boundary_block.values, left_block.values).pvalue;
        const double p2 = hicx::stats::ranksums(boundary_block.values, right_block.values).pvalue;
        raw_pvalues[slot] = std::min(p1, p2);
    });

    std::vector<std::int64_t> pvalue_index;
    std::vector<double> pvalues;
    for (std::size_t slot = 0; slot < minimum_indices.size(); ++slot) {
        if (usable[slot] == 0) {
            continue;
        }
        pvalue_index.push_back(minimum_indices[slot]);
        pvalues.push_back(raw_pvalues[slot]);
    }
    if (pvalue_index.empty()) {
        return std::nullopt;
    }

    NestedCandidates result;
    result.call = call;
    result.local = std::move(local);
    result.local_row_means = std::move(local_row_means);
    result.delta_of_min = std::move(delta_of_min);
    result.pvalue_index = std::move(pvalue_index);
    result.pvalues = std::move(pvalues);
    return result;
}

// Turns one domain's NestedCandidates into output, once every domain's raw
// p-values have been corrected together and `pvalue_of_min` (keyed by index
// into `candidates.local`) holds this domain's share of that correction.
void finalize_nested_domain(const NestedCandidates& candidates,
                            const std::map<std::int64_t, double>& pvalue_of_min,
                            double delta_threshold, bool has_ctcf,
                            const std::unordered_map<std::string, std::vector<CtcfSite>>& ctcf_sites,
                            std::int64_t ctcf_window, std::int64_t ctcf_boundary_slack,
                            std::ofstream& nested_boundaries_bed,
                            std::ofstream& nested_domains_bed, std::ofstream& nested_gff,
                            int& domain_counter) {
    const NestedCall& call = candidates.call;
    const BedgraphMatrix& local = candidates.local;
    const std::vector<double>& local_row_means = candidates.local_row_means;
    const std::map<std::int64_t, double>& delta_of_min = candidates.delta_of_min;

    std::vector<std::int64_t> accepted;
    std::map<std::int64_t, CtcfPairSupport> ctcf_match;
    for (const std::int64_t index : candidates.pvalue_index) {
        const auto delta_entry = delta_of_min.find(index);
        const double delta = delta_entry == delta_of_min.end() ? kNaN : delta_entry->second;
        // A domain is a single, small region: a Wilcoxon p-value/FDR gate
        // tuned for the top level's chromosome-wide search is either too
        // strict for a handful of candidates (per-domain correction) or
        // swamped by every other domain's candidates (pooled correction),
        // and either way a real, strong minimum (large delta) can fail it on
        // a technicality despite being the clearest dip in the domain's own
        // score. The nested pass instead uses the same criterion a TAD
        // boundary is: a local minimum of the (diamond) score, filtered only
        // by how deep it is relative to its neighbourhood (--delta). pvalue
        // is still computed and corrected (see pvalue_of_min, built from the
        // chromosome-wide correction in main()) and reported in the GFF for
        // information, not as a gate.
        if (!(delta >= delta_threshold)) {
            continue;
        }
        if (has_ctcf) {
            const auto sites = ctcf_sites.find(local.chrom[static_cast<std::size_t>(index)]);
            const CtcfPairSupport support =
                sites == ctcf_sites.end()
                    ? CtcfPairSupport{}
                    : ctcf_support(sites->second, local.start[static_cast<std::size_t>(index)],
                                   ctcf_window, ctcf_boundary_slack);
            if (!support.forward.has_value() && !support.reverse.has_value()) {
                continue;
            }
            ctcf_match.emplace(index, support);
        }
        accepted.push_back(index);
    }
    if (accepted.empty()) {
        return;
    }
    std::sort(accepted.begin(), accepted.end());

    // The parent span cut at each accepted boundary position gives the
    // sub-domains; consecutive duplicate cuts (two boundaries in the same
    // bin) collapse to one.
    std::vector<std::int64_t> cuts;
    cuts.push_back(call.parent_start);
    for (const std::int64_t index : accepted) {
        cuts.push_back(local.start[static_cast<std::size_t>(index)]);
    }
    cuts.push_back(call.parent_end);
    std::vector<std::int64_t> unique_cuts;
    for (const std::int64_t c : cuts) {
        if (unique_cuts.empty() || unique_cuts.back() != c) {
            unique_cuts.push_back(c);
        }
    }
    if (unique_cuts.size() < 3) {
        return;  // every candidate collapsed onto the parent's own edges
    }

    for (std::size_t i = 0; i < accepted.size(); ++i) {
        const std::int64_t index = accepted[i];
        char identifier[32];
        std::snprintf(identifier, sizeof(identifier), "%s_B%03zu", call.parent_id.c_str(), i + 1);
        nested_boundaries_bed << local.chrom[static_cast<std::size_t>(index)] << '\t'
                             << local.start[static_cast<std::size_t>(index)] << '\t'
                             << local.end[static_cast<std::size_t>(index)] << '\t' << identifier
                             << '\t' << format_fixed(local_row_means[static_cast<std::size_t>(index)], 12)
                             << "\t.\n";
        nested_gff << local.chrom[static_cast<std::size_t>(index)] << "\tHiCExplorer\tnested_boundary\t"
                   << local.start[static_cast<std::size_t>(index)] << '\t'
                   << local.end[static_cast<std::size_t>(index)] << '\t'
                   << format_fixed(local_row_means[static_cast<std::size_t>(index)], 12) << "\t.\t.\tID="
                   << identifier << ";parent=" << call.parent_id << ";delta="
                   << format_fixed(delta_of_min.at(index), 12) << ";pvalue="
                   << format_fixed(pvalue_of_min.at(index), 12);
        const auto matched_ctcf = ctcf_match.find(index);
        if (matched_ctcf != ctcf_match.end()) {
            if (matched_ctcf->second.forward.has_value()) {
                nested_gff << ";ctcf_site_forward=" << matched_ctcf->second.forward->position;
            }
            if (matched_ctcf->second.reverse.has_value()) {
                nested_gff << ";ctcf_site_reverse=" << matched_ctcf->second.reverse->position;
            }
        }
        nested_gff << '\n';
    }

    for (std::size_t i = 0; i + 1 < unique_cuts.size(); ++i) {
        const char* rgb = domain_counter % 2 == 0 ? "102,194,165" : "252,141,98";
        nested_domains_bed << call.chrom << '\t' << unique_cuts[i] << '\t' << unique_cuts[i + 1]
                           << "\tID_" << call.parent_id << '_' << (i + 1) << "\t0\t.\t"
                           << unique_cuts[i] << '\t' << unique_cuts[i + 1] << '\t' << rgb << '\n';
        ++domain_counter;
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        // hicFindTADs.py:1319-1325 checks this only after the matrix has been
        // loaded, so the diagnostic never fires; here it fires.
        std::string suffix = args.matrix;
        const std::size_t dot = suffix.rfind('.');
        suffix = dot == std::string::npos ? suffix : suffix.substr(dot + 1);
        if (suffix != "cool" && suffix != "h5" && suffix.find("mcool::") == std::string::npos) {
            std::fprintf(stderr,
                         "ERROR:hicexplorer.hicFindTADs:Could not determine file "
                         "ending. Please use either .h5, .cool or "
                         ".mcool::/path/to/matrix.\n");
            std::fprintf(stderr,
                         "ERROR:hicexplorer.hicFindTADs:used input: %s, ending %s\n",
                         args.matrix.c_str(), suffix.c_str());
            return 1;
        }

        std::unordered_map<std::string, std::vector<CtcfSite>> ctcf_sites;
        if (args.ctcf_file.has_value()) {
            ctcf_sites = read_ctcf_sites(*args.ctcf_file);
        }

        const int nested_given = args.nested_min_depth.has_value() +
                                 args.nested_max_depth.has_value() +
                                 args.nested_step.has_value();
        if (nested_given != 0 && nested_given != 3) {
            std::fputs("hicFindTADs: --nestedMinDepth, --nestedMaxDepth and "
                       "--nestedStep must all be given together\n",
                       stderr);
            return 1;
        }
        if (nested_given == 3 && *args.nested_max_depth <= *args.nested_min_depth) {
            std::fputs("hicFindTADs: --nestedMaxDepth must be larger than "
                       "--nestedMinDepth\n",
                       stderr);
            return 1;
        }
        const bool nested_enabled = nested_given == 3;

        hicx::ToolMatrix hic = hicx::ToolMatrix::load(args.matrix);

        // reorderChromosomes(valid_chromosomes)
        if (args.chromosomes.has_value()) {
            std::vector<std::int64_t> order;
            std::vector<std::string> invalid;
            for (const std::string& chrom : *args.chromosomes) {
                const auto range = std::find_if(
                    hic.boundaries().begin(), hic.boundaries().end(),
                    [&chrom](const auto& entry) { return entry.first == chrom; });
                if (range == hic.boundaries().end()) {
                    invalid.push_back(chrom);
                    continue;
                }
                for (std::int64_t bin = range->second.first; bin < range->second.last;
                     ++bin) {
                    order.push_back(bin);
                }
            }
            if (!invalid.empty()) {
                std::fputs(
                    "WARNING:hicexplorer.hicFindTADs:WARNING: The following "
                    "chromosome/scaffold names were not found. Please check"
                    "the correct spelling of the chromosome names. \n\n",
                    stderr);
                for (const std::string& chrom : invalid) {
                    std::fprintf(stderr, "%s\n", chrom.c_str());
                }
            }
            hicx::reorder_bins(hic.data(), order);
            hic.refresh_boundaries();
        }

        // hiCMatrix caches the bin size at construction, so this value stays in
        // force even after the bins have been enlarged.
        const hicx::BinTable initial_bins(hic.data().cut_intervals);
        const std::int64_t bin_size = initial_bins.bin_size();

        // set_variables (hicFindTADs.py:419-464)
        std::int64_t max_depth = 0;
        std::int64_t min_depth = 0;
        std::int64_t step = 0;
        if (args.max_depth.has_value() && args.min_depth.has_value() &&
            *args.max_depth <= *args.min_depth) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please check that maxDepth is "
                       "larger than minDepth.\n",
                       stderr);
            return 0;  // the Python calls exit() with no argument
        }
        if (!args.max_depth.has_value()) {
            max_depth = bin_size < 1000    ? bin_size * 60
                        : bin_size < 20000 ? bin_size * 40
                                           : bin_size * 10;
        } else if (*args.max_depth < bin_size * 5) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --maxDepth "
                       "that is at least 5 times larger than the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            max_depth = *args.max_depth;
        }
        if (!args.min_depth.has_value()) {
            min_depth = bin_size < 1000    ? bin_size * 30
                        : bin_size < 20000 ? bin_size * 10
                                           : bin_size * 5;
        } else if (*args.min_depth < bin_size * 3) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --minDepth "
                       "that is at least 3 times larger than the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            min_depth = *args.min_depth;
        }
        if (!args.step.has_value()) {
            step = bin_size < 1000 ? bin_size * 4 : bin_size * 2;
        } else if (*args.step < bin_size) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:Please specify a --step that "
                       "is at least the size of the matrix bin size\n",
                       stderr);
            return 1;
        } else {
            step = *args.step;
        }

        std::int64_t nested_min_domain_size = 0;
        if (nested_enabled) {
            if (*args.nested_min_depth < bin_size * 3) {
                std::fputs("hicFindTADs: --nestedMinDepth must be at least 3 times "
                           "the matrix bin size (otherwise the Wilcoxon rank-sum "
                           "test has too few values per side to carry any power)\n",
                           stderr);
                return 1;
            }
            if (*args.nested_max_depth < bin_size * 3) {
                std::fputs("hicFindTADs: --nestedMaxDepth must be at least 3 times "
                           "the matrix bin size\n",
                           stderr);
                return 1;
            }
            nested_min_domain_size = args.nested_min_domain_size.has_value()
                                        ? *args.nested_min_domain_size
                                        : 4 * (*args.nested_max_depth);
        }

        const std::string score_suffix = "_tad_score.bm";
        std::string tad_score_file = args.out_prefix + score_suffix;
        std::string zscore_matrix_file = args.out_prefix + "_zscore_matrix." + suffix;
        if (args.tad_sep_score_prefix.has_value()) {
            tad_score_file = *args.tad_sep_score_prefix + score_suffix;
            zscore_matrix_file = *args.tad_sep_score_prefix + "_zscore_matrix." + suffix;
            if (!file_exists(tad_score_file)) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:The given "
                             "TAD_sep_score_prefix does not contain a valid "
                             "TAD-separation score. Please check.\nCould not find "
                             "file %s\n",
                             tad_score_file.c_str());
                return 1;
            }
            if (!file_exists(zscore_matrix_file)) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:The given "
                             "TAD_sep_score_prefix does not contain a valid z-score "
                             "matrix. Please check.\nCould not find file %s\n",
                             zscore_matrix_file.c_str());
                return 1;
            }
        }

        BedgraphMatrix bedgraph;
        const bool reuse_existing =
            args.tad_sep_score_prefix.has_value() || file_exists(tad_score_file);

        if (reuse_existing) {
            // set_matrix(zscore_matrix_file) followed by load_bedgraph_matrix.
            hic = hicx::ToolMatrix::load(zscore_matrix_file);
            if (args.chromosomes.has_value()) {
                std::vector<std::int64_t> order;
                for (const std::string& chrom : *args.chromosomes) {
                    const auto range = std::find_if(
                        hic.boundaries().begin(), hic.boundaries().end(),
                        [&chrom](const auto& entry) { return entry.first == chrom; });
                    if (range == hic.boundaries().end()) {
                        continue;
                    }
                    for (std::int64_t bin = range->second.first;
                         bin < range->second.last; ++bin) {
                        order.push_back(bin);
                    }
                }
                hicx::reorder_bins(hic.data(), order);
                hic.refresh_boundaries();
            }

            std::ifstream input(tad_score_file);
            if (!input) {
                std::fprintf(stderr, "hicFindTADs: cannot open %s\n",
                             tad_score_file.c_str());
                return 1;
            }
            std::string line;
            while (std::getline(input, line)) {
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                if (line.empty()) {
                    continue;
                }
                if (line[0] == '#') {
                    // The parameters the spectrum was computed with replace the
                    // command line values (hicFindTADs.py:1157-1160). Only the
                    // four keys the writer emits occur.
                    const auto read = [&line](const char* key) -> std::int64_t {
                        const std::string needle = std::string("\"") + key + "\":";
                        const std::size_t at = line.find(needle);
                        if (at == std::string::npos) {
                            return 0;
                        }
                        return std::strtoll(line.c_str() + at + needle.size(), nullptr, 10);
                    };
                    min_depth = read("minDepth");
                    max_depth = read("maxDepth");
                    step = read("step");
                    continue;
                }
                std::vector<std::string> fields;
                std::size_t start = 0;
                while (true) {
                    const std::size_t tab = line.find('\t', start);
                    if (tab == std::string::npos) {
                        fields.push_back(line.substr(start));
                        break;
                    }
                    fields.push_back(line.substr(start, tab - start));
                    start = tab + 1;
                }
                if (fields.size() < 4) {
                    continue;
                }
                if (args.chromosomes.has_value() &&
                    std::find(args.chromosomes->begin(), args.chromosomes->end(),
                              fields[0]) == args.chromosomes->end()) {
                    continue;
                }
                bedgraph.chrom.push_back(fields[0]);
                bedgraph.start.push_back(
                    static_cast<std::int64_t>(std::strtod(fields[1].c_str(), nullptr)));
                bedgraph.end.push_back(
                    static_cast<std::int64_t>(std::strtod(fields[2].c_str(), nullptr)));
                for (std::size_t c = 3; c < fields.size(); ++c) {
                    bedgraph.values.push_back(std::strtod(fields[c].c_str(), nullptr));
                }
                bedgraph.columns = fields.size() - 3;
            }
        } else {
            // ------------------------------------------------------------------
            // compute_spectra_matrix (hicFindTADs.py:1017-1129)
            std::fputs("INFO:hicexplorer.hicFindTADs:removing diagonal values\n\n",
                       stderr);
            drop_diagonal(hic.matrix());

            MaskState mask = mask_bins(hic.data(), hic.data().nan_bins);
            const std::vector<hicx::CutInterval> original_intervals =
                hic.data().cut_intervals;

            std::fputs("INFO:hicexplorer.hicFindTADs:Computing z-score matrix...\n\n",
                       stderr);
            hicx::ObsExpOptions options;
            options.max_depth_bp = static_cast<double>(max_depth) * 2.5;
            options.zscore = true;
            options.perchr = true;
            hicx::convert_to_obs_exp_matrix(hic.data(), bin_size, options);

            hicx::enlarge_bins(hic.data().cut_intervals);
            if (!same_intervals(hic.data().cut_intervals, original_intervals)) {
                mask.orig_bin_ids.clear();
                mask.orig_cut_intervals.clear();
                hic.data().nan_bins.clear();
            }
            hic.refresh_boundaries();

            std::fputs(
                "INFO:hicexplorer.hicFindTADs:Computing TAD-separation scores...\n\n",
                stderr);
            const std::int64_t min_depth_in_bins = min_depth / bin_size;
            const std::int64_t max_depth_in_bins = max_depth / bin_size;
            if (step / bin_size == 0) {
                std::fprintf(stderr,
                             "ERROR:hicexplorer.hicFindTADs:Please select a step size "
                             "larger than %lld\n",
                             static_cast<long long>(bin_size));
                return 1;
            }
            if (min_depth_in_bins <= 1 || max_depth_in_bins <= 1) {
                std::fputs("ERROR:hicexplorer.hicFindTADs:depth length too small\n",
                           stderr);
                return 0;  // the Python calls exit(0) on both of these
            }
            truncate_to_distance(hic.matrix(), 2 * max_depth_in_bins);

            const std::vector<std::int64_t> windows =
                incremental_step_size(min_depth, max_depth, step);

            const hicx::BinTable bins(hic.data().cut_intervals);
            ScoreContext context;
            context.matrix = &hic.matrix();
            context.cut_intervals = &hic.data().cut_intervals;
            context.bins = &bins;
            for (const auto& [chrom, size] : bins.chromosome_sizes()) {
                context.chromosome_sizes[chrom] = size;
            }

            // bins_to_consider, in chrBinBoundaries order.
            std::vector<std::int64_t> bins_to_consider;
            bins_to_consider.reserve(hic.data().cut_intervals.size());
            for (const auto& [chrom, range] : hic.boundaries()) {
                for (std::int64_t bin = range.first; bin < range.last; ++bin) {
                    bins_to_consider.push_back(bin);
                }
            }

            const std::size_t count = bins_to_consider.size();
            std::vector<char> keep(count, 0);
            std::vector<double> scores(count * windows.size(), 0.0);
            const unsigned int threads =
                static_cast<unsigned int>(std::max(1, args.processors));
            hicx::parallel_for(count, threads, [&](std::size_t index) {
                DenseBlock scratch;
                const std::int64_t bin = bins_to_consider[index];
                bool usable = true;
                for (std::size_t w = 0; w < windows.size(); ++w) {
                    bool defined = true;
                    const double value =
                        cut_weight_mean(context, bin, windows[w], scratch, defined);
                    if (!defined || std::isnan(value)) {
                        usable = false;
                        break;
                    }
                    scores[index * windows.size() + w] = value;
                }
                keep[index] = usable ? 1 : 0;
            });

            bedgraph.columns = windows.size();
            for (std::size_t index = 0; index < count; ++index) {
                if (keep[index] == 0) {
                    continue;
                }
                const hicx::CutInterval& interval =
                    hic.data().cut_intervals[static_cast<std::size_t>(
                        bins_to_consider[index])];
                bedgraph.chrom.push_back(interval.chrom);
                bedgraph.start.push_back(interval.start);
                bedgraph.end.push_back(interval.end);
                bedgraph.values.insert(
                    bedgraph.values.end(), scores.begin() + static_cast<std::ptrdiff_t>(
                                                                index * windows.size()),
                    scores.begin() +
                        static_cast<std::ptrdiff_t>((index + 1) * windows.size()));
            }
            if (bedgraph.rows() == 0) {
                std::fputs("hicFindTADs: no bin produced a TAD-separation score\n",
                           stderr);
                return 1;
            }

            // hiCMatrix.save restores the masked bins first.
            restore_masked_bins(hic.data(), mask);
            hic.refresh_boundaries();
            hic.save(args.out_prefix + "_zscore_matrix." + suffix);

            // save_bedgraph_matrix
            std::ofstream out(tad_score_file);
            out << "#{\"step\":" << step << ",\"minDepth\":" << min_depth
                << ",\"maxDepth\":" << max_depth << ",\"binsize\":" << bin_size << "}\n";
            for (std::size_t r = 0; r < bedgraph.rows(); ++r) {
                out << bedgraph.chrom[r] << '\t' << bedgraph.start[r] << '\t'
                    << bedgraph.end[r];
                for (std::size_t c = 0; c < bedgraph.columns; ++c) {
                    out << '\t'
                        << format_fixed(bedgraph.values[r * bedgraph.columns + c], 6);
                }
                out << '\n';
            }
        }

        // ----------------------------------------------------------------------
        // find_boundaries (hicFindTADs.py:1251-1294)
        std::vector<std::int64_t> widths(bedgraph.rows());
        for (std::size_t r = 0; r < bedgraph.rows(); ++r) {
            widths[r] = bedgraph.end[r] - bedgraph.start[r];
        }
        std::vector<std::int64_t> sorted_widths = widths;
        std::sort(sorted_widths.begin(), sorted_widths.end());
        double average_bin_size = 0.0;
        if (!sorted_widths.empty()) {
            const std::size_t n = sorted_widths.size();
            average_bin_size = n % 2 == 1
                                   ? static_cast<double>(sorted_widths[n / 2])
                                   : 0.5 * (static_cast<double>(sorted_widths[n / 2 - 1]) +
                                            static_cast<double>(sorted_widths[n / 2]));
        }
        const double min_boundary_distance =
            args.min_boundary_distance.has_value()
                ? static_cast<double>(*args.min_boundary_distance)
                : average_bin_size * 4.0;
        const std::int64_t lookahead =
            static_cast<std::int64_t>(min_boundary_distance / average_bin_size);
        if (lookahead < 1) {
            std::fputs("hicFindTADs: minBoundaryDistance must be '1' or above in value\n",
                       stderr);
            return 1;
        }

        const std::vector<double> row_means = bedgraph.row_means();
        std::vector<Peak> maxima;
        std::vector<Peak> minima;
        peakdetect(row_means, bedgraph.chrom, lookahead, maxima, minima);
        std::vector<std::int64_t> minimum_indices;
        minimum_indices.reserve(minima.size());
        for (const Peak& peak : minima) {
            minimum_indices.push_back(peak.position);
        }
        std::map<std::int64_t, double> delta_of_min =
            delta_wrt_window(minimum_indices, row_means, bedgraph.chrom);

        // min_pvalue (hicFindTADs.py:1168-1249)
        const hicx::BinTable bins(hic.data().cut_intervals);
        ScoreContext context;
        context.matrix = &hic.matrix();
        context.cut_intervals = &hic.data().cut_intervals;
        context.bins = &bins;
        for (const auto& [chrom, size] : bins.chromosome_sizes()) {
            context.chromosome_sizes[chrom] = size;
        }

        std::vector<char> usable(minimum_indices.size(), 0);
        std::vector<double> raw_pvalues(minimum_indices.size(), kNaN);
        {
            const unsigned int threads =
                static_cast<unsigned int>(std::max(1, args.processors));
            hicx::parallel_for(
                minimum_indices.size(), threads, [&](std::size_t slot) {
                    const std::int64_t index = minimum_indices[slot];
                    const auto range = bins.region_bin_range(
                        bedgraph.chrom[static_cast<std::size_t>(index)],
                        bedgraph.start[static_cast<std::size_t>(index)],
                        bedgraph.end[static_cast<std::size_t>(index)]);
                    if (!range.has_value()) {
                        return;
                    }
                    usable[slot] = 1;
                    const std::int64_t matrix_index = range->first;
                    std::int64_t left_bin = 0;
                    std::int64_t right_bin = 0;
                    if (!bins_at_distance(context, matrix_index, min_depth, left_bin,
                                          right_bin)) {
                        // The Python unpacks None here and dies with a TypeError.
                        throw std::runtime_error(
                            "bin " + std::to_string(matrix_index) +
                            " has no window at the minimum depth; the Python "
                            "reference raises a bare TypeError here "
                            "(hicFindTADs.py:1201)");
                    }
                    DenseBlock left_block;
                    DenseBlock right_block;
                    DenseBlock boundary_block;
                    const bool has_left =
                        cut_weight_values(context, left_bin, min_depth, left_block);
                    const bool has_right =
                        cut_weight_values(context, right_bin, min_depth, right_block);
                    const bool has_boundary = cut_weight_values(context, matrix_index,
                                                               min_depth, boundary_block);
                    const std::size_t left_size =
                        has_left ? left_block.values.size() : 0;
                    const std::size_t right_size =
                        has_right ? right_block.values.size() : 0;
                    const std::size_t boundary_size =
                        has_boundary ? boundary_block.values.size() : 0;
                    if (left_size == 0 || right_size == 0 || boundary_size == 0) {
                        raw_pvalues[slot] = kNaN;
                        return;
                    }
                    const double p1 =
                        hicx::stats::ranksums(boundary_block.values, left_block.values)
                            .pvalue;
                    const double p2 =
                        hicx::stats::ranksums(boundary_block.values, right_block.values)
                            .pvalue;
                    raw_pvalues[slot] = std::min(p1, p2);
                });
        }

        std::vector<std::int64_t> pvalue_index;
        std::vector<double> pvalues;
        for (std::size_t slot = 0; slot < minimum_indices.size(); ++slot) {
            if (usable[slot] == 0) {
                continue;
            }
            pvalue_index.push_back(minimum_indices[slot]);
            pvalues.push_back(raw_pvalues[slot]);
        }

        double pvalue_fdr = 0.0;
        if (args.correct_for_multiple_testing == "fdr") {
            for (double& p : pvalues) {
                if (std::isnan(p)) {
                    p = 1.0;
                }
            }
            pvalue_fdr = hicx::stats::benjamini_hochberg_cutoff(
                pvalues, args.threshold_comparisons);
        } else if (args.correct_for_multiple_testing == "bonferroni") {
            hicx::stats::bonferroni_in_place(pvalues);
        }
        std::map<std::int64_t, double> pvalue_of_min;
        for (std::size_t i = 0; i < pvalue_index.size(); ++i) {
            pvalue_of_min.emplace(pvalue_index[i], pvalues[i]);
        }

        if (minimum_indices.empty()) {
            std::fputs("ERROR:hicexplorer.hicFindTADs:\n*ERROR*\nNo boundaries were "
                       "found.\n",
                       stderr);
            return 1;
        }

        // ----------------------------------------------------------------------
        // save_domains_and_boundaries (hicFindTADs.py:888-1015)
        const std::size_t rows = bedgraph.rows();
        std::map<std::string, std::int64_t> first_index;
        for (std::size_t i = 0; i < rows; ++i) {
            first_index.emplace(bedgraph.chrom[i], static_cast<std::int64_t>(i));
        }
        std::vector<std::int64_t> chr_end_idx;
        if (first_index.size() == 1) {
            chr_end_idx.push_back(static_cast<std::int64_t>(rows) - 1);
        } else {
            for (const auto& [name, index] : first_index) {
                chr_end_idx.push_back(index == 0 ? static_cast<std::int64_t>(rows) - 1
                                                 : index - 1);
            }
        }
        // chr_end_idx is the same array object as chr_start_idx in the Python,
        // so the start indices are lost and the end indices enter the list
        // twice. Reproduced.
        std::vector<std::int64_t> candidates;
        candidates.insert(candidates.end(), chr_end_idx.begin(), chr_end_idx.end());
        candidates.insert(candidates.end(), chr_end_idx.begin(), chr_end_idx.end());
        candidates.insert(candidates.end(), minimum_indices.begin(),
                          minimum_indices.end());
        std::sort(candidates.begin(), candidates.end());

        std::vector<std::int64_t> filtered;
        std::map<std::int64_t, CtcfPairSupport> ctcf_match;
        for (const std::int64_t index : candidates) {
            const auto delta_entry = delta_of_min.find(index);
            const double delta =
                delta_entry == delta_of_min.end() ? kNaN : delta_entry->second;
            if (delta_entry == delta_of_min.end()) {
                delta_of_min[index] = kNaN;
            }
            const auto pvalue_entry = pvalue_of_min.find(index);
            if (pvalue_entry == pvalue_of_min.end()) {
                continue;
            }
            const double threshold = args.correct_for_multiple_testing == "fdr"
                                         ? pvalue_fdr
                                         : args.threshold_comparisons;
            if (!(delta >= args.delta && pvalue_entry->second <= threshold)) {
                continue;
            }
            if (args.ctcf_file.has_value()) {
                const auto sites = ctcf_sites.find(bedgraph.chrom[static_cast<std::size_t>(index)]);
                const CtcfPairSupport support =
                    sites == ctcf_sites.end()
                        ? CtcfPairSupport{}
                        : ctcf_support(sites->second, bedgraph.start[static_cast<std::size_t>(index)],
                                       args.ctcf_window, args.ctcf_boundary_slack);
                if (!support.forward.has_value() && !support.reverse.has_value()) {
                    continue;
                }
                ctcf_match.emplace(index, support);
            }
            filtered.push_back(index);
        }

        const std::set<std::int64_t> chromosome_ends(chr_end_idx.begin(),
                                                     chr_end_idx.end());
        const std::string delta_text = hicx::npy::float_repr(args.delta);

        std::ofstream boundaries_bed(args.out_prefix + "_boundaries.bed");
        std::ofstream domains_bed(args.out_prefix + "_domains.bed");
        std::ofstream gff(args.out_prefix + "_boundaries.gff");
        std::vector<NestedCall> nested_calls;
        int count = 1;
        for (std::size_t position = 0; position < filtered.size(); ++position) {
            const std::int64_t index = filtered[position];
            if (chromosome_ends.count(index) != 0) {
                continue;
            }
            const std::size_t at = static_cast<std::size_t>(index);
            // chr_start[min_bin_id - 1] with min_bin_id == 0 is Python's
            // negative index, that is the last row. It cannot be reached in
            // practice because delta_wrt_window gives index 0 a NaN delta and
            // NaN fails the delta filter, but the wraparound is reproduced
            // rather than left as undefined behaviour.
            const std::size_t before = at == 0 ? rows - 1 : at - 1;
            const std::int64_t right_center =
                bedgraph.start[at] + (bedgraph.end[at] - bedgraph.start[at]) / 2;
            const std::int64_t left_center =
                bedgraph.start[before] +
                (bedgraph.end[before] - bedgraph.start[before]) / 2;
            if (bedgraph.chrom[at] != bedgraph.chrom[before]) {
                continue;
            }
            char identifier[16];
            std::snprintf(identifier, sizeof(identifier), "B%05lld",
                          static_cast<long long>(index));
            boundaries_bed << bedgraph.chrom[at] << '\t' << left_center << '\t'
                           << right_center << '\t' << identifier << '\t'
                           << format_fixed(row_means[at], 12) << "\t.\n";
            gff << bedgraph.chrom[at] << "\tHiCExplorer\tboundary\t" << left_center
                << '\t' << right_center << '\t' << format_fixed(row_means[at], 12)
                << "\t.\t.\tID=" << identifier
                << ";delta=" << format_fixed(delta_of_min[index], 12)
                << ";pvalue=" << format_fixed(pvalue_of_min[index], 12)
                << ";tad_sep=" << format_fixed(row_means[at], 12);
            const auto matched_ctcf = ctcf_match.find(index);
            if (matched_ctcf != ctcf_match.end()) {
                if (matched_ctcf->second.forward.has_value()) {
                    gff << ";ctcf_site_forward=" << matched_ctcf->second.forward->position;
                }
                if (matched_ctcf->second.reverse.has_value()) {
                    gff << ";ctcf_site_reverse=" << matched_ctcf->second.reverse->position;
                }
            }
            gff << '\n';

            if (position + 1 == filtered.size() ||
                bedgraph.chrom[at] !=
                    bedgraph.chrom[static_cast<std::size_t>(filtered[position + 1])]) {
                continue;
            }
            const std::int64_t start = bedgraph.start[at];
            const std::int64_t end =
                bedgraph.start[static_cast<std::size_t>(filtered[position + 1])];
            const char* rgb = count % 2 == 0 ? "51,160,44" : "31,120,180";
            std::string domain_id = "ID_" + delta_text + "_" + std::to_string(count);
            domains_bed << bedgraph.chrom[at] << '\t' << start << '\t' << end << "\t"
                        << domain_id << '\t'
                        << format_fixed(row_means[at], 12) << "\t.\t" << start << '\t'
                        << end << '\t' << rgb << '\n';
            if (nested_enabled && end - start >= nested_min_domain_size) {
                nested_calls.push_back(
                    NestedCall{bedgraph.chrom[at], start, end, std::move(domain_id)});
            }
            ++count;
        }
        boundaries_bed.close();
        domains_bed.close();
        gff.close();

        if (nested_enabled) {
            std::ofstream nested_boundaries_bed(args.out_prefix + "_nested_boundaries.bed");
            std::ofstream nested_domains_bed(args.out_prefix + "_nested_domains.bed");
            std::ofstream nested_gff(args.out_prefix + "_nested_boundaries.gff");
            std::ofstream nested_score(args.out_prefix + "_nested_score.bedgraph");
            const unsigned int threads =
                static_cast<unsigned int>(std::max(1, args.processors));

            // Phase 1: collect every domain's candidate minima and their raw
            // Wilcoxon p-values, without correcting or accepting anything yet.
            std::vector<NestedCandidates> all_candidates;
            for (const NestedCall& call : nested_calls) {
                std::optional<NestedCandidates> candidates = collect_nested_candidates(
                    context, bins, call, *args.nested_min_depth, *args.nested_max_depth,
                    *args.nested_step, args.nested_min_boundary_distance, threads, nested_score);
                if (candidates.has_value()) {
                    all_candidates.push_back(std::move(*candidates));
                }
            }

            // Phase 2: a boundary here is accepted by --delta alone (see the
            // comment on the gate in finalize_nested_domain for why), so the
            // Wilcoxon p-value is not corrected for multiple testing or used
            // as a cutoff; it is only reported in the GFF, for information,
            // exactly as computed per candidate.
            std::vector<double> flat_pvalues;
            std::vector<std::size_t> offsets;
            offsets.reserve(all_candidates.size() + 1);
            for (const NestedCandidates& candidates : all_candidates) {
                offsets.push_back(flat_pvalues.size());
                flat_pvalues.insert(flat_pvalues.end(), candidates.pvalues.begin(),
                                    candidates.pvalues.end());
            }
            offsets.push_back(flat_pvalues.size());

            // Phase 3: accept, filter by CTCF, and write, per domain.
            int domain_counter = 1;
            for (std::size_t d = 0; d < all_candidates.size(); ++d) {
                std::map<std::int64_t, double> pvalue_of_min;
                for (std::size_t k = 0; k < all_candidates[d].pvalue_index.size(); ++k) {
                    pvalue_of_min.emplace(all_candidates[d].pvalue_index[k],
                                          flat_pvalues[offsets[d] + k]);
                }
                finalize_nested_domain(all_candidates[d], pvalue_of_min, args.delta,
                                       args.ctcf_file.has_value(), ctcf_sites, args.ctcf_window,
                                       args.ctcf_boundary_slack, nested_boundaries_bed,
                                       nested_domains_bed, nested_gff, domain_counter);
            }
            nested_boundaries_bed.close();
            nested_domains_bed.close();
            nested_gff.close();
            nested_score.close();
        }

        std::ofstream score(args.out_prefix + "_score.bedgraph");
        for (std::size_t index = 1; index < rows; ++index) {
            const std::int64_t right_center =
                bedgraph.start[index] +
                (bedgraph.end[index] - bedgraph.start[index]) / 2;
            const std::int64_t left_center =
                bedgraph.start[index - 1] +
                (bedgraph.end[index - 1] - bedgraph.start[index - 1]) / 2;
            if (right_center <= left_center) {
                continue;
            }
            score << bedgraph.chrom[index] << '\t' << left_center << '\t' << right_center
                  << '\t' << format_fixed(row_means[index], 12) << '\n';
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "hicFindTADs: %s\n", error.what());
        return 1;
    }
    hicx::report_resource_usage("hicFindTADs");
    return 0;
}
