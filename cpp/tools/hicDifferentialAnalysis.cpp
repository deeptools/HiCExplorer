// hicDifferentialAnalysis: replicate-aware, count-based differential analysis
// of TADs, loops and compartments (cpp/PLAN.md 9.7, work item 2).
//
// Not a port. HiCExplorer's differential tools test a single matrix against a
// single matrix, bin by bin or block by block, with filters applied per sample
// and no correction for multiple testing; PLAN.md 9.7 measured the false
// positives that produces. This tool replaces the comparison by a count model:
//
//   * samples are grouped into two conditions with replicates; a single sample
//     per condition is refused unless --exploratory is given, and then every
//     output file says so;
//   * a bin filtered in any sample (hicCorrectMatrix's MAD rule on the cis
//     coverage) is masked in all of them;
//   * every tested unit is a raw count per sample with an offset for library
//     size and distance decay (the per-sample, per-chromosome mean count at
//     each distance) or for a local background count, tested by the negative
//     binomial quasi-likelihood model of diff_engine_impl.hpp against a minimum
//     fold change (TREAT), with Benjamini-Hochberg FDR over the units.
//
// Subcommand `tads`: for every TAD of a domain file, the intra-TAD contacts
// aggregated per distance stratum (distances of [1, 2), [2, 4), [4, 8), ...
// bins) and in total. Each stratum and the total are tested as their own
// families; the TAD's p-value is Simes' combination of its total and stratum
// p-values, and the FDR is taken over TADs. Separately, for every boundary
// between two adjacent TADs, the contacts crossing it within a window, with
// the contacts inside the two flanking windows as the offset: a change of
// boundary insulation.
//
// Calibration devices (cpp/scripts/diff_calibration.py): --splitReplicates
// divides every count of the same file into two binomial halves, the first to
// condition A and the second to condition B; --plantRegions thins the contacts
// of chosen regions in one condition by 1 / --plantFold. Both are exact count
// operations with seeded, pixel-keyed random draws (diff_contacts_impl.hpp).

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "diff_contacts_impl.hpp"
#include "diff_engine_impl.hpp"
#include "hicx/argparse.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/parallel.hpp"
#include "hicx/resource_usage.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/version.hpp"

namespace {

namespace diff = hicx::diff;
namespace diffc = hicx::diffc;

constexpr const char* kProg = "hicDifferentialAnalysis";
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t kStratumFamilyMinimum = 100;

const char* const kUsage =
    "usage: hicDifferentialAnalysis [-h] [--version] {tads,loops,compartments} ...\n";

const char* const kHelp =
    "\n"
    "Replicate-aware, count-based differential analysis of Hi-C contact matrices: TADs\n"
    "and TAD boundaries (tads), loops (loops) and A/B compartments (compartments).\n"
    "Counts are modelled with a negative binomial GLM, the\n"
    "dispersion is estimated from the replicates with quasi-likelihood empirical Bayes\n"
    "moderation, offsets take library size and distance decay, bins filtered in any\n"
    "sample are masked in all, a minimum fold change is tested (TREAT) and the\n"
    "p-values are adjusted with Benjamini-Hochberg. Not in the Python HiCExplorer.\n"
    "\n"
    "positional arguments:\n"
    "  {tads,loops,compartments}\n"
    "    tads         differential TADs and TAD boundaries\n"
    "    loops        differential loops, each against its local background\n"
    "    compartments differential compartment preference of every bin\n"
    "\n"
    "options:\n"
    "  -h, --help     show this help message and exit\n"
    "  --version      show program's version number and exit\n";

const char* const kTadsUsage =
    "usage: hicDifferentialAnalysis tads --conditionA MATRIX [MATRIX ...]\n"
    "                                    --conditionB MATRIX [MATRIX ...]\n"
    "                                    --domains DOMAINS --outFilePrefix PREFIX\n"
    "                                    [--blocks LABEL [LABEL ...]] [--exploratory]\n"
    "                                    [--fdr FDR] [--minFoldChange FOLD]\n"
    "                                    [--filterThreshold LOWER UPPER]\n"
    "                                    [--chromosomes CHROM [CHROM ...]]\n"
    "                                    [--boundaryWindow BP] [--threads THREADS]\n"
    "                                    [--splitReplicates SEED]\n"
    "                                    [--plantRegions BED] [--plantFold FOLD]\n"
    "                                    [--plantSeed SEED] [-h]\n";

const char* const kCommonRequiredHelp =
    "Required arguments:\n"
    "  --conditionA MATRIX [MATRIX ...], -a MATRIX [MATRIX ...]\n"
    "                        Cool files of condition A, one per replicate.\n"
    "  --conditionB MATRIX [MATRIX ...], -b MATRIX [MATRIX ...]\n"
    "                        Cool files of condition B, one per replicate. All files\n"
    "                        must share one bin table. Raw counts are used; a weight\n"
    "                        column is ignored.\n";

const char* const kCommonOptionalHelp =
    "  --blocks LABEL [LABEL ...]\n"
    "                        One block label per sample, condition A first, for an\n"
    "                        additive block factor (for example a batch or a paired\n"
    "                        replicate). The block must not be confounded with the\n"
    "                        condition.\n"
    "  --exploratory         Allow a condition with a single sample. The dispersion is\n"
    "                        then estimated from all samples ignoring the conditions,\n"
    "                        and every output file is labelled exploratory. Without\n"
    "                        this option such a comparison is refused.\n"
    "  --fdr FDR             Benjamini-Hochberg FDR at which a unit is called\n"
    "                        differential (Default: 0.05).\n"
    "  --minFoldChange FOLD  Minimum fold change the test is against (TREAT); 1 tests\n"
    "                        for any change (Default: 1.1).\n"
    "  --filterThreshold LOWER UPPER\n"
    "                        MAD z-score limits on the cis coverage of a bin, per\n"
    "                        sample and chromosome; a bin outside them in any sample\n"
    "                        is masked in all (Default: -1.5 5.0).\n"
    "  --chromosomes CHROM [CHROM ...]\n"
    "                        Chromosomes to analyse (Default: all).\n"
    "  --threads THREADS, -t THREADS\n"
    "                        Worker threads; the result does not depend on the number\n"
    "                        (Default: 4).\n";

const char* const kCalibrationHelp =
    "\n"
    "Calibration arguments:\n"
    "  --splitReplicates SEED\n"
    "                        Split every count of each file into two binomial halves:\n"
    "                        the first halves form condition A and the second halves\n"
    "                        condition B. --conditionA and --conditionB must list the\n"
    "                        same files in the same order.\n"
    "  --plantRegions BED    Regions whose contacts are thinned to 1 / --plantFold in\n"
    "                        one condition: chrom, start, end and an optional fourth\n"
    "                        column A or B (Default: A).\n"
    "  --plantFold FOLD      Fold difference of the planted regions (Default: 2.0).\n"
    "  --plantSeed SEED      Seed of the planted thinning (Default: 0).\n"
    "  -h, --help            show this help message and exit\n";

const char* const kTadsHelpHead =
    "\n"
    "Differential TADs and TAD boundaries. Each TAD is tested on its intra-TAD contacts\n"
    "per distance stratum and in total (Simes' combination); each boundary between two\n"
    "adjacent TADs on the contacts crossing it relative to its two flanks. Writes\n"
    "<prefix>_tads.tsv and <prefix>_boundaries.tsv.\n"
    "\n";

const char* const kTadsSpecificHelp =
    "  --domains DOMAINS, -d DOMAINS\n"
    "                        TAD domains, a BED file such as hicFindTADs' domains.bed.\n"
    "  --outFilePrefix PREFIX, -o PREFIX\n"
    "                        Prefix of the output files.\n"
    "\n"
    "Optional arguments:\n";

const char* const kTadsOptionalHelp =
    "  --boundaryWindow BP   Width of the windows on either side of a boundary, capped\n"
    "                        by the two TADs (Default: 500000).\n";

const char* const kLoopsUsage =
    "usage: hicDifferentialAnalysis loops --conditionA MATRIX [MATRIX ...]\n"
    "                                     --conditionB MATRIX [MATRIX ...]\n"
    "                                     --loops LOOPS [LOOPS ...]\n"
    "                                     --outFilePrefix PREFIX\n"
    "                                     [--blocks LABEL [LABEL ...]] [--exploratory]\n"
    "                                     [--fdr FDR] [--minFoldChange FOLD]\n"
    "                                     [--filterThreshold LOWER UPPER]\n"
    "                                     [--chromosomes CHROM [CHROM ...]]\n"
    "                                     [--peakWidth BINS] [--windowSize BINS]\n"
    "                                     [--threads THREADS]\n"
    "                                     [--splitReplicates SEED]\n"
    "                                     [--plantRegions BED] [--plantFold FOLD]\n"
    "                                     [--plantSeed SEED] [-h]\n";

const char* const kLoopsHelpHead =
    "\n"
    "Differential loops. The loop calls of all samples are united: positions mapped to\n"
    "bins of the matrices, and calls within one bin of each other (in both anchors)\n"
    "merged into one loop at their rounded mean position. Each loop is tested on the\n"
    "contacts of its peak square (the loop pixel +- --peakWidth bins) with the contacts\n"
    "of its local background as the offset: the square of +- --windowSize bins without\n"
    "the peak square and the ring of one bin around it. A change is a change of the\n"
    "loop's enrichment over its background. Writes <prefix>_loops.tsv.\n"
    "\n";

const char* const kLoopsSpecificHelp =
    "  --loops LOOPS [LOOPS ...]\n"
    "                        Loop calls, one file or more (for example one per\n"
    "                        sample): BEDPE-like, the first six columns chrom1 start1\n"
    "                        end1 chrom2 start2 end2, as hicDetectLoops writes them.\n"
    "                        Inter-chromosomal calls are ignored.\n"
    "  --outFilePrefix PREFIX, -o PREFIX\n"
    "                        Prefix of the output file.\n"
    "\n"
    "Optional arguments:\n";

const char* const kCompartmentsUsage =
    "usage: hicDifferentialAnalysis compartments --conditionA MATRIX [MATRIX ...]\n"
    "                                            --conditionB MATRIX [MATRIX ...]\n"
    "                                            --gcTrack BEDGRAPH\n"
    "                                            --outFilePrefix PREFIX\n"
    "                                            [--blocks LABEL [LABEL ...]]\n"
    "                                            [--exploratory] [--fdr FDR]\n"
    "                                            [--minFoldChange FOLD]\n"
    "                                            [--filterThreshold LOWER UPPER]\n"
    "                                            [--chromosomes CHROM [CHROM ...]]\n"
    "                                            [--minDistance BP] [--threads THREADS]\n"
    "                                            [--splitReplicates SEED]\n"
    "                                            [--plantRegions BED] [--plantFold FOLD]\n"
    "                                            [--plantSeed SEED] [-h]\n";

const char* const kCompartmentsHelpHead =
    "\n"
    "Differential compartment preference. A consensus compartment track is computed\n"
    "per chromosome from the pooled samples: the Pearson correlation of the pooled\n"
    "observed over expected matrix, of its three leading eigenvectors the one that\n"
    "correlates most with GC content, oriented so that A (positive) is GC rich. Every\n"
    "bin is then tested on its contacts with A bins, with its contacts with B bins as\n"
    "the offset (both at a distance of at least --minDistance): a change is a change of\n"
    "the bin's preference for the A over the B compartment. The per-sample score is\n"
    "log2 of observed over expected contacts with A over the same with B. Writes\n"
    "<prefix>_compartments.tsv. A planted BED region thins the region's contacts with\n"
    "the A bins, except the contacts between region bins of opposite compartments.\n"
    "Dense per chromosome: use a resolution of 50 kb or coarser.\n"
    "\n";

const char* const kCompartmentsSpecificHelp =
    "  --gcTrack BEDGRAPH    GC content: chrom, start, end, value; averaged over each\n"
    "                        bin by overlap.\n"
    "  --outFilePrefix PREFIX, -o PREFIX\n"
    "                        Prefix of the output file.\n"
    "\n"
    "Optional arguments:\n";

const char* const kCompartmentsOptionalHelp =
    "  --minDistance BP      Least distance of the contacts that are counted\n"
    "                        (Default: 200000).\n";

const char* const kLoopsOptionalHelp =
    "  --peakWidth BINS      Half width of the peak square in bins (Default: 1).\n"
    "  --windowSize BINS     Half width of the background square in bins; at least\n"
    "                        --peakWidth + 2 (Default: 5).\n";

struct Arguments {
    std::string command;
    std::vector<std::string> condition_a;
    std::vector<std::string> condition_b;
    std::vector<std::string> blocks;
    std::vector<std::string> chromosomes;
    bool exploratory = false;
    double fdr = 0.05;
    double min_fold_change = 1.1;
    double filter_lower = -1.5;
    double filter_upper = 5.0;
    std::string prefix;
    std::int64_t threads = 4;
    std::optional<std::int64_t> split_seed;
    std::optional<std::string> plant_regions;
    double plant_fold = 2.0;
    std::int64_t plant_seed = 0;
    std::string domains;
    std::int64_t boundary_window = 500000;
    std::vector<std::string> loops;
    std::int64_t peak_width = 1;
    std::int64_t window_size = 5;
    std::string gc_track;
    std::int64_t min_distance = 200000;
};

void add_common(hicx::cli::Parser& sub, hicx::cli::ArgumentGroup& required,
                hicx::cli::ArgumentGroup& optional, hicx::cli::ArgumentGroup& calibration) {
    namespace cli = hicx::cli;
    (void)sub;
    required.add({"--conditionA", "-a"})
        .nargs("+")
        .required()
        .metavar("MATRIX")
        .input({"cool", "mcool"})
        .help("Cool files of condition A, one per replicate.");
    required.add({"--conditionB", "-b"})
        .nargs("+")
        .required()
        .metavar("MATRIX")
        .input({"cool", "mcool"})
        .help("Cool files of condition B, one per replicate, on the bin table of condition A.");
    optional.add({"--blocks"})
        .nargs("+")
        .metavar("LABEL")
        .help("One block label per sample, condition A first, for an additive block factor.");
    optional.add({"--exploratory"})
        .action(cli::Action::StoreTrue)
        .help("Allow a condition with a single sample; the output is labelled exploratory.");
    optional.add({"--fdr"})
        .type("float")
        .default_value(0.05)
        .help("Benjamini-Hochberg FDR at which a unit is called differential.");
    optional.add({"--minFoldChange"})
        .type("float")
        .default_value(1.1)
        .help("Minimum fold change the test is against (TREAT); 1 tests for any change.");
    optional.add({"--filterThreshold"})
        .type("float")
        .nargs(2)
        .metavar("LIMIT")
        .help("MAD z-score limits LOWER UPPER on the cis coverage of a bin (Default: -1.5 5.0).");
    optional.add({"--chromosomes"})
        .nargs("+")
        .metavar("CHROM")
        .help("Chromosomes to analyse (Default: all).");
    calibration.add({"--splitReplicates"})
        .type("int")
        .metavar("SEED")
        .help("Split every count into binomial halves: first halves to A, second halves to B.");
    calibration.add({"--plantRegions"})
        .metavar("BED")
        .input({"bed"})
        .help("Regions whose contacts are thinned to 1 / --plantFold in one condition.");
    calibration.add({"--plantFold"})
        .type("float")
        .default_value(2.0)
        .help("Fold difference of the planted regions.");
    calibration.add({"--plantSeed"})
        .type("int")
        .default_value(0)
        .help("Seed of the planted thinning.");
}

Arguments parse_arguments(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser(kProg, "Replicate-aware, count-based differential analysis of Hi-C "
                              "contact matrices (not in the Python HiCExplorer).");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);
    cli::ArgumentGroup& options = parser.group("options");
    options.add({"-h", "--help"}).action(cli::Action::Help).help("show this help message and exit");
    options.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    parser.subcommands("command", true, std::string("{tads,loops,compartments}"));

    cli::Parser& tads = parser.add_subcommand("tads", "Differential TADs and TAD boundaries.");
    tads.set_usage(kTadsUsage)
        .set_help(std::string(kTadsHelpHead) + kCommonRequiredHelp + kTadsSpecificHelp +
                  kCommonOptionalHelp + kTadsOptionalHelp + kCalibrationHelp);
    {
        cli::ArgumentGroup& required = tads.group("Required arguments");
        cli::ArgumentGroup& optional = tads.group("Optional arguments");
        cli::ArgumentGroup& calibration = tads.group("Calibration arguments");
        add_common(tads, required, optional, calibration);
        required.add({"--domains", "-d"})
            .required()
            .input({"bed"})
            .help("TAD domains, a BED file such as hicFindTADs' domains.bed.");
        required.add({"--outFilePrefix", "-o"})
            .required()
            .metavar("PREFIX")
            .output({"tsv"}, "prefix")
            .help("Prefix of the output files <prefix>_tads.tsv and <prefix>_boundaries.tsv.");
        optional.add({"--boundaryWindow"})
            .type("int")
            .default_value(500000)
            .metavar("BP")
            .help("Width of the windows on either side of a boundary, capped by the two TADs.");
        optional.add({"--threads", "-t"})
            .type("int")
            .default_value(4)
            .help("Worker threads; the result does not depend on the number.");
        calibration.add({"-h", "--help"})
            .action(cli::Action::Help)
            .help("show this help message and exit");
    }

    cli::Parser& loops = parser.add_subcommand(
        "loops", "Differential loops, each against its local background.");
    loops.set_usage(kLoopsUsage)
        .set_help(std::string(kLoopsHelpHead) + kCommonRequiredHelp + kLoopsSpecificHelp +
                  kCommonOptionalHelp + kLoopsOptionalHelp + kCalibrationHelp);
    {
        cli::ArgumentGroup& required = loops.group("Required arguments");
        cli::ArgumentGroup& optional = loops.group("Optional arguments");
        cli::ArgumentGroup& calibration = loops.group("Calibration arguments");
        add_common(loops, required, optional, calibration);
        required.add({"--loops"})
            .nargs("+")
            .required()
            .metavar("LOOPS")
            .input({"bedgraph", "bedpe"})
            .help("Loop calls, one file or more: the first six columns chrom1 start1 end1 "
                  "chrom2 start2 end2.");
        required.add({"--outFilePrefix", "-o"})
            .required()
            .metavar("PREFIX")
            .output({"tsv"}, "prefix")
            .help("Prefix of the output file <prefix>_loops.tsv.");
        optional.add({"--peakWidth"})
            .type("int")
            .default_value(1)
            .metavar("BINS")
            .help("Half width of the peak square in bins.");
        optional.add({"--windowSize"})
            .type("int")
            .default_value(5)
            .metavar("BINS")
            .help("Half width of the background square in bins, at least --peakWidth + 2.");
        optional.add({"--threads", "-t"})
            .type("int")
            .default_value(4)
            .help("Worker threads; the result does not depend on the number.");
        calibration.add({"-h", "--help"})
            .action(cli::Action::Help)
            .help("show this help message and exit");
    }

    cli::Parser& compartments = parser.add_subcommand(
        "compartments", "Differential compartment preference of every bin.");
    compartments.set_usage(kCompartmentsUsage)
        .set_help(std::string(kCompartmentsHelpHead) + kCommonRequiredHelp +
                  kCompartmentsSpecificHelp + kCommonOptionalHelp + kCompartmentsOptionalHelp +
                  kCalibrationHelp);
    {
        cli::ArgumentGroup& required = compartments.group("Required arguments");
        cli::ArgumentGroup& optional = compartments.group("Optional arguments");
        cli::ArgumentGroup& calibration = compartments.group("Calibration arguments");
        add_common(compartments, required, optional, calibration);
        required.add({"--gcTrack"})
            .required()
            .metavar("BEDGRAPH")
            .input({"bedgraph"})
            .help("GC content per interval: chrom, start, end, value.");
        required.add({"--outFilePrefix", "-o"})
            .required()
            .metavar("PREFIX")
            .output({"tsv"}, "prefix")
            .help("Prefix of the output file <prefix>_compartments.tsv.");
        optional.add({"--minDistance"})
            .type("int")
            .default_value(200000)
            .metavar("BP")
            .help("Least distance of the contacts that are counted.");
        optional.add({"--threads", "-t"})
            .type("int")
            .default_value(4)
            .help("Worker threads; the result does not depend on the number.");
        calibration.add({"-h", "--help"})
            .action(cli::Action::Help)
            .help("show this help message and exit");
    }

    const cli::Namespace ns = parser.parse(argc, argv);
    Arguments args;
    args.command = ns.command();
    args.condition_a = ns.strs("conditionA");
    args.condition_b = ns.strs("conditionB");
    args.blocks = ns.strs("blocks");
    args.chromosomes = ns.strs("chromosomes");
    args.exploratory = ns.flag("exploratory");
    args.fdr = ns.real("fdr");
    args.min_fold_change = ns.real("minFoldChange");
    if (ns.given("filterThreshold")) {
        const std::vector<double> limits = ns.reals("filterThreshold");
        args.filter_lower = limits.at(0);
        args.filter_upper = limits.at(1);
    }
    args.threads = ns.integer("threads");
    args.split_seed = ns.opt_integer("splitReplicates");
    args.plant_regions = ns.opt_str("plantRegions");
    args.plant_fold = ns.real("plantFold");
    args.plant_seed = ns.integer("plantSeed");
    args.prefix = ns.str("outFilePrefix");
    if (args.command == "tads") {
        args.domains = ns.str("domains");
        args.boundary_window = ns.integer("boundaryWindow");
    }
    if (args.command == "loops") {
        args.loops = ns.strs("loops");
        args.peak_width = ns.integer("peakWidth");
        args.window_size = ns.integer("windowSize");
    }
    if (args.command == "compartments") {
        args.gc_track = ns.str("gcTrack");
        args.min_distance = ns.integer("minDistance");
    }
    return args;
}

// --------------------------------------------------------------------------

std::string format(double value, const char* spec) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0.0 ? "inf" : "-inf";
    }
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), spec, value);
    return buffer;
}

std::string num(double value) { return format(value, "%.6g"); }
std::string pval(double value) { return format(value, "%.6e"); }

bool write_file(const std::string& path, const std::string& content) {
    std::FILE* handle = std::fopen(path.c_str(), "wb");
    if (handle == nullptr) {
        return false;
    }
    const bool written = std::fwrite(content.data(), 1, content.size(), handle) == content.size();
    return std::fclose(handle) == 0 && written;
}

std::int64_t parse_int(const std::string& text, const std::string& what) {
    std::int64_t value = 0;
    if (!hicx::cli::python_int(text, &value)) {
        throw std::runtime_error(what + ": '" + text + "' is not an integer");
    }
    return value;
}

struct SampleSpec {
    std::string path;
    // 0 condition A, 1 condition B.
    int condition = 0;
    // Position within its condition.
    std::size_t index = 0;
};

struct Inputs {
    std::vector<SampleSpec> specs;
    std::vector<hicx::CoolFile> files;
    std::vector<std::string> labels;
    diff::DesignMatrices design;
    std::int64_t bin_size = 0;
    std::vector<std::string> chromosomes;
    std::map<std::string, std::int64_t> chrom_length;
};

// Everything that can be refused without reading a pixel.
void validate(const Arguments& args) {
    const std::size_t a = args.condition_a.size();
    const std::size_t b = args.condition_b.size();
    if ((a < 2 || b < 2) && !args.exploratory) {
        throw std::runtime_error(
            "each condition needs at least two replicates, got " + std::to_string(a) +
            " in condition A and " + std::to_string(b) +
            " in condition B. Without replicates a difference cannot be told apart from "
            "replicate variation, which is where the false positives of a one against one "
            "comparison come from (cpp/PLAN.md 9.7). Give --exploratory to run anyway: the "
            "dispersion is then estimated from all samples ignoring the conditions and the "
            "output is labelled exploratory");
    }
    if (!(args.fdr > 0.0 && args.fdr <= 1.0)) {
        throw std::runtime_error("--fdr must lie in (0, 1]");
    }
    if (!(args.min_fold_change >= 1.0)) {
        throw std::runtime_error("--minFoldChange must be at least 1");
    }
    if (!(args.filter_lower < args.filter_upper)) {
        throw std::runtime_error("--filterThreshold needs LOWER < UPPER");
    }
    if (args.threads < 1) {
        throw std::runtime_error("--threads must be at least 1");
    }
    if (args.plant_regions.has_value() && !(args.plant_fold >= 1.0)) {
        throw std::runtime_error("--plantFold must be at least 1");
    }
    if (args.split_seed.has_value() && args.condition_a != args.condition_b) {
        throw std::runtime_error("--splitReplicates needs --conditionA and --conditionB to list "
                                 "the same files in the same order");
    }
    if (args.command == "tads" && args.boundary_window < 1) {
        throw std::runtime_error("--boundaryWindow must be positive");
    }
    if (args.command == "loops" &&
        (args.peak_width < 0 || args.window_size < args.peak_width + 2 || args.window_size > 1000)) {
        throw std::runtime_error("--peakWidth must be at least 0 and --windowSize between "
                                 "--peakWidth + 2 and 1000");
    }
    if (args.command == "compartments" && args.min_distance < 1) {
        throw std::runtime_error("--minDistance must be positive");
    }
    diff::Design design;
    design.condition.assign(a, 0);
    design.condition.insert(design.condition.end(), b, 1);
    design.block = args.blocks;
    (void)diff::build_design(design);
}

Inputs open_inputs(const Arguments& args) {
    Inputs in;
    diff::Design design;
    for (int condition = 0; condition < 2; ++condition) {
        const std::vector<std::string>& paths = condition == 0 ? args.condition_a : args.condition_b;
        for (std::size_t i = 0; i < paths.size(); ++i) {
            SampleSpec spec;
            spec.path = paths[i];
            spec.condition = condition;
            spec.index = i;
            in.specs.push_back(spec);
            in.labels.push_back(std::string(condition == 0 ? "A" : "B") + std::to_string(i + 1));
            design.condition.push_back(condition);
        }
    }
    design.block = args.blocks;
    in.design = diff::build_design(design);

    for (const SampleSpec& spec : in.specs) {
        in.files.emplace_back(spec.path);
    }
    const hicx::CoolFile& reference = in.files.front();
    const auto bin_size_of = [](const hicx::CoolFile& file) -> std::int64_t {
        const hicx::json::Value* value = file.info_value("bin-size");
        return (value != nullptr && value->is_number())
                   ? static_cast<std::int64_t>(value->as_double())
                   : 0;
    };
    in.bin_size = bin_size_of(reference);
    if (in.bin_size <= 0) {
        throw std::runtime_error("'" + reference.filename() + "' has no fixed bin size");
    }
    for (std::size_t s = 1; s < in.files.size(); ++s) {
        if (bin_size_of(in.files[s]) != in.bin_size ||
            in.files[s].chrom_names() != reference.chrom_names() ||
            in.files[s].chrom_lengths() != reference.chrom_lengths()) {
            throw std::runtime_error("'" + in.specs[s].path + "' is not on the bin table of '" +
                                     in.specs[0].path + "'");
        }
    }
    for (std::size_t c = 0; c < reference.chrom_names().size(); ++c) {
        in.chrom_length[reference.chrom_names()[c]] = reference.chrom_lengths()[c];
    }
    in.chromosomes = args.chromosomes.empty() ? reference.chrom_names() : args.chromosomes;
    for (const std::string& chrom : in.chromosomes) {
        if (in.chrom_length.find(chrom) == in.chrom_length.end()) {
            throw std::runtime_error("chromosome '" + chrom + "' is not in '" + in.specs[0].path +
                                     "'");
        }
    }
    return in;
}

// A planted rectangle of pixels in chromosome-local bins: rows [row_first,
// row_last) against columns [col_first, col_last), either orientation.
struct PlantRegion {
    std::string chrom;
    std::int64_t row_first = 0;
    std::int64_t row_last = 0;
    std::int64_t col_first = 0;
    std::int64_t col_last = 0;
    int condition = 0;
    // A BED row (a region), as opposed to a BEDPE rectangle.
    bool region = false;

    [[nodiscard]] bool covers(std::int64_t i, std::int64_t j) const {
        return (i >= row_first && i < row_last && j >= col_first && j < col_last) ||
               (j >= row_first && j < row_last && i >= col_first && i < col_last);
    }
};

// BED (chrom start end [A|B]): the square of the region's contacts with
// itself. BEDPE (chrom1 start1 end1 chrom2 start2 end2 [A|B]), cis only: the
// rectangle of the first region's contacts with the second.
std::vector<PlantRegion> read_plant_regions(const std::string& path, std::int64_t bin_size) {
    const auto condition_of = [&](const std::string& text) {
        if (text == "A") {
            return 0;
        }
        if (text == "B") {
            return 1;
        }
        throw std::runtime_error("'" + path + "': the condition column must be A or B, not '" +
                                 text + "'");
    };
    std::vector<PlantRegion> regions;
    for (const std::vector<std::string>& row : diffc::read_fields(path)) {
        PlantRegion region;
        region.chrom = row[0];
        if (row.size() == 3 || row.size() == 4) {
            region.row_first = parse_int(row[1], path) / bin_size;
            region.row_last = (parse_int(row[2], path) + bin_size - 1) / bin_size;
            region.col_first = region.row_first;
            region.col_last = region.row_last;
            region.condition = row.size() == 4 ? condition_of(row[3]) : 0;
            region.region = true;
        } else if (row.size() == 6 || row.size() == 7) {
            if (row[3] != row[0]) {
                throw std::runtime_error("'" + path + "': a planted rectangle must be cis, got " +
                                         row[0] + " and " + row[3]);
            }
            region.row_first = parse_int(row[1], path) / bin_size;
            region.row_last = (parse_int(row[2], path) + bin_size - 1) / bin_size;
            region.col_first = parse_int(row[4], path) / bin_size;
            region.col_last = (parse_int(row[5], path) + bin_size - 1) / bin_size;
            region.condition = row.size() == 7 ? condition_of(row[6]) : 0;
        } else {
            throw std::runtime_error("'" + path + "': a planted region is BED (3 or 4 columns) "
                                     "or BEDPE (6 or 7 columns), got " +
                                     std::to_string(row.size()) + " columns");
        }
        regions.push_back(region);
    }
    return regions;
}

// Every sample's band of one chromosome. With --splitReplicates, each file is
// read once and its two halves go to condition A and condition B.
std::vector<diffc::ChromosomeData> load_samples(const Arguments& args, const Inputs& in,
                                                std::int64_t first, std::int64_t last,
                                                std::int64_t band) {
    const std::size_t samples = in.specs.size();
    std::vector<diffc::ChromosomeData> data(samples);
    const std::size_t in_a = args.condition_a.size();
    const auto request = [&](std::size_t s) {
        diffc::LoadRequest r;
        r.path = in.specs[s].path;
        r.index = in.specs[s].index;
        if (args.split_seed.has_value()) {
            r.split = true;
            r.split_seed = static_cast<std::uint64_t>(*args.split_seed);
        }
        return r;
    };
    if (args.split_seed.has_value()) {
        for (std::size_t i = 0; i < in_a; ++i) {
            std::vector<diffc::ChromosomeData> halves =
                diffc::load_chromosome(in.files[i], first, last, band, request(i));
            data[i] = std::move(halves[0]);
            data[in_a + i] = std::move(halves[1]);
        }
    } else {
        for (std::size_t s = 0; s < samples; ++s) {
            data[s] = std::move(diffc::load_chromosome(in.files[s], first, last, band, request(s))[0]);
        }
    }
    return data;
}

// The union of every sample's invalid bins.
std::vector<char> shared_invalid(const Arguments& args,
                                 const std::vector<diffc::ChromosomeData>& data) {
    std::vector<char> invalid(data.front().coverage.size(), 0);
    for (const diffc::ChromosomeData& sample : data) {
        const std::vector<char> outliers =
            diffc::coverage_outliers(sample.coverage, args.filter_lower, args.filter_upper);
        for (std::size_t i = 0; i < invalid.size(); ++i) {
            invalid[i] = static_cast<char>(invalid[i] | outliers[i]);
        }
    }
    return invalid;
}

// Thins the planted regions of this chromosome in the samples of their
// condition. The stream of a sample is its condition and index, so the two
// halves of a split file draw independently.
void apply_plants(const Arguments& args, const Inputs& in, const std::string& chrom,
                  const std::vector<PlantRegion>& plants, std::vector<diffc::ChromosomeData>& data,
                  unsigned threads) {
    std::vector<PlantRegion> here;
    for (const PlantRegion& region : plants) {
        if (region.chrom == chrom) {
            here.push_back(region);
        }
    }
    if (here.empty()) {
        return;
    }
    hicx::parallel_for(data.size(), threads, [&](std::size_t s) {
        const int condition = in.specs[s].condition;
        diffc::thin_pixels(
            data[s],
            [&](std::int64_t i, std::int64_t j) {
                for (const PlantRegion& region : here) {
                    if (region.condition == condition && region.covers(i, j)) {
                        return true;
                    }
                }
                return false;
            },
            1.0 / args.plant_fold, static_cast<std::uint64_t>(args.plant_seed),
            static_cast<std::uint64_t>(condition) * 65536U + in.specs[s].index);
    });
}

std::string header_common(const Arguments& args, const Inputs& in, const std::string& what) {
    std::string text = "# hicDifferentialAnalysis " + args.command + ", HiCExplorer " +
                       std::string(hicx::kVersion) + ": " + what + "\n";
    text += "# model: negative binomial GLM with offsets for library size and distance decay, "
            "trended NB dispersion, quasi-likelihood dispersion with empirical Bayes moderation, "
            "TREAT likelihood ratio test, Benjamini-Hochberg FDR; bins filtered in any sample "
            "are masked in all\n";
    text += "# samples:";
    for (std::size_t s = 0; s < in.specs.size(); ++s) {
        text += " " + in.labels[s] + "=" + in.specs[s].path;
    }
    text += "\n";
    if (!args.blocks.empty()) {
        text += "# blocks:";
        for (std::size_t s = 0; s < in.specs.size(); ++s) {
            text += " " + in.labels[s] + "=" + args.blocks[s];
        }
        text += "\n";
    }
    text += "# minimum fold change " + num(args.min_fold_change) + ", FDR " + num(args.fdr) +
            ", filterThreshold " + num(args.filter_lower) + " " + num(args.filter_upper) + "\n";
    if (in.design.exploratory) {
        text += "# EXPLORATORY: a condition has a single sample, so the dispersion is estimated "
                "from all samples ignoring the conditions; this is not a replicated test\n";
    } else {
        text += "# replicated: " + std::to_string(args.condition_a.size()) + " samples in "
                "condition A, " + std::to_string(args.condition_b.size()) + " in condition B\n";
    }
    if (args.split_seed.has_value()) {
        text += "# calibration: replicate split with seed " + std::to_string(*args.split_seed) +
                "\n";
    }
    if (args.plant_regions.has_value()) {
        text += "# calibration: planted regions " + *args.plant_regions + ", fold " +
                num(args.plant_fold) + ", seed " + std::to_string(args.plant_seed) + "\n";
    }
    return text;
}

void report_family(const std::string& command, const std::string& name,
                   const diff::FamilyResult& result) {
    std::fprintf(stderr, "%s %s: family %s: %zu units tested, prior df %s, test df %s\n", kProg,
                 command.c_str(), name.c_str(), result.tested, num(result.prior_df).c_str(),
                 num(result.test_df).c_str());
}

// --------------------------------------------------------------------------
// tads

struct Tad {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string name;
    std::int64_t first = 0;  // local bins [first, last)
    std::int64_t last = 0;
    int strata = 0;
    std::vector<double> counts;    // strata x samples
    std::vector<double> expected;  // strata x samples
};

struct Boundary {
    std::string chrom;
    std::size_t left = 0;  // indices into the TAD list
    std::size_t right = 0;
    std::int64_t position = 0;
    std::int64_t window_start = 0;
    std::int64_t window_end = 0;
    std::int64_t window_bins = 0;
    std::vector<double> cross;
    std::vector<double> cross_expected;
    std::vector<double> flank;
    std::vector<double> flank_expected;
};

int stratum_of(std::int64_t distance) {
    return static_cast<int>(std::bit_width(static_cast<std::uint64_t>(distance))) - 1;
}

void run_tads(const Arguments& args) {
    Inputs in = open_inputs(args);
    const std::size_t samples = in.specs.size();
    const auto threads = static_cast<unsigned>(std::min<std::int64_t>(args.threads, 1024));
    const std::int64_t res = in.bin_size;

    std::vector<Tad> tads;
    for (const std::vector<std::string>& row : diffc::read_fields(args.domains)) {
        if (row.size() < 3) {
            throw std::runtime_error("'" + args.domains + "': a domain needs chrom, start and end");
        }
        if (std::find(in.chromosomes.begin(), in.chromosomes.end(), row[0]) ==
            in.chromosomes.end()) {
            continue;
        }
        Tad tad;
        tad.chrom = row[0];
        tad.start = parse_int(row[1], args.domains);
        tad.end = parse_int(row[2], args.domains);
        tad.name = row.size() > 3 ? row[3] : ".";
        const std::int64_t chrom_bins = (in.chrom_length[tad.chrom] + res - 1) / res;
        tad.first = std::clamp<std::int64_t>(tad.start / res, 0, chrom_bins);
        tad.last = std::clamp<std::int64_t>((tad.end + res - 1) / res, 0, chrom_bins);
        const std::int64_t length = tad.last - tad.first;
        tad.strata = length >= 2 ? stratum_of(length - 1) + 1 : 0;
        tad.counts.assign(static_cast<std::size_t>(tad.strata) * samples, 0.0);
        tad.expected.assign(static_cast<std::size_t>(tad.strata) * samples, 0.0);
        tads.push_back(std::move(tad));
    }
    if (tads.empty()) {
        throw std::runtime_error("'" + args.domains + "' has no domain on the analysed chromosomes");
    }
    const std::vector<PlantRegion> plants = args.plant_regions.has_value()
                                                ? read_plant_regions(*args.plant_regions, res)
                                                : std::vector<PlantRegion>();
    const std::int64_t window = std::max<std::int64_t>(1, (args.boundary_window + res / 2) / res);
    std::int64_t longest = 1;
    for (const Tad& tad : tads) {
        longest = std::max(longest, tad.last - tad.first);
    }
    const std::int64_t band = std::max(longest - 1, 2 * window - 1);

    std::vector<Boundary> boundaries;
    std::size_t masked_bins = 0;
    for (const std::string& chrom : in.chromosomes) {
        std::vector<std::size_t> here;
        for (std::size_t t = 0; t < tads.size(); ++t) {
            if (tads[t].chrom == chrom) {
                here.push_back(t);
            }
        }
        if (here.empty()) {
            continue;
        }
        const auto [first, last] = in.files.front().extent(chrom);
        std::vector<diffc::ChromosomeData> data = load_samples(args, in, first, last, band);
        // The mask is fixed before any planted difference is applied.
        const std::vector<char> invalid = shared_invalid(args, data);
        masked_bins += static_cast<std::size_t>(std::count(invalid.begin(), invalid.end(), 1));
        const auto valid = [&](std::int64_t i) { return invalid[static_cast<std::size_t>(i)] == 0; };
        apply_plants(args, in, chrom, plants, data, threads);
        std::vector<std::vector<double>> decay(samples);
        hicx::parallel_for(samples, threads, [&](std::size_t s) {
            decay[s] = diffc::distance_decay(data[s], invalid, band);
        });

        hicx::parallel_for(here.size(), threads, [&](std::size_t h) {
            Tad& tad = tads[here[h]];
            const std::int64_t length = tad.last - tad.first;
            if (tad.strata == 0) {
                return;
            }
            std::vector<double> pairs(static_cast<std::size_t>(length), 0.0);
            for (std::int64_t d = 1; d < length; ++d) {
                for (std::int64_t i = tad.first; i + d < tad.last; ++i) {
                    if (valid(i) && valid(i + d)) {
                        pairs[static_cast<std::size_t>(d)] += 1.0;
                    }
                }
            }
            for (std::size_t s = 0; s < samples; ++s) {
                for (std::int64_t i = tad.first; i < tad.last; ++i) {
                    if (!valid(i)) {
                        continue;
                    }
                    const auto [k0, k1] = data[s].row_range(i, i + 1, tad.last);
                    for (std::size_t k = k0; k < k1; ++k) {
                        const std::int64_t j = data[s].col[k];
                        if (valid(j)) {
                            tad.counts[static_cast<std::size_t>(stratum_of(j - i)) * samples + s] +=
                                data[s].count[k];
                        }
                    }
                }
                for (std::int64_t d = 1; d < length; ++d) {
                    tad.expected[static_cast<std::size_t>(stratum_of(d)) * samples + s] +=
                        pairs[static_cast<std::size_t>(d)] * decay[s][static_cast<std::size_t>(d)];
                }
            }
        });

        std::vector<Boundary> found;
        for (std::size_t h = 0; h + 1 < here.size(); ++h) {
            const Tad& left = tads[here[h]];
            const Tad& right = tads[here[h + 1]];
            if (left.last != right.first || left.last <= left.first || right.last <= right.first) {
                continue;
            }
            Boundary boundary;
            boundary.chrom = chrom;
            boundary.left = here[h];
            boundary.right = here[h + 1];
            boundary.window_bins =
                std::min({window, left.last - left.first, right.last - right.first});
            const std::int64_t b = left.last;
            boundary.position = std::min(b * res, in.chrom_length[chrom]);
            boundary.window_start = (b - boundary.window_bins) * res;
            boundary.window_end = std::min((b + boundary.window_bins) * res, in.chrom_length[chrom]);
            found.push_back(std::move(boundary));
        }
        hicx::parallel_for(found.size(), threads, [&](std::size_t f) {
            Boundary& boundary = found[f];
            const std::int64_t w = boundary.window_bins;
            const std::int64_t b = tads[boundary.left].last;
            boundary.cross.assign(samples, 0.0);
            boundary.cross_expected.assign(samples, 0.0);
            boundary.flank.assign(samples, 0.0);
            boundary.flank_expected.assign(samples, 0.0);
            std::vector<double> cross_pairs(static_cast<std::size_t>(2 * w), 0.0);
            std::vector<double> flank_pairs(static_cast<std::size_t>(w), 0.0);
            for (std::int64_t i = b - w; i < b; ++i) {
                for (std::int64_t j = b; j < b + w; ++j) {
                    if (valid(i) && valid(j)) {
                        cross_pairs[static_cast<std::size_t>(j - i)] += 1.0;
                    }
                }
            }
            for (const std::int64_t p : {b - w, b}) {
                for (std::int64_t i = p; i < p + w; ++i) {
                    for (std::int64_t j = i + 1; j < p + w; ++j) {
                        if (valid(i) && valid(j)) {
                            flank_pairs[static_cast<std::size_t>(j - i)] += 1.0;
                        }
                    }
                }
            }
            for (std::size_t s = 0; s < samples; ++s) {
                for (std::int64_t i = b - w; i < b; ++i) {
                    if (!valid(i)) {
                        continue;
                    }
                    const auto [k0, k1] = data[s].row_range(i, b, b + w);
                    for (std::size_t k = k0; k < k1; ++k) {
                        if (valid(data[s].col[k])) {
                            boundary.cross[s] += data[s].count[k];
                        }
                    }
                }
                for (const std::int64_t p : {b - w, b}) {
                    for (std::int64_t i = p; i < p + w; ++i) {
                        if (!valid(i)) {
                            continue;
                        }
                        const auto [k0, k1] = data[s].row_range(i, i + 1, p + w);
                        for (std::size_t k = k0; k < k1; ++k) {
                            if (valid(data[s].col[k])) {
                                boundary.flank[s] += data[s].count[k];
                            }
                        }
                    }
                }
                for (std::size_t d = 1; d < cross_pairs.size(); ++d) {
                    boundary.cross_expected[s] += cross_pairs[d] * decay[s][d];
                }
                for (std::size_t d = 1; d < flank_pairs.size(); ++d) {
                    boundary.flank_expected[s] += flank_pairs[d] * decay[s][d];
                }
            }
        });
        for (Boundary& boundary : found) {
            boundaries.push_back(std::move(boundary));
        }
    }

    diff::FamilyOptions options;
    options.min_log_fold = std::log(args.min_fold_change);
    options.threads = threads;

    // The totals.
    diff::Family total;
    total.samples = samples;
    total.counts.assign(tads.size() * samples, 0.0);
    total.log_offsets.assign(tads.size() * samples, -std::numeric_limits<double>::infinity());
    for (std::size_t t = 0; t < tads.size(); ++t) {
        for (std::size_t s = 0; s < samples; ++s) {
            double count = 0.0;
            double expected = 0.0;
            for (int k = 0; k < tads[t].strata; ++k) {
                count += tads[t].counts[static_cast<std::size_t>(k) * samples + s];
                expected += tads[t].expected[static_cast<std::size_t>(k) * samples + s];
            }
            total.counts[t * samples + s] = count;
            total.log_offsets[t * samples + s] = expected > 0.0 ? std::log(expected) : kNaN;
        }
    }
    const diff::FamilyResult total_result = diff::test_family(total, in.design, options);
    report_family(args.command, "total", total_result);

    // The strata: stratum k is its own family while it has at least
    // kStratumFamilyMinimum units; the remaining strata form one family.
    int max_strata = 0;
    for (const Tad& tad : tads) {
        max_strata = std::max(max_strata, tad.strata);
    }
    std::vector<std::vector<double>> stratum_p(tads.size());
    for (std::size_t t = 0; t < tads.size(); ++t) {
        stratum_p[t].assign(static_cast<std::size_t>(tads[t].strata), kNaN);
    }
    for (int k = 0; k < max_strata;) {
        std::size_t members = 0;
        for (const Tad& tad : tads) {
            members += tad.strata > k ? 1 : 0;
        }
        const int k_last = members >= kStratumFamilyMinimum ? k + 1 : max_strata;
        std::vector<std::pair<std::size_t, int>> units;
        for (std::size_t t = 0; t < tads.size(); ++t) {
            for (int q = k; q < std::min(k_last, tads[t].strata); ++q) {
                units.emplace_back(t, q);
            }
        }
        diff::Family family;
        family.samples = samples;
        for (const auto& [t, q] : units) {
            for (std::size_t s = 0; s < samples; ++s) {
                const double expected = tads[t].expected[static_cast<std::size_t>(q) * samples + s];
                family.counts.push_back(tads[t].counts[static_cast<std::size_t>(q) * samples + s]);
                family.log_offsets.push_back(expected > 0.0 ? std::log(expected) : kNaN);
            }
        }
        const diff::FamilyResult result = diff::test_family(family, in.design, options);
        report_family(args.command,
                      "stratum " + std::to_string(k) +
                          (k_last > k + 1 ? "-" + std::to_string(k_last - 1) : std::string()),
                      result);
        for (std::size_t u = 0; u < units.size(); ++u) {
            stratum_p[units[u].first][static_cast<std::size_t>(units[u].second)] = result.pvalue[u];
        }
        k = k_last;
    }

    std::vector<double> tad_p(tads.size(), kNaN);
    for (std::size_t t = 0; t < tads.size(); ++t) {
        if (std::isnan(total_result.pvalue[t])) {
            continue;
        }
        std::vector<double> values = stratum_p[t];
        values.push_back(total_result.pvalue[t]);
        tad_p[t] = diff::simes(values);
    }
    const std::vector<double> tad_fdr = hicx::stats::benjamini_hochberg_adjusted(tad_p);

    // Boundaries.
    diff::Family insulation;
    insulation.samples = samples;
    for (const Boundary& boundary : boundaries) {
        double mean_cross = 0.0;
        double mean_flank = 0.0;
        for (std::size_t s = 0; s < samples; ++s) {
            insulation.counts.push_back(boundary.cross[s]);
            const double ok = boundary.flank[s] > 0.0 && boundary.flank_expected[s] > 0.0 &&
                              boundary.cross_expected[s] > 0.0;
            insulation.log_offsets.push_back(
                ok ? std::log(boundary.flank[s] * boundary.cross_expected[s] /
                              boundary.flank_expected[s])
                   : kNaN);
            mean_cross += boundary.cross[s] / static_cast<double>(samples);
            mean_flank += boundary.flank[s] / static_cast<double>(samples);
        }
        insulation.covariate.push_back(-std::log(1.0 / (mean_cross + 0.5) + 1.0 / (mean_flank + 0.5)));
    }
    const diff::FamilyResult boundary_result = diff::test_family(insulation, in.design, options);
    report_family(args.command, "boundaries", boundary_result);
    const std::vector<double> boundary_fdr =
        hicx::stats::benjamini_hochberg_adjusted(boundary_result.pvalue);

    // Output.
    std::string tad_text = header_common(args, in, "differential TADs");
    tad_text += "# pvalueTotal: the test of all intra-TAD contacts; pvalue: Simes' combination of "
                "it with the per distance stratum tests; fdr: Benjamini-Hochberg over TADs; "
                "log2FoldChange: condition B against A of the total; obsExp: observed over "
                "expected intra-TAD contacts per sample\n";
    tad_text += "# masked bins: " + std::to_string(masked_bins) + "\n";
    tad_text += "#chrom\tstart\tend\tname\tstrata\tlog2FoldChange\tpvalueTotal\tpvalue\tfdr\tdifferential";
    for (const std::string& label : in.labels) {
        tad_text += "\tobsExp_" + label;
    }
    tad_text += "\n";
    std::size_t tad_calls = 0;
    for (std::size_t t = 0; t < tads.size(); ++t) {
        const Tad& tad = tads[t];
        const bool called = !std::isnan(tad_fdr[t]) && tad_fdr[t] <= args.fdr;
        tad_calls += called ? 1 : 0;
        tad_text += tad.chrom + "\t" + std::to_string(tad.start) + "\t" + std::to_string(tad.end) +
                    "\t" + tad.name + "\t" + std::to_string(tad.strata) + "\t" +
                    num(total_result.log_fold[t] / std::log(2.0)) + "\t" +
                    pval(total_result.pvalue[t]) + "\t" + pval(tad_p[t]) + "\t" +
                    pval(tad_fdr[t]) + "\t" + (called ? "1" : "0");
        for (std::size_t s = 0; s < samples; ++s) {
            const double expected = std::exp(total.log_offsets[t * samples + s]);
            tad_text += "\t" + num(total.counts[t * samples + s] / expected);
        }
        tad_text += "\n";
    }

    std::string boundary_text = header_common(args, in, "differential TAD boundaries");
    boundary_text += "# the contacts crossing a boundary within the window, tested against the "
                     "contacts inside the two flanking windows; log2FoldChange > 0: more "
                     "crossing contacts (weaker insulation) in condition B; ratio: crossing over "
                     "flanking contacts per sample, both relative to their expectation\n";
    boundary_text += "#chrom\tposition\twindowStart\twindowEnd\tleftTad\trightTad\tlog2FoldChange\tpvalue\tfdr\tdifferential";
    for (const std::string& label : in.labels) {
        boundary_text += "\tratio_" + label;
    }
    boundary_text += "\n";
    std::size_t boundary_calls = 0;
    for (std::size_t u = 0; u < boundaries.size(); ++u) {
        const Boundary& boundary = boundaries[u];
        const bool called = !std::isnan(boundary_fdr[u]) && boundary_fdr[u] <= args.fdr;
        boundary_calls += called ? 1 : 0;
        boundary_text += boundary.chrom + "\t" + std::to_string(boundary.position) + "\t" +
                         std::to_string(boundary.window_start) + "\t" +
                         std::to_string(boundary.window_end) + "\t" + tads[boundary.left].name +
                         "\t" + tads[boundary.right].name + "\t" +
                         num(boundary_result.log_fold[u] / std::log(2.0)) + "\t" +
                         pval(boundary_result.pvalue[u]) + "\t" + pval(boundary_fdr[u]) + "\t" +
                         (called ? "1" : "0");
        for (std::size_t s = 0; s < samples; ++s) {
            const double ratio = (boundary.cross[s] / boundary.cross_expected[s]) /
                                 (boundary.flank[s] / boundary.flank_expected[s]);
            boundary_text += "\t" + num(ratio);
        }
        boundary_text += "\n";
    }

    const std::string tad_path = args.prefix + "_tads.tsv";
    const std::string boundary_path = args.prefix + "_boundaries.tsv";
    if (!write_file(tad_path, tad_text)) {
        throw std::runtime_error("cannot write '" + tad_path + "'");
    }
    if (!write_file(boundary_path, boundary_text)) {
        throw std::runtime_error("cannot write '" + boundary_path + "'");
    }
    std::fprintf(stderr, "%s tads: %zu of %zu TADs and %zu of %zu boundaries differential at FDR %s%s\n",
                 kProg, tad_calls, tads.size(), boundary_calls, boundaries.size(),
                 num(args.fdr).c_str(), in.design.exploratory ? " (EXPLORATORY)" : "");
}

// --------------------------------------------------------------------------
// loops

struct Loop {
    std::string chrom;
    std::int64_t row = 0;  // local bins of the loop pixel, row <= col
    std::int64_t col = 0;
    std::size_t calls = 0;  // input calls merged into this loop
    bool testable = false;
    std::vector<double> peak;
    std::vector<double> peak_expected;
    std::vector<double> background;
    std::vector<double> background_expected;
};

// Calls within one bin of each other in both anchors (Chebyshev distance 1)
// are linked, and every connected group becomes one loop at the rounded mean
// of its distinct positions.
std::vector<Loop> unite_calls(const std::string& chrom,
                              const std::map<std::pair<std::int64_t, std::int64_t>, std::size_t>&
                                  positions) {
    std::vector<std::pair<std::int64_t, std::int64_t>> points;
    for (const auto& [point, count] : positions) {
        points.push_back(point);
    }
    std::vector<std::size_t> parent(points.size());
    for (std::size_t k = 0; k < parent.size(); ++k) {
        parent[k] = k;
    }
    const std::function<std::size_t(std::size_t)> find = [&](std::size_t k) {
        while (parent[k] != k) {
            parent[k] = parent[parent[k]];
            k = parent[k];
        }
        return k;
    };
    for (std::size_t k = 0; k < points.size(); ++k) {
        for (std::int64_t di = -1; di <= 1; ++di) {
            for (std::int64_t dj = -1; dj <= 1; ++dj) {
                const auto it = positions.find({points[k].first + di, points[k].second + dj});
                if (it != positions.end()) {
                    const std::size_t other = static_cast<std::size_t>(
                        std::distance(positions.begin(), it));
                    const std::size_t a = find(k);
                    const std::size_t b = find(other);
                    if (a != b) {
                        parent[std::max(a, b)] = std::min(a, b);
                    }
                }
            }
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t k = 0; k < points.size(); ++k) {
        groups[find(k)].push_back(k);
    }
    std::vector<Loop> loops;
    for (const auto& [root, members] : groups) {
        double sum_row = 0.0;
        double sum_col = 0.0;
        std::size_t calls = 0;
        for (std::size_t k : members) {
            sum_row += static_cast<double>(points[k].first);
            sum_col += static_cast<double>(points[k].second);
            calls += positions.at(points[k]);
        }
        Loop loop;
        loop.chrom = chrom;
        loop.row = static_cast<std::int64_t>(std::floor(sum_row / members.size() + 0.5));
        loop.col = static_cast<std::int64_t>(std::floor(sum_col / members.size() + 0.5));
        loop.calls = calls;
        loops.push_back(std::move(loop));
    }
    std::sort(loops.begin(), loops.end(), [](const Loop& a, const Loop& b) {
        return a.row < b.row || (a.row == b.row && a.col < b.col);
    });
    return loops;
}

void run_loops(const Arguments& args) {
    Inputs in = open_inputs(args);
    const std::size_t samples = in.specs.size();
    const auto threads = static_cast<unsigned>(std::min<std::int64_t>(args.threads, 1024));
    const std::int64_t res = in.bin_size;
    const std::int64_t pw = args.peak_width;
    const std::int64_t wr = args.window_size;

    std::map<std::string, std::map<std::pair<std::int64_t, std::int64_t>, std::size_t>> positions;
    std::size_t input_calls = 0;
    for (const std::string& path : args.loops) {
        for (const std::vector<std::string>& row : diffc::read_fields(path)) {
            if (row.size() < 6) {
                throw std::runtime_error("'" + path + "': a loop needs chrom1 start1 end1 chrom2 "
                                         "start2 end2");
            }
            if (row[0] != row[3] || std::find(in.chromosomes.begin(), in.chromosomes.end(),
                                              row[0]) == in.chromosomes.end()) {
                continue;
            }
            std::int64_t i = parse_int(row[1], path) / res;
            std::int64_t j = parse_int(row[4], path) / res;
            if (i > j) {
                std::swap(i, j);
            }
            ++positions[row[0]][{i, j}];
            ++input_calls;
        }
    }
    const std::vector<PlantRegion> plants = args.plant_regions.has_value()
                                                ? read_plant_regions(*args.plant_regions, res)
                                                : std::vector<PlantRegion>();

    std::vector<Loop> loops;
    for (const std::string& chrom : in.chromosomes) {
        const auto found = positions.find(chrom);
        if (found == positions.end()) {
            continue;
        }
        std::vector<Loop> here = unite_calls(chrom, found->second);
        std::int64_t longest = 0;
        for (const Loop& loop : here) {
            longest = std::max(longest, loop.col - loop.row);
        }
        const std::int64_t band = longest + 2 * wr;
        const auto [first, last] = in.files.front().extent(chrom);
        const std::int64_t bins = last - first;
        std::vector<diffc::ChromosomeData> data = load_samples(args, in, first, last, band);
        const std::vector<char> invalid = shared_invalid(args, data);
        const auto valid = [&](std::int64_t i) { return invalid[static_cast<std::size_t>(i)] == 0; };
        apply_plants(args, in, chrom, plants, data, threads);
        std::vector<std::vector<double>> decay(samples);
        hicx::parallel_for(samples, threads, [&](std::size_t s) {
            decay[s] = diffc::distance_decay(data[s], invalid, band);
        });
        hicx::parallel_for(here.size(), threads, [&](std::size_t h) {
            Loop& loop = here[h];
            loop.peak.assign(samples, 0.0);
            loop.peak_expected.assign(samples, 0.0);
            loop.background.assign(samples, 0.0);
            loop.background_expected.assign(samples, 0.0);
            if (loop.row - wr < 0 || loop.col + wr >= bins || loop.col - loop.row < wr + 3) {
                return;
            }
            std::size_t peak_pixels = 0;
            std::size_t background_pixels = 0;
            for (std::int64_t r = loop.row - wr; r <= loop.row + wr; ++r) {
                if (!valid(r)) {
                    continue;
                }
                for (std::int64_t c = loop.col - wr; c <= loop.col + wr; ++c) {
                    if (!valid(c) || c - r < 2) {
                        continue;
                    }
                    const std::int64_t ring = std::max(std::abs(r - loop.row), std::abs(c - loop.col));
                    if (ring <= pw) {
                        ++peak_pixels;
                    } else if (ring >= pw + 2) {
                        ++background_pixels;
                    } else {
                        continue;
                    }
                    for (std::size_t s = 0; s < samples; ++s) {
                        const double e = decay[s][static_cast<std::size_t>(c - r)];
                        (ring <= pw ? loop.peak_expected : loop.background_expected)[s] += e;
                    }
                }
                for (std::size_t s = 0; s < samples; ++s) {
                    const auto [k0, k1] = data[s].row_range(r, loop.col - wr, loop.col + wr + 1);
                    for (std::size_t k = k0; k < k1; ++k) {
                        const std::int64_t c = data[s].col[k];
                        if (!valid(c) || c - r < 2) {
                            continue;
                        }
                        const std::int64_t ring =
                            std::max(std::abs(r - loop.row), std::abs(c - loop.col));
                        if (ring <= pw) {
                            loop.peak[s] += data[s].count[k];
                        } else if (ring >= pw + 2) {
                            loop.background[s] += data[s].count[k];
                        }
                    }
                }
            }
            loop.testable = peak_pixels >= 3 && background_pixels >= 20;
        });
        for (Loop& loop : here) {
            loops.push_back(std::move(loop));
        }
    }

    diff::Family family;
    family.samples = samples;
    std::vector<std::size_t> unit_of(loops.size(), static_cast<std::size_t>(-1));
    for (std::size_t l = 0; l < loops.size(); ++l) {
        const Loop& loop = loops[l];
        if (!loop.testable) {
            continue;
        }
        unit_of[l] = family.units();
        double mean_peak = 0.0;
        double mean_background = 0.0;
        for (std::size_t s = 0; s < samples; ++s) {
            family.counts.push_back(loop.peak[s]);
            const bool ok = loop.background[s] > 0.0 && loop.background_expected[s] > 0.0 &&
                            loop.peak_expected[s] > 0.0;
            family.log_offsets.push_back(
                ok ? std::log(loop.background[s] * loop.peak_expected[s] /
                              loop.background_expected[s])
                   : kNaN);
            mean_peak += loop.peak[s] / static_cast<double>(samples);
            mean_background += loop.background[s] / static_cast<double>(samples);
        }
        family.covariate.push_back(-std::log(1.0 / (mean_peak + 0.5) + 1.0 / (mean_background + 0.5)));
    }
    diff::FamilyOptions options;
    options.min_log_fold = std::log(args.min_fold_change);
    options.threads = threads;
    const diff::FamilyResult result = diff::test_family(family, in.design, options);
    report_family(args.command, "loops", result);
    const std::vector<double> fdr = hicx::stats::benjamini_hochberg_adjusted(result.pvalue);

    std::string text = header_common(args, in, "differential loops");
    text += "# " + std::to_string(input_calls) + " input calls united into " +
            std::to_string(loops.size()) + " loops; peak square +-" + std::to_string(pw) +
            " bins, background square +-" + std::to_string(wr) +
            " bins; log2FoldChange: condition B against A of the peak contacts relative to the "
            "background; enrichment: observed over expected peak contacts relative to the same "
            "for the background, per sample; nan: not testable (a masked pixel majority, the "
            "matrix edge, or anchors closer than windowSize + 3 bins)\n";
    text += "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tinputCalls\tlog2FoldChange\tpvalue\tfdr\tdifferential";
    for (const std::string& label : in.labels) {
        text += "\tenrichment_" + label;
    }
    text += "\n";
    std::size_t calls = 0;
    std::size_t tested = 0;
    for (std::size_t l = 0; l < loops.size(); ++l) {
        const Loop& loop = loops[l];
        const std::size_t u = unit_of[l];
        const double log_fold = u == static_cast<std::size_t>(-1) ? kNaN : result.log_fold[u];
        const double p = u == static_cast<std::size_t>(-1) ? kNaN : result.pvalue[u];
        const double q = u == static_cast<std::size_t>(-1) ? kNaN : fdr[u];
        const bool called = !std::isnan(q) && q <= args.fdr;
        calls += called ? 1 : 0;
        tested += std::isnan(p) ? 0 : 1;
        const std::int64_t length = in.chrom_length[loop.chrom];
        text += loop.chrom + "\t" + std::to_string(loop.row * res) + "\t" +
                std::to_string(std::min((loop.row + 1) * res, length)) + "\t" + loop.chrom + "\t" +
                std::to_string(loop.col * res) + "\t" +
                std::to_string(std::min((loop.col + 1) * res, length)) + "\t" +
                std::to_string(loop.calls) + "\t" + num(log_fold / std::log(2.0)) + "\t" + pval(p) +
                "\t" + pval(q) + "\t" + (called ? "1" : "0");
        for (std::size_t s = 0; s < samples; ++s) {
            const double enrichment =
                loop.testable ? (loop.peak[s] / loop.peak_expected[s]) /
                                    (loop.background[s] / loop.background_expected[s])
                              : kNaN;
            text += "\t" + num(enrichment);
        }
        text += "\n";
    }
    const std::string path = args.prefix + "_loops.tsv";
    if (!write_file(path, text)) {
        throw std::runtime_error("cannot write '" + path + "'");
    }
    std::fprintf(stderr, "%s loops: %zu of %zu tested loops (%zu united) differential at FDR %s%s\n",
                 kProg, calls, tested, loops.size(), num(args.fdr).c_str(),
                 in.design.exploratory ? " (EXPLORATORY)" : "");
}

// --------------------------------------------------------------------------
// compartments

// The k leading eigenpairs of a symmetric n x n matrix (row-major), by
// subspace iteration on k + 2 vectors from a fixed pseudo-random start, with
// modified Gram-Schmidt and a final Rayleigh-Ritz step solved by cyclic
// Jacobi. The matrix products reduce every row sequentially and the rows are
// independent, so the result does not depend on the thread count, which a
// threaded LAPACK would not guarantee.
std::vector<std::vector<double>> leading_eigenvectors(const std::vector<double>& a, std::size_t n,
                                                      std::size_t k, unsigned threads) {
    const std::size_t m = std::min(n, k + 2);
    std::vector<std::vector<double>> q(m, std::vector<double>(n));
    for (std::size_t t = 0; t < m; ++t) {
        for (std::size_t i = 0; i < n; ++i) {
            diffc::CounterRng rng({0x45494745ULL, t, i});
            q[t][i] = rng.uniform() - 0.5;
        }
    }
    const auto orthonormalise = [&](std::vector<std::vector<double>>& v) {
        for (std::size_t t = 0; t < v.size(); ++t) {
            for (std::size_t u = 0; u < t; ++u) {
                double dot = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    dot += v[t][i] * v[u][i];
                }
                for (std::size_t i = 0; i < n; ++i) {
                    v[t][i] -= dot * v[u][i];
                }
            }
            double norm = 0.0;
            for (double x : v[t]) {
                norm += x * x;
            }
            norm = std::sqrt(norm);
            for (double& x : v[t]) {
                x = norm > 0.0 ? x / norm : 0.0;
            }
        }
    };
    const auto multiply = [&](const std::vector<std::vector<double>>& v) {
        std::vector<std::vector<double>> z(v.size(), std::vector<double>(n, 0.0));
        hicx::parallel_for(n, threads, [&](std::size_t i) {
            const double* row = a.data() + i * n;
            for (std::size_t t = 0; t < v.size(); ++t) {
                double sum = 0.0;
                for (std::size_t j = 0; j < n; ++j) {
                    sum += row[j] * v[t][j];
                }
                z[t][i] = sum;
            }
        });
        return z;
    };
    orthonormalise(q);
    for (int iteration = 0; iteration < 500; ++iteration) {
        std::vector<std::vector<double>> z = multiply(q);
        orthonormalise(z);
        double change = 0.0;
        for (std::size_t t = 0; t < k && t < m; ++t) {
            double dot = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                dot += z[t][i] * q[t][i];
            }
            change = std::max(change, 1.0 - std::abs(dot));
        }
        q.swap(z);
        if (change < 1e-12) {
            break;
        }
    }
    // Rayleigh-Ritz.
    const std::vector<std::vector<double>> aq = multiply(q);
    std::vector<double> h(m * m);
    for (std::size_t s = 0; s < m; ++s) {
        for (std::size_t t = 0; t < m; ++t) {
            double dot = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                dot += q[s][i] * aq[t][i];
            }
            h[s * m + t] = dot;
        }
    }
    for (std::size_t s = 0; s < m; ++s) {
        for (std::size_t t = s + 1; t < m; ++t) {
            h[s * m + t] = h[t * m + s] = 0.5 * (h[s * m + t] + h[t * m + s]);
        }
    }
    std::vector<double> u(m * m, 0.0);
    for (std::size_t s = 0; s < m; ++s) {
        u[s * m + s] = 1.0;
    }
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t s = 0; s < m; ++s) {
            for (std::size_t t = s + 1; t < m; ++t) {
                off += h[s * m + t] * h[s * m + t];
            }
        }
        if (off < 1e-30) {
            break;
        }
        for (std::size_t p = 0; p < m; ++p) {
            for (std::size_t r = p + 1; r < m; ++r) {
                const double apr = h[p * m + r];
                if (std::abs(apr) < 1e-300) {
                    continue;
                }
                const double theta = (h[r * m + r] - h[p * m + p]) / (2.0 * apr);
                const double tangent = (theta >= 0.0 ? 1.0 : -1.0) /
                                       (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(tangent * tangent + 1.0);
                const double s = tangent * c;
                for (std::size_t x = 0; x < m; ++x) {
                    const double hxp = h[x * m + p];
                    const double hxr = h[x * m + r];
                    h[x * m + p] = c * hxp - s * hxr;
                    h[x * m + r] = s * hxp + c * hxr;
                }
                for (std::size_t x = 0; x < m; ++x) {
                    const double hpx = h[p * m + x];
                    const double hrx = h[r * m + x];
                    h[p * m + x] = c * hpx - s * hrx;
                    h[r * m + x] = s * hpx + c * hrx;
                }
                for (std::size_t x = 0; x < m; ++x) {
                    const double uxp = u[x * m + p];
                    const double uxr = u[x * m + r];
                    u[x * m + p] = c * uxp - s * uxr;
                    u[x * m + r] = s * uxp + c * uxr;
                }
            }
        }
    }
    std::vector<std::size_t> order(m);
    for (std::size_t s = 0; s < m; ++s) {
        order[s] = s;
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t x, std::size_t y) { return h[x * m + x] > h[y * m + y]; });
    std::vector<std::vector<double>> vectors;
    for (std::size_t e = 0; e < std::min(k, m); ++e) {
        std::vector<double> v(n, 0.0);
        for (std::size_t t = 0; t < m; ++t) {
            const double w = u[t * m + order[e]];
            for (std::size_t i = 0; i < n; ++i) {
                v[i] += w * q[t][i];
            }
        }
        vectors.push_back(std::move(v));
    }
    return vectors;
}

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
    double n = 0.0;
    double mx = 0.0;
    double my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(x[i]) && std::isfinite(y[i])) {
            n += 1.0;
            mx += x[i];
            my += y[i];
        }
    }
    if (n < 3.0) {
        return kNaN;
    }
    mx /= n;
    my /= n;
    double sxy = 0.0;
    double sxx = 0.0;
    double syy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(x[i]) && std::isfinite(y[i])) {
            sxy += (x[i] - mx) * (y[i] - my);
            sxx += (x[i] - mx) * (x[i] - mx);
            syy += (y[i] - my) * (y[i] - my);
        }
    }
    return (sxx > 0.0 && syy > 0.0) ? sxy / std::sqrt(sxx * syy) : kNaN;
}

void run_compartments(const Arguments& args) {
    Inputs in = open_inputs(args);
    const std::size_t samples = in.specs.size();
    const auto threads = static_cast<unsigned>(std::min<std::int64_t>(args.threads, 1024));
    const std::int64_t res = in.bin_size;
    const std::int64_t min_bins = std::max<std::int64_t>(1, (args.min_distance + res / 2) / res);

    // GC per bin: the overlap-weighted mean of the track.
    std::map<std::string, std::pair<std::vector<double>, std::vector<double>>> gc_sums;
    for (const std::string& chrom : in.chromosomes) {
        const std::size_t bins = static_cast<std::size_t>((in.chrom_length[chrom] + res - 1) / res);
        gc_sums[chrom] = {std::vector<double>(bins, 0.0), std::vector<double>(bins, 0.0)};
    }
    for (const std::vector<std::string>& row : diffc::read_fields(args.gc_track)) {
        if (row.size() < 4) {
            throw std::runtime_error("'" + args.gc_track + "': a GC interval needs chrom, start, "
                                     "end and value");
        }
        const auto it = gc_sums.find(row[0]);
        if (it == gc_sums.end()) {
            continue;
        }
        double value = 0.0;
        if (!hicx::cli::python_float(row[3], &value) || !std::isfinite(value)) {
            continue;
        }
        const std::int64_t start = parse_int(row[1], args.gc_track);
        const std::int64_t end = parse_int(row[2], args.gc_track);
        auto& [sums, weights] = it->second;
        for (std::int64_t b = std::max<std::int64_t>(0, start / res);
             b < static_cast<std::int64_t>(sums.size()) && b * res < end; ++b) {
            const double overlap = static_cast<double>(std::min(end, (b + 1) * res) -
                                                       std::max(start, b * res));
            if (overlap > 0.0) {
                sums[static_cast<std::size_t>(b)] += overlap * value;
                weights[static_cast<std::size_t>(b)] += overlap;
            }
        }
    }
    const std::vector<PlantRegion> plants = args.plant_regions.has_value()
                                                ? read_plant_regions(*args.plant_regions, res)
                                                : std::vector<PlantRegion>();

    struct Bin {
        std::string chrom;
        std::int64_t index = 0;
        double e1 = kNaN;
        double gc = kNaN;
        std::vector<double> with_a;
        std::vector<double> with_a_expected;
        std::vector<double> with_b;
        std::vector<double> with_b_expected;
    };
    std::vector<Bin> units;
    std::string orientation_notes;
    for (const std::string& chrom : in.chromosomes) {
        const auto [first, last] = in.files.front().extent(chrom);
        const std::int64_t bins = last - first;
        if (bins < 3) {
            continue;
        }
        std::vector<diffc::ChromosomeData> data = load_samples(args, in, first, last, bins);
        const std::vector<char> invalid = shared_invalid(args, data);
        std::vector<std::size_t> valid_bins;
        for (std::int64_t i = 0; i < bins; ++i) {
            if (invalid[static_cast<std::size_t>(i)] == 0) {
                valid_bins.push_back(static_cast<std::size_t>(i));
            }
        }
        const std::size_t nv = valid_bins.size();
        if (nv < 10) {
            continue;
        }
        if (nv > 20000) {
            throw std::runtime_error("chromosome '" + chrom + "' has " + std::to_string(nv) +
                                     " valid bins; the compartment track is dense, use a "
                                     "resolution of 50 kb or coarser");
        }
        // Consensus compartments from the pooled samples, before any plant.
        std::vector<std::vector<double>> decay(samples);
        hicx::parallel_for(samples, threads, [&](std::size_t s) {
            decay[s] = diffc::distance_decay(data[s], invalid, bins);
        });
        std::vector<std::int64_t> position(static_cast<std::size_t>(bins), -1);
        for (std::size_t v = 0; v < nv; ++v) {
            position[valid_bins[v]] = static_cast<std::int64_t>(v);
        }
        std::vector<double> oe(nv * nv, 0.0);
        std::vector<double> pooled_decay(static_cast<std::size_t>(bins) + 1, 0.0);
        for (std::size_t s = 0; s < samples; ++s) {
            for (std::size_t d = 0; d < pooled_decay.size(); ++d) {
                pooled_decay[d] += decay[s][d];
            }
            for (std::int64_t i = 0; i < bins; ++i) {
                const std::int64_t vi = position[static_cast<std::size_t>(i)];
                if (vi < 0) {
                    continue;
                }
                for (std::int64_t k = data[s].row_start[static_cast<std::size_t>(i)];
                     k < data[s].row_start[static_cast<std::size_t>(i) + 1]; ++k) {
                    const std::int64_t vj = position[static_cast<std::size_t>(data[s].col[static_cast<std::size_t>(k)])];
                    if (vj >= 0) {
                        oe[static_cast<std::size_t>(vi) * nv + static_cast<std::size_t>(vj)] +=
                            data[s].count[static_cast<std::size_t>(k)];
                    }
                }
            }
        }
        for (std::size_t x = 0; x < nv; ++x) {
            for (std::size_t y = x; y < nv; ++y) {
                const std::size_t d = valid_bins[y] - valid_bins[x];
                const double value = pooled_decay[d] > 0.0 ? oe[x * nv + y] / pooled_decay[d] : 0.0;
                oe[x * nv + y] = value;
                oe[y * nv + x] = value;
            }
        }
        // Rows standardised, correlation = Z Z^T / nv.
        for (std::size_t x = 0; x < nv; ++x) {
            double mean = 0.0;
            for (std::size_t y = 0; y < nv; ++y) {
                mean += oe[x * nv + y];
            }
            mean /= static_cast<double>(nv);
            double squares = 0.0;
            for (std::size_t y = 0; y < nv; ++y) {
                oe[x * nv + y] -= mean;
                squares += oe[x * nv + y] * oe[x * nv + y];
            }
            const double norm = std::sqrt(squares);
            for (std::size_t y = 0; y < nv; ++y) {
                oe[x * nv + y] = norm > 0.0 ? oe[x * nv + y] / norm : 0.0;
            }
        }
        std::vector<double> correlation(nv * nv, 0.0);
        hicx::parallel_for(nv, threads, [&](std::size_t x) {
            for (std::size_t y = 0; y < nv; ++y) {
                double sum = 0.0;
                for (std::size_t z = 0; z < nv; ++z) {
                    sum += oe[x * nv + z] * oe[y * nv + z];
                }
                correlation[x * nv + y] = sum;
            }
        });
        std::vector<double>().swap(oe);
        const std::vector<std::vector<double>> vectors =
            leading_eigenvectors(correlation, nv, 3, threads);
        std::vector<double>().swap(correlation);
        std::vector<double> gc(nv, kNaN);
        const auto& [gc_sum, gc_weight] = gc_sums[chrom];
        for (std::size_t v = 0; v < nv; ++v) {
            if (valid_bins[v] < gc_weight.size() && gc_weight[valid_bins[v]] > 0.0) {
                gc[v] = gc_sum[valid_bins[v]] / gc_weight[valid_bins[v]];
            }
        }
        std::size_t chosen = 0;
        double best = -1.0;
        double sign = 1.0;
        for (std::size_t e = 0; e < vectors.size(); ++e) {
            const double r = pearson(vectors[e], gc);
            if (std::isfinite(r) && std::abs(r) > best) {
                best = std::abs(r);
                chosen = e;
                sign = r < 0.0 ? -1.0 : 1.0;
            }
        }
        if (best < 0.0) {
            throw std::runtime_error("the GC track has no values on chromosome '" + chrom + "'");
        }
        orientation_notes += " " + chrom + ":eigenvector" + std::to_string(chosen + 1) + ",r=" +
                             num(best);
        std::vector<double> e1(static_cast<std::size_t>(bins), kNaN);
        std::vector<char> in_a(static_cast<std::size_t>(bins), 0);
        std::vector<char> in_b(static_cast<std::size_t>(bins), 0);
        for (std::size_t v = 0; v < nv; ++v) {
            const double value = sign * vectors[chosen][v];
            e1[valid_bins[v]] = value;
            in_a[valid_bins[v]] = value > 0.0 ? 1 : 0;
            in_b[valid_bins[v]] = value < 0.0 ? 1 : 0;
        }

        // Planted differences: a BED region's contacts with the A bins, a
        // BEDPE rectangle as it is; then the decay again, on planted counts.
        std::vector<PlantRegion> here;
        for (const PlantRegion& region : plants) {
            if (region.chrom == chrom) {
                here.push_back(region);
            }
        }
        if (!here.empty()) {
            hicx::parallel_for(samples, threads, [&](std::size_t s) {
                const int condition = in.specs[s].condition;
                diffc::thin_pixels(
                    data[s],
                    [&](std::int64_t i, std::int64_t j) {
                        for (const PlantRegion& region : here) {
                            if (region.condition != condition) {
                                continue;
                            }
                            if (!region.region) {
                                if (region.covers(i, j)) {
                                    return true;
                                }
                                continue;
                            }
                            if (j - i < min_bins) {
                                continue;
                            }
                            const bool i_in = i >= region.row_first && i < region.row_last;
                            const bool j_in = j >= region.row_first && j < region.row_last;
                            const bool i_a = in_a[static_cast<std::size_t>(i)] != 0;
                            const bool j_a = in_a[static_cast<std::size_t>(j)] != 0;
                            // A contact between two region bins of opposite
                            // compartments is a contact with A for one and
                            // with B for the other; thinning it would move
                            // both scores in opposite senses, so it is left.
                            if (i_in && j_in && i_a != j_a) {
                                continue;
                            }
                            if ((i_in && j_a) || (j_in && i_a)) {
                                return true;
                            }
                        }
                        return false;
                    },
                    1.0 / args.plant_fold, static_cast<std::uint64_t>(args.plant_seed),
                    static_cast<std::uint64_t>(condition) * 65536U + in.specs[s].index);
                decay[s] = diffc::distance_decay(data[s], invalid, bins);
            });
        }

        std::vector<Bin> chrom_units(nv);
        for (Bin& unit : chrom_units) {
            // Sized before the samples run in parallel; each writes its own slot.
            unit.with_a.assign(samples, 0.0);
            unit.with_a_expected.assign(samples, 0.0);
            unit.with_b.assign(samples, 0.0);
            unit.with_b_expected.assign(samples, 0.0);
        }
        hicx::parallel_for(samples, threads, [&](std::size_t s) {
            std::vector<double> with_a(static_cast<std::size_t>(bins), 0.0);
            std::vector<double> with_b(static_cast<std::size_t>(bins), 0.0);
            for (std::int64_t i = 0; i < bins; ++i) {
                if (invalid[static_cast<std::size_t>(i)] != 0) {
                    continue;
                }
                const auto [k0, k1] = data[s].row_range(i, i + min_bins, bins);
                for (std::size_t k = k0; k < k1; ++k) {
                    const auto j = static_cast<std::size_t>(data[s].col[k]);
                    if (invalid[j] != 0) {
                        continue;
                    }
                    const double c = data[s].count[k];
                    if (in_a[j] != 0) {
                        with_a[static_cast<std::size_t>(i)] += c;
                    } else if (in_b[j] != 0) {
                        with_b[static_cast<std::size_t>(i)] += c;
                    }
                    if (in_a[static_cast<std::size_t>(i)] != 0) {
                        with_a[j] += c;
                    } else if (in_b[static_cast<std::size_t>(i)] != 0) {
                        with_b[j] += c;
                    }
                }
            }
            for (std::size_t v = 0; v < nv; ++v) {
                const std::size_t i = valid_bins[v];
                double expected_a = 0.0;
                double expected_b = 0.0;
                for (std::size_t w = 0; w < nv; ++w) {
                    const std::size_t j = valid_bins[w];
                    const std::size_t d = i > j ? i - j : j - i;
                    if (static_cast<std::int64_t>(d) < min_bins) {
                        continue;
                    }
                    if (in_a[j] != 0) {
                        expected_a += decay[s][d];
                    } else if (in_b[j] != 0) {
                        expected_b += decay[s][d];
                    }
                }
                Bin& unit = chrom_units[v];
                unit.with_a[s] = with_a[i];
                unit.with_a_expected[s] = expected_a;
                unit.with_b[s] = with_b[i];
                unit.with_b_expected[s] = expected_b;
            }
        });
        for (std::size_t v = 0; v < nv; ++v) {
            Bin& unit = chrom_units[v];
            unit.chrom = chrom;
            unit.index = static_cast<std::int64_t>(valid_bins[v]);
            unit.e1 = e1[valid_bins[v]];
            unit.gc = gc[v];
            units.push_back(std::move(unit));
        }
    }
    if (units.empty()) {
        throw std::runtime_error("no chromosome has enough valid bins");
    }

    diff::Family family;
    family.samples = samples;
    for (const Bin& unit : units) {
        double mean_a = 0.0;
        double mean_b = 0.0;
        for (std::size_t s = 0; s < samples; ++s) {
            family.counts.push_back(unit.with_a[s]);
            const bool ok = unit.with_b[s] > 0.0 && unit.with_b_expected[s] > 0.0 &&
                            unit.with_a_expected[s] > 0.0;
            family.log_offsets.push_back(
                ok ? std::log(unit.with_b[s] * unit.with_a_expected[s] / unit.with_b_expected[s])
                   : kNaN);
            mean_a += unit.with_a[s] / static_cast<double>(samples);
            mean_b += unit.with_b[s] / static_cast<double>(samples);
        }
        family.covariate.push_back(-std::log(1.0 / (mean_a + 0.5) + 1.0 / (mean_b + 0.5)));
    }
    diff::FamilyOptions options;
    options.min_log_fold = std::log(args.min_fold_change);
    options.threads = threads;
    const diff::FamilyResult result = diff::test_family(family, in.design, options);
    report_family(args.command, "compartment bins", result);
    const std::vector<double> fdr = hicx::stats::benjamini_hochberg_adjusted(result.pvalue);

    std::string text = header_common(args, in, "differential compartment preference");
    text += "# consensus compartments from the pooled samples, oriented by GC (" + args.gc_track +
            "):" + orientation_notes + "\n";
    text += "# minimum contact distance " + std::to_string(min_bins) +
            " bins; log2FoldChange: condition B against A of the contacts with A bins relative "
            "to the contacts with B bins; score: log2 of observed over expected contacts with A "
            "over the same with B, per sample\n";
    text += "#chrom\tstart\tend\tE1\tgc\tcompartment\tlog2FoldChange\tpvalue\tfdr\tdifferential";
    for (const std::string& label : in.labels) {
        text += "\tscore_" + label;
    }
    text += "\n";
    std::size_t calls = 0;
    std::size_t tested = 0;
    for (std::size_t u = 0; u < units.size(); ++u) {
        const Bin& unit = units[u];
        const bool called = !std::isnan(fdr[u]) && fdr[u] <= args.fdr;
        calls += called ? 1 : 0;
        tested += std::isnan(result.pvalue[u]) ? 0 : 1;
        text += unit.chrom + "\t" + std::to_string(unit.index * res) + "\t" +
                std::to_string(std::min((unit.index + 1) * res, in.chrom_length[unit.chrom])) +
                "\t" + num(unit.e1) + "\t" + num(unit.gc) + "\t" +
                (unit.e1 > 0.0 ? "A" : (unit.e1 < 0.0 ? "B" : ".")) + "\t" +
                num(result.log_fold[u] / std::log(2.0)) + "\t" + pval(result.pvalue[u]) + "\t" +
                pval(fdr[u]) + "\t" + (called ? "1" : "0");
        for (std::size_t s = 0; s < samples; ++s) {
            const double score = std::log2((unit.with_a[s] / unit.with_a_expected[s]) /
                                           (unit.with_b[s] / unit.with_b_expected[s]));
            text += "\t" + num(score);
        }
        text += "\n";
    }
    const std::string path = args.prefix + "_compartments.tsv";
    if (!write_file(path, text)) {
        throw std::runtime_error("cannot write '" + path + "'");
    }
    std::fprintf(stderr, "%s compartments: %zu of %zu tested bins differential at FDR %s%s\n", kProg,
                 calls, tested, num(args.fdr).c_str(),
                 in.design.exploratory ? " (EXPLORATORY)" : "");
}

}  // namespace

int main(int argc, char** argv) {
    const Arguments args = parse_arguments(argc, argv);
    try {
        validate(args);
        if (args.command == "tads") {
            run_tads(args);
        } else if (args.command == "loops") {
            run_loops(args);
        } else if (args.command == "compartments") {
            run_compartments(args);
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s %s: error: %s\n", kProg, args.command.c_str(), error.what());
        return 1;
    }
    hicx::report_resource_usage(kProg);
    return 0;
}
