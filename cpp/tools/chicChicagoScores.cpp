// chicChicagoScores: applies a fitted CHiCAGO background model
// (chicChicagoBackgroundModel's output) to one .chinput file, writing
// CHiCAGO's own p-value and score per (bait, other-end) interaction
// (Cairns et al. 2016, Genome Biology 17:127; cpp/PLAN.md 9.15).
//
// Mirrors chicViewpoint's role (apply a fitted background model to score
// interactions) but on CHiCAGO's own design and its own getPvals/getScores
// statistics (hicx::chicago::log_pvalue, log_weight, score_from_pvalue),
// exactly reproducing R's own weightedRelative scoring given the same
// Bmean/Tmean/dispersion/eta.bar inputs (validated in cpp/tests/test_chicago.cpp
// against real R Chicago output).
//
// A (bait, other-end) pair whose other end was dropped by the background
// model's tlb pooling (its top trans-count percent, --tlbFilterTopPercent)
// is skipped, exactly as R drops it from cd@x entirely at that step; a bait
// with no fitted s_j (readSample already excluded it, or normaliseBaits
// could not normalise it) is skipped the same way.
//
// No Python HiCExplorer counterpart; C++ only (cpp/PLAN.md 9.15).

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>

#include "hicx/argparse.hpp"
#include "hicx/chicago.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: chicChicagoScores --backgroundModel BACKGROUNDMODEL\n"
    "                         (--chinput CHINPUT | --matrices MATRICES [MATRICES ...])\n"
    "                         --baitmap BAITMAP [--outFileName OUTFILENAME]\n"
    "                         [--minFragLen MINFRAGLEN] [--maxFragLen MAXFRAGLEN]\n"
    "                         [--minNPerBait MINNPERBAIT] [--maxLBrownEst MAXLBROWNEST]\n"
    "                         [--noRemoveAdjacent] [--weightAlpha WEIGHTALPHA]\n"
    "                         [--weightBeta WEIGHTBETA] [--weightGamma WEIGHTGAMMA]\n"
    "                         [--weightDelta WEIGHTDELTA] [--rmap RMAP] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicChicagoScores applies a fitted CHiCAGO background model (chicChicagoBackgroundModel's "
    "output) to one .chinput file, writing CHiCAGO's own log p-value\n"
    "and weighted score per (bait, other-end) interaction.\n"
    "\n"
    "This tool has no Python HiCExplorer counterpart; it is C++ only (cpp/PLAN.md 9.15).\n"
    "\n"
    "Required arguments:\n"
    "  --backgroundModel BACKGROUNDMODEL\n"
    "                        chicChicagoBackgroundModel's output file.\n"
    "  --chinput CHINPUT     One CHiCAGO .chinput interaction count file. Mutually\n"
    "                        exclusive with --matrices; exactly one of the two is\n"
    "                        required.\n"
    "  --matrices MATRICES [MATRICES ...], -m MATRICES [MATRICES ...]\n"
    "                        One or more Hi-C matrices (cool, h5 or .hic) to derive\n"
    "                        the per (bait, other-end) fragment N counts from\n"
    "                        directly, instead of a pre-built .chinput file (see\n"
    "                        chicChicagoBackgroundModel --help for the exact scope:\n"
    "                        cis pairs within --maxLBrownEst only). Several matrices\n"
    "                        are summed per fragment pair. Mutually exclusive with\n"
    "                        --chinput; exactly one of the two is required.\n"
    "  --baitmap BAITMAP     CHiCAGO .baitmap file (needed for eta.bar, the\n"
    "                        distance-weighting normalisation constant).\n"
    "  --rmap RMAP           CHiCAGO .rmap file (needed for eta.bar).\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The name of the scores file (Default: chicago_scores.txt).\n"
    "  --minFragLen MINFRAGLEN\n"
    "                        Minimum other-end fragment length (Default: 150).\n"
    "  --maxFragLen MAXFRAGLEN\n"
    "                        Maximum other-end fragment length (Default: 40000).\n"
    "  --minNPerBait MINNPERBAIT\n"
    "                        Minimum total read count per bait (Default: 250).\n"
    "  --maxLBrownEst MAXLBROWNEST\n"
    "                        Maximum distance the Brownian component was estimated\n"
    "                        over (Default: 1500000). Must match the value used for\n"
    "                        --backgroundModel.\n"
    "  --noRemoveAdjacent    Keep interactions with fragments immediately adjacent\n"
    "                        to their bait (Default: they are removed). Must match\n"
    "                        the value used for --backgroundModel.\n"
    "  --weightAlpha WEIGHTALPHA, --weightBeta WEIGHTBETA, --weightGamma WEIGHTGAMMA,\n"
    "  --weightDelta WEIGHTDELTA\n"
    "                        CHiCAGO's distance-weighting curve parameters\n"
    "                        (Default: R Chicago's own defaults).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct LoadedBackgroundModel {
    hicx::chicago::FilterSettings filters;
    hicx::chicago::DistFunFit dist_fun;
    double dispersion = 0.0;
    std::unordered_map<long, double> s_j;
    std::unordered_map<long, int> pool_of_other_end;
    std::unordered_map<int, double> s_i;
    std::unordered_map<long, int> tblb_of_bait;
    std::map<std::pair<int, int>, double> tmean_by_pool;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) out.push_back(tok);
    return out;
}

LoadedBackgroundModel read_background_model(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open background model file: " + path);
    LoadedBackgroundModel m;
    std::string line, section;
    bool skip_header = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line[0] == '[') {
            section = line;
            skip_header = (section == "[BaitFactors]" || section == "[OtherEndPools]" ||
                           section == "[PoolFactors]" || section == "[BaitPools]" ||
                           section == "[TechnicalNoise]");
            continue;
        }
        auto f = split_tab(line);
        if (skip_header) {
            skip_header = false;
            continue;
        }
        if (section == "[Settings]" && f.size() == 2) {
            if (f[0] == "minFragLen") m.filters.min_frag_len = std::stol(f[1]);
            else if (f[0] == "maxFragLen") m.filters.max_frag_len = std::stol(f[1]);
            else if (f[0] == "minNPerBait") m.filters.min_n_per_bait = std::stol(f[1]);
            else if (f[0] == "maxLBrownEst") m.filters.max_l_brown_est = std::stol(f[1]);
            else if (f[0] == "binsize") m.filters.binsize = std::stol(f[1]);
            else if (f[0] == "removeAdjacent") m.filters.remove_adjacent = (f[1] == "1");
        } else if (section == "[DistFun]" && f.size() == 2) {
            if (f[0] == "cubic0") m.dist_fun.cubic[0] = std::stod(f[1]);
            else if (f[0] == "cubic1") m.dist_fun.cubic[1] = std::stod(f[1]);
            else if (f[0] == "cubic2") m.dist_fun.cubic[2] = std::stod(f[1]);
            else if (f[0] == "cubic3") m.dist_fun.cubic[3] = std::stod(f[1]);
            else if (f[0] == "obsMinLog") m.dist_fun.obs_min_log = std::stod(f[1]);
            else if (f[0] == "obsMaxLog") m.dist_fun.obs_max_log = std::stod(f[1]);
            else if (f[0] == "headA") m.dist_fun.head_coef[0] = std::stod(f[1]);
            else if (f[0] == "headB") m.dist_fun.head_coef[1] = std::stod(f[1]);
            else if (f[0] == "tailA") m.dist_fun.tail_coef[0] = std::stod(f[1]);
            else if (f[0] == "tailB") m.dist_fun.tail_coef[1] = std::stod(f[1]);
        } else if (section == "[Dispersion]" && f.size() == 2 && f[0] == "alpha") {
            m.dispersion = std::stod(f[1]);
        } else if (section == "[BaitFactors]" && f.size() == 2) {
            m.s_j[std::stol(f[0])] = std::stod(f[1]);
        } else if (section == "[OtherEndPools]" && f.size() == 2) {
            m.pool_of_other_end[std::stol(f[0])] = std::stoi(f[1]);
        } else if (section == "[PoolFactors]" && f.size() == 2) {
            m.s_i[std::stoi(f[0])] = std::stod(f[1]);
        } else if (section == "[BaitPools]" && f.size() == 2) {
            m.tblb_of_bait[std::stol(f[0])] = std::stoi(f[1]);
        } else if (section == "[TechnicalNoise]" && f.size() == 3) {
            m.tmean_by_pool[{std::stoi(f[0]), std::stoi(f[1])}] = std::stod(f[2]);
        }
    }
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("chicChicagoScores",
                       "Applies a fitted CHiCAGO background model to one .chinput file.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--backgroundModel"})
        .required()
        .input({"txt"})
        .help("chicChicagoBackgroundModel's output file.");
    cli::MutuallyExclusiveGroup& input_source = parser.mutually_exclusive(required, /*required=*/true);
    input_source.add({"--chinput"}).input({"chinput", "txt"}).help("CHiCAGO .chinput file.");
    input_source.add({"--matrices", "-m"})
        .nargs("+")
        .input({"cool", "h5", "hic"})
        .help("Hi-C matrices to derive per-fragment-pair counts from directly.");
    required.add({"--baitmap"}).required().input({"baitmap", "txt"}).help("CHiCAGO .baitmap file.");
    required.add({"--rmap"}).required().input({"rmap", "txt"}).help("CHiCAGO .rmap file.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("chicago_scores.txt")
        .output({"txt"})
        .help("The name of the scores file.");
    optional.add({"--minFragLen"}).type("int").default_value(150).help("Minimum other-end length.");
    optional.add({"--maxFragLen"}).type("int").default_value(40000).help("Maximum other-end length.");
    optional.add({"--minNPerBait"}).type("int").default_value(250).help("Minimum reads per bait.");
    optional.add({"--maxLBrownEst"})
        .type("int")
        .default_value(1500000)
        .help("Brownian estimation distance cutoff.");
    optional.add({"--noRemoveAdjacent"})
        .action(cli::Action::StoreTrue)
        .help("Keep bait-adjacent interactions.");
    optional.add({"--weightAlpha"}).type("float").default_value(34.1157346557331).help("Weight curve alpha.");
    optional.add({"--weightBeta"}).type("float").default_value(-2.58688050486759).help("Weight curve beta.");
    optional.add({"--weightGamma"}).type("float").default_value(-17.1347845819659).help("Weight curve gamma.");
    optional.add({"--weightDelta"}).type("float").default_value(-7.07609217973722).help("Weight curve delta.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace args = parser.parse(argc, argv);

    try {
        using namespace hicx::chicago;
        LoadedBackgroundModel model = read_background_model(args.str("backgroundModel"));
        FilterSettings fs = model.filters;
        // Command-line filter overrides, when given, take precedence over the
        // model file's recorded settings (so a caller can re-check consistency).
        if (args.given("minFragLen")) fs.min_frag_len = args.integer("minFragLen");
        if (args.given("maxFragLen")) fs.max_frag_len = args.integer("maxFragLen");
        if (args.given("minNPerBait")) fs.min_n_per_bait = args.integer("minNPerBait");
        if (args.given("maxLBrownEst")) fs.max_l_brown_est = args.integer("maxLBrownEst");
        if (args.given("noRemoveAdjacent")) fs.remove_adjacent = false;

        auto rmap = read_rmap(args.str("rmap"));
        auto baitmap = read_baitmap(args.str("baitmap"));
        const std::vector<ChinputRecord> raw = args.given("chinput")
            ? read_chinput(args.str("chinput"))
            : chinput_from_matrices(args.strs("matrices"), rmap, baitmap, fs);
        auto x = read_sample(raw, baitmap, fs);

        WeightSettings w;
        w.alpha = args.real("weightAlpha");
        w.beta = args.real("weightBeta");
        w.gamma = args.real("weightGamma");
        w.delta = args.real("weightDelta");
        const double avg_frag_len = avg_frag_length(rmap);
        const auto nhyp = n_hypotheses(rmap, baitmap);
        const double eta_bar = eta_bar_from_design(rmap, baitmap, w, avg_frag_len, true, nhyp);

        std::ofstream out(args.str("outFileName"), std::ios::binary);
        if (!out) {
            throw std::runtime_error("[Errno 2] No such file or directory: '" +
                                     args.str("outFileName") + "'");
        }
        out.precision(15);
        out << "baitID\totherEndID\tN\tdistSign\tBmean\tTmean\tlog_p\tscore\n";

        for (const auto& c : x) {
            auto sj_it = model.s_j.find(c.bait_id);
            if (sj_it == model.s_j.end()) continue;
            auto oe_it = model.pool_of_other_end.find(c.other_end_id);
            if (oe_it == model.pool_of_other_end.end()) continue;  // dropped by tlb pooling
            const int pool = oe_it->second;
            auto si_it = model.s_i.find(pool);
            const double s_i = (si_it == model.s_i.end()) ? 1.0 : si_it->second;

            double abs_dist = std::numeric_limits<double>::infinity();
            double Bmean = 0.0;
            if (c.has_dist_sign) {
                abs_dist = static_cast<double>(std::labs(c.dist_sign));
                Bmean = sj_it->second * s_i * std::exp(eval_distance_function(model.dist_fun, abs_dist));
            }

            auto tblb_it = model.tblb_of_bait.find(c.bait_id);
            const int tblb = (tblb_it == model.tblb_of_bait.end()) ? -1 : tblb_it->second;
            auto tmean_it = model.tmean_by_pool.find({pool, tblb});
            const double Tmean = (tmean_it == model.tmean_by_pool.end()) ? 0.0 : tmean_it->second;

            const double log_p = log_pvalue(c.N, model.dispersion, Bmean, Tmean);
            const double score = score_from_pvalue(log_p, abs_dist, w, eta_bar);

            out << c.bait_id << "\t" << c.other_end_id << "\t" << c.N << "\t";
            if (c.has_dist_sign) out << c.dist_sign; else out << "NA";
            out << "\t" << Bmean << "\t" << Tmean << "\t" << log_p << "\t" << score << "\n";
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoScores: %s\n", error.what());
        return 1;
    }
    return 0;
}
