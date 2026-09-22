// chicChicagoSignificantInteractions: filters chicChicagoScores' output by
// CHiCAGO's own score threshold (Cairns et al. 2016 call interactions at
// score >= 5 by default; cpp/PLAN.md 9.15), writing the accepted calls with
// their fragment coordinates joined from the .rmap/.baitmap.
//
// Mirrors chicSignificantInteractions' role (turn a scored interaction file
// into accepted calls) but on CHiCAGO's own score, not HiCExplorer's
// negative-binomial one.
//
// No Python HiCExplorer counterpart; C++ only (cpp/PLAN.md 9.15).

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "hicx/argparse.hpp"
#include "hicx/chicago.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: chicChicagoSignificantInteractions --scores SCORES --rmap RMAP\n"
    "                                          --baitmap BAITMAP\n"
    "                                          [--outFileName OUTFILENAME]\n"
    "                                          [--scoreThreshold SCORETHRESHOLD] [--help]\n"
    "                                          [--version]\n";

const char* const kHelp =
    "\n"
    "chicChicagoSignificantInteractions filters chicChicagoScores' output by CHiCAGO's own score "
    "threshold (Default: 5, Cairns et al. 2016's own default), writing\n"
    "the accepted calls with their fragment coordinates.\n"
    "\n"
    "This tool has no Python HiCExplorer counterpart; it is C++ only (cpp/PLAN.md 9.15).\n"
    "\n"
    "Required arguments:\n"
    "  --scores SCORES       chicChicagoScores' output file.\n"
    "  --rmap RMAP           CHiCAGO .rmap file (other-end coordinates).\n"
    "  --baitmap BAITMAP     CHiCAGO .baitmap file (bait coordinates and names).\n"
    "\n"
    "Optional arguments:\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The name of the significant-interactions file\n"
    "                        (Default: chicago_significant_interactions.txt).\n"
    "  --scoreThreshold SCORETHRESHOLD\n"
    "                        Minimum CHiCAGO score to call an interaction\n"
    "                        significant (Default: 5).\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct ScoreRow {
    long bait_id = 0, other_end_id = 0;
    double N = 0, Bmean = 0, Tmean = 0, log_p = 0, score = 0;
    bool has_dist_sign = false;
    long dist_sign = 0;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) out.push_back(tok);
    return out;
}

std::vector<ScoreRow> read_scores(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open scores file: " + path);
    std::vector<ScoreRow> rows;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_tab(line);
        if (f.size() < 8) continue;
        ScoreRow r;
        r.bait_id = std::stol(f[0]);
        r.other_end_id = std::stol(f[1]);
        r.N = std::stod(f[2]);
        if (f[3] == "NA") {
            r.has_dist_sign = false;
        } else {
            r.has_dist_sign = true;
            r.dist_sign = std::stol(f[3]);
        }
        r.Bmean = std::stod(f[4]);
        r.Tmean = std::stod(f[5]);
        r.log_p = std::stod(f[6]);
        r.score = std::stod(f[7]);
        rows.push_back(r);
    }
    return rows;
}

}  // namespace

int main(int argc, char** argv) {
    namespace cli = hicx::cli;
    cli::Parser parser("chicChicagoSignificantInteractions",
                       "Filters chicChicagoScores' output by CHiCAGO's own score threshold.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--scores"}).required().input({"txt"}).help("chicChicagoScores' output file.");
    required.add({"--rmap"}).required().input({"rmap", "txt"}).help("CHiCAGO .rmap file.");
    required.add({"--baitmap"}).required().input({"baitmap", "txt"}).help("CHiCAGO .baitmap file.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--outFileName", "-o"})
        .default_value("chicago_significant_interactions.txt")
        .output({"txt"})
        .help("The name of the significant-interactions file.");
    optional.add({"--scoreThreshold"})
        .type("float")
        .default_value(5.0)
        .help("Minimum CHiCAGO score to call an interaction significant.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);

    const cli::Namespace args = parser.parse(argc, argv);

    try {
        using namespace hicx::chicago;
        auto rmap = read_rmap(args.str("rmap"));
        auto baitmap = read_baitmap(args.str("baitmap"));
        std::unordered_map<long, RmapFragment> rmap_by_id;
        for (auto& r : rmap) rmap_by_id[r.id] = r;
        std::unordered_map<long, BaitmapFragment> baitmap_by_id;
        for (auto& b : baitmap) baitmap_by_id[b.id] = b;

        auto rows = read_scores(args.str("scores"));
        const double threshold = args.real("scoreThreshold");

        std::ofstream out(args.str("outFileName"), std::ios::binary);
        if (!out) {
            throw std::runtime_error("[Errno 2] No such file or directory: '" +
                                     args.str("outFileName") + "'");
        }
        out.precision(15);
        out << "bait_chr\tbait_start\tbait_end\tbait_name\totherEnd_chr\totherEnd_start\t"
               "otherEnd_end\tN\tdistSign\tscore\n";
        for (const auto& r : rows) {
            if (!(r.score >= threshold)) continue;
            auto bait_it = baitmap_by_id.find(r.bait_id);
            auto oe_it = rmap_by_id.find(r.other_end_id);
            if (bait_it == baitmap_by_id.end() || oe_it == rmap_by_id.end()) continue;
            out << bait_it->second.chrom << "\t" << bait_it->second.start << "\t"
                << bait_it->second.end << "\t" << bait_it->second.name << "\t" << oe_it->second.chrom
                << "\t" << oe_it->second.start << "\t" << oe_it->second.end << "\t" << r.N << "\t";
            if (r.has_dist_sign) out << r.dist_sign; else out << "NA";
            out << "\t" << r.score << "\n";
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoSignificantInteractions: %s\n", error.what());
        return 1;
    }
    return 0;
}
