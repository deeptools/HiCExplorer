// chicChicagoPlotViewpoint: plots one bait's CHiCAGO viewpoint from
// chicChicagoScores' own output (cpp/PLAN.md 9.15's CHiCAGO pipeline), the
// C++-only counterpart of R Chicago::plotBaits() (no Python HiCExplorer tool
// reads chicChicagoScores' plain-text format at all: chicPlotViewpoint only
// understands chicViewpoint's own HDF5 interaction files).
//
// Reproduces plotBaits' own axes and semantics, read from its R source
// (Chicago::plotBaits, checked directly against the real installed package
// this session): x is distSign (signed distance from the bait, not
// log-scaled, not absolute), y is the raw N column (not the score: score
// only selects a point's colour); NA distSign (trans) rows are dropped, as
// are bait2bait rows by default (removeBait2bait); two significance tiers
// (plevel1 >= plevel2, R's own defaults 5 and 3) colour a point red or blue,
// everything else stays background-coloured; a grey vertical line marks the
// bait itself at x = 0. Optionally (plotBprof-equivalent, gated on
// --backgroundModel being given, since only that file carries the
// dispersion the upper band needs) the Brownian mean and its dashed upper
// 95% band are drawn too: Bmean is already a per-row column of
// chicChicagoScores' own output (it only depends on the bait's s_j, the
// other end's s_i and the fitted distance function, all baked in per row
// already), and the upper band is Bmean + 1.96 * sqrt(Bmean + Bmean^2 /
// dispersion), R's own formula, dispersion read from --backgroundModel's own
// [Dispersion] alpha (the same file chicChicagoScores itself was given).
//
// The figure is drawn by plot/hicexplorer_plot/chicChicagoPlotViewpoint.py
// (hicx::plot::draw), the same C++-computes-data / Python-draws-figure
// bridge every other plotting tool in this project uses (cpp/PLAN.md tier 7,
// option (a)); this is the first genuinely new C++-only tool to add a new
// drawing module rather than port one from the Python original.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "hicx/argparse.hpp"
#include "hicx/chicago.hpp"
#include "hicx/plot_bridge.hpp"
#include "hicx/version.hpp"

namespace {

namespace cli = hicx::cli;
namespace plot = hicx::plot;

const char* const kUsage =
    "usage: chicChicagoPlotViewpoint --scores SCORES --baitID BAITID\n"
    "                                [--baitmap BAITMAP] [--backgroundModel BACKGROUNDMODEL]\n"
    "                                [--range RANGE RANGE] [--plevel1 PLEVEL1]\n"
    "                                [--plevel2 PLEVEL2] [--keepBait2bait]\n"
    "                                [--outFileName OUTFILENAME] [--outputFormat OUTPUTFORMAT]\n"
    "                                [--dpi DPI] [--plotData FILE] [--help] [--version]\n";

const char* const kHelp =
    "\n"
    "chicChicagoPlotViewpoint plots one bait's CHiCAGO viewpoint (chicChicagoScores' own output): "
    "N against signed distance from the bait, points coloured by CHiCAGO\n"
    "score significance, matching R Chicago::plotBaits().\n"
    "\n"
    "This tool has no Python HiCExplorer counterpart; it is C++ only (cpp/PLAN.md 9.15).\n"
    "\n"
    "Required arguments:\n"
    "  --scores SCORES       chicChicagoScores' own output file.\n"
    "  --baitID BAITID       The baitID (chicChicagoScores' own baitID column) to plot.\n"
    "\n"
    "Optional arguments:\n"
    "  --baitmap BAITMAP     CHiCAGO .baitmap file. Without it the plot is titled\n"
    "                        by the bare baitID and bait2bait rows are not removed\n"
    "                        (there is no way to tell which other ends are baits).\n"
    "  --backgroundModel BACKGROUNDMODEL\n"
    "                        chicChicagoBackgroundModel's own output file (the one\n"
    "                        --scores was itself produced from). Without it the\n"
    "                        Brownian mean/upper-band overlay is not drawn (Bmean is\n"
    "                        already a column of --scores, but the upper band needs\n"
    "                        the fitted dispersion, only recorded in this file).\n"
    "  --range RANGE RANGE   Upstream and downstream distance from the bait to plot\n"
    "                        (Default: 1000000 1000000, R's own default maxD).\n"
    "  --plevel1 PLEVEL1     Score threshold for the more significant colour tier\n"
    "                        (Default: 5, R's own default).\n"
    "  --plevel2 PLEVEL2     Score threshold for the less significant colour tier\n"
    "                        (Default: 3, R's own default).\n"
    "  --keepBait2bait       Keep bait2bait rows (Default: removed, R's own\n"
    "                        removeBait2bait = TRUE default). Only has an effect\n"
    "                        with --baitmap.\n"
    "  --outFileName OUTFILENAME, -o OUTFILENAME\n"
    "                        The name of the plot file (Default: chicago_viewpoint.png).\n"
    "  --outputFormat OUTPUTFORMAT, -format OUTPUTFORMAT\n"
    "                        Output format of the plot (Default: png).\n"
    "  --dpi DPI             Resolution for the image, if the output is a raster\n"
    "                        graphics format (Default: 300).\n"
    "  --plotData FILE       Write the data the figure is drawn from as JSON to\n"
    "                        this file and do not draw it.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

// One row of chicChicagoScores' own output
// (baitID otherEndID N distSign Bmean Tmean log_p score).
struct ScoreRow {
    long bait_id = 0;
    long other_end_id = 0;
    double n = 0.0;
    bool has_dist_sign = false;
    long dist_sign = 0;
    double bmean = 0.0;
    double score = 0.0;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) out.push_back(tok);
    return out;
}

std::vector<ScoreRow> read_scores_for_bait(const std::string& path, long bait_id) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open scores file: " + path);
    std::vector<ScoreRow> rows;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_tab(line);
        if (f.size() < 8) continue;
        if (std::stol(f[0]) != bait_id) continue;
        ScoreRow r;
        r.bait_id = bait_id;
        r.other_end_id = std::stol(f[1]);
        r.n = std::stod(f[2]);
        if (f[3] == "NA") {
            r.has_dist_sign = false;
        } else {
            r.has_dist_sign = true;
            r.dist_sign = std::stol(f[3]);
        }
        r.bmean = std::stod(f[4]);
        r.score = std::stod(f[7]);
        rows.push_back(r);
    }
    return rows;
}

// Only [Dispersion] alpha is needed here.
std::optional<double> read_dispersion(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open background model file: " + path);
    std::string line, section;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line[0] == '[') {
            section = line;
            continue;
        }
        if (section != "[Dispersion]") continue;
        auto f = split_tab(line);
        if (f.size() == 2 && f[0] == "alpha") return std::stod(f[1]);
    }
    return std::nullopt;
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser("chicChicagoPlotViewpoint",
                       "Plots one bait's CHiCAGO viewpoint from chicChicagoScores' own output.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--scores"}).required().input({"txt"}).help("chicChicagoScores' own output file.");
    required.add({"--baitID"}).required().type("int").help("The baitID to plot.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--baitmap"}).input({"baitmap", "txt"}).help("CHiCAGO .baitmap file.");
    optional.add({"--backgroundModel"})
        .input({"txt"})
        .help("chicChicagoBackgroundModel's own output file.");
    optional.add({"--range"})
        .type("int")
        .nargs(2)
        .default_value(hicx::json::Value::array(
            {hicx::json::Value::integer(1000000), hicx::json::Value::integer(1000000)}))
        .help("Upstream and downstream distance from the bait to plot.");
    optional.add({"--plevel1"}).type("float").default_value(5.0).help("More significant score threshold.");
    optional.add({"--plevel2"}).type("float").default_value(3.0).help("Less significant score threshold.");
    optional.add({"--keepBait2bait"})
        .action(cli::Action::StoreTrue)
        .help("Keep bait2bait rows (Default: removed). Only has an effect with --baitmap.");
    optional.add({"--outFileName", "-o"})
        .default_value("chicago_viewpoint.png")
        .output({"png", "pdf", "svg"})
        .help("The name of the plot file.");
    optional.add({"--outputFormat", "-format"})
        .default_value("png")
        .help("Output format of the plot.");
    optional.add({"--dpi"}).type("int").default_value(std::int64_t{300}).help("Resolution for the image.");
    optional.add({"--help", "-h"}).action(cli::Action::Help).help("show this help message and exit");
    optional.add({"--version"}).version(std::string("%(prog)s ") + hicx::kVersion);
    optional.add({"--plotData"})
        .metavar("FILE")
        .output({"json"})
        .cpp_only("The data of the figure as JSON, without drawing it (cpp/PLAN.md tier 7).")
        .help("Write the data the figure is drawn from as JSON to this file and do not draw it.");

    const cli::Namespace args = parser.parse(argc, argv);
    if (const int refused = plot::preflight("chicChicagoPlotViewpoint", !args.given("plotData"));
        refused != 0) {
        return refused;
    }

    try {
        const long bait_id = args.integer("baitID");
        const auto range = args.integers("range");
        const long upstream = range[0];
        const long downstream = range[1];
        const double plevel1 = args.real("plevel1");
        const double plevel2 = args.real("plevel2");

        std::unordered_set<long> bait_ids;
        std::string bait_label = "baitID " + std::to_string(bait_id);
        if (args.given("baitmap")) {
            const auto baitmap = hicx::chicago::read_baitmap(args.str("baitmap"));
            for (const auto& b : baitmap) bait_ids.insert(b.id);
            for (const auto& b : baitmap) {
                if (b.id == bait_id) {
                    bait_label = b.name.empty() ? bait_label : b.name;
                    break;
                }
            }
        }
        const bool remove_bait2bait = args.given("baitmap") && !args.flag("keepBait2bait");

        std::optional<double> dispersion;
        if (args.given("backgroundModel")) {
            dispersion = read_dispersion(args.str("backgroundModel"));
        }
        const bool has_background = dispersion.has_value();

        auto rows = read_scores_for_bait(args.str("scores"), bait_id);

        std::vector<double> x, y, score;
        std::vector<std::pair<double, double>> bmean_points;  // (distSign, Bmean), for the overlay
        for (const auto& r : rows) {
            if (!r.has_dist_sign) continue;  // trans: no x position
            if (remove_bait2bait && bait_ids.count(r.other_end_id) > 0) continue;
            const double d = static_cast<double>(r.dist_sign);
            if (d < static_cast<double>(-upstream) || d > static_cast<double>(downstream)) continue;
            x.push_back(d);
            y.push_back(r.n);
            score.push_back(r.score);
            if (has_background) bmean_points.emplace_back(d, r.bmean);
        }

        plot::JsonObject data;
        data.add("plotFile", plot::json_string(args.str("outFileName")));
        data.add("baitLabel", plot::json_string(bait_label));
        data.add("plevel1", plot::json_number(plevel1));
        data.add("plevel2", plot::json_number(plevel2));
        data.add("x", plot::json_numbers(x));
        data.add("y", plot::json_numbers(y));
        data.add("score", plot::json_numbers(score));
        data.add("dpi", plot::json_int(args.integer("dpi")));
        data.add("hasBackground", plot::json_bool(has_background));
        if (has_background) {
            std::sort(bmean_points.begin(), bmean_points.end());
            std::vector<double> bx, by, bu;
            const double alpha = *dispersion;
            for (const auto& [d, bmean] : bmean_points) {
                bx.push_back(d);
                by.push_back(bmean);
                bu.push_back(bmean + 1.96 * std::sqrt(bmean + bmean * bmean / alpha));
            }
            data.add("bmeanX", plot::json_numbers(bx));
            data.add("bmeanY", plot::json_numbers(by));
            data.add("bmeanUpperY", plot::json_numbers(bu));
        }

        return plot::draw("chicChicagoPlotViewpoint", data.str(), args.opt_str("plotData"));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoPlotViewpoint: %s\n", error.what());
        return 1;
    }
}
