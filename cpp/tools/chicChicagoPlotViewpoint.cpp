// chicChicagoPlotViewpoint: plots one bait's CHiCAGO viewpoint, or every
// bait's interactions in a region, from chicChicagoScores' own output
// (cpp/PLAN.md 9.15's CHiCAGO pipeline). The C++-only counterpart of R
// Chicago::plotBaits() (no Python HiCExplorer tool reads chicChicagoScores'
// plain-text format at all: chicPlotViewpoint only understands
// chicViewpoint's own HDF5 interaction files), plus a second, region-wide
// view and an arc-style rendering with no R precedent to match.
//
// Two mutually exclusive scopes:
//
//   * --baitID: one bait's own viewpoint. x is distSign (signed distance
//     from the bait, not log-scaled, not absolute), matching R's own
//     plotBaits() exactly (checked directly against the real installed
//     package this session).
//   * --region CHROM START END: every bait --baitmap places in that region,
//     each of its own interactions drawn against ABSOLUTE genomic
//     coordinates (the bait's own midpoint plus distSign; no .rmap lookup
//     is needed for the other end's position, distSign already is the
//     real, signed genomic offset a real .chinput row's own coordinates
//     were computed from). A bait2bait pair between two baits both inside
//     the region would otherwise be read once from each bait's own row set
//     and drawn twice; it is deduplicated by its canonical (min id, max id)
//     pair regardless of --keepBait2bait.
//
// Two rendering styles (--style):
//
//   * scatter (--baitID only): R's own plotBaits() plot. y is the raw N
//     column (not the score: score only selects a point's colour); NA
//     distSign (trans) rows are dropped, as are bait2bait rows by default
//     (removeBait2bait); two significance tiers (plevel1 >= plevel2, R's
//     own defaults 5 and 3) colour a point red or blue, everything else
//     stays background-coloured; a grey vertical line marks the bait itself
//     at x = 0. Optionally (plotBprof-equivalent, gated on --backgroundModel
//     being given, since only that file carries the dispersion the upper
//     band needs) the Brownian mean and its dashed upper 95% band are drawn
//     too: Bmean is already a per-row column of chicChicagoScores' own
//     output, and the upper band is Bmean + 1.96 * sqrt(Bmean + Bmean^2 /
//     dispersion), R's own formula.
//   * arcs (--baitID or --region): every kept interaction becomes one arc
//     between the bait's own position and the other end's, height N (log1p
//     scaled, so one very deep interaction does not flatten the rest),
//     coloured by the same three tiers as scatter. No R precedent; this
//     project's own convention (matching the arc-style loop track already
//     used in the HiCExplorer v4 GUI's matrix browser).
//
// --onlySignificant drops every background-tier row (score < --plevel2)
// before any of the above, in both scopes and both styles: a plot of only
// the two significant colour tiers, no background clutter.
//
// --linksFile FILE writes the same kept interactions (whichever scope,
// --onlySignificant and --keepBait2bait already selected) as a
// pyGenomeTracks `links` file (chr1 start1 end1 chr2 start2 end2 score, tab
// separated, checked directly against the real installed
// pygenometracks.tracks.LinksTrack this session), so they can be combined
// with real bigwig, gene and BED region tracks (known enhancers, promoters,
// ...) in one pyGenomeTracks figure with pyGenomeTracks:
// write a tracks.ini with a [chicago links] section (file = FILE, file_type
// = links) alongside [bigwig]/[genes]/[bed] sections for the other tracks,
// then `pyGenomeTracks --tracks tracks.ini --region CHROM:START-END -o out.png`.
// Needs --baitmap (real coordinates for the bait side of every link); the
// other end's real fragment span comes from --rmap when given, else a 1 bp
// placeholder at its computed position (still a valid link for pyGenomeTracks'
// arcs style, which only draws between the two extremities).
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
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
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
    "usage: chicChicagoPlotViewpoint --scores SCORES (--baitID BAITID |\n"
    "                                --region CHROM START END)\n"
    "                                [--style {scatter,arcs}] [--baitmap BAITMAP] [--rmap RMAP]\n"
    "                                [--backgroundModel BACKGROUNDMODEL]\n"
    "                                [--range RANGE RANGE] [--plevel1 PLEVEL1]\n"
    "                                [--plevel2 PLEVEL2] [--keepBait2bait] [--onlySignificant]\n"
    "                                [--outFileName OUTFILENAME] [--outputFormat OUTPUTFORMAT]\n"
    "                                [--dpi DPI] [--linksFile FILE] [--plotData FILE] [--help]\n"
    "                                [--version]\n";

const char* const kHelp =
    "\n"
    "chicChicagoPlotViewpoint plots CHiCAGO viewpoints from chicChicagoScores' own output: "
    "either one bait (--baitID, R Chicago::plotBaits()'s own scatter of N against\n"
    "signed distance) or every bait in a region (--region, absolute genomic coordinates), "
    "either as that scatter or as an arc between each interacting pair (--style).\n"
    "\n"
    "This tool has no Python HiCExplorer counterpart; it is C++ only (cpp/PLAN.md 9.15).\n"
    "\n"
    "Required arguments:\n"
    "  --scores SCORES       chicChicagoScores' own output file.\n"
    "  --baitID BAITID       The baitID (chicChicagoScores' own baitID column) to\n"
    "                        plot. Mutually exclusive with --region; exactly one of\n"
    "                        the two is required.\n"
    "  --region CHROM START END\n"
    "                        Plot every bait --baitmap places in this region\n"
    "                        instead of one (--baitmap is then required too).\n"
    "                        Mutually exclusive with --baitID; exactly one of the\n"
    "                        two is required. --style scatter is not valid with\n"
    "                        --region (a scatter plot has one fixed bait origin).\n"
    "\n"
    "Optional arguments:\n"
    "  --style {scatter,arcs}\n"
    "                        Plot style (Default: scatter). arcs draws one arc per\n"
    "                        interaction between the bait and the other end,\n"
    "                        height log1p(N); scatter is R plotBaits()'s own plot\n"
    "                        and needs --baitID.\n"
    "  --baitmap BAITMAP     CHiCAGO .baitmap file. Required for --region. For\n"
    "                        --baitID, without it the plot is titled by the bare\n"
    "                        baitID and bait2bait rows are not removed (there is no\n"
    "                        way to tell which other ends are baits).\n"
    "  --backgroundModel BACKGROUNDMODEL\n"
    "                        chicChicagoBackgroundModel's own output file (the one\n"
    "                        --scores was itself produced from). Without it the\n"
    "                        Brownian mean/upper-band overlay is not drawn (--style\n"
    "                        scatter only; Bmean is already a column of --scores,\n"
    "                        but the upper band needs the fitted dispersion, only\n"
    "                        recorded in this file).\n"
    "  --range RANGE RANGE   Upstream and downstream distance from each bait to\n"
    "                        plot (Default: 1000000 1000000, R's own default maxD).\n"
    "  --plevel1 PLEVEL1     Score threshold for the more significant colour tier\n"
    "                        (Default: 5, R's own default).\n"
    "  --plevel2 PLEVEL2     Score threshold for the less significant colour tier\n"
    "                        (Default: 3, R's own default).\n"
    "  --keepBait2bait       Keep bait2bait rows (Default: removed, R's own\n"
    "                        removeBait2bait = TRUE default; a bait2bait pair with\n"
    "                        both baits inside --region is always deduplicated to\n"
    "                        one arc regardless of this flag). Only has an effect\n"
    "                        with --baitmap.\n"
    "  --onlySignificant     Drop background rows (score < --plevel2): only the two\n"
    "                        significant colour tiers are drawn, no background.\n"
    "  --rmap RMAP           CHiCAGO .rmap file, for --linksFile's other-end\n"
    "                        fragment spans. Without it a 1 bp placeholder at the\n"
    "                        other end's computed position is written instead.\n"
    "  --linksFile FILE      Write the kept interactions as a pyGenomeTracks\n"
    "                        `links` file (chr1 start1 end1 chr2 start2 end2\n"
    "                        score), to combine with bigwig/gene/BED tracks\n"
    "                        with pyGenomeTracks. Requires --baitmap. Written\n"
    "                        alongside the normal plot output, independent of\n"
    "                        --plotData.\n"
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

// Every row whose baitID is in `bait_ids` (one pass over the file).
std::vector<ScoreRow> read_scores_for_baits(const std::string& path,
                                            const std::unordered_set<long>& bait_ids) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open scores file: " + path);
    std::vector<ScoreRow> rows;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_tab(line);
        if (f.size() < 8) continue;
        const long bait_id = std::stol(f[0]);
        if (bait_ids.find(bait_id) == bait_ids.end()) continue;
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

std::pair<long, long> canonical_pair(long a, long b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

// One pyGenomeTracks `links` row (chr1 start1 end1 chr2 start2 end2 score),
// coordinates already 0-based half-open.
struct LinkRow {
    std::string chrom;
    long start1 = 0, end1 = 0, start2 = 0, end2 = 0;
    double score = 0.0;
};

void write_links_file(const std::string& path, const std::vector<LinkRow>& links) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot open " + path + " for writing");
    out.precision(15);
    for (const auto& l : links) {
        out << l.chrom << '\t' << l.start1 << '\t' << l.end1 << '\t' << l.chrom << '\t' << l.start2
            << '\t' << l.end2 << '\t' << l.score << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    cli::Parser parser(
        "chicChicagoPlotViewpoint",
        "Plots CHiCAGO viewpoints from chicChicagoScores' own output: one bait or a whole "
        "region, as a scatter or as arcs.");
    parser.set_usage(kUsage).set_help(kHelp).set_version_string(hicx::kVersion);

    cli::ArgumentGroup& required = parser.group("Required arguments");
    required.add({"--scores"}).required().input({"txt"}).help("chicChicagoScores' own output file.");
    cli::MutuallyExclusiveGroup& scope = parser.mutually_exclusive(required, /*required=*/true);
    scope.add({"--baitID"}).type("int").help("The baitID to plot.");
    scope.add({"--region"}).nargs(3).metavar("CHROM START END").help("Plot every bait in this region.");

    cli::ArgumentGroup& optional = parser.group("Optional arguments");
    optional.add({"--style"})
        .choices({"scatter", "arcs"})
        .default_value("scatter")
        .help("Plot style.");
    optional.add({"--baitmap"}).input({"baitmap", "txt"}).help("CHiCAGO .baitmap file.");
    optional.add({"--backgroundModel"})
        .input({"txt"})
        .help("chicChicagoBackgroundModel's own output file.");
    optional.add({"--range"})
        .type("int")
        .nargs(2)
        .default_value(hicx::json::Value::array(
            {hicx::json::Value::integer(1000000), hicx::json::Value::integer(1000000)}))
        .help("Upstream and downstream distance from each bait to plot.");
    optional.add({"--plevel1"}).type("float").default_value(5.0).help("More significant score threshold.");
    optional.add({"--plevel2"}).type("float").default_value(3.0).help("Less significant score threshold.");
    optional.add({"--keepBait2bait"})
        .action(cli::Action::StoreTrue)
        .help("Keep bait2bait rows (Default: removed). Only has an effect with --baitmap.");
    optional.add({"--onlySignificant"})
        .action(cli::Action::StoreTrue)
        .help("Drop background rows (score < --plevel2): only the two significant colour "
              "tiers are drawn, no background.");
    optional.add({"--rmap"}).input({"rmap", "txt"}).help("CHiCAGO .rmap file, for --linksFile.");
    optional.add({"--linksFile"})
        .metavar("FILE")
        .output({"links", "txt"})
        .help("Write the kept interactions as a pyGenomeTracks `links` file. Requires "
              "--baitmap.");
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
        const bool region_mode = args.given("region");
        const std::string style = args.str("style");
        if (region_mode && style == "scatter") {
            std::fprintf(stderr, "chicChicagoPlotViewpoint: --style scatter is not valid with "
                                 "--region: a scatter plot has one fixed bait origin. Use "
                                 "--style arcs.\n");
            return 2;
        }
        if (region_mode && !args.given("baitmap")) {
            std::fprintf(stderr, "chicChicagoPlotViewpoint: --region requires --baitmap.\n");
            return 2;
        }
        if (args.given("linksFile") && !args.given("baitmap")) {
            std::fprintf(stderr, "chicChicagoPlotViewpoint: --linksFile requires --baitmap.\n");
            return 2;
        }

        const auto range = args.integers("range");
        const long upstream = range[0];
        const long downstream = range[1];
        const double plevel1 = args.real("plevel1");
        const double plevel2 = args.real("plevel2");

        std::unordered_set<long> all_bait_ids;  // every id in --baitmap, for bait2bait detection
        std::unordered_map<long, hicx::chicago::BaitmapFragment> baitmap_by_id;
        if (args.given("baitmap")) {
            for (const auto& b : hicx::chicago::read_baitmap(args.str("baitmap"))) {
                all_bait_ids.insert(b.id);
                baitmap_by_id.emplace(b.id, b);
            }
        }
        const bool remove_bait2bait = !all_bait_ids.empty() && !args.flag("keepBait2bait");
        const bool only_significant = args.flag("onlySignificant");
        const bool write_links = args.given("linksFile");

        std::unordered_map<long, hicx::chicago::RmapFragment> rmap_by_id;
        if (args.given("rmap")) {
            for (const auto& f : hicx::chicago::read_rmap(args.str("rmap"))) rmap_by_id.emplace(f.id, f);
        }
        // The real fragment span when --rmap has it (the project's usual
        // 0-based half-open convention: a .rmap fragment's 1-based, inclusive
        // [start, end] becomes [start - 1, end)); otherwise a 1 bp
        // placeholder at the computed position, still a valid pyGenomeTracks
        // links row for the arcs style (only the two extremities matter).
        auto other_end_span = [&](long other_end_id, double computed_pos) {
            const auto it = rmap_by_id.find(other_end_id);
            if (it != rmap_by_id.end()) {
                return std::make_pair(it->second.start - 1, it->second.end);
            }
            const long p = static_cast<long>(std::llround(computed_pos));
            return std::make_pair(p, p + 1);
        };

        std::optional<double> dispersion;
        if (args.given("backgroundModel")) dispersion = read_dispersion(args.str("backgroundModel"));
        const bool has_background = dispersion.has_value() && style == "scatter";

        plot::JsonObject data;
        data.add("plotFile", plot::json_string(args.str("outFileName")));
        data.add("style", plot::json_string(style));
        data.add("plevel1", plot::json_number(plevel1));
        data.add("plevel2", plot::json_number(plevel2));
        data.add("dpi", plot::json_int(args.integer("dpi")));

        if (!region_mode) {
            // One bait, its own viewpoint.
            const long bait_id = args.integer("baitID");
            std::string bait_label = "baitID " + std::to_string(bait_id);
            const auto label_it = baitmap_by_id.find(bait_id);
            if (label_it != baitmap_by_id.end() && !label_it->second.name.empty()) {
                bait_label = label_it->second.name;
            }

            auto rows = read_scores_for_baits(args.str("scores"), {bait_id});
            std::vector<double> x, y, height, score;
            std::vector<std::pair<double, double>> bmean_points;
            std::vector<LinkRow> links;
            for (const auto& r : rows) {
                if (!r.has_dist_sign) continue;  // trans: no x position
                if (remove_bait2bait && all_bait_ids.count(r.other_end_id) > 0) continue;
                if (only_significant && r.score < plevel2) continue;
                const double d = static_cast<double>(r.dist_sign);
                if (d < static_cast<double>(-upstream) || d > static_cast<double>(downstream)) continue;
                score.push_back(r.score);
                if (style == "arcs") {
                    x.push_back(0.0);
                    y.push_back(d);
                    height.push_back(std::log1p(r.n));
                } else {
                    x.push_back(d);
                    y.push_back(r.n);
                }
                if (has_background) bmean_points.emplace_back(d, r.bmean);
                if (write_links && label_it != baitmap_by_id.end()) {
                    const auto& bait = label_it->second;
                    const auto [s2, e2] = other_end_span(r.other_end_id, static_cast<double>(bait.start + bait.end) / 2.0 + d);
                    links.push_back({bait.chrom, bait.start - 1, bait.end, s2, e2, r.score});
                }
            }
            if (write_links) write_links_file(args.str("linksFile"), links);

            data.add("title", plot::json_string(bait_label));
            data.add("xLabel", plot::json_string("distance from bait (bp)"));
            if (style == "arcs") {
                data.add("arcX1", plot::json_numbers(x));
                data.add("arcX2", plot::json_numbers(y));
                data.add("arcHeight", plot::json_numbers(height));
                data.add("score", plot::json_numbers(score));
                data.add("anchors", plot::json_numbers({0.0}));
            } else {
                data.add("x", plot::json_numbers(x));
                data.add("y", plot::json_numbers(y));
                data.add("score", plot::json_numbers(score));
            }
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
        } else {
            // Every bait --baitmap places in --region, arcs only, absolute
            // genomic coordinates: a bait's own midpoint, and the other
            // end's midpoint + distSign (already the real signed genomic
            // offset, no separate fragment lookup needed).
            const auto region_args = args.strs("region");
            const std::string chrom = region_args[0];
            const long region_start = std::stol(region_args[1]);
            const long region_end = std::stol(region_args[2]);

            std::unordered_set<long> selected_bait_ids;
            std::unordered_map<long, double> bait_mid;
            for (const auto& [id, b] : baitmap_by_id) {
                if (b.chrom != chrom) continue;
                if (b.end < region_start || b.start > region_end) continue;
                selected_bait_ids.insert(id);
                bait_mid[id] = static_cast<double>(b.start + b.end) / 2.0;
            }

            auto rows = read_scores_for_baits(args.str("scores"), selected_bait_ids);
            std::vector<double> x1, x2, height, score;
            std::vector<LinkRow> links;
            std::set<std::pair<long, long>> seen_pairs;  // canonical (min id, max id)
            for (const auto& r : rows) {
                if (!r.has_dist_sign) continue;
                if (only_significant && r.score < plevel2) continue;
                const double d = static_cast<double>(r.dist_sign);
                if (d < static_cast<double>(-upstream) || d > static_cast<double>(downstream)) continue;
                const bool other_is_selected_bait = selected_bait_ids.count(r.other_end_id) > 0;
                if (other_is_selected_bait) {
                    const auto key = canonical_pair(r.bait_id, r.other_end_id);
                    if (!seen_pairs.insert(key).second) continue;  // already drawn from the other side
                } else if (remove_bait2bait && all_bait_ids.count(r.other_end_id) > 0) {
                    continue;  // a bait outside the region: still a bait2bait row
                }
                const double bait_pos = bait_mid.at(r.bait_id);
                x1.push_back(bait_pos);
                x2.push_back(bait_pos + d);
                height.push_back(std::log1p(r.n));
                score.push_back(r.score);
                if (write_links) {
                    const auto& bait = baitmap_by_id.at(r.bait_id);
                    const auto [s2, e2] = other_end_span(r.other_end_id, bait_pos + d);
                    links.push_back({bait.chrom, bait.start - 1, bait.end, s2, e2, r.score});
                }
            }
            if (write_links) write_links_file(args.str("linksFile"), links);

            std::vector<double> anchors;
            anchors.reserve(bait_mid.size());
            for (const auto& [id, pos] : bait_mid) {
                (void)id;
                anchors.push_back(pos);
            }
            std::sort(anchors.begin(), anchors.end());

            data.add("title", plot::json_string(chrom + ":" + std::to_string(region_start) + "-" +
                                                std::to_string(region_end)));
            data.add("xLabel", plot::json_string("position on " + chrom + " (bp)"));
            data.add("arcX1", plot::json_numbers(x1));
            data.add("arcX2", plot::json_numbers(x2));
            data.add("arcHeight", plot::json_numbers(height));
            data.add("score", plot::json_numbers(score));
            data.add("anchors", plot::json_numbers(anchors));
            data.add("hasBackground", plot::json_bool(false));
        }

        return plot::draw("chicChicagoPlotViewpoint", data.str(), args.opt_str("plotData"));
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoPlotViewpoint: %s\n", error.what());
        return 1;
    }
}
