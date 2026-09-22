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

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>

#include "hicx/argparse.hpp"
#include "hicx/chicago.hpp"
#include "hicx/version.hpp"

namespace {

const char* const kUsage =
    "usage: chicChicagoSignificantInteractions --scores SCORES --rmap RMAP\n"
    "                                          --baitmap BAITMAP\n"
    "                                          [--outFileName OUTFILENAME]\n"
    "                                          [--scoreThreshold SCORETHRESHOLD]\n"
    "                                          [--threads THREADS] [--help] [--version]\n";

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
    "  --threads THREADS     Number of threads (Default: 1). Output row order\n"
    "                        and values do not depend on --threads.\n"
    "  --help, -h            show this help message and exit\n"
    "  --version             show program's version number and exit\n";

struct ScoreRow {
    long bait_id = 0, other_end_id = 0;
    double N = 0, Bmean = 0, Tmean = 0, log_p = 0, score = 0;
    bool has_dist_sign = false;
    long dist_sign = 0;
};

// Parses one already-known-non-empty scores data line (not the header).
// Returns false (caller skips the line) for a malformed row with fewer than
// 8 tab-separated fields, matching the old std::stringstream-based parser's
// `if (f.size() < 8) continue;` behaviour.
bool parse_score_line(std::string_view line, ScoreRow& r) {
    std::string_view fields[8];
    std::size_t n_fields = 0;
    std::string_view rest = line;
    while (n_fields < 8) {
        const auto tab = rest.find('\t');
        if (tab == std::string_view::npos) {
            fields[n_fields++] = rest;
            rest = {};
            break;
        }
        fields[n_fields++] = rest.substr(0, tab);
        rest.remove_prefix(tab + 1);
    }
    if (n_fields < 8) return false;
    auto parse_long = [](std::string_view s) {
        long v = 0;
        std::from_chars(s.data(), s.data() + s.size(), v);
        return v;
    };
    auto parse_double = [](std::string_view s) {
        double v = 0.0;
        std::from_chars(s.data(), s.data() + s.size(), v);
        return v;
    };
    r.bait_id = parse_long(fields[0]);
    r.other_end_id = parse_long(fields[1]);
    r.N = parse_double(fields[2]);
    if (fields[3] == "NA") {
        r.has_dist_sign = false;
    } else {
        r.has_dist_sign = true;
        r.dist_sign = parse_long(fields[3]);
    }
    r.Bmean = parse_double(fields[4]);
    r.Tmean = parse_double(fields[5]);
    r.log_p = parse_double(fields[6]);
    r.score = parse_double(fields[7]);
    return true;
}

// Reads the whole file once, then parses its lines across `threads` workers,
// one contiguous byte range each (moved to the next line start so no thread
// splits a line); only the first range's leading line is the header, skipped
// unconditionally there. Chunks are concatenated back in file order, so the
// result does not depend on --threads.
std::vector<ScoreRow> read_scores(const std::string& path, int threads) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open scores file: " + path);
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();

    const std::size_t total = content.size();
    const int worker_count = std::max(1, threads);
    std::vector<std::size_t> bounds(static_cast<std::size_t>(worker_count) + 1);
    bounds[0] = 0;
    bounds[static_cast<std::size_t>(worker_count)] = total;
    for (int w = 1; w < worker_count; ++w) {
        std::size_t pos = total * static_cast<std::size_t>(w) / static_cast<std::size_t>(worker_count);
        while (pos < total && content[pos] != '\n') ++pos;
        if (pos < total) ++pos;
        bounds[static_cast<std::size_t>(w)] = pos;
    }

    auto parse_range = [&](int idx) {
        std::vector<ScoreRow> local;
        std::size_t pos = bounds[static_cast<std::size_t>(idx)];
        const std::size_t end = bounds[static_cast<std::size_t>(idx) + 1];
        bool header_skipped = (idx != 0);
        while (pos < end) {
            const std::size_t eol = content.find('\n', pos);
            const std::size_t line_end = (eol == std::string::npos || eol > end) ? end : eol;
            std::string_view line(content.data() + pos, line_end - pos);
            pos = (eol == std::string::npos) ? end : eol + 1;
            if (!header_skipped) {
                header_skipped = true;
                continue;  // the header line itself, unconditionally skipped
            }
            if (line.empty()) continue;
            ScoreRow r;
            if (parse_score_line(line, r)) local.push_back(r);
        }
        return local;
    };

    if (worker_count == 1 || total < (1u << 20)) {
        return parse_range(0);
    }
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));
    std::vector<std::vector<ScoreRow>> parts(static_cast<std::size_t>(worker_count));
    for (int w = 0; w < worker_count; ++w) {
        if (bounds[static_cast<std::size_t>(w)] >= bounds[static_cast<std::size_t>(w) + 1]) continue;
        workers.emplace_back([&, w] { parts[static_cast<std::size_t>(w)] = parse_range(w); });
    }
    for (std::thread& worker : workers) worker.join();

    // Same reasoning as chicago.cpp's read_chinput: free the raw file text
    // and each chunk's own buffer as soon as it is no longer needed, instead
    // of holding the file text, the unmerged chunks and the assembled result
    // in memory all at once.
    content.clear();
    content.shrink_to_fit();

    std::size_t total_rows = 0;
    for (auto& p : parts) total_rows += p.size();
    std::vector<ScoreRow> out;
    out.reserve(total_rows);
    for (auto& p : parts) {
        out.insert(out.end(), std::make_move_iterator(p.begin()), std::make_move_iterator(p.end()));
        std::vector<ScoreRow>().swap(p);
    }
    return out;
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
    optional.add({"--threads"}).type("int").default_value(1).help("Number of threads.");
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

        const int threads = static_cast<int>(args.integer("threads"));
        auto rows = read_scores(args.str("scores"), threads);
        const double threshold = args.real("scoreThreshold");

        std::ofstream out(args.str("outFileName"), std::ios::binary);
        if (!out) {
            throw std::runtime_error("[Errno 2] No such file or directory: '" +
                                     args.str("outFileName") + "'");
        }
        out.precision(15);
        out << "bait_chr\tbait_start\tbait_end\tbait_name\totherEnd_chr\totherEnd_start\t"
               "otherEnd_end\tN\tdistSign\tscore\n";

        // Every row is filtered/formatted independently of every other (only
        // read-only lookups into rmap_by_id/baitmap_by_id), so rows are split
        // into contiguous chunks, one per thread, each formatting into its own
        // string; chunks are written out in original row order, so the output
        // file is identical no matter how many threads ran it.
        const std::size_t n = rows.size();
        const int worker_count = std::max(1, threads);
        auto format_range = [&](std::size_t first, std::size_t last) {
            std::string buf;
            for (std::size_t i = first; i < last; ++i) {
                const auto& r = rows[i];
                if (!(r.score >= threshold)) continue;
                auto bait_it = baitmap_by_id.find(r.bait_id);
                auto oe_it = rmap_by_id.find(r.other_end_id);
                if (bait_it == baitmap_by_id.end() || oe_it == rmap_by_id.end()) continue;
                std::ostringstream row;
                row.precision(15);
                row << bait_it->second.chrom << "\t" << bait_it->second.start << "\t"
                    << bait_it->second.end << "\t" << bait_it->second.name << "\t"
                    << oe_it->second.chrom << "\t" << oe_it->second.start << "\t" << oe_it->second.end
                    << "\t" << r.N << "\t";
                if (r.has_dist_sign) row << r.dist_sign; else row << "NA";
                row << "\t" << r.score << "\n";
                buf += row.str();
            }
            return buf;
        };

        if (worker_count == 1 || n < 4096) {
            out << format_range(0, n);
        } else {
            std::vector<std::string> chunks(static_cast<std::size_t>(worker_count));
            std::vector<std::thread> workers;
            workers.reserve(static_cast<std::size_t>(worker_count));
            const std::size_t chunk_size = (n + static_cast<std::size_t>(worker_count) - 1) /
                                            static_cast<std::size_t>(worker_count);
            for (int w_idx = 0; w_idx < worker_count; ++w_idx) {
                const std::size_t first = static_cast<std::size_t>(w_idx) * chunk_size;
                const std::size_t last = std::min(n, first + chunk_size);
                if (first >= last) break;
                workers.emplace_back(
                    [&, w_idx, first, last] { chunks[static_cast<std::size_t>(w_idx)] = format_range(first, last); });
            }
            for (std::thread& worker : workers) worker.join();
            for (const std::string& chunk : chunks) out << chunk;
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chicChicagoSignificantInteractions: %s\n", error.what());
        return 1;
    }
    return 0;
}
