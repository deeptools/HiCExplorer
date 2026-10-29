#include "hicx/matrix_qc_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <set>
#include <sstream>

#include "hicx/build_matrix.hpp"

namespace hicx::matrixqc {

namespace {

// --- number formats ------------------------------------------------------------------

std::string esc(const std::string& text) {
    std::string out;
    for (const char c : text) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c;
        }
    }
    return out;
}

std::string printf_string(const char* format, double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), format, value);
    return buffer;
}

// The tables: integral values without decimals, others with 6 significant digits.
std::string tsv_number(double value) {
    if (std::isnan(value)) {
        return "NA";
    }
    if (std::isinf(value)) {
        return value > 0 ? "Inf" : "-Inf";
    }
    if (std::fabs(value) < 1e15 && value == std::floor(value)) {
        return printf_string("%.0f", value);
    }
    return printf_string("%.6g", value);
}

std::string thousands(double value) {
    if (std::isnan(value)) {
        return "n/a";
    }
    std::string digits = printf_string("%.0f", std::fabs(value));
    std::string out;
    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i != 0 && (digits.size() - i) % 3 == 0) {
            out += ',';
        }
        out += digits[i];
    }
    return (value < 0 ? "-" : "") + out;
}

// A count for a narrow table cell: 12,345 or 1.23M.
std::string short_count(double value) {
    if (std::isnan(value)) {
        return "n/a";
    }
    const double magnitude = std::fabs(value);
    if (magnitude < 1e5) {
        return thousands(value);
    }
    if (magnitude < 1e6) {
        return printf_string("%.0f", value / 1e3) + "k";
    }
    if (magnitude < 1e9) {
        return printf_string("%.2f", value / 1e6) + "M";
    }
    return printf_string("%.2f", value / 1e9) + "G";
}

std::string percent(double fraction, int decimals = 2) {
    if (std::isnan(fraction)) {
        return "n/a";
    }
    const double value = 100.0 * fraction;
    if (value > 0 && value < 0.5 * std::pow(10.0, -decimals)) {
        return "<" + printf_string(decimals == 0 ? "%.0f" : (decimals == 1 ? "%.1f" : "%.2f"),
                                   std::pow(10.0, -decimals)) + "%";
    }
    return printf_string(decimals == 0 ? "%.0f" : (decimals == 1 ? "%.1f" : "%.2f"), value) + "%";
}

std::string bp_text(std::int64_t bp) {
    if (bp < 1000) {
        return std::to_string(bp) + " bp";
    }
    if (bp < 1000000) {
        return printf_string("%.3g", static_cast<double>(bp) / 1e3) + " kb";
    }
    if (bp < 1000000000) {
        return printf_string("%.3g", static_cast<double>(bp) / 1e6) + " Mb";
    }
    return printf_string("%.3g", static_cast<double>(bp) / 1e9) + " Gb";
}

std::string fixed(double value, int decimals) {
    if (std::isnan(value)) {
        return "n/a";
    }
    char format[8];
    std::snprintf(format, sizeof(format), "%%.%df", decimals);
    return printf_string(format, value);
}

std::string iso_date() {
    const std::time_t now = std::time(nullptr);
    char buffer[16];
    std::tm tm{};
    localtime_r(&now, &tm);
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d", &tm);
    return buffer;
}

// --- samples ---------------------------------------------------------------------------------

struct Sample {
    const MatrixStats* matrix = nullptr;
    const FileStats* file = nullptr;
    std::string tag;  // S1, S2, ...
};

std::vector<Sample> samples_of(const std::vector<FileStats>& files) {
    std::vector<Sample> out;
    for (const FileStats& file : files) {
        for (const MatrixStats& matrix : file.matrices) {
            out.push_back({&matrix, &file, "S" + std::to_string(out.size() + 1)});
        }
    }
    return out;
}

// Groups of at most three sample columns, so that a table stays readable at
// 400 px.
std::vector<std::vector<std::size_t>> column_blocks(std::size_t n) {
    std::vector<std::vector<std::size_t>> blocks;
    for (std::size_t i = 0; i < n; i += 3) {
        std::vector<std::size_t> block;
        for (std::size_t j = i; j < std::min(n, i + 3); ++j) {
            block.push_back(j);
        }
        blocks.push_back(std::move(block));
    }
    return blocks;
}

const MapResolution* map_of(const Sample& sample) {
    return sample.file->map_resolution.method.empty() ? nullptr : &sample.file->map_resolution;
}

}  // namespace

// --- the statistics table -----------------------------------------------------------------

std::vector<std::pair<std::string, std::string>> flatten_stats(const MatrixStats& m,
                                                               const MapResolution* map) {
    std::vector<std::pair<std::string, std::string>> rows;
    const auto add = [&rows](const std::string& key, const std::string& value) {
        rows.emplace_back(key, value);
    };
    const auto number = [&add](const std::string& key, double value) { add(key, tsv_number(value)); };
    const auto ratio = [](double part, double whole) { return whole > 0 ? part / whole : kNaN; };

    add("sample.uri", m.uri);
    number("sample.resolution_bp", static_cast<double>(m.resolution));
    number("sample.bins", static_cast<double>(m.bins));
    number("sample.pixels", static_cast<double>(m.pixels));
    number("sample.chromosomes", static_cast<double>(m.chromosomes_total));

    if (m.read_qc.has_value()) {
        const ReadQc& q = *m.read_qc;
        add("qc.tool", q.tool);
        add("qc.tool_version", q.tool_version);
        if (q.min_rest_site_distance.has_value()) {
            number("qc.min_rest_site_distance", static_cast<double>(*q.min_rest_site_distance));
        }
        number("qc.max_library_insert_size", static_cast<double>(q.max_library_insert_size));
        const auto count = [&](const std::string& key, std::int64_t value) {
            number("qc." + key, static_cast<double>(value));
        };
        count("sequenced_reads", q.sequenced_reads);
        count("pairs_mappable_unique_high_quality", q.mappable);
        count("hic_contacts", q.hic_contacts);
        count("one_mate_unmapped", q.one_mate_unmapped);
        count("one_mate_not_unique", q.one_mate_not_unique);
        count("low_mapping_quality", q.low_mapping_quality);
        for (std::size_t i = 0; i < q.dangling_end.size(); ++i) {
            const std::string restriction =
                i < q.restriction_sequences.size() ? q.restriction_sequences[i] : std::to_string(i);
            count("dangling_end." + restriction, q.dangling_end[i]);
        }
        count("self_ligation", q.self_ligation);
        count("one_mate_not_close_to_rest_site", q.not_close_to_rest_site);
        count("same_fragment", q.same_fragment);
        count("self_circle", q.self_circle);
        count("duplicated_pairs", q.duplicated_pairs);
        count("inter_chromosomal", q.inter_chromosomal);
        count("intra_short_range", q.intra_short);
        count("intra_long_range", q.intra_long);
        count("inward_pairs", q.orientation_pairs[0]);
        count("outward_pairs", q.orientation_pairs[1]);
        count("left_pairs", q.orientation_pairs[2]);
        count("right_pairs", q.orientation_pairs[3]);
        const char* names[4] = {"inward", "outward", "left", "right"};
        for (std::size_t i = 0; i < q.distance_contacts.size() && i < q.lower_edges.size(); ++i) {
            const std::string edge = std::to_string(q.lower_edges[i]);
            if (i + 1 < q.lower_edges.size() && q.lower_edges[i + 1] == q.lower_edges[i]) {
                continue;  // an empty bin of the integer grid
            }
            count("distance.contacts." + edge, q.distance_contacts[i]);
            for (std::size_t o = 0; o < 4; ++o) {
                if (i < q.orientation_distance[o].size()) {
                    count(std::string("distance.") + names[o] + "." + edge,
                          q.orientation_distance[o][i]);
                }
            }
        }
    }
    if (!m.pixel_statistics) {
        add("note", m.note);
        return rows;
    }

    number("contacts.total", m.total_contacts);
    number("contacts.cis", m.cis);
    number("contacts.trans", m.trans);
    number("contacts.cis_fraction", ratio(m.cis, m.total_contacts));
    number("contacts.trans_fraction", ratio(m.trans, m.total_contacts));
    number("contacts.short_range_cutoff_bp", static_cast<double>(m.short_range_cutoff));
    number("contacts.cis_short_range", m.intra_short);
    number("contacts.cis_long_range", m.intra_long);
    number("contacts.cis_short_range_fraction", ratio(m.intra_short, m.cis));
    number("contacts.cis_long_range_fraction", ratio(m.intra_long, m.cis));
    number("diagonal.main_fraction", m.checks.diagonal_fraction);
    number("diagonal.first_offdiagonal_fraction", m.checks.first_offdiag_fraction);

    number("decay.exponent.genome", m.decay.exponent);
    number("decay.fit_points.genome", static_cast<double>(m.decay.fit_points));
    for (const ChromosomeSummary& c : m.chromosomes) {
        number("decay.exponent." + c.name, c.decay.exponent);
    }
    for (std::size_t slot = 0; slot < m.decay.pairs.size(); ++slot) {
        if (m.decay.pairs[slot] <= 0) {
            continue;
        }
        const std::string edge = std::to_string(decay_slot_lower_edge(slot));
        number("decay.genome.contacts." + edge, m.decay.contacts[slot]);
        number("decay.genome.pairs." + edge, m.decay.pairs[slot]);
        number("decay.genome.contacts_per_pair." + edge, m.decay.contacts[slot] / m.decay.pairs[slot]);
    }

    number("coverage.empty_bins_fraction", m.empty_fraction);
    number("coverage.min", m.coverage.min);
    number("coverage.p01", m.coverage.p01);
    number("coverage.p05", m.coverage.p05);
    number("coverage.p25", m.coverage.p25);
    number("coverage.median", m.coverage.median);
    number("coverage.p75", m.coverage.p75);
    number("coverage.p95", m.coverage.p95);
    number("coverage.p99", m.coverage.p99);
    number("coverage.max", m.coverage.max);
    number("coverage.mad.median", m.mad_median);
    number("coverage.mad.mad", m.mad_value);
    number("coverage.mad.removed_bins", static_cast<double>(m.mad_removed));
    number("coverage.mad.removed_bins_empty", static_cast<double>(m.mad_removed_empty));
    number("coverage.mad.removed_fraction", m.mad_removed_fraction);
    for (const ChromosomeSummary& c : m.chromosomes) {
        number("coverage.chromosome.empty_fraction." + c.name, c.empty_fraction);
        number("coverage.chromosome.median." + c.name, c.median_coverage);
        number("coverage.chromosome.mad_removed_fraction." + c.name, c.mad_removed_fraction);
    }

    number("resolution.fraction_bins_enough_contacts", m.rao_fraction_all);
    number("resolution.fraction_nonempty_bins_enough_contacts", m.rao_fraction_nonempty);
    if (map != nullptr) {
        add("map_resolution.method", map->method);
        number("map_resolution.bp", map->bin_size > 0 ? static_cast<double>(map->bin_size) : kNaN);
        number("map_resolution.bp_nonempty_bins",
               map->bin_size_nonempty > 0 ? static_cast<double>(map->bin_size_nonempty) : kNaN);
        add("map_resolution.finest_bin_size_reached", map->finest_reached ? "yes" : "no");
        for (const LadderStep& step : map->steps) {
            number("map_resolution.fraction_all." + std::to_string(step.bin_size), step.fraction_all);
        }
    }

    const BalanceInfo& b = m.balance;
    add("balance.weight_column", b.has_weight_column ? "yes" : "no");
    add("balance.weights_usable", b.usable ? "yes" : "no");
    number("balance.nan_weight_fraction", b.nan_fraction);
    for (const auto& [name, value] : b.attributes) {
        add("balance.attribute." + name, value);
    }
    number("balance.marginal_mean", b.marginal_mean);
    number("balance.marginal_cv", b.marginal_cv);

    const TransTable& t = m.trans_table;
    number("trans.contacts_between_listed_chromosomes", t.total);
    for (std::size_t a = 0; a < t.names.size(); ++a) {
        for (std::size_t c = a + 1; c < t.names.size(); ++c) {
            const std::string pair = t.names[a] + "." + t.names[c];
            const double observed = t.observed[a * t.names.size() + c];
            const double expected = t.expected[a * t.names.size() + c];
            number("trans.observed." + pair, observed);
            number("trans.expected." + pair, expected);
            number("trans.observed_over_expected." + pair, expected > 0 ? observed / expected : kNaN);
        }
    }
    number("trans.flagged_pairs", static_cast<double>(t.flagged.size()));
    std::string flagged;
    for (const TransTable::Pair& p : t.flagged) {
        flagged += (flagged.empty() ? "" : ";") + t.names[p.a] + "-" + t.names[p.b];
    }
    add("trans.flagged_pair_list", flagged);

    const CheapChecks& k = m.checks;
    number("checks.x_autosome_coverage_ratio", k.x_autosome_ratio);
    add("checks.chrX_present", k.has_x ? "yes" : "no");
    add("checks.chrY_present", k.has_y ? "yes" : "no");
    number("checks.chrY_coverage_fraction", k.y_fraction);
    add("checks.chrM_present", k.has_m ? "yes" : "no");
    number("checks.chrM_coverage_fraction", k.m_fraction);
    return rows;
}

std::string render_table(const std::vector<FileStats>& files) {
    const std::vector<Sample> samples = samples_of(files);
    std::vector<std::vector<std::pair<std::string, std::string>>> columns;
    std::vector<std::string> keys;
    std::set<std::string> seen;
    for (const Sample& sample : samples) {
        columns.push_back(flatten_stats(*sample.matrix, map_of(sample)));
        for (const auto& [key, value] : columns.back()) {
            if (seen.insert(key).second) {
                keys.push_back(key);
            }
        }
    }
    std::string out = "statistic";
    for (const Sample& sample : samples) {
        out += "\t" + sample.matrix->label;
    }
    out += "\n";
    for (const std::string& key : keys) {
        out += key;
        for (const auto& column : columns) {
            out += "\t";
            for (const auto& [name, value] : column) {
                if (name == key) {
                    out += value;
                    break;
                }
            }
        }
        out += "\n";
    }
    return out;
}

// --- the HTML report ---------------------------------------------------------------------------

namespace {

// Colour-blind safe colours (Okabe and Ito).
const char* const kBlue = "#0072B2";
const char* const kOrange = "#E69F00";
const char* const kGreen = "#009E73";
const char* const kRed = "#D55E00";
const char* const kPurple = "#CC79A7";
const char* const kSky = "#56B4E9";
const char* const kGrey = "#999999";

std::string hex_colour(double r, double g, double b) {
    char buffer[8];
    std::snprintf(buffer, sizeof(buffer), "#%02X%02X%02X",
                  static_cast<int>(std::lround(std::clamp(r, 0.0, 255.0))),
                  static_cast<int>(std::lround(std::clamp(g, 0.0, 255.0))),
                  static_cast<int>(std::lround(std::clamp(b, 0.0, 255.0))));
    return buffer;
}

// White (t = 0) to dark blue (t = 1).
std::string sequential_colour(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return hex_colour(255 + (8 - 255) * t, 255 + (69 - 255) * t, 255 + (148 - 255) * t);
}

// Blue (t = -1), white (0), orange (t = 1).
std::string diverging_colour(double t) {
    t = std::clamp(t, -1.0, 1.0);
    if (t >= 0) {
        return hex_colour(255 + (213 - 255) * t, 255 + (94 - 255) * t, 255 + (0 - 255) * t);
    }
    return hex_colour(255 + (0 - 255) * -t, 255 + (114 - 255) * -t, 255 + (178 - 255) * -t);
}

struct Segment {
    double width;  // percent of the track
    std::string colour;
};

// One bar track: segments left to right, the rest empty.
void track_cell(std::ostream& out, const std::vector<Segment>& segments, const char* css = "track") {
    out << "<td class=\"" << css
        << "\"><table class=\"bar\" width=\"100%\" cellspacing=\"0\" cellpadding=\"0\" border=\"0\"><tr>";
    double used = 0;
    for (const Segment& segment : segments) {
        const double width = std::clamp(segment.width, 0.0, 100.0 - used);
        if (width <= 0) {
            continue;
        }
        used += width;
        out << "<td width=\"" << printf_string("%.2f", width) << "%\" bgcolor=\"" << segment.colour
            << "\" height=\"12\">&nbsp;</td>";
    }
    if (used < 99.99) {
        out << "<td height=\"12\">&nbsp;</td>";
    }
    out << "</tr></table></td>";
}

void table_open(std::ostream& out, const std::string& summary, const char* css = "data") {
    out << "<table class=\"" << css << "\" width=\"100%\" border=\"1\" cellspacing=\"0\" cellpadding=\"3\" summary=\""
        << esc(summary) << "\">\n";
}

void section_open(std::ostream& out, const std::string& id, const std::string& title) {
    out << "<h2 id=\"" << id << "\">" << esc(title) << "</h2>\n";
}

struct Row {
    std::string label;
    std::vector<std::string> cells;  // one per sample
    bool heading = false;
};

// A statistics table: statistic in the first column, one column per sample;
// more than three samples give one table per block of three.
void stat_table(std::ostream& out, const std::vector<Sample>& samples, const std::vector<Row>& rows,
                const std::string& summary) {
    for (const std::vector<std::size_t>& block : column_blocks(samples.size())) {
        table_open(out, summary);
        out << "<colgroup><col width=\"" << (block.size() == 1 ? 55 : 100 - 20 * static_cast<int>(block.size()))
            << "%\">";
        for (std::size_t i = 0; i < block.size(); ++i) {
            out << "<col>";
        }
        out << "</colgroup>\n<tr><th align=\"left\">Statistic</th>";
        for (const std::size_t i : block) {
            out << "<th align=\"right\">" << samples[i].tag << "</th>";
        }
        out << "</tr>\n";
        for (const Row& row : rows) {
            if (row.heading) {
                out << "<tr class=\"group\"><td colspan=\"" << block.size() + 1 << "\"><b>" << esc(row.label)
                    << "</b></td></tr>\n";
                continue;
            }
            out << "<tr><td>" << esc(row.label) << "</td>";
            for (const std::size_t i : block) {
                out << "<td align=\"right\">" << esc(i < row.cells.size() ? row.cells[i] : std::string()) << "</td>";
            }
            out << "</tr>\n";
        }
        out << "</table>\n";
    }
}

template <class F>
std::vector<std::string> per_sample(const std::vector<Sample>& samples, F&& f) {
    std::vector<std::string> cells;
    for (const Sample& sample : samples) {
        cells.push_back(f(*sample.matrix, sample));
    }
    return cells;
}

double fraction_of(double part, double whole) { return whole > 0 ? part / whole : kNaN; }

// ---- sections ------------------------------------------------------------------------------------

void write_samples(std::ostream& out, const std::vector<Sample>& samples) {
    section_open(out, "samples", "Samples");
    table_open(out, "The samples of the report");
    out << "<colgroup><col width=\"9%\"><col width=\"45%\"><col width=\"16%\"><col width=\"15%\"><col width=\"15%\"></colgroup>\n"
        << "<tr><th align=\"left\">ID</th><th align=\"left\">Matrix</th><th align=\"right\">Bin size</th>"
           "<th align=\"right\">Bins</th><th align=\"right\">Pixels</th></tr>\n";
    for (const Sample& sample : samples) {
        const MatrixStats& m = *sample.matrix;
        out << "<tr><td>" << sample.tag << "</td><td class=\"name\">" << esc(m.label) << "</td><td align=\"right\">"
            << (m.resolution > 0 ? bp_text(m.resolution) : std::string("variable")) << "</td><td align=\"right\">"
            << short_count(static_cast<double>(m.bins)) << "</td><td align=\"right\">"
            << short_count(static_cast<double>(m.pixels)) << "</td></tr>\n";
    }
    out << "</table>\n";
    for (const Sample& sample : samples) {
        if (!sample.matrix->note.empty()) {
            out << "<p class=\"note\">" << sample.tag << ": " << esc(sample.matrix->note) << "</p>\n";
        }
    }
}

void write_summary(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "summary", "Summary");
    std::vector<Row> rows;
    rows.push_back({"Contacts", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.pixel_statistics ? short_count(m.total_contacts) : std::string("n/a");
                    })});
    rows.push_back({"Intra-chromosomal (cis)", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.cis, m.total_contacts), 1);
                    })});
    rows.push_back({"Inter-chromosomal (trans)", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.trans, m.total_contacts), 1);
                    })});
    rows.push_back({"Cis short range (< " + bp_text(options.short_range_cutoff) + ")",
                    per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.intra_short, m.cis), 1);
                    })});
    rows.push_back({"Main diagonal", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.checks.diagonal_fraction, 1);
                    })});
    rows.push_back({"First off-diagonal", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.checks.first_offdiag_fraction, 1);
                    })});
    rows.push_back({"Decay exponent", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return fixed(m.decay.exponent, 2);
                    })});
    rows.push_back({"Empty bins", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.empty_fraction, 1);
                    })});
    rows.push_back({"Bins removed by the MAD filter", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.mad_removed_fraction, 1);
                    })});
    rows.push_back({"Map resolution", per_sample(samples, [](const MatrixStats& m, const Sample& s) {
                        const MapResolution* map = map_of(s);
                        if (map == nullptr || !m.pixel_statistics) {
                            return std::string("n/a");
                        }
                        return map->bin_size > 0 ? (map->finest_reached ? "<= " : "") + bp_text(map->bin_size)
                                                 : std::string("not reached");
                    })});
    rows.push_back({"Balancing weights", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return std::string(m.balance.usable ? "yes" : "no");
                    })});
    rows.push_back({"Flagged chromosome pairs", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.pixel_statistics ? std::to_string(m.trans_table.flagged.size()) : std::string("n/a");
                    })});
    rows.push_back({"Read pairs from the BAM files", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.read_qc.has_value() ? short_count(static_cast<double>(m.read_qc->sequenced_reads))
                                                     : std::string("n/a");
                    })});
    stat_table(out, samples, rows, "Summary of the statistics");
}

void write_read_qc(std::ostream& out, const std::vector<Sample>& samples) {
    std::set<const FileStats*> done;
    bool any = false;
    for (const Sample& sample : samples) {
        if (!sample.matrix->read_qc.has_value() || !done.insert(sample.file).second) {
            continue;
        }
        if (!any) {
            section_open(out, "library", "Library statistics stored by hicBuildMatrix");
            out << "<p>The counts below were written into the metadata of the matrix when it was built "
                   "(key hicexplorer_qc) and are the numbers of QC.log. The distance of a read pair is the "
                   "distance between the start positions of its two reads, measured before binning.</p>\n";
            any = true;
        }
        const ReadQc& q = *sample.matrix->read_qc;
        out << "<h3>" << sample.tag << ": " << esc(sample.file->path) << "</h3>\n";
        out << "<p class=\"note\">Written by " << esc(q.tool.empty() ? "hicBuildMatrix" : q.tool) << " "
            << esc(q.tool_version) << ". Maximum library insert size " << q.max_library_insert_size << " bp";
        if (q.min_rest_site_distance.has_value()) {
            out << ", minimum restriction site distance " << *q.min_rest_site_distance << " bp";
        }
        if (!q.restriction_sequences.empty()) {
            out << ", restriction sequence";
            for (std::size_t i = 0; i < q.restriction_sequences.size(); ++i) {
                out << (i == 0 ? " " : ", ") << esc(q.restriction_sequences[i]);
            }
        }
        out << ".</p>\n";

        struct Item {
            std::string label;
            std::int64_t count;
        };
        const auto group = [&](const std::string& title, double denominator, const std::vector<Item>& items) {
            table_open(out, title);
            out << "<colgroup><col width=\"46%\"><col width=\"20%\"><col width=\"14%\"><col width=\"20%\"></colgroup>\n"
                << "<tr><th align=\"left\">" << esc(title) << "</th><th align=\"right\">Count</th><th align=\"right\">%</th>"
                   "<th></th></tr>\n";
            for (const Item& item : items) {
                const double share = fraction_of(static_cast<double>(item.count), denominator);
                out << "<tr><td>" << esc(item.label) << "</td><td align=\"right\">"
                    << thousands(static_cast<double>(item.count)) << "</td><td align=\"right\">"
                    << (std::isnan(share) ? std::string("n/a") : fixed(100.0 * share, 2)) << "</td>";
                track_cell(out, {{std::isnan(share) ? 0.0 : 100.0 * share, kBlue}});
                out << "</tr>\n";
            }
            out << "</table>\n";
        };
        group("Share of sequenced reads (" + thousands(static_cast<double>(q.sequenced_reads)) + ")",
              static_cast<double>(q.sequenced_reads),
              {{"Pairs mappable, unique and high quality", q.mappable},
               {"Hi-C contacts", q.hic_contacts},
               {"One mate unmapped", q.one_mate_unmapped},
               {"One mate not unique", q.one_mate_not_unique},
               {"Low mapping quality", q.low_mapping_quality}});
        std::vector<Item> discarded;
        for (std::size_t i = 0; i < q.dangling_end.size(); ++i) {
            discarded.push_back({"Dangling end " + (i < q.dangling_sequences.size() ? q.dangling_sequences[i] : std::string()) +
                                     " (restriction sequence " +
                                     (i < q.restriction_sequences.size() ? q.restriction_sequences[i] : std::string()) + ")",
                                 q.dangling_end[i]});
        }
        if (q.restriction_cut_file) {
            discarded.push_back({q.keep_self_ligation ? "Self ligation (not removed)" : "Self ligation (removed)",
                                 q.self_ligation});
            discarded.push_back({"One mate not close to rest site", q.not_close_to_rest_site});
        }
        discarded.push_back({"Same fragment", q.same_fragment});
        if (q.restriction_cut_file) {
            discarded.push_back({"Self circle", q.self_circle});
        }
        discarded.push_back({"Duplicated pairs", q.duplicated_pairs});
        group("Share of mappable, unique and high quality pairs (" + thousands(static_cast<double>(q.mappable)) + ")",
              static_cast<double>(q.mappable), discarded);
        if (q.hic_contacts > 0) {
            group("Share of the valid pairs used (" + thousands(static_cast<double>(q.hic_contacts)) + ")",
                  static_cast<double>(q.hic_contacts),
                  {{"Inter-chromosomal", q.inter_chromosomal},
                   {"Intra-chromosomal short range (< " + bp_text(q.short_range_cutoff) + ")", q.intra_short},
                   {"Intra-chromosomal long range (>= " + bp_text(q.short_range_cutoff) + ")", q.intra_long},
                   {"Read pair type: inward", q.orientation_pairs[0]},
                   {"Read pair type: outward", q.orientation_pairs[1]},
                   {"Read pair type: left", q.orientation_pairs[2]},
                   {"Read pair type: right", q.orientation_pairs[3]}});
        }

        // The distance plots: one row per distance bin of the integer grid.
        const std::size_t bins = std::min(q.distance_contacts.size(), q.lower_edges.size());
        std::vector<std::size_t> rows;
        for (std::size_t i = 0; i < bins; ++i) {
            if (i + 1 < q.lower_edges.size() && q.lower_edges[i + 1] == q.lower_edges[i]) {
                continue;
            }
            rows.push_back(i);
        }
        const auto bin_label = [&](std::size_t i) {
            return bp_text(q.lower_edges[i]) + " to " +
                   (i + 1 < q.lower_edges.size() ? bp_text(q.lower_edges[i + 1]) : std::string("end"));
        };
        double maximum = 0;
        for (const std::size_t i : rows) {
            maximum = std::max(maximum, static_cast<double>(q.distance_contacts[i]));
        }
        if (maximum > 0) {
            out << "<h4>Intra-chromosomal read pairs by distance</h4>\n";
            table_open(out, "Read pairs per distance bin");
            out << "<colgroup><col width=\"34%\"><col width=\"46%\"><col width=\"20%\"></colgroup>\n"
                << "<tr><th align=\"left\">Distance</th><th></th><th align=\"right\">Pairs</th></tr>\n";
            for (const std::size_t i : rows) {
                const double value = static_cast<double>(q.distance_contacts[i]);
                out << "<tr><td>" << esc(bin_label(i)) << "</td>";
                track_cell(out, {{100.0 * value / maximum, kBlue}});
                out << "<td align=\"right\">" << thousands(value) << "</td></tr>\n";
            }
            out << "</table>\n";

            const char* names[4] = {"Inward", "Outward", "Left", "Right"};
            const char* colours[4] = {kBlue, kOrange, kGreen, kPurple};
            out << "<h4>Read pair type by distance</h4>\n"
                << "<p class=\"note\">The length of a bar is the number of pairs of the distance bin relative to "
                   "the largest bin; its colours are the shares of the four read pair types.</p>\n";
            out << "<table class=\"legend\" cellspacing=\"0\" cellpadding=\"2\" border=\"0\"><tr>";
            for (int o = 0; o < 4; ++o) {
                out << "<td width=\"14\" bgcolor=\"" << colours[o] << "\">&nbsp;</td><td>" << names[o] << "</td>";
            }
            out << "</tr></table>\n";
            table_open(out, "Read pair types per distance bin");
            out << "<colgroup><col width=\"34%\"><col width=\"46%\"><col width=\"20%\"></colgroup>\n"
                << "<tr><th align=\"left\">Distance</th><th></th><th align=\"right\">Pairs</th></tr>\n";
            for (const std::size_t i : rows) {
                double total = 0;
                for (int o = 0; o < 4; ++o) {
                    if (i < q.orientation_distance[static_cast<std::size_t>(o)].size()) {
                        total += static_cast<double>(q.orientation_distance[static_cast<std::size_t>(o)][i]);
                    }
                }
                out << "<tr><td>" << esc(bin_label(i)) << "</td>";
                std::vector<Segment> segments;
                for (int o = 0; o < 4 && total > 0; ++o) {
                    const auto& histogram = q.orientation_distance[static_cast<std::size_t>(o)];
                    const double value = i < histogram.size() ? static_cast<double>(histogram[i]) : 0.0;
                    segments.push_back({100.0 * total / maximum * (value / total), colours[o]});
                }
                track_cell(out, segments);
                out << "<td align=\"right\">" << thousands(total) << "</td></tr>\n";
            }
            out << "</table>\n";
        }
    }
}

void write_contacts(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "contacts", "Contacts, cis and trans");
    out << "<p>Contacts are the counts of the stored pixels. Cis contacts join two bins of one chromosome, trans "
           "contacts two chromosomes. The distance of a cis contact is the distance between the starts of its two "
           "bins, so a distance threshold is exact only up to one bin size. The contacts on the main diagonal "
           "(both reads in one bin) are closer than one bin size and count as short range. Short range is below "
        << bp_text(options.short_range_cutoff) << ".</p>\n";
    std::vector<Row> rows;
    const auto count_row = [&](const std::string& label, double MatrixStats::*member) {
        rows.push_back({label, per_sample(samples, [member](const MatrixStats& m, const Sample&) {
                            return m.pixel_statistics ? short_count(m.*member) : std::string("n/a");
                        })});
    };
    count_row("Contacts", &MatrixStats::total_contacts);
    count_row("Cis", &MatrixStats::cis);
    count_row("Trans", &MatrixStats::trans);
    count_row("Cis short range", &MatrixStats::intra_short);
    count_row("Cis long range", &MatrixStats::intra_long);
    rows.push_back({"Cis, % of contacts", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.cis, m.total_contacts));
                    })});
    rows.push_back({"Trans, % of contacts", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.trans, m.total_contacts));
                    })});
    rows.push_back({"Short range, % of cis", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.intra_short, m.cis));
                    })});
    rows.push_back({"Long range, % of cis", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(fraction_of(m.intra_long, m.cis));
                    })});
    rows.push_back({"Main diagonal, % of contacts", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.checks.diagonal_fraction);
                    })});
    rows.push_back({"First off-diagonal, % of contacts", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.checks.first_offdiag_fraction);
                    })});
    stat_table(out, samples, rows, "Contacts");

    out << "<h3>Composition of the contacts</h3>\n"
        << "<table class=\"legend\" cellspacing=\"0\" cellpadding=\"2\" border=\"0\"><tr>"
        << "<td width=\"14\" bgcolor=\"" << kBlue << "\">&nbsp;</td><td>cis short range</td>"
        << "<td width=\"14\" bgcolor=\"" << kSky << "\">&nbsp;</td><td>cis long range</td>"
        << "<td width=\"14\" bgcolor=\"" << kOrange << "\">&nbsp;</td><td>trans</td></tr></table>\n";
    table_open(out, "Composition of the contacts per sample");
    out << "<colgroup><col width=\"14%\"><col width=\"86%\"></colgroup>\n";
    for (const Sample& sample : samples) {
        const MatrixStats& m = *sample.matrix;
        if (!m.pixel_statistics || m.total_contacts <= 0) {
            continue;
        }
        out << "<tr><td>" << sample.tag << "</td>";
        track_cell(out, {{100.0 * m.intra_short / m.total_contacts, kBlue},
                         {100.0 * m.intra_long / m.total_contacts, kSky},
                         {100.0 * m.trans / m.total_contacts, kOrange}});
        out << "</tr>\n";
    }
    out << "</table>\n";
}

void write_decay(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "decay", "Contact distance curve P(s)");
    out << "<p>P(s) is the number of contacts at distance s divided by the number of bin pairs at that distance. "
           "The decay exponent is the slope of log10 P(s) against log10 s for s from "
        << bp_text(options.fit_min) << " to " << bp_text(options.fit_max)
        << ", fitted to ten groups of distances equally spaced on the log scale (ordinary least squares).</p>\n";
    std::vector<Row> rows;
    rows.push_back({"Genome", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return fixed(m.decay.exponent, 3);
                    })});
    std::vector<std::string> names;
    for (const Sample& sample : samples) {
        for (const ChromosomeSummary& c : sample.matrix->chromosomes) {
            if (std::find(names.begin(), names.end(), c.name) == names.end()) {
                names.push_back(c.name);
            }
        }
    }
    for (const std::string& name : names) {
        rows.push_back({name, per_sample(samples, [&name](const MatrixStats& m, const Sample&) {
                            for (const ChromosomeSummary& c : m.chromosomes) {
                                if (c.name == name) {
                                    return fixed(c.decay.exponent, 3);
                                }
                            }
                            return std::string("n/a");
                        })});
    }
    out << "<h3>Decay exponent</h3>\n";
    stat_table(out, samples, rows, "Decay exponents");

    for (const Sample& sample : samples) {
        const MatrixStats& m = *sample.matrix;
        if (!m.pixel_statistics) {
            continue;
        }
        double low = 1e300, high = -1e300;
        for (std::size_t slot = 0; slot < m.decay.pairs.size(); ++slot) {
            if (m.decay.pairs[slot] > 0 && m.decay.contacts[slot] > 0) {
                const double p = std::log10(m.decay.contacts[slot] / m.decay.pairs[slot]);
                low = std::min(low, p);
                high = std::max(high, p);
            }
        }
        if (!(high > low)) {
            continue;
        }
        low = std::floor(low);
        high = std::ceil(high);
        out << "<h3>" << sample.tag << ": P(s), genome-wide</h3>\n"
            << "<p class=\"note\">Bar length: log10 P(s) from " << fixed(low, 0) << " to " << fixed(high, 0)
            << ". Darker bars lie in the fitted range.</p>\n";
        table_open(out, "P(s) of " + sample.tag);
        out << "<colgroup><col width=\"36%\"><col width=\"38%\"><col width=\"26%\"></colgroup>\n"
            << "<tr><th align=\"left\">Distance</th><th></th><th align=\"right\">P(s)</th></tr>\n";
        for (std::size_t slot = 0; slot < m.decay.pairs.size(); ++slot) {
            if (m.decay.pairs[slot] <= 0) {
                continue;
            }
            const double p = m.decay.contacts[slot] / m.decay.pairs[slot];
            const std::int64_t lo = decay_slot_lower_edge(slot);
            const std::int64_t hi = slot + 1 < kDecaySlots ? decay_slot_lower_edge(slot + 1) : lo;
            const std::string label = slot == 0 ? "0 (diagonal)" : bp_text(lo) + " to " + bp_text(hi);
            const bool in_fit = m.decay.mean_distance[slot] >= static_cast<double>(options.fit_min) &&
                                m.decay.mean_distance[slot] <= static_cast<double>(options.fit_max);
            out << "<tr><td>" << esc(label) << "</td>";
            track_cell(out, {{p > 0 ? 100.0 * (std::log10(p) - low) / (high - low) : 0.0, in_fit ? kBlue : kSky}});
            out << "<td align=\"right\">" << (p > 0 ? printf_string("%.3g", p) : std::string("0")) << "</td></tr>\n";
        }
        out << "</table>\n";
    }
}

void write_map_resolution(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "resolution", "Map resolution");
    out << "<p>After Rao et al. (2014), the map resolution is the smallest bin size at which "
        << percent(options.rao_fraction, 0) << " of the bins hold at least " << short_count(options.rao_min_contacts)
        << " contacts. The contacts of a bin are its row sum of the symmetric matrix. The fraction is taken over all bins "
           "of the chromosomes in the file, which includes gaps and unmappable regions; the second column uses the "
           "bins with at least one contact only. For an mcool file the resolutions present in the file are used. For a "
           "single cool file the bins are aggregated to multiples of the bin size (within each chromosome).</p>\n";
    std::set<const FileStats*> done;
    for (const Sample& sample : samples) {
        if (!done.insert(sample.file).second || !sample.matrix->pixel_statistics) {
            continue;
        }
        const MapResolution& map = sample.file->map_resolution;
        std::string tags;
        for (const Sample& other : samples) {
            if (other.file == sample.file) {
                tags += (tags.empty() ? "" : ", ") + other.tag;
            }
        }
        out << "<h3>" << tags << ": " << esc(sample.file->path) << "</h3>\n<p>Method: "
            << (map.method == "resolutions" ? "resolutions present in the file" : "aggregation of the bins of the file")
            << ". Map resolution: "
            << (map.bin_size > 0 ? "<b>" + bp_text(map.bin_size) + "</b>" : std::string("not reached by any bin size tried"))
            << " (bins with contacts only: "
            << (map.bin_size_nonempty > 0 ? bp_text(map.bin_size_nonempty) : std::string("not reached")) << ").";
        if (map.finest_reached) {
            out << " The smallest bin size of the file already reaches the fraction, so the map resolution is that "
                   "bin size or finer.";
        }
        out << "</p>\n";
        table_open(out, "Fraction of bins with enough contacts");
        out << "<colgroup><col width=\"20%\"><col width=\"40%\"><col width=\"20%\"><col width=\"20%\"></colgroup>\n"
            << "<tr><th align=\"left\">Bin size</th><th align=\"left\">Bins with enough contacts</th>"
               "<th align=\"right\">All bins</th><th align=\"right\">With contacts</th></tr>\n";
        for (const LadderStep& step : map.steps) {
            out << "<tr><td>" << bp_text(step.bin_size) << "</td>";
            track_cell(out, {{std::isnan(step.fraction_all) ? 0.0 : 100.0 * step.fraction_all,
                              step.reached ? kGreen : kGrey}});
            out << "<td align=\"right\">" << percent(step.fraction_all, 1) << "</td><td align=\"right\">"
                << percent(step.fraction_nonempty, 1) << "</td></tr>\n";
        }
        out << "</table>\n";
    }
}

void write_coverage(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "coverage", "Coverage per bin");
    out << "<p>The coverage of a bin is its row sum of the symmetric matrix of raw counts. Quantiles are taken over "
           "the bins with at least one contact. The MAD filter is the filter of hicCorrectMatrix: the coverage "
           "without the diagonal, its median absolute deviation over all bins, and the bins whose modified "
           "z-score lies outside "
        << fixed(options.mad_lower, 1) << " to " << fixed(options.mad_upper, 1)
        << " (--filterThreshold). The bins without contacts are among the removed bins.</p>\n";
    std::vector<Row> rows;
    rows.push_back({"Empty bins", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.empty_fraction, 2);
                    })});
    const auto quantile_row = [&](const std::string& label, double Quantiles::*member) {
        rows.push_back({label, per_sample(samples, [member](const MatrixStats& m, const Sample&) {
                            return short_count(m.coverage.*member);
                        })});
    };
    quantile_row("Minimum", &Quantiles::min);
    quantile_row("1st percentile", &Quantiles::p01);
    quantile_row("5th percentile", &Quantiles::p05);
    quantile_row("25th percentile", &Quantiles::p25);
    quantile_row("Median", &Quantiles::median);
    quantile_row("75th percentile", &Quantiles::p75);
    quantile_row("95th percentile", &Quantiles::p95);
    quantile_row("99th percentile", &Quantiles::p99);
    quantile_row("Maximum", &Quantiles::max);
    rows.push_back({"MAD filter: bins removed", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.pixel_statistics ? percent(m.mad_removed_fraction, 2) : std::string("n/a");
                    })});
    rows.push_back({"of which empty bins", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.pixel_statistics
                                   ? percent(fraction_of(static_cast<double>(m.mad_removed_empty),
                                                         static_cast<double>(m.bins)), 2)
                                   : std::string("n/a");
                    })});
    stat_table(out, samples, rows, "Coverage statistics");

    for (const Sample& sample : samples) {
        const MatrixStats& m = *sample.matrix;
        if (!m.pixel_statistics || m.chromosomes.empty()) {
            continue;
        }
        out << "<h3>" << sample.tag << ": coverage along the chromosomes</h3>\n"
            << "<p class=\"note\">Each chromosome is cut into equal parts; the colour is the mean coverage of the "
               "part relative to the median coverage of the genome (white 0, dark blue 2 or more, grey no contacts).</p>\n";
        table_open(out, "Coverage along the chromosomes of " + sample.tag);
        out << "<colgroup><col width=\"15%\"><col width=\"12%\"><col width=\"12%\"><col width=\"61%\"></colgroup>\n"
            << "<tr><th align=\"left\">Chr.</th><th align=\"right\">Empty</th><th align=\"right\">MAD</th>"
               "<th align=\"left\">Coverage</th></tr>\n";
        for (const ChromosomeSummary& c : m.chromosomes) {
            out << "<tr><td>" << esc(c.name) << "</td><td align=\"right\">" << percent(c.empty_fraction, 0)
                << "</td><td align=\"right\">" << percent(c.mad_removed_fraction, 0) << "</td>"
                << "<td class=\"track\"><table class=\"bar\" width=\"100%\" cellspacing=\"0\" cellpadding=\"0\" border=\"0\"><tr>";
            for (const double value : c.strip) {
                const std::string colour = value < 0 ? "#DDDDDD" : sequential_colour(value / 2.0);
                out << "<td bgcolor=\"" << colour << "\" height=\"12\" title=\"" << fixed(std::max(value, 0.0), 2)
                    << "\">&nbsp;</td>";
            }
            out << "</tr></table></td></tr>\n";
        }
        out << "</table>\n";
    }
}

void write_balancing(std::ostream& out, const std::vector<Sample>& samples) {
    section_open(out, "balancing", "Balancing");
    out << "<p>The weights are the bin column of the cool file (cooler balance, hicCorrectMatrix). The attributes "
           "of the weight dataset are those cooler stores. The marginals are the row sums of the balanced matrix "
           "(counts multiplied by the weights of both bins) over the bins with a finite weight; a balanced matrix "
           "has equal marginals, so a large coefficient of variation shows a matrix that is not balanced or did "
           "not converge. Matrices balanced on cis contacts only have equal cis marginals, not equal total marginals.</p>\n";
    std::set<std::string> attribute_names;
    for (const Sample& sample : samples) {
        for (const auto& entry : sample.matrix->balance.attributes) {
            attribute_names.insert(entry.first);
        }
    }
    std::vector<Row> rows;
    rows.push_back({"Weight column present", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return std::string(m.balance.has_weight_column ? "yes" : "no");
                    })});
    rows.push_back({"Weights usable", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return std::string(m.balance.usable ? "yes" : "no");
                    })});
    rows.push_back({"Bins with NaN weight", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return percent(m.balance.nan_fraction, 2);
                    })});
    for (const std::string& name : attribute_names) {
        rows.push_back({name, per_sample(samples, [&name](const MatrixStats& m, const Sample&) {
                            const auto found = m.balance.attributes.find(name);
                            return found == m.balance.attributes.end() ? std::string("n/a") : found->second;
                        })});
    }
    rows.push_back({"Mean of the marginals", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.balance.has_marginals ? printf_string("%.4g", m.balance.marginal_mean) : std::string("n/a");
                    })});
    rows.push_back({"Coefficient of variation", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.balance.has_marginals ? printf_string("%.4g", m.balance.marginal_cv) : std::string("n/a");
                    })});
    stat_table(out, samples, rows, "Balancing state");
}

void write_trans(std::ostream& out, const std::vector<Sample>& samples, const Options& options) {
    section_open(out, "trans", "Chromosome pairs, observed over expected trans contacts");
    out << "<p>For a pair of chromosomes the expected number of trans contacts is the total number of trans "
           "contacts between the listed chromosomes, distributed in proportion to the product of the two "
           "chromosomes' trans marginals (their share of all trans contact ends), normalised so that the expected "
           "values of all pairs add up to the observed total. Pairs with an observed over expected ratio above "
        << fixed(options.trans_threshold, 2)
        << " are listed. Only the longest chromosomes are listed (--maxChromosomes).</p>\n";
    for (const Sample& sample : samples) {
        const MatrixStats& m = *sample.matrix;
        const TransTable& t = m.trans_table;
        if (!m.pixel_statistics || t.names.size() < 2 || t.total <= 0) {
            continue;
        }
        const std::size_t k = t.names.size();
        out << "<h3>" << sample.tag << "</h3>\n"
            << "<table class=\"legend\" cellspacing=\"0\" cellpadding=\"2\" border=\"0\"><tr>"
            << "<td width=\"14\" bgcolor=\"" << diverging_colour(-1) << "\">&nbsp;</td><td>0.25 or less</td>"
            << "<td width=\"14\" bgcolor=\"" << diverging_colour(0) << "\" >&nbsp;</td><td>1</td>"
            << "<td width=\"14\" bgcolor=\"" << diverging_colour(1) << "\">&nbsp;</td><td>4 or more</td></tr></table>\n";
        out << "<table class=\"heat\" width=\"100%\" border=\"0\" cellspacing=\"1\" cellpadding=\"0\" "
               "summary=\"Observed over expected trans contacts of "
            << sample.tag << "\">\n<tr><td width=\"8%\"></td>";
        const bool numbers = k <= 12;
        const auto short_name = [](const std::string& name) {
            return name.size() > 3 && name.compare(0, 3, "chr") == 0 ? name.substr(3) : name;
        };
        for (std::size_t b = 0; b < k; ++b) {
            out << "<th>" << esc(short_name(t.names[b])) << "</th>";
        }
        out << "</tr>\n";
        for (std::size_t a = 0; a < k; ++a) {
            out << "<tr><th>" << esc(short_name(t.names[a])) << "</th>";
            for (std::size_t b = 0; b < k; ++b) {
                if (a == b) {
                    out << "<td bgcolor=\"#EEEEEE\">&nbsp;</td>";
                    continue;
                }
                const double expected = t.expected[a * k + b];
                const double observed = t.observed[a * k + b];
                const double ratio = expected > 0 ? observed / expected : kNaN;
                const std::string colour =
                    std::isnan(ratio) ? "#EEEEEE"
                                      : diverging_colour(std::log2(std::max(ratio, 1e-6)) / 2.0);
                out << "<td bgcolor=\"" << colour << "\" align=\"center\" title=\"" << esc(t.names[a]) << " and "
                    << esc(t.names[b]) << ": " << fixed(ratio, 2) << "\">"
                    << (numbers && !std::isnan(ratio) ? fixed(ratio, 1) : std::string("&nbsp;")) << "</td>";
            }
            out << "</tr>\n";
        }
        out << "</table>\n";
        if (t.flagged.empty()) {
            out << "<p>No pair above the threshold.</p>\n";
        } else {
            out << "<p>Pairs above the threshold:</p>\n";
            table_open(out, "Flagged chromosome pairs of " + sample.tag);
            out << "<colgroup><col width=\"28%\"><col width=\"24%\"><col width=\"24%\"><col width=\"24%\"></colgroup>\n"
                << "<tr><th align=\"left\">Pair</th><th align=\"right\">Observed</th><th align=\"right\">Expected</th>"
                   "<th align=\"right\">Ratio</th></tr>\n";
            std::size_t shown = 0;
            for (const TransTable::Pair& p : t.flagged) {
                if (++shown > 40) {
                    break;
                }
                out << "<tr><td>" << esc(t.names[p.a]) << " and " << esc(t.names[p.b]) << "</td><td align=\"right\">"
                    << short_count(p.observed) << "</td><td align=\"right\">" << short_count(p.expected)
                    << "</td><td align=\"right\">" << fixed(p.ratio, 2) << "</td></tr>\n";
            }
            out << "</table>\n";
            if (t.flagged.size() > 40) {
                out << "<p class=\"note\">" << t.flagged.size() - 40 << " further pairs are in the table file.</p>\n";
            }
        }
    }
}

void write_checks(std::ostream& out, const std::vector<Sample>& samples) {
    section_open(out, "checks", "Further checks");
    out << "<p>The X over autosome ratio compares the mean coverage of the non-empty bins of chromosome X with that "
           "of the autosomes (a male genome gives about 0.5, a female genome about 1 for equal copy number and "
           "mappability). The chromosome Y and chromosome M fractions are shares of the coverage.</p>\n";
    std::vector<Row> rows;
    rows.push_back({"X over autosome coverage", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.checks.has_x ? fixed(m.checks.x_autosome_ratio, 2) : std::string("no chrX");
                    })});
    rows.push_back({"Chromosome Y present", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return std::string(m.checks.has_y ? "yes" : "no");
                    })});
    rows.push_back({"Chromosome Y share of coverage", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.checks.has_y ? percent(m.checks.y_fraction, 3) : std::string("n/a");
                    })});
    rows.push_back({"Chromosome M present", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return std::string(m.checks.has_m ? "yes" : "no");
                    })});
    rows.push_back({"Chromosome M share of coverage", per_sample(samples, [](const MatrixStats& m, const Sample&) {
                        return m.checks.has_m ? percent(m.checks.m_fraction, 3) : std::string("n/a");
                    })});
    stat_table(out, samples, rows, "Further checks");
}

const char* const kStyle =
    "body { font-family: Georgia, 'Times New Roman', serif; color: #111; background: #fff; margin: 0 auto; "
    "padding: 8px 16px 40px 16px; max-width: 920px; line-height: 1.4; font-size: 15px; }\n"
    "h1 { font-size: 1.5em; margin: 0.6em 0 0.2em 0; }\n"
    "h2 { font-size: 1.25em; margin: 1.6em 0 0.4em 0; border-bottom: 1px solid #888; }\n"
    "h3 { font-size: 1.05em; margin: 1.2em 0 0.3em 0; overflow-wrap: anywhere; }\n"
    "h4 { font-size: 1em; margin: 1em 0 0.3em 0; }\n"
    "table { border-collapse: collapse; table-layout: fixed; width: 100%; margin: 0.4em 0 0.8em 0; }\n"
    "table.data td, table.data th { border: 1px solid #bbb; padding: 2px 5px; font-size: 13px; "
    "overflow-wrap: anywhere; vertical-align: middle; }\n"
    "table.data th { background: #eee; font-weight: bold; }\n"
    "table.data tr.group td { background: #f4f4f4; }\n"
    "table.bar, table.legend { border: 0; margin: 0; width: 100%; }\n"
    "table.legend { width: auto; font-size: 13px; margin: 0.2em 0 0.4em 0; table-layout: auto; }\n"
    "table.legend td { padding: 1px 8px 1px 2px; }\n"
    "td.track { padding: 2px 5px; }\n"
    "td.track table td { padding: 0; height: 12px; font-size: 1px; line-height: 1px; border: 0; }\n"
    "td.name { overflow-wrap: anywhere; }\n"
    "table.heat { table-layout: fixed; }\n"
    "table.heat th { font-size: 10px; font-weight: normal; padding: 0; overflow: hidden; }\n"
    "table.heat td { font-size: 9px; height: 14px; padding: 0; overflow: hidden; }\n"
    "p.note { color: #444; font-size: 13px; margin: 0.2em 0 0.4em 0; }\n"
    "ul { padding-left: 1.3em; }\n";

void write_notes(std::ostream& out) {
    section_open(out, "notes", "Definitions and approximations");
    out << "<ul>\n"
        << "<li>The statistics are computed from the stored pixels of the matrix, which hold the upper triangle; "
           "every stored pixel is counted once and counts are not balanced.</li>\n"
        << "<li>Distances between bins are the distances between bin starts. A pixel at the offset d holds read "
           "pairs from distances between d - 1 and d + 1 bin sizes, so a distance threshold (short range, the "
           "fitted range of P(s)) is exact only up to one bin size, and finer bins give a closer match to "
           "the read-level numbers of QC.log, which use the distance between read starts.</li>\n"
        << "<li>The main diagonal holds pairs inside one bin; they are closer than one bin size and are counted as "
           "short range. P(s) at the offset 0 is the number of contacts on the diagonal divided by the number of bins.</li>\n"
        << "<li>The numbers of bin pairs at a distance include all bins of the chromosome, empty or not.</li>\n"
        << "<li>The statistics of a matrix with bins of different sizes (for example restriction fragments) are "
           "limited to those stored in the metadata.</li>\n"
        << "</ul>\n";
}

}  // namespace

std::string render_html(const std::vector<FileStats>& files, const Options& options,
                        const std::string& version) {
    const std::vector<Sample> samples = samples_of(files);
    std::ostringstream out;
    out << "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
           "<title>hicMatrixQC report</title>\n<style><!--\n"
        << kStyle << "--></style>\n</head>\n<body>\n";
    out << "<h1>hicMatrixQC report</h1>\n<p class=\"note\">HiCExplorer " << esc(version) << ", " << iso_date()
        << ". " << samples.size() << " sample" << (samples.size() == 1 ? "" : "s") << ".</p>\n";
    out << "<p class=\"note\">Contents: <a href=\"#samples\">samples</a>, <a href=\"#summary\">summary</a>, "
           "<a href=\"#library\">library statistics</a>, <a href=\"#contacts\">contacts</a>, "
           "<a href=\"#decay\">P(s)</a>, <a href=\"#resolution\">map resolution</a>, "
           "<a href=\"#coverage\">coverage</a>, <a href=\"#balancing\">balancing</a>, "
           "<a href=\"#trans\">trans pairs</a>, <a href=\"#checks\">checks</a>, "
           "<a href=\"#notes\">definitions</a>.</p>\n";
    write_samples(out, samples);
    write_summary(out, samples, options);
    write_read_qc(out, samples);
    write_contacts(out, samples, options);
    write_decay(out, samples, options);
    write_map_resolution(out, samples, options);
    write_coverage(out, samples, options);
    write_balancing(out, samples);
    write_trans(out, samples, options);
    write_checks(out, samples);
    write_notes(out);
    out << "</body>\n</html>\n";
    return out.str();
}

}  // namespace hicx::matrixqc
