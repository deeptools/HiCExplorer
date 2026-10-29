#include "hicx/matrix_qc.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>

#include "hicx/build_matrix.hpp"
#include "hicx/correct_ops.hpp"
#include "hicx/hdf5_util.hpp"

namespace hicx::matrixqc {

namespace {

constexpr std::size_t kGroups = 10;  // groups of distances inside the fit window
constexpr std::uint8_t kNoGroup = 255;

struct OffsetClass {
    std::uint8_t slot = 0;
    std::uint8_t group = kNoGroup;
};

// Everything one worker accumulates. All sums are doubles of contact counts;
// raw cooler counts are integers, whose sums are exact, so the merge order
// does not change a result.
struct ThreadState {
    std::vector<double> coverage_off_diagonal;  // per bin, without the diagonal pixel
    std::vector<double> marginal_off_diagonal;  // the same with balanced counts
    std::vector<std::pair<std::int64_t, double>> diagonal;           // (bin, count)
    std::vector<std::pair<std::int64_t, double>> balanced_diagonal;  // (bin, balanced count)
    double total = 0, cis = 0, trans = 0, short_range = 0, long_range = 0, diag0 = 0, diag1 = 0;
    std::int64_t pixels = 0;
    std::vector<double> cis_by_chromosome;
    std::vector<double> slots_all;  // kDecaySlots
    std::vector<double> slots_sel;  // k * kDecaySlots
    std::vector<double> groups_all;  // kGroups
    std::vector<double> groups_sel;  // k * kGroups
    std::vector<double> trans_observed;  // k * k
};

double quantile_sorted(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) {
        return kNaN;
    }
    const double position = q * static_cast<double>(sorted.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(position));
    const auto upper = static_cast<std::size_t>(std::ceil(position));
    return sorted[lower] + (sorted[upper] - sorted[lower]) * (position - static_cast<double>(lower));
}

Quantiles quantiles_of(std::vector<double> values) {
    Quantiles result;
    if (values.empty()) {
        return result;
    }
    std::sort(values.begin(), values.end());
    result.min = values.front();
    result.p01 = quantile_sorted(values, 0.01);
    result.p05 = quantile_sorted(values, 0.05);
    result.p25 = quantile_sorted(values, 0.25);
    result.median = quantile_sorted(values, 0.50);
    result.p75 = quantile_sorted(values, 0.75);
    result.p95 = quantile_sorted(values, 0.95);
    result.p99 = quantile_sorted(values, 0.99);
    result.max = values.back();
    return result;
}

std::string attribute_text(const h5::AttributeValue& value) {
    if (const auto* integer = std::get_if<std::int64_t>(&value)) {
        return std::to_string(*integer);
    }
    if (const auto* real = std::get_if<double>(&value)) {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.10g", *real);
        return buffer;
    }
    return std::get<std::string>(value);
}

std::vector<OffsetClass> classify_offsets(const Layout& layout, const Options& options,
                                          std::int64_t longest_chromosome_bins) {
    std::vector<OffsetClass> table(static_cast<std::size_t>(std::max<std::int64_t>(longest_chromosome_bins, 1)));
    const double log_min = std::log10(static_cast<double>(options.fit_min));
    const double log_max = std::log10(static_cast<double>(options.fit_max));
    for (std::size_t d = 0; d < table.size(); ++d) {
        const std::int64_t s = static_cast<std::int64_t>(d) * layout.resolution;
        OffsetClass& entry = table[d];
        entry.slot = s == 0 ? std::uint8_t{0}
                            : static_cast<std::uint8_t>(1 + qc_distance_bin(s));
        if (s >= options.fit_min && s <= options.fit_max && s > 0) {
            const double position =
                (std::log10(static_cast<double>(s)) - log_min) / (log_max - log_min);
            entry.group = static_cast<std::uint8_t>(
                std::clamp<long>(static_cast<long>(std::floor(position * kGroups)), 0,
                                 static_cast<long>(kGroups) - 1));
        }
    }
    return table;
}

}  // namespace

std::int64_t decay_slot_lower_edge(std::size_t slot) {
    return slot == 0 ? 0 : qc_distance_edge(static_cast<int>(slot) - 1);
}

bool fit_line(const std::vector<double>& x, const std::vector<double>& y, double& slope,
              double& intercept) {
    if (x.size() < 3 || x.size() != y.size()) {
        return false;
    }
    const double n = static_cast<double>(x.size());
    const double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    const double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        sxx += (x[i] - mean_x) * (x[i] - mean_x);
        sxy += (x[i] - mean_x) * (y[i] - mean_y);
    }
    if (!(sxx > 0)) {
        return false;
    }
    slope = sxy / sxx;
    intercept = mean_y - slope * mean_x;
    return true;
}

ChromosomeClass classify_chromosome(const std::string& name) {
    std::string stem = name;
    if (stem.size() > 3 && (stem.compare(0, 3, "chr") == 0 || stem.compare(0, 3, "Chr") == 0 ||
                            stem.compare(0, 3, "CHR") == 0)) {
        stem = stem.substr(3);
    }
    if (stem == "X") {
        return ChromosomeClass::X;
    }
    if (stem == "Y") {
        return ChromosomeClass::Y;
    }
    if (stem == "M" || stem == "MT" || stem == "Mt") {
        return ChromosomeClass::M;
    }
    if (!stem.empty() && std::all_of(stem.begin(), stem.end(),
                                     [](unsigned char c) { return std::isdigit(c) != 0; })) {
        return ChromosomeClass::Autosome;
    }
    if (stem.size() == 2 && std::isdigit(static_cast<unsigned char>(stem[0])) != 0 &&
        (stem[1] == 'L' || stem[1] == 'R')) {
        return ChromosomeClass::Autosome;  // Drosophila arms
    }
    return ChromosomeClass::Other;
}

std::vector<std::int64_t> ladder_factors(const Layout& layout) {
    static const std::int64_t candidates[] = {1,  2,  3,  4,   5,   6,   8,   10,  12,
                                              15, 20, 25, 30,  40,  50,  60,  80,  100,
                                              150, 200, 250, 300, 400, 500, 1000};
    std::int64_t longest = 0;
    for (const std::int64_t length : layout.lengths) {
        longest = std::max(longest, length);
    }
    std::vector<std::int64_t> factors;
    for (const std::int64_t factor : candidates) {
        if (factor == 1 || (factor * layout.resolution <= 10000000 &&
                            factor * layout.resolution <= longest)) {
            factors.push_back(factor);
        }
    }
    return factors;
}

LadderStep aggregate_ladder_step(const std::vector<double>& coverage, const Layout& layout,
                                 std::int64_t factor, const Options& options) {
    LadderStep step;
    step.bin_size = factor * layout.resolution;
    std::int64_t groups = 0, enough = 0, nonempty = 0;
    for (std::size_t c = 0; c < layout.chromosomes(); ++c) {
        const std::int64_t first = layout.first_bin[c];
        const std::int64_t last = layout.first_bin[c + 1];
        for (std::int64_t start = first; start < last; start += factor) {
            const std::int64_t stop = std::min(last, start + factor);
            double sum = 0;
            for (std::int64_t b = start; b < stop; ++b) {
                sum += coverage[static_cast<std::size_t>(b)];
            }
            ++groups;
            if (sum > 0) {
                ++nonempty;
            }
            if (sum >= options.rao_min_contacts) {
                ++enough;
            }
        }
    }
    if (groups > 0) {
        step.fraction_all = static_cast<double>(enough) / static_cast<double>(groups);
    }
    if (nonempty > 0) {
        step.fraction_nonempty = static_cast<double>(enough) / static_cast<double>(nonempty);
    }
    step.reached = step.fraction_all >= options.rao_fraction;
    step.reached_nonempty = step.fraction_nonempty >= options.rao_fraction;
    return step;
}

MapResolution pick_map_resolution(std::string method, std::vector<LadderStep> steps) {
    MapResolution result;
    result.method = std::move(method);
    std::sort(steps.begin(), steps.end(),
              [](const LadderStep& a, const LadderStep& b) { return a.bin_size < b.bin_size; });
    for (const LadderStep& step : steps) {
        if (result.bin_size == 0 && step.reached) {
            result.bin_size = step.bin_size;
        }
        if (result.bin_size_nonempty == 0 && step.reached_nonempty) {
            result.bin_size_nonempty = step.bin_size;
        }
    }
    result.finest_reached = !steps.empty() && steps.front().reached;
    result.steps = std::move(steps);
    return result;
}

MatrixStats compute_matrix_stats(const Layout& layout, const PixelSource& source,
                                 const std::optional<Correction>& correction,
                                 const Options& options) {
    if (layout.resolution <= 0 || layout.names.empty()) {
        throw std::invalid_argument("the matrix needs fixed-size bins and a chromosome table");
    }
    if (!(options.fit_min > 0 && options.fit_min < options.fit_max)) {
        throw std::invalid_argument("the decay range needs 0 < low < high");
    }
    const std::size_t bins = static_cast<std::size_t>(layout.bins());
    const std::size_t chromosomes = layout.chromosomes();
    const unsigned threads = std::max(1U, options.threads);

    // The bin -> chromosome lookup.
    std::vector<std::int32_t> chromosome_of(bins);
    std::int64_t longest_bins = 0;
    for (std::size_t c = 0; c < chromosomes; ++c) {
        std::fill(chromosome_of.begin() + layout.first_bin[c],
                  chromosome_of.begin() + layout.first_bin[c + 1], static_cast<std::int32_t>(c));
        longest_bins = std::max(longest_bins, layout.first_bin[c + 1] - layout.first_bin[c]);
    }

    // The listed chromosomes: the longest ones, in genome order.
    std::vector<std::size_t> order(chromosomes);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return layout.lengths[a] > layout.lengths[b];
    });
    order.resize(std::min(order.size(), options.max_chromosomes));
    std::sort(order.begin(), order.end());
    const std::size_t listed = order.size();
    std::vector<std::int32_t> listed_index(chromosomes, -1);
    for (std::size_t i = 0; i < listed; ++i) {
        listed_index[order[i]] = static_cast<std::int32_t>(i);
    }

    const std::vector<OffsetClass> offsets = classify_offsets(layout, options, longest_bins);
    const std::int64_t resolution = layout.resolution;
    const std::int64_t cutoff = options.short_range_cutoff;

    std::vector<ThreadState> states(threads);
    for (ThreadState& state : states) {
        state.coverage_off_diagonal.assign(bins, 0.0);
        if (correction.has_value()) {
            state.marginal_off_diagonal.assign(bins, 0.0);
        }
        state.cis_by_chromosome.assign(chromosomes, 0.0);
        state.slots_all.assign(kDecaySlots, 0.0);
        state.slots_sel.assign(listed * kDecaySlots, 0.0);
        state.groups_all.assign(kGroups, 0.0);
        state.groups_sel.assign(listed * kGroups, 0.0);
        state.trans_observed.assign(listed * listed, 0.0);
    }

    const bool divide = correction.has_value() && correction->op == '/';
    const ChunkVisitor visit = [&](unsigned thread, const PixelChunk& chunk) {
        ThreadState& st = states[thread];
        const std::size_t n = chunk.bin1.size();
        for (std::size_t k = 0; k < n; ++k) {
            std::int64_t b1 = chunk.bin1[k];
            std::int64_t b2 = chunk.bin2[k];
            const double c = chunk.count[k];
            if (!(c != 0.0) || std::isnan(c)) {
                continue;
            }
            if (b1 > b2) {
                std::swap(b1, b2);
            }
            if (b1 < 0 || b2 >= static_cast<std::int64_t>(bins)) {
                throw std::runtime_error("a pixel lies outside the bin table");
            }
            ++st.pixels;
            st.total += c;
            if (b1 == b2) {
                st.diagonal.emplace_back(b1, c);
            } else {
                st.coverage_off_diagonal[static_cast<std::size_t>(b1)] += c;
                st.coverage_off_diagonal[static_cast<std::size_t>(b2)] += c;
            }
            if (correction.has_value()) {
                const double factor = correction->factors[static_cast<std::size_t>(b1)] *
                                      correction->factors[static_cast<std::size_t>(b2)];
                const double value = divide ? c / factor : c * factor;
                if (!std::isnan(value)) {
                    if (b1 == b2) {
                        st.balanced_diagonal.emplace_back(b1, value);
                    } else {
                        st.marginal_off_diagonal[static_cast<std::size_t>(b1)] += value;
                        st.marginal_off_diagonal[static_cast<std::size_t>(b2)] += value;
                    }
                }
            }
            const std::int32_t ch1 = chromosome_of[static_cast<std::size_t>(b1)];
            const std::int32_t ch2 = chromosome_of[static_cast<std::size_t>(b2)];
            if (ch1 == ch2) {
                const std::int64_t d = b2 - b1;
                st.cis += c;
                st.cis_by_chromosome[static_cast<std::size_t>(ch1)] += c;
                const OffsetClass& cls = offsets[static_cast<std::size_t>(d)];
                st.slots_all[cls.slot] += c;
                const std::int32_t sel = listed_index[static_cast<std::size_t>(ch1)];
                if (sel >= 0) {
                    st.slots_sel[static_cast<std::size_t>(sel) * kDecaySlots + cls.slot] += c;
                }
                if (cls.group != kNoGroup) {
                    st.groups_all[cls.group] += c;
                    if (sel >= 0) {
                        st.groups_sel[static_cast<std::size_t>(sel) * kGroups + cls.group] += c;
                    }
                }
                if (d == 0) {
                    st.diag0 += c;
                } else if (d == 1) {
                    st.diag1 += c;
                }
                if (d * resolution < cutoff) {
                    st.short_range += c;
                } else {
                    st.long_range += c;
                }
            } else {
                st.trans += c;
                const std::int32_t a = listed_index[static_cast<std::size_t>(ch1)];
                const std::int32_t b = listed_index[static_cast<std::size_t>(ch2)];
                if (a >= 0 && b >= 0) {
                    st.trans_observed[static_cast<std::size_t>(a) * listed + static_cast<std::size_t>(b)] += c;
                    st.trans_observed[static_cast<std::size_t>(b) * listed + static_cast<std::size_t>(a)] += c;
                }
            }
        }
    };
    source(threads, visit);

    // --- merge the workers, in thread order --------------------------------------
    MatrixStats stats;
    stats.resolution = resolution;
    stats.bins = static_cast<std::int64_t>(bins);
    stats.chromosomes_total = chromosomes;
    stats.pixel_statistics = true;
    stats.short_range_cutoff = cutoff;
    std::vector<double> coverage_off_diagonal = std::move(states[0].coverage_off_diagonal);
    std::vector<double> marginal = std::move(states[0].marginal_off_diagonal);
    std::vector<double> diagonal(bins, 0.0);
    std::vector<double> balanced_diagonal(correction.has_value() ? bins : 0, 0.0);
    std::vector<double> cis_by_chromosome(chromosomes, 0.0);
    std::vector<double> slots_all(kDecaySlots, 0.0), slots_sel(listed * kDecaySlots, 0.0);
    std::vector<double> groups_all(kGroups, 0.0), groups_sel(listed * kGroups, 0.0);
    std::vector<double> trans_observed(listed * listed, 0.0);
    for (unsigned t = 0; t < threads; ++t) {
        ThreadState& st = states[t];
        if (t > 0) {
            for (std::size_t b = 0; b < bins; ++b) {
                coverage_off_diagonal[b] += st.coverage_off_diagonal[b];
            }
            for (std::size_t b = 0; b < marginal.size(); ++b) {
                marginal[b] += st.marginal_off_diagonal[b];
            }
            st.coverage_off_diagonal = {};
            st.marginal_off_diagonal = {};
        }
        for (const auto& [bin, value] : st.diagonal) {
            diagonal[static_cast<std::size_t>(bin)] += value;
        }
        for (const auto& [bin, value] : st.balanced_diagonal) {
            balanced_diagonal[static_cast<std::size_t>(bin)] += value;
        }
        stats.total_contacts += st.total;
        stats.cis += st.cis;
        stats.trans += st.trans;
        stats.intra_short += st.short_range;
        stats.intra_long += st.long_range;
        stats.diagonal += st.diag0;
        stats.first_offdiagonal += st.diag1;
        stats.pixels += st.pixels;
        for (std::size_t i = 0; i < chromosomes; ++i) {
            cis_by_chromosome[i] += st.cis_by_chromosome[i];
        }
        for (std::size_t i = 0; i < slots_all.size(); ++i) {
            slots_all[i] += st.slots_all[i];
        }
        for (std::size_t i = 0; i < slots_sel.size(); ++i) {
            slots_sel[i] += st.slots_sel[i];
        }
        for (std::size_t i = 0; i < kGroups; ++i) {
            groups_all[i] += st.groups_all[i];
        }
        for (std::size_t i = 0; i < groups_sel.size(); ++i) {
            groups_sel[i] += st.groups_sel[i];
        }
        for (std::size_t i = 0; i < trans_observed.size(); ++i) {
            trans_observed[i] += st.trans_observed[i];
        }
    }
    states.clear();

    // --- coverage ------------------------------------------------------------------
    std::vector<double> coverage(bins);
    for (std::size_t b = 0; b < bins; ++b) {
        coverage[b] = coverage_off_diagonal[b] + diagonal[b];
    }
    std::vector<double> nonempty;
    nonempty.reserve(bins);
    for (const double value : coverage) {
        if (value > 0) {
            nonempty.push_back(value);
        }
    }
    stats.empty_fraction =
        bins == 0 ? kNaN : static_cast<double>(bins - nonempty.size()) / static_cast<double>(bins);
    stats.coverage = quantiles_of(nonempty);
    stats.coverage_median_nonempty = stats.coverage.median;

    // hicCorrectMatrix' filter: the coverage without the diagonal, its MAD,
    // and the bins whose modified z-score leaves [mad_lower, mad_upper].
    std::vector<char> mad_removed(bins, 0);
    if (!nonempty.empty() && std::any_of(coverage_off_diagonal.begin(), coverage_off_diagonal.end(),
                                         [](double v) { return v > 0; })) {
        const correct::Mad mad(coverage_off_diagonal);
        stats.mad_median = mad.median();
        stats.mad_value = mad.median_absolute_deviation();
        const std::vector<double>& z = mad.modified_z_scores();
        for (std::size_t b = 0; b < bins; ++b) {
            if (z[b] < options.mad_lower || z[b] > options.mad_upper) {
                mad_removed[b] = 1;
                ++stats.mad_removed;
                if (coverage[b] == 0) {
                    ++stats.mad_removed_empty;
                }
            }
        }
        stats.mad_removed_fraction =
            static_cast<double>(stats.mad_removed) / static_cast<double>(bins);
    }

    // --- the decay curves ----------------------------------------------------------
    const auto build_curve = [&](const std::vector<double>& slot_contacts,
                                 const std::vector<double>& group_contacts,
                                 const std::vector<std::size_t>& members) {
        DecayCurve curve;
        curve.contacts = slot_contacts;
        curve.pairs.assign(kDecaySlots, 0.0);
        curve.mean_distance.assign(kDecaySlots, 0.0);
        std::vector<double> distance_sum(kDecaySlots, 0.0);
        std::vector<double> group_pairs(kGroups, 0.0), group_distance(kGroups, 0.0);
        for (const std::size_t c : members) {
            const std::int64_t n = layout.first_bin[c + 1] - layout.first_bin[c];
            for (std::int64_t d = 0; d < n; ++d) {
                const OffsetClass& cls = offsets[static_cast<std::size_t>(d)];
                const double pairs = static_cast<double>(n - d);
                const double distance = static_cast<double>(d) * static_cast<double>(resolution);
                curve.pairs[cls.slot] += pairs;
                distance_sum[cls.slot] += pairs * distance;
                if (cls.group != kNoGroup) {
                    group_pairs[cls.group] += pairs;
                    group_distance[cls.group] += pairs * distance;
                }
            }
        }
        for (std::size_t slot = 0; slot < kDecaySlots; ++slot) {
            if (curve.pairs[slot] > 0) {
                curve.mean_distance[slot] = distance_sum[slot] / curve.pairs[slot];
            }
        }
        std::vector<double> x, y;
        for (std::size_t g = 0; g < kGroups; ++g) {
            if (group_pairs[g] > 0 && group_contacts[g] > 0) {
                x.push_back(std::log10(group_distance[g] / group_pairs[g]));
                y.push_back(std::log10(group_contacts[g] / group_pairs[g]));
            }
        }
        curve.fit_points = x.size();
        double slope = 0, intercept = 0;
        if (fit_line(x, y, slope, intercept)) {
            curve.exponent = slope;
            curve.intercept = intercept;
        }
        return curve;
    };
    {
        std::vector<std::size_t> all(chromosomes);
        std::iota(all.begin(), all.end(), 0);
        stats.decay = build_curve(slots_all, groups_all, all);
    }

    // --- per chromosome ----------------------------------------------------------------
    const double genome_median = stats.coverage.median;
    for (std::size_t i = 0; i < listed; ++i) {
        const std::size_t c = order[i];
        ChromosomeSummary summary;
        summary.name = layout.names[c];
        summary.length = layout.lengths[c];
        const std::int64_t first = layout.first_bin[c];
        const std::int64_t last = layout.first_bin[c + 1];
        summary.bins = last - first;
        summary.cis_contacts = cis_by_chromosome[c];
        std::vector<double> values;
        std::int64_t removed = 0;
        for (std::int64_t b = first; b < last; ++b) {
            const double value = coverage[static_cast<std::size_t>(b)];
            summary.coverage_sum += value;
            if (value > 0) {
                values.push_back(value);
            }
            removed += mad_removed[static_cast<std::size_t>(b)];
        }
        if (summary.bins > 0) {
            summary.empty_fraction =
                static_cast<double>(summary.bins - static_cast<std::int64_t>(values.size())) /
                static_cast<double>(summary.bins);
            if (!std::isnan(stats.mad_removed_fraction)) {
                summary.mad_removed_fraction =
                    static_cast<double>(removed) / static_cast<double>(summary.bins);
            }
        }
        summary.median_coverage = quantiles_of(values).median;
        // The strip: the mean coverage of equal parts of the chromosome.
        const std::int64_t cells =
            std::min<std::int64_t>(static_cast<std::int64_t>(options.strip_cells), summary.bins);
        for (std::int64_t j = 0; j < cells; ++j) {
            const std::int64_t from = first + summary.bins * j / cells;
            const std::int64_t to = first + summary.bins * (j + 1) / cells;
            double sum = 0;
            for (std::int64_t b = from; b < to; ++b) {
                sum += coverage[static_cast<std::size_t>(b)];
            }
            const double mean = to > from ? sum / static_cast<double>(to - from) : 0.0;
            summary.strip.push_back(mean > 0 && genome_median > 0 ? mean / genome_median : -1.0);
        }
        std::vector<double> slot_contacts(slots_sel.begin() + static_cast<std::ptrdiff_t>(i * kDecaySlots),
                                          slots_sel.begin() + static_cast<std::ptrdiff_t>((i + 1) * kDecaySlots));
        std::vector<double> group_contacts(groups_sel.begin() + static_cast<std::ptrdiff_t>(i * kGroups),
                                           groups_sel.begin() + static_cast<std::ptrdiff_t>((i + 1) * kGroups));
        summary.decay = build_curve(slot_contacts, group_contacts, {c});
        stats.chromosomes.push_back(std::move(summary));
    }

    // --- trans enrichment ----------------------------------------------------------------
    {
        TransTable& table = stats.trans_table;
        for (const std::size_t c : order) {
            table.names.push_back(layout.names[c]);
        }
        table.observed = trans_observed;
        table.expected.assign(listed * listed, kNaN);
        std::vector<double> marginals(listed, 0.0);
        double sum = 0;
        for (std::size_t a = 0; a < listed; ++a) {
            for (std::size_t b = 0; b < listed; ++b) {
                marginals[a] += trans_observed[a * listed + b];
            }
            sum += marginals[a];
        }
        table.total = sum / 2.0;
        double squares = 0;
        for (std::size_t a = 0; a < listed; ++a) {
            const double f = sum > 0 ? marginals[a] / sum : 0.0;
            squares += f * f;
        }
        if (table.total > 0 && squares < 1.0) {
            for (std::size_t a = 0; a < listed; ++a) {
                for (std::size_t b = 0; b < listed; ++b) {
                    if (a == b) {
                        continue;
                    }
                    const double fa = marginals[a] / sum;
                    const double fb = marginals[b] / sum;
                    table.expected[a * listed + b] = table.total * 2.0 * fa * fb / (1.0 - squares);
                }
            }
            for (std::size_t a = 0; a < listed; ++a) {
                for (std::size_t b = a + 1; b < listed; ++b) {
                    const double expected = table.expected[a * listed + b];
                    const double observed = table.observed[a * listed + b];
                    if (expected > 0 && observed / expected > options.trans_threshold) {
                        table.flagged.push_back({a, b, observed, expected, observed / expected});
                    }
                }
            }
            std::sort(table.flagged.begin(), table.flagged.end(),
                      [](const TransTable::Pair& x, const TransTable::Pair& y) {
                          return x.ratio > y.ratio;
                      });
        }
    }

    // --- cheap checks ------------------------------------------------------------------------
    {
        CheapChecks& checks = stats.checks;
        if (stats.total_contacts > 0) {
            checks.diagonal_fraction = stats.diagonal / stats.total_contacts;
            checks.first_offdiag_fraction = stats.first_offdiagonal / stats.total_contacts;
        }
        double sum_all = 0, sum_x = 0, sum_y = 0, sum_m = 0, sum_auto = 0;
        std::int64_t bins_x = 0, bins_auto = 0;
        for (std::size_t c = 0; c < chromosomes; ++c) {
            double sum = 0;
            std::int64_t occupied = 0;
            for (std::int64_t b = layout.first_bin[c]; b < layout.first_bin[c + 1]; ++b) {
                const double value = coverage[static_cast<std::size_t>(b)];
                sum += value;
                occupied += value > 0 ? 1 : 0;
            }
            sum_all += sum;
            switch (classify_chromosome(layout.names[c])) {
                case ChromosomeClass::X:
                    checks.has_x = true;
                    sum_x += sum;
                    bins_x += occupied;
                    break;
                case ChromosomeClass::Y:
                    checks.has_y = true;
                    sum_y += sum;
                    break;
                case ChromosomeClass::M:
                    checks.has_m = true;
                    sum_m += sum;
                    break;
                case ChromosomeClass::Autosome:
                    sum_auto += sum;
                    bins_auto += occupied;
                    break;
                case ChromosomeClass::Other:
                    break;
            }
        }
        if (checks.has_x && bins_x > 0 && bins_auto > 0 && sum_auto > 0) {
            checks.x_autosome_ratio = (sum_x / static_cast<double>(bins_x)) /
                                      (sum_auto / static_cast<double>(bins_auto));
        }
        if (sum_all > 0) {
            checks.y_fraction = checks.has_y ? sum_y / sum_all : kNaN;
            checks.m_fraction = checks.has_m ? sum_m / sum_all : kNaN;
        }
    }

    // --- balancing --------------------------------------------------------------------------------
    if (correction.has_value()) {
        BalanceInfo& balance = stats.balance;
        balance.has_weight_column = true;
        balance.usable = true;
        std::size_t not_finite = 0;
        double sum = 0, squares = 0;
        std::size_t count = 0;
        for (std::size_t b = 0; b < bins; ++b) {
            if (!std::isfinite(correction->factors[b])) {
                ++not_finite;
                continue;
            }
            const double value = marginal[b] + balanced_diagonal[b];
            sum += value;
            squares += value * value;
            ++count;
        }
        balance.nan_fraction = static_cast<double>(not_finite) / static_cast<double>(bins);
        if (count > 0) {
            const double mean = sum / static_cast<double>(count);
            const double variance = std::max(0.0, squares / static_cast<double>(count) - mean * mean);
            balance.has_marginals = true;
            balance.marginal_mean = mean;
            balance.marginal_cv = mean != 0 ? std::sqrt(variance) / std::fabs(mean) : kNaN;
        }
    }

    // --- map resolution ladder ---------------------------------------------------------------------------
    for (const std::int64_t factor : ladder_factors(layout)) {
        stats.ladder.push_back(aggregate_ladder_step(coverage, layout, factor, options));
    }
    if (!stats.ladder.empty()) {
        stats.rao_fraction_all = stats.ladder.front().fraction_all;
        stats.rao_fraction_nonempty = stats.ladder.front().fraction_nonempty;
    }
    return stats;
}

// --- reading the QC statistics of hicBuildMatrix ----------------------------------------

namespace {

std::int64_t json_int(const json::Value& object, const char* key, std::int64_t fallback = 0) {
    const json::Value* value = object.find(key);
    return value != nullptr && value->is_number()
               ? (value->type() == json::Type::Int ? value->as_int()
                                                   : static_cast<std::int64_t>(value->as_double()))
               : fallback;
}

std::vector<std::int64_t> json_int_array(const json::Value* value) {
    std::vector<std::int64_t> out;
    if (value != nullptr && value->type() == json::Type::Array) {
        for (const json::Value& item : value->as_array()) {
            out.push_back(item.type() == json::Type::Int
                              ? item.as_int()
                              : static_cast<std::int64_t>(item.as_double()));
        }
    }
    return out;
}

std::vector<std::string> json_string_array(const json::Value* value) {
    std::vector<std::string> out;
    if (value != nullptr && value->type() == json::Type::Array) {
        for (const json::Value& item : value->as_array()) {
            out.push_back(item.is_string() ? item.as_string() : std::string());
        }
    }
    return out;
}

}  // namespace

std::optional<ReadQc> parse_read_qc(const json::Value& qc) {
    if (!qc.is_object()) {
        return std::nullopt;
    }
    const json::Value* counts = qc.find("counts");
    if (counts == nullptr || !counts->is_object()) {
        return std::nullopt;
    }
    ReadQc out;
    if (const json::Value* v = qc.find("tool"); v != nullptr && v->is_string()) {
        out.tool = v->as_string();
    }
    if (const json::Value* v = qc.find("tool_version"); v != nullptr && v->is_string()) {
        out.tool_version = v->as_string();
    }
    if (const json::Value* parameters = qc.find("parameters"); parameters != nullptr) {
        if (const json::Value* v = parameters->find("min_rest_site_distance");
            v != nullptr && v->is_number()) {
            out.min_rest_site_distance = json_int(*parameters, "min_rest_site_distance");
        }
        out.max_library_insert_size = json_int(*parameters, "max_library_insert_size");
        out.restriction_sequences = json_string_array(parameters->find("restriction_sequences"));
        out.dangling_sequences = json_string_array(parameters->find("dangling_sequences"));
        out.short_range_cutoff = json_int(*parameters, "short_range_cutoff", 20000);
        if (const json::Value* v = parameters->find("keep_self_ligation"); v != nullptr) {
            out.keep_self_ligation = v->as_bool();
        }
        if (const json::Value* v = parameters->find("restriction_cut_file"); v != nullptr) {
            out.restriction_cut_file = v->as_bool();
        }
    }
    out.sequenced_reads = json_int(*counts, "sequenced_reads");
    out.mappable = json_int(*counts, "pairs_mappable_unique_high_quality");
    out.hic_contacts = json_int(*counts, "hic_contacts");
    out.one_mate_unmapped = json_int(*counts, "one_mate_unmapped");
    out.one_mate_not_unique = json_int(*counts, "one_mate_not_unique");
    out.low_mapping_quality = json_int(*counts, "low_mapping_quality");
    if (const json::Value* list = counts->find("dangling_end");
        list != nullptr && list->type() == json::Type::Array) {
        for (const json::Value& item : list->as_array()) {
            out.dangling_end.push_back(json_int(item, "count"));
        }
    }
    out.self_ligation = json_int(*counts, "self_ligation");
    out.not_close_to_rest_site = json_int(*counts, "one_mate_not_close_to_rest_site");
    out.same_fragment = json_int(*counts, "same_fragment");
    out.self_circle = json_int(*counts, "self_circle");
    out.duplicated_pairs = json_int(*counts, "duplicated_pairs");
    out.inter_chromosomal = json_int(*counts, "inter_chromosomal");
    out.intra_short = json_int(*counts, "intra_short_range");
    out.intra_long = json_int(*counts, "intra_long_range");
    out.orientation_pairs = {json_int(*counts, "inward_pairs"), json_int(*counts, "outward_pairs"),
                             json_int(*counts, "left_pairs"), json_int(*counts, "right_pairs")};
    if (const json::Value* distance = qc.find("distance"); distance != nullptr) {
        out.lower_edges = json_int_array(distance->find("lower_edges_bp"));
        out.distance_contacts = json_int_array(distance->find("contacts"));
        if (const json::Value* orientation = distance->find("orientation"); orientation != nullptr) {
            const char* names[4] = {"inward", "outward", "left", "right"};
            for (std::size_t i = 0; i < 4; ++i) {
                out.orientation_distance[i] = json_int_array(orientation->find(names[i]));
            }
        }
    }
    return out;
}

// --- files ---------------------------------------------------------------------------------------

std::optional<Layout> layout_of(const CoolFile& cool) {
    const json::Value* bin_size = cool.info_value("bin-size");
    if (bin_size == nullptr || !bin_size->is_number() || bin_size->as_double() <= 0) {
        return std::nullopt;
    }
    Layout layout;
    layout.resolution = bin_size->type() == json::Type::Int
                            ? bin_size->as_int()
                            : static_cast<std::int64_t>(bin_size->as_double());
    layout.names = cool.chrom_names();
    layout.lengths = cool.chrom_lengths();
    layout.first_bin.push_back(0);
    for (const std::string& name : layout.names) {
        const auto [first, last] = cool.extent(name);
        if (first != layout.first_bin.back() || last < first) {
            throw std::runtime_error("the bins of '" + name +
                                     "' do not follow the previous chromosome in the bin table");
        }
        layout.first_bin.push_back(last);
    }
    if (layout.first_bin.back() != cool.nbins()) {
        throw std::runtime_error("the chromosomes do not cover the bin table");
    }
    return layout;
}

PixelSource cool_pixel_source(const CoolFile& cool, const Layout& layout) {
    return [&cool, &layout](unsigned threads, const ChunkVisitor& visit) {
        const std::int64_t bins = layout.bins();
        CoolLoadOptions raw;
        raw.apply_correction = false;
        if (threads <= 1 || bins < 4096) {
            for_each_corrected_pixel_chunk(cool, 0, bins, 0, bins, raw,
                                           [&](const PixelChunk& chunk) { visit(0, chunk); });
            return;
        }
        const std::int64_t range =
            std::max<std::int64_t>(1024, bins / (static_cast<std::int64_t>(threads) * 16));
        std::atomic<std::int64_t> next{0};
        std::exception_ptr failure;
        std::mutex failure_lock;
        std::vector<std::thread> workers;
        for (unsigned t = 0; t < threads; ++t) {
            workers.emplace_back([&, t] {
                try {
                    for (;;) {
                        const std::int64_t first = next.fetch_add(range);
                        if (first >= bins) {
                            return;
                        }
                        const std::int64_t last = std::min(bins, first + range);
                        for_each_corrected_pixel_chunk(
                            cool, first, last, 0, bins, raw,
                            [&](const PixelChunk& chunk) { visit(t, chunk); });
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> guard(failure_lock);
                    if (!failure) {
                        failure = std::current_exception();
                    }
                    next.store(bins);
                }
            });
        }
        for (std::thread& worker : workers) {
            worker.join();
        }
        if (failure) {
            std::rethrow_exception(failure);
        }
    };
}

std::vector<std::pair<std::int64_t, std::string>> cool_uris(const std::string& path) {
    std::vector<std::pair<std::int64_t, std::string>> out;
    std::vector<std::string> groups;
    try {
        const h5::File file(path);
        if (file.exists("/pixels") && file.exists("/bins")) {
            groups.push_back("/");
        } else if (file.exists("/resolutions")) {
            for (const std::string& name : file.children("/resolutions")) {
                groups.push_back("/resolutions/" + name);
            }
        } else {
            // the older layout of an mcool: one cooler per top-level group
            for (const std::string& name : file.children("/")) {
                if (file.exists("/" + name + "/pixels") && file.exists("/" + name + "/bins")) {
                    groups.push_back("/" + name);
                }
            }
        }
    } catch (const std::exception& error) {
        throw std::runtime_error("'" + path + "' is not a cool or mcool file (" + error.what() + ")");
    }
    if (groups.empty()) {
        throw std::runtime_error("'" + path + "' holds no cooler");
    }
    for (const std::string& group : groups) {
        if (group == "/") {
            out.emplace_back(0, path);
            continue;
        }
        const std::string uri = path + "::" + group;
        std::int64_t resolution = 0;
        const CoolFile cool(uri);
        if (const json::Value* bin_size = cool.info_value("bin-size");
            bin_size != nullptr && bin_size->is_number()) {
            resolution = static_cast<std::int64_t>(bin_size->as_double());
        }
        out.emplace_back(resolution, uri);
    }
    std::stable_sort(out.begin(), out.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    return out;
}

namespace {

std::string with_resolution(const std::string& label, std::int64_t resolution) {
    return resolution > 0 ? label + "@" + std::to_string(resolution) : label;
}

}  // namespace

FileStats analyse_file(const std::string& path, const std::string& label,
                       const std::vector<std::int64_t>& resolutions, const Options& options) {
    FileStats file;
    file.path = path;
    std::vector<std::pair<std::int64_t, std::string>> uris = cool_uris(path);
    file.multires = !(uris.size() == 1 && uris.front().first == 0);
    if (file.multires && !resolutions.empty()) {
        std::vector<std::pair<std::int64_t, std::string>> chosen;
        for (const std::int64_t wanted : resolutions) {
            const auto found = std::find_if(uris.begin(), uris.end(),
                                            [&](const auto& u) { return u.first == wanted; });
            if (found == uris.end()) {
                std::string present;
                for (const auto& u : uris) {
                    present += (present.empty() ? "" : ", ") + std::to_string(u.first);
                }
                throw std::runtime_error("'" + path + "' has no resolution " +
                                         std::to_string(wanted) + " (present: " + present + ")");
            }
            chosen.push_back(*found);
        }
        uris = chosen;
    }

    std::vector<LadderStep> mcool_steps;
    for (const auto& [resolution, uri] : uris) {
        const CoolFile cool(uri);
        MatrixStats stats;
        const std::optional<Layout> layout = layout_of(cool);
        if (layout.has_value()) {
            std::optional<Correction> correction;
            CoolLoadOptions load;
            load.correction_factor_table = options.weight_column;
            const bool has_column = cool.has_column(options.weight_column);
            if (const std::optional<CoolCorrection> found = cool_correction(cool, load)) {
                correction = Correction{found->factors, found->op};
            }
            stats = compute_matrix_stats(*layout, cool_pixel_source(cool, *layout), correction, options);
            stats.balance.has_weight_column = has_column;
            if (has_column) {
                try {
                    const h5::File h5(cool.filename());
                    const std::string root = cool.root() == "/" ? "" : cool.root();
                    for (const auto& [name, value] :
                         h5.attributes(root + "/bins/" + options.weight_column)) {
                        stats.balance.attributes[name] = attribute_text(value);
                    }
                } catch (const std::exception&) {
                    // a file without attributes on the weights
                }
                if (!correction.has_value() && stats.balance.nan_fraction != stats.balance.nan_fraction) {
                    stats.balance.nan_fraction = 1.0;
                }
            }
        } else {
            stats.pixel_statistics = false;
            stats.note = "The bins are not of a fixed size; the statistics that need the bin size "
                         "are not computed.";
        }
        stats.label = with_resolution(label, file.multires ? resolution : 0);
        stats.uri = uri;
        stats.resolution = layout.has_value() ? layout->resolution : 0;
        stats.bins = cool.nbins();
        stats.pixels = stats.pixel_statistics ? stats.pixels : cool.nnz();
        stats.chromosomes_total = cool.chrom_names().size();
        if (const json::Value* metadata = cool.info_value("metadata");
            metadata != nullptr && metadata->is_object()) {
            if (const json::Value* qc = metadata->find("hicexplorer_qc"); qc != nullptr) {
                stats.read_qc = parse_read_qc(*qc);
            }
        }
        if (file.multires && stats.pixel_statistics && !stats.ladder.empty()) {
            mcool_steps.push_back(stats.ladder.front());
        }
        file.matrices.push_back(std::move(stats));
    }

    if (file.multires) {
        file.map_resolution = pick_map_resolution("resolutions", std::move(mcool_steps));
    } else if (!file.matrices.empty() && file.matrices.front().pixel_statistics) {
        file.map_resolution = pick_map_resolution("aggregation", file.matrices.front().ladder);
    }
    return file;
}

}  // namespace hicx::matrixqc
