#include "stripes_impl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#include "hicx/parallel.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::stripes {

namespace {

constexpr double kEps = 1e-12;

}  // namespace

Band build_horizontal_band(const CsrMatrix& upper_triangle, std::int64_t n_bins,
                           std::int64_t max_distance, std::vector<double>* expected_out) {
    Band band;
    band.n_bins = n_bins;
    band.max_distance = max_distance;
    band.raw.assign(static_cast<std::size_t>(n_bins * max_distance), 0.0);

    const std::vector<std::int64_t>& indptr = upper_triangle.indptr();
    const std::vector<std::int32_t>& indices = upper_triangle.indices();
    const std::vector<double>& values = upper_triangle.data();
    for (std::int64_t row = 0; row < upper_triangle.rows() && row < n_bins; ++row) {
        const std::size_t begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            const std::int64_t col = static_cast<std::int64_t>(indices[k]);
            const std::int64_t distance = col - row;
            if (distance < 1 || distance > max_distance) {
                continue;
            }
            band.raw[static_cast<std::size_t>(row * max_distance + (distance - 1))] = values[k];
        }
    }

    std::vector<double> expected(static_cast<std::size_t>(max_distance), 0.0);
    for (std::int64_t d = 1; d <= max_distance; ++d) {
        const std::int64_t positions = n_bins - d;
        if (positions <= 0) {
            continue;
        }
        double sum = 0.0;
        for (std::int64_t a = 0; a < positions; ++a) {
            sum += band.raw[static_cast<std::size_t>(a * max_distance + (d - 1))];
        }
        expected[static_cast<std::size_t>(d - 1)] = sum / static_cast<double>(positions);
    }

    band.obs_exp.assign(static_cast<std::size_t>(n_bins * max_distance), 0.0);
    for (std::int64_t a = 0; a < n_bins; ++a) {
        for (std::int64_t d = 1; d <= max_distance && a + d < n_bins + 1; ++d) {
            const double e = expected[static_cast<std::size_t>(d - 1)];
            if (e <= kEps) {
                continue;
            }
            const std::size_t idx = static_cast<std::size_t>(a * max_distance + (d - 1));
            band.obs_exp[idx] = band.raw[idx] / e;
        }
    }

    if (expected_out != nullptr) {
        *expected_out = std::move(expected);
    }
    return band;
}

Band build_vertical_band(const Band& horizontal) {
    Band band;
    band.n_bins = horizontal.n_bins;
    band.max_distance = horizontal.max_distance;
    band.raw.assign(horizontal.raw.size(), 0.0);
    band.obs_exp.assign(horizontal.obs_exp.size(), 0.0);
    for (std::int64_t c = 0; c < band.n_bins; ++c) {
        for (std::int64_t d = 1; d <= band.max_distance; ++d) {
            const std::int64_t a = c - d;
            if (a < 0) {
                continue;
            }
            const std::size_t dst = static_cast<std::size_t>(c * band.max_distance + (d - 1));
            band.raw[dst] = horizontal.raw_at(a, d);
            band.obs_exp[dst] = horizontal.obs_exp_at(a, d);
        }
    }
    return band;
}

RunningSums build_running_sums(const Band& band) {
    RunningSums sums;
    sums.n_bins = band.n_bins;
    sums.max_distance = band.max_distance;
    sums.obs_exp_cum.assign(band.obs_exp.size(), 0.0);
    sums.raw_cum.assign(band.raw.size(), 0.0);
    for (std::int64_t a = 0; a < band.n_bins; ++a) {
        double oe_cum = 0.0;
        double raw_cum = 0.0;
        for (std::int64_t d = 1; d <= band.max_distance; ++d) {
            const std::size_t idx = static_cast<std::size_t>(a * band.max_distance + (d - 1));
            oe_cum += band.obs_exp[idx];
            raw_cum += band.raw[idx];
            sums.obs_exp_cum[idx] = oe_cum;
            sums.raw_cum[idx] = raw_cum;
        }
    }
    return sums;
}

namespace {

// The background anchors for a candidate at `anchor`, length `length`,
// orientation `vertical`: anchors on each side of `anchor` separated by
// `background_gap_bins`, out to `background_window_bins` further, clipped to
// [0, n_bins) and to anchors that have room for the full length (a + length
// <= n_bins for horizontal, a - length >= 0 for vertical).
std::vector<std::int64_t> background_anchors(std::int64_t anchor, std::int64_t length,
                                              bool vertical, std::int64_t n_bins,
                                              const DetectOptions& options) {
    std::vector<std::int64_t> result;
    result.reserve(static_cast<std::size_t>(options.background_window_bins) * 2);
    const auto has_room = [&](std::int64_t a) {
        return vertical ? (a - length >= 0) : (a + length <= n_bins);
    };
    for (std::int64_t offset = options.background_gap_bins + 1;
        offset <= options.background_gap_bins + options.background_window_bins; ++offset) {
        const std::int64_t lo = anchor - offset;
        const std::int64_t hi = anchor + offset;
        if (lo >= 0 && has_room(lo)) {
            result.push_back(lo);
        }
        if (hi < n_bins && has_room(hi)) {
            result.push_back(hi);
        }
    }
    return result;
}

struct MeanStd {
    double mean = 0.0;
    double std = 0.0;
    std::size_t count = 0;
};

MeanStd mean_std(const std::vector<double>& values) {
    MeanStd result;
    result.count = values.size();
    if (values.empty()) {
        return result;
    }
    double sum = 0.0;
    for (const double v : values) {
        sum += v;
    }
    result.mean = sum / static_cast<double>(values.size());
    if (values.size() < 2) {
        return result;
    }
    double sq = 0.0;
    for (const double v : values) {
        const double d = v - result.mean;
        sq += d * d;
    }
    result.std = std::sqrt(sq / static_cast<double>(values.size()));
    return result;
}

// One orientation's preselection scan, threaded over anchors. Each anchor
// writes only to its own slot.
void scan_orientation(const Band& band, const RunningSums& sums, bool vertical,
                      const DetectOptions& options, unsigned int threads,
                      std::vector<std::optional<Candidate>>& out) {
    const std::int64_t n_bins = band.n_bins;
    out.assign(static_cast<std::size_t>(n_bins), std::nullopt);
    parallel_for(static_cast<std::size_t>(n_bins), threads, [&](std::size_t index) {
        const std::int64_t anchor = static_cast<std::int64_t>(index);
        std::optional<Candidate> best;
        for (const std::int64_t length : options.length_grid_bins) {
            const bool has_room = vertical ? (anchor - length >= 0) : (anchor + length <= n_bins);
            if (!has_room) {
                continue;
            }
            const double raw_mean = sums.raw_mean(anchor, length);
            if (raw_mean < options.min_raw_count) {
                continue;
            }
            const double stripe_mean = sums.obs_exp_mean(anchor, length);
            const std::vector<std::int64_t> neighbours =
                background_anchors(anchor, length, vertical, n_bins, options);
            if (neighbours.size() < 4) {
                continue;
            }
            std::vector<double> bg_means;
            bg_means.reserve(neighbours.size());
            for (const std::int64_t other : neighbours) {
                bg_means.push_back(sums.obs_exp_mean(other, length));
            }
            const MeanStd stats = mean_std(bg_means);
            if (stats.mean <= kEps) {
                continue;
            }
            const double enrichment = stripe_mean / stats.mean;
            if (enrichment < options.min_obs_exp) {
                continue;
            }
            const double denom = stats.std > kEps ? stats.std : kEps;
            const double z = (stripe_mean - stats.mean) / denom;
            if (z < options.preselect_z) {
                continue;
            }
            if (!best.has_value() || z > best->zscore) {
                Candidate candidate;
                candidate.anchor = anchor;
                candidate.length_bins = length;
                candidate.vertical = vertical;
                candidate.enrichment = enrichment;
                candidate.zscore = z;
                best = candidate;
            }
        }
        out[index] = best;
    });
}

}  // namespace

std::vector<Candidate> preselect(const Band& horizontal, const Band& vertical,
                                 const RunningSums& horizontal_sums,
                                 const RunningSums& vertical_sums, const DetectOptions& options,
                                 unsigned int threads) {
    std::vector<std::optional<Candidate>> horizontal_best;
    std::vector<std::optional<Candidate>> vertical_best;
    scan_orientation(horizontal, horizontal_sums, false, options, threads, horizontal_best);
    scan_orientation(vertical, vertical_sums, true, options, threads, vertical_best);

    std::vector<Candidate> result;
    result.reserve(horizontal_best.size() + vertical_best.size());
    for (const auto& entry : horizontal_best) {
        if (entry.has_value()) {
            result.push_back(*entry);
        }
    }
    for (const auto& entry : vertical_best) {
        if (entry.has_value()) {
            result.push_back(*entry);
        }
    }
    return result;
}

std::vector<Candidate> suppress_non_maximal(std::vector<Candidate> candidates,
                                            std::int64_t merge_window_bins) {
    // Stable order: strongest z first, ties broken by anchor ascending, so
    // the result does not depend on the input order (which parallel_for does
    // not guarantee across a resumed run, though this port always emits
    // anchor order; this keeps the guarantee explicit).
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const Candidate& a, const Candidate& b) {
                         if (a.vertical != b.vertical) {
                             return a.vertical < b.vertical;
                         }
                         if (a.zscore != b.zscore) {
                             return a.zscore > b.zscore;
                         }
                         return a.anchor < b.anchor;
                     });

    std::vector<Candidate> kept;
    kept.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        bool suppressed = false;
        for (const Candidate& already : kept) {
            if (already.vertical != candidate.vertical) {
                continue;
            }
            if (std::llabs(static_cast<long long>(already.anchor - candidate.anchor)) <=
                merge_window_bins) {
                suppressed = true;
                break;
            }
        }
        if (!suppressed) {
            kept.push_back(candidate);
        }
    }
    // Restore a deterministic, orientation-then-anchor order for reporting.
    std::stable_sort(kept.begin(), kept.end(), [](const Candidate& a, const Candidate& b) {
        if (a.vertical != b.vertical) {
            return a.vertical < b.vertical;
        }
        return a.anchor < b.anchor;
    });
    return kept;
}

void compute_pvalues(std::vector<Candidate>& candidates, const Band& horizontal,
                     const Band& vertical, const DetectOptions& options, unsigned int threads) {
    parallel_for(candidates.size(), threads, [&](std::size_t index) {
        Candidate& candidate = candidates[index];
        const Band& band = candidate.vertical ? vertical : horizontal;
        const std::int64_t length = candidate.length_bins;
        std::vector<double> stripe_values;
        stripe_values.reserve(static_cast<std::size_t>(length));
        for (std::int64_t d = 1; d <= length; ++d) {
            stripe_values.push_back(band.obs_exp_at(candidate.anchor, d));
        }
        const std::vector<std::int64_t> neighbours = background_anchors(
            candidate.anchor, length, candidate.vertical, band.n_bins, options);
        std::vector<double> background_values;
        background_values.reserve(neighbours.size() * static_cast<std::size_t>(length));
        for (const std::int64_t other : neighbours) {
            for (std::int64_t d = 1; d <= length; ++d) {
                background_values.push_back(band.obs_exp_at(other, d));
            }
        }
        if (stripe_values.empty() || background_values.empty()) {
            candidate.pvalue = 1.0;
            return;
        }
        const stats::RanksumsResult result = stats::ranksums(stripe_values, background_values);
        candidate.pvalue = result.defined && std::isfinite(result.pvalue) ? result.pvalue : 1.0;
    });
}

std::vector<Candidate> apply_fdr(std::vector<Candidate> candidates, double fdr_q) {
    std::vector<double> pvalues;
    pvalues.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        pvalues.push_back(candidate.pvalue);
    }
    const std::vector<double> adjusted = stats::benjamini_hochberg_adjusted(pvalues);
    std::vector<Candidate> kept;
    kept.reserve(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        candidates[i].qvalue = adjusted[i];
        if (adjusted[i] <= fdr_q) {
            kept.push_back(candidates[i]);
        }
    }
    return kept;
}

}  // namespace hicx::stripes
