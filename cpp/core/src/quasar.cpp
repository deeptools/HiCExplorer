#include "hicx/quasar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace hicx::quasar {

namespace {

struct Raw {
    std::int64_t bin1 = 0;
    std::int64_t bin2 = 0;
    std::int64_t count = 0;
};

std::vector<Raw> to_raw(const ChromosomeEntries& e) {
    std::vector<Raw> raw;
    raw.reserve(e.rows.size());
    for (std::size_t i = 0; i < e.rows.size(); ++i) {
        raw.push_back({e.rows[i], e.cols[i], static_cast<std::int64_t>(std::llround(e.values[i]))});
    }
    return raw;
}

// hifive.libraries._quasar.filter_bins: counts, per bin, the number of
// other valid bins it shares a stored (nonzero) contact with within
// `width`, diagonal excluded.
void filter_bins(const std::vector<Raw>& raw, std::vector<int>& sums, const std::vector<char>& valid,
                 std::int64_t width) {
    std::fill(sums.begin(), sums.end(), 0);
    for (const Raw& r : raw) {
        const std::int64_t span = r.bin2 - r.bin1;
        if (!valid[static_cast<std::size_t>(r.bin1)] || !valid[static_cast<std::size_t>(r.bin2)] ||
            span > width || span == 0) {
            continue;
        }
        ++sums[static_cast<std::size_t>(r.bin1)];
        ++sums[static_cast<std::size_t>(r.bin2)];
    }
}

struct NormResult {
    std::int64_t n = 0;
    std::int64_t width = 0;
    std::vector<double> norm;       // n x (2*width), row-major
    std::vector<double> dist;       // n x dist_width, row-major
    std::int64_t dist_width = 0;
    std::vector<char> valid_rows;   // n
};

// Quasar._normalize, minus the downbinning step (the caller's entries are
// already at the matrix's own resolution).
NormResult normalize(const ChromosomeEntries& e, std::int64_t width) {
    NormResult out;
    const std::int64_t n = e.nbins;
    out.n = n;
    out.width = width;
    out.dist_width = std::min<std::int64_t>(width, std::max<std::int64_t>(n - 1, 0));

    const std::vector<Raw> raw = to_raw(e);
    std::vector<char> valid(static_cast<std::size_t>(n), 1);
    std::vector<int> sums(static_cast<std::size_t>(n), 0);
    std::int64_t prev_valid = n + 1;
    std::int64_t current_valid = 0;
    for (char v : valid) current_valid += v;
    while (current_valid < prev_valid) {
        prev_valid = current_valid;
        filter_bins(raw, sums, valid, width);
        current_valid = 0;
        for (std::int64_t i = 0; i < n; ++i) {
            valid[static_cast<std::size_t>(i)] = sums[static_cast<std::size_t>(i)] >= 1 ? 1 : 0;
            current_valid += valid[static_cast<std::size_t>(i)];
        }
    }
    out.valid_rows = valid;
    if (current_valid == 0) {
        return out;
    }

    // find_bg_dist_norm.
    std::vector<double> bg(static_cast<std::size_t>(width), 1.0);
    std::vector<std::int64_t> temp_count(static_cast<std::size_t>(width), 0);  // temp[:,1]
    std::vector<std::int64_t> temp_sum(static_cast<std::size_t>(width), 0);    // temp[:,0]
    out.dist.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(out.dist_width), 0.0);
    out.norm.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(2 * width), 0.0);

    auto dist_at = [&](std::int64_t i, std::int64_t span) -> double& {
        return out.dist[static_cast<std::size_t>(i) * static_cast<std::size_t>(out.dist_width) +
                        static_cast<std::size_t>(span)];
    };
    auto norm_at = [&](std::int64_t i, std::int64_t col) -> double& {
        return out.norm[static_cast<std::size_t>(i) * static_cast<std::size_t>(2 * width) +
                        static_cast<std::size_t>(col)];
    };

    for (std::int64_t i = 0; i < n; ++i) {
        if (!valid[static_cast<std::size_t>(i)]) continue;
        for (std::int64_t j = i + 1; j < std::min(n, i + width + 1); ++j) {
            if (!valid[static_cast<std::size_t>(j)]) continue;
            const std::int64_t span = j - i - 1;
            ++temp_count[static_cast<std::size_t>(span)];
            if (span < out.dist_width) dist_at(i, span) = 1.0;
        }
    }
    for (const Raw& r : raw) {
        const std::int64_t span = r.bin2 - r.bin1 - 1;
        if (!valid[static_cast<std::size_t>(r.bin1)] || !valid[static_cast<std::size_t>(r.bin2)] ||
            span >= width || span < 0) {
            continue;
        }
        temp_sum[static_cast<std::size_t>(span)] += r.count;
        if (span < out.dist_width) dist_at(r.bin1, span) = std::sqrt(1.0 + static_cast<double>(r.count));
        norm_at(r.bin1, width + span) = static_cast<double>(r.count);
        norm_at(r.bin2, width - span - 1) = static_cast<double>(r.count);
    }
    for (std::int64_t i = 0; i < width; ++i) {
        bg[static_cast<std::size_t>(i)] = (temp_count[static_cast<std::size_t>(i)] > 0 &&
                                           temp_sum[static_cast<std::size_t>(i)] > 0)
                                              ? static_cast<double>(temp_sum[static_cast<std::size_t>(i)]) /
                                                    static_cast<double>(temp_count[static_cast<std::size_t>(i)])
                                              : 1.0;
    }
    for (std::int64_t i = 0; i < n; ++i) {
        for (std::int64_t col = 0; col < width; ++col) {
            norm_at(i, col) /= bg[static_cast<std::size_t>(width - 1 - col)];
        }
        for (std::int64_t col = width; col < 2 * width; ++col) {
            norm_at(i, col) /= bg[static_cast<std::size_t>(col - width)];
        }
    }
    return out;
}

// hifive.libraries._quasar.find_correlations.
std::vector<double> find_correlations(const NormResult& norm) {
    const std::int64_t n = norm.n;
    const std::int64_t width = norm.width;
    std::vector<double> corr(static_cast<std::size_t>(n) * static_cast<std::size_t>(width),
                             std::numeric_limits<double>::infinity());
    auto norms_at = [&](std::int64_t row, std::int64_t pos) {
        return norm.norm[static_cast<std::size_t>(row) * static_cast<std::size_t>(2 * width) +
                         static_cast<std::size_t>(pos)];
    };
    auto corr_at = [&](std::int64_t row, std::int64_t col) -> double& {
        return corr[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                    static_cast<std::size_t>(col)];
    };
    const std::vector<char>& vrows = norm.valid_rows;

    for (std::int64_t x = 0; x < n; ++x) {
        if (!vrows[static_cast<std::size_t>(x)]) continue;
        for (std::int64_t y = x + 1; y < std::min(n, x + width + 1); ++y) {
            if (!vrows[static_cast<std::size_t>(y)]) continue;
            double x1 = 0, y1 = 0, x2 = 0, y2 = 0, xy = 0;
            std::int64_t count = 0;
            for (std::int64_t k = std::max<std::int64_t>(0, y - width); k < x; ++k) {
                if (!vrows[static_cast<std::size_t>(k)]) continue;
                const double tx = norms_at(x, width - (x - k));
                const double ty = norms_at(y, width - (y - k));
                x1 += tx;
                x2 += tx * tx;
                y1 += ty;
                y2 += ty * ty;
                xy += tx * ty;
                ++count;
            }
            for (std::int64_t k = x + 1; k < y; ++k) {
                if (!vrows[static_cast<std::size_t>(k)]) continue;
                const double tx = norms_at(x, width + k - x - 1);
                const double ty = norms_at(y, width - (y - k));
                x1 += tx;
                x2 += tx * tx;
                y1 += ty;
                y2 += ty * ty;
                xy += tx * ty;
                ++count;
            }
            for (std::int64_t k = y + 1; k < std::min(n, x + width + 1); ++k) {
                if (!vrows[static_cast<std::size_t>(k)]) continue;
                const double tx = norms_at(x, width + k - x - 1);
                const double ty = norms_at(y, width + k - y - 1);
                x1 += tx;
                x2 += tx * tx;
                y1 += ty;
                y2 += ty * ty;
                xy += tx * ty;
                ++count;
            }
            if (count < 3) continue;
            const double n_d = static_cast<double>(count);
            x1 /= n_d;
            x2 = x2 / n_d - x1 * x1;
            y1 /= n_d;
            y2 = y2 / n_d - y1 * y1;
            if (x2 <= 0.0 || y2 <= 0.0) continue;
            xy /= n_d;
            corr_at(x, y - x - 1) = (xy - x1 * y1) / std::sqrt(x2 * y2);
        }
    }
    return corr;
}

}  // namespace

double chromosome_score(const ChromosomeEntries& e1, const ChromosomeEntries& e2, const Options& options) {
    if (e1.nbins != e2.nbins || e1.nbins <= 0 || e1.rows.empty() || e2.rows.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const NormResult norm1 = normalize(e1, options.width);
    const NormResult norm2 = normalize(e2, options.width);
    if (norm1.norm.empty() || norm2.norm.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::vector<double> corr1 = find_correlations(norm1);
    const std::vector<double> corr2 = find_correlations(norm2);

    const std::int64_t n = norm1.n;
    const std::int64_t m1 = norm1.dist_width;
    const std::int64_t m2 = norm2.dist_width;
    const std::int64_t stop3 = std::min(m1, m2);

    auto dist1_at = [&](std::int64_t i, std::int64_t l) {
        return norm1.dist[static_cast<std::size_t>(i) * static_cast<std::size_t>(m1) + static_cast<std::size_t>(l)];
    };
    auto dist2_at = [&](std::int64_t i, std::int64_t l) {
        return norm2.dist[static_cast<std::size_t>(i) * static_cast<std::size_t>(m2) + static_cast<std::size_t>(l)];
    };
    auto corr1_at = [&](std::int64_t i, std::int64_t l) {
        return corr1[static_cast<std::size_t>(i) * static_cast<std::size_t>(norm1.width) +
                     static_cast<std::size_t>(l)];
    };
    auto corr2_at = [&](std::int64_t i, std::int64_t l) {
        return corr2[static_cast<std::size_t>(i) * static_cast<std::size_t>(norm2.width) +
                     static_cast<std::size_t>(l)];
    };

    // find_replicate_scores's own weight is (1 + dist)**0.5, where `dist`
    // here is already sqrt(1 + count) from find_bg_dist_norm: a deliberate
    // double square root, not a typo, matching the real source exactly.
    double x = 0, y = 0, x2 = 0, y2 = 0, xy = 0;
    std::int64_t n_valid = 0;
    for (std::int64_t l = 0; l < stop3; ++l) {
        const std::int64_t p = n - l - 1;
        for (std::int64_t row = 0; row < p; ++row) {
            const bool valid1 = norm1.valid_rows[static_cast<std::size_t>(row + l + 1)] &&
                                norm1.valid_rows[static_cast<std::size_t>(row)] &&
                                std::isfinite(dist1_at(row, l)) && std::isfinite(corr1_at(row, l));
            const bool valid2 = norm2.valid_rows[static_cast<std::size_t>(row + l + 1)] &&
                                norm2.valid_rows[static_cast<std::size_t>(row)] &&
                                std::isfinite(dist2_at(row, l)) && std::isfinite(corr2_at(row, l));
            if (!valid1 || !valid2) continue;
            const double d1 = std::sqrt(1.0 + dist1_at(row, l));
            const double d2 = std::sqrt(1.0 + dist2_at(row, l));
            const double t1 = corr1_at(row, l) * d1;
            const double t2 = corr2_at(row, l) * d2;
            x += t1;
            y += t2;
            x2 += t1 * t1;
            y2 += t2 * t2;
            xy += t1 * t2;
            ++n_valid;
        }
    }
    if (n_valid == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double nd = static_cast<double>(n_valid);
    const double xmu = x / nd;
    const double ymu = y / nd;
    const double x2mu = x2 / nd;
    const double y2mu = y2 / nd;
    const double xymu = xy / nd;
    if (xmu * xmu > x2mu || ymu * ymu > y2mu) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double xstd = std::sqrt(x2mu - xmu * xmu);
    const double ystd = std::sqrt(y2mu - ymu * ymu);
    if (xstd == 0.0 || ystd == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (xymu - xmu * ymu) / (xstd * ystd);
}

}  // namespace hicx::quasar
