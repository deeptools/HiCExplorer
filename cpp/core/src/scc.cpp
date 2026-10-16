#include "hicx/scc.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace hicx::scc {

namespace {

// band[d - 1][row] is the matrix value at (row, row + d), d in [1, nDiags-1],
// row in [0, nbins - d). This mirrors hicrep.utils.trimDiags(a, nDiags,
// bKeepMain=False) followed by a dense re-layout of the upper-triangle band.
using Band = std::vector<std::vector<double>>;

Band build_band(const ChromosomeEntries& e, std::int64_t n_diags) {
    const std::int64_t b = std::max<std::int64_t>(0, n_diags - 1);
    Band band(static_cast<std::size_t>(b));
    for (std::int64_t d = 1; d <= b; ++d) {
        band[static_cast<std::size_t>(d - 1)].assign(
            static_cast<std::size_t>(std::max<std::int64_t>(0, e.nbins - d)), 0.0);
    }
    for (std::size_t i = 0; i < e.rows.size(); ++i) {
        const std::int64_t row = e.rows[i];
        const std::int64_t col = e.cols[i];
        const std::int64_t d = col - row;
        if (d < 1 || d >= n_diags) {
            continue;
        }
        band[static_cast<std::size_t>(d - 1)][static_cast<std::size_t>(row)] += e.values[i];
    }
    return band;
}

// hicrep.utils.meanFilterSparse: a (2h+1) x (2h+1) box filter over the full
// nbins x nbins matrix (zero outside the band, which trimDiags already
// enforced), clipped at the matrix edges, divided by the true number of
// neighbours a cell has (not (2h+1)^2, so edge cells are not spuriously
// down-weighted). Implemented as two separable passes:
//
//  Step 1 (per row, sliding window over the diagonal offset d: this is the
//  horizontal/column direction since column = row + d): widens the band from
//  d in [1, B] to e in [1-h, B+h].
//  Step 2 (per actual column c = row + e, sliding window over row, which
//  moves along the anti-diagonal e' = c - row' of step 1's output): narrows
//  back to d in [1, B], the shape trimDiags produced.
//
// Both passes are O(nbins * (B + h)) via prefix sums; the result is
// mathematically the same 2D box filter scipy's Toeplitz-matrix matmul
// computes, just summed in a different order (so values match to floating
// point, not bit-for-bit).
Band mean_filter_band(const Band& band, std::int64_t nbins, int h) {
    const std::int64_t b = static_cast<std::int64_t>(band.size());
    if (h <= 0 || b == 0) {
        return band;
    }
    const std::int64_t m = b + 2 * h;  // widened e-range length, e in [1-h, b+h]

    // Step 1.
    std::vector<std::vector<double>> i1(static_cast<std::size_t>(nbins));
    for (std::int64_t row = 0; row < nbins; ++row) {
        const std::int64_t maxd = std::min<std::int64_t>(b, nbins - 1 - row);
        std::vector<double>& out_row = i1[static_cast<std::size_t>(row)];
        out_row.assign(static_cast<std::size_t>(m), 0.0);
        if (maxd < 1) {
            continue;
        }
        std::vector<double> prefix(static_cast<std::size_t>(maxd) + 1, 0.0);
        for (std::int64_t d = 1; d <= maxd; ++d) {
            const double v = band[static_cast<std::size_t>(d - 1)][static_cast<std::size_t>(row)];
            prefix[static_cast<std::size_t>(d)] = prefix[static_cast<std::size_t>(d - 1)] + v;
        }
        for (std::int64_t e = 1 - h; e <= b + h; ++e) {
            const std::int64_t lo = std::max<std::int64_t>(1, e - h);
            const std::int64_t hi = std::min<std::int64_t>(maxd, e + h);
            double s = 0.0;
            if (hi >= lo) {
                s = prefix[static_cast<std::size_t>(hi)] - prefix[static_cast<std::size_t>(lo - 1)];
            }
            out_row[static_cast<std::size_t>(e - (1 - h))] = s;
        }
    }

    // Step 2, producing the final band (same shape as the input band).
    Band out(static_cast<std::size_t>(b));
    for (std::int64_t d = 1; d <= b; ++d) {
        out[static_cast<std::size_t>(d - 1)].assign(
            static_cast<std::size_t>(std::max<std::int64_t>(0, nbins - d)), 0.0);
    }
    // c is the final, real column index (row + d), not the widened e-range
    // step 1 used internally: it must stay within the matrix.
    for (std::int64_t c = 1; c <= nbins - 1; ++c) {
        const std::int64_t r_lo = std::max<std::int64_t>(0, c - b - h);
        const std::int64_t r_hi = std::min<std::int64_t>(nbins - 1, c - 1 + h);
        if (r_hi < r_lo) {
            continue;
        }
        const std::int64_t len = r_hi - r_lo + 1;
        std::vector<double> prefix(static_cast<std::size_t>(len) + 1, 0.0);
        for (std::int64_t i = 0; i < len; ++i) {
            const std::int64_t rprime = r_lo + i;
            const std::int64_t eprime = c - rprime;
            double v = 0.0;
            if (eprime >= 1 - h && eprime <= b + h) {
                v = i1[static_cast<std::size_t>(rprime)][static_cast<std::size_t>(eprime - (1 - h))];
            }
            prefix[static_cast<std::size_t>(i + 1)] = prefix[static_cast<std::size_t>(i)] + v;
        }
        const std::int64_t row_lo = std::max<std::int64_t>(0, c - b);
        const std::int64_t row_hi = std::min<std::int64_t>(nbins - 1, c - 1);
        for (std::int64_t row = row_lo; row <= row_hi; ++row) {
            const std::int64_t d = c - row;
            if (d < 1 || d > b) {
                continue;
            }
            const std::int64_t win_lo = std::max<std::int64_t>(row - h, r_lo);
            const std::int64_t win_hi = std::min<std::int64_t>(row + h, r_hi);
            if (win_hi < win_lo) {
                continue;
            }
            const std::int64_t i_lo = win_lo - r_lo;
            const std::int64_t i_hi = win_hi - r_lo;
            const double s =
                prefix[static_cast<std::size_t>(i_hi + 1)] - prefix[static_cast<std::size_t>(i_lo)];
            out[static_cast<std::size_t>(d - 1)][static_cast<std::size_t>(row)] = s;
        }
    }

    // Edge correction: divide by the true number of neighbours, not (2h+1)^2.
    for (std::int64_t d = 1; d <= b; ++d) {
        std::vector<double>& row_vec = out[static_cast<std::size_t>(d - 1)];
        for (std::int64_t row = 0; row < static_cast<std::int64_t>(row_vec.size()); ++row) {
            const std::int64_t col = row + d;
            const std::int64_t row_dist = std::min<std::int64_t>(row, nbins - 1 - row);
            const std::int64_t n_dim1 = h + 1 + std::min<std::int64_t>(row_dist, h);
            const std::int64_t col_dist = std::min<std::int64_t>(col, nbins - 1 - col);
            const std::int64_t n_dim2 = h + 1 + std::min<std::int64_t>(col_dist, h);
            row_vec[static_cast<std::size_t>(row)] /= static_cast<double>(n_dim1 * n_dim2);
        }
    }
    return out;
}

// hicrep.utils.varVstran: the variance of ranked data of size n with Bessel's
// correction. Undefined (treated as zero weight downstream) for n < 2.
double var_vstran(double n) {
    return n < 2.0 ? std::numeric_limits<double>::quiet_NaN() : (1.0 + 1.0 / n) / 12.0;
}

double nan_to_zero(double v) { return std::isfinite(v) ? v : 0.0; }

// hicrep.hicrep.sccByDiag over the (already normalized/smoothed) bands.
double scc_by_diag(const Band& b1, const Band& b2) {
    const std::int64_t n_diags = static_cast<std::int64_t>(b1.size());
    double num = 0.0;
    double den = 0.0;
    for (std::int64_t d = 0; d < n_diags; ++d) {
        const std::vector<double>& r1 = b1[static_cast<std::size_t>(d)];
        const std::vector<double>& r2 = b2[static_cast<std::size_t>(d)];
        double sum1 = 0.0, sum2 = 0.0, sumsq1 = 0.0, sumsq2 = 0.0, sum12 = 0.0;
        double n_samples = 0.0;
        for (std::size_t i = 0; i < r1.size(); ++i) {
            const double a = r1[i];
            const double c = r2[i];
            if (a == 0.0 && c == 0.0) {
                continue;  // union of nonzero positions only (eliminate_zeros)
            }
            sum1 += a;
            sum2 += c;
            sumsq1 += a * a;
            sumsq2 += c * c;
            sum12 += a * c;
            n_samples += 1.0;
        }
        const double cov = sum12 - sum1 * sum2 / n_samples;
        const double var1 = sumsq1 - sum1 * sum1 / n_samples;
        const double var2 = sumsq2 - sum2 * sum2 / n_samples;
        double rho = cov / std::sqrt(var1 * var2);
        double ws = n_samples * var_vstran(n_samples);
        rho = nan_to_zero(rho);
        ws = nan_to_zero(ws);
        num += rho * ws;
        den += ws;
    }
    if (den == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return num / den;
}

// hicrep.utils.resample: multinomial downsampling of a band's flattened
// nonzero entries to `target` total count. Not numpy-equivalent (see
// scc.hpp), but the same procedure: draw `target` samples with replacement,
// weighted by each entry's share of the total.
void resample_band(Band& band, double target, std::uint64_t seed) {
    std::vector<double> weights;
    std::vector<std::pair<std::size_t, std::size_t>> positions;  // (diag, row)
    for (std::size_t d = 0; d < band.size(); ++d) {
        for (std::size_t row = 0; row < band[d].size(); ++row) {
            if (band[d][row] != 0.0) {
                weights.push_back(band[d][row]);
                positions.emplace_back(d, row);
            }
        }
    }
    if (weights.empty()) {
        return;
    }
    std::mt19937_64 rng(seed);
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    std::vector<double> counts(weights.size(), 0.0);
    const auto draws = static_cast<std::int64_t>(std::llround(target));
    for (std::int64_t i = 0; i < draws; ++i) {
        counts[dist(rng)] += 1.0;
    }
    for (std::size_t d = 0; d < band.size(); ++d) {
        std::fill(band[d].begin(), band[d].end(), 0.0);
    }
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto [d, row] = positions[i];
        band[d][row] = counts[i];
    }
}

double band_sum(const Band& band) {
    double s = 0.0;
    for (const std::vector<double>& row : band) {
        for (const double v : row) {
            s += v;
        }
    }
    return s;
}

}  // namespace

double chromosome_scc(const ChromosomeEntries& e1, const ChromosomeEntries& e2,
                      double total_contacts_1, double total_contacts_2, std::int64_t bin_size,
                      const Options& options) {
    if (e1.nbins != e2.nbins || e1.nbins <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::int64_t d_max = e1.nbins;
    if (options.max_distance_bp >= 0 && bin_size > 0) {
        d_max = options.max_distance_bp / bin_size + 1;
    }
    const std::int64_t n_diags = std::min<std::int64_t>(d_max, e1.nbins);
    if (n_diags <= 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    Band band1 = build_band(e1, n_diags);
    Band band2 = build_band(e2, n_diags);

    if (options.downsample) {
        const double size1 = band_sum(band1);
        const double size2 = band_sum(band2);
        if (size1 > size2) {
            resample_band(band1, size2, options.downsample_seed);
        } else if (size2 > size1) {
            resample_band(band2, size1, options.downsample_seed);
        }
    } else {
        for (std::vector<double>& row : band1) {
            for (double& v : row) v /= total_contacts_1;
        }
        for (std::vector<double>& row : band2) {
            for (double& v : row) v /= total_contacts_2;
        }
    }

    if (options.h > 0) {
        band1 = mean_filter_band(band1, e1.nbins, options.h);
        band2 = mean_filter_band(band2, e1.nbins, options.h);
    }

    return scc_by_diag(band1, band2);
}

}  // namespace hicx::scc
