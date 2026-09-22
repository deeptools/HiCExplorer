// See hicx/chicago.hpp.

#include "hicx/chicago.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

#include "hicx/bins.hpp"
#include "hicx/tool_matrix.hpp"

namespace hicx::chicago {

namespace {

constexpr double kDblEps = std::numeric_limits<double>::epsilon();

// log P(X > N-1) for NB(size = r, mu) via the incomplete beta identity used
// throughout this codebase (hicx::scipy::betainc mirrors it for scipy; this
// is the same identity applied directly through Boost, since R's pbeta and
// scipy's betainc both reduce to the same regularized incomplete beta and
// agree far inside the ED tolerance).
double log_nbinom_sf(double n_minus_1, double size, double mu) {
    if (size <= 0.0 || mu < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double N = n_minus_1 + 1.0;  // the smallest count with X >= N
    if (N <= 0.0) {
        return 0.0;  // P(X >= 0) = 1
    }
    // P(X > N-1) = P(X >= N) = ibeta(N, size, 1 - prob) = ibeta(N, size, mu / (size + mu)),
    // prob = size / (size + mu) being the NB "prob" parameter.
    const double x = mu / (size + mu);
    try {
        const double sf = boost::math::ibeta(N, size, x);
        return std::log(sf);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double log_poisson_sf(double n_minus_1, double lambda) {
    const double N = n_minus_1 + 1.0;
    if (N <= 0.0) {
        return 0.0;
    }
    if (lambda <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    try {
        const double sf = boost::math::gamma_p(N, lambda);  // P(X >= N)
        return std::log(sf);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// log NB pmf, k = 0..(size-1), built by the standard multiplicative
// recurrence p(k) = p(k-1) * (k-1+alpha)/k * beta/(1+beta), which is stable
// (no repeated lgamma calls, no cancellation) and lets the tail be extended
// one term at a time.
class NbPmfSeries {
   public:
    NbPmfSeries(double alpha, double beta) : alpha_(alpha), beta_(beta) {
        log_terms_.push_back(-alpha_ * std::log1p(beta_));
        log_ratio_const_ = std::log(beta_) - std::log1p(beta_);
    }
    double log_pmf(long k) {
        while (static_cast<long>(log_terms_.size()) <= k) {
            const long j = static_cast<long>(log_terms_.size());
            log_terms_.push_back(log_terms_.back() + std::log(static_cast<double>(j - 1) + alpha_) -
                                  std::log(static_cast<double>(j)) + log_ratio_const_);
        }
        return log_terms_[k];
    }

   private:
    double alpha_, beta_;
    double log_ratio_const_;
    std::vector<double> log_terms_;
};

// log Poisson pmf, j = 0.., recurrence p(j) = p(j-1) * lambda / j.
class PoissonPmfSeries {
   public:
    explicit PoissonPmfSeries(double lambda) : lambda_(lambda) {
        log_terms_.push_back(lambda_ <= 0.0 ? 0.0 : -lambda_);
    }
    double log_pmf(long j) {
        while (static_cast<long>(log_terms_.size()) <= j) {
            if (lambda_ <= 0.0) {
                log_terms_.push_back(-std::numeric_limits<double>::infinity());
                continue;
            }
            const long k = static_cast<long>(log_terms_.size());
            log_terms_.push_back(log_terms_.back() + std::log(lambda_) - std::log(static_cast<double>(k)));
        }
        return log_terms_[j];
    }

   private:
    double lambda_;
    std::vector<double> log_terms_;
};

// Delaporte survival P(X >= N) for X = NB(alpha, beta) + Poisson(lambda), N
// an integer count (n_minus_1 = N - 1, passed the way R's pdelap(N - 1, ...)
// is called; always an exact integer in this codebase, chinput counts are
// integers).
//
// Computed by summing the smaller of the two tails directly (never as 1 -
// the larger tail), because the larger tail is typically within a few ulps
// of 1 for the deep-tail p-values CHiCAGO calls on, and 1 - (that) loses
// every significant digit: measured against R's own Delaporte::pdelap, a
// naive "always sum the left tail, subtract from 1" implementation is right
// to a handful of digits for p-values around 1e-10 and pure floating-point
// noise (wrong by ~130 orders of magnitude) for the very small p-values
// CHiCAGO calls actually depend on, e.g. a real GM12878 chr20 bait/other-end
// pair with log.p = -348.5.
double log_delaporte_sf(double n_minus_1, double alpha, double beta, double lambda) {
    const long N = static_cast<long>(std::llround(n_minus_1)) + 1;
    if (N <= 0) {
        return 0.0;  // P(X >= 0) = 1
    }
    const double mean = alpha * beta + lambda;

    // Truncate the Poisson component: since lambda is CHiCAGO's technical
    // noise mean, typically well under a few hundred, its pmf is negligible
    // beyond a handful of standard deviations past the mean, and every
    // P(X = n) below sums nb_pmf(n - j) * pois_pmf(j) over that truncated j.
    PoissonPmfSeries pois(lambda);
    long j_max = static_cast<long>(std::ceil(lambda + 40.0 * std::sqrt(std::max(lambda, 1.0)) + 20.0));
    j_max = std::min(j_max, 20000L);
    while (j_max > 0 && !std::isfinite(pois.log_pmf(j_max))) {
        --j_max;  // lambda == 0: only j = 0 has mass
    }

    NbPmfSeries nb(alpha, beta);

    auto term_at = [&](long n) {
        double sum = 0.0;
        for (long j = 0; j <= std::min(j_max, n); ++j) {
            const double lp = nb.log_pmf(n - j) + pois.log_pmf(j);
            if (std::isfinite(lp)) {
                sum += std::exp(lp);
            }
        }
        return sum;
    };

    if (static_cast<double>(N - 1) <= mean) {
        // Left tail (CDF at N - 1) is the smaller one: sum it directly.
        double cdf = 0.0;
        for (long n = 0; n <= N - 1; ++n) {
            cdf += term_at(n);
        }
        const double sf = 1.0 - cdf;
        if (!(sf > 0.0)) {
            return -std::numeric_limits<double>::infinity();
        }
        return std::log(sf);
    }

    // Right tail (survival at N) is the smaller one: sum it directly,
    // extending until additional terms stop contributing.
    double sf = 0.0;
    long consecutive_negligible = 0;
    const long kHardCap = 2000000;
    for (long n = N; n < N + kHardCap; ++n) {
        const double term = term_at(n);
        sf += term;
        if (sf > 0.0 && term < sf * 1e-17) {
            if (++consecutive_negligible > 50) break;
        } else {
            consecutive_negligible = 0;
        }
        if (term == 0.0 && sf == 0.0 && n > N + 10000) {
            break;  // underflowed to nothing: genuinely (numerically) zero
        }
    }
    if (!(sf > 0.0)) {
        return -std::numeric_limits<double>::infinity();
    }
    return std::log(sf);
}

double expit(double x) { return 1.0 / (1.0 + std::exp(-x)); }

}  // namespace

double log_pvalue(double N, double alpha, double Bmean, double Tmean) {
    if (Bmean < kDblEps) {
        return log_poisson_sf(N - 1.0, Tmean);
    }
    const double beta = Bmean / alpha;
    double lp = log_delaporte_sf(N - 1.0, alpha, beta, Tmean);
    if (!std::isfinite(lp)) {
        double gamma = alpha * (1.0 + Tmean / Bmean) * (1.0 + Tmean / Bmean);
        gamma = std::min(gamma, 1e10);
        lp = log_nbinom_sf(N - 1.0, gamma, Bmean + Tmean);
    }
    return lp;
}

double log_weight(double abs_dist, const WeightSettings& w, double eta_bar) {
    const double dist = (abs_dist < 0.0) ? std::numeric_limits<double>::infinity() : abs_dist;
    const double eta = expit(w.alpha + w.beta * std::log(dist));
    const double d = expit(w.delta), g = expit(w.gamma);
    const double log_w = std::log((d - g) * eta + g) - std::log((d - g) * eta_bar + g);
    return log_w;
}

double score_from_pvalue(double log_p, double abs_dist, const WeightSettings& w, double eta_bar) {
    const double log_w = log_weight(abs_dist, w, eta_bar);
    const double log_q = log_p - log_w;
    const double minval = log_weight(0.0, w, eta_bar);
    const double score = -minval - log_q;
    return std::max(score, 0.0);
}

namespace {

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, '\t')) {
        out.push_back(tok);
    }
    return out;
}

}  // namespace

std::vector<RmapFragment> read_rmap(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open rmap file: " + path);
    }
    std::vector<RmapFragment> out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        auto f = split_tab(line);
        if (f.size() < 4) {
            continue;
        }
        RmapFragment r;
        r.chrom = f[0];
        r.start = std::stol(f[1]);
        r.end = std::stol(f[2]);
        r.id = std::stol(f[3]);
        out.push_back(std::move(r));
    }
    return out;
}

std::vector<BaitmapFragment> read_baitmap(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open baitmap file: " + path);
    }
    std::vector<BaitmapFragment> out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        auto f = split_tab(line);
        if (f.size() < 4) {
            continue;
        }
        BaitmapFragment b;
        b.chrom = f[0];
        b.start = std::stol(f[1]);
        b.end = std::stol(f[2]);
        b.id = std::stol(f[3]);
        b.name = (f.size() > 4) ? f[4] : std::string();
        out.push_back(std::move(b));
    }
    return out;
}

std::vector<ChinputRecord> read_chinput(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open chinput file: " + path);
    }
    std::vector<ChinputRecord> out;
    std::string line;
    bool header_skipped = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (!header_skipped) {
            // First non-comment line is the column header
            // (baitID otherEndID N otherEndLen distSign).
            header_skipped = true;
            if (line.rfind("baitID", 0) == 0) {
                continue;
            }
        }
        auto f = split_tab(line);
        if (f.size() < 5) {
            continue;
        }
        ChinputRecord r;
        r.bait_id = std::stol(f[0]);
        r.other_end_id = std::stol(f[1]);
        r.N = std::stod(f[2]);
        r.other_end_len = std::stol(f[3]);
        if (f[4] == "NA") {
            r.has_dist_sign = false;
        } else {
            r.has_dist_sign = true;
            r.dist_sign = std::stol(f[4]);
        }
        out.push_back(std::move(r));
    }
    return out;
}

namespace {

// floor(x + 0.5): round-half-up. Matches the .chinput distSign column's own
// midpoint rounding exactly (verified against all 270441 cis rows of
// GM_rep1.chinput: 100% agreement, where round-to-even or truncation each
// mismatch on tens of thousands of rows).
double half_up(double x) { return std::floor(x + 0.5); }

}  // namespace

std::vector<ChinputRecord> chinput_from_matrices(const std::vector<std::string>& matrix_paths,
                                                  const std::vector<RmapFragment>& rmap,
                                                  const std::vector<BaitmapFragment>& baitmap,
                                                  const FilterSettings& fs) {
    // Other-end fragments grouped by chromosome, sorted by start, so the
    // fragments within max_l_brown_est of a bait can be found with a binary
    // search rather than a scan of the whole chromosome.
    std::unordered_map<std::string, std::vector<const RmapFragment*>> by_chrom;
    for (const auto& f : rmap) by_chrom[f.chrom].push_back(&f);
    for (auto& kv : by_chrom) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const RmapFragment* a, const RmapFragment* b) { return a->start < b->start; });
    }

    // Only the chromosomes the .baitmap actually uses are ever loaded from a
    // matrix: chinput_from_matrices only derives cis pairs (see the header
    // comment), so no other chromosome's data is needed.
    std::set<std::string> bait_chroms;
    for (const auto& b : baitmap) bait_chroms.insert(b.chrom);

    std::map<std::pair<long, long>, double> n_sum;         // (baitID, otherEndID) -> N
    std::unordered_map<long, long> other_end_len;          // otherEndID -> end - start
    std::map<std::pair<long, long>, long> dist_sign_of;    // (baitID, otherEndID) -> distSign

    for (const std::string& path : matrix_paths) {
        for (const std::string& chrom : bait_chroms) {
            hicx::ToolMatrix matrix;
            try {
                matrix = hicx::ToolMatrix::load(path, chrom);
            } catch (const std::exception&) {
                continue;  // this matrix has no data for this chromosome
            }
            if (matrix.cut_intervals().empty()) continue;
            const hicx::BinTable bins(matrix.cut_intervals());

            for (const auto& bait : baitmap) {
                if (bait.chrom != chrom) continue;
                const auto bait_range = bins.region_bin_range(chrom, bait.start - 1, bait.end - 1);
                if (!bait_range.has_value()) continue;
                const auto [bait_row_first, bait_row_last] = *bait_range;

                const double bait_mid = static_cast<double>(bait.start + bait.end) / 2.0;
                const auto it = by_chrom.find(chrom);
                if (it == by_chrom.end()) continue;
                const std::vector<const RmapFragment*>& frags = it->second;

                // Binary search for the first fragment whose end could still
                // be within max_l_brown_est upstream of the bait, then walk
                // forward until the fragment starts more than max_l_brown_est
                // downstream of the bait.
                const long lo_bound = bait.start - fs.max_l_brown_est;
                auto lower = std::lower_bound(
                    frags.begin(), frags.end(), lo_bound,
                    [](const RmapFragment* f, long value) { return f->end < value; });

                for (auto fit = lower; fit != frags.end(); ++fit) {
                    const RmapFragment& oe = **fit;
                    if (oe.start - bait.end > fs.max_l_brown_est) break;
                    if (oe.id == bait.id) continue;  // self, not an other end

                    const double oe_mid = static_cast<double>(oe.start + oe.end) / 2.0;
                    const long dist_sign = static_cast<long>(half_up(oe_mid) - half_up(bait_mid));
                    if (std::labs(dist_sign) >= fs.max_l_brown_est) continue;

                    const auto oe_range = bins.region_bin_range(chrom, oe.start - 1, oe.end - 1);
                    if (!oe_range.has_value()) continue;
                    const auto [oe_col_first, oe_col_last] = *oe_range;

                    double value = 0.0;
                    for (std::int64_t row = bait_row_first; row <= bait_row_last; ++row) {
                        for (std::int64_t col = oe_col_first; col <= oe_col_last; ++col) {
                            value += matrix.matrix().at(row, col);
                        }
                    }
                    if (value == 0.0) continue;  // no observed contact, no chinput row

                    n_sum[{bait.id, oe.id}] += value;
                    other_end_len[oe.id] = oe.end - oe.start;
                    dist_sign_of[{bait.id, oe.id}] = dist_sign;
                }
            }
        }
    }

    std::vector<ChinputRecord> out;
    out.reserve(n_sum.size());
    for (const auto& [key, n] : n_sum) {
        ChinputRecord r;
        r.bait_id = key.first;
        r.other_end_id = key.second;
        r.N = n;
        r.other_end_len = other_end_len.at(key.second);
        r.has_dist_sign = true;
        r.dist_sign = dist_sign_of.at(key);
        out.push_back(std::move(r));
    }
    return out;
}

double avg_frag_length(const std::vector<RmapFragment>& rmap) {
    std::map<std::string, long> chr_max;
    for (const auto& r : rmap) {
        if (r.chrom == "MT" || r.chrom == "chrMT" || r.chrom == "M" || r.chrom == "chrM") {
            continue;
        }
        auto it = chr_max.find(r.chrom);
        if (it == chr_max.end() || r.end > it->second) {
            chr_max[r.chrom] = r.end;
        }
    }
    double sum = 0.0;
    for (const auto& kv : chr_max) {
        sum += static_cast<double>(kv.second);
    }
    return sum / static_cast<double>(rmap.size());
}

std::size_t n_hypotheses(const std::vector<RmapFragment>& rmap,
                          const std::vector<BaitmapFragment>& baitmap) {
    const double n_baits = static_cast<double>(baitmap.size());
    const double n_rmap = static_cast<double>(rmap.size());
    const double nhyp = n_baits * (2.0 * n_rmap - n_baits - 1.0) / 2.0;
    return static_cast<std::size_t>(std::llround(nhyp));
}

double eta_bar_from_design(const std::vector<RmapFragment>& rmap,
                            const std::vector<BaitmapFragment>& baitmap, const WeightSettings& w,
                            double avg_frag_len, bool /*include_trans*/,
                            std::size_t n_hypotheses) {
    std::map<std::string, long> chr_max;
    for (const auto& r : rmap) {
        auto it = chr_max.find(r.chrom);
        if (it == chr_max.end() || r.end > it->second) {
            chr_max[r.chrom] = r.end;
        }
    }
    std::map<std::string, long> n_baits;
    for (const auto& b : baitmap) {
        n_baits[b.chrom]++;
    }

    double eta_sigma = 0.0;
    for (const auto& kv : n_baits) {
        const std::string& chrom = kv.first;
        if (chrom == "MT" || chrom == "chrMT") {
            continue;
        }
        auto mit = chr_max.find(chrom);
        if (mit == chr_max.end()) {
            continue;
        }
        const double d_c = static_cast<double>(mit->second);
        const long n_c = kv.second;
        for (long i = 1; i <= n_c; ++i) {
            const double d = d_c * static_cast<double>(i) / static_cast<double>(n_c);
            const double d_near = std::min(d, d_c - d);
            // d.other <- seq(avgFragLen, max(avgFragLen, d.near), by = avgFragLen)
            const double upper1 = std::max(avg_frag_len, d_near);
            for (double x = avg_frag_len; x <= upper1 + 1e-9; x += avg_frag_len) {
                eta_sigma += 2.0 * expit(w.alpha + w.beta * std::log(x));
            }
            // d.other2 <- seq(d.near, d.c - d.near, by = avgFragLen)
            const double upper2 = d_c - d_near;
            for (double x = d_near; x <= upper2 + 1e-9; x += avg_frag_len) {
                eta_sigma += expit(w.alpha + w.beta * std::log(std::max(x, 1e-300)));
            }
        }
    }
    return eta_sigma / static_cast<double>(n_hypotheses);
}

DistFunFit fit_distance_function(const std::vector<double>& midpoints,
                                  const std::vector<double>& ref_bin_mean) {
    if (midpoints.size() != ref_bin_mean.size() || midpoints.empty()) {
        throw std::runtime_error("fit_distance_function: mismatched or empty input");
    }
    const std::size_t n = midpoints.size();
    // OLS design matrix: [1, log(m), log(m)^2, log(m)^3], normal equations
    // solved with a plain 4x4 Gaussian elimination (n is a few hundred rows
    // at most: one entry per distance bin, not per interaction).
    double XtX[4][4] = {};
    double Xty[4] = {};
    for (std::size_t i = 0; i < n; ++i) {
        const double lm = std::log(midpoints[i]);
        const double ly = std::log(ref_bin_mean[i]);
        const double row[4] = {1.0, lm, lm * lm, lm * lm * lm};
        for (int a = 0; a < 4; ++a) {
            Xty[a] += row[a] * ly;
            for (int b = 0; b < 4; ++b) {
                XtX[a][b] += row[a] * row[b];
            }
        }
    }
    // Gaussian elimination with partial pivoting.
    double A[4][5];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) A[i][j] = XtX[i][j];
        A[i][4] = Xty[i];
    }
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int r = col + 1; r < 4; ++r) {
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col])) piv = r;
        }
        std::swap(A[col], A[piv]);
        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            const double factor = A[r][col] / A[col][col];
            for (int c = col; c < 5; ++c) A[r][c] -= factor * A[col][c];
        }
    }
    double coef[4];
    for (int i = 0; i < 4; ++i) coef[i] = A[i][4] / A[i][i];

    DistFunFit fit;
    for (int i = 0; i < 4; ++i) fit.cubic[i] = coef[i];

    double obs_min = midpoints[0], obs_max = midpoints[0];
    for (double m : midpoints) {
        obs_min = std::min(obs_min, m);
        obs_max = std::max(obs_max, m);
    }
    fit.obs_min_log = std::log(obs_min);
    fit.obs_max_log = std::log(obs_max);

    const double a0 = coef[0], a1 = coef[1], a2 = coef[2], a3 = coef[3];
    auto beta_at = [&](double lx) { return a1 + 2.0 * a2 * lx + 3.0 * a3 * lx * lx; };
    auto alpha_at = [&](double lx, double beta) {
        return a0 + (a1 - beta) * lx + a2 * lx * lx + a3 * lx * lx * lx;
    };
    const double beta_min = beta_at(fit.obs_min_log);
    const double beta_max = beta_at(fit.obs_max_log);
    fit.head_coef[0] = alpha_at(fit.obs_min_log, beta_min);
    fit.head_coef[1] = beta_min;
    fit.tail_coef[0] = alpha_at(fit.obs_max_log, beta_max);
    fit.tail_coef[1] = beta_max;
    return fit;
}

double eval_distance_function(const DistFunFit& fit, double distance) {
    const double lx = std::log(distance);
    if (lx < fit.obs_min_log) {
        return fit.head_coef[0] + fit.head_coef[1] * lx;
    }
    if (lx > fit.obs_max_log) {
        return fit.tail_coef[0] + fit.tail_coef[1] * lx;
    }
    const double a0 = fit.cubic[0], a1 = fit.cubic[1], a2 = fit.cubic[2], a3 = fit.cubic[3];
    return a0 + a1 * lx + a2 * lx * lx + a3 * lx * lx * lx;
}

// ---------------------------------------------------------------------------
// Genome-wide parameter estimation (see hicx/chicago.hpp for scope notes).
// ---------------------------------------------------------------------------

namespace {

double median_of(std::vector<double> v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// R's default (type 7) sample quantile.
double quantile_type7(std::vector<double> v, double p) {
    std::sort(v.begin(), v.end());
    const long n = static_cast<long>(v.size());
    if (n == 1) return v[0];
    const double h = (static_cast<double>(n) - 1.0) * p;
    const long lo = static_cast<long>(std::floor(h));
    const long hi = std::min(lo + 1, n - 1);
    const double frac = h - static_cast<double>(lo);
    return v[lo] + frac * (v[hi] - v[lo]);
}

// geo_mean (Chicago): exp(mean(log(finite values))).
double geo_mean(const std::vector<double>& v) {
    double sum = 0.0;
    long count = 0;
    for (double x : v) {
        const double lx = std::log(x);
        if (std::isfinite(lx)) {
            sum += lx;
            ++count;
        }
    }
    if (count == 0) return std::numeric_limits<double>::quiet_NaN();
    return std::exp(sum / static_cast<double>(count));
}

std::unordered_map<long, std::vector<double>> read_bin_table(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open design table: " + path);
    }
    std::unordered_map<long, std::vector<double>> out;
    std::string line;
    bool header_skipped = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!header_skipped) {
            header_skipped = true;
            if (!line.empty() && line[0] == '#') continue;
        }
        auto f = split_tab(line);
        if (f.size() < 2) continue;
        const long id = std::stol(f[0]);
        std::vector<double> bins;
        bins.reserve(f.size() - 1);
        for (std::size_t i = 1; i < f.size(); ++i) bins.push_back(std::stod(f[i]));
        out.emplace(id, std::move(bins));
    }
    return out;
}

// distbin(abs(distSign)) for the fixed-width R cut(x, seq(0, maxL, bin)),
// right-closed intervals: 0 means undefined (out of (0, maxL]).
long distbin_of(double abs_dist, const FilterSettings& fs) {
    if (!(abs_dist > 0.0) || abs_dist > static_cast<double>(fs.max_l_brown_est)) return 0;
    return static_cast<long>(std::ceil(abs_dist / static_cast<double>(fs.binsize)));
}

}  // namespace

std::vector<ChiInteraction> read_sample(const std::vector<ChinputRecord>& raw,
                                         const std::vector<BaitmapFragment>& baitmap,
                                         const FilterSettings& fs) {
    std::unordered_set<long> bait_ids;
    for (const auto& b : baitmap) bait_ids.insert(b.id);

    std::vector<ChiInteraction> x;
    x.reserve(raw.size());
    for (const auto& r : raw) {
        if (r.other_end_len < fs.min_frag_len || r.other_end_len > fs.max_frag_len) continue;
        if (r.has_dist_sign && r.dist_sign == 0) continue;  // self-ligation
        ChiInteraction c;
        c.bait_id = r.bait_id;
        c.other_end_id = r.other_end_id;
        c.N = r.N;
        c.other_end_len = r.other_end_len;
        c.has_dist_sign = r.has_dist_sign;
        c.dist_sign = r.dist_sign;
        c.is_bait2bait = bait_ids.count(r.other_end_id) > 0;
        x.push_back(std::move(c));
    }

    // minNPerBait
    std::unordered_map<long, double> n_per_bait;
    for (const auto& c : x) n_per_bait[c.bait_id] += c.N;
    std::vector<ChiInteraction> x2;
    x2.reserve(x.size());
    for (auto& c : x) {
        if (n_per_bait[c.bait_id] >= static_cast<double>(fs.min_n_per_bait)) x2.push_back(c);
    }
    x = std::move(x2);

    if (fs.remove_adjacent) {
        std::vector<ChiInteraction> x3;
        x3.reserve(x.size());
        for (auto& c : x) {
            if (std::labs(c.bait_id - c.other_end_id) != 1) x3.push_back(c);
        }
        x = std::move(x3);
    }

    // Drop baits whose every proximal (|distSign| < maxLBrownEst, cis) row is bait2bait.
    std::unordered_map<long, bool> bait_has_prox;       // seen a proximal row at all
    std::unordered_map<long, bool> bait_all_b2b;        // all proximal rows seen so far are b2b
    for (auto& c : x) {
        if (!c.has_dist_sign) continue;
        if (!(std::labs(c.dist_sign) < fs.max_l_brown_est)) continue;
        auto it = bait_has_prox.find(c.bait_id);
        if (it == bait_has_prox.end()) {
            bait_has_prox[c.bait_id] = true;
            bait_all_b2b[c.bait_id] = c.is_bait2bait;
        } else {
            bait_all_b2b[c.bait_id] = bait_all_b2b[c.bait_id] && c.is_bait2bait;
        }
    }
    std::vector<ChiInteraction> x4;
    x4.reserve(x.size());
    for (auto& c : x) {
        auto it = bait_has_prox.find(c.bait_id);
        // isAllB2BProx is TRUE (bait dropped) when there IS at least one proximal
        // row and every one of them is bait2bait; a bait with no proximal row at
        // all keeps isAllB2BProx == TRUE too (R: `if (!length(prox)) TRUE`).
        const bool all_b2b = (it == bait_has_prox.end()) ? true : bait_all_b2b[c.bait_id];
        if (!all_b2b) x4.push_back(c);
    }
    return x4;
}

std::vector<int> cut2_native_groups(const std::vector<double>& x, long m) {
    std::vector<double> xx;
    std::vector<long> cum;
    {
        std::vector<double> sorted = x;
        std::sort(sorted.begin(), sorted.end());
        long running = 0;
        for (std::size_t i = 0; i < sorted.size();) {
            std::size_t j = i;
            while (j < sorted.size() && sorted[j] == sorted[i]) ++j;
            running += static_cast<long>(j - i);
            xx.push_back(sorted[i]);
            cum.push_back(running);
            i = j;
        }
    }
    const long nnm = static_cast<long>(x.size());
    const long g = std::max(1L, nnm / std::max(1L, m));
    std::vector<double> cuts(g);
    for (long j = 1; j <= g; ++j) {
        const double xout = static_cast<double>(j) * (static_cast<double>(nnm) / static_cast<double>(g));
        auto it = std::lower_bound(cum.begin(), cum.end(), xout,
                                    [](long c, double v) { return static_cast<double>(c) < v; });
        std::size_t idx = static_cast<std::size_t>(it - cum.begin());
        if (idx >= xx.size()) idx = xx.size() - 1;
        cuts[j - 1] = xx[idx];
    }
    cuts.back() = xx.back();

    double min_dif = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < xx.size(); ++i) min_dif = std::min(min_dif, xx[i] - xx[i - 1]);
    min_dif *= 0.5;

    std::vector<double> low(g), up(g);
    double lower = xx.front();
    double upper = 1e45;
    long i = 0;
    std::vector<int> y(x.size(), 0);
    for (long j = 1; j <= g; ++j) {
        const double cj = cuts[j - 1];
        if (cj == upper) continue;
        ++i;
        upper = cj;
        for (std::size_t k = 0; k < x.size(); ++k) {
            if (x[k] >= lower - min_dif) y[k] = static_cast<int>(i);
        }
        low[i - 1] = lower;
        if (j == g) {
            lower = upper;
        } else {
            auto it = std::upper_bound(xx.begin(), xx.end(), upper);
            lower = (it == xx.end()) ? upper : *it;
        }
        up[i - 1] = lower;
    }
    return y;
}

std::vector<double> cut2_cuts(const std::vector<double>& x, long m) {
    auto groups = cut2_native_groups(x, m);
    // Recover 'low' by re-deriving group boundaries: for each group id, the
    // minimum x value assigned to it is 'low[group]'; the breaks vector is
    // unique(c(low[1..i], max(x))), i.e. the sorted set of per-group minima
    // plus the overall maximum.
    std::map<int, double> low_of_group;
    double xmax = -std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < x.size(); ++k) {
        auto it = low_of_group.find(groups[k]);
        if (it == low_of_group.end() || x[k] < it->second) low_of_group[groups[k]] = x[k];
        xmax = std::max(xmax, x[k]);
    }
    std::vector<double> breaks;
    for (auto& [g, low] : low_of_group) breaks.push_back(low);
    breaks.push_back(xmax);
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
    return breaks;
}

int cut_with_breaks(double x, const std::vector<double>& breaks) {
    if (breaks.size() < 2) return 0;
    if (x == breaks.front()) return 1;
    if (x < breaks.front() || x > breaks.back()) return 0;
    for (std::size_t i = 1; i < breaks.size(); ++i) {
        if (x > breaks[i - 1] && x <= breaks[i]) return static_cast<int>(i);
    }
    return 0;
}

TlbResult add_tlb(const std::vector<ChiInteraction>& x, const FilterSettings& fs,
                   double tlb_filter_top_percent, long tlb_min_prox_oe_per_bin,
                   long tlb_min_prox_b2b_per_bin) {
    struct OeAgg {
        double trans_length = 0.0;
        bool is_bait2bait = false;
        double min_abs_dist = std::numeric_limits<double>::infinity();
    };
    std::map<long, OeAgg> agg;  // sorted by otherEndID, matching data.table's key order
    for (const auto& c : x) {
        OeAgg& a = agg[c.other_end_id];
        a.is_bait2bait = c.is_bait2bait;
        if (!c.has_dist_sign) {
            a.trans_length += 1.0;
        } else {
            a.min_abs_dist = std::min(a.min_abs_dist, static_cast<double>(std::labs(c.dist_sign)));
        }
    }

    std::vector<long> oe_ids;
    std::vector<double> lengths;
    for (auto& [id, a] : agg) {
        oe_ids.push_back(id);
        lengths.push_back(a.trans_length);
    }
    const double q = quantile_type7(lengths, 1.0 - tlb_filter_top_percent / 100.0);
    std::vector<long> kept_non_b2b, kept_b2b;
    std::vector<double> len_non_b2b, len_b2b;
    std::vector<double> prox_len_non_b2b, prox_len_b2b;  // subset with min_abs_dist <= maxLBrownEst
    for (long id : oe_ids) {
        const OeAgg& a = agg.at(id);
        if (a.trans_length > q) continue;  // top-percent filter: dropped entirely
        if (a.is_bait2bait) {
            kept_b2b.push_back(id);
            len_b2b.push_back(a.trans_length);
            if (a.min_abs_dist <= static_cast<double>(fs.max_l_brown_est)) prox_len_b2b.push_back(a.trans_length);
        } else {
            kept_non_b2b.push_back(id);
            len_non_b2b.push_back(a.trans_length);
            if (a.min_abs_dist <= static_cast<double>(fs.max_l_brown_est)) prox_len_non_b2b.push_back(a.trans_length);
        }
    }

    TlbResult result;
    auto assign_pools = [&](const std::vector<long>& ids, const std::vector<double>& lens,
                             const std::vector<double>& prox_lens, long min_per_bin, int id_offset) {
        if (ids.empty()) return 0;
        if (prox_lens.empty()) {
            // No proximal other end at all in this class: assign everything to one pool.
            for (long id : ids) result.pool_of_other_end[id] = id_offset + 1;
            return 1;
        }
        std::vector<double> cuts = cut2_cuts(prox_lens, min_per_bin);
        if (cuts.size() == 1) {
            for (long id : ids) result.pool_of_other_end[id] = id_offset + 1;
            return 1;
        }
        double lo = *std::min_element(lens.begin(), lens.end());
        double hi = *std::max_element(lens.begin(), lens.end());
        if (cuts.front() > lo) cuts.front() = lo;
        if (cuts.back() < hi) cuts.back() = hi;
        int max_pool = 0;
        for (std::size_t k = 0; k < ids.size(); ++k) {
            int p = cut_with_breaks(lens[k], cuts);
            if (p == 0) p = 1;  // defensive: clamped range should make this unreachable
            result.pool_of_other_end[ids[k]] = id_offset + p;
            max_pool = std::max(max_pool, p);
        }
        return max_pool;
    };

    result.n_non_b2b_pools =
        assign_pools(kept_non_b2b, len_non_b2b, prox_len_non_b2b, tlb_min_prox_oe_per_bin, 0);
    result.n_b2b_pools = assign_pools(kept_b2b, len_b2b, prox_len_b2b, tlb_min_prox_b2b_per_bin,
                                       result.n_non_b2b_pools);
    return result;
}

TechnicalNoiseResult estimate_technical_noise(const std::vector<ChiInteraction>& x, const TlbResult& tlb,
                                               const std::vector<RmapFragment>& rmap,
                                               const std::vector<BaitmapFragment>& baitmap,
                                               long min_baits_per_bin) {
    std::unordered_map<long, std::string> bait_chrom, oe_chrom;
    for (const auto& b : baitmap) bait_chrom[b.id] = b.chrom;
    for (const auto& r : rmap) oe_chrom[r.id] = r.chrom;

    std::map<long, long> trans_count_by_bait;  // sorted by baitID
    for (const auto& c : x) {
        auto it = trans_count_by_bait.find(c.bait_id);
        if (it == trans_count_by_bait.end()) trans_count_by_bait[c.bait_id] = 0;
    }
    for (const auto& c : x) {
        if (!c.has_dist_sign) trans_count_by_bait[c.bait_id] += 1;
    }
    std::vector<long> bait_ids;
    std::vector<double> counts;
    for (auto& [id, cnt] : trans_count_by_bait) {
        bait_ids.push_back(id);
        counts.push_back(static_cast<double>(cnt));
    }
    auto groups = cut2_native_groups(counts, min_baits_per_bin);

    TechnicalNoiseResult result;
    for (std::size_t i = 0; i < bait_ids.size(); ++i) result.tblb_of_bait[bait_ids[i]] = groups[i];

    struct PoolAgg {
        std::set<long> baits, oes;
        double n_trans = 0.0;
    };
    std::map<std::pair<int, int>, PoolAgg> pools;
    for (const auto& c : x) {
        auto oe_it = tlb.pool_of_other_end.find(c.other_end_id);
        if (oe_it == tlb.pool_of_other_end.end()) continue;  // dropped by addTLB's top-percent filter
        const int tlb_id = oe_it->second;
        const int tblb_id = result.tblb_of_bait.at(c.bait_id);
        PoolAgg& a = pools[{tlb_id, tblb_id}];
        a.baits.insert(c.bait_id);
        a.oes.insert(c.other_end_id);
        if (!c.has_dist_sign) a.n_trans += c.N;
    }

    for (auto& [key, a] : pools) {
        std::map<std::string, long> bait_chrom_count, oe_chrom_count;
        for (long id : a.baits) bait_chrom_count[bait_chrom.at(id)]++;
        for (long id : a.oes) oe_chrom_count[oe_chrom.at(id)]++;
        const long total_oes = static_cast<long>(a.oes.size());
        double num_pairs = 0.0;
        for (auto& [chrom, cnt] : bait_chrom_count) {
            auto it = oe_chrom_count.find(chrom);
            const long same_chrom = (it == oe_chrom_count.end()) ? 0 : it->second;
            num_pairs += static_cast<double>(cnt) * static_cast<double>(total_oes - same_chrom);
        }
        long overlap = 0;
        for (long id : a.baits) {
            if (a.oes.count(id)) ++overlap;
        }
        num_pairs -= static_cast<double>(overlap);
        result.tmean_by_pool[key] = a.n_trans / num_pairs;
    }
    return result;
}

BaitFactors normalise_baits(const std::vector<ChiInteraction>& x, const std::vector<RmapFragment>&,
                             const std::vector<BaitmapFragment>&, const std::string& npb_path,
                             const FilterSettings& fs) {
    auto npb = read_bin_table(npb_path);
    long n_bins = 0;
    for (auto& [id, v] : npb) {
        n_bins = std::max<long>(n_bins, static_cast<long>(v.size()));
    }

    // sum(N) per (baitID, distbin), non-bait2bait, distbin defined.
    std::map<std::pair<long, long>, double> sum_n;  // (baitID, distbin) -> sum N
    for (const auto& c : x) {
        if (c.is_bait2bait) continue;
        if (!c.has_dist_sign) continue;
        const long db = distbin_of(static_cast<double>(std::labs(c.dist_sign)), fs);
        if (db == 0) continue;
        sum_n[{c.bait_id, db}] += c.N;
    }

    // bbm per (baitID, distbin) = sum_n / ntot, ntot = npb[baitID][distbin].
    std::map<long, std::vector<std::pair<long, double>>> bbm_by_distbin;  // distbin -> [(baitID, bbm)]
    for (auto& [key, s] : sum_n) {
        const auto [bait_id, db] = key;
        auto it = npb.find(bait_id);
        if (it == npb.end() || db < 1 || static_cast<std::size_t>(db) > it->second.size()) continue;
        const double ntot = it->second[db - 1];
        if (!(ntot > 0.0)) continue;
        const double bbm = s / ntot;
        bbm_by_distbin[db].emplace_back(bait_id, bbm);
    }

    BaitFactors out;
    out.ref_bin_mean_by_distbin.assign(static_cast<std::size_t>(n_bins) + 1,
                                        std::numeric_limits<double>::quiet_NaN());
    std::map<long, std::vector<double>> siv_by_bait;
    for (auto& [db, entries] : bbm_by_distbin) {
        std::vector<double> vals;
        vals.reserve(entries.size());
        for (auto& [bid, bbm] : entries) vals.push_back(bbm);
        const double gm = geo_mean(vals);
        if (db >= 0 && static_cast<std::size_t>(db) < out.ref_bin_mean_by_distbin.size()) {
            out.ref_bin_mean_by_distbin[static_cast<std::size_t>(db)] = gm;
        }
        for (auto& [bid, bbm] : entries) {
            if (gm > 0.0 && std::isfinite(gm)) siv_by_bait[bid].push_back(bbm / gm);
        }
    }
    for (auto& [bid, vals] : siv_by_bait) {
        const double med = median_of(vals);
        if (std::isfinite(med)) out.s_j[bid] = med;
    }
    return out;
}

std::unordered_map<int, double> normalise_other_ends(const std::vector<ChiInteraction>& x,
                                                       const TlbResult& tlb,
                                                       const std::unordered_map<long, double>& bait_s_j,
                                                       const std::string& nbpb_path,
                                                       const FilterSettings& fs) {
    auto nbpb = read_bin_table(nbpb_path);

    // sum(NNb) per (tlb pool, distbin), restricted to |distSign| <= maxLBrownEst,
    // cis only; NNb = pmax(1, round(N / s_j)), the bait-normalised count
    // (normaliseBaits' own output column, not the raw N).
    std::map<std::pair<int, long>, double> sum_n;
    // The (otherEndID, distbin) pairs actually present in the (filtered) data:
    // R dedupes nbpb-joined-onto-x by (otherEndID, distbin) before summing,
    // so an other end contributes nbpb[otherEndID][distbin] to its pool's
    // ntot once per distinct distbin it has real interactions in, not once
    // per every bin nbpb happens to define.
    std::set<std::pair<long, long>> oe_distbin_present;
    for (const auto& c : x) {
        if (!c.has_dist_sign) continue;
        const double ad = static_cast<double>(std::labs(c.dist_sign));
        if (ad > static_cast<double>(fs.max_l_brown_est)) continue;
        auto it = tlb.pool_of_other_end.find(c.other_end_id);
        if (it == tlb.pool_of_other_end.end()) continue;
        auto sj_it = bait_s_j.find(c.bait_id);
        if (sj_it == bait_s_j.end()) continue;
        const long db = distbin_of(ad, fs);
        if (db == 0) continue;
        const double nnb = std::max(1.0, std::round(c.N / sj_it->second));
        sum_n[std::make_pair(it->second, db)] += nnb;
        oe_distbin_present.insert(std::make_pair(c.other_end_id, db));
    }

    // ntot per (tlb pool, distbin): sum over other ends in that pool, over
    // the distbins they actually appear in, of nbpb[otherEndID][distbin]
    // (matching R's nbpbSum = nbpb[, sum(bin<k>), by=(tlb,distbin)] applied
    // to the deduped, per-observed-pair nbpb join).
    std::map<std::pair<int, long>, double> ntot;
    for (const auto& [oe_id, db] : oe_distbin_present) {
        auto pool_it = tlb.pool_of_other_end.find(oe_id);
        if (pool_it == tlb.pool_of_other_end.end()) continue;
        auto it = nbpb.find(oe_id);
        if (it == nbpb.end() || db < 1 || static_cast<std::size_t>(db) > it->second.size()) continue;
        ntot[std::make_pair(pool_it->second, db)] += it->second[db - 1];
    }

    std::map<long, std::vector<std::pair<int, double>>> bbm_by_distbin;
    for (auto& [key, s] : sum_n) {
        const auto [pool, db] = key;
        auto it = ntot.find({pool, db});
        if (it == ntot.end() || !(it->second > 0.0)) continue;
        bbm_by_distbin[db].emplace_back(pool, s / it->second);
    }

    // .normaliseFragmentSets' refExcludeSuffix = "B2B" (as normaliseOtherEnds
    // always passes it): the per-distbin geometric mean that every pool's bbm
    // is normalised against is computed from the non-bait2bait pools only,
    // so a bait2bait-only other-end pool never dilutes the reference scale.
    std::map<int, std::vector<double>> siv_by_pool;
    for (auto& [db, entries] : bbm_by_distbin) {
        std::vector<double> ref_vals;
        for (auto& [pool, bbm] : entries) {
            if (pool > tlb.n_non_b2b_pools) continue;  // bait2bait pool: excluded from the reference
            ref_vals.push_back(bbm);
        }
        const double gm = geo_mean(ref_vals);
        if (!(gm > 0.0) || !std::isfinite(gm)) continue;
        for (auto& [pool, bbm] : entries) siv_by_pool[pool].push_back(bbm / gm);
    }
    std::unordered_map<int, double> s_i;
    for (auto& [pool, vals] : siv_by_pool) {
        const double med = median_of(vals);
        if (std::isfinite(med)) s_i[pool] = med;
    }
    return s_i;
}

std::vector<ProxOePair> read_poe(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open proxOE file: " + path);
    }
    std::vector<ProxOePair> out;
    std::string line;
    bool header_skipped = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!header_skipped) {
            header_skipped = true;
            if (line[0] == '#') continue;
        }
        auto f = split_tab(line);
        if (f.size() < 3) continue;
        ProxOePair p;
        p.bait_id = std::stol(f[0]);
        p.other_end_id = std::stol(f[1]);
        p.dist = std::stod(f[2]);
        out.push_back(p);
    }
    return out;
}

namespace {

// Splits [0, n) into up to `threads` contiguous chunks, runs `partial` on each
// in its own std::thread (partial(first, last) -> double), and sums the
// per-chunk results in chunk order, so the reduction order (and therefore the
// floating-point result) does not depend on how many threads ran it.
double threaded_reduce(std::size_t n, int threads,
                        const std::function<double(std::size_t, std::size_t)>& partial) {
    const int worker_count = std::max(1, threads);
    if (worker_count == 1 || n < 4096) {
        return partial(0, n);
    }
    std::vector<double> sums(static_cast<std::size_t>(worker_count), 0.0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));
    const std::size_t chunk = (n + static_cast<std::size_t>(worker_count) - 1) /
                               static_cast<std::size_t>(worker_count);
    for (int w = 0; w < worker_count; ++w) {
        const std::size_t first = static_cast<std::size_t>(w) * chunk;
        const std::size_t last = std::min(n, first + chunk);
        if (first >= last) break;
        workers.emplace_back([&, w, first, last] { sums[static_cast<std::size_t>(w)] = partial(first, last); });
    }
    for (std::thread& worker : workers) worker.join();
    double total = 0.0;
    for (double s : sums) total += s;
    return total;
}

}  // namespace

double estimate_dispersion_theta_ml(const std::vector<double>& N, const std::vector<double>& Bmean,
                                     int threads) {
    if (N.size() != Bmean.size() || N.empty()) {
        throw std::runtime_error("estimate_dispersion_theta_ml: mismatched or empty input");
    }
    const std::size_t n = N.size();
    auto score = [&](double th) {
        return threaded_reduce(n, threads, [&](std::size_t first, std::size_t last) {
            double s = 0.0;
            for (std::size_t i = first; i < last; ++i) {
                const double mu = Bmean[i];
                s += boost::math::digamma(th + N[i]) - boost::math::digamma(th) + std::log(th) + 1.0 -
                     std::log(th + mu) - (N[i] + th) / (mu + th);
            }
            return s;
        });
    };
    auto info = [&](double th) {
        return threaded_reduce(n, threads, [&](std::size_t first, std::size_t last) {
            double s = 0.0;
            for (std::size_t i = first; i < last; ++i) {
                const double mu = Bmean[i];
                s += -boost::math::trigamma(th + N[i]) + boost::math::trigamma(th) - 1.0 / th +
                     2.0 / (mu + th) - (N[i] + th) / ((mu + th) * (mu + th));
            }
            return s;
        });
    };

    const double sum_sq = threaded_reduce(n, threads, [&](std::size_t first, std::size_t last) {
        double s = 0.0;
        for (std::size_t i = first; i < last; ++i) {
            const double r = N[i] / Bmean[i] - 1.0;
            s += r * r;
        }
        return s;
    });
    double t0 = static_cast<double>(n) / sum_sq;
    long it = 0;
    double del = 1.0;
    const double eps = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
    // R's own MASS::theta.ml defaults limit = 10, which measurably (checked
    // against a real 433,821-row GM12878 fit) still trips its own "iteration
    // limit reached" warning; MASS::glm.nb's outer loop keeps refining theta
    // past that point and settles at the fully converged Newton root. This
    // reproduces that converged root directly (same score/info, same start,
    // just enough iterations to actually reach fabs(del) <= eps), which
    // matches glm.nb's reported theta to 6 significant digits on that fit.
    const long limit = 200;
    while (++it < limit && std::fabs(del) > eps) {
        t0 = std::fabs(t0);
        const double inf = info(t0);
        del = score(t0) / inf;
        t0 += del;
    }
    if (t0 < 0.0) t0 = 0.0;
    return t0;
}

double estimate_bmean(double s_j, double s_i, double abs_dist_or_nan, const DistFunFit& fit) {
    if (std::isnan(abs_dist_or_nan)) return 0.0;
    return s_j * s_i * std::exp(eval_distance_function(fit, abs_dist_or_nan));
}

BackgroundModel fit_chicago_background(const std::vector<ChinputRecord>& raw,
                                        const std::vector<RmapFragment>& rmap,
                                        const std::vector<BaitmapFragment>& baitmap,
                                        const std::string& npb_path, const std::string& nbpb_path,
                                        const std::string& poe_path, const FilterSettings& fs,
                                        long tlb_min_baits_per_bin, double tlb_filter_top_percent,
                                        long tlb_min_prox_oe_per_bin, long tlb_min_prox_b2b_per_bin,
                                        long brownian_noise_subset, int threads) {
    BackgroundModel model;
    model.filters = fs;

    const std::vector<ChiInteraction> x = read_sample(raw, baitmap, fs);
    model.tlb = add_tlb(x, fs, tlb_filter_top_percent, tlb_min_prox_oe_per_bin, tlb_min_prox_b2b_per_bin);
    model.bait_factors = normalise_baits(x, rmap, baitmap, npb_path, fs);
    model.s_i_by_tlb_pool = normalise_other_ends(x, model.tlb, model.bait_factors.s_j, nbpb_path, fs);
    model.tech_noise = estimate_technical_noise(x, model.tlb, rmap, baitmap, tlb_min_baits_per_bin);

    // estimateDistFun: sequential midpoints over the distbins with a defined
    // refBinMean, in ascending distbin order (R builds this from a sorted,
    // unique-by-distbin table, so gaps in the distbin sequence do not shift
    // the midpoint sequence: see the header comment on fit_distance_function).
    std::vector<double> midpoints, ref_bin_means;
    {
        long position = 0;
        for (std::size_t db = 1; db < model.bait_factors.ref_bin_mean_by_distbin.size(); ++db) {
            const double rbm = model.bait_factors.ref_bin_mean_by_distbin[db];
            if (!std::isfinite(rbm)) continue;
            ++position;
            midpoints.push_back(std::round(static_cast<double>(fs.binsize) / 2.0) +
                                 static_cast<double>(position - 1) * static_cast<double>(fs.binsize));
            ref_bin_means.push_back(rbm);
        }
    }
    model.dist_fun = fit_distance_function(midpoints, ref_bin_means);

    // estimateBrownianComponent / .estimateDispersion: N/Bmean pairs over the
    // full .poe design, restricted to baits that have an s_j (R's sel.baits),
    // s_i defaulting to 1 for other ends the tlb pooling excluded (top-percent
    // filter), matching normaliseOtherEnds' own NA -> 1 fallback.
    std::unordered_map<long, double> n_of_pair;  // baitID*2^32 xor otherEndID-ish key, see below
    // A real hash of the (bait, otherEnd) pair; both ids fit comfortably in 32
    // bits for every design this project ships, so a 64-bit combination key
    // is exact and collision-free here.
    auto pair_key = [](long bait, long oe) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(bait)) << 32) |
               static_cast<std::uint32_t>(oe);
    };
    std::unordered_map<std::uint64_t, double> n_by_pair;
    for (const auto& c : x) {
        if (!c.has_dist_sign) continue;
        n_by_pair[pair_key(c.bait_id, c.other_end_id)] += c.N;
    }

    const auto poe = read_poe(poe_path);
    const std::size_t n_baits_with_s_j = model.bait_factors.s_j.size();
    model.subset_would_trigger_in_r = static_cast<long>(n_baits_with_s_j) > brownian_noise_subset;

    std::vector<double> Ns, Bmeans;
    Ns.reserve(poe.size());
    Bmeans.reserve(poe.size());
    for (const auto& p : poe) {
        auto sj_it = model.bait_factors.s_j.find(p.bait_id);
        if (sj_it == model.bait_factors.s_j.end()) continue;  // not a normalised bait
        double s_i = 1.0;
        auto oe_it = model.tlb.pool_of_other_end.find(p.other_end_id);
        if (oe_it != model.tlb.pool_of_other_end.end()) {
            auto si_it = model.s_i_by_tlb_pool.find(oe_it->second);
            if (si_it != model.s_i_by_tlb_pool.end()) s_i = si_it->second;
        }
        auto n_it = n_by_pair.find(pair_key(p.bait_id, p.other_end_id));
        const double N = (n_it == n_by_pair.end()) ? 0.0 : n_it->second;
        const double bmean = sj_it->second * s_i * std::exp(eval_distance_function(model.dist_fun, p.dist));
        Ns.push_back(N);
        Bmeans.push_back(bmean);
    }
    model.dispersion_n_pairs = Ns.size();
    model.dispersion = estimate_dispersion_theta_ml(Ns, Bmeans, threads);
    return model;
}

}  // namespace hicx::chicago
