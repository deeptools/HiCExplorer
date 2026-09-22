// See hicx/chicago.hpp.

#include "hicx/chicago.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>

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

}  // namespace hicx::chicago
