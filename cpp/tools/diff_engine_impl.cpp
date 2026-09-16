// The count model of hicDifferentialAnalysis; see diff_engine_impl.hpp.

#include "diff_engine_impl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "hicx/parallel.hpp"
#include "hicx/stats_ops.hpp"

namespace hicx::diff {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kCoefficientLimit = 50.0;
constexpr double kRidge = 1e-12;
constexpr int kMaxIterations = 100;

// In-place Cholesky solve of the p x p system a x = b; false when a is not
// positive definite to within a relative 1e-10 of its diagonal.
bool cholesky_solve(std::vector<double>& a, std::vector<double>& b, std::size_t p) {
    double scale = 0.0;
    for (std::size_t i = 0; i < p; ++i) {
        scale = std::max(scale, std::abs(a[i * p + i]));
    }
    for (std::size_t j = 0; j < p; ++j) {
        double d = a[j * p + j];
        for (std::size_t k = 0; k < j; ++k) {
            d -= a[j * p + k] * a[j * p + k];
        }
        if (!(d > 1e-10 * scale)) {
            return false;
        }
        d = std::sqrt(d);
        a[j * p + j] = d;
        for (std::size_t i = j + 1; i < p; ++i) {
            double s = a[i * p + j];
            for (std::size_t k = 0; k < j; ++k) {
                s -= a[i * p + k] * a[j * p + k];
            }
            a[i * p + j] = s / d;
        }
    }
    for (std::size_t i = 0; i < p; ++i) {
        double s = b[i];
        for (std::size_t k = 0; k < i; ++k) {
            s -= a[i * p + k] * b[k];
        }
        b[i] = s / a[i * p + i];
    }
    for (std::size_t ii = p; ii-- > 0;) {
        double s = b[ii];
        for (std::size_t k = ii + 1; k < p; ++k) {
            s -= a[k * p + ii] * b[k];
        }
        b[ii] = s / a[ii * p + ii];
    }
    return true;
}

void fitted_means(std::span<const double> log_offset, std::span<const double> x,
                  std::size_t columns, const std::vector<double>& beta, std::vector<double>& eta,
                  std::vector<double>& mu) {
    const std::size_t n = log_offset.size();
    for (std::size_t s = 0; s < n; ++s) {
        double linear = 0.0;
        for (std::size_t c = 0; c < columns; ++c) {
            linear += x[s * columns + c] * beta[c];
        }
        eta[s] = linear;
        mu[s] = std::max(std::exp(std::clamp(log_offset[s] + linear, -700.0, 700.0)), 1e-300);
    }
}

// Pearson moment estimate: the phi with sum (y - mu)^2 / (mu (1 + phi mu)) = df.
double moment_dispersion(std::span<const double> y, const std::vector<double>& mu, int df) {
    const auto pearson = [&](double phi) {
        double sum = 0.0;
        for (std::size_t s = 0; s < y.size(); ++s) {
            const double r = y[s] - mu[s];
            sum += r * r / (mu[s] * (1.0 + phi * mu[s]));
        }
        return sum;
    };
    if (pearson(0.0) <= df) {
        return 0.0;
    }
    double lo = 0.0;
    double hi = 1e3;
    for (int i = 0; i < 80; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (pearson(mid) > df) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

double sorted_median(const std::vector<double>& sorted) {
    const std::size_t n = sorted.size();
    if (n == 0) {
        return kNaN;
    }
    return n % 2 == 1 ? sorted[n / 2] : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
}

}  // namespace

// --------------------------------------------------------------------------

double trigamma(double x) {
    if (!(x > 0.0)) {
        return kNaN;
    }
    double acc = 0.0;
    while (x < 6.0) {
        acc += 1.0 / (x * x);
        x += 1.0;
    }
    const double z = 1.0 / x;
    const double z2 = z * z;
    return acc + z + z2 / 2.0 +
           z * z2 * (1.0 / 6.0 - z2 * (1.0 / 30.0 - z2 * (1.0 / 42.0 - z2 * (1.0 / 30.0 - z2 * 5.0 / 66.0))));
}

double tetragamma(double x) {
    if (!(x > 0.0)) {
        return kNaN;
    }
    double acc = 0.0;
    while (x < 6.0) {
        acc -= 2.0 / (x * x * x);
        x += 1.0;
    }
    const double z = 1.0 / x;
    const double z2 = z * z;
    return acc - z2 - z * z2 - z2 * z2 / 2.0 +
           z2 * z2 * z2 * (1.0 / 6.0 - z2 * (1.0 / 6.0 - z2 * (3.0 / 10.0 - z2 * 5.0 / 6.0)));
}

double trigamma_inverse(double y) {
    if (std::isnan(y) || y < 0.0) {
        return kNaN;
    }
    if (y == 0.0) {
        return kInf;
    }
    if (y > 1e7) {
        return 1.0 / std::sqrt(y);
    }
    if (y < 1e-6) {
        return 1.0 / y;
    }
    double x = 0.5 + 1.0 / y;
    for (int i = 0; i < 50; ++i) {
        const double tri = trigamma(x);
        const double dif = tri * (1.0 - tri / y) / tetragamma(x);
        x += dif;
        if (-dif / x < 1e-8) {
            break;
        }
    }
    return x;
}

double student_t_sf(double t, double df) {
    if (std::isnan(t) || std::isnan(df) || !(df > 0.0)) {
        return kNaN;
    }
    if (df > 1e6) {
        return hicx::stats::normal_sf(t);
    }
    if (std::isinf(t)) {
        return t > 0.0 ? 0.0 : 1.0;
    }
    const double tail = 0.5 * hicx::stats::betainc(0.5 * df, 0.5, df / (df + t * t));
    return t >= 0.0 ? tail : 1.0 - tail;
}

// --------------------------------------------------------------------------

DesignMatrices build_design(const Design& design) {
    DesignMatrices out;
    const std::size_t n = design.condition.size();
    out.samples = n;
    std::size_t in_a = 0;
    std::size_t in_b = 0;
    for (int c : design.condition) {
        if (c == 0) {
            ++in_a;
        } else if (c == 1) {
            ++in_b;
        } else {
            throw std::runtime_error("a condition label must be 0 or 1");
        }
    }
    if (in_a == 0 || in_b == 0) {
        throw std::runtime_error("each condition needs at least one sample");
    }
    out.exploratory = in_a < 2 || in_b < 2;

    std::vector<std::string> levels;
    if (!design.block.empty()) {
        if (design.block.size() != n) {
            throw std::runtime_error("--blocks needs one label per sample (" + std::to_string(n) +
                                     "), got " + std::to_string(design.block.size()));
        }
        for (const std::string& label : design.block) {
            if (std::find(levels.begin(), levels.end(), label) == levels.end()) {
                levels.push_back(label);
            }
        }
    }
    out.null_columns = 1 + (levels.empty() ? 0 : levels.size() - 1);
    const std::size_t p = out.null_columns;
    out.null.assign(n * p, 0.0);
    out.full.assign(n * (p + 1), 0.0);
    out.condition.assign(n, 0.0);
    for (std::size_t s = 0; s < n; ++s) {
        out.null[s * p] = 1.0;
        out.full[s * (p + 1)] = 1.0;
        if (!levels.empty()) {
            const std::size_t level = static_cast<std::size_t>(
                std::find(levels.begin(), levels.end(), design.block[s]) - levels.begin());
            if (level > 0) {
                out.null[s * p + level] = 1.0;
                out.full[s * (p + 1) + level] = 1.0;
            }
        }
        out.condition[s] = static_cast<double>(design.condition[s]);
        out.full[s * (p + 1) + p] = out.condition[s];
    }
    const auto full_rank = [&](const std::vector<double>& x, std::size_t columns) {
        std::vector<double> xtx(columns * columns, 0.0);
        for (std::size_t s = 0; s < n; ++s) {
            for (std::size_t i = 0; i < columns; ++i) {
                for (std::size_t j = 0; j < columns; ++j) {
                    xtx[i * columns + j] += x[s * columns + i] * x[s * columns + j];
                }
            }
        }
        std::vector<double> rhs(columns, 0.0);
        return cholesky_solve(xtx, rhs, columns);
    };
    if (!full_rank(out.full, p + 1)) {
        throw std::runtime_error("the blocks are confounded with the conditions: the condition "
                                 "effect cannot be estimated");
    }
    out.full_df = static_cast<int>(n) - static_cast<int>(p + 1);
    out.null_df = static_cast<int>(n) - static_cast<int>(p);
    if (!out.exploratory && out.full_df < 1) {
        throw std::runtime_error("the design leaves no residual degrees of freedom to estimate the "
                                 "dispersion from; use fewer blocks");
    }
    if (out.exploratory && out.null_df < 1) {
        throw std::runtime_error("with a single sample in a condition the design leaves no "
                                 "degrees of freedom to estimate the dispersion from");
    }
    return out;
}

// --------------------------------------------------------------------------

double nb_deviance(std::span<const double> y, std::span<const double> mu, double phi) {
    double sum = 0.0;
    for (std::size_t s = 0; s < y.size(); ++s) {
        const double m = std::max(mu[s], 1e-300);
        const double ylog = y[s] > 0.0 ? y[s] * std::log(y[s] / m) : 0.0;
        double d = 0.0;
        if (phi > 0.0) {
            d = 2.0 * (ylog - (y[s] + 1.0 / phi) * (std::log1p(phi * y[s]) - std::log1p(phi * m)));
        } else {
            d = 2.0 * (ylog - (y[s] - m));
        }
        sum += std::max(d, 0.0);
    }
    return sum;
}

GlmFit fit_nb_glm(std::span<const double> y, std::span<const double> log_offset,
                  std::span<const double> x, std::size_t columns, double phi) {
    const std::size_t n = y.size();
    const std::size_t p = columns;
    GlmFit fit;
    fit.beta.assign(p, 0.0);

    // Start from the least squares fit of log(y + 0.5) - offset.
    {
        std::vector<double> xtx(p * p, 0.0);
        std::vector<double> xtz(p, 0.0);
        for (std::size_t s = 0; s < n; ++s) {
            const double z = std::log(y[s] + 0.5) - log_offset[s];
            for (std::size_t i = 0; i < p; ++i) {
                xtz[i] += x[s * p + i] * z;
                for (std::size_t j = 0; j < p; ++j) {
                    xtx[i * p + j] += x[s * p + i] * x[s * p + j];
                }
            }
        }
        for (std::size_t i = 0; i < p; ++i) {
            xtx[i * p + i] += kRidge;
        }
        if (cholesky_solve(xtx, xtz, p)) {
            for (std::size_t i = 0; i < p; ++i) {
                fit.beta[i] = std::clamp(xtz[i], -kCoefficientLimit, kCoefficientLimit);
            }
        }
    }

    std::vector<double> eta(n);
    std::vector<double> mu(n);
    std::vector<double> candidate(p);
    std::vector<double> eta_c(n);
    std::vector<double> mu_c(n);
    fitted_means(log_offset, x, p, fit.beta, eta, mu);
    double deviance = nb_deviance(y, mu, phi);

    for (int iteration = 1; iteration <= kMaxIterations; ++iteration) {
        std::vector<double> xtwx(p * p, 0.0);
        std::vector<double> xtwz(p, 0.0);
        for (std::size_t s = 0; s < n; ++s) {
            const double w = mu[s] / (1.0 + phi * mu[s]);
            const double z = eta[s] + (y[s] - mu[s]) / mu[s];
            for (std::size_t i = 0; i < p; ++i) {
                xtwz[i] += x[s * p + i] * w * z;
                for (std::size_t j = 0; j < p; ++j) {
                    xtwx[i * p + j] += x[s * p + i] * w * x[s * p + j];
                }
            }
        }
        for (std::size_t i = 0; i < p; ++i) {
            xtwx[i * p + i] += kRidge;
        }
        if (!cholesky_solve(xtwx, xtwz, p)) {
            break;
        }
        for (std::size_t i = 0; i < p; ++i) {
            candidate[i] = std::clamp(xtwz[i], -kCoefficientLimit, kCoefficientLimit);
        }
        fitted_means(log_offset, x, p, candidate, eta_c, mu_c);
        double deviance_c = nb_deviance(y, mu_c, phi);
        for (int halving = 0; halving < 30 && !(deviance_c <= deviance * (1.0 + 1e-12) + 1e-12);
             ++halving) {
            for (std::size_t i = 0; i < p; ++i) {
                candidate[i] = 0.5 * (candidate[i] + fit.beta[i]);
            }
            fitted_means(log_offset, x, p, candidate, eta_c, mu_c);
            deviance_c = nb_deviance(y, mu_c, phi);
        }
        double step = 0.0;
        for (std::size_t i = 0; i < p; ++i) {
            step = std::max(step, std::abs(candidate[i] - fit.beta[i]));
        }
        const double previous = deviance;
        fit.beta = candidate;
        eta.swap(eta_c);
        mu.swap(mu_c);
        deviance = deviance_c;
        fit.iterations = iteration;
        if (step < 1e-10 || std::abs(previous - deviance) <= 1e-12 * (std::abs(deviance) + 0.1)) {
            fit.converged = true;
            break;
        }
    }
    fit.deviance = deviance;
    return fit;
}

// --------------------------------------------------------------------------

std::vector<double> binned_trend(std::span<const double> covariate, std::span<const double> values,
                                 bool trimmed) {
    const std::size_t n = covariate.size();
    std::vector<double> out(n, kNaN);
    if (n == 0) {
        return out;
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return covariate[a] < covariate[b];
    });
    const std::size_t bins = std::max<std::size_t>(1, std::min<std::size_t>(20, n / 50));
    std::vector<double> centre(bins);
    std::vector<double> level(bins);
    for (std::size_t b = 0; b < bins; ++b) {
        const std::size_t begin = b * n / bins;
        const std::size_t end = (b + 1) * n / bins;
        std::vector<double> xs;
        std::vector<double> ys;
        for (std::size_t k = begin; k < end; ++k) {
            xs.push_back(covariate[order[k]]);
            ys.push_back(values[order[k]]);
        }
        centre[b] = sorted_median(xs);
        std::sort(ys.begin(), ys.end());
        std::size_t first = 0;
        std::size_t last = ys.size();
        if (trimmed) {
            const std::size_t cut = ys.size() / 10;
            if (ys.size() > 2 * cut) {
                first = cut;
                last = ys.size() - cut;
            }
        }
        double sum = 0.0;
        for (std::size_t k = first; k < last; ++k) {
            sum += ys[k];
        }
        level[b] = sum / static_cast<double>(last - first);
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double v = covariate[i];
        if (v <= centre.front()) {
            out[i] = level.front();
        } else if (v >= centre.back()) {
            out[i] = level.back();
        } else {
            const std::size_t j = static_cast<std::size_t>(
                std::upper_bound(centre.begin(), centre.end(), v) - centre.begin());
            const std::size_t i0 = j - 1;
            const double span = centre[j] - centre[i0];
            out[i] = span > 0.0 ? level[i0] + (level[j] - level[i0]) * (v - centre[i0]) / span
                                : level[i0];
        }
    }
    return out;
}

double simes(std::span<const double> pvalues) {
    std::vector<double> p;
    for (double v : pvalues) {
        if (!std::isnan(v)) {
            p.push_back(v);
        }
    }
    if (p.empty()) {
        return kNaN;
    }
    std::sort(p.begin(), p.end());
    const double m = static_cast<double>(p.size());
    double best = 1.0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        best = std::min(best, p[k] * m / static_cast<double>(k + 1));
    }
    return best;
}

// --------------------------------------------------------------------------

double quantile_sorted(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) {
        return kNaN;
    }
    const double position = q * static_cast<double>(sorted.size() - 1);
    const auto below = static_cast<std::size_t>(std::floor(position));
    const std::size_t above = std::min(below + 1, sorted.size() - 1);
    const double fraction = position - static_cast<double>(below);
    return sorted[below] + (sorted[above] - sorted[below]) * fraction;
}

namespace {

// log of the density of z = log F, F ~ F(d1, d2).
double log_f_log_density(double z, double d1, double d2) {
    const double lbeta = hicx::stats::gammaln(0.5 * d1) + hicx::stats::gammaln(0.5 * d2) -
                         hicx::stats::gammaln(0.5 * (d1 + d2));
    return 0.5 * d1 * std::log(d1 / d2) + 0.5 * d1 * z -
           0.5 * (d1 + d2) * std::log1p(d1 * std::exp(z) / d2) - lbeta;
}

// P(log F <= z).
double log_f_cdf(double z, double d1, double d2) {
    const double x = d1 * std::exp(z);
    return hicx::stats::betainc(0.5 * d1, 0.5 * d2, x / (x + d2));
}

double log_f_quantile(double q, double d1, double d2) {
    double lo = -200.0;
    double hi = 200.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (log_f_cdf(mid, d1, d2) < q) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace

TrimmedMoments trimmed_log_f_moments(double d1, double d2, double lower, double upper) {
    TrimmedMoments moments;
    const double a = log_f_quantile(lower, d1, d2);
    const double b = log_f_quantile(upper, d1, d2);
    // Composite Simpson over [a, b], normalised by the integral itself so
    // that the result is the mean and variance of z given a <= z <= b.
    const int intervals = 4000;
    const double h = (b - a) / intervals;
    double mass = 0.0;
    double first = 0.0;
    double second = 0.0;
    for (int i = 0; i <= intervals; ++i) {
        const double z = a + h * i;
        const double w = (i == 0 || i == intervals) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        const double p = std::exp(log_f_log_density(z, d1, d2));
        mass += w * p;
        first += w * p * z;
        second += w * p * z * z;
    }
    moments.mean = first / mass;
    moments.variance = second / mass - moments.mean * moments.mean;
    return moments;
}

RobustPrior robust_prior(const std::vector<double>& s2, double df,
                         std::span<const double> covariate) {
    const std::size_t n = s2.size();
    RobustPrior out;
    out.prior_s2.assign(n, kNaN);
    std::vector<double> log_s2(n);
    for (std::size_t k = 0; k < n; ++k) {
        log_s2[k] = std::log(std::max(s2[k], 1e-12));
    }
    const std::vector<double> trend = binned_trend(covariate, log_s2, true);
    std::vector<double> residual(n);
    for (std::size_t k = 0; k < n; ++k) {
        residual[k] = log_s2[k] - trend[k];
    }
    std::vector<double> sorted = residual;
    std::sort(sorted.begin(), sorted.end());
    const double lo = quantile_sorted(sorted, kPriorTrimLower);
    const double hi = quantile_sorted(sorted, kPriorTrimUpper);
    double count = 0.0;
    double sum = 0.0;
    double squares = 0.0;
    for (double r : sorted) {
        if (r >= lo && r <= hi) {
            count += 1.0;
            sum += r;
            squares += r * r;
        }
    }
    const double observed_mean = sum / count;
    const double observed_variance = squares / count - observed_mean * observed_mean;

    // The trimmed variance of log F(df, d0) falls as d0 grows; beyond
    // kPriorDfLimit the prior is taken as exact (d0 infinite).
    const auto variance_at = [&](double d0) {
        return trimmed_log_f_moments(df, d0, kPriorTrimLower, kPriorTrimUpper).variance;
    };
    double d0 = kInf;
    double reference_mean = 0.0;
    if (observed_variance > variance_at(kPriorDfLimit)) {
        double log_lo = std::log(kPriorDfFloor);
        double log_hi = std::log(kPriorDfLimit);
        if (observed_variance >= variance_at(kPriorDfFloor)) {
            log_hi = log_lo;
        } else {
            for (int i = 0; i < 60; ++i) {
                const double mid = 0.5 * (log_lo + log_hi);
                if (variance_at(std::exp(mid)) > observed_variance) {
                    log_lo = mid;
                } else {
                    log_hi = mid;
                }
            }
        }
        d0 = std::exp(0.5 * (log_lo + log_hi));
        reference_mean = trimmed_log_f_moments(df, d0, kPriorTrimLower, kPriorTrimUpper).mean;
    } else {
        reference_mean = trimmed_log_f_moments(df, kPriorDfLimit, kPriorTrimLower, kPriorTrimUpper).mean;
    }
    out.prior_df = d0;
    out.observed_trimmed_variance = observed_variance;
    for (std::size_t k = 0; k < n; ++k) {
        out.prior_s2[k] = std::exp(trend[k] + observed_mean - reference_mean);
    }
    return out;
}

// --------------------------------------------------------------------------

FamilyResult test_family(const Family& family, const DesignMatrices& design,
                         const FamilyOptions& options) {
    const std::size_t samples = family.samples;
    if (samples != design.samples) {
        throw std::runtime_error("the family and the design disagree on the number of samples");
    }
    const std::size_t units = family.units();
    FamilyResult result;
    result.log_fold.assign(units, kNaN);
    result.pvalue.assign(units, kNaN);
    result.nb_dispersion.assign(units, kNaN);
    result.ql_dispersion.assign(units, kNaN);
    result.prior_df = kNaN;
    result.test_df = kNaN;

    std::vector<std::size_t> index;
    std::vector<double> covariate;
    for (std::size_t u = 0; u < units; ++u) {
        double total = 0.0;
        bool ok = true;
        for (std::size_t s = 0; s < samples; ++s) {
            const double y = family.counts[u * samples + s];
            const double o = family.log_offsets[u * samples + s];
            if (!std::isfinite(y) || y < 0.0 || !std::isfinite(o)) {
                ok = false;
                break;
            }
            total += y;
        }
        if (!ok || !(total > 0.0)) {
            continue;
        }
        const double c = family.covariate.empty()
                             ? std::log(total / static_cast<double>(samples) + 0.5)
                             : family.covariate[u];
        if (!std::isfinite(c)) {
            continue;
        }
        index.push_back(u);
        covariate.push_back(c);
    }
    const std::size_t n = index.size();
    if (n < kMinimumFamilyUnits) {
        return result;
    }

    const bool exploratory = design.exploratory;
    const std::span<const double> x_disp = exploratory ? std::span<const double>(design.null)
                                                       : std::span<const double>(design.full);
    const std::size_t p_disp = exploratory ? design.null_columns : design.null_columns + 1;
    const int df = exploratory ? design.null_df : design.full_df;
    const unsigned threads = std::max(1U, options.threads);
    const auto y_of = [&](std::size_t k) {
        return std::span<const double>(family.counts).subspan(index[k] * samples, samples);
    };
    const auto o_of = [&](std::size_t k) {
        return std::span<const double>(family.log_offsets).subspan(index[k] * samples, samples);
    };

    // 1. NB dispersion trend from per-unit moment estimates of the Poisson fit.
    std::vector<double> moment(n);
    hicx::parallel_for(n, threads, [&](std::size_t k) {
        const GlmFit fit = fit_nb_glm(y_of(k), o_of(k), x_disp, p_disp, 0.0);
        std::vector<double> eta(samples);
        std::vector<double> mu(samples);
        fitted_means(o_of(k), x_disp, p_disp, fit.beta, eta, mu);
        moment[k] = moment_dispersion(y_of(k), mu, df);
    });
    std::vector<double> phi = binned_trend(covariate, moment, true);
    for (double& v : phi) {
        v = std::max(v, 1e-8);
    }

    // 2. Quasi-likelihood dispersion with empirical Bayes moderation.
    std::vector<double> s2(n);
    hicx::parallel_for(n, threads, [&](std::size_t k) {
        s2[k] = fit_nb_glm(y_of(k), o_of(k), x_disp, p_disp, phi[k]).deviance / df;
    });
    const RobustPrior fitted = robust_prior(s2, static_cast<double>(df), covariate);
    const double d0 = fitted.prior_df;
    const std::vector<double>& prior = fitted.prior_s2;
    std::vector<double> posterior(n);
    for (std::size_t k = 0; k < n; ++k) {
        posterior[k] = std::isinf(d0) ? prior[k] : (d0 * prior[k] + df * s2[k]) / (d0 + df);
    }
    double test_df = std::min(d0 + df, static_cast<double>(df) * static_cast<double>(n));
    if (exploratory) {
        // The unit's own residual contains the tested contrast.
        posterior = prior;
        test_df = d0;
    }

    // 3. The TREAT likelihood ratio test of the condition coefficient.
    const double tau = std::abs(options.min_log_fold);
    const std::size_t p_full = design.null_columns + 1;
    hicx::parallel_for(n, threads, [&](std::size_t k) {
        const std::span<const double> y = y_of(k);
        const std::span<const double> o = o_of(k);
        const GlmFit full = fit_nb_glm(y, o, design.full, p_full, phi[k]);
        const double beta = full.beta.back();
        std::vector<double> shifted(samples);
        const auto root = [&](double c) {
            for (std::size_t s = 0; s < samples; ++s) {
                shifted[s] = o[s] + c * design.condition[s];
            }
            const GlmFit constrained =
                fit_nb_glm(y, shifted, design.null, design.null_columns, phi[k]);
            return std::sqrt(std::max(constrained.deviance - full.deviance, 0.0) / posterior[k]);
        };
        double p = kNaN;
        if (tau == 0.0) {
            p = 2.0 * student_t_sf(root(0.0), test_df);
        } else {
            const double at_plus = root(tau);
            const double at_minus = root(-tau);
            double near = beta >= 0.0 ? at_plus : at_minus;
            const double far = beta >= 0.0 ? at_minus : at_plus;
            if (std::abs(beta) < tau) {
                near = -near;
            }
            p = student_t_sf(near, test_df) + student_t_sf(far, test_df);
        }
        const std::size_t u = index[k];
        result.log_fold[u] = beta;
        result.pvalue[u] = std::min(p, 1.0);
        result.nb_dispersion[u] = phi[k];
        result.ql_dispersion[u] = posterior[k];
    });
    result.prior_df = d0;
    result.test_df = test_df;
    result.tested = n;
    return result;
}

}  // namespace hicx::diff
