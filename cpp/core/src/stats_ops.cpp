#include "hicx/stats_ops.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "hicx/numpy_compat.hpp"

namespace hicx::stats {

namespace {

// np.finfo(np.float64).eps, the bound and the guard fit_nbinom uses.
constexpr double kMachineEpsilonStats = 2.220446049250313e-16;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// The Boost.Math policy of the special functions below. Every function is
// evaluated in double precision, never promoted to long double, so that the
// result does not depend on the width of the platform's long double. Domain,
// pole, overflow and evaluation errors return a value (NaN or infinity)
// instead of throwing; the argument ranges in which a function has a value
// are tested explicitly before the call.
using ValuePolicy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::pole_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::ignore_error>,
    boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

// The policy of the incomplete beta function: double precision, and at most
// 400 iterations of the root finders. Errors are exceptions that betainc maps
// to values.
using BetaPolicy = boost::math::policies::policy<boost::math::policies::promote_float<false>,
                                                 boost::math::policies::promote_double<false>,
                                                 boost::math::policies::max_root_iterations<400>>;

}  // namespace

double erf(double x) {
    if (std::isnan(x)) {
        return kNaN;
    }
    return boost::math::erf(x, ValuePolicy());
}

double erfc(double x) {
    if (std::isnan(x)) {
        return kNaN;
    }
    return boost::math::erfc(x, ValuePolicy());
}

// Phi(x) = erfc(-x / sqrt(2)) / 2. The complementary error function keeps the
// relative precision in the lower tail, where Phi is tiny.
double ndtr(double x) {
    if (std::isnan(x)) {
        return kNaN;
    }
    return 0.5 * erfc(-x * M_SQRT1_2);
}

double normal_sf(double x) { return ndtr(-x); }

double gammaln(double x) {
    // log|Gamma(x)|. The poles at zero and the negative integers give +inf, a
    // NaN stays NaN and an infinite argument is returned unchanged.
    if (!std::isfinite(x)) {
        return x;
    }
    if (x <= 0.0 && x == std::floor(x)) {
        return kInf;
    }
    return boost::math::lgamma(x, ValuePolicy());
}

double digamma(double x) {
    if (std::isnan(x) || x == kInf) {
        return x;
    }
    if (x == -kInf) {
        return kNaN;
    }
    if (x == 0.0) {
        // The pole at zero is approached from above for +0 and from below for -0.
        return std::copysign(kInf, -x);
    }
    if (x < 0.0 && x == std::floor(x)) {
        return kNaN;
    }
    return boost::math::digamma(x, ValuePolicy());
}

double gamma_q(double a, double x) {
    if (std::isnan(a) || std::isnan(x) || a < 0.0 || x < 0.0) {
        return kNaN;
    }
    if (a == 0.0) {
        return x > 0.0 ? 0.0 : kNaN;
    }
    if (x == 0.0) {
        return 1.0;
    }
    if (std::isinf(a)) {
        return std::isinf(x) ? kNaN : 1.0;
    }
    if (std::isinf(x)) {
        return 0.0;
    }
    return boost::math::gamma_q(a, x, ValuePolicy());
}

double gamma_p_inverse(double a, double p) {
    if (std::isnan(a) || std::isnan(p) || a < 0.0 || p < 0.0 || p > 1.0) {
        return kNaN;
    }
    if (p == 0.0) {
        return 0.0;
    }
    if (p == 1.0) {
        return kInf;
    }
    if (a == 0.0 || std::isinf(a)) {
        return kNaN;
    }
    const double x = boost::math::gamma_p_inv(a, p, ValuePolicy());
    // A root below the smallest normal number is returned as that number; as
    // a double it is 0 (for a = 1e-5 and p = 0.5 the root is 0.5 ** 1e5).
    return x <= std::numeric_limits<double>::min() ? 0.0 : x;
}

double betainc(double a, double b, double x) {
    if (std::isnan(a) || std::isnan(b) || std::isnan(x)) {
        return kNaN;
    }
    if (a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0) {
        return kNaN;
    }
    try {
        return boost::math::ibeta(a, b, x, BetaPolicy());
    } catch (const std::domain_error&) {
        return kNaN;
    } catch (const std::overflow_error&) {
        return kInf;
    } catch (const std::underflow_error&) {
        return 0.0;
    } catch (...) {
        return kNaN;
    }
}

double nbinom_sf(double x, double r, double p) {
    // 1 - cnb.cdf(x, r, p), the subtraction in float64 after the incomplete
    // beta function, as hicDetectLoops.py writes it.
    return 1.0 - betainc(r, x + 1.0, p);
}

// --------------------------------------------------------------------------

std::vector<double> rankdata_average(std::span<const double> values) {
    const std::size_t n = values.size();
    std::vector<double> ranks(n, 0.0);
    if (n == 0) {
        return ranks;
    }
    // np.argsort(kind='quicksort') is not stable, but rankdata's result does
    // not depend on the order within a tie group, so a stable sort gives the
    // same ranks and is deterministic.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&values](std::size_t a, std::size_t b) {
                         return values[a] < values[b];
                     });

    std::size_t i = 0;
    while (i < n) {
        std::size_t j = i + 1;
        while (j < n && values[order[j]] == values[order[i]]) {
            ++j;
        }
        // Ranks i+1 .. j average to (i + j + 1) / 2.
        const double rank = 0.5 * (static_cast<double>(i) + static_cast<double>(j) + 1.0);
        for (std::size_t k = i; k < j; ++k) {
            ranks[order[k]] = rank;
        }
        i = j;
    }
    return ranks;
}

RanksumsResult ranksums(std::span<const double> x, std::span<const double> y) {
    RanksumsResult result;
    const std::size_t n1 = x.size();
    const std::size_t n2 = y.size();
    if (n1 == 0 && n2 == 0) {
        result.defined = false;
        result.statistic = std::numeric_limits<double>::quiet_NaN();
        result.pvalue = std::numeric_limits<double>::quiet_NaN();
        return result;
    }

    std::vector<double> all;
    all.reserve(n1 + n2);
    all.insert(all.end(), x.begin(), x.end());
    all.insert(all.end(), y.begin(), y.end());
    const std::vector<double> ranked = rankdata_average(all);

    // np.sum over the first n1 ranks, in numpy's reduction order.
    const double s = npy::pairwise_sum(ranked.data(), n1);
    const double total = static_cast<double>(n1) + static_cast<double>(n2);
    const double expected = static_cast<double>(n1) * (total + 1.0) / 2.0;
    const double variance =
        static_cast<double>(n1) * static_cast<double>(n2) * (total + 1.0) / 12.0;
    const double z = (s - expected) / std::sqrt(variance);
    result.statistic = z;
    result.pvalue = 2.0 * normal_sf(std::abs(z));
    return result;
}

double benjamini_hochberg_cutoff(std::vector<double> pvalues, double q) {
    if (pvalues.empty()) {
        return 0.0;
    }
    std::sort(pvalues.begin(), pvalues.end());
    const double n = static_cast<double>(pvalues.size());
    double largest = 0.0;
    for (std::size_t i = 0; i < pvalues.size(); ++i) {
        const double p = pvalues[i];
        if (p <= q * (static_cast<double>(i) + 1.0) / n) {
            if (p >= largest) {
                largest = p;
            }
        }
    }
    return largest;
}

std::vector<double> benjamini_hochberg_adjusted(std::span<const double> pvalues) {
    std::vector<double> out(pvalues.begin(), pvalues.end());
    std::vector<std::size_t> order;
    order.reserve(pvalues.size());
    for (std::size_t i = 0; i < pvalues.size(); ++i) {
        if (!std::isnan(pvalues[i])) {
            order.push_back(i);
        }
    }
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return pvalues[a] < pvalues[b]; });
    const double m = static_cast<double>(order.size());
    std::vector<double> adjusted(order.size());
    for (std::size_t k = 0; k < order.size(); ++k) {
        adjusted[k] = pvalues[order[k]] * m / static_cast<double>(k + 1);
    }
    for (std::size_t k = adjusted.size(); k > 1; --k) {
        adjusted[k - 2] = std::min(adjusted[k - 2], adjusted[k - 1]);
    }
    for (std::size_t k = 0; k < order.size(); ++k) {
        out[order[k]] = std::min(adjusted[k], 1.0);
    }
    return out;
}

std::vector<double> bonferroni_adjusted(std::span<const double> pvalues) {
    std::vector<double> out(pvalues.begin(), pvalues.end());
    const double m = static_cast<double>(
        std::count_if(pvalues.begin(), pvalues.end(), [](double p) { return !std::isnan(p); }));
    for (double& p : out) {
        if (!std::isnan(p)) {
            p = std::min(p * m, 1.0);
        }
    }
    return out;
}

// --------------------------------------------------------------------------
// The float32 flavour of fit_nbinom. See the comment on NBinomPrecision for
// why this exists; the float64 flavour lives beside the optimiser in
// lbfgsb.cpp and this one is deliberately a separate function rather than a
// branch inside it, so that the tools that do not need it are untouched.

NBinomFit fit_nbinom(std::span<const double> data, NBinomPrecision precision,
                     std::optional<std::array<double, 2>> initial) {
    if (precision == NBinomPrecision::Float64) {
        return fit_nbinom(data);
    }
    NBinomFit fit;
    const std::size_t n = data.size();
    if (n == 0) {
        fit.status = 2;
        return fit;
    }

    // The array as scipy holds it. The caller has already rounded the values
    // to float32; this is the array the ufuncs actually see.
    std::vector<float> x(n, 0.0F);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = static_cast<float>(data[i]);
    }

    // np.sum(np.log(factorial(X))): scipy's factorial returns a float64 array
    // even for a float32 argument, but it evaluates gamma(x + 1) in the
    // float32 loop and widens the result, so it **overflows above about 34.6**
    // rather than above 170.6. That is not a detail: as soon as one obs/exp
    // value at a genomic distance exceeds it, this term is infinite, the whole
    // objective is infinite for every parameter pair, the forward difference
    // gradient is NaN and fmin_l_bfgs_b returns its starting point. Measured
    // on the corpus: 3 of the 199 distances of gm12878_chr1.cool are in that
    // state, and computing the term in float64 instead lets the optimiser run
    // on those three and call a loop the reference does not, which is exactly
    // how this was found.
    std::vector<double> log_factorial(n, 0.0);
    bool factorial_overflowed = false;
    for (std::size_t i = 0; i < n; ++i) {
        const float gamma =
            static_cast<float>(std::tgamma(static_cast<double>(x[i] + 1.0F)));
        if (std::isinf(gamma)) {
            factorial_overflowed = true;
            break;
        }
        log_factorial[i] = std::log(static_cast<double>(gamma));
    }

    const double log_factorial_sum =
        factorial_overflowed ? std::numeric_limits<double>::infinity()
                             : npy::pairwise_sum(log_factorial);

    // gammaln is by far the most expensive thing in the objective and it is
    // evaluated over the same array on every iteration, so the distinct values
    // are collected once and each objective evaluation looks the result up
    // instead of recomputing it. This is exact, not an approximation: gammaln
    // is a function of its argument, and the values are scattered back into
    // their original positions before the reduction, so the summation order is
    // untouched and the result is bit identical to the direct loop.
    //
    // It pays because an obs/exp value is a small integer count divided by a
    // per distance constant, so the same value recurs constantly. Measured on
    // gm12878_chr1.cool, the six sampled distance distributions hold 130,256
    // values of which 4,469 are distinct, 3.4 percent. The table is only built
    // when it at least halves the work; on a distribution of genuinely
    // distinct values, such as the 19 distances of the GSE63525 cool where the
    // ratio is 99.9 percent, it is skipped and the direct loop runs.
    std::vector<float> unique;
    std::vector<std::int32_t> slot;
    bool use_table = false;
    {
        bool finite = true;
        for (std::size_t i = 0; i < n && finite; ++i) {
            finite = std::isfinite(x[i]);
        }
        if (finite) {
            unique = x;
            std::sort(unique.begin(), unique.end());
            unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
            if (unique.size() * 2 <= n) {
                slot.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    slot[i] = static_cast<std::int32_t>(std::distance(
                        unique.begin(),
                        std::lower_bound(unique.begin(), unique.end(), x[i])));
                }
                use_table = true;
            }
        }
        if (!use_table) {
            unique.clear();
            unique.shrink_to_fit();
        }
    }
    std::vector<float> table(unique.size(), 0.0F);

    std::vector<float> scratch(n, 0.0F);
    const auto negative_log_likelihood =
        [&](std::span<const double> parameters) -> double {
        const double r = parameters[0];
        const double p = parameters[1];
        const double safe_p = p < 1.0 ? p : 1.0 - kMachineEpsilonStats;

        // gammaln(X + r): the addition runs in the float32 loop, with the
        // scalar cast to float32 first, and so does the reduction.
        const float r32 = static_cast<float>(r);
        if (use_table) {
            for (std::size_t u = 0; u < unique.size(); ++u) {
                table[u] =
                    static_cast<float>(gammaln(static_cast<double>(unique[u] + r32)));
            }
            for (std::size_t i = 0; i < n; ++i) {
                scratch[i] = table[static_cast<std::size_t>(slot[i])];
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                scratch[i] =
                    static_cast<float>(gammaln(static_cast<double>(x[i] + r32)));
            }
        }
        const float gammaln_sum = npy::pairwise_sum(scratch.data(), n);

        // X * log(1 - p), likewise float32 throughout.
        const float log1p32 = static_cast<float>(std::log(1.0 - safe_p));
        for (std::size_t i = 0; i < n; ++i) {
            scratch[i] = x[i] * log1p32;
        }
        const float tail = npy::pairwise_sum(scratch.data(), n);

        // The five terms are combined left to right, and the float32 partial
        // sums widen to float64 as they meet a float64 term.
        double value = static_cast<double>(gammaln_sum);
        value -= log_factorial_sum;
        value -= static_cast<double>(n) * gammaln(r);
        value += static_cast<double>(n) * r * std::log(p);
        value += static_cast<double>(tail);
        return -value;
    };

    // np.mean and np.var of a float32 array reduce in float32 as well, so the
    // moment estimator that seeds the optimiser is a float32 quantity.
    const float mean32 =
        npy::pairwise_sum(x.data(), n) / static_cast<float>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const float difference = x[i] - mean32;
        scratch[i] = difference * difference;
    }
    const float variance32 =
        npy::pairwise_sum(scratch.data(), n) / static_cast<float>(n);

    double start_size = 0.0;
    double start_prob = 0.0;
    if (variance32 > mean32) {
        // `size = (m ** 2) / (v - m)`. `m ** 2` with a Python int exponent
        // promotes a numpy float32 scalar to float64 in numpy 1.26, while
        // `v - m` stays float32, so the quotient and everything after it is
        // float64 even though the data is not. `m * m` would have stayed
        // float32; the difference is the eighth significant digit of the
        // starting point, which is enough to move the fit by 1e-08.
        const double mean = static_cast<double>(mean32);
        start_size = (mean * mean) / static_cast<double>(variance32 - mean32);
        const double denominator = start_size + mean;
        start_prob = denominator != 0.0 ? start_size / denominator : start_size;
    } else {
        // `size = 10` is a Python int here, and numpy 1.26 promotes
        // int + float32 to float64 for scalars, so this branch finishes in
        // double precision even though the data is float32.
        start_size = 10.0;
        const double denominator = 10.0 + static_cast<double>(mean32);
        start_prob = denominator != 0.0 ? 10.0 / denominator : 10.0;
    }

    if (initial) {
        start_size = (*initial)[0];
        start_prob = (*initial)[1];
    }
    const std::vector<double> begin{start_size, start_prob};
    const std::vector<Bound> bounds{
        Bound{kMachineEpsilonStats, std::numeric_limits<double>::infinity()},
        Bound{kMachineEpsilonStats, 1.0}};

    if (factorial_overflowed) {
        // The objective is +inf everywhere, so the forward difference gradient
        // is inf - inf = NaN, the projected gradient norm test compares
        // against NaN and std::max leaves it at zero, and minimise_lbfgsb
        // breaks out at iteration zero with the clamped starting point.
        // Returning it directly is the same answer, asserted bit for bit by
        // the unit test in cpp/tests/test_detect_loops.cpp.
        //
        // This is not an optimisation and is not claimed as one. Measured
        // interleaved on gm12878_chr1.cool, three runs each: 7.12 s of CPU
        // with it against 7.12 s without, because only 3 of the 199 distances
        // reach it. It is here because the alternative is to depend on
        // std::max(0.0, NaN) evaluating to 0.0 for the tool to terminate at
        // the right point, which is an accident of the comparison order in
        // projected_gradient_norm rather than something this code states.
        fit.size = std::min(std::max(start_size, bounds[0].lower), bounds[0].upper);
        fit.prob = std::min(std::max(start_prob, bounds[1].lower), bounds[1].upper);
        fit.status = 0;
        fit.iterations = 0;
        return fit;
    }

    const LbfgsbResult solution =
        minimise_lbfgsb(negative_log_likelihood, begin, bounds, LbfgsbOptions());
    fit.size = solution.x[0];
    fit.prob = solution.x[1];
    fit.status = solution.status;
    fit.iterations = solution.iterations;
    return fit;
}

void bonferroni_in_place(std::vector<double>& pvalues) {
    const double n = static_cast<double>(pvalues.size());
    for (double& p : pvalues) {
        p *= n;
    }
    for (double& p : pvalues) {
        if (!std::isnan(p) && p > 1.0) {
            p = 1.0;
        }
    }
}

}  // namespace hicx::stats
