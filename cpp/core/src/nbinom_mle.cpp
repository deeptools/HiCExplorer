// See hicx/nbinom_mle.hpp.

#include "hicx/nbinom_mle.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/tools/roots.hpp>

#include "hicx/numpy_compat.hpp"

namespace hicx::stats {

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

using ValuePolicy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::pole_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::ignore_error>,
    boost::math::policies::promote_double<false>>;

// psi(x + r) - psi(r) for x >= 0. For a small integer x this is the exact
// finite sum of 1 / (r + k), k < x, which keeps its relative accuracy when r
// is large and the difference of the two digamma values would cancel.
double digamma_increment(double x, double r) {
    if (x == 0.0) {
        return 0.0;
    }
    if (x <= 64.0 && x == std::floor(x)) {
        double sum = 0.0;
        for (double k = static_cast<double>(static_cast<int>(x)) - 1.0; k >= 0.0; k -= 1.0) {
            sum += 1.0 / (r + k);
        }
        return sum;
    }
    return boost::math::digamma(x + r, ValuePolicy()) - boost::math::digamma(r, ValuePolicy());
}

NBinomFit clamped_start(double size, double prob) {
    NBinomFit fit;
    // A NaN stays NaN: the comparisons below are false for it.
    fit.size = size < kEps ? kEps : size;
    fit.prob = prob < kEps ? kEps : (prob > 1.0 ? 1.0 : prob);
    return fit;
}

}  // namespace

double nbinom_profile_score(std::span<const double> distinct_values,
                            std::span<const double> multiplicities, double mean, double size) {
    double n = 0.0;
    double sum = 0.0;
    for (std::size_t j = 0; j < distinct_values.size(); ++j) {
        sum += multiplicities[j] * digamma_increment(distinct_values[j], size);
        n += multiplicities[j];
    }
    return sum - n * std::log1p(mean / size);
}

NBinomFit fit_nbinom_mle(std::span<const double> data,
                         std::optional<std::array<double, 2>> initial) {
    NBinomFit fit;
    const std::size_t n = data.size();
    if (n == 0) {
        // The moment start of an empty array is (10, NaN); a start given by the
        // caller is returned as it is, clamped into the bounds.
        fit = initial ? clamped_start((*initial)[0], (*initial)[1]) : NBinomFit{};
        if (!initial) {
            fit.size = 10.0;
            fit.prob = kNaN;
        }
        fit.status = 2;
        return fit;
    }

    // The moment estimator, sums in numpy's pairwise order as the package's
    // initial values are.
    const double count = static_cast<double>(n);
    const double mean = npy::pairwise_sum(data.data(), n) / count;
    std::vector<double> scratch(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double difference = data[i] - mean;
        scratch[i] = difference * difference;
    }
    const double variance = npy::pairwise_sum(scratch.data(), n) / count;
    const double moment_size = variance > mean ? (mean * mean) / (variance - mean) : 10.0;
    const double moment_prob =
        moment_size + mean != 0.0 ? moment_size / (moment_size + mean) : moment_size;

    // Values for which the package's objective is not finite: its optimiser
    // never moves.
    bool defined = true;
    for (const double x : data) {
        if (!std::isfinite(x) || x < 0.0 || std::isinf(std::tgamma(x + 1.0))) {
            defined = false;
            break;
        }
    }
    if (!defined) {
        if (initial) {
            return clamped_start((*initial)[0], (*initial)[1]);
        }
        return clamped_start(moment_size, moment_prob);
    }

    if (!(variance > mean)) {
        if (mean == 0.0) {
            // Every value is zero: the likelihood is flat in size at prob = 1.
            fit.size = initial ? std::max((*initial)[0], kEps) : 10.0;
            fit.prob = 1.0;
            return fit;
        }
        fit.size = kPoissonLimitSize;
        fit.prob = fit.size / (fit.size + mean);
        fit.status = 1;
        return fit;
    }

    // Distinct values with multiplicities.
    std::vector<double> sorted(data.begin(), data.end());
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> values;
    std::vector<double> weights;
    for (const double x : sorted) {
        if (!values.empty() && values.back() == x) {
            weights.back() += 1.0;
        } else {
            values.push_back(x);
            weights.push_back(1.0);
        }
    }

    int evaluations = 0;
    const auto score = [&](double r) {
        ++evaluations;
        return nbinom_profile_score(values, weights, mean, r);
    };

    // Bracket the root geometrically from the moment estimator: g > 0 below
    // the root and g < 0 above it.
    double r0 = std::isfinite(moment_size) && moment_size > 0.0 ? moment_size : 1.0;
    double g0 = score(r0);
    double lo = r0;
    double hi = r0;
    double g_lo = g0;
    double g_hi = g0;
    bool bracketed = false;
    if (g0 == 0.0) {
        fit.size = r0;
        fit.prob = r0 / (r0 + mean);
        fit.iterations = evaluations;
        return fit;
    }
    if (g0 > 0.0) {
        for (int step = 0; step < 2000; ++step) {
            lo = hi;
            g_lo = g_hi;
            hi *= 2.0;
            g_hi = score(hi);
            if (g_hi <= 0.0) {
                bracketed = true;
                break;
            }
            if (hi > 1.0e18) {
                break;
            }
        }
    } else {
        for (int step = 0; step < 2000; ++step) {
            hi = lo;
            g_hi = g_lo;
            lo *= 0.5;
            g_lo = score(lo);
            if (g_lo >= 0.0) {
                bracketed = true;
                break;
            }
            if (lo < 1.0e-300) {
                break;
            }
        }
    }
    if (!bracketed) {
        // The score does not change sign in double precision: variance and mean
        // are indistinguishable. Report the Poisson limit.
        fit.size = kPoissonLimitSize;
        fit.prob = fit.size / (fit.size + mean);
        fit.status = 1;
        fit.iterations = evaluations;
        return fit;
    }

    std::uintmax_t max_iterations = 500;
    const auto bracket = boost::math::tools::toms748_solve(
        score, lo, hi, g_lo, g_hi,
        boost::math::tools::eps_tolerance<double>(std::numeric_limits<double>::digits),
        max_iterations);
    const double r = 0.5 * (bracket.first + bracket.second);
    fit.size = r;
    fit.prob = r / (r + mean);
    fit.status = 0;
    fit.iterations = evaluations;
    return fit;
}

}  // namespace hicx::stats
