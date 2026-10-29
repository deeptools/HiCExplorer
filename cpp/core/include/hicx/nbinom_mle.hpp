// Maximum likelihood fit of the negative binomial distribution, solved
// directly instead of by a general optimiser.
//
// The model is NB(size r, prob p) with the mean r (1 - p) / p of the
// parametrisation fit_nbinom/fit_nbinom.py uses, that is
//
//     P(X = x) = Gamma(x + r) / (Gamma(r) x!) p^r (1 - p)^x.
//
// For a sample x_1 .. x_n with mean m the log likelihood is maximised over p
// in closed form, p = r / (r + m), and the remaining profile score
//
//     g(r) = sum_i [psi(x_i + r) - psi(r)] - n log(1 + m / r)
//
// (psi the digamma function) has exactly one root when the sample variance
// (divisor n) exceeds the mean: g is positive for small r and tends to zero
// from below as -n (v - m) / (2 r^2) for large r. The root is bracketed
// geometrically from the moment estimator and refined with Boost.Math's
// toms748_solve to full double precision.
//
// Where no finite maximiser exists or the objective of the Python package is
// not defined, the answer is the one its optimiser leaves, as far as that is a
// well defined value (see fit_nbinom_mle).
//
// Nothing here is derived from SciPy or from any optimiser source.

#ifndef HICX_NBINOM_MLE_HPP
#define HICX_NBINOM_MLE_HPP

#include <array>
#include <optional>
#include <span>

#include "hicx/stats_ops.hpp"

namespace hicx::stats {

// The size reported where the maximum likelihood estimate is the Poisson limit
// (sample variance not above the sample mean, mean above zero). The negative
// binomial with this size differs from the Poisson distribution of the same
// mean by a relative variance of mean / kPoissonLimitSize.
inline constexpr double kPoissonLimitSize = 1.0e6;

// fit_nbinom.fit(X) for a float64 (or integer) X.
//
//   * empty X: size 10 and prob NaN, status 2 (the moment start of an empty
//     array, as the Python package yields).
//   * a value that is negative, not finite, or whose factorial overflows a
//     double (above 170): the objective of the Python package is not finite
//     anywhere, and its optimiser returns the starting point. The same starting
//     point is returned here: the moment estimator, or `initial` where given,
//     clamped into the bounds size >= eps, eps <= prob <= 1; status 0.
//   * all values equal to zero: size 10, prob 1 (any size is a maximiser; this
//     is the start of the moment estimator).
//   * variance not above the mean, mean above zero: the likelihood increases
//     towards size = infinity, prob = 1 (the Poisson distribution). The size
//     is kPoissonLimitSize and prob = size / (size + mean); status 1. The
//     optimiser of the Python package stops at an arbitrary point of this
//     ridge (a size of order 10^2 for the cases measured), which is not a
//     value that can be derived.
//   * otherwise: the unique root of the profile score, status 0; iterations is
//     the number of score evaluations.
//
// `initial` (size, prob) replaces the moment start in the second case only; the
// fitted value does not depend on a starting point.
[[nodiscard]] NBinomFit fit_nbinom_mle(
    std::span<const double> data,
    std::optional<std::array<double, 2>> initial = std::nullopt);

// The profile score g(r) above for the sample summarised by its distinct
// values and their multiplicities. Exposed for the unit tests.
[[nodiscard]] double nbinom_profile_score(std::span<const double> distinct_values,
                                          std::span<const double> multiplicities, double mean,
                                          double size);

}  // namespace hicx::stats

#endif  // HICX_NBINOM_MLE_HPP
