// Tests on 2x2 contingency tables for chicDifferentialTest: Fisher's exact
// test, the chi-squared test of independence, and the chi-squared quantile and
// survival function.
//
// Both tests are evaluated from their textbook definitions with Boost.Math:
//
//   Fisher's exact test, two-sided: the hypergeometric distribution of the
//   upper left cell given the margins, and the p-value is the total probability
//   of all tables with the same margins whose probability does not exceed that
//   of the observed table.
//
//   Chi-squared test: expected counts row sum * column sum / total, the
//   statistic sum (O - E)^2 / E over the four cells, one degree of freedom, the
//   p-value from the chi-squared survival function. Optionally with Yates'
//   continuity correction, |O - E| reduced by min(0.5, |O - E|).
//
// Nothing here is derived from SciPy; the results agree with the SciPy
// functions of the same purpose to the precision of the arithmetic.

#ifndef HICX_CONTINGENCY_TESTS_HPP
#define HICX_CONTINGENCY_TESTS_HPP

#include <cstdint>
#include <optional>

namespace hicx::contingency {

// Two tables are tied when the probability of one is at most the probability of
// the other times (1 + kTieTolerance). The probabilities of tables that are
// equal in exact arithmetic (a table and its mirror image under symmetric
// margins) differ in the last digits of a floating point evaluation by up to
// about 1e-13 relative for the table sizes of Hi-C counts; the tolerance lies
// well above that and below any difference between distinct probabilities that
// occurs in practice.
inline constexpr double kTieTolerance = 1.0e-12;

// P(X = k) of the hypergeometric distribution for a population of `total`
// with `successes` marked items and `draws` draws (0 outside the support).
[[nodiscard]] double hypergeometric_pmf(std::int64_t k, std::int64_t total, std::int64_t successes,
                                        std::int64_t draws);

struct FisherExact {
    double pvalue = 0.0;
    // The p-value is the constant 1.0 (an empty row or column, an observed
    // table at the mode, or a total above 1), not a computed sum. The Python
    // binding returns such a value as a Python float instead of numpy.float64.
    bool constant_one = false;
};

// Two-sided p-value of [[a, b], [c, d]]; nullopt for a negative cell. A table
// with an empty row or column has p = 1.
[[nodiscard]] std::optional<FisherExact> fisher_exact_2x2(std::int64_t a, std::int64_t b,
                                                          std::int64_t c, std::int64_t d);

struct Chi2Result {
    double statistic = 0.0;
    double pvalue = 0.0;
};

// Chi-squared test of independence of [[a, b], [c, d]] with one degree of
// freedom; nullopt for a negative cell or an expected count of zero (an empty
// row or column), where the statistic is undefined.
[[nodiscard]] std::optional<Chi2Result> chi2_2x2(double a, double b, double c, double d,
                                                 bool yates_correction = false);

// Quantile (inverse cdf) of the chi-squared distribution: 0 at q = 0,
// infinity at q = 1, NaN outside [0, 1] or for df <= 0.
[[nodiscard]] double chi2_quantile(double q, double df);

// Survival function: 1 for x < 0, NaN for df <= 0 or a NaN argument.
[[nodiscard]] double chi2_sf(double x, double df);

}  // namespace hicx::contingency

#endif  // HICX_CONTINGENCY_TESTS_HPP
