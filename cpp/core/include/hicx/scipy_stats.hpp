// Statistical tests for chicDifferentialTest: a two-sided Fisher exact test on
// a 2x2 table, the chi-squared test of a 2x2 table without continuity
// correction, and the chi-squared quantile, with the results of the Python
// functions scipy.stats.fisher_exact, scipy.stats.chi2_contingency and
// scipy.stats.chi2.ppf. The special functions underneath are Boost.Math's:
//
//   hypergeometric pmf, cdf, sf   boost::math::hypergeometric_distribution,
//                                 inside the support checks of a discrete
//                                 distribution and clipped to [0, 1];
//   chi-squared sf and quantile   hicx::stats::gamma_q and gamma_p_inverse.
//
// The sums follow numpy's order for a 2x2 float64 array (((a + b) + c) + d over
// the whole array, a + b and a + c along an axis).

#ifndef HICX_SCIPY_STATS_HPP
#define HICX_SCIPY_STATS_HPP

#include <cstdint>
#include <optional>

namespace hicx::scipy {

// scipy.stats.hypergeom.pmf(k, M, n, N), .cdf and .sf for integer arguments.
[[nodiscard]] double hypergeom_pmf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N);
[[nodiscard]] double hypergeom_cdf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N);
[[nodiscard]] double hypergeom_sf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N);

// scipy.stats.fisher_exact([[a, b], [c, d]]).pvalue, two-sided; nullopt where
// scipy raises ValueError (a negative value).
[[nodiscard]] std::optional<double> fisher_exact_pvalue(std::int64_t a, std::int64_t b,
                                                        std::int64_t c, std::int64_t d);

// The same with the type scipy returns the p-value in: a Python float for the
// constant 1.0 of a table with an empty row or column, of a table at the mode
// and of min(pvalue, 1.0) when 1.0 is the smaller, numpy.float64 otherwise.
struct FisherExact {
    double pvalue = 0.0;
    bool python_float = false;
};
[[nodiscard]] std::optional<FisherExact> fisher_exact(std::int64_t a, std::int64_t b, std::int64_t c,
                                                      std::int64_t d);

struct Chi2Contingency {
    double statistic = 0.0;
    double pvalue = 0.0;
};

// scipy.stats.chi2_contingency([[a, b], [c, d]], correction=False); nullopt
// where scipy raises ValueError (a negative value, an expected frequency of 0,
// or observed and expected sums that disagree beyond sqrt(eps)).
[[nodiscard]] std::optional<Chi2Contingency> chi2_contingency_2x2(double a, double b, double c,
                                                                  double d);

// scipy.stats.chi2.ppf(q, df).
[[nodiscard]] double chi2_ppf(double q, double df);

// The chi-squared survival function of x with df degrees of freedom,
// Q(df / 2, x / 2) (hicx::stats::gamma_q): 1 for x < 0, NaN for df < 1.
[[nodiscard]] double chi2_sf(double x, double df);

}  // namespace hicx::scipy

#endif  // HICX_SCIPY_STATS_HPP
