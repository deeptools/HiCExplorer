// Statistical tests as scipy 1.14.1 computes them, for chicDifferentialTest.
//
// scipy.stats.fisher_exact on a 2x2 table (two-sided), scipy.stats.
// chi2_contingency on a 2x2 table without Yates' correction, and
// scipy.stats.chi2.ppf. Each is written against the Python it replaces,
// operation by operation, over the routines scipy/special/functions.json maps
// the special functions to at v1.14.1:
//
//   hypergeom.pmf, cdf, sf   Boost.Math hypergeometric_distribution under
//                            scipy's StatsPolicy (boost_special_functions.h,
//                            hypergeom_*_wrap), inside rv_discrete's support
//                            checks and np.clip(., 0, 1);
//   special.chdtrc           cephes_chdtrc, the Cephes igamc of scipy's C++
//                            translation (scipy/special/special/cephes);
//   special.gammaincinv      cephes_igami, the same translation.
//
// The Cephes headers are scipy's own, downloaded at the v1.14.1 commit
// (cpp/core/CMakeLists.txt), and Boost.Math is the commit scipy vendors, so
// the arithmetic is the same code; the sums follow numpy's order for a 2x2
// float64 array (((a + b) + c) + d over the whole array, a + b and a + c along
// an axis), measured.

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

// scipy.special.chdtrc(df, x) and scipy.special.gammaincinv(a, p).
[[nodiscard]] double chdtrc(double df, double x);
[[nodiscard]] double gammaincinv(double a, double p);

}  // namespace hicx::scipy

#endif  // HICX_SCIPY_STATS_HPP
