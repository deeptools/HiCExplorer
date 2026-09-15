// See hicx/scipy_stats.hpp.

#include "hicx/scipy_stats.hpp"

#include <cmath>
#include <limits>

#include <boost/math/distributions/hypergeometric.hpp>

// scipy/special/boost_special_functions.h, StatsPolicy.
using StatsPolicy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::user_error>,
    boost::math::policies::evaluation_error<boost::math::policies::user_error>,
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>,
    boost::math::policies::discrete_quantile<boost::math::policies::integer_round_up>>;

// scipy's handlers: an evaluation error warns and returns the value, an
// overflow sets OverflowError and returns 0.
namespace boost::math::policies {
template <class RealType>
RealType user_evaluation_error(const char*, const char*, const RealType& val) {
    return val;
}
template <class RealType>
RealType user_overflow_error(const char*, const char*, const RealType&) {
    return 0;
}
}  // namespace boost::math::policies

#include "special/cephes/chdtr.h"
#include "special/cephes/igami.h"

namespace hicx::scipy {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// hypergeom_pmf_wrap, hypergeom_cdf_wrap and hypergeom_sf_wrap, called with
// (k, n, N, M) as scu._hypergeom_pmf(k, n, N, M) is.
double boost_pmf(double k, double n, double N, double M) {
    if (std::isfinite(k)) {
        return boost::math::pdf(boost::math::hypergeometric_distribution<double, StatsPolicy>(
                                    static_cast<unsigned>(n), static_cast<unsigned>(N),
                                    static_cast<unsigned>(M)),
                                k);
    }
    return kNaN;
}

double boost_cdf(double k, double n, double N, double M) {
    if (std::isfinite(k)) {
        return boost::math::cdf(boost::math::hypergeometric_distribution<double, StatsPolicy>(
                                    static_cast<unsigned>(n), static_cast<unsigned>(N),
                                    static_cast<unsigned>(M)),
                                k);
    }
    return 1 - std::signbit(k);
}

double boost_sf(double k, double n, double N, double M) {
    return boost::math::cdf(boost::math::complement(
        boost::math::hypergeometric_distribution<double, StatsPolicy>(
            static_cast<unsigned>(n), static_cast<unsigned>(N), static_cast<unsigned>(M)),
        k));
}

// np.clip(x, 0, 1), NaN propagating.
double clip01(double x) {
    if (std::isnan(x)) {
        return x;
    }
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

// hypergeom_gen._argcheck
bool argcheck(std::int64_t M, std::int64_t n, std::int64_t N) {
    return M > 0 && n >= 0 && N >= 0 && n <= M && N <= M;
}

struct Support {
    std::int64_t a;
    std::int64_t b;
};

Support support(std::int64_t M, std::int64_t n, std::int64_t N) {
    const std::int64_t low = N - (M - n);
    return {low > 0 ? low : 0, n < N ? n : N};
}

// numpy int64 arithmetic wraps.
std::int64_t wrapping_multiply(std::int64_t a, std::int64_t b) {
    return static_cast<std::int64_t>(static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b));
}

// scipy.stats._binomtest._binary_search_for_binom_tst(a, d, lo, hi)
template <typename F>
std::int64_t binary_search(F a, double d, std::int64_t lo, std::int64_t hi) {
    while (lo < hi) {
        const std::int64_t mid = lo + (hi - lo) / 2;
        const double midval = a(mid);
        if (midval < d) {
            lo = mid + 1;
        } else if (midval > d) {
            hi = mid - 1;
        } else {
            return mid;
        }
    }
    if (a(lo) <= d) {
        return lo;
    }
    return lo - 1;
}

}  // namespace

double hypergeom_pmf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N) {
    if (!argcheck(M, n, N)) {
        return kNaN;
    }
    const Support s = support(M, n, N);
    if (k >= s.a && k <= s.b) {
        return clip01(boost_pmf(static_cast<double>(k), static_cast<double>(n),
                                static_cast<double>(N), static_cast<double>(M)));
    }
    return 0.0;
}

double hypergeom_cdf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N) {
    if (!argcheck(M, n, N)) {
        return kNaN;
    }
    const Support s = support(M, n, N);
    if (k >= s.b) {
        return 1.0;
    }
    if (k >= s.a) {
        return clip01(boost_cdf(static_cast<double>(k), static_cast<double>(n),
                                static_cast<double>(N), static_cast<double>(M)));
    }
    return 0.0;
}

double hypergeom_sf(std::int64_t k, std::int64_t M, std::int64_t n, std::int64_t N) {
    if (!argcheck(M, n, N)) {
        return kNaN;
    }
    const Support s = support(M, n, N);
    if (k < s.a) {
        return 1.0;
    }
    if (k < s.b) {
        return clip01(boost_sf(static_cast<double>(k), static_cast<double>(n),
                               static_cast<double>(N), static_cast<double>(M)));
    }
    return 0.0;
}

std::optional<double> fisher_exact_pvalue(std::int64_t c00, std::int64_t c01, std::int64_t c10,
                                          std::int64_t c11) {
    if (c00 < 0 || c01 < 0 || c10 < 0 || c11 < 0) {
        return std::nullopt;
    }
    if (c00 + c10 == 0 || c01 + c11 == 0 || c00 + c01 == 0 || c10 + c11 == 0) {
        return 1.0;
    }
    const std::int64_t n1 = c00 + c01;
    const std::int64_t n2 = c10 + c11;
    const std::int64_t n = c00 + c10;
    const std::int64_t total = n1 + n2;
    const auto pmf = [&](std::int64_t x) { return hypergeom_pmf(x, total, n1, n); };

    const std::int64_t mode = static_cast<std::int64_t>(
        static_cast<double>(wrapping_multiply(n + 1, n1 + 1)) / static_cast<double>(n1 + n2 + 2));
    const double pexact = pmf(c00);
    const double pmode = pmf(mode);

    const double epsilon = 1e-14;
    const double gamma = 1 + epsilon;

    // np.abs(pexact - pmode) / np.maximum(pexact, pmode) <= epsilon
    const double maximum = (std::isnan(pexact) || std::isnan(pmode))
                               ? kNaN
                               : (pexact < pmode ? pmode : pexact);
    if (std::fabs(pexact - pmode) / maximum <= epsilon) {
        return 1.0;
    }
    double pvalue = 0.0;
    if (c00 < mode) {
        const double plower = hypergeom_cdf(c00, total, n1, n);
        if (pmf(n) > pexact * gamma) {
            return plower;
        }
        const std::int64_t guess =
            binary_search([&](std::int64_t x) { return -pmf(x); }, -pexact * gamma, mode, n);
        pvalue = plower + hypergeom_sf(guess, total, n1, n);
    } else {
        const double pupper = hypergeom_sf(c00 - 1, total, n1, n);
        if (pmf(0) > pexact * gamma) {
            return pupper;
        }
        const std::int64_t guess = binary_search(pmf, pexact * gamma, 0, mode);
        pvalue = pupper + hypergeom_cdf(guess, total, n1, n);
    }
    // min(pvalue, 1.0)
    return 1.0 < pvalue ? 1.0 : pvalue;
}

std::optional<Chi2Contingency> chi2_contingency_2x2(double a, double b, double c, double d) {
    if (a < 0 || b < 0 || c < 0 || d < 0) {
        return std::nullopt;
    }
    // expected_freq: the row sums (2, 1) times the column sums (1, 2), divided
    // by the total ** 1.
    const double row0 = a + b;
    const double row1 = c + d;
    const double column0 = a + c;
    const double column1 = b + d;
    const double total = ((a + b) + c) + d;
    const double expected[4] = {row0 * column0 / total, row0 * column1 / total,
                                row1 * column0 / total, row1 * column1 / total};
    for (const double value : expected) {
        if (value == 0.0) {
            return std::nullopt;
        }
    }
    // power_divergence(observed, expected, ddof=2, axis=None, lambda_=1)
    const double observed[4] = {a, b, c, d};
    const double observed_sum = ((observed[0] + observed[1]) + observed[2]) + observed[3];
    const double expected_sum = ((expected[0] + expected[1]) + expected[2]) + expected[3];
    const double rtol = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
    const double minimum = (std::isnan(observed_sum) || std::isnan(expected_sum))
                               ? kNaN
                               : (expected_sum < observed_sum ? expected_sum : observed_sum);
    const double relative_diff = std::fabs(observed_sum - expected_sum) / minimum;
    if (relative_diff > rtol) {
        return std::nullopt;
    }
    double terms[4];
    for (int i = 0; i < 4; ++i) {
        const double difference = observed[i] - expected[i];
        terms[i] = difference * difference / expected[i];
    }
    const double statistic = ((terms[0] + terms[1]) + terms[2]) + terms[3];
    // _SimpleChi2(df = 4 - 1 - 2).sf(statistic)
    return Chi2Contingency{statistic, chdtrc(1.0, statistic)};
}

double chi2_ppf(double q, double df) {
    // rv_continuous.ppf with loc 0 and scale 1 over the support [0, inf)
    if (!(df > 0)) {
        return kNaN;
    }
    if (q == 0.0) {
        return 0.0 * 1.0 + 0.0;
    }
    if (q == 1.0) {
        return std::numeric_limits<double>::infinity() * 1.0 + 0.0;
    }
    if (0 < q && q < 1) {
        return (2 * gammaincinv(df / 2, q)) * 1.0 + 0.0;
    }
    return kNaN;
}

double chdtrc(double df, double x) {
    return special::cephes::chdtrc(df, x);
}

double gammaincinv(double a, double p) {
    return special::cephes::igami(a, p);
}

}  // namespace hicx::scipy
