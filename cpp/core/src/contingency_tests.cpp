// See hicx/contingency_tests.hpp.

#include "hicx/contingency_tests.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

namespace hicx::contingency {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Double precision throughout, and no exceptions for arguments the callers have
// already validated.
using Policy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error>,
    boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::ignore_error>,
    boost::math::policies::promote_double<false>>;

using Hypergeometric = boost::math::hypergeometric_distribution<double, Policy>;
using ChiSquared = boost::math::chi_squared_distribution<double, Policy>;

constexpr std::int64_t kMaxCount = 4294967295LL;  // Boost takes unsigned arguments

double clip01(double x) {
    if (std::isnan(x)) {
        return x;
    }
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

// The table given by its margins: `successes` is the first row sum, `draws` the
// first column sum, `total` the grand total.
struct Margins {
    std::int64_t total;
    std::int64_t successes;
    std::int64_t draws;
    std::int64_t low;   // smallest k with positive probability
    std::int64_t high;  // largest

    Hypergeometric distribution() const {
        return Hypergeometric(static_cast<unsigned>(successes), static_cast<unsigned>(draws),
                              static_cast<unsigned>(total));
    }
    double pmf(std::int64_t k) const {
        if (k < low || k > high) {
            return 0.0;
        }
        return clip01(boost::math::pdf(distribution(), static_cast<unsigned>(k)));
    }
    // P(X <= k)
    double cdf(std::int64_t k) const {
        if (k < low) {
            return 0.0;
        }
        if (k >= high) {
            return 1.0;
        }
        return clip01(boost::math::cdf(distribution(), static_cast<unsigned>(k)));
    }
    // P(X > k)
    double sf(std::int64_t k) const {
        if (k < low) {
            return 1.0;
        }
        if (k >= high) {
            return 0.0;
        }
        return clip01(boost::math::cdf(boost::math::complement(distribution(), static_cast<unsigned>(k))));
    }
};

Margins make_margins(std::int64_t total, std::int64_t successes, std::int64_t draws) {
    if (total > kMaxCount) {
        throw std::out_of_range("hicx::contingency: a table total above 2^32 - 1 is not supported");
    }
    Margins m;
    m.total = total;
    m.successes = successes;
    m.draws = draws;
    const std::int64_t low = draws - (total - successes);
    m.low = low > 0 ? low : 0;
    m.high = successes < draws ? successes : draws;
    return m;
}

}  // namespace

double hypergeometric_pmf(std::int64_t k, std::int64_t total, std::int64_t successes,
                          std::int64_t draws) {
    if (total <= 0 || successes < 0 || draws < 0 || successes > total || draws > total) {
        return kNaN;
    }
    return make_margins(total, successes, draws).pmf(k);
}

std::optional<FisherExact> fisher_exact_2x2(std::int64_t a, std::int64_t b, std::int64_t c,
                                            std::int64_t d) {
    if (a < 0 || b < 0 || c < 0 || d < 0) {
        return std::nullopt;
    }
    if (a + c == 0 || b + d == 0 || a + b == 0 || c + d == 0) {
        return FisherExact{1.0, true};
    }
    const std::int64_t row1 = a + b;
    const std::int64_t column1 = a + c;
    const std::int64_t total = row1 + c + d;
    const Margins m = make_margins(total, row1, column1);

    // The mode of the hypergeometric distribution, floor((n + 1)(K + 1) / (N + 2)).
    const std::int64_t mode = static_cast<std::int64_t>(
        (static_cast<__int128>(column1 + 1) * static_cast<__int128>(row1 + 1)) / (total + 2));
    const double observed = m.pmf(a);
    const double threshold = observed * (1.0 + kTieTolerance);

    // The pmf rises up to the mode and falls after it, so the tables that do not
    // exceed the threshold are a prefix [low, i1] and a suffix [i2, high].
    if (m.pmf(mode) <= threshold) {
        return FisherExact{1.0, true};
    }
    double pvalue = 0.0;
    {
        // largest k in [low, mode - 1] with pmf(k) <= threshold, if any
        std::int64_t first = m.low;
        std::int64_t last = mode - 1;
        std::int64_t found = m.low - 1;
        while (first <= last) {
            const std::int64_t middle = first + (last - first) / 2;
            if (m.pmf(middle) <= threshold) {
                found = middle;
                first = middle + 1;
            } else {
                last = middle - 1;
            }
        }
        if (found >= m.low) {
            pvalue += m.cdf(found);
        }
    }
    {
        // smallest k in [mode + 1, high] with pmf(k) <= threshold, if any
        std::int64_t first = mode + 1;
        std::int64_t last = m.high;
        std::int64_t found = m.high + 1;
        while (first <= last) {
            const std::int64_t middle = first + (last - first) / 2;
            if (m.pmf(middle) <= threshold) {
                found = middle;
                last = middle - 1;
            } else {
                first = middle + 1;
            }
        }
        if (found <= m.high) {
            pvalue += m.sf(found - 1);
        }
    }
    if (pvalue > 1.0) {
        return FisherExact{1.0, true};
    }
    return FisherExact{pvalue, false};
}

std::optional<Chi2Result> chi2_2x2(double a, double b, double c, double d, bool yates_correction) {
    if (a < 0 || b < 0 || c < 0 || d < 0) {
        return std::nullopt;
    }
    const double observed[4] = {a, b, c, d};
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
    double terms[4];
    for (int i = 0; i < 4; ++i) {
        double difference = observed[i] - expected[i];
        if (yates_correction) {
            const double magnitude = std::fabs(difference);
            const double shift = magnitude < 0.5 ? magnitude : 0.5;
            difference = difference > 0 ? difference - shift : difference + shift;
        }
        terms[i] = difference * difference / expected[i];
    }
    const double statistic = ((terms[0] + terms[1]) + terms[2]) + terms[3];
    return Chi2Result{statistic, chi2_sf(statistic, 1.0)};
}

double chi2_quantile(double q, double df) {
    if (!(df > 0.0) || std::isnan(q)) {
        return kNaN;
    }
    if (q == 0.0) {
        return 0.0;
    }
    if (q == 1.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (q > 0.0 && q < 1.0) {
        return boost::math::quantile(ChiSquared(df), q);
    }
    return kNaN;
}

double chi2_sf(double x, double df) {
    if (!(df > 0.0) || std::isnan(x)) {
        return kNaN;
    }
    if (x < 0.0) {
        return 1.0;
    }
    if (std::isinf(x)) {
        return 0.0;
    }
    return boost::math::cdf(boost::math::complement(ChiSquared(df), x));
}

}  // namespace hicx::contingency
