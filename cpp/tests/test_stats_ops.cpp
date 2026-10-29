// Unit tests for hicx::stats and hicx::simd.
//
// The expected values of the special functions were evaluated by mpmath at 50
// digits (cpp/scripts/special_function_reference.py), independently of
// Boost.Math, and rounded to float64; the tests bound the error in units of
// the last place. The remaining expected values were produced by the reference
// oracle, scipy 1.14.1 and numpy in ~/miniconda3/envs/__hicexplorer@3.7.6, on
// the inputs written here. The inputs are integer recurrences divided by a
// power of two, so both sides see exactly the same bits with no parsing
// question.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "hicx/numpy_compat.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/stats_ops.hpp"

namespace {

std::vector<double> ramp(std::size_t n, std::int64_t multiplier, std::int64_t modulus,
                         double divisor, double offset = 0.0) {
    std::vector<double> values(n);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = static_cast<double>((static_cast<std::int64_t>(i) * multiplier) %
                                        modulus) /
                        divisor +
                    offset;
    }
    return values;
}

}  // namespace

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// |actual - expected| <= ulps * eps * |expected|
bool within(double actual, double expected, double ulps) {
    if (expected == 0.0) {
        return actual == 0.0;
    }
    return std::fabs(actual - expected) <= ulps * kEps * std::fabs(expected);
}

struct Point1 {
    double x;
    double expected;
};
struct Point2 {
    double a;
    double b;
    double expected;
};
struct Point3 {
    double a;
    double b;
    double x;
    double expected;
};

}  // namespace

TEST_CASE("ndtr, erf and erfc match 50 digit references") {
    const Point1 ndtr[] = {
        {-37.0, 5.725571222524577e-300}, {-20.0, 2.7536241186062337e-89},
        {-6.0, 9.86587645037698e-10},    {-2.5, 0.006209665325776135},
        {-1.0, 0.15865525393145705},     {-0.5, 0.3085375387259869},
        {0.0, 0.5},                      {0.3, 0.6179114221889527},
        {0.7071, 0.7602478320101363},    {1.0, 0.8413447460685429},
        {2.5, 0.9937903346742238},       {6.0, 0.9999999990134123},
        {9.0, 1.0}};
    for (const Point1& p : ndtr) {
        // The argument x / sqrt(2) is rounded, which moves a tail value by about
        // x * x units in the last place.
        CHECK_MESSAGE(within(hicx::stats::ndtr(p.x), p.expected, 4.0 + p.x * p.x),
                      "ndtr(" << p.x << ")");
    }
    // norm.sf is ndtr(-x), the route the rank sum test takes.
    CHECK(hicx::stats::normal_sf(2.5) == hicx::stats::ndtr(-2.5));
    CHECK(within(hicx::stats::normal_sf(37.0), 5.725571222524577e-300, 4.0 + 37.0 * 37.0));

    const Point1 erf[] = {{0.1, 0.1124629160182849},  {0.5, 0.5204998778130465},
                          {1.0, 0.8427007929497149},  {2.0, 0.9953222650189527},
                          {5.0, 0.9999999999984626},  {-1.5, -0.9661051464753108}};
    for (const Point1& p : erf) {
        CHECK_MESSAGE(within(hicx::stats::erf(p.x), p.expected, 4.0), "erf(" << p.x << ")");
    }
    const Point1 erfc[] = {{0.1, 0.887537083981715},      {1.0, 0.15729920705028513},
                           {2.0, 0.004677734981047266},   {5.0, 1.537459794428035e-12},
                           {10.0, 2.088487583762545e-45}, {-1.5, 1.9661051464753108}};
    for (const Point1& p : erfc) {
        CHECK_MESSAGE(within(hicx::stats::erfc(p.x), p.expected, 4.0), "erfc(" << p.x << ")");
    }
}

TEST_CASE("ndtr, erf and erfc edges") {
    CHECK(std::isnan(hicx::stats::ndtr(kNaN)));
    CHECK(std::isnan(hicx::stats::erf(kNaN)));
    CHECK(std::isnan(hicx::stats::erfc(kNaN)));
    CHECK(hicx::stats::ndtr(kInf) == 1.0);
    CHECK(hicx::stats::ndtr(-kInf) == 0.0);
    CHECK(hicx::stats::ndtr(0.0) == 0.5);
    CHECK(hicx::stats::ndtr(-38.0) > 0.0);
    CHECK(hicx::stats::ndtr(-40.0) == 0.0);
    CHECK(hicx::stats::ndtr(-1000.0) == 0.0);
    CHECK(hicx::stats::erf(kInf) == 1.0);
    CHECK(hicx::stats::erf(-kInf) == -1.0);
    CHECK(hicx::stats::erfc(kInf) == 0.0);
    CHECK(hicx::stats::erfc(-kInf) == 2.0);
}

TEST_CASE("digamma matches 50 digit references and has the edges of psi") {
    const Point1 points[] = {{1.0, -0.5772156649015329},   {1.5, 0.03648997397857652},
                             {2.0, 0.42278433509846713},   {0.25, -4.2274535333762655},
                             {7.5, 1.9467574842460869},    {100.0, 4.600161852738087},
                             {-0.5, 0.03648997397857652},  {-2.5, 1.103156640645243},
                             {1e-05, -100000.57719921567}, {1000000.0, 13.815510057964191}};
    for (const Point1& p : points) {
        CHECK_MESSAGE(within(hicx::stats::digamma(p.x), p.expected, 4.0), "digamma(" << p.x << ")");
    }
    CHECK(std::isnan(hicx::stats::digamma(kNaN)));
    CHECK(hicx::stats::digamma(kInf) == kInf);
    CHECK(std::isnan(hicx::stats::digamma(-kInf)));
    CHECK(hicx::stats::digamma(0.0) == -kInf);
    CHECK(hicx::stats::digamma(-0.0) == kInf);
    CHECK(std::isnan(hicx::stats::digamma(-1.0)));
    CHECK(std::isnan(hicx::stats::digamma(-2.0)));
}

TEST_CASE("gammaln matches 50 digit references and has the edges of lgamma") {
    const Point1 points[] = {{0.5, 0.5723649429247001},   {7.25, 7.0521854507385395},
                             {100.5, 361.4355404677776},  {100000.0, 1051287.7089736569},
                             {-1.5, 0.860047015376481},   {1e-10, 23.025850929882736}};
    for (const Point1& p : points) {
        CHECK_MESSAGE(within(hicx::stats::gammaln(p.x), p.expected, 4.0), "gammaln(" << p.x << ")");
    }
    CHECK(hicx::stats::gammaln(1.0) == 0.0);
    CHECK(hicx::stats::gammaln(2.0) == 0.0);
    CHECK(std::isnan(hicx::stats::gammaln(kNaN)));
    CHECK(hicx::stats::gammaln(0.0) == kInf);
    CHECK(hicx::stats::gammaln(-0.0) == kInf);
    CHECK(hicx::stats::gammaln(-1.0) == kInf);
    CHECK(hicx::stats::gammaln(-7.0) == kInf);
    CHECK(hicx::stats::gammaln(kInf) == kInf);
    CHECK(hicx::stats::gammaln(-kInf) == -kInf);
    CHECK(hicx::stats::gammaln(1e308) == kInf);
}

TEST_CASE("gamma_q and gamma_p_inverse match 50 digit references") {
    const Point2 q[] = {{0.5, 0.5, 0.3173105078629141},
                        {0.5, 1.920729410347062, 0.05000000000000006},
                        {1.0, 3.0, 0.049787068367863944},
                        {2.5, 0.1, 0.9991138612111875},
                        {10.0, 5.0, 0.9681719426937951},
                        {10.0, 50.0, 1.2596084591660908e-12},
                        {0.5, 200.0, 5.5072482372124675e-89},
                        {100.0, 100.0, 0.48670120172085135}};
    for (const Point2& p : q) {
        // 1500 units of the last place: the upper tail at (0.5, 200) is 5.5e-89
        // and the evaluation is accurate to 1.5e-13 relative there.
        CHECK_MESSAGE(within(hicx::stats::gamma_q(p.a, p.b), p.expected, 1500.0),
                      "gamma_q(" << p.a << ", " << p.b << ")");
    }
    const Point2 inverse[] = {{0.5, 0.5, 0.2274682115597864},
                              {0.5, 0.95, 1.9207294103470622},
                              {0.5, 0.999, 5.4137830853313655},
                              {2.5, 1e-10, 0.00016167785731248467},
                              {10.0, 0.5, 9.668714614714132},
                              {100.0, 0.01, 78.21598305379582},
                              {1.0, 1e-20, 1e-20},
                              {0.5, 0.9999999999999999, 34.38162610583421}};
    for (const Point2& p : inverse) {
        CHECK_MESSAGE(within(hicx::stats::gamma_p_inverse(p.a, p.b), p.expected, 16.0),
                      "gamma_p_inverse(" << p.a << ", " << p.b << ")");
    }
    // The chi-squared quantile and survival function the differential test uses:
    // chi2.ppf(0.95, 1) = 3.841458820694124 and its survival function is 0.05.
    CHECK(within(2.0 * hicx::stats::gamma_p_inverse(0.5, 0.95), 3.841458820694124, 16.0));
    CHECK(within(hicx::stats::gamma_q(0.5, 3.841458820694124 / 2.0), 0.05, 1e3));
}

TEST_CASE("gamma_q and gamma_p_inverse edges") {
    CHECK(std::isnan(hicx::stats::gamma_q(kNaN, 1.0)));
    CHECK(std::isnan(hicx::stats::gamma_q(1.0, kNaN)));
    CHECK(std::isnan(hicx::stats::gamma_q(-1.0, 1.0)));
    CHECK(std::isnan(hicx::stats::gamma_q(1.0, -1.0)));
    CHECK(hicx::stats::gamma_q(0.0, 1.0) == 0.0);
    CHECK(std::isnan(hicx::stats::gamma_q(0.0, 0.0)));
    CHECK(hicx::stats::gamma_q(2.0, 0.0) == 1.0);
    CHECK(hicx::stats::gamma_q(2.0, kInf) == 0.0);
    CHECK(hicx::stats::gamma_q(kInf, 2.0) == 1.0);
    CHECK(std::isnan(hicx::stats::gamma_q(kInf, kInf)));
    CHECK(hicx::stats::gamma_q(0.5, 1e5) == 0.0);

    CHECK(std::isnan(hicx::stats::gamma_p_inverse(kNaN, 0.5)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(0.5, kNaN)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(0.5, -0.1)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(0.5, 1.1)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(-1.0, 0.5)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(0.0, 0.5)));
    CHECK(std::isnan(hicx::stats::gamma_p_inverse(kInf, 0.5)));
    CHECK(hicx::stats::gamma_p_inverse(0.5, 0.0) == 0.0);
    CHECK(hicx::stats::gamma_p_inverse(0.5, 1.0) == kInf);
    CHECK(hicx::stats::gamma_p_inverse(0.0, 0.0) == 0.0);
    CHECK(hicx::stats::gamma_p_inverse(0.0, 1.0) == kInf);
}

TEST_CASE("betainc matches 50 digit references and has the edges of the incomplete beta function") {
    const Point3 points[] = {{0.5, 0.5, 0.3, 0.36901011956554536},
                             {1.0, 2.0, 0.25, 0.4375},
                             {2.5, 3.5, 0.5, 0.669765272631355},
                             {10.0, 1.0, 0.9, 0.3486784401000001},
                             {0.1, 0.2, 0.05, 0.5095391215346399},
                             {200.0, 3.0, 0.995, 0.9183119488697338},
                             {0.001, 0.001, 0.5, 0.5},
                             {5.0, 5.0, 0.999, 0.9999999999998744},
                             {0.15883824665748106, 1.0, 0.0003361237454510874, 0.2807212068279708},
                             {100.0, 0.5, 0.2, 7.976289553984952e-72},
                             {3.0, 700.0, 0.01, 0.971359075703961},
                             {700.0, 3.0, 0.99, 0.02864092429603887},
                             {2.5, 4.0, 0.3, 0.3521975859067672}};
    for (const Point3& p : points) {
        // 64 units of the last place: the deep tail (100, 0.5, 0.2) = 8e-72 is
        // accurate to 1e-14 relative.
        CHECK_MESSAGE(within(hicx::stats::betainc(p.a, p.b, p.x), p.expected, 64.0),
                      "betainc(" << p.a << ", " << p.b << ", " << p.x << ")");
    }
    CHECK(std::isnan(hicx::stats::betainc(kNaN, 1.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(1.0, kNaN, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(1.0, 1.0, kNaN)));
    CHECK(hicx::stats::betainc(2.0, 3.0, 0.0) == 0.0);
    CHECK(hicx::stats::betainc(2.0, 3.0, 1.0) == 1.0);
    CHECK(std::isnan(hicx::stats::betainc(0.0, 3.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, -1.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, 3.0, 1.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, 3.0, -0.5)));
}

TEST_CASE("rankdata averages ties, as scipy does") {
    const std::vector<double> values{3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0};
    const std::vector<double> expected{4.0, 1.5, 5.0, 1.5, 6.0, 8.0, 3.0, 7.0};
    CHECK(hicx::stats::rankdata_average(values) == expected);

    const std::vector<double> flat{2.0, 2.0, 2.0, 2.0};
    const std::vector<double> flat_expected{2.5, 2.5, 2.5, 2.5};
    CHECK(hicx::stats::rankdata_average(flat) == flat_expected);

    CHECK(hicx::stats::rankdata_average(std::vector<double>{}).empty());
}

TEST_CASE("ranksums matches scipy.stats.ranksums") {
    // Two samples of all zeros: every rank is tied, the statistic is exactly
    // zero and the two sided p-value is exactly one. This is the case
    // hicFindTADs hits on an empty region and guards with a try/except.
    const std::vector<double> zeros_x(20, 0.0);
    const std::vector<double> zeros_y(30, 0.0);
    const hicx::stats::RanksumsResult flat =
        hicx::stats::ranksums(zeros_x, zeros_y);
    CHECK(flat.statistic == 0.0);
    CHECK(flat.pvalue == 1.0);

    // x[i] = ((i * 37) % 101) / 8, y[i] = ((i * 53) % 97) / 8 + 0.25.
    const std::vector<double> x = ramp(57, 37, 101, 8.0);
    const std::vector<double> y = ramp(83, 53, 97, 8.0, 0.25);
    const hicx::stats::RanksumsResult result = hicx::stats::ranksums(x, y);
    CHECK(result.statistic == doctest::Approx(-0.044534224258321584).epsilon(1e-15));
    CHECK(result.pvalue == doctest::Approx(0.9644785720158997).epsilon(1e-15));
}

TEST_CASE("the Benjamini-Hochberg cutoff is hicFindTADs' step up rule") {
    // The largest p with p <= q * rank / n. With q = 0.5 and five p-values
    // the thresholds are 0.1, 0.2, 0.3, 0.4 and 0.5, and the sorted p-values
    // are 0.01, 0.02, 0.35, 0.7, 0.9, so only the first two pass and the
    // cutoff is 0.02. Note that the rule stops at the largest *passing* value
    // rather than at the largest index that passes, which is the ordinary
    // step up procedure; hicFindTADs implements it this way and the port
    // follows it.
    const std::vector<double> pvalues{0.01, 0.02, 0.9, 0.35, 0.7};
    CHECK(hicx::stats::benjamini_hochberg_cutoff(pvalues, 0.5) == 0.02);
    // A value that only passes at a later rank is still taken.
    CHECK(hicx::stats::benjamini_hochberg_cutoff({0.05, 0.2, 0.25}, 1.0) == 0.25);
    // Nothing passes: the cutoff stays at zero, which rejects every boundary.
    CHECK(hicx::stats::benjamini_hochberg_cutoff({0.9, 0.95}, 0.01) == 0.0);
    CHECK(hicx::stats::benjamini_hochberg_cutoff({}, 0.5) == 0.0);
}

TEST_CASE("adjusted p-values skip NaN, keep ties equal and cap at one") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // m = 4: 0.01*4/1 = 0.04, 0.02*4/2 = 0.04, 0.03*4/3 = 0.04, 0.5*4/4 = 0.5
    const std::vector<double> p{0.5, nan, 0.02, 0.01, 0.03};
    const std::vector<double> bh = hicx::stats::benjamini_hochberg_adjusted(p);
    CHECK(bh[0] == 0.5);
    CHECK(std::isnan(bh[1]));
    CHECK(bh[2] == 0.02 * 4.0 / 2.0);
    CHECK(bh[3] == std::min(0.01 * 4.0 / 1.0, 0.02 * 4.0 / 2.0));
    CHECK(bh[4] == 0.03 * 4.0 / 3.0);
    // A later smaller ratio pulls the earlier ranks down.
    const std::vector<double> step_in{0.04, 0.045};
    const std::vector<double> step = hicx::stats::benjamini_hochberg_adjusted(step_in);
    CHECK(step[0] == 0.045);
    CHECK(step[1] == 0.045);
    // Tied p-values get the same adjusted value.
    const std::vector<double> tied_in{0.2, 0.2, 0.9};
    const std::vector<double> tied = hicx::stats::benjamini_hochberg_adjusted(tied_in);
    CHECK(tied[0] == tied[1]);
    const std::vector<double> capped_in{0.6, 0.9};
    CHECK(hicx::stats::benjamini_hochberg_adjusted(capped_in)[0] == 0.9);
    CHECK(hicx::stats::benjamini_hochberg_adjusted(std::vector<double>{}).empty());
    const std::vector<double> bonferroni_in{0.01, nan, 0.4};
    const std::vector<double> bonferroni = hicx::stats::bonferroni_adjusted(bonferroni_in);
    CHECK(bonferroni[0] == 0.01 * 2.0);
    CHECK(std::isnan(bonferroni[1]));
    CHECK(bonferroni[2] == 0.8);
    const std::vector<double> clamp_in{0.7, 0.9};
    CHECK(hicx::stats::bonferroni_adjusted(clamp_in)[1] == 1.0);
}

TEST_CASE("the Bonferroni scaling clamps at one and leaves NaN alone") {
    std::vector<double> pvalues{0.01, 0.2, 0.5,
                                std::numeric_limits<double>::quiet_NaN()};
    hicx::stats::bonferroni_in_place(pvalues);
    CHECK(pvalues[0] == doctest::Approx(0.04).epsilon(1e-15));
    CHECK(pvalues[1] == doctest::Approx(0.8).epsilon(1e-15));
    CHECK(pvalues[2] == 1.0);
    CHECK(std::isnan(pvalues[3]));
}

TEST_CASE("fit_nbinom matches the fit_nbinom package on the same data") {
    // The optimiser is not scipy's L-BFGS-B, so the agreement is measured
    // rather than exact; see cpp/core/src/lbfgsb.cpp. These tolerances are the
    // measured agreement, not a wish.
    const std::vector<double> small{1.0, 2.0, 3.0, 4.0,  5.0, 2.0, 3.0, 4.0,
                                    10.0, 0.0, 1.0, 7.0, 3.0, 2.0, 5.0};
    const hicx::stats::NBinomFit fit = hicx::stats::fit_nbinom(small);
    CHECK(fit.size == doctest::Approx(4.8614566250173565).epsilon(1e-5));
    CHECK(fit.prob == doctest::Approx(0.5837395334809616).epsilon(1e-5));

    std::vector<double> counts(400);
    for (std::size_t i = 0; i < counts.size(); ++i) {
        const std::int64_t k = static_cast<std::int64_t>(i);
        counts[i] = static_cast<double>((k * k * 7 + 3 * k) % 23);
    }
    const hicx::stats::NBinomFit wide = hicx::stats::fit_nbinom(counts);
    CHECK(wide.size == doctest::Approx(2.196162498991828).epsilon(1e-4));
    CHECK(wide.prob == doctest::Approx(0.18018029927276732).epsilon(1e-4));
}

TEST_CASE("the SIMD reduction is bit identical to the scalar reference") {
    // Every length that crosses a structural boundary of numpy's reduction:
    // below the eight way unroll, inside a 128 element block, across the
    // recursive split, and across the 8192 element ufunc buffer.
    std::vector<std::size_t> sizes;
    for (std::size_t n = 0; n <= 300; ++n) {
        sizes.push_back(n);
    }
    for (const std::size_t n : {1000u, 8191u, 8192u, 8193u, 16384u, 20000u, 65537u}) {
        sizes.push_back(n);
    }

    const std::vector<double> data = ramp(70000, 37, 1009, 64.0);
    for (const std::size_t n : sizes) {
        const double scalar = hicx::simd::pairwise_sum_scalar(data.data(), n);
        // The scalar kernel here must also agree with the one hicInfo and the
        // writers already use, or the two would drift apart.
        CHECK(scalar == hicx::npy::pairwise_sum(data.data(), n));
        CHECK(hicx::simd::pairwise_sum(data.data(), n) == scalar);
        if (hicx::simd::avx2_available()) {
            CHECK(hicx::simd::pairwise_sum_avx2(data.data(), n) == scalar);
        }
    }

    // The dispatcher picked the kernel the CPU actually supports, once.
    CHECK(std::string(hicx::simd::active_kernel()) ==
          (hicx::simd::avx2_available() ? "avx2" : "scalar"));

    // And against numpy directly on three of those lengths.
    CHECK(hicx::simd::pairwise_sum(data.data(), 129) == 973.484375);
    CHECK(hicx::simd::pairwise_sum(data.data(), 8192) == 64462.859375);
    CHECK(hicx::simd::pairwise_sum(data.data(), 20000) == 157456.875);
}

TEST_CASE("the SIMD reduction handles the awkward values too") {
    std::vector<double> data(1000, 0.0);
    for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] = static_cast<double>(i) * 1e-12 - 0.5;
    }
    data[17] = std::numeric_limits<double>::quiet_NaN();
    const double scalar = hicx::simd::pairwise_sum_scalar(data.data(), data.size());
    CHECK(std::isnan(scalar));
    CHECK(std::isnan(hicx::simd::pairwise_sum(data.data(), data.size())));

    data[17] = std::numeric_limits<double>::infinity();
    CHECK(hicx::simd::pairwise_sum(data.data(), data.size()) ==
          hicx::simd::pairwise_sum_scalar(data.data(), data.size()));
}

TEST_CASE("fit_cut_intervals snaps only when the bins are irregular") {
    std::vector<hicx::CutInterval> regular;
    for (int i = 0; i < 20; ++i) {
        regular.push_back(hicx::CutInterval{"chr1", i * 100, (i + 1) * 100, 1.0, ""});
    }
    CHECK(hicx::fit_cut_intervals(regular) == regular);

    // More than one percent of the bins deviate from the 100 bp median, so
    // every start and end is snapped to the nearest multiple of it.
    std::vector<hicx::CutInterval> irregular = regular;
    for (std::size_t i = 5; i < 15; ++i) {
        irregular[i].end += 37;
    }
    const std::vector<hicx::CutInterval> fitted = hicx::fit_cut_intervals(irregular);
    for (const hicx::CutInterval& interval : fitted) {
        CHECK(interval.start % 100 == 0);
        CHECK(interval.end % 100 == 0);
    }
    // A tie snaps downwards, because np.argmin takes the first minimum.
    std::vector<hicx::CutInterval> tied = irregular;
    tied[3].start = 350;
    CHECK(hicx::fit_cut_intervals(tied)[3].start == 300);
}

TEST_CASE("enlarge_bins closes gaps and starts every chromosome at zero") {
    // The doctest of hicexplorer.utilities.enlarge_bins.
    std::vector<hicx::CutInterval> intervals{
        hicx::CutInterval{"chr1", 10, 50, 1.0, ""},
        hicx::CutInterval{"chr1", 50, 80, 2.0, ""},
        hicx::CutInterval{"chr2", 10, 60, 3.0, ""},
        hicx::CutInterval{"chr2", 70, 90, 4.0, ""}};
    hicx::enlarge_bins(intervals);
    CHECK(intervals[0].start == 0);
    CHECK(intervals[0].end == 50);
    CHECK(intervals[1].start == 50);
    CHECK(intervals[1].end == 80);
    CHECK(intervals[2].start == 0);
    CHECK(intervals[2].end == 65);
    CHECK(intervals[3].start == 65);
    CHECK(intervals[3].end == 90);
}
