// Unit tests of the hicDifferentialAnalysis count model and its random draws
// (cpp/tools/diff_engine_impl.cpp, cpp/tools/diff_contacts_impl.cpp).
//
// Reference values: scipy 1.14.1 (special.polygamma, stats.t.sf, a brentq root
// of polygamma(1, x) = y) and statsmodels 0.15.0 (GLM with the
// NegativeBinomial(alpha) and Poisson families, offsets, tol 1e-14), computed
// on 2026-09-15. The values are small hand-made inputs; the model is validated
// on real matrices by cpp/scripts/diff_calibration.py.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/diff_contacts_impl.hpp"
#include "../tools/diff_engine_impl.hpp"

namespace {

bool relative(double a, double b, double tolerance) {
    return std::abs(a - b) <= tolerance * std::max(1.0, std::abs(b));
}

}  // namespace

TEST_CASE("diff engine: trigamma, tetragamma and the trigamma inverse against scipy") {
    const double xs[] = {0.3, 1.0, 2.5, 7.0, 40.0};
    const double tri[] = {12.245364546107734, 1.6449340668482266, 0.4903577561002349,
                          0.15354517795933756, 0.025315103841291032};
    const double tetra[] = {-75.27253658872601, -2.404113806319188, -0.23620405164172736,
                            -0.023530472985855234, -0.0006408202718352985};
    for (int i = 0; i < 5; ++i) {
        CHECK(std::abs(hicx::diff::trigamma(xs[i]) / tri[i] - 1.0) < 1e-9);
        CHECK(std::abs(hicx::diff::tetragamma(xs[i]) / tetra[i] - 1.0) < 1e-9);
    }
    const double ys[] = {0.05, 0.7, 3.0, 25.0};
    const double inverse[] = {20.49583523915461, 1.874194785205666, 0.6754781052813699,
                              0.20523750101934557};
    for (int i = 0; i < 4; ++i) {
        CHECK(std::abs(hicx::diff::trigamma_inverse(ys[i]) / inverse[i] - 1.0) < 1e-7);
    }
}

TEST_CASE("diff engine: Student t survival function against scipy") {
    CHECK(relative(hicx::diff::student_t_sf(0.5, 3.0), 0.3257239824240755, 1e-12));
    CHECK(relative(hicx::diff::student_t_sf(2.1, 4.5), 0.04796890848104311, 1e-12));
    CHECK(relative(hicx::diff::student_t_sf(-1.3, 10.0), 0.8886170913965777, 1e-12));
    CHECK(relative(hicx::diff::student_t_sf(3.7, 1.0), 0.08402226282394758, 1e-12));
    CHECK(relative(hicx::diff::student_t_sf(1.96, 250.0), 0.02555331011268659, 1e-12));
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(relative(hicx::diff::student_t_sf(1.96, inf), 0.024997895148220435, 1e-12));
}

TEST_CASE("diff engine: NB GLM with offsets and a block against statsmodels") {
    hicx::diff::Design design;
    design.condition = {0, 0, 1, 1};
    design.block = {"r1", "r2", "r1", "r2"};
    const hicx::diff::DesignMatrices x = hicx::diff::build_design(design);
    REQUIRE(x.null_columns == 2);
    CHECK(x.full_df == 1);
    CHECK_FALSE(x.exploratory);
    const std::vector<double> y = {812.0, 950.0, 1405.0, 1510.0};
    const std::vector<double> offset = {std::log(1.1e3), std::log(1.25e3), std::log(1.05e3),
                                        std::log(1.2e3)};
    const hicx::diff::GlmFit fit = hicx::diff::fit_nb_glm(y, offset, x.full, 3, 0.013);
    CHECK(fit.converged);
    CHECK(relative(fit.beta[0], -0.2801805050195082, 1e-8));
    CHECK(relative(fit.beta[1], -0.016905647739565127, 1e-8));
    CHECK(relative(fit.beta[2], 0.549333876923257, 1e-8));
    CHECK(relative(fit.deviance, 0.14742704359061776, 1e-8));

    std::vector<double> shifted = offset;
    for (std::size_t s = 0; s < 4; ++s) {
        shifted[s] += 0.2 * x.condition[s];
    }
    const hicx::diff::GlmFit constrained = hicx::diff::fit_nb_glm(y, shifted, x.null, 2, 0.013);
    CHECK(relative(constrained.beta[0], -0.08581357512712005, 1e-8));
    CHECK(relative(constrained.beta[1], -0.024268936466093377, 1e-8));
    CHECK(relative(constrained.deviance, 8.879429313724245, 1e-8));
}

TEST_CASE("diff engine: Poisson GLM with a zero count against statsmodels") {
    const std::vector<double> y = {0.0, 3.0, 7.0, 2.0};
    const std::vector<double> offset(4, 0.0);
    const std::vector<double> x = {1, 0, 1, 0, 1, 1, 1, 1};
    const hicx::diff::GlmFit fit = hicx::diff::fit_nb_glm(y, offset, x, 2, 0.0);
    CHECK(relative(fit.beta[0], 0.40546510810816405, 1e-8));
    CHECK(relative(fit.beta[1], 1.0986122886681098, 1e-8));
    CHECK(relative(fit.deviance, 7.100820750400906, 1e-8));
}

TEST_CASE("diff engine: designs that cannot be estimated are refused") {
    hicx::diff::Design confounded;
    confounded.condition = {0, 0, 1, 1};
    confounded.block = {"x", "x", "y", "y"};
    CHECK_THROWS(hicx::diff::build_design(confounded));

    hicx::diff::Design single;
    single.condition = {0, 1};
    const hicx::diff::DesignMatrices x = hicx::diff::build_design(single);
    CHECK(x.exploratory);
    CHECK(x.null_df == 1);

    hicx::diff::Design no_df;
    no_df.condition = {0, 1};
    no_df.block = {"p", "q"};
    CHECK_THROWS(hicx::diff::build_design(no_df));
}

TEST_CASE("diff engine: Simes and the binned trend") {
    const std::vector<double> p = {0.04, std::nan(""), 0.01, 0.3};
    CHECK(relative(hicx::diff::simes(p), 0.03, 1e-15));
    std::vector<double> covariate;
    std::vector<double> values;
    for (int i = 0; i < 200; ++i) {
        covariate.push_back(static_cast<double>(i));
        values.push_back(2.0);
    }
    for (double v : hicx::diff::binned_trend(covariate, values, true)) {
        CHECK(v == 2.0);
    }
}

TEST_CASE("diff engine: a family's result does not depend on the thread count") {
    // Deterministic counts with replicate noise from the counter generator.
    hicx::diff::Family family;
    family.samples = 4;
    for (std::uint64_t u = 0; u < 300; ++u) {
        hicx::diffc::CounterRng rng({7, u});
        const double level = 20.0 + static_cast<double>(u % 50) * 10.0;
        for (std::size_t s = 0; s < 4; ++s) {
            const double noise = 1.0 + 0.2 * (rng.uniform() - 0.5);
            family.counts.push_back(std::round(level * noise));
            family.log_offsets.push_back(std::log(100.0 + 5.0 * static_cast<double>(s)));
        }
    }
    hicx::diff::Design design;
    design.condition = {0, 0, 1, 1};
    const hicx::diff::DesignMatrices x = hicx::diff::build_design(design);
    hicx::diff::FamilyOptions one;
    one.threads = 1;
    one.min_log_fold = std::log(1.1);
    hicx::diff::FamilyOptions many = one;
    many.threads = 7;
    const hicx::diff::FamilyResult a = hicx::diff::test_family(family, x, one);
    const hicx::diff::FamilyResult b = hicx::diff::test_family(family, x, many);
    REQUIRE(a.tested == 300);
    for (std::size_t u = 0; u < 300; ++u) {
        CHECK(std::memcmp(&a.pvalue[u], &b.pvalue[u], sizeof(double)) == 0);
        CHECK(a.pvalue[u] >= 0.0);
        CHECK(a.pvalue[u] <= 1.0);
    }
}

TEST_CASE("diff contacts: the binomial draw is keyed, exact in its moments, and complete") {
    hicx::diffc::CounterRng first({1, 2, 3});
    hicx::diffc::CounterRng second({1, 2, 3});
    CHECK(first.next() == second.next());
    hicx::diffc::CounterRng other({1, 2, 4});
    CHECK(hicx::diffc::CounterRng({1, 2, 3}).next() != other.next());

    for (const auto& [n, p] : {std::pair<std::int64_t, double>{20, 0.5},
                               std::pair<std::int64_t, double>{1000, 0.3},
                               std::pair<std::int64_t, double>{50000, 0.5}}) {
        const int draws = 20000;
        double sum = 0.0;
        double squares = 0.0;
        for (int k = 0; k < draws; ++k) {
            hicx::diffc::CounterRng rng({11, static_cast<std::uint64_t>(k)});
            const auto v = static_cast<double>(hicx::diffc::binomial(n, p, rng));
            CHECK(v >= 0.0);
            CHECK(v <= static_cast<double>(n));
            sum += v;
            squares += v * v;
        }
        const double mean = sum / draws;
        const double variance = squares / draws - mean * mean;
        const double expected_variance = static_cast<double>(n) * p * (1.0 - p);
        // Five standard errors of the mean; 10 % on the variance.
        CHECK(std::abs(mean - static_cast<double>(n) * p) <
              5.0 * std::sqrt(expected_variance / draws));
        CHECK(std::abs(variance / expected_variance - 1.0) < 0.1);
    }
}

TEST_CASE("diff contacts: MAD coverage outliers") {
    std::vector<double> coverage = {100, 102, 98, 101, 99, 0, 10, 400, 100, 97};
    const std::vector<char> invalid = hicx::diffc::coverage_outliers(coverage, -1.5, 5.0);
    CHECK(invalid[5] == 1);
    CHECK(invalid[6] == 1);
    CHECK(invalid[7] == 1);
    CHECK(invalid[0] == 0);
    CHECK(invalid[9] == 0);
}
