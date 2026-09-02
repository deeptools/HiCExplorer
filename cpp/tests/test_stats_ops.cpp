// Unit tests for hicx::stats and hicx::simd.
//
// Every expected value in this file was produced by the reference oracle,
// scipy 1.14.1 and numpy in ~/miniconda3/envs/__hicexplorer@3.7.6, on the
// inputs written here. The inputs are integer recurrences divided by a power
// of two, so both sides see exactly the same bits with no parsing question.

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

TEST_CASE("cephes ndtr matches scipy.special.ndtr") {
    // scipy.special.ndtr, full float64 repr.
    CHECK(hicx::stats::ndtr(-6.0) == 9.865876450376946e-10);
    CHECK(hicx::stats::ndtr(-2.5) == 0.006209665325776132);
    CHECK(hicx::stats::ndtr(-1.0) == 0.15865525393145707);
    CHECK(hicx::stats::ndtr(-0.5) == 0.3085375387259869);
    CHECK(hicx::stats::ndtr(0.0) == 0.5);
    CHECK(hicx::stats::ndtr(0.3) == 0.6179114221889526);
    CHECK(hicx::stats::ndtr(0.7071) == 0.7602478320101365);
    CHECK(hicx::stats::ndtr(1.0) == 0.8413447460685429);
    CHECK(hicx::stats::ndtr(2.5) == 0.9937903346742238);
    CHECK(hicx::stats::ndtr(6.0) == 0.9999999990134123);
    CHECK(hicx::stats::ndtr(9.0) == 1.0);
    // norm.sf is ndtr(-x), which is the route scipy's rank sum test takes.
    CHECK(hicx::stats::normal_sf(2.5) == hicx::stats::ndtr(-2.5));
}

TEST_CASE("cephes erf and erfc match scipy") {
    CHECK(hicx::stats::erf(0.1) == 0.1124629160182849);
    CHECK(hicx::stats::erf(0.5) == 0.5204998778130465);
    CHECK(hicx::stats::erf(1.0) == 0.8427007929497148);
    CHECK(hicx::stats::erf(2.0) == 0.9953222650189527);
    CHECK(hicx::stats::erf(5.0) == 0.9999999999984626);
    CHECK(hicx::stats::erf(-1.5) == -0.9661051464753108);
    CHECK(hicx::stats::erfc(0.1) == 0.8875370839817152);
    CHECK(hicx::stats::erfc(1.0) == 0.15729920705028516);
    CHECK(hicx::stats::erfc(2.0) == 0.004677734981047266);
    CHECK(hicx::stats::erfc(5.0) == 1.5374597944280347e-12);
    CHECK(hicx::stats::erfc(-1.5) == 1.9661051464753108);
}

TEST_CASE("digamma matches scipy.special.psi at the branch boundaries") {
    // The four branches: a small positive integer, the [1, 2] rational
    // approximation, the downward recurrence and the asymptotic series.
    CHECK(hicx::stats::digamma(1.0) == -0.5772156649015329);
    CHECK(hicx::stats::digamma(1.5) == 0.03648997397857652);
    CHECK(hicx::stats::digamma(2.0) == 0.42278433509846713);
    CHECK(hicx::stats::digamma(0.25) == -4.2274535333762655);
    CHECK(hicx::stats::digamma(7.5) == 1.9467574842460866);
    CHECK(hicx::stats::digamma(100.0) == 4.600161852738088);
    // The reflection branch for a negative argument.
    CHECK(hicx::stats::digamma(-0.5) == doctest::Approx(0.03648997397857651).epsilon(1e-15));
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
