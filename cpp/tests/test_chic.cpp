// Unit tests for the cHi-C core components: hicx/chic_viewpoint.hpp,
// hicx/chic_hdf5.hpp and hicx/nbinom_mle.hpp.
//
// Expected values were produced by the reference environment (numpy 1.26.4,
// scipy 1.14.1, the repository's lib/viewpoint.py) on the cHi-C test data.

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <hdf5.h>

#include "doctest/doctest.h"
#include "hicx/chic_hdf5.hpp"
#include "hicx/chic_viewpoint.hpp"
#include "hicx/nbinom_mle.hpp"
#include "hicx/stats_ops.hpp"
#include "hicx/contingency_tests.hpp"

namespace {

const std::string kChic = std::string(HICX_TEST_DATA_DIR) + "/cHi-C/";

std::string temporary_path(const std::string& name) {
    return (std::filesystem::temp_directory_path() / ("hicx_test_chic_" + name)).string();
}

}  // namespace

TEST_CASE("chic: Python text semantics") {
    using namespace hicx::chic;
    CHECK(python_int(" 42\n") == 42);
    CHECK(python_int("-1_000") == -1000);
    CHECK_THROWS_AS((void)python_int("4.5"), ValueError);
    CHECK_THROWS_AS((void)python_int(""), ValueError);
    CHECK(python_float("0.000105073709\n") == 0.000105073709);
    CHECK(std::isnan(python_float("nan")));
    CHECK_THROWS_AS((void)python_float("x"), ValueError);
    CHECK(format_fixed(58.254577634279, 12) == "58.254577634279");
    CHECK(format_fixed(std::numeric_limits<double>::quiet_NaN(), 12) == "nan");
    CHECK(format_fixed(-std::numeric_limits<double>::quiet_NaN(), 12) == "nan");
    CHECK(format_fixed(2.5e-13, 12) == "0.000000000000");
    CHECK(strip(" \tchr1\t1\t2\r\n") == "chr1\t1\t2");
}

TEST_CASE("chic: reference point files") {
    const auto points = hicx::chic::read_reference_points(kChic + "referencePoints_edge_cases.bed");
    REQUIRE(points.points.size() == 7);
    CHECK(points.genes[4] == "ThreeColumns");
    // A three field line is (chrom, start, start).
    CHECK(points.points[4].start == "197195000");
    CHECK(points.points[4].end == "197195000");
    CHECK(hicx::chic::reference_point_string(points.points[1]) == "chr1_24202042_24205042");
    // A blank last line is skipped by the reader.
    const auto three = hicx::chic::read_reference_points(kChic + "referencePoints_chicViewpoint.bed");
    CHECK(three.points.size() == 3);
    CHECK(hicx::chic::read_lines(kChic + "referencePoints_chicViewpoint.bed").size() == 4);
}

TEST_CASE("chic: smoothing keeps the right border asymmetry of even windows") {
    std::vector<double> data;
    for (int i = 1; i <= 21; ++i) {
        data.push_back(std::pow(static_cast<double>(i), 1.5));
    }
    const auto mean = [&](std::size_t begin, std::size_t end) {
        double sum = 0.0;
        for (std::size_t i = begin; i < end; ++i) {
            sum += data[i];
        }
        return sum / static_cast<double>(end - begin);
    };
    for (std::int64_t window = 1; window <= 6; ++window) {
        const std::vector<double> smoothed =
            hicx::chic::smooth_interaction_values(std::span<const double>(data), window);
        const std::int64_t half = window / 2;
        const std::int64_t upstream = window % 2 == 0 ? half - 1 : half;
        const std::int64_t n = static_cast<std::int64_t>(data.size());
        for (std::int64_t i = 0; i < n; ++i) {
            std::int64_t begin = 0;
            std::int64_t end = 0;
            if (i >= n - half) {
                begin = i - half;
                end = n;
            } else {
                begin = std::max<std::int64_t>(i - upstream, 0);
                end = i + half + 1;
            }
            // Sequential and pairwise sums agree for at most seven terms.
            CHECK(smoothed[static_cast<std::size_t>(i)] ==
                  doctest::Approx(mean(static_cast<std::size_t>(begin),
                                       static_cast<std::size_t>(end))).epsilon(1e-15));
        }
    }
    // float32 means are single precision quantities widened afterwards.
    const std::vector<float> single{0.0F, 1.0F, 0.0F, 0.0F, 0.0F};
    const std::vector<double> smoothed =
        hicx::chic::smooth_interaction_values(std::span<const float>(single), 5);
    CHECK(smoothed[2] == static_cast<double>(1.0F / 5.0F));
    CHECK(smoothed[2] != 0.2);
}

TEST_CASE("chic: viewpoints on the real cHi-C matrix") {
    using namespace hicx::chic;
    const ViewpointMatrix matrix = ViewpointMatrix::load(kChic + "FL-E13-5_chr1.cool");
    CHECK(matrix.bin_size() == 1000);
    const ReferencePoint sox17{"chr1", "4487435", "4487435"};
    const ViewpointRange range = calculate_viewpoint_range(matrix, sox17, 500000, 500000);
    CHECK(range.region_start == 3987435);
    CHECK(range.region_end == 4987435);
    const ComputedViewpoint viewpoint =
        compute_viewpoint(matrix, sox17, "chr1", range.region_start, range.region_end);
    REQUIRE(viewpoint.data.size() == 1001);
    CHECK(viewpoint.index_before_viewpoint == 500);
    std::size_t nonzero = 0;
    double total = 0.0;
    for (const double value : viewpoint.data) {
        nonzero += value != 0.0 ? 1 : 0;
        total += value;
    }
    // chicQualityControl's raw filter: 0.14485514485514486
    CHECK(static_cast<double>(nonzero) / 1001.0 == 0.14485514485514486);
    // chicViewpoint's sum_of_interactions for Sox17
    CHECK(total == 810.0);

    // Clipped at the chromosome end: (197195432 - 197145200) + 1000
    const ViewpointRange near_end =
        calculate_viewpoint_range(matrix, ReferencePoint{"chr1", "197145200", "197145200"},
                                  200000, 200000);
    CHECK(near_end.region_end == 197195431);
    CHECK(near_end.downstream == 51232);

    // A position past the chromosome is a TypeError, a start after the end a
    // ValueError (numpy cannot broadcast), an unknown chromosome a ValueError.
    const ReferencePoint beyond{"chr1", "794979306", "794979306"};
    CHECK_THROWS_AS((void)compute_viewpoint(matrix, beyond, "chr1", 3987435, 4987435), TypeError);
    const ReferencePoint inverted{"chr1", "70931031", "70093803"};
    const ViewpointRange inverted_range = calculate_viewpoint_range(matrix, inverted, 500000, 500000);
    CHECK_THROWS_AS((void)compute_viewpoint(matrix, inverted, "chr1", inverted_range.region_start,
                                            inverted_range.region_end),
                    ValueError);
    CHECK_THROWS_AS((void)calculate_viewpoint_range(matrix, ReferencePoint{"chr2", "1", "1"}, 10, 10),
                    ValueError);
}

TEST_CASE("chic: relative values treat a zero denominator as absent") {
    const std::vector<double> data{1.0, 3.0};
    CHECK(hicx::chic::compute_relative_values(data, 8.0)[1] == 0.375);
    CHECK(hicx::chic::compute_relative_values(data, 0.0)[1] == 0.75);
}

TEST_CASE("chic: background model extension and the p-value lookup rule") {
    using namespace hicx::chic;
    const BackgroundModel model =
        read_background_model(kChic + "background.txt", 600000, 700000, 500000, false);
    CHECK(model.min_key() == -600000);
    CHECK(model.max_key() == 700000);
    CHECK(model.at(700000) == model.at(500000));
    CHECK(model.at(-500000)[0] == 58.254577634279);
    const BackgroundModel means =
        read_background_model(kChic + "background.txt", 200000, 200000, 500000, true);
    CHECK(interaction_background_data(means, 200000, 200000).size() == 401);

    const BackgroundModel plain =
        read_background_model(kChic + "background.txt", 200000, 200000, 500000, false);
    const std::vector<double> raw{9.8, 0.0, 26.6};
    const std::vector<double> p = p_values(plain, raw, 1);
    // Index 0 is one bin upstream, which is not a key in base pairs, so the
    // model at -500 kb is used; scipy gives 1 - betainc = 0 exactly here.
    CHECK(p[0] == 1.0 - hicx::stats::betainc(58.254577634279, 10.8, 0.998604709282));
    CHECK(p[0] == 0.0);
    CHECK(p[1] == 1.0);
    const std::vector<double>& upstream = plain.at(plain.max_key());
    CHECK(p[2] == 1.0 - hicx::stats::betainc(upstream[0], 27.6, upstream[1]));
}

TEST_CASE("stats::betainc is exactly 1.0 in the upper tail where an evaluation one ulp short would not be") {
    // The reference scipy.special.betainc(58.254577634279, x + 1, 0.998604709282)
    // is 1.0 for all four; chicViewpoint stores 1 - betainc as the p-value, so a
    // result of 1 - 1.1e-16 would be stored as a p-value of 1.1e-16 instead of 0.
    for (const double x : {9.8, 26.6, 75.0, 90.0}) {
        CHECK(hicx::stats::betainc(58.254577634279, x + 1.0, 0.998604709282) == 1.0);
    }
    // Values of scipy 1.14.1, compared bit for bit.
    CHECK(hicx::stats::betainc(2.5, 3.5, 0.5) == 0.6697652726313549);
    CHECK(hicx::stats::betainc(0.10000000000000001, 0.20000000000000001, 0.5) ==
          0.6705707961028995);
    CHECK(std::isnan(hicx::stats::betainc(0.0, 1.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(1.0, 1.0, 1.5)));
    CHECK(hicx::stats::betainc(2.0, 3.0, 0.0) == 0.0);
    CHECK(hicx::stats::betainc(2.0, 3.0, 1.0) == 1.0);
}

TEST_CASE("fit_nbinom_mle solves the profile score to full precision") {
    // Reference values: the root of the profile score
    //   sum_i [psi(x_i + r) - psi(r)] - n log(1 + mean / r)
    // found by bisection in mpmath with 60 digits (scratch generator
    // oracle_exact.py), prob = r / (r + mean).
    struct Case {
        std::vector<double> data;
        double size;
        double prob;
    };
    const std::vector<Case> cases{
        {{1, 2, 3, 5, 8, 2, 1, 0, 4}, 3.1809295484709001, 0.52405678708480785},
        {{1, 5}, 6.6042876143208726, 0.68763950847049549},
        {{0, 0, 1, 3, 2, 0, 7, 12, 5, 1, 0, 0, 2, 9, 4, 0, 1, 1, 3, 25}, 0.54325659891490552,
         0.12508047510953638},
        {{3, 4, 5, 6, 5, 4, 7, 2, 9, 12, 3, 1, 0, 6, 5, 4, 8, 5, 3, 11, 6, 4, 2, 7, 18},
         4.1582275080392357, 0.4261252880826476},
        {{0.5, 1.7, 3.2, 0.1, 9.9, 4.4, 0.0, 2.5}, 1.0862892035212945, 0.28042031882732593},
    };
    for (const Case& c : cases) {
        const auto fit = hicx::stats::fit_nbinom_mle(c.data);
        CHECK(fit.status == 0);
        CHECK(fit.size == doctest::Approx(c.size).epsilon(1e-12));
        CHECK(fit.prob == doctest::Approx(c.prob).epsilon(1e-12));
        // Order and repetition of the values do not change the answer.
        std::vector<double> reversed(c.data.rbegin(), c.data.rend());
        const auto again = hicx::stats::fit_nbinom_mle(reversed);
        CHECK(again.size == doctest::Approx(fit.size).epsilon(1e-13));
    }
}

TEST_CASE("fit_nbinom_mle: the score vanishes at the fit and the likelihood is maximal there") {
    const std::vector<double> data{3, 4, 5, 6, 5, 4, 7, 2, 9, 12, 3, 1, 0, 6, 5, 4, 8, 5, 3, 11, 6, 4, 2, 7, 18};
    const auto fit = hicx::stats::fit_nbinom_mle(data);
    double mean = 0.0;
    for (const double x : data) {
        mean += x;
    }
    mean /= static_cast<double>(data.size());
    const std::vector<double> values(data.begin(), data.end());
    const std::vector<double> ones(data.size(), 1.0);
    const double step = 1e-4 * fit.size;
    CHECK(hicx::stats::nbinom_profile_score(values, ones, mean, fit.size - step) > 0.0);
    CHECK(hicx::stats::nbinom_profile_score(values, ones, mean, fit.size + step) < 0.0);
    // prob is the closed form maximiser for the fitted size
    CHECK(fit.prob == doctest::Approx(fit.size / (fit.size + mean)).epsilon(1e-15));
}

TEST_CASE("fit_nbinom_mle returns the Python's answer where there is nothing to fit") {
    const std::vector<double> empty;
    const auto none = hicx::stats::fit_nbinom_mle(empty);
    CHECK(none.size == 10.0);
    CHECK(std::isnan(none.prob));
    CHECK(none.status == 2);
    const auto given_empty = hicx::stats::fit_nbinom_mle(empty, std::array<double, 2>{2.0, 0.5});
    CHECK(given_empty.size == 2.0);
    CHECK(given_empty.prob == 0.5);
    // factorial overflows above 170: the objective is infinite everywhere and
    // the starting point, the moment estimator, is returned.
    const std::vector<double> overflow{177.0, 0.0, 2.0, 5.0};
    const auto start = hicx::stats::fit_nbinom_mle(overflow);
    const double mean = 46.0;
    const double variance = (131.0 * 131.0 + 46.0 * 46.0 + 44.0 * 44.0 + 41.0 * 41.0) / 4.0;
    const double size = (mean * mean) / (variance - mean);
    CHECK(start.size == doctest::Approx(size).epsilon(1e-15));
    CHECK(start.prob == doctest::Approx(size / (size + mean)).epsilon(1e-15));
    // A start given by the caller replaces the moment start there, clamped.
    const auto given = hicx::stats::fit_nbinom_mle(overflow, std::array<double, 2>{-1.0, 5.0});
    CHECK(given.size == std::numeric_limits<double>::epsilon());
    CHECK(given.prob == 1.0);
    // all zeros: size 10, prob 1
    const auto zeros = hicx::stats::fit_nbinom_mle(std::vector<double>(7, 0.0));
    CHECK(zeros.size == 10.0);
    CHECK(zeros.prob == 1.0);
}

TEST_CASE("fit_nbinom_mle: variance not above the mean is the Poisson limit") {
    const std::vector<double> underdispersed{2.0, 3.0, 2.0, 3.0, 2.0, 3.0};
    const auto fit = hicx::stats::fit_nbinom_mle(underdispersed);
    CHECK(fit.status == 1);
    CHECK(fit.size == hicx::stats::kPoissonLimitSize);
    CHECK(fit.prob == doctest::Approx(fit.size / (fit.size + 2.5)).epsilon(1e-15));
    // the mean of NB(size, prob) is size (1 - prob) / prob: the sample mean
    CHECK(fit.size * (1.0 - fit.prob) / fit.prob == doctest::Approx(2.5).epsilon(1e-12));
    const auto single = hicx::stats::fit_nbinom_mle(std::vector<double>{4.0});
    CHECK(single.status == 1);
}

TEST_CASE("chic: contingency tests agree with exact rational arithmetic and scipy 1.14.1") {
    namespace ct = hicx::contingency;
    // Exact Fisher p-values (Python fractions, tables whose probability does not
    // exceed that of the observed one, tie tolerance 1e-12); scipy 1.14.1
    // returns the same to 1e-15 relative.
    struct Table {
        std::int64_t a, b, c, d;
        double p;
    };
    const std::vector<Table> tables{
        {673, 6, 832, 10, 0.62176500764164044},
        {810, 298, 978, 401, 0.24331170304805214},
        {5, 0, 0, 5, 0.0079365079365079361},
        {20, 3, 900, 150, 1.0},
        {3, 7, 7, 3, 0.17889540799757517},
        {1, 9, 11, 3, 0.0027594561852200836},
        {2, 2, 2, 2, 1.0},
        {0, 5, 5, 0, 0.0079365079365079361},
        {12, 5, 29, 3, 0.10642628903375535},
        {100, 150, 120, 90, 0.00026209889973324846},
        {7, 0, 1, 6, 0.004662004662004662},
        {40, 20, 10, 30, 8.3086587062391521e-05},
        {1, 0, 0, 1, 1.0},
        {250, 230, 260, 300, 0.07147456273562916},
    };
    for (const Table& t : tables) {
        const auto result = ct::fisher_exact_2x2(t.a, t.b, t.c, t.d);
        REQUIRE(result.has_value());
        CHECK(result->pvalue == doctest::Approx(t.p).epsilon(1e-12));
        // the mirror image tables have the same p-value
        const auto swapped = ct::fisher_exact_2x2(t.d, t.c, t.b, t.a);
        CHECK(swapped->pvalue == doctest::Approx(result->pvalue).epsilon(1e-12));
        const auto transposed = ct::fisher_exact_2x2(t.a, t.c, t.b, t.d);
        CHECK(transposed->pvalue == doctest::Approx(result->pvalue).epsilon(1e-12));
    }
    CHECK_FALSE(ct::fisher_exact_2x2(-1, 3, 4, 5).has_value());
    // an empty row or column, the mode and a total above one are the constant 1.0
    CHECK(ct::fisher_exact_2x2(0, 0, 3, 4)->constant_one);
    CHECK(ct::fisher_exact_2x2(2, 2, 2, 2)->constant_one);
    CHECK_FALSE(ct::fisher_exact_2x2(5, 0, 0, 5)->constant_one);

    // closed form: the hypergeometric pmf of a small table
    CHECK(ct::hypergeometric_pmf(2, 10, 4, 5) == doctest::Approx(6.0 * 20.0 / 252.0).epsilon(1e-14));
    CHECK(ct::hypergeometric_pmf(5, 10, 4, 5) == 0.0);
    // sum of the pmf over the support is one
    double total = 0.0;
    for (int k = 0; k <= 4; ++k) {
        total += ct::hypergeometric_pmf(k, 10, 4, 5);
    }
    CHECK(total == doctest::Approx(1.0).epsilon(1e-14));

    // chi2_contingency of scipy 1.14.1 with and without Yates' correction
    struct Chi {
        double a, b, c, d;
        bool yates;
        double statistic;
        double pvalue;
    };
    const std::vector<Chi> chis{
        {673.0, 5.6, 832.0, 0.0, false, 6.891447919451126, 0.00866090872217125},
        {673.0, 5.6, 832.0, 0.0, true, 4.83821474664194, 0.02783569741540067},
        {810.0, 297.8, 978.0, 400.4, false, 1.4268145854722973, 0.23228500408601602},
        {810.0, 297.8, 978.0, 400.4, true, 1.3215803522506986, 0.2503086237278468},
        {12.0, 5.0, 29.0, 3.0, false, 3.2626939114060254, 0.07087260587160922},
        {12.0, 5.0, 29.0, 3.0, true, 1.9608156048242464, 0.16142611786143468},
        {3.0, 7.0, 7.0, 3.0, false, 3.2, 0.07363827012030258},
        {3.0, 7.0, 7.0, 3.0, true, 1.8, 0.17971249487899593},
        {100.5, 150.2, 120.7, 90.1, false, 13.527390785252017, 0.00023510668972423121},
        {1.0, 1.0, 5.0, 0.3, true, 0.09734216255442668, 0.7550427612004261},
    };
    for (const Chi& c : chis) {
        const auto r = ct::chi2_2x2(c.a, c.b, c.c, c.d, c.yates);
        REQUIRE(r.has_value());
        CHECK(r->statistic == doctest::Approx(c.statistic).epsilon(1e-13));
        CHECK(r->pvalue == doctest::Approx(c.pvalue).epsilon(1e-11));
    }
    // closed form: the statistic is N (ad - bc)^2 / (row1 row2 col1 col2)
    const auto closed = ct::chi2_2x2(3.0, 7.0, 7.0, 3.0);
    CHECK(closed->statistic == doctest::Approx(20.0 * 40.0 * 40.0 / (10.0 * 10.0 * 10.0 * 10.0)).epsilon(1e-15));
    // an expected frequency of zero (an empty row) and a negative cell are refused
    CHECK_FALSE(ct::chi2_2x2(673.0, 0.0, 832.0, 0.0).has_value());
    CHECK_FALSE(ct::chi2_2x2(-1.0, 1.0, 2.0, 3.0).has_value());

    // chi2.ppf(q, 1) and chi2.sf(x, 1) of scipy 1.14.1
    CHECK(ct::chi2_quantile(0.5, 1.0) == doctest::Approx(0.454936423119572).epsilon(1e-13));
    CHECK(ct::chi2_quantile(0.9, 1.0) == doctest::Approx(2.705543454095404).epsilon(1e-13));
    CHECK(ct::chi2_quantile(0.95, 1.0) == doctest::Approx(3.841458820694124).epsilon(1e-13));
    CHECK(ct::chi2_quantile(0.99, 1.0) == doctest::Approx(6.6348966010212145).epsilon(1e-13));
    CHECK(ct::chi2_quantile(0.999, 1.0) == doctest::Approx(10.827566170662733).epsilon(1e-13));
    CHECK(ct::chi2_quantile(0.999999999, 1.0) == doctest::Approx(37.32489310651872).epsilon(1e-9));
    CHECK(ct::chi2_quantile(0.0, 1.0) == 0.0);
    CHECK(ct::chi2_quantile(1.0, 1.0) == std::numeric_limits<double>::infinity());
    CHECK(std::isnan(ct::chi2_quantile(1.5, 1.0)));
    CHECK(ct::chi2_sf(0.0, 1.0) == 1.0);
    CHECK(ct::chi2_sf(3.84, 1.0) == doctest::Approx(0.05004352124870519).epsilon(1e-12));
    CHECK(ct::chi2_sf(10.0, 1.0) == doctest::Approx(0.001565402258002549).epsilon(1e-12));
    CHECK(ct::chi2_sf(40.0, 1.0) == doctest::Approx(2.5396285894708634e-10).epsilon(1e-10));
    CHECK(ct::chi2_sf(-1.0, 1.0) == 1.0);
}

TEST_CASE("chic_hdf5 writes h5py's layout") {
    const std::string path = temporary_path("layout.hdf5");
    {
        hicx::chic::Hdf5Writer writer(path);
        writer.set_attribute("/", "type", std::string("interactions"));
        const std::vector<std::int64_t> range{200000, 100000};
        writer.set_attribute("/", "range", std::span<const std::int64_t>(range));
        writer.create_group("m/chr1/g");
        CHECK_THROWS((void)writer.create_group("m/chr1/g"));
        hicx::chic::InteractionFileData data;
        data.chromosome = "chr1";
        data.gene = "g";
        data.sum_of_interactions = 3.0;
        data.starts = {0, 1000};
        data.ends = {1000, 2000};
        data.relative_positions = {0, 1000};
        data.interaction_data = {0.5, 0.25};
        data.raw = {1.5, 0.75};
        data.pvalues = {1.0, 0.5};
        data.xfold = {2.0, 4.0};
        hicx::chic::write_interaction_datasets(writer, "m/chr1/g", data, 5, 6);
        CHECK(writer.hard_link("m/chr1/g", "m/genes/g"));
        CHECK_FALSE(writer.hard_link("m/chr1/g", "m/genes/g"));
    }
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file >= 0);
    const hid_t dataset = H5Dopen2(file, "m/genes/g/pvalue", H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    const hid_t plist = H5Dget_create_plist(dataset);
    CHECK(H5Pget_layout(plist) == H5D_CHUNKED);
    CHECK(H5Pget_nfilters(plist) == 1);
    unsigned int flags = 0;
    std::size_t count = 1;
    unsigned int level = 0;
    char name[32];
    CHECK(H5Pget_filter2(plist, 0, &flags, &count, &level, sizeof(name), name, nullptr) ==
          H5Z_FILTER_DEFLATE);
    CHECK(level == 9);
    H5Pclose(plist);
    H5Dclose(dataset);
    const hid_t gene = H5Dopen2(file, "m/chr1/g/gene", H5P_DEFAULT);
    const hid_t type = H5Dget_type(gene);
    CHECK(H5Tis_variable_str(type) > 0);
    CHECK(H5Tget_cset(type) == H5T_CSET_UTF8);
    H5Tclose(type);
    H5Dclose(gene);
    H5Fclose(file);
    std::remove(path.c_str());
}

TEST_CASE("chic: collapse_viewpoint is the second half of compute_viewpoint") {
    using namespace hicx::chic;
    // A wide reference point (three bins) on the real matrix: the binding sums
    // the rows itself and collapses them with collapse_viewpoint.
    const ViewpointMatrix matrix = ViewpointMatrix::load(kChic + "MB-E10-5_chr1.cool");
    const ReferencePoint wide{"chr1", "24202042", "24205042"};
    const ComputedViewpoint expected = compute_viewpoint(matrix, wide, "chr1", 24102042, 24305042);
    const auto [first, last] = matrix.reference_point_indices(wide);
    const auto range = matrix.region_bin_range("chr1", 24102042, 24305042);
    REQUIRE(range.has_value());
    std::vector<double> rows(static_cast<std::size_t>(range->second + 1 - range->first), 0.0);
    for (std::int64_t row = first; row <= last; ++row) {
        for (std::size_t j = 0; j < rows.size(); ++j) {
            rows[j] += matrix.value(row, range->first + static_cast<std::int64_t>(j));
        }
    }
    const ComputedViewpoint collapsed = collapse_viewpoint(rows, first, last, range->first);
    CHECK(collapsed.data == expected.data);
    CHECK(collapsed.index_before_viewpoint == expected.index_before_viewpoint);
    CHECK(collapsed.data.size() == 201);
    // A start after the end lengthens the result (numpy's slices), and an
    // index before the region is numpy's IndexError.
    CHECK(collapse_viewpoint(std::vector<double>(10, 1.0), 5, 2, 0).data.size() == 13);
    CHECK_THROWS_AS((void)collapse_viewpoint(std::vector<double>(3, 1.0), 0, 0, 9), IndexError);
}

TEST_CASE("chic: H5Reader gives h5py's view of a cHi-C file") {
    using namespace hicx::chic;
    const H5Reader reader(kChic + "chicViewpoint/two_matrices.hdf5");
    CHECK(reader.kind("FL-E13-5_chr1/chr1/Eya1") == H5Reader::Kind::Group);
    CHECK(reader.kind("FL-E13-5_chr1/chr1/Eya1/pvalue") == H5Reader::Kind::Dataset);
    CHECK(reader.kind("FL-E13-5_chr1/chr1/Eya1/pvalue/x") == H5Reader::Kind::Missing);
    CHECK(reader.kind("FL-E13-5_chr1/chrX/Eya1") == H5Reader::Kind::Missing);
    CHECK(reader.keys("FL-E13-5_chr1/chr1") == std::vector<std::string>{"Eya1", "Sox17", "Tfap2d"});
    const H5Value positions = reader.read("FL-E13-5_chr1/chr1/Eya1/relative_position_list");
    CHECK(positions.dtype_str() == "<i8");
    CHECK(positions.length == 401);
    CHECK(positions.as_int64(0) == -200000);
    const H5Value chromosome = reader.read("FL-E13-5_chr1/chr1/Eya1/chromosome");
    CHECK(chromosome.scalar);
    CHECK(chromosome.variable_length);
    CHECK(chromosome.strings == std::vector<std::string>{"chr1"});
    const H5Value start = reader.read("FL-E13-5_chr1/chr1/Eya1/reference_point_start");
    CHECK(start.scalar);
    CHECK(start.as_int64(0) == 14300280);
    // h5py's KeyError messages
    CHECK(reader.missing_message("FL-E13-5_chr1/x") ==
          "Unable to synchronously open object (object 'x' doesn't exist)");
    CHECK(reader.missing_message("FL-E13-5_chr1/x/y") ==
          "Unable to synchronously open object (component not found)");
    try {
        (void)H5Reader(kChic + "no_such_file.hdf5");
        FAIL("no exception");
    } catch (const hicx::PythonError& error) {
        CHECK(error.type() == "FileNotFoundError");
    }
    try {
        (void)H5Reader(kChic + "referencePoints.bed");
        FAIL("no exception");
    } catch (const hicx::PythonError& error) {
        CHECK(error.type() == "OSError");
        CHECK(std::string(error.what()) == "Unable to synchronously open file (file signature not found)");
    }
    // The target file stores its coordinates as strings.
    const H5Reader target(kChic + "chicSignificantInteractions/targetFile_dual.hdf5");
    const H5Value starts = target.read("FL-E13-5_chr1/MB-E10-5_chr1/chr1/Eya1/start_list");
    CHECK(starts.dtype_str() == "|O");
    CHECK(starts.strings.at(0) == "14274000");
}

TEST_CASE("chic: target regions are sorted and merged as bedtools merge -d does") {
    using hicx::chic::TargetRegion;
    const std::vector<TargetRegion> merged = hicx::chic::sort_merge_regions(
        {{"chr2", 50, 60}, {"chr1", 30, 40}, {"chr1", 0, 10}, {"chr1", 11, 20}, {"chr1", 41, 45}}, 1);
    REQUIRE(merged.size() == 3);
    CHECK(merged[0].start == 0);
    CHECK(merged[0].end == 20);
    CHECK(merged[1].start == 30);
    CHECK(merged[1].end == 45);
    CHECK(merged[2].chromosome == "chr2");
}

TEST_CASE("chic: fisher_exact_2x2 marks the p-values that are the constant 1.0") {
    // [[10, 4], [5, 7]]: hypergeometric tails, scipy 1.14.1 gives 0.2328577197469377
    const auto computed = hicx::contingency::fisher_exact_2x2(10, 4, 5, 7);
    REQUIRE(computed.has_value());
    CHECK_FALSE(computed->constant_one);
    CHECK(computed->pvalue == doctest::Approx(0.2328577197469377).epsilon(1e-12));
    // an empty row, and a table at the mode: the constant 1.0
    CHECK(hicx::contingency::fisher_exact_2x2(0, 0, 3, 0)->constant_one);
    CHECK(hicx::contingency::fisher_exact_2x2(3, 3, 3, 3)->constant_one);
    CHECK(hicx::contingency::fisher_exact_2x2(3, 3, 3, 3)->pvalue == 1.0);
    CHECK_FALSE(hicx::contingency::fisher_exact_2x2(-1, 2, 3, 4).has_value());
}
