// Unit tests for the computational core of hicCompartmentalization
// (tools/compartmentalization_impl.hpp).
//
// Every expected value was produced by the reference environment (numpy 1.26.4,
// scipy 1.14.1), either by numpy directly or by running
// hicexplorer.hicCompartmentalization.count_interactions and
// within_vs_between_compartments on the same fixture, and is quoted as a hex
// float so that the comparison is bit for bit. Nothing is asserted against a
// number this code computed itself.

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/compartmentalization_impl.hpp"
#include "hicx/sparse_matrix.hpp"

namespace cm = hicx::compartments;

namespace {

const double kInf = std::numeric_limits<double>::infinity();
const double kNan = std::numeric_limits<double>::quiet_NaN();

void check_bits(const std::vector<double>& actual, const std::vector<double>& expected) {
    REQUIRE(actual.size() == expected.size());
    for (std::size_t k = 0; k < expected.size(); ++k) {
        CAPTURE(k);
        CAPTURE(actual[k]);
        CAPTURE(expected[k]);
        if (std::isnan(expected[k])) {
            CHECK(std::isnan(actual[k]));
        } else {
            CHECK(actual[k] == expected[k]);
            CHECK(std::signbit(actual[k]) == std::signbit(expected[k]));
        }
    }
}

// The 6 bin fixture: an upper triangle as hicmatrix stores it, with a NaN at
// (1, 5) and an infinity at (3, 4).
//
//   U[0,0]=4    U[0,1]=1.5  U[0,3]=0.25  U[1,1]=3   U[1,2]=2   U[1,5]=nan
//   U[2,2]=5.5  U[2,4]=0.125 U[3,3]=1    U[3,4]=inf U[4,4]=2   U[4,5]=0.75
//   U[5,5]=6
hicx::CsrMatrix fixture_upper() {
    std::vector<std::int64_t> indptr = {0, 3, 6, 8, 10, 12, 13};
    std::vector<std::int32_t> indices = {0, 1, 3, 1, 2, 5, 2, 4, 3, 4, 4, 5, 5};
    std::vector<double> data = {4.0, 1.5, 0.25, 3.0, 2.0, kNan, 5.5,
                                0.125, 1.0, kInf, 2.0, 0.75, 6.0};
    hicx::CsrMatrix matrix(6, 6, std::move(indptr), std::move(indices), std::move(data),
                           "float64");
    matrix.set_symmetry(hicx::Symmetry::UpperTriangle);
    return matrix;
}

// pc1 rows: bin_id, quantile and chromosome as the Python fixture had them.
const std::vector<std::vector<std::int64_t>> kBinIds = {{0, 1}, {2}, {3, 4}, {5}, {1, 2}};
const std::vector<std::int64_t> kQuantiles = {1, 0, 2, 1, 0};
const std::int64_t kChromosomes = 3;  // "a", "b", "c"

const std::vector<double> kCountNone = {
    0x1.d555555555555p+1, 0x1.1000000000000p+0, 0x1.5555555555555p-5, 0x0.0p+0,
    0x1.1000000000000p+0, 0x1.2492492492492p+1, 0x1.5555555555555p-3, 0x0.0p+0,
    0x1.5555555555555p-5, 0x1.5555555555555p-3, 0x1.8000000000000p+0, 0x0.0p+0,
    0x0.0p+0,             0x0.0p+0,             0x0.0p+0,             0x0.0p+0};
const std::vector<double> kRatioNone = {kInf, 0x1.6fc57c57c57c6p+4, 0x1.59daf304859dap+0};

const std::vector<double> kCountOffset0 = {
    0x1.0000000000000p+1, 0x1.9249249249249p-1, 0x1.5555555555555p-5, 0x0.0p+0,
    0x1.9249249249249p-1, 0x1.8000000000000p-1, 0x1.5555555555555p-3, 0x0.0p+0,
    0x1.5555555555555p-5, 0x1.5555555555555p-3, 0x0.0p+0,             0x0.0p+0,
    0x0.0p+0,             0x0.0p+0,             0x0.0p+0,             0x0.0p+0};
const std::vector<double> kRatioOffset0 = {kInf, 0x1.4be2be2be2be3p+3, 0x1.85fb37072d754p+0};

const std::vector<double> kCountOffset02 = {
    0x1.0000000000000p+1, 0x1.199999999999ap+0, 0x0.0p+0, 0x0.0p+0,
    0x1.199999999999ap+0, 0x1.8000000000000p-1, 0x1.0000000000000p-2, 0x0.0p+0,
    0x0.0p+0,             0x1.0000000000000p-2, 0x0.0p+0, 0x0.0p+0,
    0x0.0p+0,             0x0.0p+0,             0x0.0p+0, 0x0.0p+0};
const std::vector<double> kRatioOffset02 = {kInf, 0x1.3cccccccccccdp+3, 0x1.6cefa8d9df51bp+0};

}  // namespace

TEST_CASE("nanquantile_linear matches numpy on float32 widened values with NaN") {
    std::vector<double> values;
    for (const double literal : {0.491006273550, kNan, -0.25, 3.5, 1e-3, -7.125,
                                 0.3333333333, kNan, 2.0}) {
        values.push_back(static_cast<double>(static_cast<float>(literal)));
    }
    std::vector<double> quantiles;
    for (int j = 0; j < 8; ++j) {
        quantiles.push_back(static_cast<double>(j) / 7.0);
    }
    check_bits(cm::nanquantile_linear(values, quantiles),
               {-0x1.c800000000000p+2, -0x1.3b6db6db6db6fp+0, -0x1.21a54d85b6db9p-4,
                0x1.86f9142bffffdp-3, 0x1.9a877c4924924p-2, 0x1.d823b26db6db5p-1,
                0x1.1b6db6db6db6cp+1, 0x1.c000000000000p+1});
    check_bits(cm::nanquantile_linear(values, {5.0 / 100, (100 - 5.0) / 100}),
               {-0x1.4400000000000p+2, 0x1.8666666666664p+1});
}

TEST_CASE("nanquantile_linear rejects quantiles outside [0, 1] and NaN") {
    CHECK_THROWS_AS((void)cm::nanquantile_linear({1.0, 2.0}, {1.5}), std::invalid_argument);
    CHECK_THROWS_AS((void)cm::nanquantile_linear({1.0, 2.0}, {kNan}), std::invalid_argument);
    const auto all_nan = cm::nanquantile_linear({kNan, kNan}, {0.0, 1.0});
    REQUIRE(all_nan.size() == 2);
    CHECK(std::isnan(all_nan[0]));
    CHECK(std::isnan(all_nan[1]));
}

TEST_CASE("linspace matches numpy, including the step == 0 branch") {
    check_bits(cm::linspace(-0.3, 0.7, 7),
               {-0x1.3333333333333p-2, -0x1.1111111111111p-3, 0x1.1111111111110p-5,
                0x1.999999999999ap-3, 0x1.7777777777777p-2, 0x1.1111111111110p-1,
                0x1.6666666666666p-1});
    check_bits(cm::linspace(5e-324, 1e-323, 3),
               {0x0.0000000000001p-1022, 0x0.0000000000001p-1022, 0x0.0000000000002p-1022});
    check_bits(cm::linspace(2.5, 9.0, 1), {0x1.4000000000000p+1});
    CHECK(cm::linspace(1.0, 2.0, 0).empty());
    CHECK_THROWS_AS((void)cm::linspace(1.0, 2.0, -2), std::invalid_argument);
}

TEST_CASE("searchsorted_right sorts NaN last and seeds from the previous key") {
    CHECK(cm::searchsorted_right({-1.0, 0.0, 0.0, 2.0, kNan},
                                 {0.0, -5.0, kNan, 2.0, 1.5, -1.0, 3.0}) ==
          std::vector<std::int64_t>{3, 0, 5, 4, 3, 1, 4});
    // An unsorted array exposes the seeding. Found by a search against numpy:
    // a plain binary search answers {3, 3, 3} here, numpy {3, 5, 5}, because
    // the second key starts from the bounds the first one left.
    CHECK(cm::searchsorted_right({1.0, 2.0, 1.0, 4.0, 1.0, 5.0}, {2.5, 3.5, 3.5}) ==
          std::vector<std::int64_t>{3, 5, 5});
}

TEST_CASE("SymmetricRows answers the rows fillLowerTriangle would produce") {
    const hicx::CsrMatrix upper = fixture_upper();
    hicx::CsrMatrix full = fixture_upper();
    full.materialize_full();
    const cm::SymmetricRows rows(upper);
    for (std::int64_t r = 0; r < 6; ++r) {
        std::vector<double> dense(6, 0.0);
        std::vector<int> seen(6, 0);
        rows.for_each_in_row(r, [&](std::int64_t c, double v) {
            dense[static_cast<std::size_t>(c)] = v;
            ++seen[static_cast<std::size_t>(c)];
        });
        for (std::int64_t c = 0; c < 6; ++c) {
            CAPTURE(r);
            CAPTURE(c);
            CHECK(seen[static_cast<std::size_t>(c)] <= 1);
            const double expected = full.at(r, c);
            if (std::isnan(expected)) {
                CHECK(std::isnan(dense[static_cast<std::size_t>(c)]));
            } else {
                CHECK(dense[static_cast<std::size_t>(c)] == expected);
            }
        }
    }
}

TEST_CASE("normalised_sum_per_quantile and within_vs_between match the Python functions") {
    const hicx::CsrMatrix upper = fixture_upper();
    hicx::CsrMatrix full = fixture_upper();
    full.materialize_full();

    struct Variant {
        std::vector<std::int64_t> offsets;
        const std::vector<double>* count;
        const std::vector<double>* ratio;
    };
    const std::vector<Variant> variants = {{{}, &kCountNone, &kRatioNone},
                                           {{0}, &kCountOffset0, &kRatioOffset0},
                                           {{0, 2}, &kCountOffset02, &kRatioOffset02}};
    for (const auto& variant : variants) {
        CAPTURE(variant.offsets.size());
        const std::array<const hicx::CsrMatrix*, 2> matrices = {&upper, &full};
        for (const hicx::CsrMatrix* matrix : matrices) {
            const auto normalised = cm::normalised_sum_per_quantile(
                *matrix, kBinIds, kQuantiles, 4, variant.offsets, kChromosomes);
            check_bits(normalised, *variant.count);
            check_bits(cm::within_vs_between(normalised, 4), *variant.ratio);
        }
    }
}

TEST_CASE("the chromosome count changes the last bits exactly as the repeated additions do") {
    // The same fixture with non dyadic values, so that s + s + s is not 3 * s.
    // count_interactions was run with one, three and five distinct chromosome
    // labels; the results differ from each other in the last bit of up to 5
    // of the 16 cells, and the replay has to land on each of them.
    std::vector<std::int64_t> indptr = {0, 3, 6, 8, 10, 12, 13};
    std::vector<std::int32_t> indices = {0, 1, 3, 1, 2, 5, 2, 4, 3, 4, 4, 5, 5};
    std::vector<double> data = {0.1, 0.7, 1.0 / 3.0, 0.3, 2.2, kNan, 5.5,
                                0.125, 1.1, kInf, 0.9, 0.75, 6.0};
    hicx::CsrMatrix upper(6, 6, std::move(indptr), std::move(indices), std::move(data),
                          "float64");
    upper.set_symmetry(hicx::Symmetry::UpperTriangle);

    const std::vector<double> one = {
        0x1.ba4fa4fa4fa50p+1, 0x1.599999999999ap-1, 0x1.5555555555555p-5, 0x0.0p+0,
        0x1.599999999999ap-1, 0x1.1d41d41d41d42p+0, 0x1.71c71c71c71c7p-3, 0x0.0p+0,
        0x1.5555555555555p-5, 0x1.71c71c71c71c7p-3, 0x1.0000000000000p+0, 0x0.0p+0,
        0x0.0p+0,             0x0.0p+0,             0x0.0p+0,             0x0.0p+0};
    const std::vector<double> three = {
        0x1.ba4fa4fa4fa4fp+1, 0x1.5999999999999p-1, 0x1.5555555555555p-5, 0x0.0p+0,
        0x1.5999999999999p-1, 0x1.1d41d41d41d42p+0, 0x1.71c71c71c71c6p-3, 0x0.0p+0,
        0x1.5555555555555p-5, 0x1.71c71c71c71c6p-3, 0x1.0000000000000p+0, 0x0.0p+0,
        0x0.0p+0,             0x0.0p+0,             0x0.0p+0,             0x0.0p+0};
    const std::vector<double> five = {
        0x1.ba4fa4fa4fa50p+1, 0x1.5999999999999p-1, 0x1.5555555555555p-5, 0x0.0p+0,
        0x1.5999999999999p-1, 0x1.1d41d41d41d41p+0, 0x1.71c71c71c71c7p-3, 0x0.0p+0,
        0x1.5555555555555p-5, 0x1.71c71c71c71c7p-3, 0x1.0000000000000p+0, 0x0.0p+0,
        0x0.0p+0,             0x0.0p+0,             0x0.0p+0,             0x0.0p+0};
    check_bits(cm::normalised_sum_per_quantile(upper, kBinIds, kQuantiles, 4, {}, 1), one);
    check_bits(cm::normalised_sum_per_quantile(upper, kBinIds, kQuantiles, 4, {}, 3), three);
    check_bits(cm::normalised_sum_per_quantile(upper, kBinIds, kQuantiles, 4, {}, 5), five);
}

TEST_CASE("a quantile outside the range, as for NaN or the maximum pc1, is never counted") {
    const hicx::CsrMatrix upper = fixture_upper();
    std::vector<std::int64_t> quantiles = kQuantiles;
    quantiles.push_back(4);  // == Q, what searchsorted returns for NaN
    std::vector<std::vector<std::int64_t>> bin_ids = kBinIds;
    bin_ids.push_back({0, 1, 2, 3, 4, 5});
    check_bits(cm::normalised_sum_per_quantile(upper, bin_ids, quantiles, 4, {},
                                               kChromosomes),
               kCountNone);
}

TEST_CASE("nan_to_num and savetxt_line spell values the way numpy does") {
    std::vector<double> values = {kNan, kInf, -kInf, 1.0};
    cm::nan_to_num_in_place(values);
    check_bits(values, {0x0.0p+0, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
                        0x1.0000000000000p+0});
    CHECK(cm::savetxt_line({kInf, -kNan, -0.0, 1.0 / 3.0, 1e-300, -2.5e17}) ==
          "inf nan -0.000000000000000000e+00 3.333333333333333148e-01 "
          "1.000000000000000025e-300 -2.500000000000000000e+17\n");
    CHECK(cm::savetxt_line({}) == "\n");
}

TEST_CASE("quantile_boundaries follows the two branches of main:186-194") {
    std::vector<cm::PcaRow> rows;
    for (const double value : {0.5, -1.0, 2.0, 0.25}) {
        cm::PcaRow row;
        row.chrom = "chrX";
        row.pc1 = static_cast<float>(value);
        rows.push_back(row);
    }
    check_bits(cm::quantile_boundaries(rows, 3, 0.0), {-1.0, 0.375, 2.0});
    CHECK_THROWS_AS((void)cm::quantile_boundaries(rows, 1, 0.0), std::domain_error);
    CHECK(cm::quantile_boundaries(rows, 1, 5.0).size() == 1);
    CHECK(cm::quantile_boundaries(rows, 0, 0.0).empty());
    CHECK_THROWS_AS((void)cm::quantile_boundaries(rows, 30, 150.0), std::invalid_argument);
}
