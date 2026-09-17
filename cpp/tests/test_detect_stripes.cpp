// Unit tests for hicDetectStripes' computational core
// (tools/stripes_impl.hpp), the only tier 9 tool (PLAN.md 9.3) with no Python
// reference to pin values against. These tests check the band construction
// against hand-computed values and the end-to-end pipeline against a
// synthetic matrix with stripes planted by construction, not against a
// reference implementation.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/stripes_impl.hpp"
#include "hicx/sparse_matrix.hpp"

using hicx::stripes::Band;
using hicx::stripes::Candidate;
using hicx::stripes::DetectOptions;
using hicx::stripes::RunningSums;

namespace {

hicx::CsrMatrix upper_triangle_from_dense(const std::vector<std::vector<double>>& dense) {
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    for (std::size_t r = 0; r < dense.size(); ++r) {
        for (std::size_t c = r + 1; c < dense[r].size(); ++c) {
            if (dense[r][c] != 0.0) {
                rows.push_back(static_cast<std::int32_t>(r));
                cols.push_back(static_cast<std::int32_t>(c));
                values.push_back(dense[r][c]);
            }
        }
    }
    hicx::CsrMatrix matrix = hicx::CsrMatrix::from_coo(
        static_cast<std::int64_t>(dense.size()), static_cast<std::int64_t>(dense[0].size()),
        rows, cols, std::move(values), "float64");
    matrix.set_symmetry(hicx::Symmetry::Full);
    return matrix;
}

}  // namespace

TEST_CASE("build_horizontal_band reads the right pixels and normalises by distance") {
    //      0  1  2  3
    // 0    .  4  2  0
    // 1    .  .  6  0
    // 2    .  .  .  8
    // 3    .  .  .  .
    // distance 1: (0,1)=4, (1,2)=6, (2,3)=8 -> expected = 6
    // distance 2: (0,2)=2, (1,3)=0          -> expected = 1
    // distance 3: (0,3)=0                    -> expected = 0
    const hicx::CsrMatrix matrix =
        upper_triangle_from_dense({{0, 4, 2, 0}, {0, 0, 6, 0}, {0, 0, 0, 8}, {0, 0, 0, 0}});
    std::vector<double> expected;
    const Band band = hicx::stripes::build_horizontal_band(matrix, 4, 3, &expected);

    CHECK(expected[0] == doctest::Approx(6.0));
    CHECK(expected[1] == doctest::Approx(1.0));
    CHECK(expected[2] == doctest::Approx(0.0));

    CHECK(band.raw_at(0, 1) == doctest::Approx(4.0));
    CHECK(band.obs_exp_at(0, 1) == doctest::Approx(4.0 / 6.0));
    CHECK(band.raw_at(1, 1) == doctest::Approx(6.0));
    CHECK(band.obs_exp_at(1, 1) == doctest::Approx(1.0));
    CHECK(band.raw_at(0, 2) == doctest::Approx(2.0));
    CHECK(band.obs_exp_at(0, 2) == doctest::Approx(2.0));
    // distance 3's expected value is 0 (only one position, unstored), so
    // obs_exp is defined as 0 rather than dividing by zero.
    CHECK(band.obs_exp_at(0, 3) == doctest::Approx(0.0));
}

TEST_CASE("build_vertical_band mirrors the horizontal band by the (anchor, distance) shift") {
    const hicx::CsrMatrix matrix =
        upper_triangle_from_dense({{0, 4, 2, 1}, {0, 0, 6, 3}, {0, 0, 0, 8}, {0, 0, 0, 0}});
    const Band horizontal = hicx::stripes::build_horizontal_band(matrix, 4, 3);
    const Band vertical = hicx::stripes::build_vertical_band(horizontal);

    // vertical.raw_at(c, d) == horizontal.raw_at(c - d, d)
    CHECK(vertical.raw_at(3, 1) == doctest::Approx(horizontal.raw_at(2, 1)));
    CHECK(vertical.raw_at(3, 2) == doctest::Approx(horizontal.raw_at(1, 2)));
    CHECK(vertical.raw_at(3, 3) == doctest::Approx(horizontal.raw_at(0, 3)));
    // No pixel exists above the top row: anchor 1 at distance 2 would need
    // row -1.
    CHECK(vertical.raw_at(1, 2) == doctest::Approx(0.0));
}

TEST_CASE("build_running_sums accumulates obs_exp and raw along the distance axis") {
    const hicx::CsrMatrix matrix =
        upper_triangle_from_dense({{0, 2, 2, 2}, {0, 0, 2, 2}, {0, 0, 0, 2}, {0, 0, 0, 0}});
    const Band band = hicx::stripes::build_horizontal_band(matrix, 4, 3);
    const RunningSums sums = hicx::stripes::build_running_sums(band);

    // A uniform matrix has obs_exp 1.0 everywhere it is defined, so the mean
    // at any length is 1.0 for anchor 0 (which has all three distances).
    CHECK(sums.obs_exp_mean(0, 1) == doctest::Approx(1.0));
    CHECK(sums.obs_exp_mean(0, 3) == doctest::Approx(1.0));
    CHECK(sums.raw_mean(0, 1) == doctest::Approx(2.0));
}

namespace {

// A 60-bin synthetic chromosome: a flat background of 3 counts at every
// distance up to 35 bins, with a horizontal stripe planted at anchor row 10
// (columns 11..30, 20 bins long) and a vertical stripe at anchor column 50
// (rows 30..49, 20 bins long), each about 5-fold enriched over the
// background. The two stripes do not overlap each other's footprint or
// background window.
hicx::CsrMatrix synthetic_matrix(std::int64_t n, std::int64_t max_distance,
                                 double background, double stripe_value,
                                 std::int64_t horizontal_anchor, std::int64_t horizontal_length,
                                 std::int64_t vertical_anchor, std::int64_t vertical_length) {
    std::vector<std::vector<double>> dense(static_cast<std::size_t>(n),
                                           std::vector<double>(static_cast<std::size_t>(n), 0.0));
    // Small Poisson-like noise around `background`, seeded for
    // reproducibility, so the background window has nonzero variance: a
    // perfectly flat background (std 0) makes every candidate length tie on
    // z-score, which is not representative of real, noisy contact data. A
    // fixed pattern (for example a modulo function of row and column) was
    // tried first and rejected: its period aliases with the length grid and
    // produces the same kind of degenerate ties.
    std::mt19937 rng(20260915);
    std::poisson_distribution<int> noise(background);
    for (std::int64_t row = 0; row < n; ++row) {
        for (std::int64_t col = row + 1; col < n && col - row <= max_distance; ++col) {
            dense[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] =
                static_cast<double>(noise(rng));
        }
    }
    for (std::int64_t d = 1; d <= horizontal_length; ++d) {
        const std::int64_t col = horizontal_anchor + d;
        if (col < n) {
            dense[static_cast<std::size_t>(horizontal_anchor)][static_cast<std::size_t>(col)] =
                stripe_value;
        }
    }
    for (std::int64_t d = 1; d <= vertical_length; ++d) {
        const std::int64_t row = vertical_anchor - d;
        if (row >= 0) {
            dense[static_cast<std::size_t>(row)][static_cast<std::size_t>(vertical_anchor)] =
                stripe_value;
        }
    }
    return upper_triangle_from_dense(dense);
}

}  // namespace

TEST_CASE("the end-to-end pipeline recovers a planted horizontal and vertical stripe") {
    constexpr std::int64_t n = 60;
    constexpr std::int64_t max_distance = 35;
    const hicx::CsrMatrix matrix =
        synthetic_matrix(n, max_distance, /*background=*/3.0, /*stripe_value=*/15.0,
                         /*horizontal_anchor=*/10, /*horizontal_length=*/20,
                         /*vertical_anchor=*/50, /*vertical_length=*/20);

    const Band horizontal = hicx::stripes::build_horizontal_band(matrix, n, max_distance);
    const Band vertical = hicx::stripes::build_vertical_band(horizontal);
    const RunningSums horizontal_sums = hicx::stripes::build_running_sums(horizontal);
    const RunningSums vertical_sums = hicx::stripes::build_running_sums(vertical);

    DetectOptions options;
    options.length_grid_bins = {5, 10, 15, 20, 25};
    options.background_window_bins = 5;
    options.background_gap_bins = 1;
    options.min_obs_exp = 1.5;
    options.preselect_z = 1.5;
    options.min_raw_count = 1.0;
    options.merge_window_bins = 2;
    options.fdr_q = 0.2;

    std::vector<Candidate> candidates = hicx::stripes::preselect(
        horizontal, vertical, horizontal_sums, vertical_sums, options, /*threads=*/1);
    REQUIRE(!candidates.empty());
    candidates = hicx::stripes::suppress_non_maximal(std::move(candidates), options.merge_window_bins);
    hicx::stripes::compute_pvalues(candidates, horizontal, vertical, options, /*threads=*/1);
    const std::vector<Candidate> kept = hicx::stripes::apply_fdr(candidates, options.fdr_q);

    bool found_horizontal = false;
    bool found_vertical = false;
    for (const Candidate& candidate : kept) {
        if (!candidate.vertical && candidate.anchor == 10) {
            found_horizontal = true;
            CHECK(candidate.length_bins >= 15);
            CHECK(candidate.enrichment > 2.0);
            CHECK(candidate.qvalue <= options.fdr_q);
        }
        if (candidate.vertical && candidate.anchor == 50) {
            found_vertical = true;
            CHECK(candidate.length_bins >= 15);
            CHECK(candidate.enrichment > 2.0);
        }
    }
    CHECK(found_horizontal);
    CHECK(found_vertical);
}

TEST_CASE("a flat matrix with no planted stripe yields no calls after FDR") {
    // The preselection z-test is deliberately cheap and uncorrected
    // (hicx::stripes::preselect docstring), so pure noise routinely clears a
    // z of 1.5 at a few of the many (anchor, orientation, length) triples
    // tested; that is what stage 3's rank-sum test and stage 4's
    // Benjamini-Hochberg FDR exist to remove. This test therefore checks the
    // full pipeline's output, not the raw preselection.
    constexpr std::int64_t n = 40;
    constexpr std::int64_t max_distance = 25;
    const hicx::CsrMatrix matrix =
        synthetic_matrix(n, max_distance, /*background=*/3.0, /*stripe_value=*/3.0,
                         /*horizontal_anchor=*/0, /*horizontal_length=*/0,
                         /*vertical_anchor=*/0, /*vertical_length=*/0);
    const Band horizontal = hicx::stripes::build_horizontal_band(matrix, n, max_distance);
    const Band vertical = hicx::stripes::build_vertical_band(horizontal);
    const RunningSums horizontal_sums = hicx::stripes::build_running_sums(horizontal);
    const RunningSums vertical_sums = hicx::stripes::build_running_sums(vertical);

    DetectOptions options;
    options.length_grid_bins = {5, 10, 15, 20};
    options.background_window_bins = 5;
    options.background_gap_bins = 1;
    options.min_obs_exp = 1.5;
    options.preselect_z = 1.5;
    options.min_raw_count = 1.0;
    options.merge_window_bins = 2;
    options.fdr_q = 0.05;

    std::vector<Candidate> candidates = hicx::stripes::preselect(
        horizontal, vertical, horizontal_sums, vertical_sums, options, /*threads=*/1);
    candidates = hicx::stripes::suppress_non_maximal(std::move(candidates), options.merge_window_bins);
    hicx::stripes::compute_pvalues(candidates, horizontal, vertical, options, /*threads=*/1);
    const std::vector<Candidate> kept = hicx::stripes::apply_fdr(candidates, options.fdr_q);
    CHECK(kept.empty());
}

TEST_CASE("suppress_non_maximal keeps only the strongest anchor within the merge window") {
    std::vector<Candidate> candidates;
    Candidate a;
    a.anchor = 10;
    a.vertical = false;
    a.zscore = 3.0;
    a.pvalue = 0.01;
    candidates.push_back(a);
    Candidate b = a;
    b.anchor = 12;
    b.zscore = 5.0;
    b.pvalue = 0.001;
    candidates.push_back(b);
    Candidate c = a;
    c.anchor = 40;
    c.zscore = 1.0;
    candidates.push_back(c);

    const std::vector<Candidate> kept = hicx::stripes::suppress_non_maximal(candidates, 5);
    REQUIRE(kept.size() == 2);
    CHECK(kept[0].anchor == 12);
    CHECK(kept[1].anchor == 40);
}

TEST_CASE("apply_fdr keeps only candidates at or below the requested q-value") {
    std::vector<Candidate> candidates;
    for (const double p : {0.001, 0.01, 0.5, 0.9}) {
        Candidate candidate;
        candidate.pvalue = p;
        candidates.push_back(candidate);
    }
    const std::vector<Candidate> kept = hicx::stripes::apply_fdr(candidates, 0.05);
    CHECK(kept.size() <= candidates.size());
    for (const Candidate& candidate : kept) {
        CHECK(candidate.qvalue <= 0.05);
    }
}
