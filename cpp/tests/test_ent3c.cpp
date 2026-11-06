// Unit tests for hicx::ent3c. The expected values of the edge-case windows
// are ENT3C 2.2.2's own entropies of those windows (the 2.5 Mb GM12878 test
// matrix, chromosome Y); the agreement on whole matrices is checked by
// tests/cli/hicCorrelate_ent3c_cli.sh.

#include <doctest/doctest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "hicx/ent3c.hpp"

using hicx::ent3c::ChromosomeEntries;
using hicx::ent3c::Options;

TEST_CASE("windows: n from the chromosome split, Phi from the shift") {
    Options options;  // split 7, shift 1, at most 1000 windows
    const hicx::ent3c::Windows w = hicx::ent3c::windows(457, options);
    CHECK(w.size == 65);  // floor(457 / 7 + 0.5)
    CHECK(w.count == 393);
    CHECK(w.shift == 1);
}

TEST_CASE("windows: the shift grows until there are at most max_windows") {
    Options options;
    options.max_windows = 50;
    const hicx::ent3c::Windows w = hicx::ent3c::windows(457, options);
    CHECK(w.shift == 8);  // shift 7 would give 57
    CHECK(w.count == 50);
    options.max_windows = 0;  // no limit
    CHECK(hicx::ent3c::windows(457, options).count == 393);
}

TEST_CASE("windows: a fixed size, a size above the bin count, invalid options") {
    Options options;
    options.submatrix_size = 40;
    CHECK(hicx::ent3c::windows(457, options).count == 418);
    options.submatrix_size = 500;
    CHECK(hicx::ent3c::windows(457, options).count == 0);
    options.submatrix_size = -1;
    CHECK_THROWS_AS(hicx::ent3c::windows(457, options), std::invalid_argument);
    options = Options();
    options.window_shift = 0;
    CHECK_THROWS_AS(hicx::ent3c::windows(457, options), std::invalid_argument);
    options = Options();
    options.chromosome_split = 0;
    CHECK_THROWS_AS(hicx::ent3c::windows(457, options), std::invalid_argument);
}

TEST_CASE("window_entropy: ENT3C's values for degenerate windows") {
    // Both columns constant after the fill but kept (their logs do not sum
    // to 0): P is the identity, rho = I / 2, entropy ln 2.
    CHECK(hicx::ent3c::window_entropy({2.0, 0.0, 0.0, 181.0}, 2) ==
          doctest::Approx(0.6931471805599453).epsilon(1e-15));
    // The second column holds only 1 and 0: its logs sum to 0, it leaves,
    // and one bin cannot be correlated.
    CHECK(std::isnan(hicx::ent3c::window_entropy({4.0, 1.0, 1.0, 0.0}, 2)));
    // The first column is empty.
    CHECK(std::isnan(hicx::ent3c::window_entropy({0.0, 0.0, 0.0, 322.0}, 2)));
    CHECK(std::isnan(hicx::ent3c::window_entropy({0.0, 0.0, 0.0, 0.0}, 2)));
}

TEST_CASE("window_entropy: a window of equal counts has entropy ln k") {
    // Every column constant and kept: P = I, rho = I / 4.
    const std::vector<double> counts(16, 3.0);
    CHECK(hicx::ent3c::window_entropy(counts, 4) == doctest::Approx(std::log(4.0)).epsilon(1e-14));
}

TEST_CASE("similarity: windows without entropy are left out") {
    const double nan = std::nan("");
    CHECK(hicx::ent3c::similarity({1.0, 2.0, nan, 3.0}, {2.0, 4.0, 5.0, 6.0}) ==
          doctest::Approx(1.0).epsilon(1e-15));
    CHECK(hicx::ent3c::similarity({1.0, 2.0, 3.0}, {3.0, 2.0, 1.0}) ==
          doctest::Approx(-1.0).epsilon(1e-15));
    CHECK(std::isnan(hicx::ent3c::similarity({1.0, 1.0, 1.0}, {1.0, 2.0, 3.0})));  // constant
    CHECK(std::isnan(hicx::ent3c::similarity({1.0, nan}, {1.0, 2.0})));            // one window
    CHECK(std::isnan(hicx::ent3c::similarity({1.0, 2.0}, {1.0, 2.0, 3.0})));       // lengths
}

TEST_CASE("empty_bins and entropy_signal: both triangles, thread independent") {
    // A 30-bin chromosome with a banded, varying signal; bin 7 is empty.
    ChromosomeEntries e;
    e.nbins = 30;
    for (std::int32_t i = 0; i < 30; ++i) {
        for (std::int32_t d = 0; d < 8 && i + d < 30; ++d) {
            if (i == 7 || i + d == 7) {
                continue;
            }
            e.rows.push_back(i);
            e.cols.push_back(i + d);
            e.values.push_back(static_cast<double>(1 + ((i * 7 + d * 3) % 11) * (8 - d)));
        }
    }
    const std::vector<bool> empty = hicx::ent3c::empty_bins(e);
    for (std::size_t i = 0; i < empty.size(); ++i) {
        CHECK(empty[i] == (i == 7));
    }
    std::vector<std::int64_t> keep;
    for (std::int64_t i = 0; i < 30; ++i) {
        if (!empty[static_cast<std::size_t>(i)]) {
            keep.push_back(i);
        }
    }
    Options options;
    options.submatrix_size = 10;
    const hicx::ent3c::Windows w = hicx::ent3c::windows(static_cast<std::int64_t>(keep.size()), options);
    REQUIRE(w.count == 20);
    const std::vector<double> one = hicx::ent3c::entropy_signal(e, keep, w, 1);
    const std::vector<double> four = hicx::ent3c::entropy_signal(e, keep, w, 4);
    REQUIRE(one.size() == 20);
    for (std::size_t i = 0; i < one.size(); ++i) {
        CHECK(std::isfinite(one[i]));
        CHECK(one[i] == four[i]);  // bitwise
    }
    CHECK(hicx::ent3c::similarity(one, four) == doctest::Approx(1.0).epsilon(1e-15));
}
