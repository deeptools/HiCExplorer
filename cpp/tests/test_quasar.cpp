// Unit tests for hicx::quasar. These cases are small enough to verify by
// hand; the algorithm's fidelity to the real hifive 1.5.7 source on
// realistic data is checked against the real compiled
// hifive.libraries._quasar Cython extension, called directly under Python
// 2.7 (not read and guessed at) on real GM12878/A549 data (see
// quasar.hpp and the commit message for the exact reference values).

#include <doctest/doctest.h>

#include <cmath>

#include "hicx/quasar.hpp"

using hicx::quasar::ChromosomeEntries;
using hicx::quasar::Options;
using hicx::quasar::chromosome_score;

TEST_CASE("chromosome_score: a matrix against itself scores (numerically) exactly 1.0") {
    // A 12-bin chain dense enough, with a small width, that every bin has
    // several valid neighbours and the correlation/weight machinery has
    // real variance to work with.
    ChromosomeEntries e;
    e.nbins = 12;
    for (std::int32_t i = 0; i < 11; ++i) {
        e.rows.push_back(i);
        e.cols.push_back(i + 1);
        e.values.push_back(static_cast<double>(5 + (i % 4)));
    }
    for (std::int32_t i = 0; i < 10; ++i) {
        e.rows.push_back(i);
        e.cols.push_back(i + 2);
        e.values.push_back(static_cast<double>(3 + (i % 3)));
    }

    Options opts;
    opts.width = 4;
    const double score = chromosome_score(e, e, opts);
    REQUIRE(!std::isnan(score));
    CHECK(score == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("chromosome_score: mismatched chromosome sizes are refused with NaN") {
    ChromosomeEntries e1;
    e1.nbins = 4;
    ChromosomeEntries e2;
    e2.nbins = 5;
    Options opts;
    CHECK(std::isnan(chromosome_score(e1, e2, opts)));
}

TEST_CASE("chromosome_score: an empty chromosome is refused with NaN") {
    ChromosomeEntries e1;
    e1.nbins = 4;
    ChromosomeEntries e2;
    e2.nbins = 4;
    e2.rows = {0};
    e2.cols = {1};
    e2.values = {1.0};
    Options opts;
    CHECK(std::isnan(chromosome_score(e1, e2, opts)));
}

TEST_CASE("chromosome_score: a chromosome with no contacts within width is refused with NaN") {
    // Every bin's only neighbour is more than --width away, so the
    // iterative filter_bins fixed point drops every bin to invalid.
    ChromosomeEntries e;
    e.nbins = 10;
    e.rows = {0};
    e.cols = {9};
    e.values = {5.0};
    Options opts;
    opts.width = 2;
    CHECK(std::isnan(chromosome_score(e, e, opts)));
}

