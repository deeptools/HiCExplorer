// Unit tests for hicx::genomedisco. These cases are small enough to verify
// by hand; the algorithm's fidelity to the real kundajelab/genomedisco
// source on realistic data is instead checked against a hand-transcribed
// python3 reproduction of that source on real GM12878/A549 data (see
// genomedisco.hpp and the commit message for the exact values), since the
// real package is Python 2 only and does not install in this environment.

#include <doctest/doctest.h>

#include <cmath>

#include "hicx/genomedisco.hpp"

using hicx::genomedisco::ChromosomeEntries;
using hicx::genomedisco::Options;
using hicx::genomedisco::chromosome_score;

TEST_CASE("chromosome_score: a matrix against itself scores a perfect 1.0") {
    // A 4-node ring with varied weights, rich enough that sqrtvc/to_transition
    // leave no degenerate (all-equal) rows, so self-comparison truly reaches
    // the maximum score rather than being diluted the way hicx::scc's
    // self-vs-self can be on a sparse/fragmented fixture.
    ChromosomeEntries e;
    e.nbins = 4;
    e.rows = {0, 0, 1, 2};
    e.cols = {1, 2, 2, 3};
    e.values = {3.0, 1.0, 2.0, 5.0};

    Options opts;
    opts.subsample = false;
    const double score = chromosome_score(e, e, opts);
    CHECK(score == doctest::Approx(1.0));
}

TEST_CASE("chromosome_score: disjoint contact patterns score low") {
    ChromosomeEntries e1;
    e1.nbins = 5;
    e1.rows = {0, 1, 2};
    e1.cols = {1, 2, 3};
    e1.values = {5.0, 5.0, 5.0};

    ChromosomeEntries e2;
    e2.nbins = 5;
    e2.rows = {0, 2, 3};
    e2.cols = {4, 4, 4};
    e2.values = {5.0, 5.0, 5.0};

    Options opts;
    opts.subsample = false;
    const double score = chromosome_score(e1, e2, opts);
    REQUIRE(!std::isnan(score));
    CHECK(score < 0.5);
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
    opts.subsample = false;
    CHECK(std::isnan(chromosome_score(e1, e2, opts)));
}

TEST_CASE("chromosome_score: a bin with exactly one contact is not silently dropped") {
    // Regression test for an off-by-one in the internal CSR builder that
    // could drop a row's only entry when that row immediately followed a
    // row that already had entries (hicx::genomedisco's own CSR assembly,
    // not hicx::CsrMatrix). Comparing a matrix against a copy with that one
    // entry zeroed out must show a real difference; a dropped row would
    // instead make the two sides look identical where they are not.
    ChromosomeEntries full;
    full.nbins = 5;
    full.rows = {0, 0, 1, 2, 3};
    full.cols = {1, 2, 2, 3, 4};
    full.values = {4.0, 3.0, 2.0, 6.0, 5.0};  // row 2 has exactly one entry

    ChromosomeEntries missing_row2;
    missing_row2.nbins = 5;
    missing_row2.rows = {0, 0, 3};
    missing_row2.cols = {1, 2, 4};
    missing_row2.values = {4.0, 3.0, 5.0};

    Options opts;
    opts.subsample = false;
    const double self_score = chromosome_score(full, full, opts);
    const double diff_score = chromosome_score(full, missing_row2, opts);
    CHECK(self_score == doctest::Approx(1.0));
    REQUIRE(!std::isnan(diff_score));
    CHECK(diff_score < self_score);
}
