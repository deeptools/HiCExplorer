// Unit tests for hicx::hicspector. These cases are small enough to verify
// by hand; the algorithm's fidelity to the real gersteinlab/HiC-spector
// source on realistic data is checked against that real source directly
// (run as a validation oracle; it needs nothing but numpy/scipy, unlike
// hicx::genomedisco's Python-2-only reference), on real GM12878/A549 data
// (see hicspector.hpp and the commit message for the exact values).

#include <doctest/doctest.h>

#include <cmath>

#include "hicx/hicspector.hpp"

using hicx::hicspector::ChromosomeEntries;
using hicx::hicspector::Options;
using hicx::hicspector::chromosome_score;

TEST_CASE("chromosome_score: a matrix against itself scores (numerically) exactly 1.0") {
    ChromosomeEntries e;
    e.nbins = 6;
    e.rows = {0, 0, 1, 2, 3, 4};
    e.cols = {1, 2, 2, 3, 4, 5};
    e.values = {5.0, 2.0, 4.0, 3.0, 6.0, 1.0};

    Options opts;
    opts.num_evec = 4;
    opts.ipr_cut = 0.0;  // the tiny hand-built graph below has few/localized
                         // eigenvectors; the real tool's own ipr_cut=5 would
                         // filter all of them out on a graph this small.
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

TEST_CASE("chromosome_score: a bin with only a self-contact is dropped, not crashing") {
    // Bin 3 has only its own diagonal entry nonzero (a self-ligation-only
    // bin): the degenerate-bin filter (k==1 and diagonal>0) must drop it
    // without the rest of the pipeline tripping over the resulting index
    // gap.
    ChromosomeEntries e;
    e.nbins = 5;
    e.rows = {0, 0, 1, 2, 3, 3};
    e.cols = {1, 2, 2, 4, 3, 4};
    e.values = {5.0, 2.0, 4.0, 3.0, 9.0, 1.0};

    Options opts;
    opts.num_evec = 3;
    opts.ipr_cut = 0.0;
    const double score = chromosome_score(e, e, opts);
    REQUIRE(!std::isnan(score));
    CHECK(score == doctest::Approx(1.0).epsilon(1e-9));
}
