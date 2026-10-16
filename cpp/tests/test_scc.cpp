// Unit tests for hicx::scc, the HiCRep stratum-adjusted correlation
// coefficient. The core cases here are small enough to check by hand; the
// smoothing (--sccSmoothWindow) and downsampling paths are instead validated
// end to end against the real hicrep 0.2.6 Python package on real GM12878
// replicate .cool files (see the hicCorrelate --method scc CLI test), since
// there is no convenient way to hand-verify a 2D edge-corrected box filter.

#include <doctest/doctest.h>

#include <cmath>

#include "hicx/scc.hpp"

using hicx::scc::ChromosomeEntries;
using hicx::scc::Options;
using hicx::scc::chromosome_scc;

TEST_CASE("chromosome_scc: identical shape scaled by a constant factor is a perfect match") {
    // Both matrices have only diagonal d=1 populated, with matrix 2 being
    // exactly half of matrix 1. After dividing each by its own total
    // contacts, the two normalized diagonals are bit-identical, so the
    // single diagonal with data has a weighted correlation of 1, and it is
    // the only diagonal that contributes any weight.
    ChromosomeEntries e1;
    e1.nbins = 4;
    e1.rows = {0, 1, 2};
    e1.cols = {1, 2, 3};
    e1.values = {2.0, 4.0, 6.0};

    ChromosomeEntries e2;
    e2.nbins = 4;
    e2.rows = {0, 1, 2};
    e2.cols = {1, 2, 3};
    e2.values = {1.0, 2.0, 3.0};

    Options opts;
    const double result = chromosome_scc(e1, e2, 12.0, 6.0, 1000, opts);
    CHECK(result == doctest::Approx(1.0));
}

TEST_CASE("chromosome_scc: anti-correlated diagonals score negative") {
    ChromosomeEntries e1;
    e1.nbins = 5;
    e1.rows = {0, 1, 2, 3};
    e1.cols = {1, 2, 3, 4};
    e1.values = {1.0, 2.0, 3.0, 4.0};

    ChromosomeEntries e2;
    e2.nbins = 5;
    e2.rows = {0, 1, 2, 3};
    e2.cols = {1, 2, 3, 4};
    e2.values = {4.0, 3.0, 2.0, 1.0};

    Options opts;
    const double result = chromosome_scc(e1, e2, 10.0, 10.0, 1000, opts);
    CHECK(result == doctest::Approx(-1.0));
}

TEST_CASE("chromosome_scc: max_distance_bp drops diagonals beyond the cutoff") {
    // A single far-away diagonal (d=3) outside a 2-bin cutoff must not
    // contribute: with it excluded, every remaining diagonal is empty on
    // both sides, so the result is NaN (no diagonal has two or more samples
    // with nonzero weight).
    ChromosomeEntries e1;
    e1.nbins = 5;
    e1.rows = {0};
    e1.cols = {3};
    e1.values = {1.0};

    ChromosomeEntries e2 = e1;
    e2.values = {5.0};

    Options opts;
    opts.max_distance_bp = 2000;  // bin_size 1000 -> d_max = 2000/1000+1 = 3, n_diags=min(3,5)=3
    const double result = chromosome_scc(e1, e2, 1.0, 5.0, 1000, opts);
    CHECK(std::isnan(result));
}

TEST_CASE("chromosome_scc: mismatched chromosome sizes are refused with NaN") {
    ChromosomeEntries e1;
    e1.nbins = 4;
    ChromosomeEntries e2;
    e2.nbins = 5;
    Options opts;
    CHECK(std::isnan(chromosome_scc(e1, e2, 1.0, 1.0, 1000, opts)));
}
