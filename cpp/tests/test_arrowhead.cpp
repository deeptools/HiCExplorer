// Unit tests for hicArrowhead's computational core (tools/arrowhead_impl.hpp):
// a faithful C++ port of Juicer/Juicebox's Arrowhead contact-domain caller.
// There is no Python reference to pin values against (and no local Java/
// juicer_tools environment to generate one from the real tool), so these
// check the low-level primitive against a hand-computed value and the
// end-to-end pipeline against synthetic matrices with a known structure
// (uniform: no domains; two strong diagonal blocks: at least one domain).

#include <cmath>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/arrowhead_impl.hpp"

using hicx::arrowhead::Matrix;
using hicx::arrowhead::Options;

TEST_CASE("directionality_index_upstream matches a hand-computed row") {
    // n = 30, gap = 7, i = 15: window = min(30 - 22, 15 - 7) = 8, so
    // off = 0, 1 are computed (off <= window - gap = 1).
    //   off=0: left = i-gap-off = 8,  right = i+gap+off = 22
    //   off=1: left = i-gap-off = 7,  right = i+gap+off = 23
    Matrix observed(30);
    observed.at(15, 8) = 6.0;
    observed.at(15, 22) = 2.0;   // (6-2)/(6+2) = 0.5
    observed.at(15, 7) = 3.0;
    observed.at(15, 23) = 9.0;   // (3-9)/(3+9) = -0.5

    const Matrix d = hicx::arrowhead::directionality_index_upstream(observed, 7);
    CHECK(d.at(15, 22) == doctest::Approx(0.5));
    CHECK(d.at(15, 23) == doctest::Approx(-0.5));
    // Every other entry on row 15 outside those two columns stays zero
    // (no data was placed there, so a=b=0 would be the degenerate 0/0 case,
    // but off only reaches those two columns here).
    CHECK(d.at(15, 21) == doctest::Approx(0.0));
}

TEST_CASE("directionality_index_upstream: rows too close to the matrix edge are untouched") {
    Matrix observed(10);
    for (double& v : observed.v) v = 1.0;
    // n=10, gap=7: window = min(10-(i+7), i-7); for every i in [0,10), at
    // least one side has less than gap=7 room, so no row is touched.
    const Matrix d = hicx::arrowhead::directionality_index_upstream(observed, 7);
    for (const double v : d.v) {
        CHECK(v == doctest::Approx(0.0));
    }
}

TEST_CASE("call_domains: a uniform matrix has no detectable contact domains") {
    const int n = 200;
    const auto observed_window = [&](int lim_start, int lim_end) {
        Matrix m(lim_end - lim_start);
        for (double& v : m.v) v = 5.0;
        return m;
    };
    Options opts;
    opts.matrix_width = 100;
    const std::vector<hicx::arrowhead::HighScore> domains =
        hicx::arrowhead::call_domains(n, opts, observed_window);
    CHECK(domains.empty());
}

TEST_CASE("call_domains: two strong diagonal blocks yield at least one domain") {
    // Two 40-bin blocks (rows/cols [0,40) and [40,80)) with strong intra-
    // block contacts and weak inter-block contacts, a textbook planted
    // contact-domain pattern.
    const int n = 80;
    std::vector<double> dense(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 1.0);
    auto set = [&](int i, int j, double v) { dense[static_cast<std::size_t>(i) * n + j] = v; };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const bool same_block = (i < 40) == (j < 40);
            // A little deterministic variation on top of the flat blocks:
            // a perfectly uniform block makes the directionality index's
            // upstream/downstream sign sums collapse to an exact zero on
            // one side (an artifact of the step-function symmetry, not
            // something real contact data would produce), which a sign
            // threshold's strict "<" then rejects at any positive
            // threshold including floating-point noise near zero.
            const double noise = same_block ? static_cast<double>((i + j) % 7)
                                            : static_cast<double>((i * 3 + j) % 3) * 0.1;
            set(i, j, (same_block ? 50.0 : 1.0) + noise);
        }
    }

    const auto observed_window = [&](int lim_start, int lim_end) {
        const int size = lim_end - lim_start;
        Matrix m(size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                m.at(i, j) = dense[static_cast<std::size_t>(lim_start + i) * n + (lim_start + j)];
            }
        }
        return m;
    };

    Options opts;
    opts.matrix_width = 80;  // one window covers the whole synthetic chromosome
    // The default --minBlockSize (60 bins) assumes real TAD-scale windows
    // (Juicer's own default --matrixWidth is 2000); the planted blocks here
    // are 40 bins wide, so the size filter is relaxed to match this test's
    // own scale.
    opts.min_block_size = 10;
    // The algorithm's sign gate requires BOTH the upper- and lower-triangle
    // sign means at a candidate corner to clear the threshold. A flat,
    // two-value planted block is a clean step function, so its two
    // strongest corners (the ones actually detected) each have one side's
    // sign mean collapse near zero by construction -- real, noisier contact
    // data does not do this as cleanly. The sign gate is relaxed here so
    // the test exercises the variance/connected-component/merge machinery
    // without depending on a synthetic pattern's own sign-symmetry quirks.
    opts.max_low_sign_threshold = -1.0;
    opts.min_low_sign_threshold = -1.0;
    const std::vector<hicx::arrowhead::HighScore> domains =
        hicx::arrowhead::call_domains(n, opts, observed_window);
    REQUIRE(!domains.empty());
    // Sorted descending by sort_value() (uVarScore + lVarScore): the top
    // call should be a non-trivial block, not a single-bin degenerate one.
    CHECK(domains.front().width() > 1);
}
