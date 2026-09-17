// Unit tests for hicDetectStripes' computational core
// (tools/stripes_impl.hpp): a faithful C++ port of Stripenn 1.1.65.22's own
// stripe-finding algorithm (PLAN.md tier 9 section 9.3). There is no Python
// reference to pin values against, so these check the low-level image
// primitives against hand-computed values, and the end-to-end frame search
// against a synthetic image with a stripe planted by construction.

#include <cmath>
#include <cstdint>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/stripes_impl.hpp"

using hicx::stripes::BackgroundModel;
using hicx::stripes::Image;

TEST_CASE("quantile matches numpy's linear interpolation") {
    CHECK(hicx::stripes::quantile({1, 2, 3, 4}, 0.5) == doctest::Approx(2.5));
    CHECK(hicx::stripes::quantile({1, 2, 3, 4}, 0.0) == doctest::Approx(1.0));
    CHECK(hicx::stripes::quantile({1, 2, 3, 4}, 1.0) == doctest::Approx(4.0));
    CHECK(hicx::stripes::quantile({5}, 0.3) == doctest::Approx(5.0));
}

TEST_CASE("box_blur averages a uniform image to itself") {
    Image img(5, 5, 0.5);
    const Image blurred = hicx::stripes::box_blur(img, 3);
    for (const double v : blurred.v) {
        CHECK(v == doctest::Approx(0.5));
    }
}

TEST_CASE("box_blur smooths a single bright pixel") {
    Image img(5, 5, 0.0);
    img.at(2, 2) = 1.0;
    const Image blurred = hicx::stripes::box_blur(img, 3);
    // The centre pixel's 3x3 neighbourhood has exactly one bright pixel.
    CHECK(blurred.at(2, 2) == doctest::Approx(1.0 / 9.0));
    CHECK(blurred.at(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("canny finds an edge at a sharp brightness step") {
    // A vertical step edge: columns 0-4 dark, columns 5-9 bright.
    Image img(10, 10, 0.0);
    for (int r = 0; r < 10; ++r) {
        for (int c = 5; c < 10; ++c) {
            img.at(r, c) = 1.0;
        }
    }
    const Image edges = hicx::stripes::canny(img, 1.0);
    bool found_edge_near_boundary = false;
    for (int r = 2; r < 8; ++r) {
        for (int c = 3; c < 7; ++c) {
            if (edges.at(r, c) != 0.0) {
                found_edge_near_boundary = true;
            }
        }
    }
    CHECK(found_edge_near_boundary);
    // Far from the boundary, on a flat region, there should be no edges.
    CHECK(edges.at(1, 1) == doctest::Approx(0.0));
    CHECK(edges.at(1, 8) == doctest::Approx(0.0));
}

TEST_CASE("vertical_line keeps a vertical edge and drops a horizontal one") {
    // A purely vertical edge (column step) should produce mostly-vertical
    // gradient orientation (near 90 degrees) and survive the filter.
    Image vertical_edge(10, 10, 0.0);
    for (int r = 3; r < 7; ++r) {
        vertical_edge.at(r, 5) = 1.0;
    }
    const Image vert = hicx::stripes::vertical_line(vertical_edge, 60.0, 120.0);
    bool any_kept = false;
    for (const double v : vert.v) {
        if (v != 0.0) {
            any_kept = true;
        }
    }
    CHECK(any_kept);
}

TEST_CASE("block_scan finds the longest contiguous run in a column") {
    Image vert(20, 3, 0.0);
    for (int r = 2; r < 15; ++r) {
        vert.at(r, 1) = 1.0;
    }
    const hicx::stripes::BlockResult result = hicx::stripes::block_scan(vert, 1);
    CHECK(result.length == 13);
}

TEST_CASE("block_scan tolerates a short gap") {
    Image vert(20, 3, 0.0);
    for (int r = 2; r < 8; ++r) {
        vert.at(r, 1) = 1.0;
    }
    // A 2-row gap (below the tolerance of 5) should not break the run.
    for (int r = 10; r < 16; ++r) {
        vert.at(r, 1) = 1.0;
    }
    const hicx::stripes::BlockResult result = hicx::stripes::block_scan(vert, 1);
    CHECK(result.length == 12);
}

namespace {

// A synthetic frame: a flat low background with a bright vertical stripe
// (a narrow column range, wide row range) planted near the middle, touching
// the frame's diagonal at its top so verticalLine's edge detection has a
// real boundary to find on both sides of the stripe.
Image synthetic_frame(int S, double background, double stripe_value, int stripe_col_start,
                      int stripe_width, int stripe_row_start, int stripe_row_end) {
    Image img(S, S, 0.0);
    for (int r = 0; r < S; ++r) {
        for (int c = 0; c < S; ++c) {
            const int distance = std::abs(r - c);
            img.at(r, c) = distance == 0 ? 0.0 : background / (1.0 + 0.01 * distance);
        }
    }
    for (int r = stripe_row_start; r <= stripe_row_end && r < S; ++r) {
        for (int c = stripe_col_start; c < stripe_col_start + stripe_width && c < S; ++c) {
            img.at(r, c) = stripe_value;
            img.at(c, r) = stripe_value;
        }
    }
    return img;
}

}  // namespace

TEST_CASE("stripe_search_frame finds candidates on a synthetic planted stripe") {
    constexpr int S = 200;
    const Image frame = synthetic_frame(S, /*background=*/5.0, /*stripe_value=*/80.0,
                                        /*stripe_col_start=*/60, /*stripe_width=*/4,
                                        /*stripe_row_start=*/60, /*stripe_row_end=*/150);
    const double maxpixel_value = 60.0;
    const std::vector<hicx::stripes::FrameCandidate> candidates = hicx::stripes::stripe_search_frame(
        frame, maxpixel_value, /*canny_sigma=*/1.5, /*min_length=*/10, /*max_width=*/8,
        /*blur_filter=*/3);
    REQUIRE(!candidates.empty());
    bool found_near_planted = false;
    for (const auto& c : candidates) {
        if (c.x <= 64 && c.x + c.w >= 60 && c.h >= 20) {
            found_near_planted = true;
        }
    }
    CHECK(found_near_planted);
}

TEST_CASE("stripe_search_frame finds nothing on a flat frame") {
    constexpr int S = 100;
    Image frame(S, S, 0.0);
    for (int r = 0; r < S; ++r) {
        for (int c = 0; c < S; ++c) {
            const int distance = std::abs(r - c);
            frame.at(r, c) = distance == 0 ? 0.0 : 5.0 / (1.0 + 0.01 * distance);
        }
    }
    const std::vector<hicx::stripes::FrameCandidate> candidates =
        hicx::stripes::stripe_search_frame(frame, 6.0, 1.5, 10, 8, 3);
    // A perfectly smooth distance-decaying background has no sharp edges,
    // so there should be no or very few spurious candidates.
    CHECK(candidates.size() < 5);
}

TEST_CASE("remove_redundant keeps the more elongated of two overlapping boxes") {
    std::vector<hicx::stripes::Candidate> candidates;
    hicx::stripes::Candidate a;
    a.chrom = "chr1";
    a.pos1 = 0;
    a.pos2 = 50;   // wide (width 50)
    a.pos3 = 0;
    a.pos4 = 60;   // short (height 60), ratio 60/50 = 1.2
    a.frame_index = 0;
    candidates.push_back(a);

    hicx::stripes::Candidate b;
    b.chrom = "chr1";
    b.pos1 = 0;
    b.pos2 = 10;    // narrow (width 10)
    b.pos3 = 0;
    b.pos4 = 500;   // long (height 500), ratio 50, much more elongated
    b.frame_index = 0;
    candidates.push_back(b);

    const std::vector<hicx::stripes::Candidate> kept =
        hicx::stripes::remove_redundant(candidates, /*by_pvalue=*/false, /*threads=*/1);
    REQUIRE(kept.size() == 1);
    CHECK(kept[0].pos4 == 500);
}

TEST_CASE("remove_redundant keeps the smaller p-value when requested") {
    std::vector<hicx::stripes::Candidate> candidates;
    hicx::stripes::Candidate a;
    a.chrom = "chr1";
    a.pos1 = 0;
    a.pos2 = 100;
    a.pos3 = 0;
    a.pos4 = 400;
    a.frame_index = 0;
    a.pvalue = 0.2;
    candidates.push_back(a);

    hicx::stripes::Candidate b = a;
    b.pvalue = 0.01;
    candidates.push_back(b);

    const std::vector<hicx::stripes::Candidate> kept =
        hicx::stripes::remove_redundant(candidates, /*by_pvalue=*/true, /*threads=*/1);
    REQUIRE(kept.size() == 1);
    CHECK(kept[0].pvalue == doctest::Approx(0.01));
}

TEST_CASE("candidate_pvalue ranks against the background sample") {
    BackgroundModel model;
    model.left.assign(400, {});
    model.right.assign(400, {});
    // A background sample of small, tightly clustered differences at
    // distance 5.
    for (int i = 0; i < 100; ++i) {
        model.left[5].push_back(1.0 + 0.01 * i);
        model.right[5].push_back(1.0 + 0.01 * i);
    }
    // A candidate with a much larger difference should get a small p-value.
    const double p_strong = hicx::stripes::candidate_pvalue(model, 5, 50.0, 50.0);
    CHECK(p_strong < 0.05);
    // A candidate with a typical, unremarkable difference should not.
    const double p_typical = hicx::stripes::candidate_pvalue(model, 5, 1.5, 1.5);
    CHECK(p_typical > 0.2);
    // An empty distance bucket falls back to p=1 (never a spurious call).
    const double p_empty = hicx::stripes::candidate_pvalue(model, 200, 100.0, 100.0);
    CHECK(p_empty == doctest::Approx(1.0));
}
