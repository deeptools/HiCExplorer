// Unit tests of the core code behind the Python API of hicAverageRegions,
// hicCompartmentalization, hicPCA, hicPlotDistVsCounts, hicPlotSVL,
// hicPlotMatrix and hicAggregateContacts (cpp/python/API_BINDINGS.md).
// The equivalence with 3.7.6 is tested by hicexplorer/test/api/cases_*.py;
// these cases pin the arithmetic and the clamping rules.

#include <doctest/doctest.h>

#include <cstdint>

#include "hicx/average_regions.hpp"

namespace ar = hicx::average_regions;

TEST_CASE("average_regions: centre truncates towards zero") {
    CHECK(ar::center_of(10.0, 21.0) == 15);
    CHECK(ar::center_of(2500001.5, 3500000.0) == 3000000);
    CHECK(ar::center_of(-3.0, 0.0) == -1);
}

TEST_CASE("average_regions: genomic window clamps to 0 and to the chromosome end minus one") {
    ar::Window w = ar::viewpoint_range(100, 30, 50, 1000);
    CHECK(w.start == 70);
    CHECK(w.end == 150);
    CHECK_FALSE(w.start_out_of_range);
    CHECK_FALSE(w.end_out_of_range);
    w = ar::viewpoint_range(10, 30, 5000, 1000);
    CHECK(w.start == 0);
    CHECK(w.end == 999);
    CHECK(w.start_out_of_range);
    CHECK(w.end_out_of_range);
    // reaching exactly the end is not out of range
    w = ar::viewpoint_range(500, 0, 500, 1000);
    CHECK(w.end == 1000);
    CHECK_FALSE(w.end_out_of_range);
}

TEST_CASE("average_regions: bin window clamps to the chromosome's bins, last bin inclusive") {
    ar::Window w = ar::viewpoint_range_bins(105, 10, 20, 100, 120);
    CHECK(w.start == 100);
    CHECK(w.end == 120);
    CHECK(w.start_out_of_range);
    CHECK(w.end_out_of_range);
    w = ar::viewpoint_range_bins(110, 10, 10, 100, 120);
    CHECK(w.start == 100);
    CHECK(w.end == 120);
    CHECK_FALSE(w.start_out_of_range);
    CHECK_FALSE(w.end_out_of_range);
}
