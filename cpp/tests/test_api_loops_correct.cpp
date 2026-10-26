// The core code behind the public Python API of hicexplorer.lib.cnb,
// hicexplorer.iterativeCorrection, hicexplorer.reduceMatrix,
// hicexplorer.hicMergeMatrixBins, hicexplorer.hicMergeLoops,
// hicexplorer.hicCorrectMatrix, hicexplorer.hicAdjustMatrix and
// hicexplorer.hicDetectLoops. Their equivalence with HiCExplorer 3.7.6 is
// tested function by function from Python (hicexplorer/test/api); these
// tests pin the C++ contracts and check the API kernels against the tool
// kernels that implement the same Python functions.

#include <doctest/doctest.h>

#include <cmath>
#include <limits>

#include "hicx/cnb.hpp"
#include "hicx/stats_ops.hpp"

TEST_CASE("cnb: scipy 1.14 values") {
    // scipy.special.betainc(2.5, 4.0, 0.3), gammaln(7.25), xlog1py(3, -0.3)
    CHECK(hicx::cnb::cdf(3.0, 2.5, 0.3) == 0.3521975859067672);
    CHECK(hicx::cnb::gammaln(7.25) == 7.0521854507385395);
    CHECK(hicx::cnb::xlog1py(3.0, -0.3) == -1.070024831816197);
    CHECK(hicx::cnb::xlog1py(0.0, -1.0) == 0.0);
    CHECK(std::isnan(hicx::cnb::xlog1py(0.0, std::numeric_limits<double>::quiet_NaN())));
    // the value where cephes and Boost differ (hicx/scipy_special.hpp)
    CHECK(hicx::cnb::cdf(9.8, 58.254577634279, 0.998604709282) == 1.0);
}

TEST_CASE("cnb: cdf agrees with the tool's survival function") {
    for (const double x : {0.0, 1.0, 3.5, 12.0}) {
        const double sf = hicx::stats::nbinom_sf(x, 2.5, 0.4);
        CHECK(std::fabs((1.0 - hicx::cnb::cdf(x, 2.5, 0.4)) - sf) <= 1e-14);
    }
}

TEST_CASE("cnb: pdf is the derivative of cdf in x for integer x") {
    // for integer x the continuous pdf is the negative binomial pmf, and
    // cdf(x) - cdf(x - 1) is the same probability
    for (int x = 1; x < 10; ++x) {
        const double pmf = hicx::cnb::pdf(x, 3.0, 0.4);
        const double difference = hicx::cnb::cdf(x, 3.0, 0.4) - hicx::cnb::cdf(x - 1, 3.0, 0.4);
        CHECK(std::fabs(pmf - difference) <= 1e-12);
    }
}
