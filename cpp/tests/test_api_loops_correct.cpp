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
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "hicx/cnb.hpp"
#include "hicx/iterative_correction.hpp"
#include "hicx/math/ice.hpp"
#include "hicx/math/sparse_kernels.hpp"
#include "hicx/sparse_matrix.hpp"
#include "hicx/stats_ops.hpp"

namespace {

// A random symmetric matrix as the tools hold it (upper triangle CSR) and as
// the API receives it (the full COO of its csr, row major).
struct SymmetricPair {
    hicx::CsrMatrix upper;
    hicx::iterative_correction::Coo full;
};

SymmetricPair make_pair(std::int64_t n, double density, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> value(0.5, 20.0);
    std::uniform_real_distribution<double> chance(0.0, 1.0);
    std::vector<double> dense(static_cast<std::size_t>(n * n), 0.0);
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> data;
    for (std::int64_t row = 0; row < n; ++row) {
        for (std::int64_t column = row; column < n; ++column) {
            if (row != column && chance(rng) > density) {
                continue;
            }
            const double entry = value(rng);
            indices.push_back(static_cast<std::int32_t>(column));
            data.push_back(entry);
            dense[static_cast<std::size_t>(row * n + column)] = entry;
            dense[static_cast<std::size_t>(column * n + row)] = entry;
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(data.size());
    }
    SymmetricPair pair;
    pair.upper = hicx::CsrMatrix(n, n, std::move(indptr), std::move(indices), std::move(data), "float64");
    pair.upper.set_symmetry(hicx::Symmetry::UpperTriangle);
    pair.full.rows = n;
    pair.full.cols = n;
    for (std::int64_t row = 0; row < n; ++row) {
        for (std::int64_t column = 0; column < n; ++column) {
            const double v = dense[static_cast<std::size_t>(row * n + column)];
            if (v != 0.0) {
                pair.full.row.push_back(row);
                pair.full.col.push_back(column);
                pair.full.data.push_back(v);
            }
        }
    }
    return pair;
}

}  // namespace

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

// Not bit for bit: the Python scales the two mirrored entries of a pair
// separately, v * s[i] * s[j] and v * s[j] * s[i], which round differently,
// while the tool's triangle holds one value per pair. The API follows the
// Python; the two agree to rounding.
TEST_CASE("iterative_correction: agrees with the tool's ICE on upper triangle storage") {
    SymmetricPair pair = make_pair(120, 0.25, 5);
    hicx::ice::Options options;
    options.max_iterations = 300;
    options.threads = 1;
    hicx::kernels::DiagonalBlock block(pair.upper);
    const hicx::ice::Result tool = hicx::ice::correct(block, options);
    const hicx::iterative_correction::Result api =
        hicx::iterative_correction::correct(pair.full, 300, 1e-5);
    REQUIRE(api.outcome == hicx::iterative_correction::Outcome::Finished);
    REQUIRE(tool.correction_factors.size() == api.total_bias.size());
    CHECK(tool.converged);
    CHECK(api.converged_at == tool.iterations);
    for (std::size_t i = 0; i < api.total_bias.size(); ++i) {
        CHECK(api.total_bias[i] == doctest::Approx(tool.correction_factors[i]).epsilon(1e-12));
    }
    // the corrected values: every stored upper entry against the full COO
    for (std::size_t k = 0; k < pair.full.data.size(); ++k) {
        const std::int64_t r = pair.full.row[k];
        const std::int64_t c = pair.full.col[k];
        CHECK(pair.full.data[k] == doctest::Approx(pair.upper.at(std::min(r, c), std::max(r, c))).epsilon(1e-12));
    }
}

TEST_CASE("iterative_correction: symmetry measure and guards") {
    hicx::iterative_correction::Coo m;
    m.rows = m.cols = 2;
    m.row = {0, 0, 1, 1};
    m.col = {0, 1, 0, 1};
    m.data = {1.0, 2.0, 2.0, 3.0};
    CHECK(hicx::iterative_correction::asymmetry(m) == 0.0);
    m.data[1] = 4.0;  // |A - A.T| has two entries of 2, mean 1, matrix mean 2.5
    CHECK(hicx::iterative_correction::asymmetry(m) == doctest::Approx(0.4));
    hicx::iterative_correction::Coo tiny;
    tiny.rows = tiny.cols = 2;
    tiny.row = {0, 1};
    tiny.col = {0, 1};
    tiny.data = {1e-200, 1.0};
    CHECK(hicx::iterative_correction::correct(tiny, 50, 1e-5).outcome ==
          hicx::iterative_correction::Outcome::LargeValuesInLoop);
}
