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
#include "hicx/correct_api.hpp"
#include "hicx/correct_ops.hpp"
#include "hicx/iterative_correction.hpp"
#include "hicx/math/ice.hpp"
#include "hicx/math/sparse_kernels.hpp"
#include "hicx/merge_loops.hpp"
#include "hicx/reduce_matrix.hpp"
#include "hicx/reduce_matrix_entries.hpp"
#include "hicx/scipy_coo.hpp"
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

TEST_CASE("reduce_entries: agrees with the tools' reduce_matrix on a CSR") {
    SymmetricPair pair = make_pair(90, 0.3, 9);
    std::vector<std::int32_t> row(pair.full.row.begin(), pair.full.row.end());
    std::vector<std::int32_t> col(pair.full.col.begin(), pair.full.col.end());
    const hicx::CsrMatrix full = hicx::CsrMatrix::from_coo(90, 90, row, col, pair.full.data, "float64");
    std::vector<std::vector<std::int64_t>> groups;
    std::vector<std::int64_t> group_of(90, -1);
    for (std::int64_t i = 0; i < 90; i += 3) {
        if (i == 30) {
            continue;  // a dropped group
        }
        groups.push_back({i, i + 1, i + 2});
        for (std::int64_t j = i; j < i + 3; ++j) {
            group_of[static_cast<std::size_t>(j)] = static_cast<std::int64_t>(groups.size()) - 1;
        }
    }
    hicx::reduce_entries::Arithmetic a;
    a.cast = [](double v) { return v; };
    a.add = [](double x, double y) { return x + y; };
    a.subtract = [](double x, double y) { return x - y; };
    a.twice = [](double x) { return x * 2.0; };
    hicx::reduce_entries::Entries entries;
    entries.rows = entries.cols = 90;
    entries.row = pair.full.row;
    entries.col = pair.full.col;
    entries.data = pair.full.data;
    for (const bool diagonal : {true, false}) {
        const hicx::CsrMatrix tool = hicx::reduce_matrix(full, groups, true, diagonal);
        const hicx::reduce_entries::Cells api = hicx::reduce_entries::reduce(
            hicx::reduce_entries::upper_triangle(entries), group_of, static_cast<std::int64_t>(groups.size()), true,
            diagonal, a);
        std::size_t k = 0;
        for (std::int64_t r = 0; r < tool.rows(); ++r) {
            for (std::int64_t c = 0; c < tool.cols(); ++c) {
                const double expected = tool.at(r, c);
                if (expected == 0.0) {
                    continue;
                }
                REQUIRE(k < api.data.size());
                CHECK(api.row[k] == r);
                CHECK(api.col[k] == c);
                CHECK(api.data[k] == expected);
                ++k;
            }
        }
        CHECK(k == api.data.size());
    }
}

TEST_CASE("scipy_coo: running_window agrees with the tools' running_window on counts") {
    // integer counts: every order of summation gives the same values
    SymmetricPair pair = make_pair(60, 0.3, 21);
    std::vector<std::int32_t> row(pair.full.row.begin(), pair.full.row.end());
    std::vector<std::int32_t> col(pair.full.col.begin(), pair.full.col.end());
    std::vector<double> counts;
    for (const double v : pair.full.data) {
        counts.push_back(std::floor(v));
    }
    const hicx::CsrMatrix full = hicx::CsrMatrix::from_coo(60, 60, row, col, counts, "int64");
    hicx::scipy_coo::Csr<std::int64_t> csr;
    csr.rows = csr.cols = 60;
    csr.indptr = full.indptr();
    csr.indices.assign(full.indices().begin(), full.indices().end());
    for (const double v : full.data()) {
        csr.data.push_back(static_cast<std::int64_t>(v));
    }
    for (const std::int64_t bins : {3, 5}) {
        const hicx::CsrMatrix tool = hicx::running_window(full, bins);
        const hicx::scipy_coo::Csr<std::int64_t> api = hicx::scipy_coo::running_window(csr, bins);
        for (std::int64_t r = 0; r < 60; ++r) {
            for (std::int64_t k = api.indptr[static_cast<std::size_t>(r)];
                 k < api.indptr[static_cast<std::size_t>(r) + 1]; ++k) {
                CHECK(static_cast<double>(api.data[static_cast<std::size_t>(k)]) ==
                      tool.at(r, api.indices[static_cast<std::size_t>(k)]));
            }
        }
        CHECK(static_cast<std::size_t>(api.indptr.back()) == tool.nnz());
    }
}

TEST_CASE("scipy_coo: tocsr sums duplicates and keeps zeros") {
    hicx::scipy_coo::Coo<double> coo;
    coo.rows = coo.cols = 2;
    coo.row = {1, 0, 1, 0};
    coo.col = {0, 1, 0, 0};
    coo.data = {1.5, 0.0, 2.5, 3.0};
    const hicx::scipy_coo::Csr<double> csr = hicx::scipy_coo::tocsr(coo);
    CHECK(csr.indptr == std::vector<std::int64_t>{0, 2, 3});
    CHECK(csr.indices == std::vector<std::int64_t>{0, 1, 0});
    CHECK(csr.data == std::vector<double>{3.0, 0.0, 4.0});
}

TEST_CASE("merge_loops: the constructed frame of cases_hicMergeLoops.py") {
    // ROWS of hicexplorer/test/api/cases_hicMergeLoops.py at resolution
    // 10000; 3.7.6 returns [0, 3]
    const std::vector<std::int64_t> chrom = {0, 0, 0, 1, 1};
    const std::vector<std::int64_t> xb = {10000, 10000, 12000, 50000, 52000};
    const std::vector<std::int64_t> xe = {20000, 15000, 17000, 60000, 57000};
    const std::vector<std::int64_t> yb = {250000, 250000, 252000, 90000, 92000};
    const std::vector<std::int64_t> ye = {260000, 255000, 257000, 100000, 97000};
    hicx::merge_loops::Intervals x;
    hicx::merge_loops::Intervals y;
    x.tree = y.tree = chrom;
    x.begin = xb;
    x.end = xe;
    y.begin = yb;
    y.end = ye;
    x.tree_count = y.tree_count = 2;
    x.dictionary = y.dictionary = {0, 1};
    std::vector<hicx::merge_loops::Loop> loops;
    for (std::size_t r = 0; r < chrom.size(); ++r) {
        hicx::merge_loops::Loop loop;
        loop.x_tree = loop.y_tree = chrom[r];
        loop.y_guard = true;
        const std::int64_t fx = 10000 - (xe[r] - xb[r]);
        const std::int64_t fy = 10000 - (ye[r] - yb[r]);
        loop.x_begin = xb[r] - fx - 1;
        loop.x_end = xe[r] + fx + 1;
        loop.y_begin = yb[r] - fy - 1;
        loop.y_end = ye[r] + fy + 1;
        loops.push_back(loop);
    }
    const hicx::merge_loops::Result result = hicx::merge_loops::merge(x, y, loops);
    CHECK(result.error == hicx::merge_loops::Error::None);
    CHECK(result.ids == std::vector<std::int64_t>{0, 3});
    loops[0].x_tree = -1;
    CHECK(hicx::merge_loops::merge(x, y, loops).error == hicx::merge_loops::Error::UnboundX);
}

TEST_CASE("correct_api: MAD agrees with the tool's MAD") {
    std::mt19937 rng(17);
    std::gamma_distribution<double> gamma(2.0, 30.0);
    std::vector<double> points;
    for (int i = 0; i < 501; ++i) {
        points.push_back(i % 37 == 0 ? 0.0 : std::floor(gamma(rng)));
    }
    const hicx::correct_api::MadStats api = hicx::correct_api::mad(points, 501, 1);
    const hicx::correct::Mad tool(points);
    CHECK(api.median == tool.median());
    CHECK(api.med_abs_deviation == tool.median_absolute_deviation());
    REQUIRE(api.modified_z_score.size() == tool.modified_z_scores().size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        CHECK(api.modified_z_score[i] == doctest::Approx(tool.modified_z_scores()[i]).epsilon(1e-15));
    }
    CHECK(std::isnan(hicx::correct_api::median({1.0, std::nan(""), 3.0})));
}

TEST_CASE("correct_api: fill_gaps fills an isolated bin with the mean of its neighbours") {
    // a 5 by 5 matrix of ones but for an empty bin 2
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    std::vector<double> data;
    for (std::int64_t i = 0; i < 5; ++i) {
        for (std::int64_t j = 0; j < 5; ++j) {
            if (i != 2 && j != 2) {
                row.push_back(i);
                col.push_back(j);
                data.push_back(1.0);
            }
        }
    }
    const auto filled = hicx::scipy_coo::Csr<double>(hicx::correct_api::fill_gaps<double>(5, 5, row, col, data, {2}));
    // row 2, columns 1..3: (M[1, j - 1] + M[3, j + 1]) / 2; the four neighbours of M[2, 2] get the mean of
    // the diagonal neighbours, 1
    std::vector<double> dense(25, 0.0);
    for (std::int64_t i = 0; i < 5; ++i) {
        for (std::int64_t k = filled.indptr[static_cast<std::size_t>(i)]; k < filled.indptr[static_cast<std::size_t>(i) + 1]; ++k) {
            dense[static_cast<std::size_t>(i * 5 + filled.indices[static_cast<std::size_t>(k)])] =
                filled.data[static_cast<std::size_t>(k)];
        }
    }
    CHECK(dense[2 * 5 + 1] == 1.0);
    CHECK(dense[2 * 5 + 2] == 1.0);
    CHECK(dense[1 * 5 + 2] == 1.0);
    CHECK(dense[2 * 5 + 0] == 0.0);
    CHECK_THROWS_AS(hicx::correct_api::fill_gaps<double>(5, 5, row, col, data, {4}), hicx::PythonError);
}

TEST_CASE("scipy_coo: assign_zero_block inserts explicit zeros as scipy's _insert_many") {
    hicx::scipy_coo::Csr<double> m;
    m.rows = m.cols = 3;
    m.indptr = {0, 2, 3, 4};
    m.indices = {2, 0, 1, 2};  // row 0 unsorted
    m.data = {5.0, 1.0, 2.0, 3.0};
    const auto unsorted = hicx::scipy_coo::assign_zero_block(m, false, 0, 2, 1, 3);
    CHECK(unsorted.inserted);
    // row 0: the inserted column 1 first, then the old entries (2 zeroed, 0 kept)
    CHECK(unsorted.matrix.indices == std::vector<std::int64_t>{1, 2, 0, 2, 1, 2});
    CHECK(unsorted.matrix.data == std::vector<double>{0.0, 0.0, 1.0, 0.0, 0.0, 3.0});
    const auto sorted = hicx::scipy_coo::assign_zero_block(m, true, 0, 2, 1, 3);
    CHECK(sorted.matrix.indices == std::vector<std::int64_t>{0, 1, 2, 1, 2, 2});
    const auto overwrite = hicx::scipy_coo::assign_zero_block(m, true, 2, 3, 2, 3);
    CHECK_FALSE(overwrite.inserted);
    CHECK(overwrite.matrix.data == std::vector<double>{5.0, 1.0, 2.0, 0.0});
}
