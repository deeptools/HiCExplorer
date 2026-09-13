// Unit tests for hicDetectLoops: the statistics added to hicx::stats and the
// computational core in tools/detect_loops_impl.hpp.
//
// Every expected value here was produced by the reference environment, either
// by scipy directly or by running the hicexplorer helper the code replaces on
// the same fixture, and is quoted to seventeen significant digits. Nothing is
// asserted against a number this code computed itself.

#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "../tools/detect_loops_impl.hpp"
#include "hicx/adjust_ops.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/numpy_compat.hpp"
#include "hicx/simd_reduce.hpp"
#include "hicx/stats_ops.hpp"

namespace {

// The 6 by 6 fixture the hicexplorer helpers were run on to produce the
// expected values below: the strict upper triangle, no diagonal, which is the
// shape hicDetectLoops hands to obs_exp_matrix.
//
//   (0,1)=5 (0,2)=3 (0,4)=10 (1,2)=7 (1,3)=2 (2,5)=4 (3,4)=6 (4,5)=1
hicx::CsrMatrix fixture(const std::string& dtype) {
    const std::vector<std::int32_t> rows{0, 0, 0, 1, 1, 2, 3, 4};
    const std::vector<std::int32_t> cols{1, 2, 4, 2, 3, 5, 4, 5};
    std::vector<double> values{5, 3, 10, 7, 2, 4, 6, 1};
    hicx::CsrMatrix matrix =
        hicx::CsrMatrix::from_coo(6, 6, rows, cols, std::move(values), dtype);
    matrix.set_symmetry(hicx::Symmetry::Full);
    return matrix;
}

hicx::CsrMatrix from_dense(const std::vector<std::vector<double>>& dense) {
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    for (std::size_t r = 0; r < dense.size(); ++r) {
        for (std::size_t c = 0; c < dense[r].size(); ++c) {
            if (dense[r][c] != 0.0) {
                rows.push_back(static_cast<std::int32_t>(r));
                cols.push_back(static_cast<std::int32_t>(c));
                values.push_back(dense[r][c]);
            }
        }
    }
    hicx::CsrMatrix matrix = hicx::CsrMatrix::from_coo(
        static_cast<std::int64_t>(dense.size()),
        static_cast<std::int64_t>(dense[0].size()), rows, cols, std::move(values),
        "float64");
    matrix.set_symmetry(hicx::Symmetry::Full);
    return matrix;
}

const char* const kDataDir = HICX_TEST_DATA_DIR;

}  // namespace

TEST_CASE("betainc reproduces scipy.special.betainc") {
    // scipy 1.x in the reference environment, twelve to seventeen digits.
    CHECK(hicx::stats::betainc(0.5, 0.5, 0.29999999999999999) ==
          doctest::Approx(0.36901011956554536).epsilon(1e-13));
    CHECK(hicx::stats::betainc(1, 2, 0.25) ==
          doctest::Approx(0.4375).epsilon(1e-13));
    CHECK(hicx::stats::betainc(2.5, 3.5, 0.5) ==
          doctest::Approx(0.66976527263135488).epsilon(1e-13));
    CHECK(hicx::stats::betainc(10, 1, 0.90000000000000002) ==
          doctest::Approx(0.34867844010000004).epsilon(1e-13));
    CHECK(hicx::stats::betainc(0.10000000000000001, 0.20000000000000001,
                               0.050000000000000003) ==
          doctest::Approx(0.50953912153463998).epsilon(1e-13));
    CHECK(hicx::stats::betainc(200, 3, 0.995) ==
          doctest::Approx(0.91831194886973366).epsilon(1e-13));
    CHECK(hicx::stats::betainc(0.001, 0.001, 0.5) ==
          doctest::Approx(0.49999999999999994).epsilon(1e-13));
    CHECK(hicx::stats::betainc(5, 5, 0.999) ==
          doctest::Approx(0.99999999999987443).epsilon(1e-13));
    // The parameters of a real degenerate fit from small_test_matrix.h5.
    CHECK(hicx::stats::betainc(0.15883824665748106, 1, 0.00033612374545108738) ==
          doctest::Approx(0.2807212068279708).epsilon(1e-13));
    // The logarithmic branch, where the gamma functions would overflow.
    CHECK(hicx::stats::betainc(100, 0.5, 0.20000000000000001) ==
          doctest::Approx(7.9762895539848858e-72).epsilon(1e-12));
    CHECK(hicx::stats::betainc(3, 700, 0.01) ==
          doctest::Approx(0.97135907570396118).epsilon(1e-13));
    CHECK(hicx::stats::betainc(700, 3, 0.98999999999999999) ==
          doctest::Approx(0.028640924296038869).epsilon(1e-13));
}

TEST_CASE("betainc edge cases follow cephes") {
    CHECK(hicx::stats::betainc(2.0, 3.0, 0.0) == 0.0);
    CHECK(hicx::stats::betainc(2.0, 3.0, 1.0) == 1.0);
    CHECK(std::isnan(hicx::stats::betainc(0.0, 3.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, -1.0, 0.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, 3.0, 1.5)));
    CHECK(std::isnan(hicx::stats::betainc(2.0, 3.0, -0.5)));
}

TEST_CASE("nbinom_sf reproduces 1 - hicexplorer.lib.cnb.cdf") {
    CHECK(hicx::stats::nbinom_sf(2, 10, 0.90000000000000002) ==
          doctest::Approx(0.11086997774499996).epsilon(1e-12));
    CHECK(hicx::stats::nbinom_sf(4411, 0.15883824665748106, 0.00033612374545108738) ==
          doctest::Approx(0.019744103473751995).epsilon(1e-12));
    CHECK(hicx::stats::nbinom_sf(7, 207.07414628962866, 0.99466269511004068) ==
          doctest::Approx(2.3835777908431943e-05).epsilon(1e-11));
    CHECK(hicx::stats::nbinom_sf(1.5, 3, 0.5) ==
          doctest::Approx(0.59109707489812935).epsilon(1e-12));
    CHECK(hicx::stats::nbinom_sf(0, 2, 0.40000000000000002) ==
          doctest::Approx(0.83999999999999997).epsilon(1e-12));
}

TEST_CASE("fit_nbinom returns its starting point when factorial overflows") {
    // scipy.special.factorial is infinite above 170, so the objective
    // np.sum(np.log(factorial(X))) is infinite for every parameter pair, the
    // forward difference gradient is inf - inf = NaN and fmin_l_bfgs_b stops
    // at iteration zero. The fit is then exactly R's fitdistr moment
    // estimator, and it is bit reproducible. This is a real distance
    // distribution from small_test_matrix.h5 chromosome chr2L.
    const std::vector<double> data{24, 7, 15, 17, 31, 172, 4411, 2221,
                                   7,  15, 19, 19, 31, 28,  69};
    const hicx::stats::NBinomFit fit = hicx::stats::fit_nbinom(data);
    // fit_nbinom.fit on the same array in the reference environment.
    CHECK(fit.size == doctest::Approx(0.15883824665748106).epsilon(1e-15));
    CHECK(fit.prob == doctest::Approx(0.00033612374545108738).epsilon(1e-15));
    CHECK(fit.iterations == 0);

    // Without the overflow the optimiser really runs.
    const std::vector<double> small{2.5, 3.1, 1.7, 4.0, 2.2, 9.9, 1.1, 0.5};
    const hicx::stats::NBinomFit fit_small = hicx::stats::fit_nbinom(small);
    CHECK(fit_small.iterations > 0);
    CHECK(fit_small.size > 2.0);
}

TEST_CASE("fit_nbinom in float32 reproduces the reference's failed line search") {
    // A float32 obs/exp matrix makes two of the five terms of the log
    // likelihood single precision, which drowns the 1e-8 forward difference
    // step in rounding noise. scipy then aborts its line search after three
    // iterations a few parts in 1e-6 from the start, at size 10 rather than at
    // the maximum likelihood 207. The distribution is distance 1 of
    // GSE63525_GM12878_insitu_primary_2_5mb.cool chromosome 1.
    std::vector<double> data(90, 1.0);
    for (std::size_t i = 0; i < 10; ++i) {
        data[i] = 2.0;
    }
    data[42] = 2.0020203590393066;
    const hicx::stats::NBinomFit narrow =
        hicx::stats::fit_nbinom(data, hicx::stats::NBinomPrecision::Float32);
    const hicx::stats::NBinomFit wide =
        hicx::stats::fit_nbinom(data, hicx::stats::NBinomPrecision::Float64);
    // The float32 objective leaves the fit at its start; the float64 one
    // converges to a size an order of magnitude larger. The point of the test
    // is that the two are different and that the float32 one is the one that
    // stays put, which is what the reference does.
    CHECK(narrow.size < 11.0);
    CHECK(wide.size > 50.0);
    CHECK(hicx::stats::fit_nbinom(data).size == wide.size);
}

TEST_CASE("in float32 the factorial term overflows at 34.6, not at 170.6") {
    // scipy.special.factorial returns a float64 array whatever it is given,
    // but for a float32 argument it evaluates gamma(x + 1) in the float32 loop
    // and widens the result, so it overflows at the float32 limit. The
    // objective is then infinite, the fit is degenerate and fmin_l_bfgs_b
    // returns its starting point; the same data in float64 converges instead.
    // On gm12878_chr1.cool 168 of the 199 genomic distances are in this state,
    // and evaluating the term in float64 there calls a loop the reference does
    // not.
    const std::vector<double> data{0.5, 1.5, 40.0, 2.25, 0.125, 3.75, 60.0, 0.0625};
    const hicx::stats::NBinomFit narrow =
        hicx::stats::fit_nbinom(data, hicx::stats::NBinomPrecision::Float32);
    const hicx::stats::NBinomFit wide =
        hicx::stats::fit_nbinom(data, hicx::stats::NBinomPrecision::Float64);
    // fit_nbinom.fit on the same values as a float32 array. The fit is its own
    // starting point, so this also pins that the starting point is right,
    // including the fact that `size = (m ** 2) / (v - m)` promotes a float32
    // mean to float64 through the integer exponent while `v - m` stays
    // float32. Computing it wholly in float32 gives 0.40079742670059204.
    CHECK(narrow.size == doctest::Approx(0.40079743194169981).epsilon(1e-15));
    CHECK(narrow.prob == doctest::Approx(0.028784161851670913).epsilon(1e-15));
    CHECK(narrow.iterations == 0);
    CHECK(narrow.size != doctest::Approx(0.40079742670059204).epsilon(1e-12));
    // The same values in float64 are inside the float64 factorial range, so
    // the optimiser runs and lands somewhere else entirely. The tolerance is
    // 1e-4 rather than the 1e-15 above because this is the only branch where
    // the optimiser actually runs, and the projected L-BFGS of core/lbfgsb.cpp
    // stops at a different point inside the same tolerance from the Fortran
    // L-BFGS-B scipy calls; that difference is measured in cpp/STATUS.md, and
    // it is precisely what makes the degenerate cases above worth reproducing.
    CHECK(wide.size == doctest::Approx(0.32561896056641176).epsilon(1e-4));
    CHECK(wide.iterations > 0);
    CHECK(wide.size < narrow.size * 0.9);
}

TEST_CASE("numpy_slice keeps Python's negative index semantics") {
    using hicx::loops::numpy_slice;
    // a[:3] of a length 10 axis
    CHECK(numpy_slice(10, 0, 3).begin == 0);
    CHECK(numpy_slice(10, 0, 3).end == 3);
    // a[:-2] is everything but the last two, not the empty slice. This is the
    // case candidate_region_test hits when the peak sits within --peakWidth of
    // the edge of its window.
    CHECK(numpy_slice(10, 0, -2).begin == 0);
    CHECK(numpy_slice(10, 0, -2).end == 8);
    // a[-3:] is the last three.
    CHECK(numpy_slice(10, -3, 10).begin == 7);
    CHECK(numpy_slice(10, -3, 10).end == 10);
    // Out of range bounds clamp rather than throw.
    CHECK(numpy_slice(4, -100, 100).begin == 0);
    CHECK(numpy_slice(4, -100, 100).end == 4);
    // An inverted range is empty, never negative.
    CHECK(numpy_slice(10, 7, 3).size() == 0);
}

TEST_CASE("dense_block and flatten_slice match numpy basic slicing") {
    const hicx::CsrMatrix matrix = from_dense({{1, 2, 0, 4},
                                               {0, 5, 6, 0},
                                               {7, 0, 8, 9},
                                               {0, 1, 0, 2}});
    const hicx::loops::DenseBlock block =
        hicx::loops::dense_block(matrix, 0, 3, 1, 4);
    CHECK(block.rows == 3);
    CHECK(block.cols == 3);
    CHECK(block.at(0, 0) == 2);
    CHECK(block.at(0, 2) == 4);
    CHECK(block.at(2, 1) == 8);
    CHECK(block.at(1, 0) == 5);

    // block[:-1, :] is the first two rows.
    const std::vector<double> top =
        hicx::loops::flatten_slice(block, 0, -1, 0, block.cols);
    CHECK(top == std::vector<double>{2, 0, 4, 5, 6, 0});
    // block[2:, 1:] is one row of two.
    const std::vector<double> corner =
        hicx::loops::flatten_slice(block, 2, block.rows, 1, block.cols);
    CHECK(corner == std::vector<double>{8, 9});
    CHECK(hicx::loops::element(matrix, 2, 3) == 9);
    CHECK(hicx::loops::element(matrix, 0, 2) == 0);
}

TEST_CASE("group_by_distance buckets the stored entries in CSR order") {
    const hicx::CsrMatrix matrix = fixture("float64");
    const hicx::loops::DistanceGroups groups =
        hicx::loops::group_by_distance(matrix);
    CHECK(groups.min_distance == 1);
    CHECK(groups.max_distance == 4);
    // distance 1: (0,1), (1,2), (3,4), (4,5)
    CHECK(groups.bucket(0).size() == 4);
    // distance 2: (0,2), (1,3)
    CHECK(groups.bucket(1).size() == 2);
    // distance 3: (2,5)
    CHECK(groups.bucket(2).size() == 1);
    // distance 4: (0,4)
    CHECK(groups.bucket(3).size() == 1);
    const std::vector<double>& values = matrix.data();
    CHECK(values[static_cast<std::size_t>(groups.bucket(0)[0])] == 5);
    CHECK(values[static_cast<std::size_t>(groups.bucket(0)[1])] == 7);
    CHECK(values[static_cast<std::size_t>(groups.bucket(0)[2])] == 6);
    CHECK(values[static_cast<std::size_t>(groups.bucket(0)[3])] == 1);
}

TEST_CASE("expected_interactions divides by n + 1 - d, not by n - d") {
    const hicx::CsrMatrix matrix = fixture("float64");
    const hicx::loops::DistanceGroups groups =
        hicx::loops::group_by_distance(matrix);
    for (const unsigned int threads : {1U, 2U, 8U}) {
        const std::vector<double> expected =
            hicx::loops::expected_interactions(matrix, groups, threads);
        REQUIRE(expected.size() == 6);
        // hicexplorer.utilities.expected_interactions on the same fixture.
        CHECK(expected[0] == 0.0);
        CHECK(expected[1] == doctest::Approx(3.1666666666666665).epsilon(1e-15));
        CHECK(expected[2] == doctest::Approx(1.0).epsilon(1e-15));
        CHECK(expected[3] == doctest::Approx(1.0).epsilon(1e-15));
        CHECK(expected[4] == doctest::Approx(3.3333333333333335).epsilon(1e-15));
        CHECK(expected[5] == 0.0);
        // 19 / 6, the diagonal length plus one, rather than 19 / 5.
        CHECK(expected[1] != doctest::Approx(19.0 / 5.0));
    }
}

TEST_CASE("expected_interactions_non_zero divides by the stored count") {
    const std::vector<double> expected =
        hicx::loops::expected_interactions_non_zero(fixture("float64"));
    REQUIRE(expected.size() == 6);
    CHECK(expected[0] == 0.0);
    CHECK(expected[1] == doctest::Approx(4.75).epsilon(1e-15));
    CHECK(expected[2] == doctest::Approx(2.5).epsilon(1e-15));
    CHECK(expected[3] == doctest::Approx(4.0).epsilon(1e-15));
    CHECK(expected[4] == doctest::Approx(10.0).epsilon(1e-15));
    CHECK(expected[5] == 0.0);
}

TEST_CASE("obs_exp_matrix looks the expected value up at ceil(d / 2)") {
    const hicx::CsrMatrix matrix = fixture("float64");
    const hicx::loops::DistanceGroups groups =
        hicx::loops::group_by_distance(matrix);
    const hicx::CsrMatrix obs_exp =
        hicx::loops::obs_exp_matrix(matrix, groups, 1);
    const std::vector<double>& v = obs_exp.data();
    REQUIRE(v.size() == 8);
    // hicexplorer.utilities.obs_exp_matrix on the same fixture. Note the
    // entries at distance 2 and 4 are divided by the expected value of
    // distance 1 and 2 respectively.
    CHECK(v[0] == doctest::Approx(1.5789473684210527).epsilon(1e-15));
    CHECK(v[1] == doctest::Approx(0.94736842105263164).epsilon(1e-15));
    CHECK(v[2] == doctest::Approx(10.0).epsilon(1e-15));
    CHECK(v[3] == doctest::Approx(2.2105263157894739).epsilon(1e-15));
    CHECK(v[4] == doctest::Approx(0.63157894736842113).epsilon(1e-15));
    CHECK(v[5] == doctest::Approx(4.0).epsilon(1e-15));
    CHECK(v[6] == doctest::Approx(1.8947368421052633).epsilon(1e-15));
    CHECK(v[7] == doctest::Approx(0.31578947368421056).epsilon(1e-15));
    // The distance 4 entry divided by expected[2] = 1, not by expected[4].
    CHECK(v[2] == 10.0);
}

TEST_CASE("obs_exp_matrix truncates to the input dtype") {
    const hicx::CsrMatrix matrix = fixture("int64");
    const hicx::loops::DistanceGroups groups =
        hicx::loops::group_by_distance(matrix);
    hicx::CsrMatrix obs_exp = hicx::loops::obs_exp_matrix(matrix, groups, 1);
    CHECK(obs_exp.data() == std::vector<double>{1, 0, 10, 2, 0, 4, 1, 0});
    // Three of the eight ratios truncate to zero, so eliminate_zeros leaves a
    // different sparsity pattern from the raw matrix and hicDetectLoops.py:888
    // abandons the chromosome.
    obs_exp.eliminate_zeros();
    CHECK(obs_exp.stored_nnz() == 5);
    CHECK(obs_exp.stored_nnz() != matrix.stored_nnz());
}

TEST_CASE("obs_exp_matrix_non_zero stays in float32 and uses a 1e-9 epsilon") {
    const hicx::CsrMatrix plain =
        hicx::loops::obs_exp_matrix_non_zero(fixture("float64"), false);
    const std::vector<double>& v = plain.data();
    REQUIRE(v.size() == 8);
    CHECK(plain.dtype() == "float32");
    CHECK(v[0] == doctest::Approx(1.0526316165924072).epsilon(1e-15));
    CHECK(v[1] == doctest::Approx(1.2000000476837158).epsilon(1e-15));
    CHECK(v[2] == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(v[3] == doctest::Approx(1.4736841917037964).epsilon(1e-15));
    CHECK(v[4] == doctest::Approx(0.80000001192092896).epsilon(1e-15));
    CHECK(v[5] == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(v[6] == doctest::Approx(1.263157844543457).epsilon(1e-15));
    CHECK(v[7] == doctest::Approx(0.21052631735801697).epsilon(1e-15));
    // An integer input goes through the same float32 path, so the two agree.
    const hicx::CsrMatrix from_int =
        hicx::loops::obs_exp_matrix_non_zero(fixture("int64"), false);
    CHECK(from_int.data() == v);

    const hicx::CsrMatrix ligation =
        hicx::loops::obs_exp_matrix_non_zero(fixture("float64"), true);
    const std::vector<double>& w = ligation.data();
    REQUIRE(w.size() == 8);
    CHECK(w[0] == doctest::Approx(0.24691358208656311).epsilon(1e-15));
    CHECK(w[1] == doctest::Approx(0.63333332538604736).epsilon(1e-15));
    CHECK(w[2] == doctest::Approx(2.1111111640930176).epsilon(1e-15));
    CHECK(w[3] == doctest::Approx(1.5555555820465088).epsilon(1e-15));
    CHECK(w[4] == doctest::Approx(0.56296294927597046).epsilon(1e-15));
    // Bin 5 has an empty row, so its ligation factor is zero and the quotient
    // is infinite; pToEpsilon replaces it with 1e-9, in float32.
    CHECK(w[5] == doctest::Approx(9.9999997171806854e-10).epsilon(1e-15));
    CHECK(w[6] == doctest::Approx(8.0).epsilon(1e-15));
    CHECK(w[7] == doctest::Approx(9.9999997171806854e-10).epsilon(1e-15));
}

TEST_CASE("neighborhood_merge keeps only the maximum of its own window") {
    //  a 7 by 7 matrix with two peaks four bins apart
    std::vector<std::vector<double>> dense(7, std::vector<double>(7, 0.0));
    dense[1][2] = 3.0;
    dense[1][3] = 9.0;   // the local maximum
    dense[2][3] = 4.0;
    dense[5][6] = 7.0;   // isolated, so also a maximum
    const hicx::CsrMatrix matrix = from_dense(dense);
    const std::vector<hicx::loops::Candidate> candidates{
        {1, 2}, {1, 3}, {2, 3}, {5, 6}};
    for (const unsigned int threads : {1U, 4U}) {
        const std::vector<hicx::loops::Candidate> kept =
            hicx::loops::neighborhood_merge(candidates, 1, matrix, threads);
        REQUIRE(kept.size() == 2);
        CHECK(kept[0].row == 1);
        CHECK(kept[0].col == 3);
        CHECK(kept[1].row == 5);
        CHECK(kept[1].col == 6);
    }
}

TEST_CASE("candidate_region_test is independent of the worker count") {
    // A block with a clear central peak, big enough for a window of 3 and a
    // peak width of 1.
    std::vector<std::vector<double>> dense(20, std::vector<double>(20, 0.0));
    for (std::size_t r = 0; r < 20; ++r) {
        for (std::size_t c = r + 1; c < 20; ++c) {
            dense[r][c] = 1.0 + static_cast<double>((r * 7 + c * 3) % 5);
        }
    }
    dense[8][12] = 40.0;
    const hicx::CsrMatrix matrix = from_dense(dense);
    const std::vector<hicx::loops::Candidate> candidates{
        {8, 12}, {4, 9}, {6, 15}, {2, 5}, {11, 17}};
    const hicx::loops::RegionTestResult one =
        hicx::loops::candidate_region_test(matrix, candidates, 3, 0.5, 1, 1);
    for (const unsigned int threads : {2U, 5U, 16U}) {
        const hicx::loops::RegionTestResult many =
            hicx::loops::candidate_region_test(matrix, candidates, 3, 0.5, 1, threads);
        REQUIRE(many.candidates.size() == one.candidates.size());
        for (std::size_t i = 0; i < one.candidates.size(); ++i) {
            CHECK(many.candidates[i].row == one.candidates[i].row);
            CHECK(many.candidates[i].col == one.candidates[i].col);
            // Byte identical, not merely close.
            CHECK(many.pvalues[i] == one.pvalues[i]);
        }
    }
}

TEST_CASE("the whole preselection is independent of the worker count") {
    const hicx::MatrixData data =
        hicx::read_hicexplorer_h5(std::string(kDataDir) + "/small_test_matrix.h5");
    // chr2L is bins 3462 to 8065 of that file; take a slice of it that is
    // large enough to have hundreds of distinct genomic distances.
    std::vector<std::int64_t> selection;
    for (std::int64_t bin = 3462; bin < 4462; ++bin) {
        selection.push_back(bin);
    }
    hicx::CsrMatrix block = hicx::select_bins(data.matrix, selection);
    // triu without the diagonal, which is what the tool works on.
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> cols;
    std::vector<double> values;
    for (std::int64_t row = 0; row < block.rows(); ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(block.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end =
            static_cast<std::size_t>(block.indptr()[static_cast<std::size_t>(row) + 1]);
        for (std::size_t k = begin; k < end; ++k) {
            if (block.indices()[k] <= row || block.data()[k] == 0.0) {
                continue;
            }
            rows.push_back(static_cast<std::int32_t>(row));
            cols.push_back(block.indices()[k]);
            values.push_back(block.data()[k]);
        }
    }
    REQUIRE(!values.empty());
    hicx::CsrMatrix upper = hicx::CsrMatrix::from_coo(
        block.rows(), block.cols(), rows, cols, std::move(values), block.dtype());
    upper.set_symmetry(hicx::Symmetry::Full);

    const hicx::loops::DistanceGroups groups =
        hicx::loops::group_by_distance(upper);
    const hicx::CsrMatrix obs_exp =
        hicx::loops::obs_exp_matrix_non_zero(upper, false);
    const hicx::loops::DistanceGroups obs_exp_groups =
        hicx::loops::group_by_distance(obs_exp);
    const std::map<std::int64_t, double> no_table;
    const hicx::loops::PreselectionResult one = hicx::loops::preselect_candidates(
        obs_exp, obs_exp_groups, 0.5, no_table, 5000, 1.5, 1);
    CHECK(std::accumulate(one.mask.begin(), one.mask.end(), 0) > 0);
    for (const unsigned int threads : {2U, 7U, 32U}) {
        const hicx::loops::PreselectionResult many =
            hicx::loops::preselect_candidates(obs_exp, obs_exp_groups, 0.5, no_table,
                                              5000, 1.5, threads);
        CHECK(many.mask == one.mask);
        REQUIRE(many.fits.size() == one.fits.size());
        for (std::size_t i = 0; i < one.fits.size(); ++i) {
            CHECK(many.fits[i].size == one.fits[i].size);
            CHECK(many.fits[i].prob == one.fits[i].prob);
        }
    }
    (void)groups;
}

TEST_CASE("the AVX2 reduction agrees with the scalar one on this tool's arrays") {
    // cpp/OPTIMIZATION.md 2: every SIMD kernel has a scalar reference and a
    // test asserting they agree. expected_interactions and numpy_mean are the
    // two places hicDetectLoops reduces, and both go through
    // hicx::simd::pairwise_sum, so the shapes tested here are the shapes it
    // actually feeds it: a per distance bucket of a few thousand values and a
    // peak or background window of a few dozen.
    if (!hicx::simd::avx2_available()) {
        MESSAGE("no AVX2 on this CPU, only the scalar kernel is exercised");
    }
    for (const std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{7},
                                std::size_t{25}, std::size_t{121}, std::size_t{128},
                                std::size_t{129}, std::size_t{1000},
                                std::size_t{4603}, std::size_t{8192},
                                std::size_t{8193}, std::size_t{20000}}) {
        std::vector<double> values(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            // Obs/exp values span several orders of magnitude, which is the
            // regime where a reordered reduction would show up.
            values[i] = 1e-6 + static_cast<double>((i * 37) % 9973) * 0.5;
        }
        const double scalar = hicx::simd::pairwise_sum_scalar(values.data(), n);
        const double dispatched = hicx::simd::pairwise_sum(values.data(), n);
        CHECK(dispatched == scalar);
        if (hicx::simd::avx2_available()) {
            CHECK(hicx::simd::pairwise_sum_avx2(values.data(), n) == scalar);
        }
        // numpy_mean is the reduction candidate_region_test compares on.
        if (n > 0) {
            CHECK(hicx::loops::numpy_mean(values) ==
                  scalar / static_cast<double>(n));
        }
    }
}
