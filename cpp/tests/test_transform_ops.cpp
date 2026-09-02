// Tests for the obs/exp, covariance, Pearson and eigenvector kernels that
// hicPCA and hicTransform share.
//
// The small matrix cases assert against numbers produced by the Python
// reference itself (numpy 1.26.4 and scipy 1.14.1 in the contract's
// environment), not against this implementation, so a change in either side is
// a failure rather than a silent drift. The real-matrix cases check the two
// properties no small example can: that the covariance is independent of the
// thread count, and that the SIMD and the scalar path of the row finalisation
// agree.

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "hicx/sparse_matrix.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"

using hicx::CsrMatrix;

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

// The symmetric integer matrix every small case below starts from:
//
//   10  5  2  0  1
//    5  8  4  3  0
//    2  4  9  6  2
//    0  3  6  7  5
//    1  0  2  5  6
//
// 21 stored entries, four exact zeros not stored.
CsrMatrix sample_matrix(const std::string& dtype = "float64") {
    const std::vector<std::vector<double>> dense{{10, 5, 2, 0, 1},
                                                 {5, 8, 4, 3, 0},
                                                 {2, 4, 9, 6, 2},
                                                 {0, 3, 6, 7, 5},
                                                 {1, 0, 2, 5, 6}};
    std::vector<std::int32_t> row;
    std::vector<std::int32_t> col;
    std::vector<double> data;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != 0.0) {
                row.push_back(i);
                col.push_back(j);
                data.push_back(dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
            }
        }
    }
    CsrMatrix matrix = CsrMatrix::from_coo(5, 5, row, col, std::move(data), dtype);
    // The kernels want the explicit symmetric form, which from_coo already
    // produces here because every mirrored entry was supplied.
    return matrix;
}

void check_close(const std::vector<double>& got, const std::vector<double>& want,
                 double tolerance = 1e-15) {
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < got.size(); ++i) {
        const double scale = std::max(1.0, std::fabs(want[i]));
        CHECK(std::fabs(got[i] - want[i]) / scale <= tolerance);
    }
}

}  // namespace

TEST_CASE("expected_interactions_in_distance matches utilities.py:293") {
    const CsrMatrix matrix = sample_matrix();
    // Python: expected_interactions_in_distance(5, 1, csr_matrix(D))
    check_close(hicx::expected_interactions_in_distance(matrix, 5, 1),
                {8.0, 10.0, 4.666666666666667, 0.0, 2.0});
}

TEST_CASE("expected_interactions_non_zero matches utilities.py:317") {
    const CsrMatrix matrix = sample_matrix();
    check_close(hicx::expected_interactions_non_zero(matrix),
                {8.0, 5.0, 2.3333333333333335, 0.0, 1.0});
}

TEST_CASE("expected_interactions matches utilities.py:356") {
    const CsrMatrix matrix = sample_matrix();
    // The divisor is np.arange(n + 1, 1, -1), that is n + 1 - i, not the
    // number of cells on the diagonal. Reproduced as it stands.
    check_close(hicx::expected_interactions(matrix),
                {6.666666666666667, 8.0, 3.5, 0.0, 1.0});
}

TEST_CASE("obs_exp_matrix_lieberman reads the expected array at ceil(d/2)") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    // Entry (0, 2) sits at distance 2 but is divided by expected[1], not by
    // expected[2], because utilities.py:499 indexes with ceil(|i-j| / 2). The
    // value is 2 / 10 = 0.2 and not 2 / 4.667 = 0.4286.
    CHECK(matrix.at(0, 2) == 0.2);
    // Entry (0, 4), distance 4, is divided by expected[2] = 4.6667.
    CHECK(matrix.at(0, 4) == doctest::Approx(0.21428571428571427).epsilon(1e-15));
    check_close(matrix.data(),
                {1.25, 0.5, 0.2, 0.21428571428571427,
                 0.5, 1.0, 0.4, 0.3,
                 0.2, 0.4, 1.125, 0.6, 0.2,
                 0.3, 0.6, 0.875, 0.5,
                 0.21428571428571427, 0.2, 0.5, 0.75});
    // astype(data_type) casts back to what the matrix was, so a float64
    // matrix stays float64 even though the division ran through float32.
    CHECK(matrix.dtype() == "float64");
}

TEST_CASE("obs_exp_matrix_non_zero keeps its result in float32") {
    SUBCASE("without the ligation factor") {
        CsrMatrix matrix = sample_matrix();
        hicx::obs_exp_non_zero_in_place(matrix, false);
        CHECK(matrix.dtype() == "float32");
        // Every value is the float32 rounding of the quotient, because the
        // Python assigns each one back into a float32 array one at a time.
        check_close(matrix.data(),
                    {1.25, 1.0, 0.8571428656578064, 1.0,
                     1.0, 1.0, 0.800000011920929, 1.2857142686843872,
                     0.8571428656578064, 0.800000011920929, 1.125, 1.2000000476837158,
                     0.8571428656578064, 1.2857142686843872, 1.2000000476837158, 0.875,
                     1.0, 1.0, 0.8571428656578064, 1.0, 0.75},
                    0.0);
    }
    SUBCASE("with the ligation factor") {
        CsrMatrix matrix = sample_matrix();
        hicx::obs_exp_non_zero_in_place(matrix, true);
        check_close(matrix.data(),
                    {0.37037035822868347, 0.2666666805744171, 0.19875776767730713,
                     0.380952388048172, 0.2666666805744171, 0.23999999463558197,
                     0.16695652902126312, 0.2938775420188904, 0.19875776767730713,
                     0.16695652902126312, 0.20415878295898438, 0.23850931227207184,
                     0.2555457055568695, 0.2938775420188904, 0.23850931227207184,
                     0.190476194024086, 0.3265306055545807, 0.380952388048172,
                     0.2555457055568695, 0.3265306055545807, 0.36734694242477417},
                    0.0);
    }
}

TEST_CASE("obs_exp_matrix matches utilities.py:554") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_in_place(matrix);
    check_close(matrix.data(),
                {1.5, 0.625, 0.25, 0.2857142857142857,
                 0.625, 1.2, 0.5, 0.375,
                 0.25, 0.5, 1.3499999999999999, 0.75, 0.25,
                 0.375, 0.75, 1.05, 0.625,
                 0.2857142857142857, 0.25, 0.625, 0.8999999999999999});
}

TEST_CASE("covariance_of_symmetric reproduces np.cov of the dense block") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    const hicx::DenseSymmetric covariance = hicx::covariance_of_symmetric(matrix, 1);
    // np.cov(csr_matrix(obs_exp).todense()) from the reference environment.
    const std::vector<double> expected{
        0.24039795918367346,  0.06317857142857143,   -0.09377678571428572,
        -0.1519017857142857,  -0.06295663265306123,  0.06317857142857143,
        0.133,                0.004750000000000005,  -0.049625,
        -0.09878571428571428, -0.09377678571428572,  0.004750000000000005,
        0.14762499999999998,  0.06778124999999999,   -0.030651785714285715,
        -0.1519017857142857,  -0.049625,             0.06778124999999999,
        0.107625,             0.0438125,             -0.06295663265306123,
        -0.09878571428571428, -0.030651785714285715, 0.0438125,
        0.08611224489795918};
    std::vector<double> got(covariance.data(), covariance.data() + 25);
    // The identity (A A^T - n m m^T) / (n - 1) is algebraically np.cov's
    // centred product but not the same arithmetic, so the agreement is at the
    // level of the cancellation it introduces, not bit for bit. Measured on
    // the real validation inputs: 2.1e-12 pure relative on the covariance of
    // small_test_matrix.h5 chrX and 4.3e-14 relative to a floor of 1 on
    // mm9_reduced_chr1.cool.
    check_close(got, expected, 1e-13);
}

TEST_CASE("pearson_row reproduces np.corrcoef of the dense block") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    const hicx::DenseSymmetric covariance = hicx::covariance_of_symmetric(matrix, 1);
    const std::vector<double> scaling = hicx::pearson_scaling(covariance);
    std::vector<double> got;
    std::vector<double> row(5);
    for (std::int64_t i = 0; i < 5; ++i) {
        hicx::pearson_row(covariance, scaling, i, row.data());
        got.insert(got.end(), row.begin(), row.end());
    }
    const std::vector<double> expected{
        1.0,
        0.35332843748576254,
        -0.4977943915337237,
        -0.9443672309251671,
        -0.43756607887516746,
        0.3533284374857626,
        1.0,
        0.0338990739077925,
        -0.41478029587074233,
        -0.9230726181963105,
        -0.4977943915337236,
        0.0338990739077925,
        1.0,
        0.5377408584507045,
        -0.2718587839009439,
        -0.9443672309251671,
        -0.4147802958707423,
        0.5377408584507045,
        1.0,
        0.45510226555565353,
        -0.43756607887516746,
        -0.9230726181963105,
        -0.2718587839009439,
        0.45510226555565353,
        1.0};
    check_close(got, expected, 1e-12);
}

TEST_CASE("dsyevr returns the requested eigenvectors in descending order") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    hicx::DenseSymmetric covariance = hicx::covariance_of_symmetric(matrix, 1);
    hicx::pin_blas_to_one_thread();
    const hicx::EigenResult result =
        hicx::leading_eigenvectors(covariance, {1, 2, 3}, hicx::EigenSolver::Dsyevr);
    // numpy.linalg.eigh of the same covariance, largest first.
    CHECK(result.values[0] == doctest::Approx(0.4331053441544388).epsilon(1e-11));
    CHECK(result.values[1] == doctest::Approx(0.20799007375990344).epsilon(1e-11));
    CHECK(result.values[2] == doctest::Approx(0.0664921673316918).epsilon(1e-11));
    for (const std::vector<double>& vector : result.vectors) {
        REQUIRE(vector.size() == 5);
        // The sign convention: the largest magnitude component is positive.
        std::size_t extreme = 0;
        for (std::size_t j = 1; j < vector.size(); ++j) {
            if (std::fabs(vector[j]) > std::fabs(vector[extreme])) {
                extreme = j;
            }
        }
        CHECK(vector[extreme] > 0.0);
        double norm = 0.0;
        for (double value : vector) {
            norm += value * value;
        }
        CHECK(std::sqrt(norm) == doctest::Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("dgeev and dsyevr agree on the eigenpairs up to order and sign") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    hicx::pin_blas_to_one_thread();

    hicx::DenseSymmetric a = hicx::covariance_of_symmetric(matrix, 1);
    hicx::DenseSymmetric b = hicx::covariance_of_symmetric(matrix, 1);
    const hicx::EigenResult general =
        hicx::leading_eigenvectors(a, {1, 2, 3, 4, 5}, hicx::EigenSolver::Dgeev);
    const hicx::EigenResult symmetric =
        hicx::leading_eigenvectors(b, {1, 2, 3, 4, 5}, hicx::EigenSolver::Dsyevr);

    // Every eigenvalue the symmetric solver returns is one the general solver
    // returned too, whatever order dgeev happened to use.
    for (double wanted : symmetric.values) {
        bool found = false;
        for (double got : general.values) {
            if (std::fabs(got - wanted) <= 1e-10 * std::max(1.0, std::fabs(wanted))) {
                found = true;
            }
        }
        CHECK(found);
    }
}

TEST_CASE("a requested eigenvector beyond the matrix size comes back empty") {
    CsrMatrix matrix = sample_matrix();
    hicx::obs_exp_lieberman_in_place(matrix, 5, 1);
    hicx::pin_blas_to_one_thread();
    hicx::DenseSymmetric covariance = hicx::covariance_of_symmetric(matrix, 1);
    const hicx::EigenResult result =
        hicx::leading_eigenvectors(covariance, {1, 9}, hicx::EigenSolver::Dgeev);
    CHECK(result.vectors[0].size() == 5);
    // numpy's out-of-range column slice is empty, which shortens the row and
    // makes hicPCA skip the whole chromosome when writing.
    CHECK(result.vectors[1].empty());
}

TEST_CASE("the SIMD row finalisation agrees with the scalar reference") {
    // cpp/OPTIMIZATION.md section 2 requires a test that every dispatched path
    // agrees with the always-compiled scalar one, and section 3 says the
    // comparison is at the ED tolerance rather than at equality, because the
    // fused multiply-add rounds once where the scalar expression rounds twice.
    const std::size_t n = 1031;  // deliberately not a multiple of 4 or 8
    std::vector<double> accumulator(n);
    std::vector<double> means(n);
    for (std::size_t j = 0; j < n; ++j) {
        // A spread of magnitudes, deterministic, no random number generator.
        accumulator[j] = std::sin(static_cast<double>(j) * 0.37) *
                         std::pow(10.0, static_cast<double>(j % 9) - 4.0);
        means[j] = std::cos(static_cast<double>(j) * 0.11) * 1e-3;
    }
    const double mean_i = 7.25e-4;
    const double observation_count = static_cast<double>(n);
    const double inverse = 1.0 / static_cast<double>(n - 1);

    std::vector<double> scalar(n);
    std::vector<double> dispatched(n);
    hicx::covariance_row_finalise_scalar(scalar.data(), accumulator.data(), means.data(),
                                         mean_i, observation_count, inverse, n);
    hicx::covariance_row_finalise(dispatched.data(), accumulator.data(), means.data(),
                                  mean_i, observation_count, inverse, n);
    INFO("SIMD path: " << hicx::simd_path_name());
    for (std::size_t j = 0; j < n; ++j) {
        if (scalar[j] == 0.0) {
            CHECK(dispatched[j] == 0.0);
            continue;
        }
        CHECK(std::fabs(dispatched[j] - scalar[j]) / std::fabs(scalar[j]) <= 1e-3);
    }
    // In practice they are far closer than the gate: the only difference is
    // the fused multiply-add, so the agreement is at the ulp level.
    for (std::size_t j = 0; j < n; ++j) {
        const double scale = std::max(1e-30, std::fabs(scalar[j]));
        CHECK(std::fabs(dispatched[j] - scalar[j]) / scale <= 1e-14);
    }
    // The in-place form, which is how covariance_of_symmetric calls it.
    std::vector<double> in_place = accumulator;
    hicx::covariance_row_finalise(in_place.data(), in_place.data(), means.data(), mean_i,
                                  observation_count, inverse, n);
    for (std::size_t j = 0; j < n; ++j) {
        CHECK(in_place[j] == dispatched[j]);
    }
}

TEST_CASE("the covariance does not depend on the thread count") {
    // The determinism rule of cpp/OPTIMIZATION.md section 3, on a real
    // chromosome block rather than on a toy: 449 bins of chrX from the 50 kb
    // matrix, which is small enough for a unit test and large enough that the
    // row partition really is split across the workers.
    hicx::ToolMatrix hic =
        hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    hicx::keep_only_chromosomes(hic.data(), {"chrX"});
    CsrMatrix block = hic.matrix();
    block.materialize_full();
    hicx::obs_exp_lieberman_in_place(block, block.rows(), 1);

    const hicx::DenseSymmetric one = hicx::covariance_of_symmetric(block, 1);
    const hicx::DenseSymmetric four = hicx::covariance_of_symmetric(block, 4);
    const hicx::DenseSymmetric sixteen = hicx::covariance_of_symmetric(block, 16);
    REQUIRE(one.size() == four.size());
    REQUIRE(one.size() == sixteen.size());
    const std::size_t count =
        static_cast<std::size_t>(one.size()) * static_cast<std::size_t>(one.size());
    bool identical_four = true;
    bool identical_sixteen = true;
    for (std::size_t k = 0; k < count; ++k) {
        // Byte identical, not merely close. A parallel reduction whose result
        // depends on the thread count is a defect even when it passes ED.
        identical_four = identical_four && (one.data()[k] == four.data()[k]);
        identical_sixteen = identical_sixteen && (one.data()[k] == sixteen.data()[k]);
    }
    CHECK(identical_four);
    CHECK(identical_sixteen);
    // And it really is symmetric, which is what lets the same buffer be handed
    // to LAPACK as a column major matrix without a transpose.
    bool symmetric = true;
    for (std::int64_t i = 0; i < one.size(); ++i) {
        for (std::int64_t j = 0; j < one.size(); ++j) {
            symmetric = symmetric && (one.at(i, j) == one.at(j, i));
        }
    }
    CHECK(symmetric);
}

TEST_CASE("keep_only_chromosomes selects in ascending bin order") {
    hicx::ToolMatrix hic =
        hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    const std::size_t before = hic.data().cut_intervals.size();
    // Named in the reverse of their order in the file: the selection must
    // still come out in bin order, because the Python takes np.flatnonzero of
    // a boolean mask.
    hicx::keep_only_chromosomes(hic.data(), {"chrX", "chr2RHet"});
    const std::vector<hicx::CutInterval>& intervals = hic.data().cut_intervals;
    CHECK(intervals.size() < before);
    CHECK(intervals.front().chrom == "chr2RHet");
    CHECK(intervals.back().chrom == "chrX");
    CHECK(hic.matrix().rows() == static_cast<std::int64_t>(intervals.size()));
    // keepOnlyTheseChr clears distance_counts unconditionally (HiCMatrix.py:677).
    CHECK(!hic.data().distance_counts.has_value());
}

TEST_CASE("keep_only_chromosomes rejects a name that is not in the matrix") {
    hicx::ToolMatrix hic =
        hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    CHECK_THROWS(hicx::keep_only_chromosomes(hic.data(), {"chrNotThere"}));
}
