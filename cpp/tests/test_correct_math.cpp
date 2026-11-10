// Unit tests for the hicCorrectMatrix numerics: the sparse kernels and their
// SIMD paths, the modified z-score filter, the bin masking and the two
// balancing algorithms.
//
// The kernels are checked against dense reference implementations written
// straight out of the definitions, which is the only way to be sure that the
// upper triangle storage and the partitioned reductions really produce the
// symmetric matrix's values. The SIMD paths are checked against the scalar path
// on the same input, bit for bit where the operation is elementwise.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <doctest/doctest.h>

#include "hicx/correct_ops.hpp"
#include "hicx/math/ice.hpp"
#include "hicx/math/kr_balancing.hpp"
#include "hicx/math/sparse_kernels.hpp"
#include "hicx/sparse_matrix.hpp"

namespace {

using hicx::CsrMatrix;
using hicx::Symmetry;

// A random symmetric matrix, held as the upper triangle, plus its dense form.
struct Fixture {
    CsrMatrix matrix;
    std::vector<double> dense;  // n * n, symmetric
    std::int64_t n = 0;
};

Fixture make_symmetric(std::int64_t n, double density, unsigned seed,
                       bool full_diagonal = true) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> value(0.5, 20.0);
    std::uniform_real_distribution<double> chance(0.0, 1.0);

    Fixture fixture;
    fixture.n = n;
    fixture.dense.assign(static_cast<std::size_t>(n * n), 0.0);
    std::vector<std::int64_t> indptr(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::int32_t> indices;
    std::vector<double> data;
    for (std::int64_t row = 0; row < n; ++row) {
        for (std::int64_t column = row; column < n; ++column) {
            const bool diagonal = row == column;
            if (!(diagonal && full_diagonal) && chance(rng) > density) {
                continue;
            }
            const double entry = value(rng);
            indices.push_back(static_cast<std::int32_t>(column));
            data.push_back(entry);
            fixture.dense[static_cast<std::size_t>(row * n + column)] = entry;
            fixture.dense[static_cast<std::size_t>(column * n + row)] = entry;
        }
        indptr[static_cast<std::size_t>(row) + 1] = static_cast<std::int64_t>(data.size());
    }
    fixture.matrix = CsrMatrix(n, n, std::move(indptr), std::move(indices), std::move(data),
                               "float64");
    fixture.matrix.set_symmetry(Symmetry::UpperTriangle);
    return fixture;
}

std::vector<double> dense_marginals(const Fixture& fixture) {
    std::vector<double> out(static_cast<std::size_t>(fixture.n), 0.0);
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        double sum = 0.0;
        for (std::int64_t column = 0; column < fixture.n; ++column) {
            sum += fixture.dense[static_cast<std::size_t>(row * fixture.n + column)];
        }
        out[static_cast<std::size_t>(row)] = sum;
    }
    return out;
}

// Guard so that a test that changes the split threshold cannot leak it.
struct SplitThresholdGuard {
    std::size_t previous = hicx::kernels::block_split_threshold();
    explicit SplitThresholdGuard(std::size_t entries) {
        hicx::kernels::set_block_split_threshold(entries);
    }
    ~SplitThresholdGuard() { hicx::kernels::set_block_split_threshold(previous); }
};

struct SimdGuard {
    hicx::kernels::SimdPath previous = hicx::kernels::active_simd_path();
    explicit SimdGuard(hicx::kernels::SimdPath path) {
        hicx::kernels::force_simd_path(path);
    }
    ~SimdGuard() { hicx::kernels::force_simd_path(previous); }
};

}  // namespace

TEST_CASE("numpy_median matches numpy on odd, even and empty input") {
    CHECK(hicx::correct::numpy_median({3.0, 1.0, 2.0}) == doctest::Approx(2.0));
    CHECK(hicx::correct::numpy_median({4.0, 1.0, 3.0, 2.0}) == doctest::Approx(2.5));
    CHECK(hicx::correct::numpy_median({7.0}) == doctest::Approx(7.0));
    CHECK(std::isnan(hicx::correct::numpy_median({})));
    // The even case is the mean of the two middle elements, not the lower one.
    CHECK(hicx::correct::numpy_median({1.0, 1.0, 5.0, 5.0}) == doctest::Approx(3.0));
}

TEST_CASE("MAD reproduces the modified z-score of hicCorrectMatrix") {
    // median over the positive entries is 4, diff is {-4, -3, -1, 0, 2, 6},
    // |diff| has median (1 + 2) / 2 = 1.5.
    const std::vector<double> points = {0.0, 1.0, 3.0, 4.0, 6.0, 10.0};
    const hicx::correct::Mad mad(points);
    // The median is taken over the strictly positive entries only, so the zero
    // is excluded and the median of 1, 3, 4, 6, 10 is 4.
    CHECK(mad.median() == doctest::Approx(4.0));
    CHECK(mad.median_absolute_deviation() == doctest::Approx(2.5));
    const hicx::correct::Mad recomputed(points);
    const std::vector<double>& z = recomputed.modified_z_scores();
    REQUIRE(z.size() == points.size());
    // Whatever the median turns out to be, the definition has to hold exactly.
    const double median = recomputed.median();
    const double deviation = recomputed.median_absolute_deviation();
    for (std::size_t i = 0; i < points.size(); ++i) {
        CHECK(z[i] == doctest::Approx(0.6745 * ((points[i] - median) / deviation)));
    }
}

TEST_CASE("symmetric marginals equal the dense row sums") {
    Fixture fixture = make_symmetric(37, 0.3, 11);
    hicx::kernels::DiagonalBlock block(fixture.matrix);
    std::vector<double> out(static_cast<std::size_t>(fixture.n), 0.0);
    hicx::kernels::symmetric_marginals(block, out.data(), 1);
    const std::vector<double> expected = dense_marginals(fixture);
    for (std::size_t i = 0; i < out.size(); ++i) {
        CHECK(out[i] == doctest::Approx(expected[i]).epsilon(1e-12));
    }
}

TEST_CASE("a diagonal block sees only its own chromosome") {
    Fixture fixture = make_symmetric(24, 0.5, 5);
    const std::int64_t first = 8;
    const std::int64_t last = 17;
    hicx::kernels::DiagonalBlock block(fixture.matrix, first, last);
    REQUIRE(block.size() == last - first);
    std::vector<double> out(static_cast<std::size_t>(block.size()), 0.0);
    hicx::kernels::symmetric_marginals(block, out.data(), 1);
    for (std::int64_t local = 0; local < block.size(); ++local) {
        double expected = 0.0;
        for (std::int64_t column = first; column < last; ++column) {
            expected += fixture.dense[static_cast<std::size_t>((first + local) * fixture.n +
                                                               column)];
        }
        CHECK(out[static_cast<std::size_t>(local)] == doctest::Approx(expected).epsilon(1e-12));
    }
}

TEST_CASE("partitioning does not depend on the thread count") {
    SplitThresholdGuard guard(1);  // force the partitioned path
    Fixture fixture = make_symmetric(200, 0.2, 7);
    hicx::kernels::DiagonalBlock block(fixture.matrix);
    REQUIRE(block.partitions() > 1);
    // Partition boundaries are monotone and cover every row exactly once.
    CHECK(block.partition().front() == 0);
    CHECK(block.partition().back() == fixture.n);
    for (std::size_t p = 1; p < block.partition().size(); ++p) {
        CHECK(block.partition()[p] >= block.partition()[p - 1]);
    }

    std::vector<double> one(static_cast<std::size_t>(fixture.n), 0.0);
    std::vector<double> many(static_cast<std::size_t>(fixture.n), 0.0);
    hicx::kernels::symmetric_marginals(block, one.data(), 1);
    hicx::kernels::symmetric_marginals(block, many.data(), 8);
    for (std::size_t i = 0; i < one.size(); ++i) {
        // Bit identical, not merely close: the partitions are the same and are
        // combined in the same order whatever the thread count.
        CHECK(one[i] == many[i]);
    }
}

TEST_CASE("symmetric spmv equals the dense product, identity term included") {
    Fixture fixture = make_symmetric(41, 0.25, 3);
    std::mt19937 rng(99);
    std::uniform_real_distribution<double> spread(0.1, 2.0);
    std::vector<double> x(static_cast<std::size_t>(fixture.n));
    for (double& value : x) {
        value = spread(rng);
    }
    const double addend = 0.00001;

    hicx::kernels::DiagonalBlock block(fixture.matrix);
    std::vector<double> y(static_cast<std::size_t>(fixture.n), 0.0);
    hicx::kernels::symmetric_spmv(block, x.data(), y.data(), addend, 1);
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        double expected = addend * x[static_cast<std::size_t>(row)];
        for (std::int64_t column = 0; column < fixture.n; ++column) {
            expected += fixture.dense[static_cast<std::size_t>(row * fixture.n + column)] *
                        x[static_cast<std::size_t>(column)];
        }
        CHECK(y[static_cast<std::size_t>(row)] ==
              doctest::Approx(expected).epsilon(1e-12));
    }
}

TEST_CASE("the fused ICE pass agrees with scaling and then reducing") {
    SplitThresholdGuard guard(1);
    Fixture reference = make_symmetric(120, 0.2, 21);
    Fixture fused = make_symmetric(120, 0.2, 21);
    std::mt19937 rng(4);
    std::uniform_real_distribution<double> spread(0.2, 3.0);
    std::vector<double> scale(static_cast<std::size_t>(reference.n));
    for (double& value : scale) {
        value = spread(rng);
    }

    hicx::kernels::DiagonalBlock reference_block(reference.matrix);
    const double reference_largest =
        hicx::kernels::scale_rows_and_cols(reference_block, scale.data(), 1);
    std::vector<double> reference_marginals(static_cast<std::size_t>(reference.n), 0.0);
    hicx::kernels::symmetric_marginals(reference_block, reference_marginals.data(), 1);

    hicx::kernels::DiagonalBlock fused_block(fused.matrix);
    std::vector<double> fused_marginals(static_cast<std::size_t>(fused.n), 0.0);
    const double fused_largest = hicx::kernels::scale_and_marginals(
        fused_block, scale.data(), fused_marginals.data(), 1);

    CHECK(fused_largest == reference_largest);
    for (std::size_t k = 0; k < reference.matrix.data().size(); ++k) {
        CHECK(fused.matrix.data()[k] == reference.matrix.data()[k]);
    }
    for (std::size_t i = 0; i < fused_marginals.size(); ++i) {
        CHECK(fused_marginals[i] == reference_marginals[i]);
    }
}

TEST_CASE("the AVX2 scaling kernels agree with the scalar reference") {
    if (hicx::kernels::active_simd_path() == hicx::kernels::SimdPath::Scalar) {
        MESSAGE("no AVX2 on this machine, the scalar path is the only one");
        return;
    }
    Fixture scalar = make_symmetric(97, 0.4, 31);
    Fixture vector = make_symmetric(97, 0.4, 31);
    std::mt19937 rng(8);
    std::uniform_real_distribution<double> spread(0.3, 4.0);
    std::vector<double> scale(static_cast<std::size_t>(scalar.n));
    for (double& value : scale) {
        value = spread(rng);
    }

    double scalar_largest = 0.0;
    {
        SimdGuard guard(hicx::kernels::SimdPath::Scalar);
        hicx::kernels::DiagonalBlock block(scalar.matrix);
        scalar_largest = hicx::kernels::scale_rows_and_cols(block, scale.data(), 1);
    }
    double vector_largest = 0.0;
    {
        SimdGuard guard(hicx::kernels::SimdPath::Avx2);
        hicx::kernels::DiagonalBlock block(vector.matrix);
        vector_largest = hicx::kernels::scale_rows_and_cols(block, scale.data(), 1);
    }
    // Elementwise, so the vectorised path performs the same two multiplications
    // in the same pairing and the results are bit identical.
    CHECK(scalar_largest == vector_largest);
    for (std::size_t k = 0; k < scalar.matrix.data().size(); ++k) {
        CHECK(vector.matrix.data()[k] == scalar.matrix.data()[k]);
    }

    Fixture scalar_symmetric = make_symmetric(83, 0.35, 55);
    Fixture vector_symmetric = make_symmetric(83, 0.35, 55);
    std::vector<double> x(static_cast<std::size_t>(scalar_symmetric.n));
    for (double& value : x) {
        value = spread(rng);
    }
    {
        SimdGuard guard(hicx::kernels::SimdPath::Scalar);
        hicx::kernels::DiagonalBlock block(scalar_symmetric.matrix);
        hicx::kernels::scale_symmetric(block, x.data(), 1);
    }
    {
        SimdGuard guard(hicx::kernels::SimdPath::Avx2);
        hicx::kernels::DiagonalBlock block(vector_symmetric.matrix);
        hicx::kernels::scale_symmetric(block, x.data(), 1);
    }
    for (std::size_t k = 0; k < scalar_symmetric.matrix.data().size(); ++k) {
        CHECK(vector_symmetric.matrix.data()[k] == scalar_symmetric.matrix.data()[k]);
    }
}

TEST_CASE("ICE drives the marginals to one and is thread count invariant") {
    SplitThresholdGuard guard(1);
    Fixture one = make_symmetric(150, 0.3, 77);
    Fixture many = make_symmetric(150, 0.3, 77);

    hicx::ice::Options options;
    options.max_iterations = 200;
    options.threads = 1;
    hicx::kernels::DiagonalBlock block_one(one.matrix);
    const hicx::ice::Result result_one = hicx::ice::correct(block_one, options);
    options.threads = 8;
    hicx::kernels::DiagonalBlock block_many(many.matrix);
    const hicx::ice::Result result_many = hicx::ice::correct(block_many, options);

    CHECK(result_one.converged);
    CHECK(result_one.iterations == result_many.iterations);
    for (std::size_t i = 0; i < result_one.correction_factors.size(); ++i) {
        CHECK(result_one.correction_factors[i] == result_many.correction_factors[i]);
    }
    for (std::size_t k = 0; k < one.matrix.data().size(); ++k) {
        CHECK(many.matrix.data()[k] == one.matrix.data()[k]);
    }

    // The defining property: after correction every marginal is the same.
    std::vector<double> marginals(static_cast<std::size_t>(one.n), 0.0);
    hicx::kernels::symmetric_marginals(block_one, marginals.data(), 1);
    const double first = marginals.front();
    for (const double value : marginals) {
        CHECK(value == doctest::Approx(first).epsilon(1e-4));
    }
}

TEST_CASE("Knight-Ruiz balances the matrix and is thread count invariant") {
    SplitThresholdGuard guard(1);
    Fixture one = make_symmetric(120, 0.35, 13);
    Fixture many = make_symmetric(120, 0.35, 13);

    hicx::kernels::DiagonalBlock block_one(one.matrix);
    hicx::kr::Options options;
    options.threads = 1;
    hicx::kr::Balancer balancer_one(block_one, options);
    REQUIRE(balancer_one.compute());

    hicx::kernels::DiagonalBlock block_many(many.matrix);
    options.threads = 8;
    hicx::kr::Balancer balancer_many(block_many, options);
    REQUIRE(balancer_many.compute());

    const std::vector<double> x_one = balancer_one.normalisation_vector(false);
    const std::vector<double> x_many = balancer_many.normalisation_vector(false);
    for (std::size_t i = 0; i < x_one.size(); ++i) {
        CHECK(x_one[i] == x_many[i]);
    }

    // The balancing condition: diag(x) A diag(x) has unit row sums.
    for (std::int64_t row = 0; row < one.n; ++row) {
        double sum = 0.00001 * x_one[static_cast<std::size_t>(row)] *
                     x_one[static_cast<std::size_t>(row)];
        for (std::int64_t column = 0; column < one.n; ++column) {
            sum += one.dense[static_cast<std::size_t>(row * one.n + column)] *
                   x_one[static_cast<std::size_t>(row)] *
                   x_one[static_cast<std::size_t>(column)];
        }
        CHECK(sum == doctest::Approx(1.0).epsilon(1e-5));
    }
}

TEST_CASE("the rescaled normalisation vector preserves the matrix sum") {
    Fixture fixture = make_symmetric(60, 0.4, 91);
    hicx::kernels::DiagonalBlock block(fixture.matrix);
    hicx::kr::Options options;
    hicx::kr::Balancer balancer(block, options);
    REQUIRE(balancer.compute());
    CHECK_FALSE(balancer.rescaled());
    const std::vector<double> vector = balancer.normalisation_vector(true);
    CHECK(balancer.rescaled());
    CHECK(std::isfinite(balancer.normalisation_factor()));

    // rescale_norm_vector divides x by sqrt(sum(A_ij x_i x_j) / sum(A_ij)) over
    // the upper triangle with the off diagonal counted twice, so the balanced
    // matrix has the same total as the raw one.
    double original = 0.0;
    double balanced = 0.0;
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        for (std::int64_t column = 0; column < fixture.n; ++column) {
            const double value =
                fixture.dense[static_cast<std::size_t>(row * fixture.n + column)] +
                (row == column ? 0.00001 : 0.0);
            original += value;
            balanced += value * vector[static_cast<std::size_t>(row)] *
                        vector[static_cast<std::size_t>(column)];
        }
    }
    CHECK(balanced == doctest::Approx(original).epsilon(1e-6));
}

TEST_CASE("v3 float32 accumulators move the normalisation factor") {
    Fixture wide = make_symmetric(70, 0.5, 17);
    Fixture narrow = make_symmetric(70, 0.5, 17);

    hicx::kernels::DiagonalBlock wide_block(wide.matrix);
    hicx::kr::Options wide_options;
    hicx::kr::Balancer wide_balancer(wide_block, wide_options);
    REQUIRE(wide_balancer.compute());
    wide_balancer.normalisation_vector(true);

    hicx::kernels::DiagonalBlock narrow_block(narrow.matrix);
    hicx::kr::Options narrow_options;
    narrow_options.float32_input = true;
    narrow_options.float32_rescale = true;
    hicx::kr::Balancer narrow_balancer(narrow_block, narrow_options);
    REQUIRE(narrow_balancer.compute());
    narrow_balancer.normalisation_vector(true);

    // Same algorithm, different precision: close but not equal.
    CHECK(narrow_balancer.normalisation_factor() ==
          doctest::Approx(wide_balancer.normalisation_factor()).epsilon(1e-3));
    CHECK(narrow_balancer.normalisation_factor() != wide_balancer.normalisation_factor());
}

TEST_CASE("KR reports non convergence instead of exiting the process") {
    // krbalancing calls exit(0) after 300 outer iterations. One iteration is
    // never enough for a real matrix, so the guard is easy to reach.
    Fixture fixture = make_symmetric(80, 0.3, 23);
    hicx::kernels::DiagonalBlock block(fixture.matrix);
    hicx::kr::Options options;
    options.max_outer_iterations = 1;
    hicx::kr::Balancer balancer(block, options);
    CHECK_FALSE(balancer.compute());
    CHECK(balancer.message().find("did not converge") != std::string::npos);
}

TEST_CASE("add_missing_diagonal completes the diagonal and nothing else") {
    Fixture fixture = make_symmetric(30, 0.2, 43, /*full_diagonal=*/false);
    const std::size_t before = fixture.matrix.data().size();
    const std::size_t inserted = hicx::correct::add_missing_diagonal(fixture.matrix, 7.0);
    CHECK(fixture.matrix.data().size() == before + inserted);
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        const std::size_t begin =
            static_cast<std::size_t>(fixture.matrix.indptr()[static_cast<std::size_t>(row)]);
        const std::size_t end = static_cast<std::size_t>(
            fixture.matrix.indptr()[static_cast<std::size_t>(row) + 1]);
        REQUIRE(begin < end);
        CHECK(fixture.matrix.indices()[begin] == static_cast<std::int32_t>(row));
        const double expected =
            fixture.dense[static_cast<std::size_t>(row * fixture.n + row)];
        CHECK(fixture.matrix.data()[begin] == (expected == 0.0 ? 7.0 : expected));
        // The columns of a row stay sorted.
        for (std::size_t k = begin + 1; k < end; ++k) {
            CHECK(fixture.matrix.indices()[k] > fixture.matrix.indices()[k - 1]);
        }
    }
    // Idempotent.
    CHECK(hicx::correct::add_missing_diagonal(fixture.matrix, 7.0) == 0);
}

TEST_CASE("remove_diagonal drops exactly the diagonal") {
    Fixture fixture = make_symmetric(25, 0.3, 61);
    const std::size_t before = fixture.matrix.data().size();
    hicx::correct::remove_diagonal(fixture.matrix);
    CHECK(fixture.matrix.data().size() == before - static_cast<std::size_t>(fixture.n));
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        for (std::size_t k =
                 static_cast<std::size_t>(fixture.matrix.indptr()[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(
                     fixture.matrix.indptr()[static_cast<std::size_t>(row) + 1]);
             ++k) {
            CHECK(fixture.matrix.indices()[k] != static_cast<std::int32_t>(row));
        }
    }
}

TEST_CASE("masking and restoring bins is a round trip") {
    Fixture fixture = make_symmetric(20, 0.4, 71);
    hicx::MatrixData data;
    data.matrix = fixture.matrix;
    for (std::int64_t bin = 0; bin < fixture.n; ++bin) {
        data.cut_intervals.push_back(
            hicx::CutInterval{"chr1", bin * 10, bin * 10 + 10, 1.0, ""});
    }

    std::vector<char> masked(static_cast<std::size_t>(fixture.n), 0);
    masked[3] = 1;
    masked[4] = 1;
    masked[17] = 1;
    const hicx::correct::MaskState state =
        hicx::correct::mask_bins_in_place(data, masked);
    CHECK(data.matrix.rows() == fixture.n - 3);
    CHECK(state.removed == std::vector<std::int64_t>{3, 4, 17});
    CHECK(data.cut_intervals.size() == static_cast<std::size_t>(fixture.n - 3));

    data.correction_factors = std::vector<double>(
        static_cast<std::size_t>(data.matrix.rows()), 2.0);
    hicx::correct::restore_masked_bins(data, state);
    CHECK(data.matrix.rows() == fixture.n);
    CHECK(data.nan_bins == std::vector<std::int64_t>{3, 4, 17});
    REQUIRE(data.correction_factors.has_value());
    CHECK(std::isnan((*data.correction_factors)[3]));
    CHECK(std::isnan((*data.correction_factors)[17]));
    CHECK((*data.correction_factors)[0] == 2.0);

    // Every surviving entry is back at its original coordinates with its
    // original value; every masked row and column is empty.
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        for (std::size_t k =
                 static_cast<std::size_t>(data.matrix.indptr()[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(
                     data.matrix.indptr()[static_cast<std::size_t>(row) + 1]);
             ++k) {
            const std::int64_t column = data.matrix.indices()[k];
            CHECK(masked[static_cast<std::size_t>(row)] == 0);
            CHECK(masked[static_cast<std::size_t>(column)] == 0);
            CHECK(data.matrix.data()[k] ==
                  fixture.dense[static_cast<std::size_t>(row * fixture.n + column)]);
        }
    }
}

TEST_CASE("keep_only_diagonal_blocks drops the inter chromosomal contacts") {
    Fixture fixture = make_symmetric(12, 0.9, 83);
    std::vector<std::pair<std::string, hicx::BinRange>> boundaries = {
        {"chrA", hicx::BinRange{0, 5}}, {"chrB", hicx::BinRange{5, 12}}};
    hicx::correct::keep_only_diagonal_blocks(fixture.matrix, boundaries);
    for (std::int64_t row = 0; row < fixture.n; ++row) {
        const std::int64_t limit = row < 5 ? 5 : 12;
        for (std::size_t k =
                 static_cast<std::size_t>(fixture.matrix.indptr()[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(
                     fixture.matrix.indptr()[static_cast<std::size_t>(row) + 1]);
             ++k) {
            CHECK(fixture.matrix.indices()[k] < limit);
        }
    }
}

TEST_CASE("CsrMatrix release and adopt move the arrays without copying") {
    Fixture fixture = make_symmetric(10, 0.5, 101);
    const double* values = fixture.matrix.data().data();
    CsrMatrix::Arrays arrays = fixture.matrix.release();
    CHECK(fixture.matrix.rows() == 0);
    CHECK(fixture.matrix.data().empty());
    CHECK(arrays.data.data() == values);
    CHECK(arrays.symmetry == Symmetry::UpperTriangle);
    CsrMatrix restored = CsrMatrix::adopt(std::move(arrays));
    CHECK(restored.rows() == 10);
    CHECK(restored.data().data() == values);
    CHECK(restored.symmetry() == Symmetry::UpperTriangle);
}

TEST_CASE("cpython_set_order reproduces CPython's set iteration order") {
    // The 64 MAD outliers of small_test_matrix.h5 restricted to chrUextra and
    // chr3LHet, in the ascending order sorted() hands to set(), and the order
    // CPython 3.12 then iterates them in. Both were read off the reference
    // implementation; the second is the line order of the checked in
    // hicCorrectMatrix/filtered.bed.
    const std::vector<std::int64_t> ascending = {
        3, 11, 23, 27, 28, 29, 39, 41, 43, 45, 46, 48, 49, 52, 56, 57, 59, 63,
        64, 65, 67, 68, 73, 74, 78, 80, 81, 82, 87, 89, 90, 91, 94, 96, 99, 100,
        101, 102, 104, 109, 110, 111, 113, 114, 117, 120, 121, 128, 136, 141,
        150, 153, 155, 156, 157, 159, 165, 167, 168, 169, 171, 174, 177, 179};
    const std::vector<std::int64_t> expected = {
        128, 3, 136, 11, 141, 150, 23, 153, 27, 28, 29, 155, 156, 157, 159, 165,
        39, 167, 41, 168, 43, 169, 45, 46, 171, 48, 49, 174, 177, 52, 179, 56,
        57, 59, 63, 64, 65, 67, 68, 73, 74, 78, 80, 81, 82, 87, 89, 90, 91, 94,
        96, 99, 100, 101, 102, 104, 109, 110, 111, 113, 114, 117, 120, 121};
    CHECK(hicx::correct::cpython_set_order(ascending) == expected);

    // Small sets: with eight slots and one resize the order is still not the
    // insertion order.
    CHECK(hicx::correct::cpython_set_order({1, 2, 3}) ==
          std::vector<std::int64_t>{1, 2, 3});
    CHECK(hicx::correct::cpython_set_order({8, 1}) == std::vector<std::int64_t>{8, 1});
    CHECK(hicx::correct::cpython_set_order({}).empty());
    // A duplicate is absorbed, as a set does.
    CHECK(hicx::correct::cpython_set_order({5, 5, 5}) == std::vector<std::int64_t>{5});
}
