// Unit tests for the core hicDifferentialTAD and hicInterIntraTAD share,
// core/include/hicx/tad_contacts.hpp.
//
// The numeric expectations were produced by scipy 1.14.1 and numpy 1.26.4 in
// the reference environment on the same fixture: the matrix is built as
// csr_matrix((v, (r, c))) + triu(., 1).T, which is fillLowerTriangle, and then
// sliced, summed and passed to scipy.stats.ranksums. The block geometry is
// derived by hand from hicDifferentialTAD.py:116-190 and cooler's
// region_to_extent, with the line of the Python each expectation follows.

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include <doctest/doctest.h>

#include "hicx/tad_contacts.hpp"

namespace {

using hicx::tads::Block;
using hicx::tads::ContactMatrix;
using hicx::tads::Domain;
using hicx::tads::PyNumber;

// Six 10 bp bins on chrA, upper triangle:
//   (0,0)=1 (0,1)=2 (0,3)=3 (1,1)=4 (1,2)=5 (2,4)=6 (3,3)=7 (3,5)=8 (4,5)=9 (5,5)=10
ContactMatrix fixture(hicx::tads::Format format) {
    const std::vector<std::int32_t> rows{0, 0, 0, 1, 1, 2, 3, 3, 4, 5};
    const std::vector<std::int32_t> cols{0, 1, 3, 1, 2, 4, 3, 5, 5, 5};
    std::vector<double> values{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    ContactMatrix matrix;
    matrix.format = format;
    matrix.matrix = hicx::CsrMatrix::from_coo(6, 6, rows, cols, std::move(values), "float64");
    matrix.matrix.set_symmetry(hicx::Symmetry::UpperTriangle);
    std::vector<hicx::CutInterval> bins;
    for (std::int64_t i = 0; i < 6; ++i) {
        bins.push_back(hicx::CutInterval{"chrA", i * 10, i * 10 + 10, 1.0, ""});
    }
    matrix.bins = hicx::BinTable(std::move(bins));
    matrix.chrom_names = {"chrA"};
    matrix.chrom_lengths = {60};
    matrix.bin_size = 10;
    return matrix;
}

Domain domain(const std::string& chrom, std::int64_t start, std::int64_t end) {
    Domain d;
    d.chrom = chrom;
    d.start = start;
    d.end = end;
    d.text = {chrom, std::to_string(start), std::to_string(end), "ID", "0.5", "."};
    return d;
}

// Three TADs: bins [0,2), [2,3), [3,5). The last TAD ends at 50, inside the
// chromosome, which the h5 path needs because it looks up the bin containing
// the end coordinate.
std::vector<Domain> three_tads() {
    return {domain("chrA", 0, 20), domain("chrA", 20, 30), domain("chrA", 30, 50)};
}

}  // namespace

TEST_CASE("normalise_slice follows slice.indices and scipy's inverted slice") {
    using hicx::tads::normalise_slice;
    CHECK(normalise_slice(0, -1, 5) == std::pair<std::int64_t, std::int64_t>{0, 4});
    CHECK(normalise_slice(-7, 10, 5) == std::pair<std::int64_t, std::int64_t>{0, 5});
    CHECK(normalise_slice(3, 2, 5) == std::pair<std::int64_t, std::int64_t>{3, 3});
    CHECK(normalise_slice(2, 9, 5) == std::pair<std::int64_t, std::int64_t>{2, 5});
}

TEST_CASE("extract_block mirrors the upper triangle as fillLowerTriangle does") {
    const ContactMatrix matrix = fixture(hicx::tads::Format::H5);
    // F[1:4, 2:6]
    const hicx::tads::DenseBlock block = hicx::tads::extract_block(matrix, Block{1, 4, 2, 6});
    CHECK(block.values == std::vector<double>{5, 0, 0, 0, 0, 0, 6, 0, 0, 7, 0, 8});
    CHECK(block.nnz == 4);
    const PyNumber sum = hicx::tads::block_sum(block, hicx::DType::Float64);
    CHECK(sum.str() == "26.0");
    CHECK(hicx::tads::block_sum(block, hicx::DType::Integer).str() == "26");

    // The same answer from the explicitly materialised matrix.
    ContactMatrix full = fixture(hicx::tads::Format::H5);
    full.matrix.materialize_full();
    const hicx::tads::DenseBlock same = hicx::tads::extract_block(full, Block{1, 4, 2, 6});
    CHECK(same.values == block.values);
    CHECK(same.nnz == block.nnz);

    // An empty block sums to 0.0, like (M @ ones((0, 1))).sum().
    const hicx::tads::DenseBlock empty = hicx::tads::extract_block(matrix, Block{3, 3, 0, 6});
    CHECK(empty.values.empty());
    CHECK(hicx::tads::block_sum(empty, hicx::DType::Float64).str() == "0.0");
}

TEST_CASE("rank_sum_test matches scipy.stats.ranksums on dense blocks") {
    const ContactMatrix matrix = fixture(hicx::tads::Format::H5);
    const auto values = [&](const Block& block) {
        return hicx::tads::extract_block(matrix, block).values;
    };
    // ranksums(F[0:2, 2:4].flatten(), F[2:4, 4:6].flatten())
    const hicx::tads::RankTest first =
        hicx::tads::rank_sum_test(values(Block{0, 2, 2, 4}), values(Block{2, 4, 4, 6}));
    CHECK(first.statistic == doctest::Approx(-0.5773502691896258).epsilon(1e-15));
    CHECK(first.pvalue == doctest::Approx(0.5637028616507731).epsilon(1e-15));
    // ranksums(F[1:4, 1:4].flatten(), F[2:5, 2:5].flatten())
    const hicx::tads::RankTest second =
        hicx::tads::rank_sum_test(values(Block{1, 4, 1, 4}), values(Block{2, 5, 2, 5}));
    CHECK(second.statistic == doctest::Approx(0.08830215713766959).epsilon(1e-15));
    CHECK(second.pvalue == doctest::Approx(0.9296365245070707).epsilon(1e-15));

    // NaN anywhere propagates, and two empty samples give NaN, as scipy does.
    const std::vector<double> with_nan{1.0, std::nan("")};
    const std::vector<double> plain{2.0, 3.0};
    CHECK(std::isnan(hicx::tads::rank_sum_test(with_nan, plain).pvalue));
    CHECK(std::isnan(hicx::tads::rank_sum_test(with_nan, plain).statistic));
    const std::vector<double> none;
    CHECK(std::isnan(hicx::tads::rank_sum_test(none, none).pvalue));
}

// The three TADs of three_tads() are bins [0, 2), [2, 3) and [3, 5). Every TAD
// but the first has a left block (previous TAD x this TAD) and every TAD but
// the last a right block (this TAD x next TAD), the same for h5 and cool.
TEST_CASE("geometry: neighbouring TADs tile the matrix, the same for h5 and cool") {
    for (const auto format : {hicx::tads::Format::H5, hicx::tads::Format::Cool}) {
        CAPTURE(static_cast<int>(format));
        const ContactMatrix matrix = fixture(format);
        const std::vector<Domain> tads = three_tads();

        const auto g0 = hicx::tads::tad_geometry(matrix, tads, 0);
        CHECK(g0.intra == Block{0, 2, 0, 2, false});
        CHECK_FALSE(g0.left.has_value());
        REQUIRE(g0.right.has_value());
        CHECK(*g0.right == Block{0, 2, 2, 3, false});

        const auto g1 = hicx::tads::tad_geometry(matrix, tads, 1);
        CHECK(g1.intra == Block{2, 3, 2, 3, false});
        REQUIRE(g1.left.has_value());
        CHECK(*g1.left == Block{0, 2, 2, 3, false});
        REQUIRE(g1.right.has_value());
        CHECK(*g1.right == Block{2, 3, 3, 5, false});  // keeps the next TAD's last bin

        // The last TAD: its left block is TAD 1 x TAD 2 (the Python reused the
        // previous iteration's index: an empty block on h5, a wrong one on cool).
        const auto g2 = hicx::tads::tad_geometry(matrix, tads, 2);
        CHECK(g2.intra == Block{3, 5, 3, 5, false});
        CHECK_FALSE(g2.right.has_value());
        REQUIRE(g2.left.has_value());
        CHECK(*g2.left == Block{2, 3, 3, 5, false});

        // A chromosome of two TADs: the second one has its left block (the
        // Python's `i - 1 > 0` skipped it).
        const std::vector<Domain> two{domain("chrA", 0, 20), domain("chrA", 20, 50)};
        const auto t1 = hicx::tads::tad_geometry(matrix, two, 1);
        REQUIRE(t1.left.has_value());
        CHECK(*t1.left == Block{0, 2, 2, 5, false});
    }
}

TEST_CASE("cool geometry: a region past the chromosome end is rejected") {
    const ContactMatrix matrix = fixture(hicx::tads::Format::Cool);
    const std::vector<Domain> beyond{domain("chrA", 0, 20), domain("chrA", 20, 70)};
    CHECK_THROWS((void)hicx::tads::tad_geometry(matrix, beyond, 1));
}

TEST_CASE("cool weights are applied per region, not per matrix") {
    ContactMatrix matrix = fixture(hicx::tads::Format::Cool);
    matrix.weights = std::vector<double>{1.0, 2.0, std::nan(""), 1.0, 1.0, 1.0};
    const std::vector<Domain> single{domain("chrA", 0, 10)};

    // chrA:0-10 holds one stored pixel, (0,0), so len(data) > 1 fails and the
    // raw count stays.
    const auto one = hicx::tads::tad_geometry(matrix, single, 0);
    CHECK_FALSE(one.intra.apply_correction);
    CHECK(hicx::tads::extract_block(matrix, one.intra).values == std::vector<double>{1});
    CHECK(hicx::tads::block_dtype(matrix, one.intra) == hicx::DType::Float64);

    // chrA:0-20 holds (0,0), (0,1), (1,0), (1,1): corrected, count * w_i * w_j.
    const std::vector<Domain> pair{domain("chrA", 0, 20)};
    const auto two = hicx::tads::tad_geometry(matrix, pair, 0);
    CHECK(two.intra.apply_correction);
    CHECK(hicx::tads::extract_block(matrix, two.intra).values ==
          std::vector<double>{1, 4, 4, 16});

    // chrA:10-30 includes the NaN weight of bin 2: (1,2) becomes NaN and then 0,
    // and it no longer counts as stored.
    const std::vector<Domain> nan_region{domain("chrA", 10, 30)};
    const auto with_nan = hicx::tads::tad_geometry(matrix, nan_region, 0);
    const auto block = hicx::tads::extract_block(matrix, with_nan.intra);
    CHECK(block.values == std::vector<double>{16, 0, 0, 0});
    CHECK(block.nnz == 1);
}

TEST_CASE("group_by_chromosome splits on every change of name") {
    const std::vector<Domain> domains{domain("chr1", 0, 10), domain("chr1", 10, 20),
                                      domain("chr2", 0, 10), domain("chr1", 20, 30)};
    const auto groups = hicx::tads::group_by_chromosome(domains);
    REQUIRE(groups.size() == 3);
    CHECK(groups[0].size() == 2);
    CHECK(groups[1].size() == 1);
    CHECK(groups[2].size() == 1);
}

TEST_CASE("PyNumber follows numpy 1.26 scalar promotion and str()") {
    using hicx::tads::py_add;
    using hicx::tads::py_divide;
    // 0 / np.float64(0.0) -> nan, np.int64(3) / np.int64(0) -> inf
    CHECK(py_divide(PyNumber::py_int(0), PyNumber::float64(0.0)).str() == "nan");
    CHECK(py_divide(PyNumber::int64(3), PyNumber::int64(0)).str() == "inf");
    CHECK(py_divide(PyNumber::py_int(0), PyNumber::float64(2.0)).str() == "0.0");
    const PyNumber added = py_add(PyNumber::py_int(0), PyNumber::int64(3));
    CHECK(added.kind == PyNumber::Kind::Int64);
    CHECK(added.str() == "3");
    CHECK(py_add(PyNumber::py_int(0), PyNumber::py_int(0)).str() == "0");
    CHECK(PyNumber::float64(0.1 + 0.2).str() == "0.30000000000000004");
    CHECK_THROWS((void)py_divide(PyNumber::py_int(1), PyNumber::py_int(0)));
}
