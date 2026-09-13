// read_cool with a chromosome name cuts the chromosome's block out of the
// pixel table in place (core/src/cool_file.cpp). Before that change the cut
// was select_bins over the chromosome's contiguous bin range, and every tool
// that loads one chromosome of a cooler (hicAdjustMatrix, hicCorrectMatrix,
// hicPlotSVL) was validated against that. This pins that the in-place cut
// produces exactly the matrix select_bins produces, on real coolers: every
// chromosome of a many-contig file, and a file with a weight column.

#include <doctest/doctest.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "hicx/adjust_ops.hpp"
#include "hicx/cool_file.hpp"

namespace {

const std::string kData = HICX_TEST_DATA_DIR;

void check_every_chromosome(const std::string& path) {
    const hicx::CoolFile cool(path);
    hicx::CoolLoadOptions whole_options;
    const hicx::CoolLoadResult whole = hicx::read_cool(path, whole_options);
    for (const std::string& chromosome : cool.chrom_names()) {
        CAPTURE(path);
        CAPTURE(chromosome);
        std::int64_t first = -1;
        std::int64_t last = -1;
        for (std::size_t bin = 0; bin < whole.data.cut_intervals.size(); ++bin) {
            if (whole.data.cut_intervals[bin].chrom == chromosome) {
                if (first < 0) {
                    first = static_cast<std::int64_t>(bin);
                }
                last = static_cast<std::int64_t>(bin) + 1;
            }
        }
        std::vector<std::int64_t> selection;
        for (std::int64_t bin = first; bin < last; ++bin) {
            selection.push_back(bin);
        }
        const hicx::CsrMatrix expected = hicx::select_bins(whole.data.matrix, selection);

        hicx::CoolLoadOptions options;
        options.chrom_name = chromosome;
        const hicx::CoolLoadResult cut = hicx::read_cool(path, options);
        const hicx::CsrMatrix& got = cut.data.matrix;

        CHECK(got.rows() == expected.rows());
        CHECK(got.cols() == expected.cols());
        CHECK(got.dtype() == expected.dtype());
        CHECK(got.symmetry() == expected.symmetry());
        CHECK(got.indptr() == expected.indptr());
        CHECK(got.indices() == expected.indices());
        REQUIRE(got.data().size() == expected.data().size());
        CHECK(std::memcmp(got.data().data(), expected.data().data(),
                          got.data().size() * sizeof(double)) == 0);
        CHECK(cut.data.cut_intervals.size() == selection.size());
    }
}

}  // namespace

TEST_CASE("the in-place chromosome cut of read_cool equals select_bins, 15 contigs") {
    check_every_chromosome(kData + "/small_test_matrix_50kb_res.cool");
}

TEST_CASE("the in-place chromosome cut of read_cool equals select_bins, weight column") {
    check_every_chromosome(kData + "/hicCorrectMatrix/kr_full.cool");
}

TEST_CASE("the in-place chromosome cut of read_cool equals select_bins, two chromosomes") {
    check_every_chromosome(kData +
                           "/hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1_chr2.cool");
}
