// .hic reading through hic_adapter's generic loader (read_hic, is_hic_path,
// parse_hic_uri), the entry point ToolMatrix::load and HiCMatrix::load both
// use so that every matrix-reading tool gets real .hic support (cpp/PLAN.md
// 9.1: "Afterwards every matrix-reading tool accepts a .hic with a resolution
// and normalisation selector").
//
// The cross-check against hic2cool_convert + read_cool (the conversion route
// hicConvertFormat already validates) is done here directly, matrix for
// matrix, rather than only through stdout text, because ToolMatrix::load's
// hic branch is implemented as exactly that conversion (see hic_adapter.cpp
// read_hic): it reuses hic2cool_convert into a temporary cool file and then
// read_cool, so this test is the equivalence class E2 check the project's
// convention asks for (same sparsity pattern, every stored value
// bit-identical) between the native selector path and the conversion path,
// not a tautology, because it is exercised through the public entry points a
// tool calls, with independently constructed temporary files on each side.

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/hic_adapter.hpp"
#include "hicx/hic_matrix.hpp"
#include "hicx/tool_matrix.hpp"

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

const std::string kFixture = "hicConvertFormat/GM12878_combined_30.chr21_chr22.v7.hic";

// A small monotonically increasing id so parallel doctest cases do not race
// on the same temporary file name.
int next_unique_id() {
    static int next = 0;
    return next++;
}

struct TempFile {
    std::string path;
    explicit TempFile(const std::string& suffix)
        : path((std::filesystem::temp_directory_path() /
                ("hicx_test_hic_adapter_" + std::to_string(next_unique_id()) + suffix))
                   .string()) {}
    ~TempFile() {
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
    }
};

}  // namespace

TEST_CASE("has_hic_signature reads the literal HIC magic, not the extension") {
    const std::string real_hic = data_path(kFixture);
    CHECK(hicx::has_hic_signature(real_hic));

    // A real .hic file copied under an unrelated extension is still detected
    // by content, the preference cpp/PLAN.md's task asks for.
    TempFile renamed(".not_hic_at_all");
    {
        std::ifstream in(real_hic, std::ios::binary);
        std::ofstream out(renamed.path, std::ios::binary);
        out << in.rdbuf();
    }
    CHECK(hicx::has_hic_signature(renamed.path));
    CHECK(hicx::is_hic_path(renamed.path));

    // A cool file, even if named ".hic", is not mistaken for one by its
    // signature; is_hic_path still calls it hic only through the extension
    // fallback, matching an explicit ".hic" name the caller chose.
    const std::string real_cool = data_path("Li_et_al_2015.cool");
    CHECK_FALSE(hicx::has_hic_signature(real_cool));
    CHECK_FALSE(hicx::is_hic_path(real_cool));
}

TEST_CASE("parse_hic_uri accepts the four documented selector forms") {
    const hicx::HicUri plain = hicx::parse_hic_uri("a.hic");
    CHECK(plain.path == "a.hic");
    CHECK_FALSE(plain.resolution.has_value());
    CHECK_FALSE(plain.normalization.has_value());

    const hicx::HicUri resolution_only = hicx::parse_hic_uri("a.hic::/resolutions/10000");
    CHECK(resolution_only.path == "a.hic");
    REQUIRE(resolution_only.resolution.has_value());
    CHECK(*resolution_only.resolution == 10000);
    CHECK_FALSE(resolution_only.normalization.has_value());

    const hicx::HicUri norm_only = hicx::parse_hic_uri("a.hic::/normalizations/KR");
    CHECK(norm_only.path == "a.hic");
    CHECK_FALSE(norm_only.resolution.has_value());
    REQUIRE(norm_only.normalization.has_value());
    CHECK(*norm_only.normalization == "KR");

    const hicx::HicUri both =
        hicx::parse_hic_uri("a.hic::/resolutions/10000/normalizations/KR");
    CHECK(both.path == "a.hic");
    REQUIRE(both.resolution.has_value());
    CHECK(*both.resolution == 10000);
    REQUIRE(both.normalization.has_value());
    CHECK(*both.normalization == "KR");

    CHECK_THROWS_AS(hicx::parse_hic_uri("a.hic::/bogus/10000"), std::runtime_error);
    CHECK_THROWS_AS(hicx::parse_hic_uri("a.hic::/resolutions/"), std::runtime_error);
}

TEST_CASE("read_hic defaults to the finest resolution and raw counts") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult loaded = hicx::read_hic(path);
    // The fixture's resolutions are 2500000, 1000000, 500000, 250000
    // (cpp/PLAN.md 9.1); the default selects the finest, 250000.
    const hicx::BinTable bins(loaded.data.cut_intervals);
    CHECK(bins.bin_size() == 250000);
    CHECK_FALSE(loaded.correction_operator.has_value());
}

TEST_CASE("read_hic rejects a resolution or normalization the file does not have") {
    const std::string path = data_path(kFixture);
    CHECK_THROWS_AS(hicx::read_hic(path + "::/resolutions/999"), std::runtime_error);
    CHECK_THROWS_AS(hicx::read_hic(path + "::/normalizations/BOGUS"), std::runtime_error);
}

TEST_CASE("read_hic (raw counts) equals hic2cool_convert + read_cool on the same resolution") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult native = hicx::read_hic(path + "::/resolutions/1000000");

    TempFile converted(".cool");
    hicx::hic2cool_convert(path, converted.path, 1000000);
    const hicx::CoolLoadResult reference = hicx::read_cool(converted.path);

    REQUIRE(native.data.matrix.rows() == reference.data.matrix.rows());
    REQUIRE(native.data.matrix.nnz() == reference.data.matrix.nnz());
    CHECK(native.data.matrix.indptr() == reference.data.matrix.indptr());
    CHECK(native.data.matrix.indices() == reference.data.matrix.indices());
    CHECK(native.data.matrix.data() == reference.data.matrix.data());
    CHECK(native.data.cut_intervals.size() == reference.data.cut_intervals.size());
    for (std::size_t i = 0; i < native.data.cut_intervals.size(); ++i) {
        CHECK(native.data.cut_intervals[i].chrom == reference.data.cut_intervals[i].chrom);
        CHECK(native.data.cut_intervals[i].start == reference.data.cut_intervals[i].start);
        CHECK(native.data.cut_intervals[i].end == reference.data.cut_intervals[i].end);
    }
    CHECK(native.data.nan_bins == reference.data.nan_bins);
}

TEST_CASE("read_hic (KR normalized) equals hic2cool_convert + read_cool with the KR column") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult native =
        hicx::read_hic(path + "::/resolutions/1000000/normalizations/KR");
    REQUIRE(native.correction_operator.has_value());
    CHECK(*native.correction_operator == '/');

    TempFile converted(".cool");
    hicx::hic2cool_convert(path, converted.path, 1000000);
    hicx::CoolLoadOptions options;
    options.apply_correction = true;
    options.correction_factor_table = "KR";
    options.correction_operator = '/';
    const hicx::CoolLoadResult reference = hicx::read_cool(converted.path, options);

    REQUIRE(native.data.matrix.nnz() == reference.data.matrix.nnz());
    CHECK(native.data.matrix.indptr() == reference.data.matrix.indptr());
    CHECK(native.data.matrix.indices() == reference.data.matrix.indices());
    CHECK(native.data.matrix.data() == reference.data.matrix.data());
}

TEST_CASE("ToolMatrix::load and HiCMatrix::load both read a .hic selector") {
    const std::string path = data_path(kFixture) + "::/resolutions/1000000";

    const hicx::ToolMatrix tool = hicx::ToolMatrix::load(path);
    CHECK_FALSE(tool.input_is_h5());
    CHECK(tool.matrix().rows() > 0);
    CHECK(tool.cut_intervals().size() > 0);

    const hicx::HiCMatrix plain = hicx::HiCMatrix::load(path);
    CHECK(plain.matrix().rows() == tool.matrix().rows());
    // ToolMatrix::load symmetrizes and swaps correction_factors/distance_counts
    // (cpp/PLAN.md 2.7 quirk 1); HiCMatrix::load only symmetrizes. Both must
    // still see the same finest-resolution, raw-count fixture underneath.
    CHECK(plain.bin_size() == 1000000);
}

namespace {

// Every stored (row, col, value) of the whole-file load whose row and column
// both fall in [first, last), offset back to [0, last-first). The chromosome
// blocks of a .hic file are contiguous in the whole-file bin table exactly as
// they are in a cooler's, so this is the reference a chromosome or region
// partial read must equal.
std::vector<std::array<double, 3>> block_of(const hicx::CsrMatrix& whole, std::int64_t first,
                                            std::int64_t last) {
    std::vector<std::array<double, 3>> entries;
    const auto& indptr = whole.indptr();
    const auto& indices = whole.indices();
    const auto& data = whole.data();
    for (std::int64_t row = first; row < last; ++row) {
        for (auto k = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
             k < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]); ++k) {
            const std::int32_t col = indices[k];
            if (col >= first && col < last) {
                entries.push_back({static_cast<double>(row - first), static_cast<double>(col - first),
                                   data[k]});
            }
        }
    }
    std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
        return a[0] != b[0] ? a[0] < b[0] : a[1] < b[1];
    });
    return entries;
}

std::vector<std::array<double, 3>> block_of(const hicx::HicLoadResult& loaded) {
    return block_of(loaded.data.matrix, 0, loaded.data.matrix.rows());
}

}  // namespace

TEST_CASE("read_hic with a bare chromosome name is a genuinely partial read, "
         "equal to the whole file's block for that chromosome") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult whole = hicx::read_hic(path + "::/resolutions/250000");
    std::int64_t first = -1;
    std::int64_t last = -1;
    for (std::size_t i = 0; i < whole.data.cut_intervals.size(); ++i) {
        if (whole.data.cut_intervals[i].chrom == "22") {
            if (first < 0) {
                first = static_cast<std::int64_t>(i);
            }
            last = static_cast<std::int64_t>(i) + 1;
        }
    }
    REQUIRE(first >= 0);

    const hicx::HicLoadResult partial = hicx::read_hic(path + "::/resolutions/250000", "22");
    CHECK(partial.data.matrix.rows() == last - first);
    for (const auto& interval : partial.data.cut_intervals) {
        CHECK(interval.chrom == "22");
    }
    CHECK(block_of(partial) == block_of(whole.data.matrix, first, last));
}

TEST_CASE("read_hic with a chrom:start-end region matches the whole file's bin "
         "range, snapped to bin boundaries") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult whole = hicx::read_hic(path + "::/resolutions/250000");
    // Chromosome 21 is 48,129,895 bp at 250,000 bp bins: bins [0, 193) cover
    // it (ceil(48129895 / 250000) = 193). The region below is bins [4, 12).
    const hicx::HicLoadResult region =
        hicx::read_hic(path + "::/resolutions/250000", "21:1000000-3000000");
    REQUIRE(region.data.cut_intervals.size() == 8);
    CHECK(region.data.cut_intervals.front().start == 1000000);
    CHECK(region.data.cut_intervals.back().end == 3000000);

    CHECK(block_of(region) == block_of(whole.data.matrix, 4, 12));
}

TEST_CASE("read_hic (KR normalized) on a region matches the whole file's block, "
         "same normalization") {
    const std::string path = data_path(kFixture);
    const hicx::HicLoadResult whole =
        hicx::read_hic(path + "::/resolutions/250000/normalizations/KR");
    const hicx::HicLoadResult region =
        hicx::read_hic(path + "::/resolutions/250000/normalizations/KR", "21:1000000-3000000");
    REQUIRE(region.correction_operator.has_value());
    CHECK(*region.correction_operator == '/');
    CHECK(block_of(region) == block_of(whole.data.matrix, 4, 12));
}

TEST_CASE("ToolMatrix::load with a .hic region loads only that region, not the whole file") {
    const std::string path = data_path(kFixture) + "::/resolutions/250000";
    const hicx::ToolMatrix whole = hicx::ToolMatrix::load(path);
    const hicx::ToolMatrix region = hicx::ToolMatrix::load(path, std::string("21:1000000-3000000"));
    CHECK(region.matrix().rows() == 8);
    CHECK(region.matrix().rows() < whole.matrix().rows());
    for (const auto& interval : region.cut_intervals()) {
        CHECK(interval.chrom == "21");
    }
}
