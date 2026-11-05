// The core code behind the public Python API: hicx/py_text.hpp,
// hicx/utilities.hpp (hicexplorer.utilities) and hicx/read_bed.hpp
// (hicexplorer.readBed). Their equivalence with HiCExplorer 3.7.6 is tested
// function by function from Python (hicexplorer/test/api); these tests pin
// the C++ contracts, and check the API kernels against the tool kernels that
// already implement the same Python functions.

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <limits>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/build_matrix.hpp"
#include "hicx/obsexp_ops.hpp"
#include "hicx/py_text.hpp"
#include "hicx/python_error.hpp"
#include "hicx/read_bed.hpp"
#include "hicx/tool_matrix.hpp"
#include "hicx/transform_ops.hpp"
#include "hicx/utilities.hpp"

namespace pt = hicx::pytext;
namespace util = hicx::utilities;

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

std::string python_error_type(const std::function<void()>& f) {
    try {
        f();
    } catch (const hicx::PythonError& error) {
        return error.type() + ": " + error.what();
    }
    return "";
}

// chrX of the 50 kb test matrix in its explicit symmetric form
hicx::CsrMatrix real_block() {
    hicx::ToolMatrix hic = hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    hicx::keep_only_chromosomes(hic.data(), {"chrX"});
    hicx::CsrMatrix block = hic.matrix();
    block.materialize_full();
    return block;
}

util::SparseEntries entries_of(const hicx::CsrMatrix& m) {
    util::SparseEntries e;
    e.rows = m.rows();
    e.cols = m.cols();
    for (std::int64_t r = 0; r < m.rows(); ++r) {
        for (std::int64_t k = m.indptr()[static_cast<std::size_t>(r)];
             k < m.indptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            e.row.push_back(r);
            e.col.push_back(m.indices()[static_cast<std::size_t>(k)]);
        }
    }
    return e;
}

}  // namespace

TEST_CASE("pytext follows str.isspace, strip and split") {
    CHECK(pt::whitespace_length("\x1c", 0) == 1);
    CHECK(pt::whitespace_length("　", 0) == 3);
    CHECK(pt::whitespace_length("\u0085", 0) == 2);
    CHECK(pt::whitespace_length("​", 0) == 0);  // zero width space is not whitespace
    CHECK(pt::whitespace_length("᠎", 0) == 0);
    CHECK(pt::strip("　 a b\t\n ") == "a b");
    CHECK(pt::strip(" \t ").empty());
    const std::vector<std::string_view> ws = pt::split_whitespace("  a b  c\n");
    REQUIRE(ws.size() == 3);
    CHECK(ws[1] == "b");
    const std::vector<std::string_view> tabs = pt::split("a\t\tb\t", "\t");
    REQUIRE(tabs.size() == 4);
    CHECK(tabs[1].empty());
    CHECK(tabs[3].empty());
    CHECK(pt::remove_whitespace(" c h　r ") == "chr");
    CHECK(pt::code_point_count("aé北") == 3);
}

TEST_CASE("pytext reads text files as open(path, 'r') does") {
    const std::vector<std::string> lines = pt::universal_newline_lines("a\r\nb\rc\n\nd");
    REQUIRE(lines.size() == 5);
    CHECK(lines[0] == "a\n");
    CHECK(lines[1] == "b\n");
    CHECK(lines[3] == "\n");
    CHECK(lines[4] == "d");
    CHECK(pt::universal_newline_lines("").empty());
    CHECK(!pt::first_invalid_utf8("ok é € \U0001f600"));
    CHECK(pt::first_invalid_utf8("ab\xff") == 2u);
    CHECK(pt::first_invalid_utf8("a\xc0\x80") == 1u);          // overlong
    CHECK(pt::first_invalid_utf8("\xed\xa0\x80") == 0u);       // surrogate
    CHECK(pt::first_invalid_utf8("x\xe2\x82") == 1u);          // cut off
    CHECK(pt::first_invalid_utf8("\xf4\x90\x80\x80") == 0u);   // beyond U+10FFFF
}

TEST_CASE("Python arithmetic and numpy casts") {
    CHECK(util::python_mod(-7, 3) == 2);
    CHECK(util::python_mod(7, -3) == -2);
    CHECK(util::python_mod(10, 3) == 1);
    CHECK(util::python_float_mod(-7.5, 2.0) == doctest::Approx(0.5));
    CHECK(std::signbit(util::python_float_mod(0.0, -2.0)));
    CHECK(util::python_greater(std::int64_t{9007199254740993}, 9007199254740992.0));
    CHECK(!util::python_greater(1, std::nan("")));
    CHECK(util::cast_from_double<std::int64_t>(std::nan("")) == std::numeric_limits<std::int64_t>::min());
    CHECK(util::cast_from_double<std::int32_t>(3e9) == std::numeric_limits<std::int32_t>::min());
    CHECK(util::cast_from_double<std::int32_t>(-2.9) == -2);
}

TEST_CASE("readBed and readTargetBed") {
    const auto viewpoints = util::read_bed({"#c\n", "chr1\t1\t2\textra\n", "chr2 3 4\n"});
    REQUIRE(viewpoints.size() == 2);
    CHECK(viewpoints[1][0] == "chr2");
    CHECK(viewpoints[1][2] == "4");
    CHECK(python_error_type([] { (void)util::read_bed({"\n"}); }) ==
          "ValueError: not enough values to unpack (expected 3, got 0)");

    const util::TargetBed target =
        util::read_target_bed({"chr1\t1\t2\tg1\n", "chr1\t3\t4\n", "chr2 5 6 g2\n", "chr3\t7\t8\tg1\n"});
    REQUIRE(target.genes.size() == 2);
    CHECK(target.genes[0] == "g1");
    CHECK(target.positions[0].chromosome == "chr3");
    CHECK(target.positions[0].starts == std::vector<std::string>{"1", "3", "7"});
    CHECK(python_error_type([] { (void)util::read_target_bed({"chr1\t1\t2\n"}); }).rfind("UnboundLocalError", 0) == 0);
    CHECK(python_error_type([] { (void)util::read_target_bed({"chr1 1\n"}); }) ==
          "ValueError: not enough values to unpack (expected 4, got 2)");
}

TEST_CASE("genomicRegion and getUserRegion") {
    CHECK(util::genomic_region(" chr1 : 1,000 - 2,000 ") == std::optional<std::string>("chr1:1000:2000"));
    CHECK(!util::genomic_region(" \t　"));
    CHECK(python_error_type([] { (void)util::genomic_region(",;"); }) ==
          "argparse.ArgumentTypeError: ,; is not a valid region");
    // agrees with the tools' normalise_region on ASCII input
    CHECK(*util::genomic_region("chr1:{1,0}-(5)") == hicx::normalise_region("chr1:{1,0}-(5)"));

    // the two doctests of getUserRegion
    util::UserRegion a = util::user_region_bounds(1000, 10, 1001, "chr2:10:1001");
    util::user_region_finish(a, std::nullopt, util::ChunkLimit{});
    CHECK(a.start == 10);
    CHECK(a.end == 1000);
    CHECK(a.end_is_chrom_size);
    CHECK(a.chunk_size == 990);
    util::UserRegion b = util::user_region_bounds(200000, 10, 123344, "chr2:10:123344:3");
    util::user_region_finish(b, 3, util::ChunkLimit{});
    CHECK(b.start == 9);
    CHECK(b.end == 123345);
    CHECK(b.chunk_size == 123336);
    CHECK(python_error_type([] { (void)util::user_region_bounds(100, 50, 10, "c:50:10"); }).rfind("NameError", 0) == 0);
}

TEST_CASE("enlarge_bins agrees with the tools' implementation") {
    hicx::ToolMatrix hic = hicx::ToolMatrix::load(data_path("small_test_matrix_50kb_res.h5"));
    std::vector<hicx::CutInterval> intervals;
    for (std::size_t i = 0; i < hic.data().cut_intervals.size(); ++i) {
        if (i % 5 != 2) {
            intervals.push_back(hic.data().cut_intervals[i]);
        }
    }
    std::vector<std::int64_t> starts;
    std::vector<std::int64_t> ends;
    std::vector<bool> same;
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        starts.push_back(intervals[i].start);
        ends.push_back(intervals[i].end);
        if (i + 1 < intervals.size()) {
            same.push_back(intervals[i].chrom == intervals[i + 1].chrom);
        }
    }
    const std::vector<util::EnlargedBin> bins = util::enlarge_bins(starts, ends, same);
    hicx::enlarge_bins(intervals);
    bool agree = true;
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        agree = agree && bins[i].start == intervals[i].start && bins[i].end == intervals[i].end;
    }
    CHECK(agree);
    // the doctest
    const std::vector<util::EnlargedBin> doc = util::enlarge_bins({10, 50, 10, 70}, {50, 80, 60, 90},
                                                                  {true, false, true});
    CHECK(doc[0].start == 0);
    CHECK(doc[2].start == 0);
    CHECK(doc[2].end == 65);
    CHECK(doc[3].start == 65);
    CHECK(doc[3].start_kind == util::Coordinate::Middle);
}

TEST_CASE("statistics follow numpy") {
    CHECK(util::median(std::vector<double>{3, 1, 2, 10}) == 2.5);
    CHECK(std::isnan(util::median(std::vector<double>{1, std::nan("")})));
    CHECK(std::isnan(util::median(std::vector<double>{})));
    CHECK(util::median(std::vector<std::int64_t>{5, 1, 3}) == 3.0);
    CHECK(util::mean(std::vector<float>{1.0F, 2.0F}) == 1.5);
    const auto keep = util::remove_outliers_keep(std::vector<double>{1, 2, 3, 2, 1, 100}, 3.5);
    REQUIRE(keep);
    CHECK(keep->size() == 5);
    CHECK(!util::remove_outliers_keep(std::vector<double>{4, 4, 4}, 3.5));
}

TEST_CASE("the API's expected and obs/exp kernels agree with the tools' on a real block") {
    const hicx::CsrMatrix block = real_block();
    const util::SparseEntries entries = entries_of(block);
    const std::vector<double>& data = block.data();
    const util::Nonzero nz = util::nonzero(entries, data.data());
    const std::vector<std::int64_t> distance = nz.distances();

    const std::optional<std::vector<double>> expected =
        util::expected_interactions(block.rows(), distance, data.data(), data.size());
    REQUIRE(expected);
    CHECK(*expected == hicx::expected_interactions(block));
    CHECK(util::expected_interactions_non_zero(block.rows(), distance, data.data()) ==
          hicx::expected_interactions_non_zero(block));
    const std::vector<double> lieberman = util::expected_interactions_in_distance(
        block.rows(), distance, data.data(), static_cast<std::int32_t>(block.rows()), 1);
    const std::vector<double> tools_lieberman =
        hicx::expected_interactions_in_distance(block, block.rows(), 1);
    bool same = lieberman.size() == tools_lieberman.size();
    for (std::size_t i = 0; same && i < lieberman.size(); ++i) {
        same = (lieberman[i] == tools_lieberman[i]) || (std::isnan(lieberman[i]) && std::isnan(tools_lieberman[i]));
    }
    CHECK(same);

    hicx::CsrMatrix tools = block;
    hicx::obs_exp_in_place(tools);
    // the test matrix is int64, so the quotient is truncated back to integers
    REQUIRE(block.dtype() == "int64");
    const std::vector<std::int64_t> counts(data.begin(), data.end());
    const std::vector<std::int64_t> divided =
        util::obs_exp_divide(counts.data(), counts.size(), *expected, distance, 0.0);
    CHECK(std::vector<double>(divided.begin(), divided.end()) == tools.data());
    const std::vector<double> as_float =
        util::obs_exp_divide(data.data(), data.size(), *expected, distance, 0.0);
    CHECK(std::count_if(as_float.begin(), as_float.end(), [](double v) { return v != std::floor(v); }) > 0);

    // the threaded Python splits the distances into ranges; the sums are the same
    const auto [lo, hi] = std::minmax_element(distance.begin(), distance.end());
    std::vector<double> parts(static_cast<std::size_t>(block.rows()), 0.0);
    const std::int64_t mid = (*lo + *hi) / 2;
    for (const auto& range : {std::pair{*lo, mid}, std::pair{mid, *hi + 1}}) {
        const std::vector<double> part =
            util::distance_sums(block.rows(), distance, data.data(), data.size(), range.first, range.second);
        for (std::size_t k = 0; k < part.size(); ++k) {
            parts[k] += part[k];
        }
    }
    CHECK(parts == util::distance_sums(block.rows(), distance, data.data(), data.size(), *lo, *hi + 1));
}

TEST_CASE("explicit zeros reproduce the Python's misalignment and errors") {
    // csr [[1, 0*, 2], [3, 0*], ...]: nonzero() drops the stored zeros
    util::SparseEntries e;
    e.rows = 3;
    e.cols = 3;
    e.row = {0, 0, 0, 1, 1};
    e.col = {0, 1, 2, 1, 2};
    const std::vector<double> data = {1.0, 0.0, 2.0, 3.0, 0.0};
    const util::Nonzero nz = util::nonzero(e, data.data());
    REQUIRE(nz.size() == 3);
    const std::vector<double> mean = util::expected_interactions_non_zero(3, nz.distances(), data.data());
    // distances 0, 2, 0 take data[0], data[1], data[2] = 1, 0, 2
    CHECK(mean[0] == 1.5);
    CHECK(mean[2] == 0.0);
    CHECK(python_error_type([&] { (void)util::expected_interactions(3, nz.distances(), data.data(), data.size()); }) ==
          "IndexError: boolean index did not match indexed array along dimension 0; dimension is 5 but "
          "corresponding boolean dimension is 3");
    std::vector<std::int64_t> indptr = {0, 3, 5, 5};
    const std::vector<std::size_t> keep = util::eliminate_zeros(indptr, data.data());
    CHECK(keep == std::vector<std::size_t>{0, 2, 3});
    CHECK(indptr == std::vector<std::int64_t>{0, 2, 3, 3});
}

TEST_CASE("text helpers of utilities") {
    CHECK(util::change_chrom_name("chr1") == "1");
    CHECK(util::change_chrom_name("X") == "chrX");
    CHECK(util::in_units(2e6) == "2.00  Mbp");
    CHECK(util::in_units(1501) == "2 Kbp");
    CHECK(util::in_units(1500) == "1500.00  bp");
    CHECK(util::in_units(std::nan("")) == "nan  bp");
    CHECK(util::unidecode("Knösós") == "Knosos");
    CHECK(util::unidecode("北亰") == "Bei Jing ");
    std::vector<char32_t> surrogates;
    CHECK(util::unidecode("a\xed\xb2\x80" "b", &surrogates) == "ab");
    CHECK(surrogates.size() == 1);
    CHECK(util::gzip_magic("\x1f\x8b"));
    CHECK(!util::gzip_magic("\x1f"));
}

TEST_CASE("ReadBed rules") {
    std::string assumed;
    CHECK(hicx::readbed::guess_file_type(4) == "bedgraph");
    CHECK(hicx::readbed::guess_file_type(7, &assumed) == "bed6");
    CHECK(assumed == "bed6");
    CHECK(hicx::readbed::guess_file_type(5, &assumed) == "bed3");
    CHECK(hicx::readbed::interval_fields("bed9").size() == 9);
    CHECK(hicx::readbed::interval_fields("bedgraph").size() == 6);
    CHECK(!hicx::readbed::required_field_count("bed12"));
    bool invalid = false;
    CHECK(hicx::readbed::normalize_strand("-1", invalid) == "-");
    CHECK(!invalid);
    CHECK(hicx::readbed::normalize_strand("x", invalid) == ".");
    CHECK(invalid);
    CHECK(hicx::readbed::is_comment_line("browser x"));
    CHECK(hicx::readbed::is_comment_line(" 　\n"));
    CHECK(!hicx::readbed::is_comment_line(" chr1"));
    CHECK(hicx::readbed::split_commas("10,20,", true) == std::vector<std::string>{"10", "20"});
    CHECK(hicx::readbed::split_line("a\tb\n", false).back() == "b\n");
}
