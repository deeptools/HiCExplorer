// Unit tests of the .pairs reader behind hicBuildMatrix --pairsFile
// (cpp/PLAN.md tier 9, item 9.2). The whole-file tests read the committed
// inputs of hicexplorer/test/test_data/hicBuildMatrix/pairs, which
// cpp/scripts/make_bam_route_pairs.py wrote from the test BAM files.

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "hicx/pairs_file.hpp"

using hicx::PairKind;
using hicx::PairsHeader;
using hicx::PairsReader;
using hicx::PairsRecord;

namespace {

const std::string kPairsDir = std::string(HICX_TEST_DATA_DIR) + "/hicBuildMatrix/pairs/";

PairsRecord parse_line(const std::string& line, const PairsHeader& header) {
    std::vector<std::string_view> fields;
    PairsRecord record;
    const auto error = hicx::parse_pairs_record(line, header, fields, record);
    if (error.has_value()) {
        throw std::runtime_error(*error);
    }
    return record;
}

PairKind kind_of(const std::string& type, const PairsHeader& header) {
    return parse_line(".\tchr1\t10\tchr2\t20\t+\t-\t" + type + "\t30\t40", header).kind;
}

std::int64_t count_lines(const std::string& path, int threads) {
    PairsReader reader(path, threads);
    std::string_view line;
    std::int64_t count = 0;
    while (reader.next(line)) {
        ++count;
    }
    return count;
}

}  // namespace

TEST_CASE("pairs header: the reserved columns without a #columns line") {
    const PairsHeader header = hicx::parse_pairs_header(
        {"## pairs format v1.0", "#sorted: none", "#chromsize: chrX 1000"});
    CHECK(header.format_version == "1.0");
    CHECK_FALSE(header.columns_declared);
    CHECK(header.chrom1 == 1);
    CHECK(header.pos1 == 2);
    CHECK(header.chrom2 == 3);
    CHECK(header.pos2 == 4);
    CHECK(header.strand1 == 5);
    CHECK(header.strand2 == 6);
    CHECK(header.mapq1 == -1);
    CHECK(header.pair_type == -1);
    CHECK(header.fields_needed == 7);
    CHECK(header.sorted_keys == std::vector<std::string>{"none"});
    CHECK_FALSE(header.sorted_by_chroms_then_positions());
    REQUIRE(header.chrom_sizes.size() == 1);
    CHECK(header.chrom_sizes[0].first == "chrX");
    CHECK(header.chrom_sizes[0].second == 1000);
}

TEST_CASE("pairs header: chr1 and chr2 are chrom1 and chrom2, optional columns found") {
    const PairsHeader header = hicx::parse_pairs_header(
        {"## pairs format v1.0.0", "#sorted: chr1-chr2-pos1-pos2", "#shape: upper triangle",
         "#genome_assembly: mm10", "#samheader: @SQ\tSN:chr1\tLN:10",
         "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 frag1 frag2 pair_type mapq1 mapq2"});
    CHECK(header.sorted_by_chroms_then_positions());
    CHECK(header.upper_triangle());
    CHECK(header.genome_assembly == "mm10");
    CHECK(header.chrom1 == 1);
    CHECK(header.chrom2 == 3);
    CHECK(header.pair_type == 9);
    CHECK(header.mapq1 == 10);
    CHECK(header.mapq2 == 11);
    CHECK(header.fields_needed == 12);
    CHECK(header.chrom_sizes.empty());
}

TEST_CASE("pairs header: refused headers") {
    CHECK_THROWS_WITH_AS(hicx::parse_pairs_header({"#columns: readID chrom1 pos1 chrom2 pos2"}),
                         doctest::Contains("## pairs format v1.0"), std::runtime_error);
    CHECK_THROWS_AS(hicx::parse_pairs_header({}), std::runtime_error);
    CHECK_THROWS_WITH_AS(
        hicx::parse_pairs_header({"## pairs format v1.0", "#chromsize: chr1 10", "#chromsize: chr1 12"}),
        doctest::Contains("more than one"), std::runtime_error);
    CHECK_THROWS_WITH_AS(hicx::parse_pairs_header({"## pairs format v1.0", "#chromsize: chr1"}),
                         doctest::Contains("malformed #chromsize"), std::runtime_error);
    CHECK_THROWS_WITH_AS(
        hicx::parse_pairs_header({"## pairs format v1.0", "#columns: readID chrom1 pos1 chrom2"}),
        doctest::Contains("pos2"), std::runtime_error);
    CHECK_THROWS_WITH_AS(
        hicx::parse_pairs_header({"## pairs format v1.0", "#columns: chrom1 pos1 chr1 pos2 chrom2"}),
        doctest::Contains("twice"), std::runtime_error);
}

TEST_CASE("pairs record: fields, strands, mapping qualities and pair types") {
    const PairsHeader header = hicx::parse_pairs_header(
        {"## pairs format v1.0",
         "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2"});
    // The record's chromosome names are views into the line, which must outlive them.
    const std::string text = ".\tchr2L\t8815\tchr3R\t53053\t-\t+\tUU\t44\t7";
    const PairsRecord record = parse_line(text, header);
    CHECK(record.chrom1 == "chr2L");
    CHECK(record.chrom2 == "chr3R");
    CHECK(record.pos1 == 8815);
    CHECK(record.pos2 == 53053);
    CHECK(record.strand1 == -1);
    CHECK(record.strand2 == 1);
    CHECK(record.mapq1 == 44);
    CHECK(record.mapq2 == 7);
    CHECK(record.kind == PairKind::Mapped);

    CHECK(kind_of("UR", header) == PairKind::Mapped);
    CHECK(kind_of("RU", header) == PairKind::Mapped);
    CHECK(kind_of(".", header) == PairKind::Mapped);
    CHECK(kind_of("MU", header) == PairKind::NotUnique);
    CHECK(kind_of("MM", header) == PairKind::NotUnique);
    CHECK(kind_of("NU", header) == PairKind::Unmapped);
    CHECK(kind_of("MN", header) == PairKind::Unmapped);
    CHECK(kind_of("WW", header) == PairKind::Unmapped);
    CHECK(kind_of("XX", header) == PairKind::Unmapped);
    CHECK(kind_of("DD", header) == PairKind::MarkedDuplicate);

    // '!' is unmapped when the pair type does not say more.
    CHECK(parse_line(".\t!\t0\tchr1\t20\t-\t+\tUU\t0\t40", header).kind == PairKind::Unmapped);
    CHECK(parse_line(".\t!\t0\tchr1\t20\t-\t+\tMU\t0\t40", header).kind == PairKind::NotUnique);
    // A strand other than + and - is unknown.
    CHECK(parse_line(".\tchr1\t1\tchr1\t2\t.\t+\tUU\t1\t1", header).strand1 == 0);

    std::vector<std::string_view> fields;
    PairsRecord scratch;
    CHECK(hicx::parse_pairs_record(".\tchr1\t1\tchr1\t2\t+\t+\tUU\t1", header, fields, scratch)
              ->find("needs at least 10") != std::string::npos);
    CHECK(hicx::parse_pairs_record(".\tchr1\tx1\tchr1\t2\t+\t+\tUU\t1\t1", header, fields, scratch)
              ->find("pos1 'x1'") != std::string::npos);
    CHECK(hicx::parse_pairs_record(".\tchr1\t-5\tchr1\t2\t+\t+\tUU\t1\t1", header, fields, scratch)
              .has_value());
    CHECK(hicx::parse_pairs_record(".\tchr1\t1\tchr1\t2\t+\t+\tUU\t1\t60x", header, fields, scratch)
              ->find("mapq2") != std::string::npos);
}

TEST_CASE("pairs record: the reserved layout needs seven fields only") {
    const PairsHeader header = hicx::parse_pairs_header({"## pairs format v1.0"});
    const PairsRecord record = parse_line(".\tchr1\t5\tchr1\t9\t+\t-\textra\tcolumns", header);
    CHECK(record.pos2 == 9);
    CHECK(record.mapq1 == -1);
    CHECK(record.kind == PairKind::Mapped);
}

TEST_CASE("pairs reader: the committed plain and bgzip inputs") {
    {
        PairsReader reader(kPairsDir + "R1_1000_all.pairs", 1);
        const PairsHeader& header = reader.header();
        CHECK(header.sorted_keys == std::vector<std::string>{"none"});
        CHECK(header.pair_type == 7);
        CHECK(header.mapq2 == 9);
        std::string_view line;
        REQUIRE(reader.next(line));
        CHECK(line.find('\n') == std::string_view::npos);
        CHECK(reader.line_number() == static_cast<std::int64_t>(header.chrom_sizes.size()) + 6);
    }
    CHECK(count_lines(kPairsDir + "R1_1000_all.pairs", 1) == 983);

    const std::string sorted = kPairsDir + "small_test_valid_5kb.pairs.gz";
    {
        PairsReader reader(sorted, 1);
        CHECK(reader.header().sorted_by_chroms_then_positions());
        CHECK(reader.header().upper_triangle());
    }
    CHECK(count_lines(sorted, 1) == 37321);
    CHECK(count_lines(sorted, 4) == 37321);

    const std::string reserved = kPairsDir + "small_test_valid_rf.pairs.gz";
    {
        PairsReader reader(reserved, 2);
        CHECK_FALSE(reader.header().columns_declared);
        CHECK(reader.header().fields_needed == 7);
    }
    CHECK(count_lines(reserved, 1) == 36627);

    // chr1 and chr2 column names, as ENCODE writes them.
    const std::string all = kPairsDir + "small_test_all.pairs.gz";
    PairsReader reader(all, 3);
    CHECK(reader.header().chrom1 == 1);
    CHECK(reader.header().sorted_by_chroms_then_positions());
    CHECK(count_lines(all, 3) == 99983);

    CHECK_THROWS_AS(PairsReader(kPairsDir + "does_not_exist.pairs", 1), std::runtime_error);
    CHECK_THROWS_WITH_AS(PairsReader(std::string(HICX_TEST_DATA_DIR) + "/DpnII.bed", 1),
                         doctest::Contains("not a .pairs file"), std::runtime_error);
}
