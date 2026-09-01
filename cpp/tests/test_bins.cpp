#include <doctest/doctest.h>

#include <vector>

#include "hicx/bins.hpp"

using hicx::BinTable;
using hicx::CutInterval;

namespace {

std::vector<CutInterval> fixed_bins(const std::string& chrom, std::int64_t count,
                                    std::int64_t width, std::int64_t offset = 0) {
    std::vector<CutInterval> intervals;
    for (std::int64_t i = 0; i < count; ++i) {
        intervals.push_back(CutInterval{chrom, offset + i * width,
                                        offset + (i + 1) * width, 1.0, ""});
    }
    return intervals;
}

}  // namespace

TEST_CASE("chromosome boundaries follow intervalListToIntervalTree") {
    std::vector<CutInterval> intervals = fixed_bins("chr1", 3, 10);
    for (const CutInterval& interval : fixed_bins("chr2", 2, 10)) {
        intervals.push_back(interval);
    }
    const BinTable table(intervals);

    REQUIRE(table.chrom_bin_boundaries().size() == 2);
    CHECK(table.chrom_bin_boundaries()[0].first == "chr1");
    CHECK(table.chrom_bin_boundaries()[0].second.first == 0);
    CHECK(table.chrom_bin_boundaries()[0].second.last == 3);
    CHECK(table.chrom_bin_boundaries()[1].first == "chr2");
    CHECK(table.chrom_bin_boundaries()[1].second.first == 3);
    CHECK(table.chrom_bin_boundaries()[1].second.last == 5);

    const std::vector<std::string> names = table.chrom_names();
    CHECK(names == std::vector<std::string>{"chr1", "chr2"});
}

TEST_CASE("chromosome sizes use the end of the last bin") {
    std::vector<CutInterval> intervals{
        CutInterval{"chr1", 0, 10, 1.0, ""},
        CutInterval{"chr1", 10, 20, 1.0, ""},
        CutInterval{"chr1", 20, 25, 1.0, ""},  // truncated last bin
        CutInterval{"chr2", 5, 15, 1.0, ""},
    };
    const BinTable table(intervals);

    const auto sizes = table.chromosome_sizes();
    REQUIRE(sizes.size() == 2);
    CHECK(sizes[0].first == "chr1");
    CHECK(sizes[0].second == 25);
    CHECK(sizes[1].second == 15);

    // get_chromosome_sizes_real subtracts the start of the first bin.
    const auto real_sizes = table.chromosome_sizes_real();
    CHECK(real_sizes[0].second == 26);
    CHECK(real_sizes[1].second == 11);
}

TEST_CASE("bin size is the median difference between consecutive starts") {
    const BinTable uniform(fixed_bins("chr1", 10, 1000));
    CHECK(uniform.bin_size() == 1000);
    CHECK(uniform.bin_size_homogeneous());

    // A single bin returns its own width.
    const BinTable single({CutInterval{"chr1", 0, 4242, 1.0, ""}});
    CHECK(single.bin_size() == 4242);

    // Restriction fragment style bins: the median is used and the table is
    // reported as not homogeneous.
    std::vector<CutInterval> variable{
        CutInterval{"chrX", 0, 100, 1.0, ""},   CutInterval{"chrX", 100, 350, 1.0, ""},
        CutInterval{"chrX", 350, 400, 1.0, ""}, CutInterval{"chrX", 400, 900, 1.0, ""},
        CutInterval{"chrX", 900, 1000, 1.0, ""}};
    const BinTable fragments(variable);
    // Differences between starts: 100, 250, 50, 500 -> median (100+250)/2 = 175
    CHECK(fragments.bin_size() == 175);
    CHECK_FALSE(fragments.bin_size_homogeneous());
}

TEST_CASE("region lookup returns the bins overlapping a position") {
    const BinTable table(fixed_bins("chr1", 5, 100));

    CHECK(table.bin_at("chr1", 0).value() == 0);
    CHECK(table.bin_at("chr1", 99).value() == 0);
    CHECK(table.bin_at("chr1", 100).value() == 1);
    CHECK(table.bin_at("chr1", 450).value() == 4);
    CHECK_FALSE(table.bin_at("chr1", 500).has_value());
    CHECK_FALSE(table.bin_at("chr2", 0).has_value());

    const auto range = table.region_bin_range("chr1", 150, 320);
    REQUIRE(range.has_value());
    CHECK(range->first == 1);
    CHECK(range->second == 3);
    // A position past the last bin gives no range, like the Python IndexError
    // path that returns None.
    CHECK_FALSE(table.region_bin_range("chr1", 150, 5000).has_value());
}
