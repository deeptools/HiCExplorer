// Unit tests for the hicBuildMatrix machinery.
//
// The expected values of the binning tests are the ones the Python doctests in
// hicexplorer/lib/buildMatrixMethods.py declare, so a divergence in either
// direction is caught here rather than in a whole tool run.

#include <doctest/doctest.h>

#include <cmath>
#include <string>
#include <vector>

#include "hicx/bam_file.hpp"
#include "hicx/build_matrix.hpp"

using hicx::BinSearchIndex;
using hicx::ChromNames;
using hicx::ChromSizes;
using hicx::GenomeInterval;
using hicx::RestrictionSiteIndex;

namespace {

std::vector<std::string> as_text(const std::vector<GenomeInterval>& intervals,
                                 const ChromNames& names) {
    std::vector<std::string> out;
    out.reserve(intervals.size());
    for (const auto& interval : intervals) {
        out.push_back(names.name(interval.chrom) + ":" +
                      std::to_string(interval.start) + "-" +
                      std::to_string(interval.end));
    }
    return out;
}

std::vector<GenomeInterval> sites(ChromNames& names,
                                  const std::vector<std::tuple<const char*, std::int64_t,
                                                               std::int64_t>>& items) {
    std::vector<GenomeInterval> out;
    for (const auto& [chrom, start, end] : items) {
        out.push_back({names.intern(chrom), start, end});
    }
    return out;
}

}  // namespace

TEST_CASE("get_bins matches the buildMatrixMethods doctest") {
    ChromNames names;
    const ChromSizes chrom_sizes{{"contig-1", 7125}, {"contig-2", 3345}};
    const auto bins = hicx::get_bins(50000, chrom_sizes, "", names);
    CHECK(as_text(bins, names) ==
          std::vector<std::string>{"contig-1:0-7125", "contig-2:0-3345"});

    const auto one = hicx::get_bins(50000, chrom_sizes, "contig-1", names);
    CHECK(as_text(one, names) == std::vector<std::string>{"contig-1:0-7125"});
}

TEST_CASE("get_bins cuts a chromosome that is not a multiple of the bin size") {
    ChromNames names;
    const ChromSizes chrom_sizes{{"chr1", 25}};
    const auto bins = hicx::get_bins(10, chrom_sizes, "", names);
    CHECK(as_text(bins, names) ==
          std::vector<std::string>{"chr1:0-10", "chr1:10-20", "chr1:20-25"});
}

TEST_CASE("get_rf_bins matches the buildMatrixMethods doctest") {
    ChromNames names;
    const auto cut_sites = sites(names, {{"chr1", 10, 20},
                                         {"chr1", 60, 70},
                                         {"chr2", 20, 30},
                                         {"chr2", 40, 50},
                                         {"chr2", 70, 80}});
    const auto bins = hicx::get_rf_bins(cut_sites, 10, 20);
    CHECK(as_text(bins, names) == std::vector<std::string>{"chr1:0-40", "chr1:40-90",
                                                           "chr2:0-60", "chr2:60-100"});
}

TEST_CASE("enlarge_bins matches the buildMatrixMethods doctest") {
    ChromNames names;
    auto bins = sites(names, {{"chr1", 10, 30},
                              {"chr1", 50, 80},
                              {"chr2", 10, 60},
                              {"chr2", 60, 90}});
    const ChromSizes chrom_sizes{{"chr1", 100}, {"chr2", 100}};
    hicx::enlarge_bins(bins, chrom_sizes, names);
    CHECK(as_text(bins, names) == std::vector<std::string>{"chr1:0-40", "chr1:40-100",
                                                           "chr2:0-60", "chr2:60-100"});
}

TEST_CASE("getUserRegion accepts a bare chromosome and rejects an unknown one") {
    const ChromSizes chrom_sizes{{"chr2", 1000}};
    const auto region = hicx::user_region(chrom_sizes, "chr2:10:1001");
    CHECK(region.start == 10);
    CHECK(region.end == 1000);  // clipped to the chromosome length
    CHECK(region.chrom_sizes.size() == 1);
    CHECK(region.chrom_sizes[0].second == 1000);

    // The failure mode that makes test_build_matrix_restrictionCutFile_two fail
    // on master: a chromosome name that is not in the BAM header.
    CHECK_THROWS_AS(hicx::user_region(chrom_sizes, "Chr2"), std::runtime_error);
}

TEST_CASE("genomicRegion strips punctuation and turns a dash into a colon") {
    CHECK(hicx::normalise_region("chr1:1,000-2,000") == "chr1:1000:2000");
    CHECK(hicx::normalise_region(" chr X ") == "chrX");
    CHECK(hicx::normalise_region("") == "");
}

TEST_CASE("the bin lookup reproduces the shared array binary search") {
    ChromNames names;
    const auto bins = sites(names, {{"chr1", 0, 10}, {"chr1", 10, 20}, {"chr1", 20, 30},
                                    {"chr2", 0, 5}});
    const BinSearchIndex index(bins, names);
    const std::uint32_t chr1 = *names.lookup("chr1");
    const std::uint32_t chr2 = *names.lookup("chr2");

    CHECK(index.bin_at(chr1, 0) == 0);
    CHECK(index.bin_at(chr1, 5) == 0);
    CHECK(index.bin_at(chr1, 15) == 1);
    CHECK(index.bin_at(chr1, 25) == 2);
    CHECK(index.bin_at(chr2, 3) == 3);
    // Past the last bin of the chromosome.
    CHECK(!index.bin_at(chr1, 100).has_value());
    // A chromosome with no bins at all is the KeyError branch.
    CHECK(!index.has_chrom(names.intern("chr3")));
    CHECK(!index.bin_at(names.intern("chr3"), 0).has_value());
}

TEST_CASE("the bin lookup treats both interval ends as inclusive") {
    // buildMatrixMethods.py:663 compares begin <= x <= end, so a position on a
    // shared boundary can land in either neighbour depending on which one the
    // binary search probes first. With three bins the search starts at the
    // middle one, so position 10 and position 20 both resolve to bin 1.
    ChromNames names;
    const auto bins = sites(names, {{"chr1", 0, 10}, {"chr1", 10, 20}, {"chr1", 20, 30}});
    const BinSearchIndex index(bins, names);
    const std::uint32_t chr1 = *names.lookup("chr1");
    CHECK(index.bin_at(chr1, 10) == 1);
    CHECK(index.bin_at(chr1, 20) == 1);
}

TEST_CASE("the restriction site index answers the intervaltree slice query") {
    ChromNames names;
    const auto cut_sites = sites(names, {{"chr1", 100, 104}, {"chr1", 200, 204},
                                         {"chr2", 50, 54}});
    const RestrictionSiteIndex index(cut_sites, names.size());
    const std::uint32_t chr1 = *names.lookup("chr1");
    const std::uint32_t chr2 = *names.lookup("chr2");

    CHECK(index.has_chrom(chr1));
    CHECK(index.overlaps(chr1, 90, 110));
    CHECK(index.overlaps(chr1, 103, 150));
    // tree[a:b] needs begin < b and end > a.
    CHECK(!index.overlaps(chr1, 104, 200));
    CHECK(!index.overlaps(chr1, 0, 100));
    // A null slice returns an empty set.
    CHECK(!index.overlaps(chr1, 150, 150));
    CHECK(!index.overlaps(chr1, 200, 150));
    CHECK(!index.overlaps(chr2, 0, 40));
    CHECK(index.overlaps(chr2, 0, 51));
    CHECK(!index.has_chrom(names.intern("chr3")));
}

TEST_CASE("reverse_complement follows Bio.Seq over the IUPAC alphabet") {
    CHECK(hicx::reverse_complement("GATC") == "GATC");
    CHECK(hicx::reverse_complement("AGCT") == "AGCT");
    CHECK(hicx::reverse_complement("AAGCTT") == "AAGCTT");
    CHECK(hicx::reverse_complement("AACCGGTT") == "AACCGGTT");
    CHECK(hicx::reverse_complement("ACGTN") == "NACGT");
    CHECK(hicx::reverse_complement("ATTG") == "CAAT");
}

TEST_CASE("the QC counters combine as an order independent integer sum") {
    hicx::QcCounters a;
    a.dangling_end = {3, 0};
    a.self_circle = 5;
    a.pair_added = 11;
    a.iter_num = 100;
    hicx::QcCounters b;
    b.dangling_end = {0, 7};
    b.self_circle = 2;
    b.pair_added = 4;
    b.iter_num = 50;

    hicx::QcCounters forwards;
    forwards.add(a);
    forwards.add(b);
    hicx::QcCounters backwards;
    backwards.add(b);
    backwards.add(a);

    CHECK(forwards.dangling_end == backwards.dangling_end);
    CHECK(forwards.dangling_end == std::vector<std::int64_t>{3, 7});
    CHECK(forwards.self_circle == backwards.self_circle);
    CHECK(forwards.pair_added == 15);
    CHECK(forwards.iter_num == 150);
}

TEST_CASE("the QC log reproduces the reference layout") {
    hicx::QcLogInputs inputs;
    inputs.out_file_name = "/tmp/example.h5";
    inputs.min_distance = 300;
    inputs.max_library_insert_size = 1000;
    inputs.keep_self_ligation = false;
    inputs.has_restriction_cut_file = true;
    inputs.dangling_sequences = {{"GATC", "GATC"}};

    hicx::QcCounters counters;
    counters.iter_num = 99983;
    counters.one_mate_unmapped = 8777;
    counters.one_mate_not_unique = 3603;
    counters.one_mate_low_quality = 34877;
    counters.dangling_end = {54};
    counters.self_ligation = 5056;
    counters.mate_not_close_to_rf = 0;
    counters.same_fragment = 10283;
    counters.self_circle = 651;
    counters.duplicated_pairs = 12;
    counters.pair_added = 37321;
    counters.inter_chromosomal = 5955;
    counters.short_range = 8853;
    counters.long_range = 22513;
    counters.count_inward = 7145;
    counters.count_outward = 9731;
    counters.count_left = 7156;
    counters.count_right = 7334;

    const std::string log = hicx::format_qc_log(inputs, counters);
    CHECK(log.rfind("\nFile\t/tmp/example.h5\t\t\n", 0) == 0);
    CHECK(log.find("Sequenced reads\t99983\t\t\n") != std::string::npos);
    CHECK(log.find("Min rest. site distance\t300\t\t\n") != std::string::npos);
    CHECK(log.find("Pairs mappable, unique and high quality\t52726\t(52.73)\n") !=
          std::string::npos);
    CHECK(log.find("Hi-C contacts\t37321\t(37.33)\n") != std::string::npos);
    CHECK(log.find("dangling end GATC (restriction sequence GATC)\t54\t(0.10)\n") !=
          std::string::npos);
    CHECK(log.find("self ligation (removed)\t5056\t(9.59)\n") != std::string::npos);
    CHECK(log.find("self circle\t651\t(1.23)\n") != std::string::npos);
    CHECK(log.find("Read pair type: right pairs\t7334\t(19.65)\n") !=
          std::string::npos);
    // Without a restriction cut file three lines disappear from the log.
    inputs.has_restriction_cut_file = false;
    const std::string without = hicx::format_qc_log(inputs, counters);
    CHECK(without.find("self ligation") == std::string::npos);
    CHECK(without.find("self circle") == std::string::npos);
    CHECK(without.find("One mate not close to rest site") == std::string::npos);
    CHECK(without.find("same fragment\t10283\t(19.50)\n") != std::string::npos);
    // --keepSelfLigation only changes the label.
    inputs.has_restriction_cut_file = true;
    inputs.keep_self_ligation = true;
    CHECK(hicx::format_qc_log(inputs, counters).find("self ligation (not removed)") !=
          std::string::npos);
}

TEST_CASE("the BAM reader reports pysam's fields on a real file") {
    const std::string path = std::string(HICX_TEST_DATA_DIR) + "/R1_1000.bam";
    hicx::BamReader reader(path);
    CHECK(reader.references().size() == 15);
    CHECK(reader.references().front() == "chr2L");
    CHECK(reader.lengths().front() == 23011544);

    const auto sizes = hicx::chrom_sizes_of(reader);
    CHECK(sizes.size() == 15);
    CHECK(sizes.front().first == "chr2L");
    CHECK(sizes.back().first == "chrYHet");

    hicx::BamRecord record;
    std::int64_t records = 0;
    std::int64_t mapped = 0;
    while (reader.read(record)) {
        const hicx::ReadFields fields = record.fields();
        CHECK(fields.seq_len == 50);
        CHECK(fields.qlen <= fields.seq_len);
        if (!fields.is_unmapped()) {
            ++mapped;
            CHECK(fields.tid >= 0);
            CHECK(fields.pos >= 0);
        }
        ++records;
    }
    // Fixed by the file, so a change in the reader shows up immediately.
    CHECK(records == 983);
    CHECK(mapped == 542);
}
