// Unit tests for the C++ cores behind the public Python names of
// hicexplorer.lib.buildMatrixMethods, hicexplorer.hicFindRestSite,
// hicexplorer.parserCommon, hicexplorer.hicBuildMatrixMicroC and
// hicexplorer.hicValidateLocations (cpp/python/API_BINDINGS.md). The expected
// values are the ones the 3.7.6 doctests declare.

#include <doctest/doctest.h>

#include <fstream>
#include <string>
#include <vector>

#include "hicx/fasta_reader.hpp"
#include "hicx/parser_common.hpp"

namespace {

std::string write_file(const std::string& name, const std::string& content) {
    const std::string path = std::string(HICX_TEST_OUTPUT_DIR) + "/" + name;
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << content;
    return path;
}

}  // namespace

TEST_CASE("find_restriction_sites: the doctests of hicFindRestSite.find_pattern") {
    const std::string fasta = write_file("api_find_pattern.fa", ">chr1\nCTACGGTACGAACGTACGGTACGcgtaCGNAGTCATG\n");
    CHECK(hicx::fasta::restriction_sites_bed(hicx::fasta::find_restriction_sites({"GTAC"}, fasta)) ==
          "chr1\t5\t9\t.\t0\t+\nchr1\t13\t17\t.\t0\t+\nchr1\t18\t22\t.\t0\t+\nchr1\t24\t28\t.\t0\t+\n");
    CHECK(hicx::fasta::restriction_sites_bed(hicx::fasta::find_restriction_sites({"CG.AG"}, fasta)) ==
          "chr1\t0\t5\t.\t0\t-\nchr1\t27\t32\t.\t0\t+\n");
    CHECK(hicx::fasta::find_restriction_sites({}, fasta).empty());
}

TEST_CASE("find_restriction_sites: sort -u keeps the first site of a key") {
    const std::string fasta = write_file("api_find_pattern_u.fa", ">b\nAAGCTT\n>a\nAAGCTT\n");
    const auto sites = hicx::fasta::find_restriction_sites({"AAGCTT", "AAGC"}, fasta);
    REQUIRE(sites.size() == 4);
    CHECK(sites[0].chrom == "a");
    CHECK(sites[0].end == 6);  // AAGCTT came first; AAGC at the same start is dropped
    CHECK(sites[1].start == 2);
    CHECK(sites[1].strand == '-');
    CHECK(sites[2].chrom == "b");
    CHECK(sites[3].start == 2);
}

TEST_CASE("parser_common: open(path, 'w').close()") {
    const std::string path = write_file("api_writable.txt", "content");
    CHECK(hicx::parser_common::truncate_for_writing(path) == 0);
    std::ifstream in(path);
    std::string text;
    std::getline(in, text);
    CHECK(text.empty());
    CHECK(hicx::parser_common::truncate_for_writing(std::string(HICX_TEST_OUTPUT_DIR)) == EISDIR);
    CHECK(hicx::parser_common::truncate_for_writing(std::string(HICX_TEST_OUTPUT_DIR) + "/no/such/file") == ENOENT);
    CHECK(hicx::parser_common::close_for_writing(987654) == EBADF);
    CHECK(hicx::parser_common::writable_file_message("x") == "x file can be opened for writting");
}
