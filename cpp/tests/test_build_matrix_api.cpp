// Unit tests for the C++ cores behind the public Python names of
// hicexplorer.lib.buildMatrixMethods, hicexplorer.hicFindRestSite,
// hicexplorer.parserCommon, hicexplorer.hicBuildMatrixMicroC and
// hicexplorer.hicValidateLocations (cpp/python/API_BINDINGS.md). The expected
// values are the ones the 3.7.6 doctests declare.

#include <doctest/doctest.h>

#include <fstream>
#include <string>
#include <vector>

#include "hicx/build_matrix_api.hpp"
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

TEST_CASE("build_matrix_api: posixpath and Tester paths") {
    namespace api = hicx::build_matrix_api;
    CHECK(api::normpath("") == ".");
    CHECK(api::normpath("/a/./b/../c//d/") == "/a/c/d");
    CHECK(api::normpath("//a/b") == "//a/b");
    CHECK(api::normpath("///a/../..") == "/");
    CHECK(api::normpath("../x/../../y") == "../../y");
    CHECK(api::abspath("x/y.py", "/cwd/dir") == "/cwd/dir/x/y.py");
    CHECK(api::dirname("/a/b/c.py") == "/a/b");
    CHECK(api::dirname("/c.py") == "/");
    CHECK(api::dirname("c.py").empty());
    CHECK(api::dirname("/a//b//c") == "/a//b");
    CHECK(api::join("/root/", "hic.bam") == "/root/hic.bam");
    CHECK(api::join("/root", "hic.bam") == "/root/hic.bam");
    CHECK(api::join("/root", "/abs") == "/abs");
    const auto lib = api::tester_paths(std::nullopt, "/pkg/hicexplorer/lib/buildMatrixMethods.py", "/", 2);
    CHECK(lib.root == "/pkg/hicexplorer/test/test_data/");
    CHECK(lib.bam_file_1 == "/pkg/hicexplorer/test/test_data/hic.bam");
    const auto micro_c = api::tester_paths(std::string(), "/pkg/hicexplorer/hicBuildMatrixMicroC.py", "/", 1);
    CHECK(micro_c.root == "/pkg/hicexplorer/test/test_data/");
    const auto env = api::tester_paths(std::string("/data"), "/x.py", "/", 1);
    CHECK(env.root == "/data");
    CHECK(env.bam_file_1 == "/data/hic.bam");
}
