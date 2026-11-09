// Unit tests for core/src/merge_domains.cpp, the core of hicMergeDomains.
//
// The central test compares the DOT source of every relation tree against
// graphviz's own Digraph.source, byte for byte. The reference sources in
// hicexplorer/test/test_data/hicMergeDomains/tree_sources/ were recorded from
// the Python tool on the real domain files, with graphviz.Digraph.render
// replaced by a recorder that appends '### <filename>' and self.source for
// every render call (no Python source file was modified; the recorder patched
// the class at run time). two_files.txt and three_files_protein.txt were
// recorded again from the C++ tool when it repaired the lone nodes of the last
// chromosome and the relation that opens a chromosome (#790, #803, #839), both
// of which 3.7.6 dropped from the trees; the lists are unchanged. The recorded
// invocations, run in an empty directory:
//
//   two_files.txt            -d 10kbtad_domains.bed 50kbtad_domains.bed
//                            -om m1 -or r1 -ot tree
//   three_files_protein.txt  -d 10kbtad_domains.bed 50kbtad_domains.bed
//                            100kbtad_domains.bed -p ctcf_sorted_nochr.bed
//                            -om m2 -or r2 -ot tree
//
// hicexplorer/test/general/test_hicMergeDomains.py checks the same files
// against the Python, so a change in either implementation fails a test.
// Equivalence of the text outputs on the real corpus is the harness's job
// (cpp/scripts/cases/hicMergeDomains.json).

#include <doctest/doctest.h>

#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "hicx/merge_domains.hpp"

namespace md = hicx::merge_domains;

namespace {

const std::string kData = std::string(HICX_TEST_DATA_DIR) + "/hicMergeDomains/";

struct RecordedGraph {
    std::string filename;
    std::string source;
};

std::vector<RecordedGraph> read_recording(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());
    const std::string content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    std::vector<RecordedGraph> graphs;
    std::size_t position = 0;
    while (position < content.size()) {
        REQUIRE(content.compare(position, 4, "### ") == 0);
        const std::size_t name_end = content.find('\n', position);
        std::size_t next = content.find("\n### ", name_end);
        next = next == std::string::npos ? content.size() : next + 1;
        graphs.push_back({content.substr(position + 4, name_end - position - 4),
                          content.substr(name_end + 1, next - name_end - 1)});
        position = next;
    }
    return graphs;
}

std::vector<md::TreeGraph> trees_for(const std::vector<std::string>& files,
                                     const std::string& protein_file) {
    md::ProteinList proteins;
    if (!protein_file.empty()) {
        proteins = md::read_protein(protein_file);
    }
    md::MergedDomains domains = md::merge_domain_files(
        files, protein_file.empty() ? nullptr : &proteins, 1, 5000);
    const std::vector<md::Relation> relations =
        md::create_relationship_list(domains.pool, domains.merged, 0.5);
    std::vector<md::TreeGraph> graphs;
    md::create_tree(relations, domains.pool, domains.merged, "tree",
                    [&graphs](const md::TreeGraph& graph) { graphs.push_back(graph); });
    return graphs;
}

void check_against_recording(const std::vector<md::TreeGraph>& graphs,
                             const std::string& recording) {
    const std::vector<RecordedGraph> expected = read_recording(recording);
    REQUIRE(graphs.size() == expected.size());
    for (std::size_t i = 0; i < graphs.size(); ++i) {
        CAPTURE(expected[i].filename);
        CHECK(graphs[i].filename == expected[i].filename);
        CHECK(graphs[i].chrom == expected[i].filename.substr(5));
        // Byte identity; a CHECK on the strings would print megabytes on
        // failure, so report the first differing offset instead.
        std::size_t first_difference = 0;
        while (first_difference < graphs[i].source.size() &&
               first_difference < expected[i].source.size() &&
               graphs[i].source[first_difference] == expected[i].source[first_difference]) {
            ++first_difference;
        }
        CAPTURE(first_difference);
        CHECK(graphs[i].source.size() == expected[i].source.size());
        CHECK(graphs[i].source == expected[i].source);
    }
}

std::string temp_file_with(const std::string& content) {
    static int counter = 0;
    const std::string path = std::string(P_tmpdir) + "/hicx-merge-domains-" +
                             std::to_string(++counter) + "-" +
                             std::to_string(static_cast<long long>(::getpid())) + ".bed";
    std::ofstream out(path, std::ios::binary);
    out << content;
    return path;
}

}  // namespace

TEST_CASE("relation tree DOT sources are byte identical to graphviz, two domain files") {
    const std::vector<md::TreeGraph> graphs =
        trees_for({kData + "10kbtad_domains.bed", kData + "50kbtad_domains.bed"}, "");
    check_against_recording(graphs, kData + "tree_sources/two_files.txt");
    // hicMergeDomains.py:271 against :297: only the first graph is strict.
    REQUIRE(graphs.size() == 23);
    CHECK(graphs[0].source.rfind("strict digraph {\n", 0) == 0);
    for (std::size_t i = 1; i < graphs.size(); ++i) {
        CHECK(graphs[i].source.rfind("digraph {\n", 0) == 0);
    }
}

TEST_CASE("relation tree DOT sources are byte identical to graphviz, three files and "
          "protein peaks") {
    check_against_recording(trees_for({kData + "10kbtad_domains.bed",
                                       kData + "50kbtad_domains.bed",
                                       kData + "100kbtad_domains.bed"},
                                      kData + "ctcf_sorted_nochr.bed"),
                            kData + "tree_sources/three_files_protein.txt");
}

TEST_CASE("the domain files are merged from the finest to the coarsest bin size") {
    // 3.7.6 merged in command line order and failed with an IndexError in
    // create_tree for 100, 50, 10 kb (#790, #803, #839).
    const std::string ten = kData + "10kbtad_domains.bed";
    const std::string fifty = kData + "50kbtad_domains.bed";
    const std::string hundred = kData + "100kbtad_domains.bed";
    check_against_recording(trees_for({fifty, ten}, ""), kData + "tree_sources/two_files.txt");
    const md::MergedDomains sorted = md::merge_domain_files({ten, fifty, hundred}, nullptr, 1, 5000);
    for (const auto& order : std::vector<std::vector<std::string>>{
             {hundred, fifty, ten}, {fifty, hundred, ten}, {ten, hundred, fifty}}) {
        const md::MergedDomains other = md::merge_domain_files(order, nullptr, 1, 5000);
        REQUIRE(other.merged.size() == sorted.merged.size());
        for (std::size_t i = 0; i < sorted.merged.size(); ++i) {
            CHECK(other.pool.rows[other.merged[i]] == sorted.pool.rows[sorted.merged[i]]);
        }
        const auto relations = md::create_relationship_list(other.pool, other.merged, 0.5);
        std::size_t graphs = 0;
        md::create_tree(relations, other.pool, other.merged, "tree",
                        [&graphs](const md::TreeGraph&) { ++graphs; });
        CHECK(graphs == 23);
    }
}

TEST_CASE("the protein filter only acts when the chromosome names match") {
    // ctcf_sorted.bed names chromosomes chr1..chrX while the domain files use
    // 1..X, so compare_boundaries_protein never finds a chromosome and removes
    // nothing; the derived file without the prefix removes 1,168 TADs.
    const std::vector<std::string> files = {kData + "10kbtad_domains.bed"};
    const md::ProteinList with_prefix = md::read_protein(kData + "ctcf_sorted.bed");
    const md::ProteinList without_prefix = md::read_protein(kData + "ctcf_sorted_nochr.bed");
    CHECK(md::merge_domain_files(files, &with_prefix, 1, 5000).merged.size() == 9175);
    CHECK(md::merge_domain_files(files, &without_prefix, 1, 5000).merged.size() == 8007);
}

TEST_CASE("create_list_of_file derives the bin size from the last non zero digit") {
    md::RowPool pool;
    const std::string path = temp_file_with(
        "1\t710000\t1250000\tx\t0.1\n1\t1250000\t1480000\tx\t0.2\n1\t30\t40\tx\t0.3\n");
    const md::DomainList list = md::create_list_of_file(pool, path);
    std::remove(path.c_str());
    CHECK(list.rows.size() == 3);
    // 710000 -> 10000, 1250000 -> 50000, 30 -> 30: the smallest wins.
    CHECK(list.bin_size == 30);
}

TEST_CASE("merge_protein with a bin size of zero is refused instead of looping forever") {
    const md::ProteinList proteins = {{{"1", "5", "10"}, {"1", "50", "60"}}};
    CHECK_THROWS_AS(static_cast<void>(md::merge_protein(proteins, 0, 1)),
                    md::ReferenceNeverTerminates);
}

TEST_CASE("Python int() and float() semantics") {
    CHECK(md::py_int(" 42 ") == 42);
    CHECK(md::py_int("-7") == -7);
    CHECK(md::py_int("1_000") == 1000);
    CHECK(md::py_int("007") == 7);
    CHECK_THROWS_AS(static_cast<void>(md::py_int("")), md::PythonError);
    CHECK_THROWS_AS(static_cast<void>(md::py_int("1.0")), md::PythonError);
    CHECK_THROWS_AS(static_cast<void>(md::py_int("1__0")), md::PythonError);
    CHECK(md::py_float("-1.092317509315") == -1.092317509315);
    CHECK(md::py_float(" 1e3 ") == 1000.0);
    CHECK(md::py_float(".5") == 0.5);
    CHECK_THROWS_AS(static_cast<void>(md::py_float("0x10")), md::PythonError);
    CHECK_THROWS_AS(static_cast<void>(md::py_float(".")), md::PythonError);
    CHECK(md::is_graphviz_format("PNG"));
    CHECK_FALSE(md::is_graphviz_format("jpeg2000"));
}

TEST_CASE("read_lines uses universal newlines and rstrip") {
    const std::string path = temp_file_with("a\tb  \r\nc\rd\n\ne \t");
    const std::vector<std::string> lines = md::read_lines(path);
    std::remove(path.c_str());
    CHECK(lines == std::vector<std::string>{"a\tb", "c", "d", "", "e"});
}
