// Unit tests for the core code the tier 2 tools brought in:
//
//   * hicx::TextTable, pandas.read_csv/to_csv with its dtype inference and its
//     precise_xstrtod float parser
//   * hicx::bedtools, the sort, merge and intersect -c that replace pybedtools
//   * hicx::fasta, the FASTA reader and the IUPAC reverse complement
//   * hicx::npz, the .npy and .npz reader and writer
//
// Contract rule 2 applies: these are mechanics, and the equivalence against
// the Python is checked on the real corpus by cpp/scripts/equiv.py. What is
// pinned here is the behaviour a hand written five line file can show and a
// 12,000 line one cannot: the exact digits precise_xstrtod produces, the
// permutation std::sort leaves behind, and the byte layout of a .npy header.

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "hicx/bedtools_ops.hpp"
#include "hicx/fasta_reader.hpp"
#include "hicx/npz_file.hpp"
#include "hicx/text_table.hpp"

namespace {

class TempFile {
  public:
    explicit TempFile(const std::string& suffix) {
        static int counter = 0;
        path_ = (std::filesystem::temp_directory_path() /
                 ("hicx-tier2-" + std::to_string(++counter) + "-" +
                  std::to_string(::getpid()) + suffix))
                    .string();
    }
    ~TempFile() { std::remove(path_.c_str()); }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }
    void write(const std::string& content) const {
        std::ofstream out(path_, std::ios::binary);
        out << content;
    }

  private:
    std::string path_;
};

}  // namespace

TEST_CASE("pandas_strtod keeps 17 significant digits and then scales") {
    bool ok = false;

    // The value that made this function necessary. pandas' precise_xstrtod
    // accumulates the first 17 significant digits into a double and divides by
    // an exact power of ten, so it does not produce the correctly rounded
    // conversion strtod gives; hicValidateLocations writes the shorter repr
    // that results.
    const double parsed = hicx::pandas_strtod("0.0019430592210407135", &ok);
    CHECK(ok);
    CHECK(parsed == 19430592210407.0 / 1e16);
    CHECK(parsed != std::strtod("0.0019430592210407135", nullptr));

    // Values with at most 17 significant digits are unaffected.
    CHECK(hicx::pandas_strtod("0.001", &ok) == 0.001);
    CHECK(ok);
    CHECK(hicx::pandas_strtod("1e-05", &ok) == 1e-05);
    CHECK(ok);
    CHECK(hicx::pandas_strtod("-12.5", &ok) == -12.5);
    CHECK(ok);
    CHECK(hicx::pandas_strtod("  3.5  ", &ok) == 3.5);
    CHECK(ok);

    // Not numbers.
    CHECK(hicx::pandas_strtod("chr1", &ok) == 0.0);
    CHECK_FALSE(ok);
    CHECK(hicx::pandas_strtod("1.5x", &ok) == 1.5);
    CHECK_FALSE(ok);
    CHECK(hicx::pandas_strtod("", &ok) == 0.0);
    CHECK_FALSE(ok);
}

TEST_CASE("TextTable infers one dtype per column and writes it back") {
    TempFile file(".bed");
    file.write("chr1\t10\t20\t0.5\n"
               "chrX\t30\t40\t2\n");

    const hicx::TextTable table = hicx::TextTable::read_tsv(file.path());
    REQUIRE(table.rows() == 2);
    REQUIRE(table.cols() == 4);
    CHECK(table.column(0).type == hicx::ColumnType::String);
    CHECK(table.column(1).type == hicx::ColumnType::Int64);
    CHECK(table.column(2).type == hicx::ColumnType::Int64);
    // One float in the column makes the whole column float, and the integer
    // 2 is then written back as 2.0.
    CHECK(table.column(3).type == hicx::ColumnType::Float64);
    CHECK(table.to_tsv() == "chr1\t10\t20\t0.5\nchrX\t30\t40\t2.0\n");

    // A column of whole numbers stays integral and keeps its spelling.
    TempFile integers(".bed");
    integers.write("chr1\t10\t20\t3\n");
    CHECK(hicx::TextTable::read_tsv(integers.path()).to_tsv() == "chr1\t10\t20\t3\n");
}

TEST_CASE("TextTable reproduces the pandas round trip and the chr prefix rules") {
    TempFile file(".bed");
    file.write("1\t10\t20\n"
               "X\t30\t40\n");
    hicx::TextTable table = hicx::TextTable::read_tsv(file.path());
    // The chromosome column holds 'X', so it is an object column already.
    CHECK(table.column(0).type == hicx::ColumnType::String);

    table.add_chr_prefix(0);
    CHECK(table.to_tsv() == "chr1\t10\t20\nchrX\t30\t40\n");

    // lstrip('chr') strips a character set, not a prefix: 'chrX' loses the
    // leading c, h and r and keeps the X, and a name made only of those
    // characters would be emptied.
    table.remove_chr_prefix(0);
    CHECK(table.to_tsv() == "1\t10\t20\nX\t30\t40\n");

    // reparse re-infers the dtypes from the text, which is what
    // BedTool.from_dataframe(...).to_dataframe() does. An all-numeric
    // chromosome column comes back as int64 and loses nothing here.
    TempFile numeric(".bed");
    numeric.write("1\t10\t20\n2\t30\t40\n");
    hicx::TextTable plain = hicx::TextTable::read_tsv(numeric.path());
    CHECK(plain.column(0).type == hicx::ColumnType::Int64);
    plain.add_chr_prefix(0);
    CHECK(plain.column(0).type == hicx::ColumnType::String);
    CHECK(plain.reparse().column(0).type == hicx::ColumnType::String);
}

TEST_CASE("TextTable drops duplicates in both pandas modes") {
    TempFile file(".bed");
    file.write("chr1\t10\t20\n"
               "chr1\t10\t20\n"
               "chr2\t30\t40\n");

    hicx::TextTable first = hicx::TextTable::read_tsv(file.path());
    first.drop_duplicates();
    CHECK(first.to_tsv() == "chr1\t10\t20\nchr2\t30\t40\n");

    // keep=False removes every copy, which is what hicMergeLoops asks for.
    hicx::TextTable all = hicx::TextTable::read_tsv(file.path());
    all.drop_duplicates_drop_all();
    CHECK(all.to_tsv() == "chr2\t30\t40\n");
}

TEST_CASE("bedtools sort orders chromosomes by byte and starts numerically") {
    const std::vector<std::string> chrom = {"chr2", "chr10", "chr1", "chr1"};
    const std::vector<std::int64_t> start = {5, 5, 10000, 402};

    const std::vector<std::size_t> order = hicx::bedtools::sort_order(chrom, start);
    // std::map order over the names is chr1, chr10, chr2, which is bedtools'
    // lexicographic order and not a natural one; and within chr1 the start is
    // compared as a number, so 402 comes before 10000.
    CHECK(order == std::vector<std::size_t>{3, 2, 1, 0});
}

TEST_CASE("bedtools sort is deterministic and keeps every record") {
    // A start tie, where bedtools' comparator does not look at any other
    // field. The permutation is whatever std::sort leaves; what is asserted
    // here is that it is a permutation, that it is stable for a two element
    // run, and that repeating the call gives the same answer, which is the
    // property the tools depend on.
    const std::vector<std::string> chrom = {"c", "c", "c"};
    const std::vector<std::int64_t> start = {7, 7, 3};
    const std::vector<std::size_t> first = hicx::bedtools::sort_order(chrom, start);
    const std::vector<std::size_t> second = hicx::bedtools::sort_order(chrom, start);
    CHECK(first == second);
    CHECK(first.size() == 3);
    CHECK(first[0] == 2);
    std::vector<std::size_t> sorted = first;
    std::sort(sorted.begin(), sorted.end());
    CHECK(sorted == std::vector<std::size_t>{0, 1, 2});
}

TEST_CASE("bedtools merge joins overlapping and book-ended intervals") {
    hicx::bedtools::Intervals input;
    input.push_back("chr1", 0, 10);
    input.push_back("chr1", 10, 20);   // book-ended, merged at -d 0
    input.push_back("chr1", 25, 30);
    input.push_back("chr1", 28, 40);   // overlapping
    input.push_back("chr2", 0, 5);

    const hicx::bedtools::Intervals merged = hicx::bedtools::merge(input);
    REQUIRE(merged.size() == 3);
    CHECK(merged.chrom[0] == "chr1");
    CHECK(merged.start[0] == 0);
    CHECK(merged.end[0] == 20);
    CHECK(merged.start[1] == 25);
    CHECK(merged.end[1] == 40);
    CHECK(merged.chrom[2] == "chr2");
}

TEST_CASE("bedtools intersect counts overlaps of at least one base") {
    hicx::bedtools::Intervals peaks;
    peaks.push_back("chr1", 10, 20);
    peaks.push_back("chr1", 30, 40);
    peaks.push_back("chr2", 0, 100);

    hicx::bedtools::Intervals queries;
    queries.push_back("chr1", 0, 10);    // touching only, no overlap
    queries.push_back("chr1", 19, 21);   // one base of the first peak
    queries.push_back("chr1", 15, 35);   // both peaks
    queries.push_back("chr3", 0, 100);   // chromosome absent from B
    queries.push_back("chr2", 99, 200);

    const std::vector<std::int64_t> counts =
        hicx::bedtools::intersect_count(queries, peaks);
    CHECK(counts == std::vector<std::int64_t>{0, 1, 2, 0, 1});
}

TEST_CASE("fasta reverse complement follows the Bio.Seq table") {
    CHECK(hicx::fasta::reverse_complement("AAGCTT") == "AAGCTT");
    CHECK(hicx::fasta::reverse_complement("GATC") == "GATC");
    CHECK(hicx::fasta::reverse_complement("GGATC") == "GATCC");
    // The lowercase ambiguity codes are complemented too, which is where the
    // table differs from a switch over the four bases.
    CHECK(hicx::fasta::reverse_complement("rywskmhbvdn") == "nhbvdkmswry");
    // hicFindRestSite hands regular expressions to this function, so the
    // characters outside the table have to pass through and end up mirrored.
    CHECK(hicx::fasta::reverse_complement("CG..GC") == "GC..CG");
    CHECK(hicx::fasta::reverse_complement("CG.AG") == "CT.CG");
    CHECK(hicx::fasta::reverse_complement("A[CG]T") == "A]CG[T");
}

TEST_CASE("fasta reader takes the name up to the first space and folds lines") {
    TempFile file(".fasta");
    file.write(">chrM some description here\n"
               "ACGT acgt\n"
               "NNNN\n"
               ">chr2\n"
               "TTTT\n");

    std::vector<std::pair<std::string, std::string>> records;
    hicx::fasta::read_fasta(file.path(), false,
                            [&](const std::string& name, const std::string& sequence) {
                                records.emplace_back(name, sequence);
                            });
    REQUIRE(records.size() == 2);
    CHECK(records[0].first == "chrM");
    // Spaces inside a sequence line are removed, as SimpleFastaParser does.
    CHECK(records[0].second == "ACGTacgtNNNN");
    CHECK(records[1].first == "chr2");
    CHECK(records[1].second == "TTTT");

    // The gzip decision is made from the file name alone, never from the
    // magic bytes, because that is what mimetypes.guess_type does.
    CHECK(hicx::fasta::looks_gzipped("genome.fa.gz"));
    CHECK_FALSE(hicx::fasta::looks_gzipped("genome.fa.bz2"));
    CHECK_FALSE(hicx::fasta::looks_gzipped("genome.fa"));
}

TEST_CASE("npy header matches numpy's byte for byte") {
    // Captured from numpy 1.26.4: the header is padded with spaces so that the
    // whole thing is a multiple of 64 bytes and ends in a newline, and the
    // length field counts the padded body including that newline.
    const std::string header = hicx::npz::npy_header("<i4", {191});
    CHECK(header.size() == 128);
    CHECK(header.compare(0, 8, std::string("\x93NUMPY\x01\x00", 8)) == 0);
    CHECK(static_cast<unsigned char>(header[8]) == 118);
    CHECK(static_cast<unsigned char>(header[9]) == 0);
    CHECK(header.find("{'descr': '<i4', 'fortran_order': False, 'shape': (191,), }") ==
          10);
    CHECK(header.back() == '\n');

    // A zero dimensional array, which is how the csr format tag is stored.
    const std::string scalar = hicx::npz::npy_header("|S3", {});
    CHECK(scalar.size() == 128);
    CHECK(scalar.find("{'descr': '|S3', 'fortran_order': False, 'shape': (), }") ==
          10);

    // A two element shape keeps ", " between the entries and no trailing comma.
    const std::string pair = hicx::npz::npy_header("<i8", {2, 3});
    CHECK(pair.find("{'descr': '<i8', 'fortran_order': False, 'shape': (2, 3), }") ==
          10);
}

TEST_CASE("npz round trips a csr matrix through its own reader") {
    TempFile file(".npz");
    const std::vector<std::int32_t> indptr = {0, 2, 3};
    const std::vector<std::int32_t> indices = {0, 2, 1};
    const std::vector<double> data = {1.5, 2.5, 3.5};
    hicx::npz::save_csr_npz(file.path(), 2, 3, indptr, indices, data);

    const std::vector<hicx::npz::Array> arrays = hicx::npz::read_npz(file.path());
    REQUIRE(arrays.size() == 5);
    // The entry order is scipy.sparse.save_npz's dictionary order.
    CHECK(arrays[0].name == "indices");
    CHECK(arrays[1].name == "indptr");
    CHECK(arrays[2].name == "format");
    CHECK(arrays[3].name == "shape");
    CHECK(arrays[4].name == "data");

    CHECK(arrays[2].dtype == "|S3");
    CHECK(arrays[2].data == "csr");
    CHECK(arrays[2].shape.empty());

    CHECK(arrays[3].dtype == "<i8");
    std::int64_t shape[2] = {0, 0};
    std::memcpy(shape, arrays[3].data.data(), sizeof(shape));
    CHECK(shape[0] == 2);
    CHECK(shape[1] == 3);

    CHECK(arrays[4].dtype == "<f8");
    REQUIRE(arrays[4].shape == std::vector<std::int64_t>{3});
    std::vector<double> read_back(3);
    std::memcpy(read_back.data(), arrays[4].data.data(), 3 * sizeof(double));
    CHECK(read_back == data);

    CHECK(arrays[0].dtype == "<i4");
    std::vector<std::int32_t> read_indices(3);
    std::memcpy(read_indices.data(), arrays[0].data.data(), 3 * sizeof(std::int32_t));
    CHECK(read_indices == indices);
}

TEST_CASE("npz reads an uncompressed archive as well") {
    TempFile file(".npz");
    std::vector<hicx::npz::Array> arrays;
    arrays.push_back(hicx::npz::Array{"x", "<i4", {2}, std::string("\1\0\0\0\2\0\0\0", 8)});
    hicx::npz::write_npz(file.path(), arrays, false);

    const std::vector<hicx::npz::Array> back = hicx::npz::read_npz(file.path());
    REQUIRE(back.size() == 1);
    CHECK(back[0].name == "x");
    CHECK(back[0].dtype == "<i4");
    CHECK(back[0].shape == std::vector<std::int64_t>{2});
    CHECK(back[0].data.size() == 8);
}
