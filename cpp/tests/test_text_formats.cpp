// The homer, ginteractions, hicpro and 2D-text readers and writers.
//
// Small hand written files pin the mechanics: which cells a line touches, how
// a value is spelled, where a newline is and is not, and which of the two
// triangles a writer emits. The equivalence against Python written files is
// checked on the real matrices of the corpus by cpp/scripts/equiv.py, which is
// where contract rule 2 applies.

#include <doctest/doctest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "hicx/text_formats.hpp"

namespace {

class TempFile {
  public:
    explicit TempFile(const std::string& suffix) {
        path_ = (std::filesystem::temp_directory_path() /
                 ("hicx-text-" + std::to_string(++counter_) + "-" +
                  std::to_string(::getpid()) + suffix))
                    .string();
    }
    ~TempFile() {
        std::remove(path_.c_str());
        std::remove((path_ + ".tsv").c_str());
    }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;
    [[nodiscard]] const std::string& path() const { return path_; }

  private:
    static inline int counter_ = 0;
    std::string path_;
};

void write_text(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    out << content;
}

std::string read_text(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// Three bins on two chromosomes, upper triangle of
//   5 2 0
//   2 0 3
//   0 3 7
hicx::MatrixData toy_matrix(const std::string& dtype) {
    hicx::MatrixData data;
    const std::vector<std::int32_t> row{0, 0, 1, 2};
    const std::vector<std::int32_t> col{0, 1, 2, 2};
    std::vector<double> values{5.0, 2.0, 3.0, 7.0};
    data.matrix = hicx::CsrMatrix::from_coo(3, 3, row, col, std::move(values), dtype);
    data.cut_intervals = {
        hicx::CutInterval{"chr1", 0, 10, 1.0, ""},
        hicx::CutInterval{"chr1", 10, 20, 1.0, ""},
        hicx::CutInterval{"chr2", 0, 10, 1.0, ""},
    };
    return data;
}

}  // namespace

TEST_CASE("value_repr spells a value the way str() of a numpy scalar does") {
    CHECK(hicx::value_repr(0.0, "int32") == "0");
    CHECK(hicx::value_repr(7.0, "int64") == "7");
    CHECK(hicx::value_repr(0.0, "float64") == "0.0");
    CHECK(hicx::value_repr(1.3e-05, "float64") == "1.3e-05");
    CHECK(hicx::value_repr(1054.4748299674175, "float64") == "1054.4748299674175");
    CHECK(hicx::value_repr(0.0, "float32") == "0.0");
    CHECK(hicx::value_repr(0.1, "float32") == "0.1");
}

TEST_CASE("maximum_with_transpose mirrors without doubling the diagonal") {
    const hicx::CsrMatrix matrix = toy_matrix("int32").matrix;
    const hicx::CsrMatrix symmetric = hicx::maximum_with_transpose(matrix);
    CHECK(symmetric.at(0, 0) == 5.0);  // a sum would give 10
    CHECK(symmetric.at(0, 1) == 2.0);
    CHECK(symmetric.at(1, 0) == 2.0);
    CHECK(symmetric.at(1, 2) == 3.0);
    CHECK(symmetric.at(2, 1) == 3.0);
    CHECK(symmetric.at(2, 2) == 7.0);
    CHECK(symmetric.stored_nnz() == 6);
}

TEST_CASE("maximum_with_transpose clamps a negative off diagonal value") {
    // scipy's maximum against the structural zero of the mirror cell keeps the
    // larger of the two, which for a negative value is the zero.
    const std::vector<std::int32_t> row{0, 0, 1};
    const std::vector<std::int32_t> col{0, 1, 1};
    std::vector<double> values{-4.0, -3.0, 2.0};
    const hicx::CsrMatrix matrix =
        hicx::CsrMatrix::from_coo(2, 2, row, col, std::move(values), "float64");
    const hicx::CsrMatrix symmetric = hicx::maximum_with_transpose(matrix);
    CHECK(symmetric.at(0, 0) == -4.0);  // the diagonal keeps its sign
    CHECK(symmetric.at(0, 1) == 0.0);
    CHECK(symmetric.at(1, 0) == 0.0);
    CHECK(symmetric.stored_nnz() == 2);
}

TEST_CASE("upper_triangle_after_maximum equals the triangle of the mirror") {
    const hicx::CsrMatrix matrix = toy_matrix("int32").matrix;
    const hicx::CsrMatrix full = hicx::maximum_with_transpose(matrix);
    const hicx::CsrMatrix upper = hicx::upper_triangle_after_maximum(matrix);
    for (std::int64_t row = 0; row < matrix.rows(); ++row) {
        for (std::int64_t col = row; col < matrix.cols(); ++col) {
            CHECK(upper.at(row, col) == full.at(row, col));
        }
        for (std::int64_t col = 0; col < row; ++col) {
            CHECK(upper.at(row, col) == 0.0);
        }
    }
}

TEST_CASE("the homer writer emits the header, no trailing newline and str()") {
    TempFile file(".homer");
    hicx::write_homer(file.path(), toy_matrix("int32"));
    // The writer always compresses; the reader sniffs the magic bytes.
    const hicx::MatrixData back = hicx::read_homer(file.path());
    REQUIRE(back.cut_intervals.size() == 3);
    CHECK(back.cut_intervals[0].chrom == "chr1");
    CHECK(back.cut_intervals[0].start == 0);
    // The bin size comes from the first two names, so every bin gets width 10
    // even though the third one starts a new chromosome at 0.
    CHECK(back.cut_intervals[0].end == 10);
    CHECK(back.cut_intervals[2].chrom == "chr2");
    CHECK(back.cut_intervals[2].start == 0);
    // The dense table holds the stored triangle only, because write_homer was
    // handed the matrix as it is; hicConvertFormat symmetrises first.
    CHECK(back.matrix.at(0, 0) == 5.0);
    CHECK(back.matrix.at(0, 1) == 2.0);
    CHECK(back.matrix.at(1, 0) == 0.0);
    CHECK(back.matrix.at(2, 2) == 7.0);
}

TEST_CASE("the homer reader accepts an uncompressed table") {
    TempFile file(".homer");
    write_text(file.path(),
               "HiCMatrix (directory=.)\tRegions\tchr1-0\tchr1-10\t\n"
               "chr1-0\tchr1-0\t1.5\t2.5\n"
               "chr1-10\tchr1-10\t2.5\t0.0");
    const hicx::MatrixData data = hicx::read_homer(file.path());
    REQUIRE(data.cut_intervals.size() == 2);
    CHECK(data.cut_intervals[1].start == 10);
    CHECK(data.cut_intervals[1].end == 20);
    CHECK(data.matrix.at(0, 0) == 1.5);
    CHECK(data.matrix.at(0, 1) == 2.5);
    CHECK(data.matrix.at(1, 0) == 2.5);
    // csr_matrix over a dense table stores no zero.
    CHECK(data.matrix.stored_nnz() == 3);
}

TEST_CASE("a homer bin name with a dash in the chromosome is rejected") {
    TempFile file(".homer");
    write_text(file.path(),
               "HiCMatrix (directory=.)\tRegions\tchr-1-0\tchr-1-10\t\n"
               "chr-1-0\tchr-1-0\t1\t2\n"
               "chr-1-10\tchr-1-10\t2\t3");
    // The Python unpacks value.split(b'-') into exactly two names and raises.
    CHECK_THROWS(hicx::read_homer(file.path()));
}

TEST_CASE("ginteractions writes the sidecar and leaves the named file alone") {
    TempFile file(".ginteractions");
    write_text(file.path(), "");
    hicx::write_ginteractions(file.path(), toy_matrix("int32"));
    CHECK(read_text(file.path()).empty());
    CHECK(read_text(file.path() + ".tsv") ==
          "chr1\t0\t10\tchr1\t0\t10\t5\n"
          "chr1\t0\t10\tchr1\t10\t20\t2\n"
          "chr1\t10\t20\tchr2\t0\t10\t3\n"
          "chr2\t0\t10\tchr2\t0\t10\t7\n");
}

TEST_CASE("ginteractions prints float values with the Python float repr") {
    TempFile file(".ginteractions");
    hicx::MatrixData data = toy_matrix("float64");
    data.matrix.mutable_data()[0] = 1.3e-05;
    hicx::write_ginteractions(file.path(), data);
    const std::string text = read_text(file.path() + ".tsv");
    CHECK(text.rfind("chr1\t0\t10\tchr1\t0\t10\t1.3e-05\n", 0) == 0);
}

TEST_CASE("hicpro writes one based ids for every stored entry") {
    TempFile matrix_file(".hicpro");
    TempFile bed_file(".bed");
    hicx::write_hicpro(matrix_file.path(), bed_file.path(), toy_matrix("int32"));
    CHECK(read_text(matrix_file.path()) == "1\t1\t5\n1\t2\t2\n2\t3\t3\n3\t3\t7\n");
    CHECK(read_text(bed_file.path()) ==
          "chr1\t0\t10\t1\nchr1\t10\t20\t2\nchr2\t0\t10\t3\n");

    const hicx::MatrixData back =
        hicx::read_hicpro(matrix_file.path(), bed_file.path());
    REQUIRE(back.cut_intervals.size() == 3);
    CHECK(back.cut_intervals[2] == hicx::CutInterval{"chr2", 0, 10, 3.0, ""});
    CHECK(back.matrix.at(0, 0) == 5.0);
    CHECK(back.matrix.at(1, 2) == 3.0);
    // The reader is float64 whatever the writer's dtype was.
    CHECK(back.matrix.dtype() == "float64");
}

TEST_CASE("the hicpro reader sums duplicate coordinates") {
    TempFile matrix_file(".hicpro");
    TempFile bed_file(".bed");
    write_text(bed_file.path(), "chr1\t0\t10\t1\nchr1\t10\t20\t2\n");
    write_text(matrix_file.path(), "1\t2\t1.5\n1\t2\t2.5\n");
    const hicx::MatrixData data =
        hicx::read_hicpro(matrix_file.path(), bed_file.path());
    CHECK(data.matrix.at(0, 1) == 4.0);
    CHECK(data.matrix.stored_nnz() == 1);
}

TEST_CASE("the chromosome sizes reader keeps file order") {
    TempFile file(".sizes");
    write_text(file.path(), "chr2\t200\nchr1\t100\n");
    const auto sizes = hicx::read_chromosome_sizes(file.path());
    REQUIRE(sizes.size() == 2);
    CHECK(sizes[0].first == "chr2");
    CHECK(sizes[0].second == 200);
    CHECK(sizes[1].first == "chr1");
}

TEST_CASE("2D text assigns, does not accumulate, and keeps the last value") {
    TempFile file(".txt");
    write_text(file.path(),
               "chr1\t0\t5\tchr1\t0\t5\t1.0\n"
               "chr1\t0\t5\tchr1\t0\t5\t3.0\n");
    const hicx::MatrixData data =
        hicx::read_two_dimensional_text(file.path(), {{"chr1", 30}}, 10);
    REQUIRE(data.matrix.rows() == 3);
    CHECK(data.matrix.at(0, 0) == 3.0);
    CHECK(data.matrix.stored_nnz() == 1);
}

TEST_CASE("2D text sets two cells when an interval straddles a bin boundary") {
    TempFile file(".txt");
    // The first interval ends at 10, which is the first position of bin 1, so
    // getRegionBinRange returns (0, 1) and the pair indexes two cells.
    write_text(file.path(), "chr1\t0\t10\tchr1\t20\t25\t2.0\n");
    const hicx::MatrixData data =
        hicx::read_two_dimensional_text(file.path(), {{"chr1", 30}}, 10);
    CHECK(data.matrix.at(0, 2) == 2.0);
    CHECK(data.matrix.at(1, 2) == 2.0);
    CHECK(data.matrix.stored_nnz() == 2);
}

TEST_CASE("2D text drops a line whose position is past the chromosome end") {
    TempFile file(".txt");
    write_text(file.path(),
               "chr1\t0\t5\tchr1\t100\t105\t2.0\n"
               "chr1\t0\t5\tchr1\t10\t15\t4.0\n");
    const hicx::MatrixData data =
        hicx::read_two_dimensional_text(file.path(), {{"chr1", 30}}, 10);
    CHECK(data.matrix.stored_nnz() == 1);
    CHECK(data.matrix.at(0, 1) == 4.0);
}

TEST_CASE("2D text generates a bin table that ends the chromosome exactly") {
    TempFile file(".txt");
    write_text(file.path(), "");
    const hicx::MatrixData data =
        hicx::read_two_dimensional_text(file.path(), {{"chr1", 25}, {"chr2", 10}}, 10);
    REQUIRE(data.cut_intervals.size() == 4);
    CHECK(data.cut_intervals[2] == hicx::CutInterval{"chr1", 20, 25, 1.0, ""});
    CHECK(data.cut_intervals[3] == hicx::CutInterval{"chr2", 0, 10, 1.0, ""});
}
