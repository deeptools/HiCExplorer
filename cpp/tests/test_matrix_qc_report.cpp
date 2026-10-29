// Tests of the hicMatrixQC table and HTML report.

#include <doctest/doctest.h>

#include <string>
#include <tuple>
#include <vector>

#include "hicx/matrix_qc.hpp"
#include "hicx/matrix_qc_report.hpp"

namespace qc = hicx::matrixqc;

namespace {

using Pixel = std::tuple<std::int64_t, std::int64_t, double>;

qc::PixelSource source_of(const std::vector<Pixel>& pixels) {
    return [pixels](unsigned, const qc::ChunkVisitor& visit) {
        hicx::PixelChunk chunk;
        for (const auto& [a, b, c] : pixels) {
            chunk.bin1.push_back(a);
            chunk.bin2.push_back(b);
            chunk.count.push_back(c);
        }
        visit(0, chunk);
    };
}

qc::Options small_options() {
    qc::Options options;
    options.threads = 1;
    options.rao_min_contacts = 7;
    options.fit_min = 20000;
    options.fit_max = 80000;
    return options;
}

qc::Layout make_layout() {
    qc::Layout layout;
    layout.resolution = 10000;
    layout.names = {"chrA", "chrB"};
    layout.lengths = {100000, 60000};
    layout.first_bin = {0, 10, 16};
    return layout;
}

const std::vector<Pixel> kTwoChromosomes = {{0, 0, 5}, {0, 1, 3}, {1, 1, 4}, {2, 5, 2},
                                            {9, 9, 1}, {3, 12, 7}, {10, 10, 6}, {10, 12, 2}};

}  // namespace

TEST_CASE("the statistics table has one column per sample") {
    const qc::Layout layout = make_layout();
    qc::FileStats file;
    file.path = "synthetic.cool";
    qc::MatrixStats a =
        qc::compute_matrix_stats(layout, source_of(kTwoChromosomes), std::nullopt, small_options());
    a.label = "first";
    qc::MatrixStats b = a;
    b.label = "second";
    file.matrices = {a, b};
    const std::string table = qc::render_table({file});
    CHECK(table.rfind("statistic\tfirst\tsecond\n", 0) == 0);
    CHECK(table.find("\ncontacts.total\t30\t30\n") != std::string::npos);
    CHECK(table.find("\ncontacts.cis_fraction\t") != std::string::npos);
    const std::string html = qc::render_html({file}, small_options(), "test");
    CHECK(html.find("<script") == std::string::npos);
    CHECK(html.find("<svg") == std::string::npos);
    CHECK(html.find("<style><!--") != std::string::npos);
    CHECK(html.find("bgcolor=") != std::string::npos);
    CHECK(html.find("\xE2\x80\x94") == std::string::npos);
}

