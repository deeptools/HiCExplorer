// Tests for the cool and h5 writers.
//
// Small hand written matrices pin the mechanics: which nodes appear, what the
// derived tables and attributes contain, that the upper triangle is selected
// from a full symmetric matrix without materialising anything, and that the
// correction factors are inverted and divided back out the way hicmatrix does
// it. The equivalence against a Python written file is checked on the real
// matrices by cpp/scripts/equiv.py, which is where contract rule 2 applies;
// here the point is the behaviour, not the corpus.

#include <doctest/doctest.h>

#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/matrix_data.hpp"

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

// A temporary file that removes itself, so a failing assertion cannot leave
// output behind in the build tree.
class TempFile {
  public:
    explicit TempFile(const std::string& suffix) {
        path_ = (std::filesystem::temp_directory_path() /
                 ("hicx-test-" + std::to_string(++counter_) + "-" +
                  std::to_string(::getpid()) + suffix))
                    .string();
    }
    ~TempFile() { std::remove(path_.c_str()); }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }

  private:
    static inline int counter_ = 0;
    std::string path_;
};

// Three bins on two chromosomes and the upper triangle of
//   5 2 0
//   2 0 3
//   0 3 7
hicx::MatrixData toy_matrix() {
    hicx::MatrixData data;
    const std::vector<std::int32_t> row{0, 0, 1, 2};
    const std::vector<std::int32_t> col{0, 1, 2, 2};
    std::vector<double> values{5.0, 2.0, 3.0, 7.0};
    data.matrix = hicx::CsrMatrix::from_coo(3, 3, row, col, std::move(values), "int32");
    data.cut_intervals = {
        hicx::CutInterval{"chr1", 0, 10, 1.0, ""},
        hicx::CutInterval{"chr1", 10, 20, 1.0, ""},
        hicx::CutInterval{"chr2", 0, 10, 1.0, ""},
    };
    return data;
}

// cut_intervals compare with == on the extra value, which is NaN for some
// matrices in the corpus, so a round trip needs an equality that treats two
// NaNs as equal.
bool same_bins(const std::vector<hicx::CutInterval>& a,
               const std::vector<hicx::CutInterval>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        const bool extra_equal = a[i].extra == b[i].extra ||
                                 (std::isnan(a[i].extra) && std::isnan(b[i].extra));
        if (a[i].chrom != b[i].chrom || a[i].start != b[i].start ||
            a[i].end != b[i].end || !extra_equal ||
            a[i].extra_text != b[i].extra_text) {
            return false;
        }
    }
    return true;
}

std::string attribute_string(const hicx::h5::File& file, const std::string& name) {
    const auto attributes = file.attributes("/");
    const auto found = attributes.find(name);
    REQUIRE(found != attributes.end());
    return std::get<std::string>(found->second);
}

std::int64_t attribute_int(const hicx::h5::File& file, const std::string& name) {
    const auto attributes = file.attributes("/");
    const auto found = attributes.find(name);
    REQUIRE(found != attributes.end());
    return std::get<std::int64_t>(found->second);
}

}  // namespace

TEST_CASE("guess_chunk reproduces h5py's automatic chunk layout") {
    // Measured on the files cooler writes through hicmatrix, see the dataset
    // listings in cpp/PLAN.md 2.4: a 33,754 bin int32 column gets 2,110, the
    // 168,770 pixel int32 columns get 3,470 and the float64 count column of
    // the 11,104 bin matrix gets 1,735.
    CHECK(hicx::h5::guess_chunk(33754, 4) == 2110);
    CHECK(hicx::h5::guess_chunk(33755, 8) == 1055);
    CHECK(hicx::h5::guess_chunk(168770, 4) == 5275);
    CHECK(hicx::h5::guess_chunk(55520, 4) == 3470);
    CHECK(hicx::h5::guess_chunk(55520, 8) == 1735);
    CHECK(hicx::h5::guess_chunk(11104, 8) == 1388);
    // Small datasets fit into a single chunk.
    CHECK(hicx::h5::guess_chunk(15, 9) == 15);
    CHECK(hicx::h5::guess_chunk(1, 4) == 1);
}

TEST_CASE("the cool writer produces the cooler object tree") {
    hicx::MatrixData data = toy_matrix();
    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.creation_date = "2026-09-01T00:00:00.000000";
    hicx::write_cool(output.path(), data, options);

    const hicx::h5::File file(output.path());
    for (const char* path : {"/chroms/name", "/chroms/length", "/bins/chrom",
                             "/bins/start", "/bins/end", "/pixels/bin1_id",
                             "/pixels/bin2_id", "/pixels/count",
                             "/indexes/chrom_offset", "/indexes/bin1_offset"}) {
        CHECK_MESSAGE(file.exists(path), path);
    }
    CHECK_FALSE(file.exists("/bins/weight"));

    CHECK(file.read_strings("/chroms/name") == std::vector<std::string>{"chr1", "chr2"});
    // The chromosome length is the end of its last bin.
    CHECK(file.read_int64("/chroms/length") == std::vector<std::int64_t>{20, 10});
    CHECK(file.read_int32("/bins/chrom") == std::vector<std::int32_t>{0, 0, 1});
    CHECK(file.read_int64("/bins/start") == std::vector<std::int64_t>{0, 10, 0});
    CHECK(file.read_int64("/indexes/chrom_offset") == std::vector<std::int64_t>{0, 2, 3});
    CHECK(file.read_int64("/indexes/bin1_offset") ==
          std::vector<std::int64_t>{0, 2, 3, 4});

    CHECK(file.read_int32("/pixels/bin1_id") == std::vector<std::int32_t>{0, 0, 1, 2});
    CHECK(file.read_int32("/pixels/bin2_id") == std::vector<std::int32_t>{0, 1, 2, 2});
    CHECK(file.read_int64("/pixels/count") == std::vector<std::int64_t>{5, 2, 3, 7});
    CHECK(file.dataset_dtype("/pixels/count") == "int32");
    CHECK(file.dataset_dtype("/pixels/bin1_id") == "int32");
    CHECK(file.dataset_dtype("/indexes/bin1_offset") == "int64");

    CHECK(attribute_int(file, "nbins") == 3);
    CHECK(attribute_int(file, "nchroms") == 2);
    CHECK(attribute_int(file, "nnz") == 4);
    CHECK(attribute_int(file, "sum") == 17);
    CHECK(attribute_int(file, "format-version") == 3);
    CHECK(attribute_int(file, "bin-size") == 10);
    CHECK(attribute_string(file, "bin-type") == "fixed");
    CHECK(attribute_string(file, "storage-mode") == "symmetric-upper");
    CHECK(attribute_string(file, "format") == "HDF5::Cooler");
    CHECK(attribute_string(file, "generated-by") == "HiCMatrix-17.2");
    CHECK(attribute_string(file, "genome-assembly") == "unknown");
    CHECK(attribute_string(file, "creation-date") == "2026-09-01T00:00:00.000000");
    CHECK(attribute_string(file, "metadata") ==
          "{\"format\": \"HDF5::Cooler\", \"format-url\": "
          "\"https://github.com/mirnylab/cooler\", \"generated-by\": "
          "\"HiCMatrix-17.2\", \"generated-by-cooler-lib\": \"cooler-0.10.2\", "
          "\"tool-url\": \"https://github.com/deeptools/HiCMatrix\"}");
}

TEST_CASE("a variable bin size is written as the string null") {
    // cooler infers the bin size from every bin except the last of each
    // chromosome, so two different widths among those make it variable.
    hicx::MatrixData data = toy_matrix();
    data.cut_intervals = {
        hicx::CutInterval{"chr1", 0, 10, 1.0, ""},
        hicx::CutInterval{"chr1", 10, 17, 1.0, ""},
        hicx::CutInterval{"chr1", 17, 30, 1.0, ""},
    };
    const TempFile output(".cool");
    hicx::write_cool(output.path(), data);

    const hicx::h5::File file(output.path());
    CHECK(attribute_string(file, "bin-type") == "variable");
    CHECK(attribute_string(file, "bin-size") == "null");
}

TEST_CASE("only the upper triangle of a symmetric matrix is written") {
    hicx::MatrixData data = toy_matrix();
    data.matrix.materialize_full();
    REQUIRE(data.matrix.symmetry() == hicx::Symmetry::Full);
    REQUIRE(data.matrix.stored_nnz() == 6);

    const TempFile output(".cool");
    hicx::write_cool(output.path(), data);
    const hicx::h5::File file(output.path());
    CHECK(file.read_int32("/pixels/bin1_id") == std::vector<std::int32_t>{0, 0, 1, 2});
    CHECK(file.read_int32("/pixels/bin2_id") == std::vector<std::int32_t>{0, 1, 2, 2});
    CHECK(attribute_int(file, "nnz") == 4);
}

TEST_CASE("correction factors are inverted and divided back out for an h5 source") {
    hicx::MatrixData data = toy_matrix();
    data.matrix.set_dtype("float64");
    data.correction_factors = std::vector<double>{2.0, 4.0, 0.0};

    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.file_was_h5 = true;
    hicx::write_cool(output.path(), data, options);

    const hicx::h5::File file(output.path());
    // 1/2, 1/4 and 1/0, the last of which becomes 0 rather than infinity.
    const std::vector<double> weight = file.read_doubles("/bins/weight");
    REQUIRE(weight.size() == 3);
    CHECK(weight[0] == doctest::Approx(0.5));
    CHECK(weight[1] == doctest::Approx(0.25));
    CHECK(weight[2] == 0.0);

    // The counts are divided by the product of the inverted factors, so
    // 5 / (0.5 * 0.5) = 20 and 2 / (0.5 * 0.25) = 16.
    const std::vector<double> counts = file.read_doubles("/pixels/count");
    REQUIRE(counts.size() == 4);
    CHECK(counts[0] == doctest::Approx(20.0));
    CHECK(counts[1] == doctest::Approx(16.0));
    // A zero factor turns its entries into infinity, which cool stores as it
    // is. That is what hicmatrix does and the port does not improve on it.
    CHECK(std::isinf(counts[2]));
    CHECK(file.dataset_dtype("/pixels/count") == "float64");
}

TEST_CASE("a nan bin pair is dropped when the input was an h5 file") {
    hicx::MatrixData data = toy_matrix();
    // Bins 1 and 2 are NaN bins, so both (1, 2) and the diagonal entry (2, 2)
    // have both of their ends in the set and are removed, while (0, 1) keeps
    // one good end and survives.
    data.nan_bins = {1, 2};
    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.file_was_h5 = true;
    hicx::write_cool(output.path(), data, options);

    const hicx::h5::File file(output.path());
    CHECK(file.read_int32("/pixels/bin1_id") == std::vector<std::int32_t>{0, 0});
    CHECK(file.read_int32("/pixels/bin2_id") == std::vector<std::int32_t>{0, 1});
    CHECK(attribute_int(file, "nnz") == 2);
}

TEST_CASE("without pApplyCorrection the counts keep their correction") {
    hicx::MatrixData data = toy_matrix();
    data.matrix.set_dtype("float64");
    data.correction_factors = std::vector<double>{2.0, 4.0, 8.0};

    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.apply_correction = false;
    options.file_was_h5 = true;
    hicx::write_cool(output.path(), data, options);

    const hicx::h5::File file(output.path());
    // The factors are neither inverted nor divided out; they are only stored.
    CHECK(file.read_doubles("/bins/weight") == std::vector<double>{2.0, 4.0, 8.0});
    CHECK(file.read_int64("/pixels/count") == std::vector<std::int64_t>{5, 2, 3, 7});
}

TEST_CASE("a cooler group URI writes into that group and leaves the root alone") {
    const TempFile output(".mcool");
    {
        hicx::MatrixData data = toy_matrix();
        hicx::write_cool(output.path() + "::/resolutions/10000", data);
    }
    {
        // The second resolution appends, so the first one survives and the
        // root provenance is not written again.
        hicx::MatrixData data = toy_matrix();
        hicx::CoolSaveOptions options;
        options.append = true;
        options.generated_by = "should-not-reach-the-root";
        hicx::write_cool(output.path() + "::/resolutions/20000", data, options);
    }

    const hicx::h5::File file(output.path());
    CHECK(file.exists("/resolutions/10000/pixels/count"));
    CHECK(file.exists("/resolutions/20000/pixels/count"));
    CHECK(file.children("/") == std::vector<std::string>{"resolutions"});

    // cooler's own provenance stays on the group; hicmatrix overwrites the
    // root, and only in mode 'w'.
    const auto group_attrs = file.attributes("/resolutions/10000");
    CHECK(std::get<std::string>(group_attrs.at("generated-by")) == "cooler-0.10.2");
    CHECK(std::get<std::string>(group_attrs.at("format-url")) ==
          "https://github.com/open2c/cooler");
    const auto root_attrs = file.attributes("/");
    CHECK(std::get<std::string>(root_attrs.at("generated-by")) == "HiCMatrix-17.2");
    CHECK(std::get<std::string>(root_attrs.at("format-url")) ==
          "https://github.com/mirnylab/cooler");
    CHECK(root_attrs.count("nbins") == 0);
}

TEST_CASE("enforce_integer rounds half to even into an int32 count column") {
    hicx::MatrixData data = toy_matrix();
    data.matrix.set_dtype("float64");
    data.matrix.mutable_data() = {0.5, 1.5, 2.5, -0.5};

    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.enforce_integer = true;
    hicx::write_cool(output.path(), data, options);

    const hicx::h5::File file(output.path());
    CHECK(file.dataset_dtype("/pixels/count") == "int32");
    CHECK(file.read_int64("/pixels/count") == std::vector<std::int64_t>{0, 2, 2, 0});
    // The four entries are still stored even though two of them round to zero,
    // because hicmatrix does not eliminate zeros after np.rint.
    CHECK(attribute_int(file, "nnz") == 4);
}

TEST_CASE("the h5 writer produces the HiCExplorer node set") {
    hicx::MatrixData data = toy_matrix();
    data.nan_bins = {2};
    data.correction_factors = std::vector<double>{1.0, 2.0, 3.0};

    const TempFile output(".h5");
    hicx::write_hicexplorer_h5(output.path(), data);

    const hicx::h5::File file(output.path());
    for (const char* path : {"/matrix/data", "/matrix/indices", "/matrix/indptr",
                             "/matrix/shape", "/intervals/chr_list",
                             "/intervals/start_list", "/intervals/end_list",
                             "/intervals/extra_list", "/nan_bins",
                             "/correction_factors"}) {
        CHECK_MESSAGE(file.exists(path), path);
    }
    CHECK_FALSE(file.exists("/distance_counts"));

    CHECK(file.read_int64("/matrix/shape") == std::vector<std::int64_t>{3, 3});
    CHECK(file.read_int64("/matrix/indptr") == std::vector<std::int64_t>{0, 2, 3, 4});
    CHECK(file.read_int32("/matrix/indices") == std::vector<std::int32_t>{0, 1, 2, 2});
    CHECK(file.read_int64("/matrix/data") == std::vector<std::int64_t>{5, 2, 3, 7});
    CHECK(file.dataset_dtype("/matrix/data") == "int32");
    CHECK(file.dataset_dtype("/matrix/indices") == "int32");
    CHECK(file.read_strings("/intervals/chr_list") ==
          std::vector<std::string>{"chr1", "chr1", "chr2"});
    CHECK(file.read_int64("/intervals/end_list") == std::vector<std::int64_t>{10, 20, 10});
    CHECK(file.read_doubles("/intervals/extra_list") ==
          std::vector<double>{1.0, 1.0, 1.0});
    CHECK(file.read_int64("/nan_bins") == std::vector<std::int64_t>{2});
}

TEST_CASE("a NaN correction factor is written as zero") {
    hicx::MatrixData data = toy_matrix();
    data.correction_factors = std::vector<double>{
        1.0, std::numeric_limits<double>::quiet_NaN(), 3.0};
    const TempFile output(".h5");
    hicx::write_hicexplorer_h5(output.path(), data);

    const hicx::h5::File file(output.path());
    CHECK(file.read_doubles("/correction_factors") ==
          std::vector<double>{1.0, 0.0, 3.0});
}

TEST_CASE("the h5 writer keeps a blosc compressed file readable by the reader") {
    // The round trip through the compression half of the filter, which is what
    // open question 1 of cpp/STATUS.md turned on.
    const hicx::H5MatrixData source = hicx::read_hicexplorer_h5(
        data_path("small_test_matrix.h5"));
    const TempFile output(".h5");
    hicx::MatrixData copy = source;
    hicx::write_hicexplorer_h5(output.path(), copy);

    const hicx::H5MatrixData written = hicx::read_hicexplorer_h5(output.path());
    CHECK(written.matrix.rows() == source.matrix.rows());
    CHECK(written.matrix.stored_nnz() == source.matrix.stored_nnz());
    CHECK(written.matrix.data() == source.matrix.data());
    CHECK(written.matrix.indices() == source.matrix.indices());
    CHECK(written.matrix.indptr() == source.matrix.indptr());
    CHECK(same_bins(written.cut_intervals, source.cut_intervals));
    CHECK(written.matrix.dtype() == source.matrix.dtype());
}

TEST_CASE("a cool file written from a real matrix reads back unchanged") {
    hicx::CoolLoadResult loaded = hicx::read_cool(data_path("small_test_matrix.cool"));
    const std::size_t nnz = loaded.data.matrix.upper_triangle_nnz();
    const std::vector<hicx::CutInterval> bins = loaded.data.cut_intervals;

    const TempFile output(".cool");
    hicx::CoolSaveOptions options;
    options.has_hic_metadata = true;
    options.hic_metadata = loaded.metadata;
    hicx::write_cool(output.path(), loaded.data, options);

    const hicx::CoolLoadResult again = hicx::read_cool(output.path());
    CHECK(again.data.matrix.upper_triangle_nnz() == nnz);
    CHECK(same_bins(again.data.cut_intervals, bins));
    CHECK(again.data.matrix.data() == loaded.data.matrix.data());
    CHECK(again.data.matrix.indices() == loaded.data.matrix.indices());
}

namespace {

// Identifiers still open on any HDF5 file, with their names, for the message
// of a failing check. An identifier left open keeps its file open until the
// library shuts down, so a tool that execs the drawing process afterwards
// would leave an unflushed, unreadable file behind.
std::string open_hdf5_objects() {
    const ssize_t count = H5Fget_obj_count(static_cast<hid_t>(H5F_OBJ_ALL), H5F_OBJ_ALL);
    std::string text = std::to_string(count);
    if (count <= 0) {
        return text;
    }
    std::vector<hid_t> ids(static_cast<std::size_t>(count));
    const ssize_t listed =
        H5Fget_obj_ids(static_cast<hid_t>(H5F_OBJ_ALL), H5F_OBJ_ALL, ids.size(), ids.data());
    for (ssize_t i = 0; i < listed; ++i) {
        char name[512] = {0};
        H5Iget_name(ids[static_cast<std::size_t>(i)], name, sizeof(name));
        text += " [type " + std::to_string(H5Iget_type(ids[static_cast<std::size_t>(i)])) +
                " " + name + "]";
    }
    return text;
}

}  // namespace

TEST_CASE("the h5 and cool writers leave no HDF5 identifier open") {
    // Compared against the count before the write, so an identifier another
    // test left open does not fail this one, and without REQUIRE, which would
    // abort doctest's re-entry for the second subcase.
    const std::string before = open_hdf5_objects();
    SUBCASE("h5") {
        hicx::MatrixData data = toy_matrix();
        data.correction_factors = std::vector<double>{1.0, 2.0, 4.0};
        const TempFile output(".h5");
        hicx::H5SaveOptions options;
        options.symmetric = true;
        hicx::write_hicexplorer_h5(output.path(), data, options);
        CHECK(open_hdf5_objects() == before);
    }
    SUBCASE("cool") {
        hicx::MatrixData data = toy_matrix();
        const TempFile output(".cool");
        hicx::CoolSaveOptions options;
        options.symmetric = true;
        hicx::write_cool(output.path(), data, options);
        CHECK(open_hdf5_objects() == before);
    }
}
