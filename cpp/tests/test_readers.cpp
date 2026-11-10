// Reader tests against the real matrices of the Python test suite. The
// expected values come from the Python implementation, see
// hicexplorer/test/general/test_hicInfo.py.

#include <doctest/doctest.h>

#include <string>

#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/hic_matrix.hpp"
#include "hicx/numpy_compat.hpp"

namespace {

std::string data_path(const std::string& relative) {
    return std::string(HICX_TEST_DATA_DIR) + "/" + relative;
}

}  // namespace

TEST_CASE("cool metadata of Li_et_al_2015.cool") {
    const hicx::CoolFile cool(data_path("Li_et_al_2015.cool"));
    CHECK(cool.nbins() == 11104);
    CHECK(cool.nnz() == 1661678);
    CHECK(cool.chrom_names() == std::vector<std::string>{"X"});
    CHECK(cool.chrom_lengths() == std::vector<std::int64_t>{22422827});
    CHECK(cool.bin_columns() ==
          std::vector<std::string>{"chrom", "start", "end"});

    // bin-size is the literal string "null" in this file because the bins are
    // restriction fragments, so cooler reports None and hicInfo prints nothing.
    const hicx::json::Value* bin_size = cool.info_value("bin-size");
    REQUIRE(bin_size != nullptr);
    CHECK(bin_size->is_null());
    CHECK(cool.info_value("bin-type")->as_string() == "variable");
    CHECK(cool.info_value("generated-by")->as_string() == "cooler-0.7.11");
    CHECK(cool.info_value("creation-date")->as_string() == "2018-11-08T20:03:12.737522");
    // Older coolers carry no storage-mode attribute; the upper triangle layout
    // is then implied.
    CHECK(cool.info_value("storage-mode") == nullptr);
}

TEST_CASE("cool metadata of a 100 kb matrix") {
    const hicx::CoolFile cool(
        data_path("hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.cool"));
    CHECK(cool.nbins() == 1972);
    CHECK(cool.nnz() == 1252980);
    CHECK(cool.info_value("bin-size")->as_int() == 100000);
    CHECK(cool.info_value("storage-mode")->as_string() == "symmetric-upper");
    CHECK(cool.chrom_names() == std::vector<std::string>{"chr1"});

    const std::vector<hicx::CutInterval> bins = cool.read_bins();
    REQUIRE(bins.size() == 1972);
    CHECK(bins.front().chrom == "chr1");
    CHECK(bins.front().start == 0);
    CHECK(bins.front().end == 100000);
    CHECK(bins.back().end == 197195432);
}

TEST_CASE("h5 reader on Li_et_al_2015.h5") {
    const hicx::H5MatrixData data = hicx::read_hicexplorer_h5(
        data_path("Li_et_al_2015.h5"));
    CHECK(data.matrix.rows() == 11104);
    CHECK(data.matrix.cols() == 11104);
    CHECK(data.matrix.nnz() == 1661678);  // upper triangle as stored
    CHECK(data.matrix.dtype() == "float64");
    CHECK(data.nan_bins.size() == 855);
    REQUIRE(data.correction_factors.has_value());
    CHECK(data.correction_factors->size() == 11104);
    REQUIRE(data.cut_intervals.size() == 11104);
    CHECK(data.cut_intervals.front().chrom == "X");
    CHECK(data.cut_intervals.back().end == 22422827);
}

TEST_CASE("HiCMatrix reproduces the numbers hicInfo prints for the h5 file") {
    const hicx::HiCMatrix hic = hicx::HiCMatrix::load(data_path("Li_et_al_2015.h5"));
    CHECK(hic.matrix().rows() == 11104);
    CHECK(hic.matrix().nnz() == 3313107);
    CHECK(hic.bin_size() == 1843);
    CHECK(hic.nan_bins().size() == 855);
    CHECK(hicx::npy::float_repr(hic.matrix().data_min().as_double()) ==
          "0.1702102389180371");
    CHECK(hicx::npy::float_repr(hic.matrix().data_max().as_double()) ==
          "1914.0146478655674");

    const hicx::Scalar total = hic.matrix().sum();
    const hicx::Scalar diagonal = hic.matrix().diagonal_sum();
    const double sum_elements =
        ((total.as_double() - diagonal.as_double()) / 2.0) + diagonal.as_double();
    CHECK(hicx::npy::float_repr(sum_elements) == "17548966.536917936");

    const auto sizes = hic.chromosome_sizes();
    REQUIRE(sizes.size() == 1);
    CHECK(sizes[0].first == "X");
    CHECK(sizes[0].second == 22422827);
}

TEST_CASE("the cool and the h5 version of Li_et_al_2015 agree") {
    const hicx::HiCMatrix from_cool =
        hicx::HiCMatrix::load(data_path("Li_et_al_2015.cool"));
    CHECK(from_cool.matrix().nnz() == 3313107);
    CHECK(from_cool.bin_size() == 1843);
    CHECK(from_cool.nan_bins().size() == 855);
    const hicx::Scalar total = from_cool.matrix().sum();
    const hicx::Scalar diagonal = from_cool.matrix().diagonal_sum();
    const double sum_elements =
        ((total.as_double() - diagonal.as_double()) / 2.0) + diagonal.as_double();
    CHECK(hicx::npy::float_repr(sum_elements) == "17548966.536917936");
}

TEST_CASE("check_cooler follows hicexplorer.utilities") {
    CHECK(hicx::check_cooler(data_path("Li_et_al_2015.cool")));
    CHECK_FALSE(hicx::check_cooler(data_path("Li_et_al_2015.h5")));
    CHECK(hicx::is_cooler(data_path("Li_et_al_2015.cool")));
    CHECK_FALSE(hicx::is_cooler(data_path("Li_et_al_2015.h5")));
}
