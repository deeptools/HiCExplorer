// Reader for the cool format (HDF5, https://github.com/open2c/cooler).
//
// A cooler is an HDF5 group holding
//   chroms/name, chroms/length
//   bins/chrom, bins/start, bins/end and optional weight columns
//   pixels/bin1_id, pixels/bin2_id, pixels/count
//   indexes/bin1_offset, indexes/chrom_offset
// plus file level attributes: nbins, nchroms, nnz, sum, bin-size, bin-type,
// storage-mode, generated-by, metadata and others.
//
// Only reading is implemented in this milestone.

#ifndef HICX_COOL_FILE_HPP
#define HICX_COOL_FILE_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/hdf5_util.hpp"
#include "hicx/json_lite.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// cooler.fileops.is_cooler: an HDF5 group carrying the four cooler groups.
[[nodiscard]] bool is_cooler(const std::string& path);

// hicexplorer.utilities.check_cooler.
[[nodiscard]] bool check_cooler(const std::string& path);

class CoolFile {
  public:
    // Accepts a plain path or a cooler URI "file.mcool::/resolutions/10000".
    explicit CoolFile(const std::string& uri);

    [[nodiscard]] const std::string& filename() const noexcept { return filename_; }
    [[nodiscard]] const std::string& root() const noexcept { return root_; }

    // Cooler.info: file attributes with cooler's json coercion applied.
    [[nodiscard]] const std::map<std::string, json::Value>& info() const noexcept {
        return info_;
    }
    [[nodiscard]] const json::Value* info_value(const std::string& key) const;

    [[nodiscard]] const std::vector<std::string>& chrom_names() const noexcept {
        return chrom_names_;
    }
    [[nodiscard]] const std::vector<std::int64_t>& chrom_lengths() const noexcept {
        return chrom_lengths_;
    }

    // Cooler.bins().columns.values: chrom, start, end first, then any further
    // column of the bins group in HDF5 name order.
    [[nodiscard]] std::vector<std::string> bin_columns() const;

    [[nodiscard]] std::int64_t nbins() const;
    [[nodiscard]] std::int64_t nnz() const;

    // The bin table as cut intervals, with the constant extra value 1.0 that
    // hicmatrix uses for cool files.
    [[nodiscard]] std::vector<CutInterval> read_bins() const;

    [[nodiscard]] bool has_column(const std::string& column) const;
    [[nodiscard]] std::vector<double> read_column(const std::string& column) const;

    // The raw pixel table as a CSR matrix, without any balancing applied.
    [[nodiscard]] CsrMatrix read_matrix() const;

  private:
    std::string filename_;
    std::string root_ = "/";
    std::shared_ptr<h5::File> file_;
    std::map<std::string, json::Value> info_;
    std::vector<std::string> chrom_names_;
    std::vector<std::int64_t> chrom_lengths_;

    [[nodiscard]] std::string path_of(const std::string& relative) const;
};

}  // namespace hicx

#endif  // HICX_COOL_FILE_HPP
