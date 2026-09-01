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
// The writer reproduces what cooler.create_cooler produces when hicmatrix
// calls it (hicmatrix/lib/cool.py:406-426), including the dataset layout,
// because cool output is compared structurally (class E1 of cpp/PLAN.md).

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
#include "hicx/matrix_data.hpp"
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

struct CoolLoadOptions {
    // Cool.applyCorrectionLoad.
    bool apply_correction = true;
    // Cool.correctionFactorTable, the bin column holding the weights.
    std::string correction_factor_table = "weight";
    // Cool.correctionOperator when the caller sets it before the load, which
    // is what hicConvertFormat --correction_division does. The whole block
    // that derives the operator from the column name is then skipped
    // (cool.py:195-207), so the hic2cool and hicmatrix versions are not read
    // out of 'generated-by' either.
    std::optional<char> correction_operator;
    // Cool.chrnameList holding exactly one name (cool.py:120-131, :156-161):
    // only that chromosome's bins, pixels and weights are loaded, and the NaN
    // bin heuristic runs on that block alone rather than on the whole matrix.
    // hicAdjustMatrix takes this path for a single --chromosomes with
    // --action keep on a cool input (hicAdjustMatrix.py:73-76), and it is
    // observable: the block has more empty rows on its own than it has inside
    // the whole matrix, so the NaN bin list is not the same one a whole file
    // load followed by a selection would produce.
    std::optional<std::string> chrom_name;
};

struct CoolLoadResult {
    MatrixData data;
    // The state the loader leaves behind on the Cool object and that a
    // following save reads back (hicmatrix/lib/cool.py:195-207): the operator
    // the correction was applied with, and the version of whatever wrote the
    // file. hicConvertFormat carries both from the input handler to the output
    // handler, so they belong to the load result and not to the matrix.
    std::optional<char> correction_operator;
    std::optional<std::string> hic2cool_version;
    std::optional<std::string> hicmatrix_version;
    // Cooler.info with every value rendered the way Python's str() would.
    std::map<std::string, std::string> metadata;
};

// Port of hicmatrix.lib.Cool.load for the whole matrix case, which is the only
// one the tools use when no chromosome is preselected.
[[nodiscard]] CoolLoadResult read_cool(const std::string& uri,
                                       const CoolLoadOptions& options = CoolLoadOptions());

// The state hicmatrix.lib.Cool carries into save(), one field per Python
// attribute that changes what is written.
struct CoolSaveOptions {
    // pSymmetric: store triu(matrix, k=0). Always true in HiCExplorer.
    bool symmetric = true;
    // pApplyCorrection: divide the correction factors back out of the counts
    // so that 'count' on disk is raw and 'weight' carries the correction.
    bool apply_correction = true;
    // Cool.enforceInteger: round the counts to int32 with numpy's
    // round-half-to-even.
    bool enforce_integer = false;
    // Cool.appendData: cooler.create_cooler is called with mode 'a' instead of
    // 'w'. hicConvertFormat sets it for every resolution of an mcool after the
    // first, and it is also what decides whether the hicmatrix provenance
    // attributes are written onto the file root (cool.py:422-426).
    bool append = false;
    // Cool.fileWasH5: the input was an h5 file, which both triggers the
    // nan-bin masking and forces the correction factors to be inverted.
    bool file_was_h5 = false;
    // Cool.hic2cool_version, compared as a string against "0.5".
    std::optional<std::string> hic2cool_version;
    // Cool.correctionOperator, '*' or '/', set by the loader.
    std::optional<char> correction_operator;
    // Cool.hic_metadata, the info dictionary of the source cooler. Only
    // 'genome-assembly', 'matrix-generated-by' and 'matrix-generated-by-url'
    // are read out of it.
    std::map<std::string, std::string> hic_metadata;
    bool has_hic_metadata = false;
    // The provenance strings. cpp/PLAN.md 2.4 decides that the port emits the
    // hicmatrix identity verbatim until the whole suite is green, so that the
    // four normalised attributes are not a free pass in the comparator.
    std::string generated_by = "HiCMatrix-17.2";
    std::string generated_by_cooler_lib = "cooler-0.10.2";
    std::string tool_url = "https://github.com/deeptools/HiCMatrix";
    std::string format_url = "https://github.com/mirnylab/cooler";
    // cooler's own format-url, which survives on a cooler written into a group
    // because hicmatrix overwrites the provenance attributes on the file root
    // only. Note that hicmatrix's format_url above is the old mirnylab one.
    std::string cooler_format_url = "https://github.com/open2c/cooler";
    // ISO 8601 local time, like datetime.now().isoformat(). Empty means "take
    // the current time"; a fixed value makes a test reproducible.
    std::string creation_date;
};

// Port of hicmatrix.lib.Cool.create_cooler_input and Cool.save.
//
// `path` is a cooler URI: either a plain file name, in which case the cooler
// occupies the file root, or "file::/group/path", in which case it is written
// into that group and the file keeps whatever else it holds. The second form
// is how hicConvertFormat produces an mcool, one group per resolution. A
// cooler written into a group keeps cooler's own provenance attributes
// (format-url open2c, generated-by cooler-<version>) because hicmatrix only
// overwrites them on the file root, and it overwrites them there only in mode
// 'w', that is for the first resolution.
//
// `data` is modified in place exactly where the Python modifies it: NaN counts
// become zero, the pairs of NaN bins are dropped for an h5 source, the
// correction factors are inverted and the counts are divided by them. With
// pSymmetric the reverted counts are the upper triangle only, because at that
// point the Python has already replaced its matrix with that triangle. Nothing
// is copied; the pixel table is streamed out of the CSR one block at a time,
// so the writer adds a bounded buffer and the O(nbins) row offset array to the
// resident set and nothing that scales with the pixel count.
void write_cool(const std::string& path, MatrixData& data,
                const CoolSaveOptions& options = CoolSaveOptions());

}  // namespace hicx

#endif  // HICX_COOL_FILE_HPP
