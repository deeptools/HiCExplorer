// cool and mcool files for HiCExplorer v4: the hicmatrix layer over coolercpp.
//
// Every read and write of the cool format goes through coolercpp, the
// cooler-compatible library (cmake/HicxCoolercpp.cmake). This adapter keeps
// only what hicmatrix.lib.Cool adds on top of cooler, and hands the tools the
// v4 matrix types:
//
//   * the correction handling of Cool.load: weights applied on load,
//     multiplicative for 'weight' and divisive for the hic2cool tables, the
//     operator and the hic2cool/hicmatrix versions read out of generated-by;
//   * the NaN bin list rebuilt on every load as the bins with an empty row
//     and an empty column;
//   * the single chromosome load, which reads only that chromosome's block;
//   * Cool.create_cooler_input and Cool.save on write: NaN bin pair masking
//     for an h5 source, the inversion and reversion of the correction
//     factors, the weight column, int32 bin IDs, the count dtype, the pixel
//     table split into 10,000 parts above 10^7 pixels (which fixes the order
//     the 'sum' attribute is accumulated in), hicmatrix's metadata dictionary
//     and its key order, and hicmatrix's provenance written onto the file root
//     only in mode 'w'.
//
// The load time correction_factors/distance_counts swap of hiCMatrix
// (STATUS.md F9, F10, F21) stays with hicx::ToolMatrix, as before.

#ifndef HICX_COOL_ADAPTER_HPP
#define HICX_COOL_ADAPTER_HPP

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/json_lite.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace coolercpp {
class Cooler;
}

namespace hicx {

// The cooler of a path or URI, as cooler.fileops.is_cooler decides; a missing
// group is false rather than an error.
[[nodiscard]] bool is_cooler(const std::string& path);

// hicexplorer.utilities.check_cooler.
[[nodiscard]] bool check_cooler(const std::string& path);

// Sorts the columns of every row of `arrays` that is stored out of order,
// stably (so duplicate columns keep their relative order), in place. The
// cooler format requires each row's pixels to be stored in ascending column
// order, but a real file need not honour it (cpp/STATUS.md 2026-10-01: a 4D
// Nucleome GM12878 map has out-of-order rows); every reader that assembles a
// CsrMatrix::Arrays from pixel rows read directly off disk, rather than
// through CsrMatrix::from_coo (which always sorts), must call this before
// CsrMatrix::adopt, because CsrMatrix code downstream (dense_block's
// lower_bound window search in detect_loops_impl.cpp, for one) assumes the
// per-row order the format promises.
void sort_row_columns(CsrMatrix::Arrays& arrays);

// sort_row_columns, then every run of equal columns in a row merged into one
// entry holding their sum, in place. A cooler need not store each (bin1, bin2)
// pixel once (the IMR90 4DN mcool's 250 kb group stores 34,561 pixels twice),
// and hicmatrix builds its matrix as scipy's csr_matrix((data, (bin1, bin2))),
// whose COO to CSR conversion sums duplicates (cpp/PLAN.md 2.3 rule 1). Every
// read of a cooler's pixels into a CsrMatrix goes through this, so every tool
// sees the summed matrix. Rows already strictly increasing are left as they
// are, so a well formed file costs one scan.
void sort_and_sum_row_columns(CsrMatrix::Arrays& arrays);

// One chunk of stored pixels, in storage order.
struct PixelChunk {
    std::vector<std::int64_t> bin1;
    std::vector<std::int64_t> bin2;
    std::vector<double> count;
};

// A cooler opened through coolercpp, with its metadata in the form the tools
// consume.
class CoolFile {
  public:
    // A plain path or a cooler URI "file.mcool::/resolutions/10000".
    explicit CoolFile(const std::string& uri);

    [[nodiscard]] const std::string& filename() const noexcept { return filename_; }
    [[nodiscard]] const std::string& root() const noexcept { return root_; }

    // Cooler.info: string attributes JSON decoded, the rest as numbers.
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

    // Cooler.bins().columns: chrom, start, end, then the other bin columns.
    [[nodiscard]] std::vector<std::string> bin_columns() const;
    [[nodiscard]] std::int64_t nbins() const;
    [[nodiscard]] std::int64_t nnz() const;
    // Cooler.extent(region): the bin range [first, last) of a region such as
    // "chrX:3000000-3500000" (cooler.util.parse_region against the chromosome
    // sizes, then region_to_extent).
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> extent(const std::string& region) const;

    // The bin table as cut intervals with hicmatrix's constant extra 1.0.
    [[nodiscard]] std::vector<CutInterval> read_bins() const;
    [[nodiscard]] bool has_column(const std::string& column) const;
    [[nodiscard]] std::vector<double> read_column(const std::string& column) const;
    // The numpy dtype name of the pixel count column.
    [[nodiscard]] std::string count_dtype() const;

    // The stored pixel table as a CSR matrix, without balancing.
    [[nodiscard]] CsrMatrix read_matrix() const;
    // The stored pixels with first <= bin1, bin2 < last, renumbered from
    // `first`, columns sorted within each row: the chromosome block of a
    // contiguous bin range, read without touching the other rows.
    [[nodiscard]] CsrMatrix read_block(std::int64_t first, std::int64_t last) const;
    // The stored pixels of rows [row_first, row_last) whose column lies in
    // [col_first, col_last), in storage order, one chunk at a time.
    void for_each_pixel_chunk(std::int64_t row_first, std::int64_t row_last,
                              std::int64_t col_first, std::int64_t col_last,
                              const std::function<void(const PixelChunk&)>& visit) const;

  private:
    // read_matrix (first = 0, last = nbins) or read_block through
    // indexes/bin1_offset, or nullopt when the index is missing or does not
    // describe the pixel table.
    [[nodiscard]] std::optional<CsrMatrix> read_from_index(std::int64_t first, std::int64_t last,
                                                           const std::string& dtype) const;

    std::shared_ptr<const coolercpp::Cooler> cooler_;
    std::string filename_;
    std::string root_ = "/";
    std::map<std::string, json::Value> info_;
    std::vector<std::string> chrom_names_;
    std::vector<std::int64_t> chrom_lengths_;
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
    // Cool.distance (hiCMatrix pDistance), used together with chrom_name on a
    // fixed bin size cooler (cool.py:122-145): only pixels with
    // bin2 - bin1 < distance / binsize are kept, their counts as float32.
    // The weights are applied afterwards like on any other load.
    std::optional<std::int64_t> distance;
};

struct CoolLoadResult {
    MatrixData data;
    // The state the loader leaves behind on the Cool object and that a
    // following save reads back (hicmatrix/lib/cool.py:195-207).
    std::optional<char> correction_operator;
    std::optional<std::string> hic2cool_version;
    std::optional<std::string> hicmatrix_version;
    // Cooler.info with every value rendered the way Python's str() would.
    std::map<std::string, std::string> metadata;
};

// The correction a load applies: Cool.load's weight step (cool.py:165-200),
// the factors of every bin of the file and the operator. nullopt when the load
// applies none (apply_correction off, no such column, every factor NaN, or at
// most one stored pixel).
struct CoolCorrection {
    std::vector<double> factors;
    char op = '*';
};
[[nodiscard]] std::optional<CoolCorrection> cool_correction(const CoolFile& cool,
                                                            const CoolLoadOptions& options);

// CoolFile::for_each_pixel_chunk with the load's correction applied to every
// count; a pixel whose corrected value is NaN or 0 is dropped, as Cool.load
// drops it. For tools that stream a band instead of holding the matrix.
void for_each_corrected_pixel_chunk(const CoolFile& cool, std::int64_t row_first,
                                    std::int64_t row_last, std::int64_t col_first,
                                    std::int64_t col_last, const CoolLoadOptions& options,
                                    const std::function<void(const PixelChunk&)>& visit);

// hicmatrix.lib.Cool.load.
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
    // 'w', and hicmatrix's provenance is not written onto the file root.
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
    // The provenance strings, passed to coolercpp and written onto the root.
    // cpp/PLAN.md 5.0.1 decides when they change; until then the port emits
    // the hicmatrix identity verbatim.
    std::string generated_by = "HiCMatrix-17.2";
    std::string generated_by_cooler_lib = "cooler-0.10.2";
    std::string tool_url = "https://github.com/deeptools/HiCMatrix";
    std::string format_url = "https://github.com/mirnylab/cooler";
    // ISO 8601 local time, like datetime.now().isoformat(). Empty means "take
    // the current time"; a fixed value makes a test reproducible.
    std::string creation_date;
};

// hicmatrix.lib.Cool.create_cooler_input and Cool.save, through
// coolercpp::create_cooler. `path` is a plain file name or "file::/group".
// `data` is modified in place exactly where the Python modifies it.
void write_cool(const std::string& path, MatrixData& data,
                const CoolSaveOptions& options = CoolSaveOptions());

}  // namespace hicx

#endif  // HICX_COOL_ADAPTER_HPP
