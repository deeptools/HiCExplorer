// Reading and writing Juicer .hic files for HiCExplorer v4.
//
// The file format itself is hicfilecpp's (the Juicer .hic reader and writer
// the build obtains through cmake/HicxHicfilecpp.cmake). This adapter keeps
// what HiCExplorer adds on top of it:
//
//   * hic2cool_convert: the Python reference converts .hic to cool with
//     hic2cool (hicConvertFormat.py:124-138). The port writes the same cool
//     and mcool files, dataset for dataset: hic2cool's layout, dtypes, chunks
//     and filters, its attributes, its norm columns and its pixel order.
//   * write_hic: HiCExplorer matrices (fixed bins, as loaded from h5, cool or
//     mcool) handed to hicfilecpp as pixels at one or several resolutions.

#ifndef HICX_HIC_ADAPTER_HPP
#define HICX_HIC_ADAPTER_HPP

#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "hicx/matrix_data.hpp"

namespace hicx {

// A condition under which hic2cool prints a message and calls sys.exit(1)
// (hic2cool_utils.py force_exit). The message is what it prints.
class Hic2coolExit : public std::runtime_error {
  public:
    explicit Hic2coolExit(const std::string& message) : std::runtime_error(message) {}
};

// hic2cool.hic2cool_convert(infile, outfile, resolution): resolution 0 means
// every base pair resolution of the file. hic2cool renames outfile to .cool
// or .mcool by the number of resolutions written; the path written is
// returned. An existing file at that path is replaced.
std::string hic2cool_convert(const std::string& hic_path, const std::string& outfile,
                             std::int64_t resolution);

// hic2cool's multi resolution layout for a chosen list of resolutions,
// written to exactly `outfile` (no renaming), even for a single resolution.
void hic2cool_convert_mcool(const std::string& hic_path, const std::string& outfile,
                            const std::vector<std::int64_t>& resolutions);

// The base pair resolutions of a .hic file, in file order.
[[nodiscard]] std::vector<std::int64_t> hic_resolutions(const std::string& hic_path);

struct HicWriteOptions {
    // .hic version, 8 or 9.
    int version = 8;
    std::vector<std::string> normalizations{"VC", "VC_SQRT", "KR", "SCALE"};
    std::string genome = "unknown";
    int threads = 1;
};

// Writes HiCExplorer matrices as a .hic file. With one matrix, its bin size
// and every value in extra_resolutions (each a multiple of that bin size) are
// written, the coarser ones binned from it. With several matrices, as loaded
// from the resolutions of an mcool file, each is written as its own
// resolution and extra_resolutions must be empty. All matrices need fixed
// size bins starting at 0 and the same chromosomes. The upper triangle's
// finite, non-zero values become the pixels.
void write_hic(const std::string& path, const std::vector<const MatrixData*>& matrices,
               const std::vector<std::int64_t>& extra_resolutions,
               const HicWriteOptions& options);

// ----------------------------------------------------------- generic loading
//
// hicmatrix.HiCMatrix.__init__ (hicmatrix.py:47-51) has never known about
// .hic; every matrix-reading tool in v4 gains it here, through ToolMatrix,
// with a resolution and normalisation selector modelled on the mcool
// convention the tools already use ("file.mcool::/resolutions/10000",
// hicConvertFormat --help). The .hic selector nests the same way:
//
//   file.hic                                        finest resolution, raw counts
//   file.hic::/resolutions/10000                     10000 bp, raw counts
//   file.hic::/resolutions/10000/normalizations/KR    10000 bp, KR applied
//   file.hic::/normalizations/KR                     finest resolution, KR applied
//
// "raw counts" and "finest resolution" are hic2cool_convert's own defaults
// (resolution 0 means every resolution, and hicConvertFormat's own default
// --correction_name "weight" names a column no hic2cool cool file has, so by
// default nothing gets applied); the selector follows them rather than
// choosing new defaults independently.

// The literal bytes "HIC" hicfilecpp's reader itself requires at the start of
// a file (src/reader.cpp: "Hi-C magic string is missing"). Cheap enough to
// call before deciding a file's format, and preferred over the ".hic"
// extension where the two disagree. False when the file cannot be opened or
// is shorter than 3 bytes.
[[nodiscard]] bool has_hic_signature(const std::string& path);

// True when `path`, with any "::" selector suffix stripped, is a .hic file:
// its content signature first, the ".hic" extension as a fallback.
[[nodiscard]] bool is_hic_path(const std::string& path);

// One parsed selector.
struct HicUri {
    std::string path;
    std::optional<std::int64_t> resolution;
    std::optional<std::string> normalization;
};

// Parses the selector documented above. Throws std::runtime_error on a
// suffix that is not one of the four forms.
[[nodiscard]] HicUri parse_hic_uri(const std::string& uri);

// What ToolMatrix::load needs from a .hic source, matching hicx::CoolLoadResult
// closely enough that the caller can treat it the same way.
struct HicLoadResult {
    MatrixData data;
    std::optional<char> correction_operator;
    std::map<std::string, std::string> metadata;
};

// hicmatrix.lib.Cool.load for a .hic source, in the same two shapes
// CoolLoadOptions.chrom_name gives the cool loader (cool_adapter.hpp):
//
//  * `chrom_name` empty: the whole file, at the selected resolution and
//    normalisation. Goes through hic2cool_convert (the same, already
//    validated conversion hicConvertFormat runs) into a temporary single
//    resolution cool file, removed before this returns, and then
//    hicx::read_cool on it; not a reimplementation of .hic parsing.
//  * `chrom_name` a bare chromosome name or a "chrom:start-end" region
//    (coolercpp::parse_region_string, the same parser cool's own
//    CoolFile::extent uses, so both formats accept the same region syntax):
//    only that block is read, through hicfilecpp's MatrixZoomData::getRecords,
//    which itself decodes only the blocks the requested genomic range's block
//    index selects (hicfilecpp/src/reader.cpp), not the whole chromosome or
//    file. This is the .hic side of cpp/PLAN.md's tier 10 note that a matrix
//    is never loaded whole for an on screen region, and it is what lets
//    hicPlotMatrix's existing "a cool input without --region2 ... is loaded
//    as the region or chromosome only" fast path (hicPlotMatrix.cpp's header
//    comment) carry over to .hic without loading the file hic2cool_convert
//    would have to touch in full.
//
// Either way the requested normalisation, when not "NONE", is applied as
// hic2cool's own tables always are: divisively (observed / normI / normJ),
// matching what hicConvertFormat's --correction_division names explicitly
// rather than leaving it to the column-name heuristic cool files use for
// their own 'weight' column.
[[nodiscard]] HicLoadResult read_hic(const std::string& uri,
                                     const std::optional<std::string>& chrom_name = std::nullopt);

}  // namespace hicx

#endif  // HICX_HIC_ADAPTER_HPP
