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

}  // namespace hicx

#endif  // HICX_HIC_ADAPTER_HPP
