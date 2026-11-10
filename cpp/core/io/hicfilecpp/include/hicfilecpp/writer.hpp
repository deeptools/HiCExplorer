#ifndef HICFILECPP_WRITER_HPP
#define HICFILECPP_WRITER_HPP

// Writing .hic files, versions 8 and 9. writeHicFile corresponds to Juicer
// tools "pre" (followed by its normalization step) and addNorm to "addNorm":
// version 8 files follow the layout and arithmetic of Juicer tools 1.22.01,
// version 9 files those of 2.20.00. Expected values and normalization vectors
// agree with Juicer's to within floating point rounding; blocks may be
// encoded differently, so the bytes are not Juicer's. The output depends only
// on the input and the options, never on the number of threads.

#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "hicfilecpp/errors.hpp"

namespace hicfilecpp {

// A pixel handed to the writer: bins of the requested resolution, local to
// their chromosomes, and the count.
struct Pixel {
    int32_t bin1 = 0;
    int32_t bin2 = 0;
    float count = 0.0f;
};

class PixelSource {
public:
    virtual ~PixelSource() = default;

    // Hands every pixel of the chromosome pair (chr1, chr2) at `resolution`
    // to `consume`, in any number of batches. chr1 <= chr2 are 0-based
    // indexes into WriteOptions::chromosomes; bin1 lies on chr1 and bin2 on
    // chr2. An intra-chromosomal pixel is given once, in either triangle;
    // repeated pixels are summed. The writer asks for one pair and resolution
    // at a time, in pair order, and may ask for a pair once per resolution.
    virtual void pixels(int32_t resolution, int32_t chr1, int32_t chr2,
                        const std::function<void(const Pixel* pixels, size_t count)>& consume) = 0;
};

struct WriteOptions {
    // 8 or 9. Versions 6 and 7 are refused with their own message: no Juicer
    // tools release that writes them can be obtained to validate against
    // (docs/PROVENANCE.md).
    int32_t version = 9;
    std::string genomeId;
    // Name and length in base pairs, in file order. The whole-genome
    // pseudo-chromosome "All" is added in front, as Juicer tools does.
    std::vector<std::pair<std::string, int64_t>> chromosomes;
    // Base pair resolutions, written from the coarsest to the finest.
    std::vector<int32_t> resolutions;
    // The resolution of the pixels the source hands over. Each pixel counts
    // as one contact at the start of its bins, so every resolution, which
    // must be a multiple of this one, holds what Juicer tools pre makes of
    // the same contacts.
    int32_t sourceResolution = 0;
    // Ask the source for pixels at every resolution instead of binning the
    // source resolution. The whole-genome matrix then uses the finest one.
    bool sourceProvidesEveryResolution = false;
    // Computed as Juicer tools addNorm does; any of VC, VC_SQRT, KR, SCALE.
    std::vector<std::string> normalizations{"VC", "VC_SQRT", "KR", "SCALE"};
    // Normalization vectors the caller supplies instead of having them
    // computed, such as the vectors of another .hic file kept in a cool file.
    // A name is any normalization type label (KR, GW_KR, a custom one) that
    // `normalizations` does not also list. For every chromosome with
    // intra-chromosomal pixels at a resolution, `normVector` is asked for the
    // vector of each name, chrIndex being the 0-based index into
    // `chromosomes`; an empty vector writes none for that chromosome. The
    // values are stored as given, without Juicer's scaling to the matrix sum,
    // and the normalized expected values are computed from them. A vector
    // shorter than the matrix grid is padded as Juicer pads its own: with 0
    // for VC types, NaN otherwise; a longer one is cut.
    std::vector<std::string> providedNormalizations;
    std::function<std::vector<double>(const std::string& name, int32_t chrIndex, int32_t resolution)> normVector;
    // Header attributes written after "software", in order.
    std::vector<std::pair<std::string, std::string>> attributes;
    // The "software" attribute; empty means "hicfilecpp <version>".
    std::string software;
    // Threads for block encoding and compression.
    int threads = 1;
    // zlib compression level; -1 is zlib's default, which Juicer tools uses.
    int compressionLevel = -1;
};

// Writes a .hic file from the source's pixels.
void writeHicFile(const std::string& path, const WriteOptions& options, PixelSource& source);

// Computes the normalization vectors and normalized expected values of an
// existing version 8 or 9 file with base pair resolutions only, replacing any
// it holds, as Juicer tools addNorm does.
void addNorm(const std::string& path, const std::vector<std::string>& normalizations, int threads = 1);

}  // namespace hicfilecpp

#endif  // HICFILECPP_WRITER_HPP
