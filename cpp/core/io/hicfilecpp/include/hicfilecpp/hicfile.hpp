#ifndef HICFILECPP_HICFILE_HPP
#define HICFILECPP_HICFILE_HPP

// Reading .hic files, versions 6 to 9. HiCFile and MatrixZoomData reproduce
// hicstraw 1.3.1 (the pybind11 module of aidenlab/straw): the same names, the
// same arguments and the same values, including hicstraw's float arithmetic.
// The members below the first block of each class go beyond hicstraw: they
// expose the file's metadata, its vectors and its blocks for bulk reading.

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "hicfilecpp/errors.hpp"

namespace hicfilecpp {

class HiCFile;

namespace detail {
struct FileState;
struct ZoomState;
const FileState& stateOf(const HiCFile& file);
}  // namespace detail

// hicstraw.chromosome
struct Chromosome {
    std::string name;
    int32_t index = 0;
    int64_t length = 0;
};

// hicstraw.contactRecord. getRecords reports binX and binY as genomic start
// positions (bin times resolution), as hicstraw does; readBlock reports bins.
struct ContactRecord {
    int32_t binX = 0;
    int32_t binY = 0;
    float counts = 0.0f;
};

// The float32 array getRecordsAsMatrix returns, row major.
struct FloatMatrix {
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<float> values;

    float at(int64_t row, int64_t col) const { return values[static_cast<size_t>(row * cols + col)]; }
};

// The header of one resolution of one chromosome pair matrix.
struct ZoomHeader {
    std::string unit;
    int32_t zoomIndex = 0;
    float sumCounts = 0.0f;
    float occupiedCellCount = 0.0f;
    float percent5 = 0.0f;
    float percent95 = 0.0f;
    int32_t binSize = 0;
    int32_t blockBinCount = 0;
    int32_t blockColumnCount = 0;
};

// One entry of a zoom level's block index.
struct BlockIndexEntry {
    int32_t number = 0;
    int64_t position = 0;
    int32_t size = 0;
};

// The key of one entry of the footer's expected-value sections. The
// normalization of the raw section is "NONE".
struct ExpectedValuesKey {
    std::string normalization;
    std::string unit;
    int32_t binSize = 0;
};

struct ExpectedValues {
    ExpectedValuesKey key;
    std::vector<double> values;
    // chromosome index, normalization factor, in file order
    std::vector<std::pair<int32_t, double>> normalizationFactors;
};

// One entry of the footer's normalization vector index.
struct NormVectorIndexEntry {
    std::string normalization;
    int32_t chrIndex = 0;
    std::string unit;
    int32_t resolution = 0;
    int64_t position = 0;
    int64_t sizeInBytes = 0;
};

// hicstraw.MatrixZoomData: one chromosome pair at one resolution, with a
// matrix type ("observed", "oe" or "expected") and a normalization.
class MatrixZoomData {
public:
    std::vector<ContactRecord> getRecords(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) const;
    FloatMatrix getRecordsAsMatrix(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) const;
    std::vector<double> getNormVector(int32_t index) const;
    std::vector<double> getExpectedValues() const;
    int64_t getNumberOfTotalRecords() const;

    // False where hicstraw prints a message and returns no records: the pair is
    // absent from the file, or the expected values an "oe" query needs are.
    bool found() const;
    const std::string& message() const;
    int32_t chr1Index() const;
    int32_t chr2Index() const;
    bool isIntra() const;
    int32_t resolution() const;
    const ZoomHeader& zoomHeader() const;
    // In file order.
    const std::vector<BlockIndexEntry>& blockIndex() const;
    // The block's records as stored: bins and raw counts, in block order.
    std::vector<ContactRecord> readBlock(const BlockIndexEntry& entry) const;
    // Decodes every block, in ascending block number, on up to `threads`
    // threads, and hands each to `visit` in that order on the calling thread.
    void forEachBlock(const std::function<void(const BlockIndexEntry&, std::vector<ContactRecord>&)>& visit,
                      int threads = 1) const;

private:
    friend class HiCFile;
    explicit MatrixZoomData(std::shared_ptr<const detail::ZoomState> state);
    std::shared_ptr<const detail::ZoomState> state_;
};

// hicstraw.HiCFile
class HiCFile {
public:
    explicit HiCFile(const std::string& fileName);

    std::string getGenomeID() const;
    std::vector<int32_t> getResolutions() const;
    std::vector<Chromosome> getChromosomes() const;
    MatrixZoomData getMatrixZoomData(const std::string& chr1, const std::string& chr2,
                                     const std::string& matrixType, const std::string& norm,
                                     const std::string& unit, int32_t resolution) const;

    int32_t version() const;
    const std::string& fileName() const;
    int64_t masterIndexPosition() const;
    // The normalization vector index position and length of a version 9
    // header; (0, 0) in version 8.
    std::pair<int64_t, int64_t> normVectorIndexHeader() const;
    const std::vector<std::pair<std::string, std::string>>& attributes() const;
    std::vector<int32_t> getFragResolutions() const;
    // Restriction sites per chromosome, present when FRAG resolutions are.
    std::vector<int32_t> fragmentSiteCounts() const;
    bool hasMatrix(int32_t chr1Index, int32_t chr2Index) const;
    std::vector<ZoomHeader> matrixZoomHeaders(int32_t chr1Index, int32_t chr2Index) const;
    // False when the footer ends after the raw expected values, with no
    // section for normalized ones (hic2cool warns about such files).
    bool hasNormalizedExpectedSection() const;
    // Every expected-value entry, raw section first, in file order.
    std::vector<ExpectedValuesKey> expectedValuesKeys() const;
    std::optional<ExpectedValues> readExpectedValues(const ExpectedValuesKey& key) const;
    const std::vector<NormVectorIndexEntry>& normVectorIndex() const;
    // Normalization types with at least one vector, in index order.
    std::vector<std::string> getNormalizationTypes() const;
    std::optional<std::vector<double>> readNormVector(const std::string& norm, int32_t chrIndex,
                                                      const std::string& unit, int32_t resolution) const;

private:
    friend const detail::FileState& detail::stateOf(const HiCFile& file);
    std::shared_ptr<const detail::FileState> state_;
};

}  // namespace hicfilecpp

#endif  // HICFILECPP_HICFILE_HPP
