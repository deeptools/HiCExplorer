#ifndef HICFILECPP_DETAIL_STATE_HPP
#define HICFILECPP_DETAIL_STATE_HPP

// The parsed index of an open .hic file and of one MatrixZoomData query,
// shared by the reader and by addNorm.

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "hicfilecpp/hicfile.hpp"
#include "io.hpp"

namespace hicfilecpp::detail {

struct VectorLocation {
    int64_t valuesOffset = 0;
    int64_t nValues = 0;
    int64_t factorsOffset = 0;
    int32_t nFactors = 0;
};

struct ExpectedEntry {
    ExpectedValuesKey key;
    VectorLocation location;
};

struct MatrixEntry {
    int64_t position = 0;
    int32_t size = 0;
};

struct FileState {
    explicit FileState(const std::string& fileName);

    std::vector<double> readValues(int64_t offset, int64_t count) const;
    std::vector<std::pair<int32_t, double>> readFactors(int64_t offset, int32_t count) const;
    std::vector<double> readNormVectorAt(int64_t position) const;
    void decodeBlock(const BlockIndexEntry& entry, std::vector<char>& scratch,
                     std::vector<ContactRecord>& out) const;

    std::string path;
    File file;
    int32_t version = 0;
    int64_t master = 0;
    std::string genome;
    int64_t nviPosition = 0;
    int64_t nviLength = 0;
    std::vector<std::pair<std::string, std::string>> attributes;
    std::vector<Chromosome> chromosomes;
    std::map<std::string, Chromosome> chromosomeByName;
    std::vector<int32_t> bpResolutions;
    std::vector<int32_t> fragResolutions;
    std::vector<int32_t> fragmentSiteCounts;
    std::map<std::string, MatrixEntry> matrices;
    std::vector<ExpectedEntry> expectedNone;
    std::vector<ExpectedEntry> expectedNormalized;
    std::vector<NormVectorIndexEntry> normIndex;
    // Where the normalized expected values begin: the end of the raw
    // expected-value section, which is where Juicer tools addNorm writes.
    int64_t normalizedSectionPosition = 0;
    // False when the file ends right after the raw expected values, as files
    // written by old Juicer releases without normalizations do.
    bool normalizedSectionPresent = false;
};

struct ZoomState {
    std::set<int32_t> blockNumbers(const int64_t* regionIndices) const;

    std::shared_ptr<const FileState> file;
    bool isIntra = false;
    int32_t c1 = 0;
    int32_t c2 = 0;
    std::string matrixType;
    std::string norm;
    std::string unit;
    int32_t resolution = 0;
    int32_t numBins1 = 0;
    int32_t numBins2 = 0;
    bool foundFooter = false;
    std::string message;
    std::vector<double> expectedValues;
    std::vector<double> c1Norm;
    std::vector<double> c2Norm;
    ZoomHeader header;
    std::map<int32_t, BlockIndexEntry> blockMap;
    std::vector<BlockIndexEntry> blockList;
    double avgCount = 0.0;
};

}  // namespace hicfilecpp::detail

#endif  // HICFILECPP_DETAIL_STATE_HPP
