// Reading .hic files. The query semantics follow hicstraw 1.3.1
// (pybind11_python/src/straw.cpp of aidenlab/straw at 82fba9cf); the notes
// name the hicstraw function each part reproduces.

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <thread>

#include "codec.hpp"
#include "hicfilecpp/hicfile.hpp"
#include "io.hpp"
#include "state.hpp"

namespace hicfilecpp {

namespace detail {

namespace {

const std::string kNone = "NONE";

std::string pairKey(int32_t c1, int32_t c2) {
    return std::to_string(c1) + "_" + std::to_string(c2);
}

int64_t checkedCount(int64_t count, const std::string& what, const std::string& path) {
    if (count < 0) {
        throw HicError("negative " + what + " count in " + path);
    }
    return count;
}

}  // namespace

// readHeader, readResolutionsFromHeader and readFooter of hicstraw, reading the
// whole footer index once instead of once per query.
FileState::FileState(const std::string& fileName) : path(fileName), file(fileName) {
    SequentialReader in(file, 0);
    std::string magic;
    try {
        magic = in.cstr(64);
    } catch (const HicError&) {
        magic.clear();
    }
    if (magic.size() < 3 || magic.compare(0, 3, "HIC") != 0) {
        throw HicError("Hi-C magic string is missing, does not appear to be a hic file");
    }
    version = in.get<int32_t>();
    // hicstraw refuses versions below 6 (readHeader); versions 6 and 7 share
    // version 8's header, footer, zoom, expected value and vector layout, and
    // version 6 differs only in its block records (decodeBlock).
    if (version < 6 || version > 9) {
        throw HicError("Version " + std::to_string(version) +
                       " is not supported: hicfilecpp reads .hic versions 6 to 9");
    }
    master = in.get<int64_t>();
    genome = in.cstr();
    if (version > 8) {
        nviPosition = in.get<int64_t>();
        nviLength = in.get<int64_t>();
    }
    const int64_t nAttributes = checkedCount(in.get<int32_t>(), "attribute", path);
    for (int64_t i = 0; i < nAttributes; ++i) {
        std::string key = in.cstr();
        std::string value = in.cstr();
        attributes.emplace_back(std::move(key), std::move(value));
    }
    const int64_t nChromosomes = checkedCount(in.get<int32_t>(), "chromosome", path);
    for (int64_t i = 0; i < nChromosomes; ++i) {
        Chromosome chromosome;
        chromosome.name = in.cstr();
        chromosome.length = version > 8 ? in.get<int64_t>() : static_cast<int64_t>(in.get<int32_t>());
        chromosome.index = static_cast<int32_t>(i);
        chromosomeByName[chromosome.name] = chromosome;
        chromosomes.push_back(std::move(chromosome));
    }
    const int64_t nBp = checkedCount(in.get<int32_t>(), "resolution", path);
    for (int64_t i = 0; i < nBp; ++i) {
        bpResolutions.push_back(in.get<int32_t>());
    }
    const int64_t nFrag = checkedCount(in.get<int32_t>(), "resolution", path);
    for (int64_t i = 0; i < nFrag; ++i) {
        fragResolutions.push_back(in.get<int32_t>());
    }
    if (nFrag > 0) {
        for (int64_t i = 0; i < nChromosomes; ++i) {
            const int32_t nSites = in.get<int32_t>();
            fragmentSiteCounts.push_back(nSites);
            in.skip(4 * static_cast<int64_t>(std::max(nSites, 0)));
        }
    }

    if (master <= 0 || master >= file.size()) {
        throw HicError("master index position " + std::to_string(master) + " lies outside " + path);
    }
    SequentialReader footer(file, master);
    if (version > 8) {
        footer.get<int64_t>();
    } else {
        footer.get<int32_t>();
    }
    const int64_t nEntries = checkedCount(footer.get<int32_t>(), "matrix", path);
    for (int64_t i = 0; i < nEntries; ++i) {
        std::string key = footer.cstr();
        MatrixEntry entry;
        entry.position = footer.get<int64_t>();
        entry.size = footer.get<int32_t>();
        matrices[key] = entry;
    }
    const int64_t valueSize = version > 8 ? 4 : 8;
    auto readExpected = [&](bool normalized, std::vector<ExpectedEntry>& out) {
        const int64_t count = checkedCount(footer.get<int32_t>(), "expected value", path);
        for (int64_t i = 0; i < count; ++i) {
            ExpectedEntry entry;
            entry.key.normalization = normalized ? footer.cstr() : kNone;
            entry.key.unit = footer.cstr();
            entry.key.binSize = footer.get<int32_t>();
            entry.location.nValues = checkedCount(
                version > 8 ? footer.get<int64_t>() : static_cast<int64_t>(footer.get<int32_t>()), "value", path);
            entry.location.valuesOffset = footer.tell();
            footer.skip(entry.location.nValues * valueSize);
            entry.location.nFactors = static_cast<int32_t>(checkedCount(footer.get<int32_t>(), "factor", path));
            entry.location.factorsOffset = footer.tell();
            footer.skip(static_cast<int64_t>(entry.location.nFactors) * (4 + valueSize));
            if (footer.tell() > file.size()) {
                throw HicError("truncated expected values in " + path);
            }
            out.push_back(std::move(entry));
        }
    };
    normalizedSectionPosition = footer.tell();
    if (footer.eof()) {
        return;
    }
    readExpected(false, expectedNone);
    normalizedSectionPosition = footer.tell();
    if (footer.eof()) {
        return;
    }
    normalizedSectionPresent = true;
    readExpected(true, expectedNormalized);
    if (footer.eof()) {
        return;
    }
    const int64_t nNorm = checkedCount(footer.get<int32_t>(), "normalization vector", path);
    for (int64_t i = 0; i < nNorm; ++i) {
        NormVectorIndexEntry entry;
        entry.normalization = footer.cstr();
        entry.chrIndex = footer.get<int32_t>();
        entry.unit = footer.cstr();
        entry.resolution = footer.get<int32_t>();
        entry.position = footer.get<int64_t>();
        entry.sizeInBytes = version > 8 ? footer.get<int64_t>() : static_cast<int64_t>(footer.get<int32_t>());
        normIndex.push_back(std::move(entry));
    }
}

std::vector<double> FileState::readValues(int64_t offset, int64_t count) const {
    std::vector<double> values(static_cast<size_t>(count));
    if (count == 0) {
        return values;
    }
    if (version > 8) {
        std::vector<float> raw(static_cast<size_t>(count));
        file.readAt(offset, reinterpret_cast<char*>(raw.data()), raw.size() * sizeof(float));
        std::copy(raw.begin(), raw.end(), values.begin());
    } else {
        file.readAt(offset, reinterpret_cast<char*>(values.data()), values.size() * sizeof(double));
    }
    return values;
}

std::vector<std::pair<int32_t, double>> FileState::readFactors(int64_t offset, int32_t count) const {
    const size_t valueSize = version > 8 ? 4 : 8;
    std::vector<char> raw(static_cast<size_t>(count) * (4 + valueSize));
    if (!raw.empty()) {
        file.readAt(offset, raw.data(), raw.size());
    }
    std::vector<std::pair<int32_t, double>> factors;
    factors.reserve(static_cast<size_t>(count));
    MemReader in(raw.data(), raw.size());
    for (int32_t i = 0; i < count; ++i) {
        const int32_t chrIndex = in.get<int32_t>();
        const double factor = version > 8 ? static_cast<double>(in.get<float>()) : in.get<double>();
        factors.emplace_back(chrIndex, factor);
    }
    return factors;
}

// readNormalizationVector
std::vector<double> FileState::readNormVectorAt(int64_t position) const {
    int64_t count = 0;
    if (version > 8) {
        file.readAt(position, reinterpret_cast<char*>(&count), sizeof(int64_t));
        position += 8;
    } else {
        int32_t count32 = 0;
        file.readAt(position, reinterpret_cast<char*>(&count32), sizeof(int32_t));
        count = count32;
        position += 4;
    }
    return readValues(position, checkedCount(count, "normalization vector value", path));
}

namespace {

inline void appendRecord(std::vector<ContactRecord>& out, int32_t binX, int32_t binY, float counts) {
    out.push_back(ContactRecord{binX, binY, counts});
}

template <class RowIndex, class ColIndex, class Value>
void decodeRows(MemReader& in, int32_t binXOffset, int32_t binYOffset, std::vector<ContactRecord>& out) {
    const auto rowCount = static_cast<int64_t>(in.get<RowIndex>());
    for (int64_t i = 0; i < rowCount; ++i) {
        in.require(sizeof(RowIndex) + sizeof(ColIndex));
        const int32_t binY = binYOffset + static_cast<int32_t>(in.getUnchecked<RowIndex>());
        const auto colCount = static_cast<int64_t>(in.getUnchecked<ColIndex>());
        if (colCount <= 0) {
            continue;
        }
        in.require(static_cast<size_t>(colCount) * (sizeof(ColIndex) + sizeof(Value)));
        for (int64_t j = 0; j < colCount; ++j) {
            const int32_t binX = binXOffset + static_cast<int32_t>(in.getUnchecked<ColIndex>());
            const auto counts = static_cast<float>(in.getUnchecked<Value>());
            appendRecord(out, binX, binY, counts);
        }
    }
}

}  // namespace

// readBlock
void FileState::decodeBlock(const BlockIndexEntry& entry, std::vector<char>& scratch,
                            std::vector<ContactRecord>& out) const {
    if (entry.size <= 0) {
        return;
    }
    std::vector<char> compressed(static_cast<size_t>(entry.size));
    file.readAt(entry.position, compressed.data(), compressed.size());
    inflateZlib(compressed.data(), compressed.size(), scratch);
    MemReader in(scratch.data(), scratch.size());
    const int32_t nRecords = in.get<int32_t>();
    out.reserve(out.size() + static_cast<size_t>(std::max(nRecords, 0)));
    if (version < 7) {
        // Version 6: nRecords plain records of int32 binX, int32 binY and
        // float32 counts.
        in.require(static_cast<size_t>(std::max(nRecords, 0)) * 12);
        for (int32_t i = 0; i < nRecords; ++i) {
            const auto binX = in.getUnchecked<int32_t>();
            const auto binY = in.getUnchecked<int32_t>();
            const auto counts = in.getUnchecked<float>();
            appendRecord(out, binX, binY, counts);
        }
        return;
    }
    const int32_t binXOffset = in.get<int32_t>();
    const int32_t binYOffset = in.get<int32_t>();
    const bool useShort = in.get<char>() == 0;  // yes, 0 means short counts
    bool useShortBinX = true;
    bool useShortBinY = true;
    if (version > 8) {
        useShortBinX = in.get<char>() == 0;
        useShortBinY = in.get<char>() == 0;
    }
    const char type = in.get<char>();
    if (type == 1) {
        if (useShortBinX && useShortBinY) {
            if (useShort) {
                decodeRows<int16_t, int16_t, int16_t>(in, binXOffset, binYOffset, out);
            } else {
                decodeRows<int16_t, int16_t, float>(in, binXOffset, binYOffset, out);
            }
        } else if (useShortBinX && !useShortBinY) {
            if (useShort) {
                decodeRows<int32_t, int16_t, int16_t>(in, binXOffset, binYOffset, out);
            } else {
                decodeRows<int32_t, int16_t, float>(in, binXOffset, binYOffset, out);
            }
        } else if (!useShortBinX && useShortBinY) {
            if (useShort) {
                decodeRows<int16_t, int32_t, int16_t>(in, binXOffset, binYOffset, out);
            } else {
                decodeRows<int16_t, int32_t, float>(in, binXOffset, binYOffset, out);
            }
        } else {
            if (useShort) {
                decodeRows<int32_t, int32_t, int16_t>(in, binXOffset, binYOffset, out);
            } else {
                decodeRows<int32_t, int32_t, float>(in, binXOffset, binYOffset, out);
            }
        }
    } else if (type == 2) {
        const int32_t nPts = in.get<int32_t>();
        const int16_t w = in.get<int16_t>();
        if (nPts > 0 && w == 0) {
            throw HicError("dense block with zero width in " + path);
        }
        in.require(static_cast<size_t>(std::max(nPts, 0)) * (useShort ? 2 : 4));
        for (int32_t i = 0; i < nPts; ++i) {
            const int32_t row = i / w;
            const int32_t col = i - row * w;
            const int32_t bin1 = binXOffset + col;
            const int32_t bin2 = binYOffset + row;
            if (useShort) {
                const auto c = in.getUnchecked<int16_t>();
                if (c != -32768) {
                    appendRecord(out, bin1, bin2, static_cast<float>(c));
                }
            } else {
                const auto counts = in.getUnchecked<float>();
                if (!std::isnan(counts)) {
                    appendRecord(out, bin1, bin2, counts);
                }
            }
        }
    }
}

// getBlockNumbersForRegionFromBinPosition and
// getBlockNumbersForRegionFromBinPositionV9Intra
std::set<int32_t> ZoomState::blockNumbers(const int64_t* r) const {
    std::set<int32_t> blocks;
    const int32_t bbc = header.blockBinCount;
    const int32_t bcc = header.blockColumnCount;
    if (bbc <= 0) {
        return blocks;
    }
    auto number = [](int64_t row, int64_t col, int64_t columns) {
        return static_cast<int32_t>(row * columns + col);
    };
    if (file->version > 8 && isIntra) {
        const auto lowerPad = static_cast<int32_t>((r[0] + r[2]) / 2 / bbc);
        const auto higherPad = static_cast<int32_t>((r[1] + r[3]) / 2 / bbc + 1);
        const auto nearer =
            static_cast<int32_t>(std::log2(1 + static_cast<double>(std::abs(r[0] - r[3])) / std::sqrt(2) / bbc));
        const auto further =
            static_cast<int32_t>(std::log2(1 + static_cast<double>(std::abs(r[1] - r[2])) / std::sqrt(2) / bbc));
        int32_t nearerDepth = std::min(nearer, further);
        if ((r[0] > r[3] && r[1] < r[2]) || (r[1] > r[2] && r[0] < r[3])) {
            nearerDepth = 0;
        }
        const int32_t furtherDepth = std::max(nearer, further) + 1;
        for (int64_t depth = nearerDepth; depth <= furtherDepth; ++depth) {
            for (int64_t pad = lowerPad; pad <= higherPad; ++pad) {
                blocks.insert(number(depth, pad, bcc));
            }
        }
        return blocks;
    }
    const auto col1 = static_cast<int32_t>(r[0] / bbc);
    const auto col2 = static_cast<int32_t>((r[1] + 1) / bbc);
    const auto row1 = static_cast<int32_t>(r[2] / bbc);
    const auto row2 = static_cast<int32_t>((r[3] + 1) / bbc);
    for (int64_t row = row1; row <= row2; ++row) {
        for (int64_t col = col1; col <= col2; ++col) {
            blocks.insert(number(row, col, bcc));
        }
    }
    if (isIntra) {
        for (int64_t row = col1; row <= col2; ++row) {
            for (int64_t col = row1; col <= row2; ++col) {
                blocks.insert(number(row, col, bcc));
            }
        }
    }
    return blocks;
}

const FileState& stateOf(const HiCFile& file) {
    return *file.state_;
}

}  // namespace detail

using detail::FileState;
using detail::ZoomState;

namespace {

double normAt(const std::vector<double>& vector, int32_t index) {
    if (index < 0 || static_cast<size_t>(index) >= vector.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return vector[static_cast<size_t>(index)];
}

std::string describe(int32_t resolution, const std::string& unit) {
    return std::to_string(resolution) + " " + unit;
}

}  // namespace

// ---------------------------------------------------------------- HiCFile

HiCFile::HiCFile(const std::string& fileName) : state_(std::make_shared<const FileState>(fileName)) {}

std::string HiCFile::getGenomeID() const {
    return state_->genome;
}

std::vector<int32_t> HiCFile::getResolutions() const {
    return state_->bpResolutions;
}

std::vector<Chromosome> HiCFile::getChromosomes() const {
    // hicstraw keys chromosomes by name, so a repeated name keeps the last one.
    std::vector<Chromosome> result(state_->chromosomes.size());
    for (const auto& entry : state_->chromosomeByName) {
        result[static_cast<size_t>(entry.second.index)] = entry.second;
    }
    return result;
}

int32_t HiCFile::version() const {
    return state_->version;
}

const std::string& HiCFile::fileName() const {
    return state_->path;
}

int64_t HiCFile::masterIndexPosition() const {
    return state_->master;
}

std::pair<int64_t, int64_t> HiCFile::normVectorIndexHeader() const {
    return {state_->nviPosition, state_->nviLength};
}

const std::vector<std::pair<std::string, std::string>>& HiCFile::attributes() const {
    return state_->attributes;
}

std::vector<int32_t> HiCFile::getFragResolutions() const {
    return state_->fragResolutions;
}

std::vector<int32_t> HiCFile::fragmentSiteCounts() const {
    return state_->fragmentSiteCounts;
}

bool HiCFile::hasMatrix(int32_t chr1Index, int32_t chr2Index) const {
    return state_->matrices.count(detail::pairKey(chr1Index, chr2Index)) > 0;
}

std::vector<ZoomHeader> HiCFile::matrixZoomHeaders(int32_t chr1Index, int32_t chr2Index) const {
    std::vector<ZoomHeader> headers;
    const auto it = state_->matrices.find(detail::pairKey(chr1Index, chr2Index));
    if (it == state_->matrices.end()) {
        return headers;
    }
    detail::SequentialReader in(state_->file, it->second.position, 1 << 16);
    in.get<int32_t>();
    in.get<int32_t>();
    const int32_t nRes = in.get<int32_t>();
    for (int32_t i = 0; i < nRes; ++i) {
        ZoomHeader header;
        header.unit = in.cstr();
        header.zoomIndex = in.get<int32_t>();
        header.sumCounts = in.get<float>();
        header.occupiedCellCount = in.get<float>();
        header.percent5 = in.get<float>();
        header.percent95 = in.get<float>();
        header.binSize = in.get<int32_t>();
        header.blockBinCount = in.get<int32_t>();
        header.blockColumnCount = in.get<int32_t>();
        const int32_t nBlocks = in.get<int32_t>();
        in.skip(16 * static_cast<int64_t>(std::max(nBlocks, 0)));
        headers.push_back(std::move(header));
    }
    return headers;
}

bool HiCFile::hasNormalizedExpectedSection() const {
    return state_->normalizedSectionPresent;
}

std::vector<ExpectedValuesKey> HiCFile::expectedValuesKeys() const {
    std::vector<ExpectedValuesKey> keys;
    for (const auto& entry : state_->expectedNone) {
        keys.push_back(entry.key);
    }
    for (const auto& entry : state_->expectedNormalized) {
        keys.push_back(entry.key);
    }
    return keys;
}

std::optional<ExpectedValues> HiCFile::readExpectedValues(const ExpectedValuesKey& key) const {
    const auto& section = key.normalization == detail::kNone ? state_->expectedNone : state_->expectedNormalized;
    for (const auto& entry : section) {
        if (entry.key.normalization == key.normalization && entry.key.unit == key.unit &&
            entry.key.binSize == key.binSize) {
            ExpectedValues result;
            result.key = entry.key;
            result.values = state_->readValues(entry.location.valuesOffset, entry.location.nValues);
            result.normalizationFactors = state_->readFactors(entry.location.factorsOffset, entry.location.nFactors);
            return result;
        }
    }
    return std::nullopt;
}

const std::vector<NormVectorIndexEntry>& HiCFile::normVectorIndex() const {
    return state_->normIndex;
}

std::vector<std::string> HiCFile::getNormalizationTypes() const {
    std::vector<std::string> types;
    for (const auto& entry : state_->normIndex) {
        if (std::find(types.begin(), types.end(), entry.normalization) == types.end()) {
            types.push_back(entry.normalization);
        }
    }
    return types;
}

std::optional<std::vector<double>> HiCFile::readNormVector(const std::string& norm, int32_t chrIndex,
                                                           const std::string& unit, int32_t resolution) const {
    std::optional<int64_t> position;
    for (const auto& entry : state_->normIndex) {
        if (entry.normalization == norm && entry.chrIndex == chrIndex && entry.unit == unit &&
            entry.resolution == resolution) {
            position = entry.position;
        }
    }
    if (!position) {
        return std::nullopt;
    }
    return state_->readNormVectorAt(*position);
}

// MatrixZoomData's constructor in hicstraw: readFooter, the normalization
// vectors and readMatrix.
MatrixZoomData HiCFile::getMatrixZoomData(const std::string& chr1, const std::string& chr2,
                                          const std::string& matrixType, const std::string& norm,
                                          const std::string& unit, int32_t resolution) const {
    const FileState& file = *state_;
    const auto chromosome = [&](const std::string& name) -> const Chromosome& {
        const auto it = file.chromosomeByName.find(name);
        if (it == file.chromosomeByName.end()) {
            throw HicError(name + " not found in the file.");
        }
        return it->second;
    };
    const Chromosome& chrom1 = chromosome(chr1);
    const Chromosome& chrom2 = chromosome(chr2);
    if (resolution <= 0) {
        throw HicError("resolution must be positive, got " + std::to_string(resolution));
    }
    auto zoom = std::make_shared<ZoomState>();
    ZoomState& z = *zoom;
    z.file = state_;
    if (chrom1.index <= chrom2.index) {
        z.c1 = chrom1.index;
        z.c2 = chrom2.index;
        z.numBins1 = static_cast<int32_t>(chrom1.length / resolution);
        z.numBins2 = static_cast<int32_t>(chrom2.length / resolution);
    } else {
        z.c1 = chrom2.index;
        z.c2 = chrom1.index;
        z.numBins1 = static_cast<int32_t>(chrom2.length / resolution);
        z.numBins2 = static_cast<int32_t>(chrom1.length / resolution);
    }
    z.isIntra = z.c1 == z.c2;
    z.matrixType = matrixType;
    z.norm = norm;
    z.unit = unit;
    z.resolution = resolution;
    const bool wantsExpected = matrixType == "oe" || matrixType == "expected";

    const std::string key = detail::pairKey(z.c1, z.c2);
    const auto matrix = file.matrices.find(key);
    if (matrix == file.matrices.end()) {
        z.message = "File doesn't have the given chr_chr map " + key;
        return MatrixZoomData(zoom);
    }

    const bool needsFooter = !((matrixType == "observed" && norm == detail::kNone) ||
                               (wantsExpected && norm == detail::kNone && z.c1 != z.c2));
    if (needsFooter) {
        const auto collect = [&](const std::vector<detail::ExpectedEntry>& section, bool normalized) {
            for (const auto& entry : section) {
                const bool store = z.isIntra && wantsExpected &&
                                   (normalized ? entry.key.normalization == norm : norm == detail::kNone) &&
                                   entry.key.unit == unit && entry.key.binSize == resolution;
                if (!store) {
                    continue;
                }
                z.expectedValues = file.readValues(entry.location.valuesOffset, entry.location.nValues);
                for (const auto& factor : file.readFactors(entry.location.factorsOffset, entry.location.nFactors)) {
                    if (factor.first == z.c1) {
                        for (double& value : z.expectedValues) {
                            value = value / factor.second;
                        }
                    }
                }
            }
        };
        collect(file.expectedNone, false);
        if (z.isIntra && wantsExpected && norm == detail::kNone) {
            if (z.expectedValues.empty()) {
                z.message = "File did not contain expected values vectors at " + describe(resolution, unit);
                return MatrixZoomData(zoom);
            }
        } else {
            collect(file.expectedNormalized, true);
            if (z.isIntra && wantsExpected && z.expectedValues.empty()) {
                z.message =
                    "File did not contain normalized expected values vectors at " + describe(resolution, unit);
                return MatrixZoomData(zoom);
            }
            if (norm != detail::kNone) {
                std::optional<int64_t> position1;
                std::optional<int64_t> position2;
                for (const auto& entry : file.normIndex) {
                    if (entry.normalization == norm && entry.unit == unit && entry.resolution == resolution) {
                        if (entry.chrIndex == z.c1) {
                            position1 = entry.position;
                        }
                        if (entry.chrIndex == z.c2) {
                            position2 = entry.position;
                        }
                    }
                }
                if (!position1 || !position2) {
                    throw HicError("File did not contain " + norm +
                                   " normalization vectors for one or both chromosomes at " +
                                   describe(resolution, unit));
                }
                z.c1Norm = file.readNormVectorAt(*position1);
                z.c2Norm = z.isIntra ? z.c1Norm : file.readNormVectorAt(*position2);
            }
        }
    }
    z.foundFooter = true;

    // readMatrix and readMatrixZoomData
    detail::SequentialReader in(file.file, matrix->second.position, 1 << 16);
    in.get<int32_t>();
    in.get<int32_t>();
    const int32_t nRes = in.get<int32_t>();
    bool found = false;
    for (int32_t i = 0; i < nRes && !found; ++i) {
        ZoomHeader header;
        header.unit = in.cstr();
        header.zoomIndex = in.get<int32_t>();
        header.sumCounts = in.get<float>();
        header.occupiedCellCount = in.get<float>();
        header.percent5 = in.get<float>();
        header.percent95 = in.get<float>();
        header.binSize = in.get<int32_t>();
        header.blockBinCount = in.get<int32_t>();
        header.blockColumnCount = in.get<int32_t>();
        const int32_t nBlocks = in.get<int32_t>();
        if (header.unit == unit && header.binSize == resolution) {
            found = true;
            z.header = header;
            z.blockList.reserve(static_cast<size_t>(std::max(nBlocks, 0)));
            for (int32_t b = 0; b < nBlocks; ++b) {
                BlockIndexEntry entry;
                entry.number = in.get<int32_t>();
                entry.position = in.get<int64_t>();
                entry.size = in.get<int32_t>();
                z.blockMap[entry.number] = entry;
                z.blockList.push_back(entry);
            }
        } else {
            in.skip(16 * static_cast<int64_t>(std::max(nBlocks, 0)));
        }
    }
    if (!found) {
        throw HicError("Error finding block data: no " + describe(resolution, unit) + " zoom for chromosome pair " +
                       key);
    }
    if (!z.isIntra) {
        // hicstraw divides the float sum by int bin counts in float arithmetic.
        const float average = (z.header.sumCounts / static_cast<float>(z.numBins1)) / static_cast<float>(z.numBins2);
        z.avgCount = static_cast<double>(average);
    }
    return MatrixZoomData(zoom);
}

// --------------------------------------------------------- MatrixZoomData

MatrixZoomData::MatrixZoomData(std::shared_ptr<const ZoomState> state) : state_(std::move(state)) {}

bool MatrixZoomData::found() const {
    return state_->foundFooter;
}

const std::string& MatrixZoomData::message() const {
    return state_->message;
}

int32_t MatrixZoomData::chr1Index() const {
    return state_->c1;
}

int32_t MatrixZoomData::chr2Index() const {
    return state_->c2;
}

bool MatrixZoomData::isIntra() const {
    return state_->isIntra;
}

int32_t MatrixZoomData::resolution() const {
    return state_->resolution;
}

const ZoomHeader& MatrixZoomData::zoomHeader() const {
    return state_->header;
}

const std::vector<BlockIndexEntry>& MatrixZoomData::blockIndex() const {
    return state_->blockList;
}

std::vector<ContactRecord> MatrixZoomData::readBlock(const BlockIndexEntry& entry) const {
    std::vector<char> scratch;
    std::vector<ContactRecord> records;
    state_->file->decodeBlock(entry, scratch, records);
    return records;
}

void MatrixZoomData::forEachBlock(
    const std::function<void(const BlockIndexEntry&, std::vector<ContactRecord>&)>& visit, int threads) const {
    std::vector<BlockIndexEntry> order = state_->blockList;
    std::sort(order.begin(), order.end(),
              [](const BlockIndexEntry& a, const BlockIndexEntry& b) { return a.number < b.number; });
    if (threads <= 1 || order.size() <= 1) {
        std::vector<char> scratch;
        std::vector<ContactRecord> records;
        for (const auto& entry : order) {
            records.clear();
            state_->file->decodeBlock(entry, scratch, records);
            visit(entry, records);
        }
        return;
    }
    const size_t window = static_cast<size_t>(threads) * 8;
    for (size_t begin = 0; begin < order.size(); begin += window) {
        const size_t end = std::min(order.size(), begin + window);
        std::vector<std::vector<ContactRecord>> decoded(end - begin);
        std::vector<std::exception_ptr> errors(static_cast<size_t>(threads));
        std::mutex lock;
        size_t next = begin;
        auto work = [&](size_t worker) {
            std::vector<char> scratch;
            try {
                while (true) {
                    size_t item = 0;
                    {
                        std::lock_guard<std::mutex> guard(lock);
                        if (next >= end) {
                            return;
                        }
                        item = next++;
                    }
                    state_->file->decodeBlock(order[item], scratch, decoded[item - begin]);
                }
            } catch (...) {
                errors[worker] = std::current_exception();
            }
        };
        std::vector<std::thread> pool;
        const size_t workers = std::min(static_cast<size_t>(threads), end - begin);
        for (size_t t = 1; t < workers; ++t) {
            pool.emplace_back(work, t);
        }
        work(0);
        for (auto& thread : pool) {
            thread.join();
        }
        for (const auto& error : errors) {
            if (error) {
                std::rethrow_exception(error);
            }
        }
        for (size_t item = begin; item < end; ++item) {
            visit(order[item], decoded[item - begin]);
        }
    }
}

// getRecords
std::vector<ContactRecord> MatrixZoomData::getRecords(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) const {
    const ZoomState& z = *state_;
    std::vector<ContactRecord> records;
    if (!z.foundFooter) {
        return records;
    }
    const int64_t orig[4] = {gx0, gx1, gy0, gy1};
    int64_t region[4];
    for (int q = 0; q < 4; ++q) {
        region[q] = orig[q] / z.resolution;
    }
    const bool normalized = z.norm != "NONE";
    const bool oe = z.matrixType == "oe";
    const bool expected = z.matrixType == "expected";
    std::vector<char> scratch;
    std::vector<ContactRecord> block;
    for (const int32_t number : z.blockNumbers(region)) {
        const auto it = z.blockMap.find(number);
        if (it == z.blockMap.end() || it->second.size <= 0) {
            continue;
        }
        block.clear();
        z.file->decodeBlock(it->second, scratch, block);
        for (const ContactRecord& rec : block) {
            const int64_t x = static_cast<int64_t>(rec.binX) * z.resolution;
            const int64_t y = static_cast<int64_t>(rec.binY) * z.resolution;
            const bool inside = (x >= orig[0] && x <= orig[1] && y >= orig[2] && y <= orig[3]) ||
                                (z.isIntra && y >= orig[0] && y <= orig[1] && x >= orig[2] && x <= orig[3]);
            if (!inside) {
                continue;
            }
            float c = rec.counts;
            if (normalized) {
                c = static_cast<float>(c / (normAt(z.c1Norm, rec.binX) * normAt(z.c2Norm, rec.binY)));
            }
            if (oe || expected) {
                double denominator = z.avgCount;
                if (z.isIntra) {
                    const auto distance = static_cast<size_t>(std::abs(y - x) / z.resolution);
                    denominator = z.expectedValues[std::min(z.expectedValues.size() - 1, distance)];
                }
                c = oe ? static_cast<float>(c / denominator) : static_cast<float>(denominator);
            }
            if (!std::isnan(c) && !std::isinf(c)) {
                records.push_back(ContactRecord{static_cast<int32_t>(x), static_cast<int32_t>(y), c});
            }
        }
    }
    return records;
}

// getRecordsAsMatrix
FloatMatrix MatrixZoomData::getRecordsAsMatrix(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1) const {
    const ZoomState& z = *state_;
    const std::vector<ContactRecord> records = getRecords(gx0, gx1, gy0, gy1);
    FloatMatrix matrix;
    if (records.empty()) {
        matrix.rows = 1;
        matrix.cols = 1;
        matrix.values.assign(1, 0.0f);
        return matrix;
    }
    const int64_t originR = gx0 / z.resolution;
    const int64_t endR = gx1 / z.resolution;
    const int64_t originC = gy0 / z.resolution;
    const int64_t endC = gy1 / z.resolution;
    const auto numRows = static_cast<int32_t>(endR - originR + 1);
    const auto numCols = static_cast<int32_t>(endC - originC + 1);
    matrix.rows = std::max(numRows, 0);
    matrix.cols = std::max(numCols, 0);
    matrix.values.assign(static_cast<size_t>(matrix.rows * matrix.cols), 0.0f);
    const auto inRange = [&](int32_t r, int32_t c) { return 0 <= r && r < numRows && 0 <= c && c < numCols; };
    for (const ContactRecord& cr : records) {
        if (std::isnan(cr.counts) || std::isinf(cr.counts)) {
            continue;
        }
        auto r = static_cast<int32_t>(cr.binX / z.resolution - originR);
        auto c = static_cast<int32_t>(cr.binY / z.resolution - originC);
        if (inRange(r, c)) {
            matrix.values[static_cast<size_t>(static_cast<int64_t>(r) * matrix.cols + c)] = cr.counts;
        }
        if (z.isIntra) {
            r = static_cast<int32_t>(cr.binY / z.resolution - originR);
            c = static_cast<int32_t>(cr.binX / z.resolution - originC);
            if (inRange(r, c)) {
                matrix.values[static_cast<size_t>(static_cast<int64_t>(r) * matrix.cols + c)] = cr.counts;
            }
        }
    }
    return matrix;
}

// getNormVector
std::vector<double> MatrixZoomData::getNormVector(int32_t index) const {
    if (index == state_->c1) {
        return state_->c1Norm;
    }
    if (index == state_->c2) {
        return state_->c2Norm;
    }
    return {};
}

std::vector<double> MatrixZoomData::getExpectedValues() const {
    return state_->expectedValues;
}

// getNumberOfTotalRecords
int64_t MatrixZoomData::getNumberOfTotalRecords() const {
    const ZoomState& z = *state_;
    if (!z.foundFooter) {
        return 0;
    }
    const int64_t region[4] = {0, z.numBins1, 0, z.numBins2};
    int64_t total = 0;
    std::vector<char> compressed;
    std::vector<char> scratch;
    for (const int32_t number : z.blockNumbers(region)) {
        const auto it = z.blockMap.find(number);
        if (it == z.blockMap.end() || it->second.size <= 0) {
            continue;
        }
        compressed.resize(static_cast<size_t>(it->second.size));
        z.file->file.readAt(it->second.position, compressed.data(), compressed.size());
        detail::inflateZlib(compressed.data(), compressed.size(), scratch);
        detail::MemReader in(scratch.data(), scratch.size());
        total += in.get<int32_t>();
    }
    return total;
}

}  // namespace hicfilecpp
