// Writing .hic files: Juicer tools pre (juicebox.tools.utils.original.
// Preprocessor, MatrixPP, MatrixZoomDataPP) and addNorm
// (juicebox.tools.utils.norm.NormalizationVectorUpdater, NormVectorUpdater),
// following release 1.22.01 for version 8 and 2.20.00 for version 9.
//
// Layout of a written file: the header; for every chromosome pair with
// pixels, the blocks of every resolution followed by the pair's matrix header
// and block indexes; the whole-genome matrix the same way; the footer with the
// matrix index, the raw expected values, the normalized expected values, the
// normalization vector index and the vectors. Juicer places each matrix
// header before its blocks; readers locate both through the index, so the
// order does not matter to them.

#include "hicfilecpp/writer.hpp"

#include <algorithm>
#include <climits>
#include <functional>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "codec.hpp"
#include "hicfilecpp/hicfile.hpp"
#include "hicfilecpp/version.hpp"
#include "norms.hpp"
#include "output.hpp"
#include "state.hpp"

namespace hicfilecpp {

namespace detail {

namespace {

const std::vector<std::string> kNormOrder{"VC", "VC_SQRT", "KR", "SCALE"};
constexpr int32_t kBlockSize = 1000;  // Preprocessor.BLOCK_SIZE
constexpr int32_t kMaxSqrt = 46340;   // (int) Math.sqrt(Integer.MAX_VALUE)

struct Layout {
    int32_t binSize = 0;
    int32_t blockBinCount = 0;
    int32_t blockColumnCount = 0;
    bool intra = false;
    int32_t version = 9;
};

// MatrixPP.getNumColumnsFromNumBins
int32_t numColumns(int32_t version, int32_t nBins, int32_t binSize, bool intra) {
    int32_t nColumns = nBins / kBlockSize + 1;
    if (version == 8) {
        if (nColumns > std::sqrt(static_cast<double>(INT32_MAX))) {
            nColumns = kMaxSqrt - 1;
        }
        return nColumns;
    }
    const int32_t cutoff = intra ? 500 : 5000;
    if (binSize < cutoff) {
        const int64_t numerator = static_cast<int64_t>(nBins) * binSize;
        const int64_t denominator = static_cast<int64_t>(kBlockSize) * cutoff;
        nColumns = static_cast<int32_t>(numerator / denominator) + 1;
    }
    return std::min(nColumns, kMaxSqrt - 1);
}

// The MatrixPP and MatrixZoomDataPP constructors.
Layout pairLayout(int32_t version, int64_t length1, int64_t length2, int32_t binSize, bool intra) {
    Layout layout;
    layout.binSize = binSize;
    layout.intra = intra;
    layout.version = version;
    const int64_t length = std::max(length1, length2);
    const auto nBins = static_cast<int32_t>(length / binSize + 1);
    layout.blockColumnCount = numColumns(version, nBins, binSize, intra);
    layout.blockBinCount = nBins / layout.blockColumnCount + 1;
    return layout;
}

// Preprocessor.getInitialGenomeWideMatrixPP
Layout wholeGenomeLayout(int32_t version, int32_t genomeLengthKb) {
    Layout layout;
    layout.binSize = std::max(genomeLengthKb / 500, 1);
    layout.intra = true;
    layout.version = version;
    const int32_t nBinsX = genomeLengthKb / layout.binSize + 1;
    layout.blockColumnCount = nBinsX / kBlockSize + 1;
    layout.blockBinCount = nBinsX / layout.blockColumnCount + 1;
    return layout;
}

// MatrixZoomDataPP.incrementCount; version 9 numbers intra-chromosomal blocks
// by their distance from the diagonal (LogDepth, base 2).
int32_t blockNumber(const Layout& layout, int32_t x, int32_t y) {
    if (layout.intra && layout.version > 8) {
        const double v = std::abs(x - y) / std::sqrt(2) / layout.blockBinCount;
        const auto depth = static_cast<int32_t>(std::log(1 + v) / std::log(2));
        const int32_t position = (x + y) / 2 / layout.blockBinCount;
        return depth * layout.blockColumnCount + position;
    }
    return layout.blockColumnCount * (y / layout.blockBinCount) + x / layout.blockBinCount;
}

struct Entry {
    int32_t block;
    int32_t x;
    int32_t y;
    float count;
};

// A stable LSD radix sort on the block number.
void sortByBlock(std::vector<Entry>& entries) {
    std::vector<Entry> scratch(entries.size());
    std::vector<size_t> counts(size_t{1} << 16);
    for (int pass = 0; pass < 2; ++pass) {
        const int shift = pass * 16;
        std::fill(counts.begin(), counts.end(), 0);
        for (const Entry& e : entries) {
            counts[(static_cast<uint32_t>(e.block) >> shift) & 0xFFFF]++;
        }
        size_t total = 0;
        for (auto& count : counts) {
            const size_t n = count;
            count = total;
            total += n;
        }
        for (const Entry& e : entries) {
            scratch[counts[(static_cast<uint32_t>(e.block) >> shift) & 0xFFFF]++] = e;
        }
        entries.swap(scratch);
    }
}

// MatrixZoomDataPP.writeBlock, before compression. Records are sorted by y,
// then x. Version 8 picks the smaller of the list-of-rows and dense
// representations as 1.22.01 does; version 9 always writes list of rows, as
// 2.20.00 does.
void encodeBlock(int32_t version, const std::vector<NormRecord>& records, ByteWriter& out) {
    out.bytes.clear();
    out.bytes.reserve(records.size() * 6 + 32);
    int32_t binXOffset = INT32_MAX;
    int32_t binYOffset = INT32_MAX;
    int32_t binXMax = 0;
    int32_t binYMax = 0;
    bool isInteger = true;
    float maxCounts = 0;
    size_t rows = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        const NormRecord& r = records[i];
        binXOffset = std::min(binXOffset, r.x);
        binYOffset = std::min(binYOffset, r.y);
        binXMax = std::max(binXMax, r.x);
        binYMax = std::max(binYMax, r.y);
        isInteger = isInteger && std::floor(r.counts) == r.counts;
        maxCounts = std::max(r.counts, maxCounts);
        if (i == 0 || r.y != records[i - 1].y) {
            rows++;
        }
    }
    const bool useShort = isInteger && maxCounts < 32767;
    const int64_t valueSize = useShort ? 2 : 4;
    const int64_t lorSize = static_cast<int64_t>(rows) * 4 + static_cast<int64_t>(records.size()) * valueSize;
    auto putValue = [&](float counts) {
        if (useShort) {
            out.put<int16_t>(static_cast<int16_t>(counts));
        } else {
            out.put<float>(counts);
        }
    };
    out.put<int32_t>(static_cast<int32_t>(records.size()));
    out.put<int32_t>(binXOffset);
    out.put<int32_t>(binYOffset);
    out.put<int8_t>(useShort ? 0 : 1);

    bool rowShort = true;
    bool colShort = true;
    if (version > 8) {
        colShort = binXMax - binXOffset + 1 < 32767;
        rowShort = binYMax - binYOffset + 1 < 32767;
        out.put<int8_t>(colShort ? 0 : 1);
        out.put<int8_t>(rowShort ? 0 : 1);
    } else {
        const auto w = static_cast<int16_t>(binXMax - binXOffset + 1);
        const NormRecord& last = records.back();
        const int64_t nDensePts = static_cast<int64_t>(last.y - binYOffset) * w + (last.x - binXOffset) + 1;
        if (lorSize >= nDensePts * valueSize) {
            out.put<int8_t>(2);
            out.put<int32_t>(static_cast<int32_t>(nDensePts));
            out.put<int16_t>(w);
            int64_t lastIdx = 0;
            for (const NormRecord& r : records) {
                const int64_t idx = static_cast<int64_t>(r.y - binYOffset) * w + (r.x - binXOffset);
                for (int64_t i = lastIdx; i < idx; ++i) {
                    if (useShort) {
                        out.put<int16_t>(INT16_MIN);
                    } else {
                        out.put<float>(std::numeric_limits<float>::quiet_NaN());
                    }
                }
                putValue(r.counts);
                lastIdx = idx + 1;
            }
            return;
        }
    }
    out.put<int8_t>(1);
    auto putRowIndex = [&](int64_t value) {
        if (rowShort) {
            out.put<int16_t>(static_cast<int16_t>(value));
        } else {
            out.put<int32_t>(static_cast<int32_t>(value));
        }
    };
    auto putColIndex = [&](int64_t value) {
        if (colShort) {
            out.put<int16_t>(static_cast<int16_t>(value));
        } else {
            out.put<int32_t>(static_cast<int32_t>(value));
        }
    };
    putRowIndex(static_cast<int64_t>(rows));
    for (size_t begin = 0; begin < records.size();) {
        size_t end = begin + 1;
        while (end < records.size() && records[end].y == records[begin].y) {
            end++;
        }
        putRowIndex(records[begin].y - binYOffset);
        putColIndex(static_cast<int64_t>(end - begin));
        for (size_t i = begin; i < end; ++i) {
            putColIndex(records[i].x - binXOffset);
            putValue(records[i].counts);
        }
        begin = end;
    }
}

struct EncodedBlock {
    int32_t number = 0;
    std::vector<char> bytes;
    std::vector<NormRecord> records;
};

// Sorts the entries of one zoom level into blocks, merges repeated pixels in
// input order, encodes and compresses the blocks on the pool and writes them
// in ascending block number. With `keep`, also returns the merged records in
// that order, which is the order Juicer reads them back for normalization.
std::vector<BlockIndexEntry> writeBlocks(OutputFile& out, ThreadPool& pool, int32_t version, int level,
                                         std::vector<Entry>& entries, std::vector<NormRecord>* keep) {
    sortByBlock(entries);
    std::vector<std::pair<size_t, size_t>> ranges;
    for (size_t i = 0; i < entries.size();) {
        size_t j = i + 1;
        while (j < entries.size() && entries[j].block == entries[i].block) {
            j++;
        }
        ranges.emplace_back(i, j);
        i = j;
    }
    std::vector<BlockIndexEntry> index;
    index.reserve(ranges.size());
    const size_t window = static_cast<size_t>(pool.size()) * 32;
    std::vector<EncodedBlock> encoded;
    for (size_t begin = 0; begin < ranges.size(); begin += window) {
        const size_t end = std::min(ranges.size(), begin + window);
        encoded.assign(end - begin, EncodedBlock{});
        pool.parallelFor(end - begin, [&](size_t i) {
            const auto [first, last] = ranges[begin + i];
            std::stable_sort(entries.begin() + static_cast<std::ptrdiff_t>(first),
                             entries.begin() + static_cast<std::ptrdiff_t>(last),
                             [](const Entry& a, const Entry& b) { return a.y != b.y ? a.y < b.y : a.x < b.x; });
            std::vector<NormRecord> records;
            records.reserve(last - first);
            for (size_t k = first; k < last; ++k) {
                const Entry& e = entries[k];
                if (!records.empty() && records.back().x == e.x && records.back().y == e.y) {
                    records.back().counts += e.count;
                } else {
                    records.push_back(NormRecord{e.x, e.y, e.count});
                }
            }
            ByteWriter raw;
            encodeBlock(version, records, raw);
            EncodedBlock& block = encoded[i];
            block.number = entries[first].block;
            deflateZlib(raw.bytes.data(), raw.bytes.size(), level, block.bytes);
            if (keep != nullptr) {
                block.records = std::move(records);
            }
        });
        for (auto& block : encoded) {
            const int64_t position = out.position();
            out.write(block.bytes);
            index.push_back(BlockIndexEntry{block.number, position, static_cast<int32_t>(block.bytes.size())});
            if (keep != nullptr) {
                keep->insert(keep->end(), block.records.begin(), block.records.end());
            }
        }
    }
    return index;
}

struct ZoomWritten {
    int32_t zoomIndex = 0;
    double sum = 0;
    Layout layout;
    std::vector<BlockIndexEntry> index;
};

// Preprocessor.writeMatrix and writeZoomHeader. Juicer writes the header
// before the blocks, so the stored sum holds the contacts (off-diagonal
// intra-chromosomal ones twice) and the cell count and percentiles are 0.
MatrixEntry writeMatrixHeader(OutputFile& out, int32_t chr1, int32_t chr2, const std::vector<ZoomWritten>& zooms) {
    ByteWriter m;
    m.put<int32_t>(chr1);
    m.put<int32_t>(chr2);
    m.put<int32_t>(static_cast<int32_t>(zooms.size()));
    for (const auto& zoom : zooms) {
        m.cstr("BP");
        m.put<int32_t>(zoom.zoomIndex);
        m.put<float>(static_cast<float>(zoom.sum));
        m.put<float>(0.0f);
        m.put<float>(0.0f);
        m.put<float>(0.0f);
        m.put<int32_t>(zoom.layout.binSize);
        m.put<int32_t>(zoom.layout.blockBinCount);
        m.put<int32_t>(zoom.layout.blockColumnCount);
        m.put<int32_t>(static_cast<int32_t>(zoom.index.size()));
        for (const auto& entry : zoom.index) {
            m.put<int32_t>(entry.number);
            m.put<int64_t>(entry.position);
            m.put<int32_t>(entry.size);
        }
    }
    MatrixEntry entry;
    entry.position = out.position();
    entry.size = static_cast<int32_t>(m.size());
    out.write(m);
    return entry;
}

void putValues(ByteWriter& w, int32_t version, const std::vector<double>& values) {
    if (version > 8) {
        w.put<int64_t>(static_cast<int64_t>(values.size()));
        for (const double v : values) {
            w.put<float>(static_cast<float>(v));
        }
    } else {
        w.put<int32_t>(static_cast<int32_t>(values.size()));
        for (const double v : values) {
            w.put<double>(v);
        }
    }
}

void putExpected(ByteWriter& w, int32_t version, ExpectedValueCalculation& ev, bool withType) {
    ev.computeDensity();
    if (withType) {
        w.cstr(ev.type());
    }
    w.cstr("BP");
    w.put<int32_t>(ev.gridSize());
    putValues(w, version, ev.densityAvg());
    w.put<int32_t>(static_cast<int32_t>(ev.chrScaleFactors().size()));
    for (const auto& [chr, factor] : ev.chrScaleFactors()) {
        w.put<int32_t>(chr);
        if (version > 8) {
            w.put<float>(static_cast<float>(factor));
        } else {
            w.put<double>(factor);
        }
    }
}

void validateNormalizations(const std::vector<std::string>& normalizations) {
    std::set<std::string> seen;
    for (const auto& norm : normalizations) {
        if (std::find(kNormOrder.begin(), kNormOrder.end(), norm) == kNormOrder.end()) {
            throw HicError("unsupported normalization " + norm + "; hicfilecpp computes VC, VC_SQRT, KR and SCALE");
        }
        if (!seen.insert(norm).second) {
            throw HicError("normalization " + norm + " is listed twice");
        }
    }
}

// NormalizationVectorUpdater: the vectors of every chromosome at every
// resolution and the normalized expected values they give.
class NormBuilder {
public:
    NormBuilder(int32_t version, std::vector<int64_t> lengths, std::vector<int32_t> resolutions,
                const std::vector<std::string>& requested, std::vector<std::string> provided = {},
                std::function<std::vector<double>(const std::string&, int32_t, int32_t)> provider = {})
        : version_(version),
          lengths_(std::move(lengths)),
          resolutions_(std::move(resolutions)),
          names_(kNormOrder),
          provided_(std::move(provided)),
          provider_(std::move(provider)) {
        for (const auto& type : kNormOrder) {
            if (std::find(requested.begin(), requested.end(), type) != requested.end()) {
                types_.push_back(type);
            }
        }
        // Provided types are indexed after the computed ones, in the order given.
        for (const auto& name : provided_) {
            if (std::find(names_.begin(), names_.end(), name) == names_.end()) {
                names_.push_back(name);
            }
        }
    }

    bool active() const { return !types_.empty() || !provided_.empty(); }

    // vectorLength is the bin count of the matrix's grid axis, which Juicer
    // takes from the zoom level's block layout (MatrixZoomData: blockBinCount
    // times blockColumnCount), not from the chromosome length.
    void addChromosome(int32_t chrIndex, int32_t resolution, const std::vector<NormRecord>& records,
                       int64_t vectorLength) {
        if (records.empty() || !active()) {
            return;
        }
        if (version_ > 8) {
            compute<float>(chrIndex, resolution, records, vectorLength);
            provide<float>(chrIndex, resolution, records, vectorLength);
        } else {
            compute<double>(chrIndex, resolution, records, vectorLength);
            provide<double>(chrIndex, resolution, records, vectorLength);
        }
    }

    // Writes the normalized expected values, the vector index and the vectors
    // at the output's position; returns the index's position and size.
    std::pair<int64_t, int64_t> write(OutputFile& out) {
        ByteWriter e;
        int32_t count = 0;
        for (auto& [key, ev] : expected_) {
            count += ev.hasData() ? 1 : 0;
        }
        e.put<int32_t>(count);
        for (auto& [key, ev] : expected_) {
            if (ev.hasData()) {
                putExpected(e, version_, ev, true);
            }
        }
        out.write(e);
        std::stable_sort(vectors_.begin(), vectors_.end(), [](const StoredVector& a, const StoredVector& b) {
            return std::tie(a.zoomOrder, a.chr, a.typeOrder) < std::tie(b.zoomOrder, b.chr, b.typeOrder);
        });
        const int64_t indexPosition = out.position();
        int64_t indexSize = 4;
        for (const auto& v : vectors_) {
            indexSize += static_cast<int64_t>(names_[static_cast<size_t>(v.typeOrder)].size()) + 1 + 4 + 3 + 4 +
                         8 + (version_ > 8 ? 8 : 4);
        }
        ByteWriter index;
        index.put<int32_t>(static_cast<int32_t>(vectors_.size()));
        int64_t position = indexPosition + indexSize;
        for (const auto& v : vectors_) {
            const auto n = static_cast<int64_t>(v.values.size());
            const int64_t size = version_ > 8 ? 8 + 4 * n : 4 + 8 * n;
            index.cstr(names_[static_cast<size_t>(v.typeOrder)]);
            index.put<int32_t>(v.chr);
            index.cstr("BP");
            index.put<int32_t>(v.resolution);
            index.put<int64_t>(position);
            if (version_ > 8) {
                index.put<int64_t>(size);
            } else {
                index.put<int32_t>(static_cast<int32_t>(size));
            }
            position += size;
        }
        out.write(index);
        for (const auto& v : vectors_) {
            ByteWriter w;
            putValues(w, version_, v.values);
            out.write(w);
        }
        return {indexPosition, indexSize};
    }

private:
    struct StoredVector {
        int32_t zoomOrder = 0;
        int32_t chr = 0;
        int32_t typeOrder = 0;
        int32_t resolution = 0;
        std::vector<double> values;
    };

    int32_t zoomOrder(int32_t resolution) const {
        return static_cast<int32_t>(std::find(resolutions_.begin(), resolutions_.end(), resolution) -
                                    resolutions_.begin());
    }
    int32_t typeOrder(const std::string& type) const {
        return static_cast<int32_t>(std::find(names_.begin(), names_.end(), type) - names_.begin());
    }
    bool wants(const std::string& type) const {
        return std::find(types_.begin(), types_.end(), type) != types_.end();
    }
    ExpectedValueCalculation& expected(int32_t resolution, const std::string& type) {
        const auto key = std::make_pair(zoomOrder(resolution), typeOrder(type));
        auto it = expected_.find(key);
        if (it == expected_.end()) {
            it = expected_.emplace(key, ExpectedValueCalculation(lengths_, resolution, type)).first;
        }
        return it->second;
    }

    template <class T>
    void compute(int32_t chr, int32_t resolution, const std::vector<NormRecord>& records, int64_t size) {
        const bool vc = wants("VC");
        const bool vcSqrt = wants("VC_SQRT");
        if (vc || vcSqrt) {
            std::vector<T> values = computeVC<T>(records, size);
            std::vector<T> roots(values.size(), T(0));
            if (vcSqrt) {
                for (size_t i = 0; i < values.size(); ++i) {
                    roots[i] = static_cast<T>(std::sqrt(static_cast<double>(values[i])));
                }
            }
            if (vc) {
                update<T>(chr, resolution, "VC", std::move(values), records);
            }
            if (vcSqrt) {
                update<T>(chr, resolution, "VC_SQRT", std::move(roots), records);
            }
        }
        if (wants("KR")) {
            update<T>(chr, resolution, "KR", computeKR<T>(records, size, version_ > 8), records);
        }
        if (wants("SCALE")) {
            std::vector<T> scale = computeScale<T>(records, size);
            if (!scale.empty()) {
                update<T>(chr, resolution, "SCALE", std::move(scale), records);
            }
        }
    }

    // The caller's vectors, stored as given (writeHicFile's chromosome indexes
    // start at 1 for the whole-genome pseudo-chromosome; the caller's at 0).
    template <class T>
    void provide(int32_t chr, int32_t resolution, const std::vector<NormRecord>& records, int64_t size) {
        for (const auto& name : provided_) {
            const std::vector<double> given = provider_(name, chr - 1, resolution);
            if (given.empty()) {
                continue;
            }
            const bool vcType = name == "VC" || name == "VC_SQRT" ||
                                (name.size() > 3 && name.compare(name.size() - 3, 3, "_VC") == 0) ||
                                (name.size() > 8 && name.compare(name.size() - 8, 8, "_VC_SQRT") == 0);
            const T pad = vcType ? T(0) : std::numeric_limits<T>::quiet_NaN();
            std::vector<T> vector(static_cast<size_t>(std::max<int64_t>(size, 0)), pad);
            const size_t n = std::min(vector.size(), given.size());
            for (size_t i = 0; i < n; ++i) {
                vector[i] = static_cast<T>(given[i]);
            }
            store<T>(chr, resolution, name, std::move(vector), records);
        }
    }

    // NormalizationVectorUpdater.updateExpectedValueCalculationForChr
    template <class T>
    void update(int32_t chr, int32_t resolution, const std::string& type, std::vector<T> vector,
                const std::vector<NormRecord>& records) {
        const double factor = sumFactor(records, vector);
        for (auto& value : vector) {
            value = static_cast<T>(value * factor);
        }
        store<T>(chr, resolution, type, std::move(vector), records);
    }

    template <class T>
    void store(int32_t chr, int32_t resolution, const std::string& type, std::vector<T> vector,
               const std::vector<NormRecord>& records) {
        addDistancesFromRecords(expected(resolution, type), chr, records, vector);
        StoredVector stored;
        stored.zoomOrder = zoomOrder(resolution);
        stored.chr = chr;
        stored.typeOrder = typeOrder(type);
        stored.resolution = resolution;
        stored.values.assign(vector.begin(), vector.end());
        vectors_.push_back(std::move(stored));
    }

    int32_t version_;
    std::vector<int64_t> lengths_;
    std::vector<int32_t> resolutions_;
    std::vector<std::string> types_;
    std::vector<std::string> names_;
    std::vector<std::string> provided_;
    std::function<std::vector<double>(const std::string&, int32_t, int32_t)> provider_;
    std::map<std::pair<int32_t, int32_t>, ExpectedValueCalculation> expected_;
    std::vector<StoredVector> vectors_;
};

void validate(const WriteOptions& options) {
    if (options.version == 6 || options.version == 7) {
        throw HicError("writing .hic version " + std::to_string(options.version) +
                       " is not supported: no Juicer tools release that writes it can be obtained to "
                       "validate against; hicfilecpp reads versions 6 to 9 and writes 8 and 9");
    }
    if (options.version != 8 && options.version != 9) {
        throw HicError("version must be 8 or 9, got " + std::to_string(options.version));
    }
    if (options.chromosomes.empty()) {
        throw HicError("no chromosomes given");
    }
    std::set<std::string> names;
    for (const auto& [name, length] : options.chromosomes) {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
        if (name.empty() || lower == "all") {
            throw HicError("invalid chromosome name '" + name + "'");
        }
        if (!names.insert(name).second) {
            throw HicError("chromosome " + name + " is listed twice");
        }
        if (length <= 0 || length > INT32_MAX) {
            throw HicError("chromosome " + name + " has length " + std::to_string(length) +
                           "; lengths must lie in 1 to 2147483647");
        }
    }
    if (options.resolutions.empty()) {
        throw HicError("no resolutions given");
    }
    std::set<int32_t> resolutions;
    for (const int32_t resolution : options.resolutions) {
        if (resolution <= 0 || !resolutions.insert(resolution).second) {
            throw HicError("resolutions must be positive and distinct");
        }
    }
    if (!options.sourceProvidesEveryResolution) {
        if (options.sourceResolution <= 0) {
            throw HicError("sourceResolution must be positive");
        }
        for (const int32_t resolution : options.resolutions) {
            if (resolution % options.sourceResolution != 0) {
                throw HicError("resolution " + std::to_string(resolution) + " is not a multiple of the source resolution " +
                               std::to_string(options.sourceResolution));
            }
        }
    }
    validateNormalizations(options.normalizations);
    std::set<std::string> provided;
    for (const auto& name : options.providedNormalizations) {
        if (name.empty() || name == "NONE" ||
            std::any_of(name.begin(), name.end(), [](unsigned char c) { return c <= ' ' || c >= 127; })) {
            throw HicError("invalid provided normalization name '" + name + "'");
        }
        if (!provided.insert(name).second) {
            throw HicError("provided normalization " + name + " is listed twice");
        }
        if (std::find(options.normalizations.begin(), options.normalizations.end(), name) !=
            options.normalizations.end()) {
            throw HicError("normalization " + name + " is both computed and provided");
        }
    }
    if (!options.providedNormalizations.empty() && !options.normVector) {
        throw HicError("providedNormalizations needs a normVector function");
    }
    if (options.threads < 1) {
        throw HicError("threads must be at least 1");
    }
    if (options.compressionLevel < -1 || options.compressionLevel > 9) {
        throw HicError("compressionLevel must lie in -1 to 9");
    }
}

}  // namespace

}  // namespace detail

void writeHicFile(const std::string& path, const WriteOptions& options, PixelSource& source) {
    using namespace detail;
    validate(options);
    const int32_t version = options.version;
    std::vector<int32_t> resolutions = options.resolutions;
    std::sort(resolutions.begin(), resolutions.end(), std::greater<>());
    const auto nChromosomes = static_cast<int32_t>(options.chromosomes.size());

    std::vector<int64_t> lengths(static_cast<size_t>(nChromosomes) + 1);
    std::vector<int64_t> offsets(static_cast<size_t>(nChromosomes) + 1, 0);
    int64_t genomeLength = 0;
    for (int32_t i = 0; i < nChromosomes; ++i) {
        offsets[static_cast<size_t>(i) + 1] = genomeLength;
        lengths[static_cast<size_t>(i) + 1] = options.chromosomes[static_cast<size_t>(i)].second;
        genomeLength += options.chromosomes[static_cast<size_t>(i)].second;
    }
    const auto genomeLengthKb = static_cast<int32_t>(genomeLength / 1000);
    lengths[0] = genomeLengthKb;

    OutputFile out(path, std::nullopt);
    ByteWriter header;
    header.cstr("HIC");
    header.put<int32_t>(version);
    const auto masterOffset = static_cast<int64_t>(header.size());
    header.put<int64_t>(0);
    header.cstr(options.genomeId);
    const auto nviOffset = static_cast<int64_t>(header.size());
    if (version > 8) {
        header.put<int64_t>(0);
        header.put<int64_t>(0);
    }
    header.put<int32_t>(static_cast<int32_t>(options.attributes.size()) + 1);
    header.cstr("software");
    header.cstr(options.software.empty() ? std::string("hicfilecpp ") + kVersion : options.software);
    for (const auto& [key, value] : options.attributes) {
        header.cstr(key);
        header.cstr(value);
    }
    header.put<int32_t>(nChromosomes + 1);
    header.cstr("All");
    if (version > 8) {
        header.put<int64_t>(genomeLengthKb);
    } else {
        header.put<int32_t>(genomeLengthKb);
    }
    for (const auto& [name, length] : options.chromosomes) {
        header.cstr(name);
        if (version > 8) {
            header.put<int64_t>(length);
        } else {
            header.put<int32_t>(static_cast<int32_t>(length));
        }
    }
    header.put<int32_t>(static_cast<int32_t>(resolutions.size()));
    for (const int32_t resolution : resolutions) {
        header.put<int32_t>(resolution);
    }
    header.put<int32_t>(0);
    out.write(header);

    ThreadPool pool(options.threads);
    std::vector<ExpectedValueCalculation> expectedNone;
    for (const int32_t resolution : resolutions) {
        expectedNone.emplace_back(lengths, resolution, "NONE");
    }
    NormBuilder norms(version, lengths, resolutions, options.normalizations, options.providedNormalizations,
                      options.normVector);

    const Layout wholeLayout = wholeGenomeLayout(version, genomeLengthKb);
    const int64_t wholeBins = genomeLengthKb / wholeLayout.binSize + 1;
    std::vector<float> wholeCells(static_cast<size_t>(wholeBins * wholeBins), 0.0f);
    std::vector<uint8_t> wholePresent(wholeCells.size(), 0);
    double wholeSum = 0;

    std::vector<std::pair<std::string, MatrixEntry>> matrices;
    std::vector<Entry> entries;
    std::vector<NormRecord> records;
    for (int32_t c1 = 0; c1 < nChromosomes; ++c1) {
        for (int32_t c2 = c1; c2 < nChromosomes; ++c2) {
            const bool intra = c1 == c2;
            const int64_t length1 = lengths[static_cast<size_t>(c1) + 1];
            const int64_t length2 = lengths[static_cast<size_t>(c2) + 1];
            std::vector<ZoomWritten> zooms;
            for (size_t z = 0; z < resolutions.size(); ++z) {
                const int32_t resolution = resolutions[z];
                const int32_t request =
                    options.sourceProvidesEveryResolution ? resolution : options.sourceResolution;
                const Layout layout = pairLayout(version, length1, length2, resolution, intra);
                const bool wholeGenomePass = z + 1 == resolutions.size();
                ExpectedValueCalculation& ev = expectedNone[z];
                entries.clear();
                double sum = 0;
                source.pixels(request, c1, c2, [&](const Pixel* pixels, size_t count) {
                    for (size_t k = 0; k < count; ++k) {
                        const Pixel& p = pixels[k];
                        const int64_t pos1 = static_cast<int64_t>(p.bin1) * request;
                        const int64_t pos2 = static_cast<int64_t>(p.bin2) * request;
                        if (p.bin1 < 0 || p.bin2 < 0 || pos1 >= length1 || pos2 >= length2) {
                            throw HicError("pixel (" + std::to_string(p.bin1) + ", " + std::to_string(p.bin2) +
                                           ") at " + std::to_string(request) + " bp lies outside chromosomes " +
                                           options.chromosomes[static_cast<size_t>(c1)].first + " and " +
                                           options.chromosomes[static_cast<size_t>(c2)].first);
                        }
                        if (!std::isfinite(p.count)) {
                            throw HicError("pixel count is not finite");
                        }
                        const float score = p.count;
                        auto x = static_cast<int32_t>(pos1 / resolution);
                        auto y = static_cast<int32_t>(pos2 / resolution);
                        sum += score;
                        if (intra) {
                            if (x > y) {
                                std::swap(x, y);
                            }
                            if (x != y) {
                                sum += score;
                            }
                            ev.addDistance(c1 + 1, x, y, score);
                        }
                        entries.push_back(Entry{blockNumber(layout, x, y), x, y, score});
                        if (wholeGenomePass) {
                            // Preprocessor.getGenomicPosition, in kilobases
                            const auto g1 = static_cast<int32_t>((offsets[static_cast<size_t>(c1) + 1] + pos1) / 1000);
                            const auto g2 = static_cast<int32_t>((offsets[static_cast<size_t>(c2) + 1] + pos2) / 1000);
                            wholeSum += score;
                            int32_t wx = g1 / wholeLayout.binSize;
                            int32_t wy = g2 / wholeLayout.binSize;
                            if (wx > wy) {
                                std::swap(wx, wy);
                            }
                            if (wx != wy) {
                                wholeSum += score;
                            }
                            const auto cell = static_cast<size_t>(wx * wholeBins + wy);
                            if (wholePresent[cell] != 0) {
                                wholeCells[cell] += score;
                            } else {
                                wholePresent[cell] = 1;
                                wholeCells[cell] = score;
                            }
                        }
                    }
                });
                if (entries.empty()) {
                    if (z == 0) {
                        break;
                    }
                    throw HicError("the source gave pixels for " + options.chromosomes[static_cast<size_t>(c1)].first +
                                   " and " + options.chromosomes[static_cast<size_t>(c2)].first +
                                   " at some resolutions only");
                }
                const bool keep = intra && norms.active();
                records.clear();
                ZoomWritten zoom;
                zoom.zoomIndex = static_cast<int32_t>(z);
                zoom.sum = sum;
                zoom.layout = layout;
                zoom.index = writeBlocks(out, pool, version, options.compressionLevel, entries,
                                         keep ? &records : nullptr);
                zooms.push_back(std::move(zoom));
                if (keep) {
                    norms.addChromosome(c1 + 1, resolution, records,
                                        static_cast<int64_t>(layout.blockBinCount) * layout.blockColumnCount);
                }
            }
            if (!zooms.empty()) {
                matrices.emplace_back(std::to_string(c1 + 1) + "_" + std::to_string(c2 + 1),
                                      writeMatrixHeader(out, c1 + 1, c2 + 1, zooms));
            }
        }
    }

    // The whole-genome matrix, Preprocessor.computeWholeGenomeMatrix.
    entries.clear();
    for (int64_t x = 0; x < wholeBins; ++x) {
        for (int64_t y = x; y < wholeBins; ++y) {
            const auto cell = static_cast<size_t>(x * wholeBins + y);
            if (wholePresent[cell] != 0) {
                entries.push_back(Entry{blockNumber(wholeLayout, static_cast<int32_t>(x), static_cast<int32_t>(y)),
                                        static_cast<int32_t>(x), static_cast<int32_t>(y), wholeCells[cell]});
            }
        }
    }
    if (!entries.empty()) {
        ZoomWritten zoom;
        zoom.sum = wholeSum;
        zoom.layout = wholeLayout;
        zoom.index = writeBlocks(out, pool, version, options.compressionLevel, entries, nullptr);
        matrices.insert(matrices.begin(), {"0_0", writeMatrixHeader(out, 0, 0, {zoom})});
    }
    std::vector<Entry>().swap(entries);

    // Preprocessor.writeFooter
    const int64_t masterPosition = out.position();
    ByteWriter footer;
    if (version > 8) {
        footer.put<int64_t>(0);
    } else {
        footer.put<int32_t>(0);
    }
    const size_t counted = footer.size();
    footer.put<int32_t>(static_cast<int32_t>(matrices.size()));
    for (const auto& [key, entry] : matrices) {
        footer.cstr(key);
        footer.put<int64_t>(entry.position);
        footer.put<int32_t>(entry.size);
    }
    footer.put<int32_t>(static_cast<int32_t>(expectedNone.size()));
    for (auto& ev : expectedNone) {
        putExpected(footer, version, ev, false);
    }
    const int64_t nBytes = static_cast<int64_t>(footer.size() - counted);
    if (version > 8) {
        std::memcpy(footer.bytes.data(), &nBytes, sizeof(int64_t));
    } else {
        const auto nBytes32 = static_cast<int32_t>(nBytes);
        std::memcpy(footer.bytes.data(), &nBytes32, sizeof(int32_t));
    }
    out.write(footer);
    const auto [indexPosition, indexSize] = norms.write(out);
    if (version > 8) {
        out.patchValue<int64_t>(nviOffset, indexPosition);
        out.patchValue<int64_t>(nviOffset + 8, indexSize);
    }
    out.patchValue<int64_t>(masterOffset, masterPosition);
    out.close();
}

void addNorm(const std::string& path, const std::vector<std::string>& normalizations, int threads) {
    using namespace detail;
    validateNormalizations(normalizations);
    if (threads < 1) {
        throw HicError("threads must be at least 1");
    }
    int64_t truncateAt = 0;
    int32_t version = 0;
    int64_t nviOffset = 0;
    std::optional<NormBuilder> builder;
    {
        const HiCFile file(path);
        const FileState& state = stateOf(file);
        if (!state.fragResolutions.empty()) {
            throw HicError("addNorm supports base pair resolutions only; " + path + " has FRAG resolutions");
        }
        version = state.version;
        truncateAt = state.normalizedSectionPosition;
        nviOffset = 4 + 4 + 8 + static_cast<int64_t>(state.genome.size()) + 1;
        std::vector<int64_t> lengths;
        for (const auto& chromosome : state.chromosomes) {
            lengths.push_back(chromosome.length);
        }
        builder.emplace(version, lengths, state.bpResolutions, normalizations);
        for (const int32_t resolution : state.bpResolutions) {
            for (size_t chr = 1; chr < state.chromosomes.size(); ++chr) {
                const auto index = static_cast<int32_t>(chr);
                const auto headers = file.matrixZoomHeaders(index, index);
                const auto header = std::find_if(headers.begin(), headers.end(), [&](const ZoomHeader& h) {
                    return h.unit == "BP" && h.binSize == resolution;
                });
                if (header == headers.end()) {
                    continue;
                }
                const auto mzd = file.getMatrixZoomData(state.chromosomes[chr].name, state.chromosomes[chr].name,
                                                        "observed", "NONE", "BP", resolution);
                std::vector<NormRecord> records;
                mzd.forEachBlock(
                    [&](const BlockIndexEntry&, std::vector<ContactRecord>& block) {
                        for (const auto& r : block) {
                            records.push_back(NormRecord{r.binX, r.binY, r.counts});
                        }
                    },
                    threads);
                builder->addChromosome(index, resolution, records,
                                       static_cast<int64_t>(header->blockBinCount) * header->blockColumnCount);
            }
        }
    }
    OutputFile out(path, truncateAt);
    const auto [indexPosition, indexSize] = builder->write(out);
    if (version > 8) {
        out.patchValue<int64_t>(nviOffset, indexPosition);
        out.patchValue<int64_t>(nviOffset + 8, indexSize);
    }
    out.close();
}

}  // namespace hicfilecpp
