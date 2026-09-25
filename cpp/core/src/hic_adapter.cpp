#include "hicx/hic_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <thread>
#include <utility>

#include <coolercpp/region.hpp>
#include <hic2cool/hic2cool.hpp>
#include <hicfilecpp/hicfilecpp.hpp>

#include "hicx/cool_adapter.hpp"
#include "hicx/hdf5_util.hpp"

namespace hicx {

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// hic2cool prints its progress; hicConvertFormat does not.
hic2cool::Console quiet_console() {
    hic2cool::Console console;
    console.out = [](const std::string&) {};
    console.err = [](const std::string& line) { std::fprintf(stderr, "%s\n", line.c_str()); };
    console.input = [](const std::string&) -> std::string {
        throw hic2cool::ExitError("EOFError: EOF when reading a line");
    };
    return console;
}

}  // namespace

std::vector<std::int64_t> hic_resolutions(const std::string& hic_path) {
    const hicfilecpp::HiCFile hic(hic_path);
    std::vector<std::int64_t> resolutions;
    for (const std::int32_t resolution : hic.getResolutions()) {
        resolutions.push_back(resolution);
    }
    return resolutions;
}

std::string hic2cool_convert(const std::string& hic_path, const std::string& outfile,
                             std::int64_t resolution) {
    try {
        return hic2cool::hic2cool_convert(hic_path, outfile, resolution, 0, false, true, "symmetric-upper",
                                          quiet_console());
    } catch (const hic2cool::ExitError& error) {
        throw Hic2coolExit(error.what());
    }
}

void hic2cool_convert_mcool(const std::string& hic_path, const std::string& outfile,
                            const std::vector<std::int64_t>& resolutions) {
    try {
        hic2cool::hic2cool_convert_mcool(hic_path, outfile, resolutions, 0, false, true, quiet_console());
    } catch (const hic2cool::ExitError& error) {
        throw Hic2coolExit(error.what());
    }
}

// ------------------------------------------------------------- write_hic

namespace {

struct ChromBins {
    std::string name;
    std::int64_t length = 0;
    std::int64_t first = 0;
    std::int64_t end = 0;
};

struct MatrixLayout {
    std::int64_t bin_size = 0;
    std::vector<ChromBins> chromosomes;
};

MatrixLayout layout_of(const MatrixData& matrix) {
    MatrixLayout layout;
    const auto& intervals = matrix.cut_intervals;
    for (const auto& interval : intervals) {
        layout.bin_size = std::max(layout.bin_size, interval.end - interval.start);
    }
    if (intervals.empty() || layout.bin_size <= 0) {
        throw std::runtime_error("the matrix has no bins to write as .hic");
    }
    for (std::size_t i = 0; i < intervals.size(); ++i) {
        if (layout.chromosomes.empty() || layout.chromosomes.back().name != intervals[i].chrom) {
            for (const auto& seen : layout.chromosomes) {
                if (seen.name == intervals[i].chrom) {
                    throw std::runtime_error("chromosome " + seen.name +
                                             " is not contiguous in the bin table");
                }
            }
            layout.chromosomes.push_back(ChromBins{intervals[i].chrom, 0,
                                                   static_cast<std::int64_t>(i),
                                                   static_cast<std::int64_t>(i)});
        }
        ChromBins& chrom = layout.chromosomes.back();
        const std::int64_t local = static_cast<std::int64_t>(i) - chrom.first;
        if (intervals[i].start != local * layout.bin_size) {
            throw std::runtime_error(
                "the .hic format needs fixed size bins starting at 0; bin " + std::to_string(i) +
                " of " + chrom.name + " starts at " + std::to_string(intervals[i].start));
        }
        chrom.end = static_cast<std::int64_t>(i) + 1;
        chrom.length = intervals[i].end;
    }
    return layout;
}

class MatrixSource : public hicfilecpp::PixelSource {
  public:
    MatrixSource(std::map<std::int32_t, std::pair<const MatrixData*, MatrixLayout>> matrices,
                 std::int32_t base)
        : matrices_(std::move(matrices)), base_(base) {}

    void pixels(std::int32_t resolution, std::int32_t chr1, std::int32_t chr2,
                const std::function<void(const hicfilecpp::Pixel*, std::size_t)>& consume) override {
        auto it = matrices_.find(resolution);
        if (it == matrices_.end()) {
            it = matrices_.find(base_);
        }
        const MatrixData& matrix = *it->second.first;
        const MatrixLayout& layout = it->second.second;
        const ChromBins& a = layout.chromosomes[static_cast<std::size_t>(chr1)];
        const ChromBins& b = layout.chromosomes[static_cast<std::size_t>(chr2)];
        const auto& indptr = matrix.matrix.indptr();
        const auto& indices = matrix.matrix.indices();
        const auto& data = matrix.matrix.data();
        std::vector<hicfilecpp::Pixel> batch;
        batch.reserve(1 << 16);
        for (std::int64_t row = a.first; row < a.end; ++row) {
            for (std::int64_t k = indptr[static_cast<std::size_t>(row)];
                 k < indptr[static_cast<std::size_t>(row) + 1]; ++k) {
                const std::int64_t col = indices[static_cast<std::size_t>(k)];
                if (col < b.first || col >= b.end || (chr1 == chr2 && col < row)) {
                    continue;
                }
                const double value = data[static_cast<std::size_t>(k)];
                if (!std::isfinite(value) || value == 0) {
                    continue;
                }
                batch.push_back(hicfilecpp::Pixel{static_cast<std::int32_t>(row - a.first),
                                                  static_cast<std::int32_t>(col - b.first),
                                                  static_cast<float>(value)});
                if (batch.size() == batch.capacity()) {
                    consume(batch.data(), batch.size());
                    batch.clear();
                }
            }
        }
        if (!batch.empty()) {
            consume(batch.data(), batch.size());
        }
    }

  private:
    std::map<std::int32_t, std::pair<const MatrixData*, MatrixLayout>> matrices_;
    std::int32_t base_;
};

}  // namespace

void write_hic(const std::string& path, const std::vector<const MatrixData*>& matrices,
               const std::vector<std::int64_t>& extra_resolutions,
               const HicWriteOptions& options) {
    if (matrices.empty()) {
        throw std::runtime_error("no matrix to write as .hic");
    }
    if (matrices.size() > 1 && !extra_resolutions.empty()) {
        throw std::runtime_error("extra resolutions need a single input matrix");
    }
    // A chromosome's length is the end of its last bin. The resolutions of an
    // mcool file need not clip that bin at the same position, so the header
    // takes the largest end, which every resolution's bins fit in.
    std::map<std::int32_t, std::pair<const MatrixData*, MatrixLayout>> by_resolution;
    std::vector<std::pair<std::string, std::int64_t>> chromosomes;
    for (const MatrixData* matrix : matrices) {
        MatrixLayout layout = layout_of(*matrix);
        if (chromosomes.empty()) {
            for (const auto& chrom : layout.chromosomes) {
                chromosomes.emplace_back(chrom.name, chrom.length);
            }
        } else {
            if (layout.chromosomes.size() != chromosomes.size()) {
                throw std::runtime_error("the resolutions do not share the same chromosomes");
            }
            for (std::size_t c = 0; c < chromosomes.size(); ++c) {
                if (layout.chromosomes[c].name != chromosomes[c].first) {
                    throw std::runtime_error("the resolutions do not share the same chromosomes");
                }
                chromosomes[c].second = std::max(chromosomes[c].second, layout.chromosomes[c].length);
            }
        }
        const auto resolution = static_cast<std::int32_t>(layout.bin_size);
        if (!by_resolution.emplace(resolution, std::make_pair(matrix, std::move(layout))).second) {
            throw std::runtime_error("resolution " + std::to_string(resolution) + " is given twice");
        }
    }
    hicfilecpp::WriteOptions write;
    write.version = options.version;
    write.genomeId = options.genome;
    write.chromosomes = chromosomes;
    write.normalizations = options.normalizations;
    write.threads = options.threads;
    const std::int32_t base = by_resolution.begin()->first;
    for (const auto& entry : by_resolution) {
        write.resolutions.push_back(entry.first);
    }
    if (matrices.size() > 1) {
        write.sourceProvidesEveryResolution = true;
    } else {
        write.sourceResolution = base;
        for (const std::int64_t resolution : extra_resolutions) {
            if (resolution != base) {
                write.resolutions.push_back(static_cast<std::int32_t>(resolution));
            }
        }
    }
    MatrixSource source(std::move(by_resolution), base);
    hicfilecpp::writeHicFile(path, write, source);
}

// --------------------------------------------------------- generic loading

namespace {

std::vector<std::string> split_slash(const std::string& text) {
    std::vector<std::string> parts;
    std::string current;
    for (const char c : text) {
        if (c == '/') {
            if (!current.empty()) {
                parts.push_back(current);
                current.clear();
            }
        } else {
            current += c;
        }
    }
    if (!current.empty()) {
        parts.push_back(current);
    }
    return parts;
}

// Removes a file when it goes out of scope, the same pattern
// hicConvertFormat.cpp uses for its own hic2cool_convert temporaries.
struct TemporaryFile {
    std::string path;
    ~TemporaryFile() {
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
    }
};

std::string unique_temp_cool_path() {
    const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    std::ostringstream name;
    name << "hicx_hic_load_" << std::this_thread::get_id() << "_" << now << ".cool";
    return (std::filesystem::temp_directory_path() / name.str()).string();
}

}  // namespace

bool has_hic_signature(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }
    char magic[3] = {0, 0, 0};
    file.read(magic, 3);
    return file.gcount() == 3 && magic[0] == 'H' && magic[1] == 'I' && magic[2] == 'C';
}

bool is_hic_path(const std::string& path) {
    const std::string base = path.substr(0, path.find("::"));
    if (has_hic_signature(base)) {
        return true;
    }
    return ends_with(base, ".hic");
}

HicUri parse_hic_uri(const std::string& uri) {
    HicUri parsed;
    const std::size_t separator = uri.find("::");
    parsed.path = uri.substr(0, separator);
    if (separator == std::string::npos) {
        return parsed;
    }
    const std::vector<std::string> tokens = split_slash(uri.substr(separator + 2));
    const auto syntax_error = [&]() -> HicUri {
        throw std::runtime_error(
            "invalid .hic selector '" + uri +
            "'; expected 'file.hic', 'file.hic::/resolutions/<N>', "
            "'file.hic::/normalizations/<NAME>' or "
            "'file.hic::/resolutions/<N>/normalizations/<NAME>'");
    };
    if (tokens.size() == 2 && tokens[0] == "resolutions") {
        parsed.resolution = std::stoll(tokens[1]);
    } else if (tokens.size() == 2 && tokens[0] == "normalizations") {
        parsed.normalization = tokens[1];
    } else if (tokens.size() == 4 && tokens[0] == "resolutions" && tokens[2] == "normalizations") {
        parsed.resolution = std::stoll(tokens[1]);
        parsed.normalization = tokens[3];
    } else {
        return syntax_error();
    }
    return parsed;
}

namespace {

std::int64_t resolve_resolution(const hicfilecpp::HiCFile& hic, const std::string& path,
                                const std::optional<std::int64_t>& requested) {
    const std::vector<std::int32_t> available = hic.getResolutions();
    if (available.empty()) {
        throw std::runtime_error(".hic file " + path + " has no resolutions");
    }
    if (requested.has_value()) {
        if (std::find(available.begin(), available.end(), *requested) == available.end()) {
            throw std::runtime_error("resolution " + std::to_string(*requested) +
                                     " is not one of " + path + "'s: " + python_list(available));
        }
        return *requested;
    }
    // hic2cool_convert's own resolution 0 ("every resolution") is not useful
    // for a single matrix load; the finest one is what a caller giving no
    // selector expects, matching cpp/PLAN.md's stated default.
    return *std::min_element(available.begin(), available.end());
}

std::string resolve_normalization(const hicfilecpp::HiCFile& hic, const std::string& path,
                                  const std::optional<std::string>& requested) {
    const std::string normalization = requested.value_or("NONE");
    if (normalization == "NONE") {
        return normalization;
    }
    const std::vector<std::string> types = hic.getNormalizationTypes();
    if (std::find(types.begin(), types.end(), normalization) == types.end()) {
        std::string list = "[";
        for (std::size_t i = 0; i < types.size(); ++i) {
            list += (i ? ", " : "") + types[i];
        }
        list += "]";
        throw std::runtime_error("normalization " + normalization + " is not one of " + path +
                                 "'s: " + list);
    }
    return normalization;
}

// The whole file, through hic2cool_convert (already validated) into a
// temporary cool file and hicx::read_cool. Used when no chrom_name narrows
// the load.
HicLoadResult read_hic_whole(const std::string& path, std::int64_t resolution,
                             const std::string& normalization) {
    TemporaryFile temporary{unique_temp_cool_path()};
    hic2cool_convert(path, temporary.path, resolution);

    CoolLoadOptions options;
    options.apply_correction = normalization != "NONE";
    options.correction_factor_table = normalization;
    if (normalization != "NONE") {
        // Every hic2cool normalization table (KR, VC, VC_SQRT, SCALE, and the
        // GW_/INTER_ variants hic2cool also writes) is a Juicer normalization
        // vector: divisive, observed / normI / normJ. cool's own auto
        // detection (cool_adapter.cpp) only recognises KR, VC and SQRT_VC by
        // column name, so it is set explicitly here instead of relied on.
        options.correction_operator = '/';
    }
    CoolLoadResult loaded = read_cool(temporary.path, options);

    HicLoadResult result;
    result.data = std::move(loaded.data);
    result.correction_operator = loaded.correction_operator;
    result.metadata = std::move(loaded.metadata);
    return result;
}

// cool.py's own NaN bin heuristic (cool_adapter.cpp read_cool, "the bins that
// hold no interaction at all"), reproduced on a matrix built directly from
// .hic records instead of a cooler's pixel table.
std::vector<std::int64_t> nan_bins_of(const CsrMatrix& matrix) {
    std::vector<std::int64_t> nan_bins;
    const std::int64_t shape = std::min(matrix.rows(), matrix.cols());
    std::vector<char> used_as_column(static_cast<std::size_t>(shape), 0);
    for (const std::int32_t column : matrix.indices()) {
        if (column >= 0 && column < shape) {
            used_as_column[static_cast<std::size_t>(column)] = 1;
        }
    }
    for (std::int64_t bin = 0; bin < shape; ++bin) {
        if (used_as_column[static_cast<std::size_t>(bin)] != 0) {
            continue;
        }
        const bool empty_row = matrix.indptr()[static_cast<std::size_t>(bin)] ==
                               matrix.indptr()[static_cast<std::size_t>(bin) + 1];
        if (empty_row) {
            nan_bins.push_back(bin);
        }
    }
    return nan_bins;
}

// One chromosome, or one genomic sub-range of it, through hicfilecpp's block
// indexed MatrixZoomData::getRecords, which decodes only the blocks the
// query's genomic range selects (hicfilecpp/src/reader.cpp: z.blockNumbers(region)
// picks the block numbers before anything is decoded). Neither the whole
// file nor the whole chromosome is read when the region is narrower than it;
// this is the .hic side of cool's existing CoolFile::extent + read_block fast
// path (cool_adapter.cpp read_cool, options.chrom_name).
HicLoadResult read_hic_region(const hicfilecpp::HiCFile& hic, const std::string& path,
                              std::int64_t resolution, const std::string& normalization,
                              const std::string& chrom_name) {
    const coolercpp::GenomicRange range = coolercpp::parse_region_string(chrom_name);
    const std::vector<hicfilecpp::Chromosome> chromosomes = hic.getChromosomes();
    const auto found = std::find_if(chromosomes.begin(), chromosomes.end(),
                                    [&](const hicfilecpp::Chromosome& c) { return c.name == range.chrom; });
    if (found == chromosomes.end()) {
        throw std::runtime_error("chromosome " + range.chrom + " is not one of " + path + "'s");
    }
    const std::int64_t length = found->length;
    const std::int64_t start = std::max<std::int64_t>(0, range.start.value_or(0));
    const std::int64_t end = std::min(length, range.end.value_or(length));
    if (start >= end) {
        throw std::runtime_error("region " + chrom_name + " of " + path + " is empty");
    }

    const std::int64_t nbins_chrom = (length + resolution - 1) / resolution;
    const std::int64_t bin_first = start / resolution;
    const std::int64_t bin_last = std::min(nbins_chrom, (end + resolution - 1) / resolution);
    const std::int64_t nbins_local = bin_last - bin_first;
    const std::int64_t gx0 = bin_first * resolution;
    const std::int64_t gx1 = std::min(bin_last * resolution, length);

    hicfilecpp::MatrixZoomData mzd =
        hic.getMatrixZoomData(range.chrom, range.chrom, "observed", normalization, "BP", resolution);
    if (!mzd.found()) {
        throw std::runtime_error("no matrix for " + range.chrom + " x " + range.chrom + " in " + path +
                                 ": " + mzd.message());
    }
    const std::vector<hicfilecpp::ContactRecord> records = mzd.getRecords(gx0, gx1, gx0, gx1);

    // COO triplets, canonicalised to the upper triangle (row <= col) the way
    // convert()'s MatrixSource assumes on write, so symmetrize_in_place (every
    // ToolMatrix::load caller runs it) sees exactly what a cool block gives it.
    std::vector<std::int32_t> row;
    std::vector<std::int32_t> col;
    std::vector<double> data;
    row.reserve(records.size());
    col.reserve(records.size());
    data.reserve(records.size());
    for (const hicfilecpp::ContactRecord& record : records) {
        const std::int64_t x = (static_cast<std::int64_t>(record.binX) - gx0) / resolution;
        const std::int64_t y = (static_cast<std::int64_t>(record.binY) - gx0) / resolution;
        const std::int64_t r = std::min(x, y);
        const std::int64_t c = std::max(x, y);
        if (r < 0 || c >= nbins_local) {
            continue;
        }
        row.push_back(static_cast<std::int32_t>(r));
        col.push_back(static_cast<std::int32_t>(c));
        data.push_back(static_cast<double>(record.counts));
    }

    HicLoadResult result;
    result.data.matrix =
        CsrMatrix::from_coo(nbins_local, nbins_local, row, col, std::move(data), "float64");

    result.data.cut_intervals.reserve(static_cast<std::size_t>(nbins_local));
    for (std::int64_t i = 0; i < nbins_local; ++i) {
        CutInterval interval;
        interval.chrom = range.chrom;
        interval.start = (bin_first + i) * resolution;
        interval.end = std::min((bin_first + i + 1) * resolution, length);
        result.data.cut_intervals.push_back(std::move(interval));
    }

    if (normalization != "NONE") {
        result.correction_operator = '/';
        const std::optional<std::vector<double>> vector =
            hic.readNormVector(normalization, found->index, "BP", static_cast<std::int32_t>(resolution));
        std::vector<double> factors(static_cast<std::size_t>(nbins_local),
                                    std::numeric_limits<double>::quiet_NaN());
        if (vector.has_value()) {
            for (std::int64_t i = 0; i < nbins_local; ++i) {
                const auto source = static_cast<std::size_t>(bin_first + i);
                if (source < vector->size()) {
                    factors[static_cast<std::size_t>(i)] = (*vector)[source];
                }
            }
        }
        result.data.correction_factors = std::move(factors);
    }

    result.data.nan_bins = nan_bins_of(result.data.matrix);
    return result;
}

}  // namespace

HicLoadResult read_hic(const std::string& uri, const std::optional<std::string>& chrom_name) {
    const HicUri parsed = parse_hic_uri(uri);
    const hicfilecpp::HiCFile hic(parsed.path);
    const std::int64_t resolution = resolve_resolution(hic, parsed.path, parsed.resolution);
    const std::string normalization = resolve_normalization(hic, parsed.path, parsed.normalization);

    if (chrom_name.has_value()) {
        return read_hic_region(hic, parsed.path, resolution, normalization, *chrom_name);
    }
    return read_hic_whole(parsed.path, resolution, normalization);
}

}  // namespace hicx
