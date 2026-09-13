#include "hicx/hic_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <utility>

#include <hicfilecpp/hicfilecpp.hpp>

#include "hicx/hdf5_util.hpp"

namespace hicx {

namespace {

// hic2cool/_version.py of the reference environment.
constexpr const char* kHic2coolVersion = "0.8.3";
// hic2cool_config.py URL, COOLER_FORMAT, COOLER_FORMAT_VERSION, MCOOL_FORMAT
// and MCOOL_FORMAT_VERSION.
constexpr const char* kHic2coolUrl = "https://github.com/4dn-dcic/hic2cool";
constexpr std::size_t kChromNameWidth = 32;  // CHROM_DTYPE S32

std::string lowercase(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// datetime.utcnow().isoformat()
std::string utc_now_isoformat() {
    const auto now = std::chrono::system_clock::now();
    const auto micros =
        std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count() %
        1000000;
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
    gmtime_r(&seconds, &utc);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &utc);
    std::string text = buffer;
    if (micros != 0) {
        std::snprintf(buffer, sizeof(buffer), ".%06lld", static_cast<long long>(micros));
        text += buffer;
    }
    return text;
}

// The repr of a Python list of ints.
std::string python_list(const std::vector<std::int32_t>& values) {
    std::string text = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        text += (i ? ", " : "") + std::to_string(values[i]);
    }
    return text + "]";
}

// The outfile renaming of hic2cool_convert.
std::string hic2cool_output_name(const std::string& outfile, bool multi_res) {
    if (ends_with(outfile, ".multi.cool")) {
        return multi_res ? outfile : outfile.substr(0, outfile.size() - 11) + ".cool";
    }
    if (ends_with(outfile, ".mcool")) {
        return multi_res ? outfile : outfile.substr(0, outfile.size() - 6) + ".cool";
    }
    if (ends_with(outfile, ".cool")) {
        return multi_res ? outfile.substr(0, outfile.size() - 5) + ".mcool" : outfile;
    }
    return outfile + (multi_res ? ".mcool" : ".cool");
}

// Assigning a Python float to an element of an int32 numpy array: truncation,
// with the platform's INT32_MIN for values it cannot represent.
std::int32_t numpy_int32(float value) {
    if (!std::isfinite(value) || value >= 2147483648.0f || value < -2147483648.0f) {
        return std::numeric_limits<std::int32_t>::min();
    }
    return static_cast<std::int32_t>(value);
}

struct Chrom {
    std::int32_t index = 0;
    std::string name;
    std::int64_t length = 0;
};

h5::Handle write_column(h5::FileWriter& writer, const std::string& path, hid_t file_type,
                        hid_t memory_type, std::size_t length, const void* data) {
    h5::Handle dataset =
        writer.create_dataset(path, file_type, length, length, h5::Filter::Hic2cool);
    if (length > 0) {
        h5::FileWriter::write_block(dataset.get(), memory_type, 0, length, data);
    }
    return dataset;
}

// hic2cool_convert after the header and footer are read: initialize_res,
// the pixel loop of parse_hic and write_pixels_chunk, and
// finalize_resolution_cool, for every resolution in use.
void convert(const hicfilecpp::HiCFile& hic, const std::string& outfile,
             const std::vector<std::int32_t>& resolutions, bool multi_res) {
    // read_header keeps a chromosome only when its name and length are truthy.
    std::vector<Chrom> used;
    for (const auto& chromosome : hic.getChromosomes()) {
        if (!chromosome.name.empty() && chromosome.length != 0) {
            used.push_back(Chrom{chromosome.index, chromosome.name, chromosome.length});
        }
    }
    // read_footer: NORMS in the order the normalized expected values name them.
    std::vector<std::string> norms;
    for (const auto& key : hic.expectedValuesKeys()) {
        if (key.normalization != "NONE" &&
            std::find(norms.begin(), norms.end(), key.normalization) == norms.end()) {
            norms.push_back(key.normalization);
        }
    }
    std::vector<const Chrom*> chroms;
    for (const auto& chrom : used) {
        if (lowercase(chrom.name) != "all") {
            chroms.push_back(&chrom);
        }
    }
    std::map<std::int32_t, std::string> name_of;
    for (const auto& chrom : used) {
        name_of[chrom.index] = chrom.name;
    }
    const std::string genome = hic.getGenomeID();
    std::map<std::string, std::string> metadata_seen;
    std::vector<std::pair<std::string, std::string>> metadata;
    for (const auto& [key, value] : hic.attributes()) {
        if (metadata_seen.count(key) == 0) {
            metadata.emplace_back(key, value);
        } else {
            for (auto& entry : metadata) {
                if (entry.first == key) {
                    entry.second = value;
                }
            }
        }
        metadata_seen[key] = value;
    }

    for (std::size_t r = 0; r < resolutions.size(); ++r) {
        const std::int32_t binsize = resolutions[r];
        h5::FileWriter writer(outfile, h5::WriteMode::Append);
        std::string group = "/";
        if (multi_res) {
            if (r == 0) {
                writer.set_attribute("/", "format", std::string("HDF5::MCOOL"));
                writer.set_attribute("/", "format-version", static_cast<std::int64_t>(2));
                writer.create_group("/resolutions");
            }
            group = "/resolutions/" + std::to_string(binsize);
            writer.create_group(group);
        }
        const std::string prefix = group == "/" ? "/" : group + "/";

        // write_chroms
        writer.create_group(prefix + "chroms");
        std::vector<char> names(chroms.size() * kChromNameWidth, '\0');
        std::vector<std::int32_t> lengths;
        std::vector<std::string> enum_names;
        for (std::size_t i = 0; i < chroms.size(); ++i) {
            const std::string& name = chroms[i]->name;
            std::copy_n(name.begin(), std::min(name.size(), kChromNameWidth),
                        names.begin() + static_cast<std::ptrdiff_t>(i * kChromNameWidth));
            enum_names.push_back(name.substr(0, std::min(name.size(), kChromNameWidth)));
            lengths.push_back(static_cast<std::int32_t>(chroms[i]->length));
        }
        {
            const h5::Handle string_type = h5::fixed_string_type(kChromNameWidth);
            write_column(writer, prefix + "chroms/name", string_type.get(), string_type.get(),
                         chroms.size(), names.data());
        }
        write_column(writer, prefix + "chroms/length", H5T_STD_I32LE, H5T_NATIVE_INT32,
                     lengths.size(), lengths.data());

        // create_bins
        std::vector<std::int32_t> chrom_ids;
        std::vector<std::int32_t> starts;
        std::vector<std::int32_t> ends;
        std::vector<std::int64_t> chrom_offsets{0};
        std::map<std::int32_t, std::int64_t> bins_of;
        std::map<std::int32_t, std::int64_t> offset_of;
        std::int32_t position = 0;
        for (const auto& chrom : used) {
            if (lowercase(chrom.name) == "all") {
                continue;
            }
            std::int64_t start = 0;
            std::int64_t count = 0;
            while (start < chrom.length) {
                const std::int64_t end = std::min(start + binsize, chrom.length);
                chrom_ids.push_back(position);
                starts.push_back(static_cast<std::int32_t>(start));
                ends.push_back(static_cast<std::int32_t>(end));
                count++;
                start = end;
            }
            offset_of[chrom.index] = chrom_offsets.back();
            bins_of[chrom.index] = count;
            if (chrom_offsets.size() < used.size()) {
                chrom_offsets.push_back(
                    static_cast<std::int64_t>(std::ceil(static_cast<double>(chrom.length) / binsize)) +
                    chrom_offsets.back());
            }
            position++;
        }
        const std::size_t n_bins = chrom_ids.size();

        // write_bins
        writer.create_group(prefix + "bins");
        {
            const h5::Handle file_enum = h5::enum_type(enum_names);
            const h5::Handle memory_enum = h5::enum_type(enum_names, H5T_NATIVE_INT32);
            write_column(writer, prefix + "bins/chrom", file_enum.get(), memory_enum.get(), n_bins,
                         chrom_ids.data());
        }
        write_column(writer, prefix + "bins/start", H5T_STD_I32LE, H5T_NATIVE_INT32, n_bins,
                     starts.data());
        write_column(writer, prefix + "bins/end", H5T_STD_I32LE, H5T_NATIVE_INT32, n_bins,
                     ends.data());
        for (const auto& norm : norms) {
            std::vector<double> column;
            column.reserve(n_bins);
            for (const Chrom* chrom : chroms) {
                const auto chr_bins = static_cast<std::size_t>(bins_of[chrom->index]);
                const auto vector = hic.readNormVector(norm, chrom->index, "BP", binsize);
                if (!vector) {
                    column.insert(column.end(), chr_bins, std::numeric_limits<double>::quiet_NaN());
                } else {
                    column.insert(column.end(), vector->begin(),
                                  vector->begin() + static_cast<std::ptrdiff_t>(
                                                        std::min(vector->size(), chr_bins)));
                }
            }
            if (column.size() != n_bins) {
                throw Hic2coolExit("!!! ERROR. Length of normalization vector " + norm +
                                   " does not match the number of bins.\nThis is likely a "
                                   "problem with the hic file");
            }
            write_column(writer, prefix + "bins/" + norm, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, n_bins,
                         column.data());
        }

        // write_chrom_offset
        writer.create_group(prefix + "indexes");
        write_column(writer, prefix + "indexes/chrom_offset", H5T_STD_I64LE, H5T_NATIVE_INT64,
                     chrom_offsets.size(), chrom_offsets.data());

        // initialize_pixels
        writer.create_group(prefix + "pixels");
        h5::Handle bin1_dataset =
            writer.create_dataset(prefix + "pixels/bin1_id", H5T_STD_I64LE, 0,
                                  h5::FileWriter::kUnlimited, h5::Filter::Hic2cool);
        h5::Handle bin2_dataset =
            writer.create_dataset(prefix + "pixels/bin2_id", H5T_STD_I64LE, 0,
                                  h5::FileWriter::kUnlimited, h5::Filter::Hic2cool);
        h5::Handle count_dataset =
            writer.create_dataset(prefix + "pixels/count", H5T_STD_I32LE, 0,
                                  h5::FileWriter::kUnlimited, h5::Filter::Hic2cool);

        // The resolution's attributes: the .hic header's attributes, then
        // hic2cool's own.
        for (const auto& [key, value] : metadata) {
            writer.set_attribute(group, key, value);
        }
        writer.set_attribute(group, "nchroms", static_cast<std::int64_t>(chroms.size()));
        writer.set_attribute(group, "nbins", static_cast<std::int64_t>(n_bins));
        writer.set_attribute(group, "bin-type", std::string("fixed"));
        writer.set_attribute(group, "bin-size", static_cast<std::int64_t>(binsize));
        writer.set_attribute(group, "format", std::string("HDF5::Cooler"));
        writer.set_attribute(group, "format-url", std::string(kHic2coolUrl));
        writer.set_attribute(group, "format-version", static_cast<std::int64_t>(3));
        writer.set_attribute(group, "storage-mode", std::string("symmetric-upper"));
        writer.set_attribute(group, "generated-by", std::string("hic2cool-") + kHic2coolVersion);
        writer.set_attribute(group, "genome-assembly", genome);
        writer.set_attribute(group, "creation-date", utc_now_isoformat());

        // The pixel loop: every chromosome pair once, grouped by the first
        // chromosome in file order, each group sorted by (bin1, bin2).
        struct Pixel {
            std::int64_t bin1;
            std::int64_t bin2;
            std::int32_t count;
        };
        std::set<std::pair<std::int32_t, std::int32_t>> covered;
        std::vector<std::int64_t> bin1_counts(n_bins, 0);
        std::size_t nnz = 0;
        std::vector<std::int64_t> column64;
        std::vector<std::int32_t> column32;
        for (const auto& chr_a : used) {
            if (lowercase(chr_a.name) == "all") {
                continue;
            }
            std::vector<Pixel> total;
            for (const auto& chr_b : used) {
                if (lowercase(chr_b.name) == "all") {
                    continue;
                }
                const std::int32_t c1 = std::min(chr_a.index, chr_b.index);
                const std::int32_t c2 = std::max(chr_a.index, chr_b.index);
                if (!covered.insert({c1, c2}).second) {
                    continue;
                }
                if (!hic.hasMatrix(c1, c2)) {
                    continue;
                }
                const auto headers = hic.matrixZoomHeaders(c1, c2);
                const bool present =
                    std::any_of(headers.begin(), headers.end(), [&](const hicfilecpp::ZoomHeader& h) {
                        return h.unit == "BP" && h.binSize == binsize;
                    });
                if (!present) {
                    continue;
                }
                const auto mzd = hic.getMatrixZoomData(name_of[c1], name_of[c2], "observed", "NONE",
                                                       "BP", binsize);
                const std::int64_t bins1 = bins_of[c1];
                const std::int64_t bins2 = bins_of[c2];
                const std::int64_t offset1 = offset_of[c1];
                const std::int64_t offset2 = offset_of[c2];
                for (const auto& entry : mzd.blockIndex()) {
                    for (const auto& record : mzd.readBlock(entry)) {
                        if (record.binX >= 0 && record.binX < bins1 && record.binY >= 0 &&
                            record.binY < bins2) {
                            total.push_back(Pixel{record.binX + offset1, record.binY + offset2,
                                                  numpy_int32(record.counts)});
                        }
                    }
                }
            }
            std::sort(total.begin(), total.end(), [](const Pixel& a, const Pixel& b) {
                return a.bin1 != b.bin1 ? a.bin1 < b.bin1 : a.bin2 < b.bin2;
            });
            if (total.empty()) {
                continue;
            }
            const std::size_t grown = nnz + total.size();
            h5::FileWriter::resize(bin1_dataset.get(), grown);
            h5::FileWriter::resize(bin2_dataset.get(), grown);
            h5::FileWriter::resize(count_dataset.get(), grown);
            column64.resize(total.size());
            for (std::size_t i = 0; i < total.size(); ++i) {
                column64[i] = total[i].bin1;
                bin1_counts[static_cast<std::size_t>(total[i].bin1)]++;
            }
            h5::FileWriter::write_block(bin1_dataset.get(), H5T_NATIVE_INT64, nnz, total.size(),
                                        column64.data());
            for (std::size_t i = 0; i < total.size(); ++i) {
                column64[i] = total[i].bin2;
            }
            h5::FileWriter::write_block(bin2_dataset.get(), H5T_NATIVE_INT64, nnz, total.size(),
                                        column64.data());
            column32.resize(total.size());
            for (std::size_t i = 0; i < total.size(); ++i) {
                column32[i] = total[i].count;
            }
            h5::FileWriter::write_block(count_dataset.get(), H5T_NATIVE_INT32, nnz, total.size(),
                                        column32.data());
            nnz = grown;
        }

        // finalize_resolution_cool
        std::vector<std::int64_t> bin1_offset(n_bins + 1, 0);
        for (std::size_t i = 0; i < n_bins; ++i) {
            bin1_offset[i + 1] = bin1_offset[i] + bin1_counts[i];
        }
        write_column(writer, prefix + "indexes/bin1_offset", H5T_STD_I64LE, H5T_NATIVE_INT64,
                     bin1_offset.size(), bin1_offset.data());
        writer.set_attribute(group, "nnz", static_cast<std::int64_t>(nnz));
    }
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
    const hicfilecpp::HiCFile hic(hic_path);
    const std::vector<std::int32_t> resolutions = hic.getResolutions();
    if (resolution != 0 &&
        std::find(resolutions.begin(), resolutions.end(), resolution) == resolutions.end()) {
        throw Hic2coolExit("!!! ERROR. Given binsize (in bp) is not a supported resolution in "
                           "this file.\nPlease use 0 (all resolutions) or use one of: " +
                           python_list(resolutions));
    }
    const std::vector<std::int32_t> use =
        resolution == 0 ? resolutions
                        : std::vector<std::int32_t>{static_cast<std::int32_t>(resolution)};
    const bool multi_res = use.size() > 1;
    const std::string written = hic2cool_output_name(outfile, multi_res);
    std::filesystem::remove(written);
    convert(hic, written, use, multi_res);
    return written;
}

void hic2cool_convert_mcool(const std::string& hic_path, const std::string& outfile,
                            const std::vector<std::int64_t>& resolutions) {
    const hicfilecpp::HiCFile hic(hic_path);
    const std::vector<std::int32_t> available = hic.getResolutions();
    std::vector<std::int32_t> use;
    if (resolutions.empty()) {
        use = available;
    } else {
        for (const std::int64_t resolution : resolutions) {
            if (std::find(available.begin(), available.end(), resolution) == available.end()) {
                throw Hic2coolExit(
                    "!!! ERROR. Given binsize (in bp) is not a supported resolution in "
                    "this file.\nPlease use 0 (all resolutions) or use one of: " +
                    python_list(available));
            }
            use.push_back(static_cast<std::int32_t>(resolution));
        }
    }
    std::filesystem::remove(outfile);
    convert(hic, outfile, use, true);
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

}  // namespace hicx
