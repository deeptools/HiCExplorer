// hic2cool_convert (hic2cool_utils.py of hic2cool 1.0.1): initialize_res,
// write_chroms, create_bins, write_bins, parse_hic, write_pixels_chunk and
// finalize_resolution_cool, with the same datasets, attributes, messages and
// pixel order, on several threads.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <map>
#include <set>

#include <hicfilecpp/hicfilecpp.hpp>

#include "io.hpp"
#include "hic2cool/hic2cool.hpp"
#include "hic_header.hpp"
#include "pyformat.hpp"

namespace hic2cool {

namespace {

using namespace detail;

constexpr const char* kFormatUrl = "https://github.com/4dn-dcic/hic2cool";
constexpr std::size_t kChromNameWidth = 32;  // CHROM_DTYPE S32
constexpr std::size_t kFlushChunks = 1024;     // pixel chunks buffered per column before writing
constexpr std::size_t kAppendSlice = std::size_t{1} << 20;  // rows converted per append

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string output_name(const std::string& outfile, bool multi_res) {
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

// counts[k] = c: a Python float stored into an int32 numpy array.
std::int32_t numpy_int32(float value) {
    if (!std::isfinite(value) || value >= 2147483648.0f || value < -2147483648.0f) {
        return std::numeric_limits<std::int32_t>::min();
    }
    return static_cast<std::int32_t>(value);
}

// The smallest bin1 a block can hold, from its number: column * blockBinCount
// for versions 6 to 8 and inter-chromosomal matrices; for version 9
// intra-chromosomal matrices, whose blocks run along the diagonal, the
// position along the diagonal less half the widest distance the block's depth
// allows (one depth level of margin).
std::int64_t block_bin1_floor(std::int32_t number, std::int32_t block_bin_count, std::int32_t block_column_count,
                              bool v9_intra) {
    if (block_bin_count <= 0 || block_column_count <= 0 || number < 0) {
        return 0;
    }
    const std::int64_t column = number % block_column_count;
    if (!v9_intra) {
        return column * block_bin_count;
    }
    const std::int64_t depth = std::min<std::int64_t>(number / block_column_count, 40);
    const double half_span =
        std::sqrt(2.0) * block_bin_count * (std::ldexp(1.0, static_cast<int>(depth) + 2) - 1.0) / 2.0;
    const double floor = static_cast<double>(column * block_bin_count) - half_span;
    return floor <= 0.0 ? 0 : static_cast<std::int64_t>(floor);
}

struct Pixel {
    std::int32_t bin1;
    std::int32_t bin2;
    std::int32_t count;
};

// numpy's sort of the structured chunk with order=['bin1_id', 'bin2_id']:
// the remaining field breaks ties.
bool pixel_less(const Pixel& a, const Pixel& b) {
    if (a.bin1 != b.bin1) {
        return a.bin1 < b.bin1;
    }
    if (a.bin2 != b.bin2) {
        return a.bin2 < b.bin2;
    }
    return a.count < b.count;
}

void parallel_sort(ThreadPool& pool, Pixel* begin, std::size_t count) {
    const std::size_t parts = static_cast<std::size_t>(pool.threads());
    if (parts <= 1 || count < (std::size_t{1} << 16)) {
        std::sort(begin, begin + count, pixel_less);
        return;
    }
    std::vector<std::size_t> bounds(parts + 1);
    for (std::size_t p = 0; p <= parts; ++p) {
        bounds[p] = count * p / parts;
    }
    pool.run(parts, [&](std::size_t p) { std::sort(begin + bounds[p], begin + bounds[p + 1], pixel_less); });
    for (std::size_t width = 1; width < parts; width *= 2) {
        std::vector<std::size_t> merges;
        for (std::size_t p = 0; p + width < parts; p += 2 * width) {
            merges.push_back(p);
        }
        pool.run(merges.size(), [&](std::size_t m) {
            const std::size_t p = merges[m];
            std::inplace_merge(begin + bounds[p], begin + bounds[p + width],
                               begin + bounds[std::min(parts, p + 2 * width)], pixel_less);
        });
    }
}

class PixelWriter {
  public:
    PixelWriter(ThreadPool& pool, h5::File& file, const std::string& prefix, std::size_t n_bins)
        : pool_(pool),
          bin1_(file.create_dataset(prefix + "pixels/bin1_id", h5::ColumnType::int64(), 0, true), 8, true),
          bin2_(file.create_dataset(prefix + "pixels/bin2_id", h5::ColumnType::int64(), 0, true), 8, true),
          count_(file.create_dataset(prefix + "pixels/count", h5::ColumnType::int32(), 0, true), 4, true),
          bin1_counts_(n_bins, 0) {}

    // Appends pixels already in file order, a slice at a time, so that a
    // large sorted group is never held in the column buffers as a whole.
    void append(const Pixel* pixels, std::size_t n) {
        for (std::size_t lo = 0; lo < n; lo += kAppendSlice) {
            const std::size_t count = std::min(kAppendSlice, n - lo);
            b1_.resize(count);
            b2_.resize(count);
            c_.resize(count);
            for (std::size_t i = 0; i < count; ++i) {
                const Pixel& p = pixels[lo + i];
                b1_[i] = p.bin1;
                b2_[i] = p.bin2;
                c_[i] = p.count;
                bin1_counts_[static_cast<std::size_t>(p.bin1)]++;
            }
            bin1_.append(b1_.data(), count);
            bin2_.append(b2_.data(), count);
            count_.append(c_.data(), count);
            nnz_ += count;
            // A fixed flush size, not one derived from the thread count: where
            // HDF5 places the chunk index nodes follows the flushes, and the
            // file is to be the same byte for byte on any number of threads.
            if (bin1_.buffered_chunks() >= kFlushChunks) {
                h5::flush(pool_, {&bin1_, &bin2_, &count_}, false);
            }
        }
    }

    void finish() { h5::flush(pool_, {&bin1_, &bin2_, &count_}, true); }

    [[nodiscard]] std::size_t nnz() const noexcept { return nnz_; }
    [[nodiscard]] const std::vector<std::int64_t>& bin1_counts() const noexcept { return bin1_counts_; }

  private:
    ThreadPool& pool_;
    h5::ChunkedWriter bin1_;
    h5::ChunkedWriter bin2_;
    h5::ChunkedWriter count_;
    std::vector<std::int64_t> bin1_counts_;
    std::vector<std::int64_t> b1_;
    std::vector<std::int64_t> b2_;
    std::vector<std::int32_t> c_;
    std::size_t nnz_ = 0;
};

void convert_resolution(ThreadPool& pool, const hicfilecpp::HiCFile& hic, const HicHeader& header, h5::File& file,
                        std::int32_t binsize, bool multi_res, bool square, bool show_warnings,
                        const Console& console) {
    std::string group = "/";
    if (multi_res) {
        if (!file.exists("/resolutions")) {
            file.set_attribute("/", "format", std::string("HDF5::MCOOL"));
            file.set_attribute("/", "format-version", std::int64_t{2});
            file.create_group("/resolutions");
        }
        group = "/resolutions/" + std::to_string(binsize);
        file.create_group(group);
    }
    const std::string prefix = group == "/" ? "/" : group + "/";
    const auto& used = header.used;
    std::vector<const HicChrom*> chroms;
    std::map<std::int32_t, std::string> name_of;
    for (const auto& chrom : used) {
        name_of[chrom.index] = chrom.name;
        if (lowercase(chrom.name) != "all") {
            chroms.push_back(&chrom);
        }
    }

    // write_chroms
    file.create_group(prefix + "chroms");
    std::vector<char> names(chroms.size() * kChromNameWidth, '\0');
    std::vector<std::int32_t> lengths;
    std::vector<std::string> enum_names;
    for (std::size_t i = 0; i < chroms.size(); ++i) {
        const std::string name = chroms[i]->name.substr(0, kChromNameWidth);
        std::copy(name.begin(), name.end(), names.begin() + static_cast<std::ptrdiff_t>(i * kChromNameWidth));
        enum_names.push_back(name);
        lengths.push_back(static_cast<std::int32_t>(chroms[i]->length));
    }
    {
        const h5::ColumnType string_type = h5::ColumnType::fixed_string(kChromNameWidth);
        h5::write_dataset(pool, file, prefix + "chroms/name", string_type, names.data(), chroms.size());
    }
    h5::write_dataset(pool, file, prefix + "chroms/length", h5::ColumnType::int32(), lengths.data(), lengths.size());

    // create_bins
    std::vector<std::int32_t> chrom_ids;
    std::vector<std::int32_t> starts;
    std::vector<std::int32_t> ends;
    std::vector<std::int64_t> chrom_offsets{0};
    std::map<std::int32_t, std::int64_t> bins_of;
    std::map<std::int32_t, std::int64_t> offset_of;
    std::int32_t position = 0;
    for (const HicChrom* chrom : chroms) {
        std::int64_t count = 0;
        for (std::int64_t start = 0; start < chrom->length;) {
            const std::int64_t end = std::min(start + binsize, chrom->length);
            chrom_ids.push_back(position);
            starts.push_back(static_cast<std::int32_t>(start));
            ends.push_back(static_cast<std::int32_t>(end));
            count++;
            start = end;
        }
        offset_of[chrom->index] = chrom_offsets.back();
        bins_of[chrom->index] = count;
        if (chrom_offsets.size() < used.size()) {
            chrom_offsets.push_back(
                static_cast<std::int64_t>(std::ceil(static_cast<double>(chrom->length) / binsize)) +
                chrom_offsets.back());
        }
        position++;
    }
    const std::size_t n_bins = chrom_ids.size();
    if (n_bins > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw ExitError("!!! ERROR. " + std::to_string(n_bins) + " bins exceed what this converter indexes");
    }

    // write_bins
    file.create_group(prefix + "bins");
    {
        const h5::ColumnType chrom_type = h5::ColumnType::enumeration(enum_names);
        h5::write_dataset(pool, file, prefix + "bins/chrom", chrom_type, chrom_ids.data(), n_bins);
    }
    h5::write_dataset(pool, file, prefix + "bins/start", h5::ColumnType::int32(), starts.data(), n_bins);
    h5::write_dataset(pool, file, prefix + "bins/end", h5::ColumnType::int32(), ends.data(), n_bins);
    for (const auto& norm : header.norms) {
        std::vector<double> column;
        column.reserve(n_bins);
        for (const HicChrom* chrom : chroms) {
            const auto chr_bins = static_cast<std::size_t>(bins_of[chrom->index]);
            const auto vector = hic.readNormVector(norm, chrom->index, "BP", binsize);
            if (!vector) {
                if (show_warnings) {
                    console.err("!!! WARNING. Normalization vector " + norm + " does not exist for chr idx " +
                                std::to_string(chrom->index) + ".");
                }
                column.insert(column.end(), chr_bins, std::numeric_limits<double>::quiet_NaN());
            } else {
                column.insert(column.end(), vector->begin(),
                              vector->begin() + static_cast<std::ptrdiff_t>(std::min(vector->size(), chr_bins)));
            }
        }
        if (column.size() != n_bins) {
            throw ExitError("!!! ERROR. Length of normalization vector " + norm +
                            " does not match the number of bins.\nThis is likely a problem with the hic file");
        }
        h5::write_dataset(pool, file, prefix + "bins/" + norm, h5::ColumnType::float64(), column.data(), n_bins);
    }

    // write_chrom_offset and initialize_pixels
    file.create_group(prefix + "indexes");
    h5::write_dataset(pool, file, prefix + "indexes/chrom_offset", h5::ColumnType::int64(), chrom_offsets.data(),
                      chrom_offsets.size());
    file.create_group(prefix + "pixels");
    PixelWriter writer(pool, file, prefix, n_bins);

    for (const auto& [key, value] : header.metadata) {
        file.set_attribute(group, key, value);
    }
    file.set_attribute(group, "nchroms", static_cast<std::int64_t>(chroms.size()));
    file.set_attribute(group, "nbins", static_cast<std::int64_t>(n_bins));
    file.set_attribute(group, "bin-type", std::string("fixed"));
    file.set_attribute(group, "bin-size", static_cast<std::int64_t>(binsize));
    file.set_attribute(group, "format", std::string("HDF5::Cooler"));
    file.set_attribute(group, "format-url", std::string(kFormatUrl));
    file.set_attribute(group, "format-version", std::int64_t{3});
    file.set_attribute(group, "storage-mode", std::string(square ? "square" : "symmetric-upper"));
    file.set_attribute(group, "generated-by", std::string("hic2cool-") + kVersion);
    file.set_attribute(group, "genome-assembly", header.genome);
    file.set_attribute(group, "creation-date", utcnow_isoformat());

    // The pixel loop: every chromosome pair once, grouped by the first
    // chromosome of the pair in file order, each group in (bin1, bin2) order.
    // hic2cool concatenates a group and sorts it in memory. Here the group is
    // streamed: blocks are visited in the order of the smallest bin1 each can
    // hold, and pending pixels below the next block's bound are final, so they
    // are sorted and written. Blocks are decoded on the pool in batches.
    std::set<std::pair<std::int32_t, std::int32_t>> covered;
    std::vector<Pixel> pending;
    const std::size_t batch = static_cast<std::size_t>(pool.threads());
    const auto emit = [&](std::size_t count) {
        if (count == 0) {
            return;
        }
        parallel_sort(pool, pending.data(), count);
        writer.append(pending.data(), count);
        pending.erase(pending.begin(), pending.begin() + static_cast<std::ptrdiff_t>(count));
    };
    for (const HicChrom* chr_a : chroms) {
        struct PairMatrix {
            hicfilecpp::MatrixZoomData mzd;
            std::int64_t bins1;
            std::int64_t bins2;
            std::int64_t offset1;
            std::int64_t offset2;
            bool mirror;  // square storage: the pair transposed, its pixels in chr_a's rows as bin1
        };
        struct BlockRef {
            std::int64_t floor;
            std::size_t pair;
            hicfilecpp::BlockIndexEntry entry;
        };
        std::vector<PairMatrix> pairs;
        std::vector<BlockRef> blocks;
        const auto add_pair = [&](std::int32_t c1, std::int32_t c2, bool mirror) {
            const auto headers = hic.hasMatrix(c1, c2) ? hic.matrixZoomHeaders(c1, c2)
                                                        : std::vector<hicfilecpp::ZoomHeader>{};
            const auto zoom = std::find_if(headers.begin(), headers.end(), [&](const hicfilecpp::ZoomHeader& h) {
                return h.unit == "BP" && h.binSize == binsize;
            });
            if (zoom == headers.end()) {
                if (show_warnings && !mirror) {
                    console.err("... The intersection between " + name_of[c1] + " and " + name_of[c2] +
                                " cannot be found in the hic file.");
                }
                return;
            }
            pairs.push_back(PairMatrix{hic.getMatrixZoomData(name_of[c1], name_of[c2], "observed", "NONE", "BP", binsize),
                                       bins_of[c1], bins_of[c2], offset_of[c1], offset_of[c2], mirror});
            const bool v9_intra = hic.version() > 8 && c1 == c2;
            for (const auto& entry : pairs.back().mzd.blockIndex()) {
                std::int64_t floor = 0;
                if (!mirror) {
                    floor = pairs.back().offset1 +
                            block_bin1_floor(entry.number, zoom->blockBinCount, zoom->blockColumnCount, v9_intra);
                } else if (v9_intra) {
                    // binY >= (binX + binY) / 2, so the bound of binX holds for binY
                    floor = pairs.back().offset2 +
                            block_bin1_floor(entry.number, zoom->blockBinCount, zoom->blockColumnCount, true);
                } else {
                    // the block's row: binY / blockBinCount
                    floor = pairs.back().offset2 + static_cast<std::int64_t>(entry.number / zoom->blockColumnCount) *
                                                       zoom->blockBinCount;
                }
                blocks.push_back(BlockRef{floor, pairs.size() - 1, entry});
            }
        };
        for (const HicChrom* chr_b : chroms) {
            const std::int32_t c1 = std::min(chr_a->index, chr_b->index);
            const std::int32_t c2 = std::max(chr_a->index, chr_b->index);
            if (!covered.insert({c1, c2}).second) {
                continue;
            }
            add_pair(c1, c2, false);
        }
        if (square) {
            // The lower triangle of chr_a's rows: every pair (c, chr_a) with c
            // up to chr_a, transposed, the diagonal left out.
            for (const HicChrom* chr_b : chroms) {
                if (chr_b->index <= chr_a->index) {
                    add_pair(chr_b->index, chr_a->index, true);
                }
            }
        }
        std::sort(blocks.begin(), blocks.end(), [](const BlockRef& a, const BlockRef& b) {
            if (a.floor != b.floor) {
                return a.floor < b.floor;
            }
            return a.pair != b.pair ? a.pair < b.pair : a.entry.number < b.entry.number;
        });
        std::int64_t final_below = std::numeric_limits<std::int64_t>::min();
        for (std::size_t first = 0; first < blocks.size(); first += batch) {
            const std::size_t n = std::min(batch, blocks.size() - first);
            std::vector<std::vector<Pixel>> decoded(n);
            std::vector<std::int64_t> lowest(n, std::numeric_limits<std::int64_t>::max());
            pool.run(n, [&](std::size_t k) {
                const BlockRef& block = blocks[first + k];
                const PairMatrix& pair = pairs[block.pair];
                for (const auto& record : pair.mzd.readBlock(block.entry)) {
                    if (record.binX >= 0 && record.binX < pair.bins1 && record.binY >= 0 &&
                        record.binY < pair.bins2) {
                        std::int64_t bin1 = record.binX + pair.offset1;
                        std::int64_t bin2 = record.binY + pair.offset2;
                        if (pair.mirror) {
                            if (bin1 == bin2) {
                                continue;
                            }
                            std::swap(bin1, bin2);
                        }
                        lowest[k] = std::min(lowest[k], bin1);
                        decoded[k].push_back(Pixel{static_cast<std::int32_t>(bin1), static_cast<std::int32_t>(bin2),
                                                   numpy_int32(record.counts)});
                    }
                }
            });
            for (std::size_t k = 0; k < n; ++k) {
                const BlockRef& block = blocks[first + k];
                if (block.floor > final_below) {
                    const auto split = std::partition(pending.begin(), pending.end(),
                                                      [&](const Pixel& p) { return p.bin1 < block.floor; });
                    emit(static_cast<std::size_t>(split - pending.begin()));
                    final_below = block.floor;
                }
                if (lowest[k] < final_below) {
                    throw ExitError("!!! ERROR. Block " + std::to_string(block.entry.number) +
                                    " holds bins outside the range its block number allows; the hic file is not "
                                    "laid out as Juicer tools writes it");
                }
                pending.insert(pending.end(), decoded[k].begin(), decoded[k].end());
                std::vector<Pixel>().swap(decoded[k]);
            }
        }
        emit(pending.size());
        std::vector<Pixel>().swap(pending);
    }
    writer.finish();

    // finalize_resolution_cool
    std::vector<std::int64_t> bin1_offset(n_bins + 1, 0);
    for (std::size_t i = 0; i < n_bins; ++i) {
        bin1_offset[i + 1] = bin1_offset[i] + writer.bin1_counts()[i];
    }
    h5::write_dataset(pool, file, prefix + "indexes/bin1_offset", h5::ColumnType::int64(), bin1_offset.data(),
                      bin1_offset.size());
    file.set_attribute(group, "nnz", static_cast<std::int64_t>(writer.nnz()));
}

}  // namespace

namespace {

// resolutions empty: every resolution of the file. exact_multi: a multi
// resolution layout written to exactly `outfile`, even for one resolution.
std::string convert_impl(const std::string& infile, const std::string& outfile,
                         const std::vector<std::int64_t>& resolutions, bool all_when_empty, bool exact_multi,
                         int nproc, bool show_warnings, bool silent, const std::string& storage_mode,
                         const Console& console) {
    if (storage_mode != "symmetric-upper" && storage_mode != "square") {
        throw ExitError("!!! ERROR. storage mode must be symmetric-upper or square, not " + storage_mode);
    }
    const bool square = storage_mode == "square";
    check_hic_magic(infile);
    const hicfilecpp::HiCFile hic(infile);
    const HicHeader header = read_hic_header(hic, console);
    if (!silent) {
        std::vector<std::string> chr_names;
        for (const auto& chrom : header.used) {
            chr_names.push_back(chrom.name);
        }
        console.out("##########################");
        console.out("### hic2cool / convert ###");
        console.out("##########################");
        console.out("### Header info from hic");
        console.out("... Chromosomes:  " + py_repr(chr_names));
        console.out("... Resolutions:  " + py_repr(header.resolutions));
        console.out("... Normalizations:  " + py_repr(header.norms));
        console.out("... Genome:  " + header.genome);
    }
    std::vector<std::int32_t> use;
    if (resolutions.empty() && all_when_empty) {
        use = header.resolutions;
    } else {
        for (const std::int64_t resolution : resolutions) {
            if (std::find(header.resolutions.begin(), header.resolutions.end(), resolution) ==
                header.resolutions.end()) {
                throw ExitError("!!! ERROR. Given binsize (in bp) is not a supported resolution in this file.\n"
                                "Please use 0 (all resolutions) or use one of: " +
                                py_repr(header.resolutions));
            }
            use.push_back(static_cast<std::int32_t>(resolution));
        }
    }
    const bool multi_res = exact_multi || use.size() > 1;
    const std::string written = exact_multi ? outfile : output_name(outfile, multi_res);
    std::error_code error;
    if (std::filesystem::exists(written, error)) {
        if (!std::filesystem::remove(written, error) || error) {
            throw ExitError("!!! ERROR. Output file path " + written +
                            " already exists. This can cause issues with the hdf5 structure. Please remove that "
                            "file or choose a different output name.");
        }
    }
    console.out("### Converting");
    ThreadPool pool(nproc > 0 ? nproc : available_threads());
    h5::File file(written, h5::Mode::Create);
    for (const std::int32_t binsize : use) {
        const auto start = std::chrono::steady_clock::now();
        convert_resolution(pool, hic, header, file, binsize, multi_res, square, show_warnings, console);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (!silent) {
            console.out("... Resolution " + std::to_string(binsize) + " took: " + py_repr(elapsed) + " seconds.");
        }
    }
    file.close();
    if (!silent) {
        console.out("### Finished! Output written to: " + written);
        if (multi_res) {
            console.out("... This file is higlass compatible.");
        } else {
            console.out("... This file is single resolution and NOT higlass compatible. Run with `-r 0` for "
                        "multi-resolution.");
        }
    }
    return written;
}

}  // namespace

std::string hic2cool_convert(const std::string& infile, const std::string& outfile, std::int64_t resolution,
                             int nproc, bool show_warnings, bool silent, const std::string& storage_mode,
                             const Console& console) {
    const std::vector<std::int64_t> resolutions =
        resolution == 0 ? std::vector<std::int64_t>{} : std::vector<std::int64_t>{resolution};
    return convert_impl(infile, outfile, resolutions, true, false, nproc, show_warnings, silent, storage_mode,
                        console);
}

std::string hic2cool_convert_mcool(const std::string& infile, const std::string& outfile,
                                   const std::vector<std::int64_t>& resolutions, int nproc, bool show_warnings,
                                   bool silent, const Console& console) {
    return convert_impl(infile, outfile, resolutions, true, true, nproc, show_warnings, silent,
                        "symmetric-upper", console);
}

}  // namespace hic2cool
