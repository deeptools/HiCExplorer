// cool2hic: a .cool or .mcool file to a .hic file, the opposite of
// hic2cool_convert. Not part of the Python hic2cool package.
//
// The pixels of every chosen resolution are handed to hicfilecpp's writer,
// which lays the file out as Juicer tools pre does (version 8 as release
// 1.22.01, version 9 as 2.x). Normalization vectors that hic2cool stored as
// bins table columns are written back as they are, so that a .hic file taken
// through hic2cool and cool2hic keeps its vectors; otherwise Juicer's VC,
// VC_SQRT, KR and SCALE are computed.

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <cstring>
#include <sstream>

#include <hicfilecpp/hicfilecpp.hpp>

#include <coolercpp/coolercpp.hpp>

#include "io.hpp"
#include "hic2cool/hic2cool.hpp"
#include "pyformat.hpp"

namespace hic2cool {

namespace {

using namespace detail;

constexpr std::size_t kReadSlice = std::size_t{1} << 24;

const std::set<std::string> kJuicerNormNames{"VC",     "VC_SQRT",  "KR",       "SCALE",      "GW_KR",
                                             "GW_VC",  "GW_SCALE", "INTER_KR", "INTER_VC", "INTER_SCALE"};

// Attributes cooler and hic2cool write themselves; the others of a cool file
// written by hic2cool come from the .hic header and go back into it.
const std::set<std::string> kCoolerAttributes{
    "nchroms",      "nbins",        "nnz",          "bin-type",        "bin-size",      "format",
    "format-url",   "format-version", "storage-mode", "generated-by",  "genome-assembly", "creation-date",
    "update-date",  "metadata",     "sum",          "cis",             "software"};

struct Resolution {
    std::string path;
    std::string prefix;
    std::int32_t binsize = 0;
    std::vector<std::int64_t> chrom_offset;  // bins, per chromosome, n + 1
    std::vector<std::int64_t> bin1_offset;
    std::vector<std::string> columns;       // norm columns carried from this group
    bool square = false;                    // storage-mode square: both triangles stored
};

std::string text_of(const std::optional<h5::Value>& value) {
    if (!value) {
        return {};
    }
    if (const auto* s = std::get_if<std::string>(&*value)) {
        return *s;
    }
    return {};
}

std::vector<std::string> split_list(const std::string& text) {
    std::vector<std::string> out;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if (!item.empty()) {
            out.push_back(item);
        }
    }
    return out;
}

class CoolSource : public hicfilecpp::PixelSource {
  public:
    // lower: take the lower triangle of square coolers, transposed.
    CoolSource(ThreadPool& pool, const h5::File& file, std::vector<Resolution>& resolutions, bool lower)
        : pool_(pool), file_(file), resolutions_(resolutions), lower_(lower) {}

    // The group whose pixels give `resolution`: its own, or the finest one,
    // binned.
    const Resolution& group_for(std::int32_t resolution) const {
        for (const auto& r : resolutions_) {
            if (r.binsize == resolution) {
                return r;
            }
        }
        return *std::min_element(resolutions_.begin(), resolutions_.end(),
                                 [](const Resolution& a, const Resolution& b) { return a.binsize < b.binsize; });
    }

    void pixels(std::int32_t resolution, std::int32_t chr1, std::int32_t chr2,
                const std::function<void(const hicfilecpp::Pixel*, std::size_t)>& consume) override {
        const Resolution& group = group_for(resolution);
        // The writer asks for one pair at a time, once per resolution; a
        // pair's pixels are kept until it moves on to the next pair.
        if (chr1 != last_chr1_ || chr2 != last_chr2_) {
            release(last_chr1_, last_chr2_);
            last_chr1_ = chr1;
            last_chr2_ = chr2;
        }
        // Upper: the rows of chr1, pixels whose bin2 lies on chr2. Lower: the
        // rows of chr2, pixels whose bin2 lies on chr1, transposed.
        const std::int32_t row_chrom = lower_ ? chr2 : chr1;
        const std::int32_t other = lower_ ? chr1 : chr2;
        auto& loaded = caches_[group.path];
        auto it = loaded.find(row_chrom);
        if (it == loaded.end()) {
            if (!lower_) {
                loaded.clear();  // upper: one row chromosome at a time
            }
            it = loaded.emplace(row_chrom, load(group, row_chrom)).first;
        }
        const auto& bucket = it->second[static_cast<std::size_t>(other)];
        if (bucket.empty()) {
            return;
        }
        if (group.binsize == resolution || resolution % group.binsize != 0) {
            consume(bucket.data(), bucket.size());
            return;
        }
        const std::int32_t factor = resolution / group.binsize;
        std::vector<hicfilecpp::Pixel> binned(bucket.begin(), bucket.end());
        for (auto& p : binned) {
            p.bin1 /= factor;
            p.bin2 /= factor;
        }
        consume(binned.data(), binned.size());
    }

    [[nodiscard]] std::size_t skipped() const noexcept { return skipped_; }

    // Whether a square cooler's lower triangle mirrors its upper triangle:
    // the number of off-diagonal pixels and an order-independent 64-bit
    // hash of (smaller bin, larger bin, count) agree between the triangles.
    bool symmetric(const Resolution& group) {
        std::uint64_t upper_hash = 0;
        std::uint64_t lower_hash = 0;
        std::uint64_t upper_count = 0;
        std::uint64_t lower_count = 0;
        const std::size_t n_chroms = group.chrom_offset.size() - 1;
        for (std::size_t chrom = 0; chrom < n_chroms; ++chrom) {
            scan_rows(group, static_cast<std::int32_t>(chrom), [&](std::int64_t b1, std::int64_t b2, double value) {
                if (b1 == b2 || value == 0) {
                    return;
                }
                std::uint64_t h = mix(static_cast<std::uint64_t>(std::min(b1, b2)));
                h = mix(h ^ static_cast<std::uint64_t>(std::max(b1, b2)));
                std::uint64_t bits = 0;
                std::memcpy(&bits, &value, sizeof(bits));
                h = mix(h ^ bits);
                if (b1 < b2) {
                    upper_hash += h;
                    ++upper_count;
                } else {
                    lower_hash += h;
                    ++lower_count;
                }
            });
        }
        return upper_hash == lower_hash && upper_count == lower_count;
    }

  private:
    using Buckets = std::vector<std::vector<hicfilecpp::Pixel>>;

    static std::uint64_t mix(std::uint64_t x) {  // splitmix64
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    void release(std::int32_t chr1, std::int32_t chr2) {
        if (chr1 < 0) {
            return;
        }
        const std::int32_t row_chrom = lower_ ? chr2 : chr1;
        const std::int32_t other = lower_ ? chr1 : chr2;
        for (auto& [path, loaded] : caches_) {
            const auto it = loaded.find(row_chrom);
            if (it != loaded.end()) {
                std::vector<hicfilecpp::Pixel>().swap(it->second[static_cast<std::size_t>(other)]);
            }
        }
    }

    // Visits (bin1, bin2, count) of every pixel in the rows of `chrom`.
    template <class Visit>
    void scan_rows(const Resolution& group, std::int32_t chrom, Visit visit) {
        const std::int64_t first_bin = group.chrom_offset[static_cast<std::size_t>(chrom)];
        const std::int64_t end_bin = group.chrom_offset[static_cast<std::size_t>(chrom) + 1];
        const auto lo = static_cast<std::size_t>(group.bin1_offset[static_cast<std::size_t>(first_bin)]);
        const auto hi = static_cast<std::size_t>(group.bin1_offset[static_cast<std::size_t>(end_bin)]);
        if (hi <= lo) {
            return;
        }
        auto& readers = readers_[group.path];
        if (!readers.first) {
            readers.first = std::make_unique<h5::ColumnReader>(file_, group.prefix + "pixels/bin2_id");
            readers.second = std::make_unique<h5::ColumnReader>(file_, group.prefix + "pixels/count");
        }
        std::vector<std::int64_t> bin2(std::min(kReadSlice, hi - lo));
        std::vector<double> count(bin2.size());
        std::int64_t bin1 = first_bin;
        for (std::size_t start = lo; start < hi; start += kReadSlice) {
            const std::size_t n = std::min(kReadSlice, hi - start);
            readers.first->read(pool_, start, start + n, bin2.data());
            readers.second->read(pool_, start, start + n, count.data());
            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t row = start + k;
                while (bin1 < end_bin &&
                       static_cast<std::size_t>(group.bin1_offset[static_cast<std::size_t>(bin1) + 1]) <= row) {
                    ++bin1;
                }
                visit(bin1, bin2[k], count[k]);
            }
        }
    }

    // The pixels of the rows of `chrom`, by the chromosome of their other bin:
    // for the upper triangle bin2 on the same or a later chromosome, for the
    // lower triangle bin2 on the same or an earlier one, transposed.
    Buckets load(const Resolution& group, std::int32_t chrom) {
        const std::size_t n_chroms = group.chrom_offset.size() - 1;
        const std::int64_t first_bin = group.chrom_offset[static_cast<std::size_t>(chrom)];
        Buckets buckets(n_chroms);
        scan_rows(group, chrom, [&](std::int64_t b1, std::int64_t b2, double value) {
            const auto it = std::upper_bound(group.chrom_offset.begin(), group.chrom_offset.end(), b2);
            const auto chr2 = static_cast<std::int32_t>(it - group.chrom_offset.begin()) - 1;
            if (chr2 < 0 || chr2 >= static_cast<std::int32_t>(n_chroms)) {
                return;
            }
            const bool upper = chr2 > chrom || (chr2 == chrom && b2 >= b1);
            const bool lower = chr2 < chrom || (chr2 == chrom && b2 <= b1);
            if (lower_ ? !lower : !upper) {
                return;
            }
            if (value == 0) {
                return;
            }
            if (!std::isfinite(value)) {
                ++skipped_;
                return;
            }
            const auto local1 = static_cast<std::int32_t>(b1 - first_bin);
            const auto local2 = static_cast<std::int32_t>(b2 - group.chrom_offset[static_cast<std::size_t>(chr2)]);
            buckets[static_cast<std::size_t>(chr2)].push_back(
                lower_ ? hicfilecpp::Pixel{local2, local1, static_cast<float>(value)}
                       : hicfilecpp::Pixel{local1, local2, static_cast<float>(value)});
        });
        return buckets;
    }

    ThreadPool& pool_;
    const h5::File& file_;
    std::vector<Resolution>& resolutions_;
    bool lower_;
    std::map<std::string, std::map<std::int32_t, Buckets>> caches_;
    std::map<std::string, std::pair<std::unique_ptr<h5::ColumnReader>, std::unique_ptr<h5::ColumnReader>>> readers_;
    std::int32_t last_chr1_ = -1;
    std::int32_t last_chr2_ = -1;
    std::size_t skipped_ = 0;
};

}  // namespace

std::string cool2hic_convert(const std::string& infile, const std::string& outfile, const Cool2hicOptions& options,
                             const Console& console) {
    if (options.hic_version != 8 && options.hic_version != 9) {
        throw ExitError("!!! ERROR. The .hic version must be 8 or 9, not " + std::to_string(options.hic_version));
    }
    const auto separator = infile.find("::");
    const std::string path = infile.substr(0, separator);
    std::string root = separator == std::string::npos ? "" : infile.substr(separator + 2);
    if (!root.empty() && root[0] != '/') {
        root = "/" + root;
    }
    if (!h5::is_hdf5(path)) {
        throw ExitError("!!! ERROR. " + path + " is not a cool or mcool file");
    }
    const std::vector<std::string> cooler_paths = coolercpp::list_coolers(path);
    const h5::File file(path, h5::Mode::Read);
    std::vector<std::string> groups;
    for (const auto& group : cooler_paths) {
        if (root.empty() || root == "/" || group == root || group.rfind(root + "/", 0) == 0) {
            groups.push_back(group);
        }
    }
    if (!root.empty() && root != "/" && std::find(groups.begin(), groups.end(), root) != groups.end()) {
        groups = {root};
    }
    if (groups.empty()) {
        // cooler's own check compares the format attribute with a str, so a
        // cool file whose attributes were written as bytes (HiCMatrix up to
        // version 10, h5py 2) lists no cooler although its group is complete.
        const std::string candidate = root.empty() ? "/" : root;
        const std::string prefix = candidate == "/" ? "" : candidate;
        bool complete = file.is_group(candidate);
        for (const char* name : {"chroms", "bins", "pixels", "indexes"}) {
            complete = complete && file.exists(prefix + "/" + name);
        }
        if (complete) {
            groups.push_back(candidate);
        }
    }
    if (groups.empty()) {
        throw ExitError("!!! ERROR. No cooler found in " + infile);
    }

    // The resolutions of the input.
    std::vector<Resolution> available;
    for (const auto& group : groups) {
        Resolution r;
        r.path = group;
        r.prefix = group == "/" ? "/" : group + "/";
        const auto size = file.attribute(group, "bin-size");
        const std::string bin_type = text_of(file.attribute(group, "bin-type"));
        if (!size || !std::holds_alternative<std::int64_t>(*size) || (!bin_type.empty() && bin_type != "fixed")) {
            throw ExitError("!!! ERROR. The cooler " + group + " in " + path +
                            " does not have fixed size bins, which the .hic format needs");
        }
        const std::int64_t binsize = std::get<std::int64_t>(*size);
        if (binsize <= 0 || binsize > std::numeric_limits<std::int32_t>::max()) {
            throw ExitError("!!! ERROR. Unusable bin size " + std::to_string(binsize) + " in " + group);
        }
        r.binsize = static_cast<std::int32_t>(binsize);
        r.square = text_of(file.attribute(group, "storage-mode")) == "square";
        available.push_back(std::move(r));
    }
    std::sort(available.begin(), available.end(),
              [](const Resolution& a, const Resolution& b) { return a.binsize > b.binsize; });

    std::vector<Resolution> chosen;
    if (options.resolution == 0) {
        chosen = available;
    } else {
        for (const auto& r : available) {
            if (r.binsize == options.resolution) {
                chosen = {r};
                break;
            }
        }
        if (chosen.empty()) {
            std::vector<std::int64_t> sizes;
            for (const auto& r : available) {
                sizes.push_back(r.binsize);
            }
            throw ExitError("!!! ERROR. Given binsize (in bp) is not a resolution of this file.\n"
                            "Please use 0 (all resolutions) or use one of: " +
                            py_repr(sizes));
        }
    }
    // Duplicate resolutions (two groups with the same bin size): the last one
    // in cooler's listing order wins, as in hic2cool extract-norms.
    for (auto it = chosen.begin(); it != chosen.end();) {
        const auto later = std::find_if(it + 1, chosen.end(), [&](const Resolution& r) { return r.binsize == it->binsize; });
        it = later != chosen.end() ? chosen.erase(it) : it + 1;
    }
    const std::int32_t finest = chosen.back().binsize;

    // Chromosomes: shared by every chosen resolution.
    std::vector<std::string> names = file.read_strings(chosen.front().prefix + "chroms/name");
    std::vector<std::int64_t> lengths = file.read_int64(chosen.front().prefix + "chroms/length");
    for (auto& r : chosen) {
        if (file.read_strings(r.prefix + "chroms/name") != names) {
            throw ExitError("!!! ERROR. The resolutions of " + path + " do not have the same chromosomes");
        }
        const auto other = file.read_int64(r.prefix + "chroms/length");
        for (std::size_t i = 0; i < lengths.size(); ++i) {
            lengths[i] = std::max(lengths[i], other[i]);
        }
        r.chrom_offset = file.read_int64(r.prefix + "indexes/chrom_offset");
        r.bin1_offset = file.read_int64(r.prefix + "indexes/bin1_offset");
        if (r.chrom_offset.size() != names.size() + 1) {
            throw ExitError("!!! ERROR. Malformed chrom_offset index in " + r.path);
        }
    }

    // Normalizations.
    std::vector<std::string> computed;
    bool carry = false;
    std::string mode = options.normalizations;
    if (mode == "auto") {
        carry = true;
    } else if (mode != "none") {
        computed = split_list(mode);
    }
    std::vector<std::string> carried;
    for (auto& r : chosen) {
        if (!carry) {
            break;
        }
        const bool from_hic2cool = text_of(file.attribute(r.path, "generated-by")).rfind("hic2cool", 0) == 0;
        const std::size_t n_bins = static_cast<std::size_t>(r.chrom_offset.back());
        for (const auto& column : file.children(r.prefix + "bins")) {
            if (column == "chrom" || column == "start" || column == "end" || column == "weight") {
                continue;
            }
            if (!from_hic2cool && kJuicerNormNames.count(column) == 0) {
                continue;
            }
            if (!file.is_dataset(r.prefix + "bins/" + column) || file.length(r.prefix + "bins/" + column) != n_bins) {
                if (options.show_warnings && !options.silent) {
                    console.err("!!! WARNING. Bins column " + column + " of " + r.path +
                                " is not a normalization vector; not written.");
                }
                continue;
            }
            r.columns.push_back(column);
            if (std::find(carried.begin(), carried.end(), column) == carried.end()) {
                carried.push_back(column);
            }
        }
    }
    if (carry && carried.empty()) {
        computed = {"VC", "VC_SQRT", "KR", "SCALE"};
    }
    std::vector<std::string> provided = carried;
    if (options.cooler_weight) {
        const bool any = std::any_of(chosen.begin(), chosen.end(),
                                     [&](const Resolution& r) { return file.is_dataset(r.prefix + "bins/weight"); });
        if (any) {
            provided.push_back(*options.cooler_weight);
        } else if (options.show_warnings && !options.silent) {
            console.err("!!! WARNING. No resolution of " + infile + " has a weight column; " +
                        *options.cooler_weight + " is not written.");
        }
    }

    // Resolutions to write.
    std::vector<std::int32_t> write_resolutions;
    for (const auto& r : chosen) {
        write_resolutions.push_back(r.binsize);
    }
    for (const std::int64_t extra : options.add_resolutions) {
        if (extra <= 0 || extra % finest != 0 || extra > std::numeric_limits<std::int32_t>::max()) {
            throw ExitError("!!! ERROR. Added resolution " + std::to_string(extra) +
                            " is not a multiple of the finest resolution written, " + std::to_string(finest));
        }
        if (std::find(write_resolutions.begin(), write_resolutions.end(), extra) == write_resolutions.end()) {
            write_resolutions.push_back(static_cast<std::int32_t>(extra));
        }
    }

    const std::string genome = !options.genome.empty() ? options.genome
                               : !text_of(file.attribute(chosen.back().path, "genome-assembly")).empty()
                                   ? text_of(file.attribute(chosen.back().path, "genome-assembly"))
                                   : std::string("unknown");
    if (!options.silent) {
        std::vector<std::int32_t> sorted = write_resolutions;
        std::sort(sorted.begin(), sorted.end(), std::greater<>());
        console.out("##########################");
        console.out("### cool2hic / convert ###");
        console.out("##########################");
        console.out("### Header info from cool");
        console.out("... Chromosomes:  " + py_repr(names));
        console.out("... Resolutions:  " + py_repr(sorted));
        console.out("... Normalizations carried:  " + py_repr(provided));
        console.out("... Normalizations computed:  " + py_repr(computed));
        console.out("... Genome:  " + genome);
        console.out("... hic version:  " + std::to_string(options.hic_version));
        if (std::any_of(chosen.begin(), chosen.end(), [](const Resolution& r) { return r.square; })) {
            console.out("... Square cooler, triangle:  " + options.triangle);
        }
        console.out("### Converting");
    }

    const int threads = options.nproc > 0 ? options.nproc : available_threads();
    ThreadPool pool(threads);
    hicfilecpp::WriteOptions write;
    write.version = options.hic_version;
    write.genomeId = genome;
    for (std::size_t i = 0; i < names.size(); ++i) {
        write.chromosomes.emplace_back(names[i], lengths[i]);
    }
    write.resolutions = write_resolutions;
    write.threads = threads;
    write.normalizations = computed;
    write.providedNormalizations = provided;
    write.software = std::string("cool2hic (hic2cool ") + kVersion + ")";
    for (const auto& [key, value] : file.attributes(chosen.back().path)) {
        if (kCoolerAttributes.count(key) == 0 && std::holds_alternative<std::string>(value)) {
            write.attributes.emplace_back(key, std::get<std::string>(value));
        }
    }
    if (chosen.size() == 1) {
        write.sourceResolution = finest;
    } else {
        write.sourceProvidesEveryResolution = true;
    }

    // Norm columns, read once per group and column.
    std::map<std::pair<std::string, std::string>, std::vector<double>> columns;
    write.normVector = [&](const std::string& name, std::int32_t chr, std::int32_t resolution) -> std::vector<double> {
        const auto group = std::find_if(chosen.begin(), chosen.end(),
                                        [&](const Resolution& r) { return r.binsize == resolution; });
        if (group == chosen.end()) {
            return {};
        }
        const bool weight = options.cooler_weight && name == *options.cooler_weight;
        const std::string column = weight ? "weight" : name;
        if (!weight && std::find(group->columns.begin(), group->columns.end(), column) == group->columns.end()) {
            return {};
        }
        const auto key = std::make_pair(group->path, column);
        auto it = columns.find(key);
        if (it == columns.end()) {
            const std::string dataset = group->prefix + "bins/" + column;
            std::vector<double> values;
            if (file.is_dataset(dataset)) {
                values = file.read_doubles(dataset);
                if (weight) {
                    for (double& v : values) {
                        v = 1.0 / v;  // cooler multiplies by weights, Juicer divides by vectors
                    }
                }
            }
            it = columns.emplace(key, std::move(values)).first;
        }
        const auto& values = it->second;
        const auto lo = static_cast<std::size_t>(group->chrom_offset[static_cast<std::size_t>(chr)]);
        const auto hi = static_cast<std::size_t>(group->chrom_offset[static_cast<std::size_t>(chr) + 1]);
        if (values.size() < hi) {
            return {};
        }
        std::vector<double> vector(values.begin() + static_cast<std::ptrdiff_t>(lo),
                                   values.begin() + static_cast<std::ptrdiff_t>(hi));
        // hic2cool fills a chromosome without a vector with NaN.
        if (std::all_of(vector.begin(), vector.end(), [](double v) { return std::isnan(v); })) {
            return {};
        }
        return vector;
    };

    // Square coolers: .hic stores one triangle.
    const std::string triangle = options.triangle;
    if (triangle != "auto" && triangle != "upper" && triangle != "lower") {
        throw ExitError("!!! ERROR. triangle must be auto, upper or lower, not " + triangle);
    }
    const bool any_square = std::any_of(chosen.begin(), chosen.end(), [](const Resolution& r) { return r.square; });
    if (triangle == "lower" && std::any_of(chosen.begin(), chosen.end(), [](const Resolution& r) { return !r.square; })) {
        throw ExitError("!!! ERROR. --triangle lower needs square coolers; " + path +
                        " stores the upper triangle only");
    }
    CoolSource source(pool, file, chosen, triangle == "lower");
    if (any_square && triangle == "auto") {
        for (const auto& r : chosen) {
            if (r.square && !source.symmetric(r)) {
                throw ExitError("!!! ERROR. The square cooler " + r.path + " in " + path +
                                " is not symmetric, and .hic files store one triangle. Choose it with "
                                "--triangle upper or --triangle lower.");
            }
        }
    }
    try {
        hicfilecpp::writeHicFile(outfile, write, source);
    } catch (const hicfilecpp::HicError& e) {
        throw ExitError(std::string("!!! ERROR. ") + e.what());
    }
    if (!options.silent) {
        if (source.skipped() > 0 && options.show_warnings) {
            console.err("!!! WARNING. " + std::to_string(source.skipped()) + " pixels with a count that is not finite were left out.");
        }
        console.out("### Finished! Output written to: " + outfile);
    }
    return outfile;
}

}  // namespace hic2cool
