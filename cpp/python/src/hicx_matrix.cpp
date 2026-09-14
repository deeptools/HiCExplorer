// hicx_matrix: one region of a Hi-C matrix as a numpy array.
//
//   f = hicx_matrix.open(path)   # .cool, .mcool, "x.mcool::/resolutions/N", .h5, .hic
//   f.fetch(region1, region2=None, resolution=None, normalization="none", transform="none")
//
// Every read touches only the rows the region needs: cool and mcool through
// coolercpp's matrix selector (the bin1_offset index and a range query over
// the pixel table), h5 through hyperslabs of the CSR rows (h5_rows.hpp), and
// .hic through hicfilecpp's block index. The values follow the reference
// readers: cooler's Cooler.matrix().fetch for cool and mcool, hicmatrix's
// hiCMatrix for h5 (lower triangle filled as on load), and hicstraw 1.3.1's
// getRecordsAsMatrix for .hic.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <coolercpp/api.hpp>
#include <coolercpp/errors.hpp>
#include <coolercpp/rangequery.hpp>
#include <coolercpp/region.hpp>
#include <hicfilecpp/errors.hpp>
#include <hicfilecpp/hicfile.hpp>

#include "h5_rows.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/hdf5_util.hpp"

namespace py = pybind11;

namespace hicx::python {

namespace {

// HDF5 is not built thread safe, and the GIL is released while reading, so
// all reads are serialised here.
std::mutex& io_mutex() {
    static std::mutex mutex;
    return mutex;
}

struct Dense {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<double> values;
};

enum class Transform { None, Log1p, Log, ObsExp };

struct Region {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

std::string join(const std::vector<std::string>& items, std::size_t limit = 25) {
    std::string out;
    for (std::size_t i = 0; i < items.size() && i < limit; ++i) {
        out += (i == 0 ? "" : ", ") + items[i];
    }
    if (items.size() > limit) {
        out += ", ... (" + std::to_string(items.size()) + " in total)";
    }
    return out;
}

// cooler's region grammar and bounds for every format: "chr1",
// "chr1:1,000,000-3,000,000", "chr1:1.5M-2M".
Region parse_region(const std::string& text, const coolercpp::ChromSizes& sizes) {
    coolercpp::GenomicRange range;
    try {
        range = coolercpp::parse_region_string(text);
    } catch (const coolercpp::Error& e) {
        throw py::value_error("malformed region '" + text + "': " + e.what());
    }
    if (!sizes.find(range.chrom).has_value()) {
        throw py::key_error("unknown chromosome '" + range.chrom + "' in region '" + text +
                            "'; the file has: " + join(sizes.names()));
    }
    try {
        const coolercpp::RegionTuple tuple = coolercpp::parse_region(coolercpp::Region(text), &sizes);
        return Region{tuple.chrom, tuple.start, tuple.end};
    } catch (const coolercpp::Error& e) {
        throw py::value_error("invalid region '" + text + "': " + e.what());
    }
}

class Backend {
  public:
    virtual ~Backend() = default;
    [[nodiscard]] virtual std::string format() const = 0;
    [[nodiscard]] virtual std::vector<std::pair<std::string, std::int64_t>> chromosomes() const = 0;
    [[nodiscard]] virtual std::vector<std::int64_t> resolutions() const = 0;
    [[nodiscard]] virtual std::vector<std::string> normalizations() const = 0;
    // Checks the arguments that need no reading, before the GIL is released.
    virtual void check(const std::optional<std::int64_t>& resolution, const std::string& normalization,
                       Transform transform) const = 0;
    [[nodiscard]] virtual Dense fetch(const std::string& region1, const std::string& region2,
                                      const std::optional<std::int64_t>& resolution,
                                      const std::string& normalization, Transform transform) = 0;
};

void check_normalization(const Backend& backend, const std::string& normalization) {
    const std::vector<std::string> names = backend.normalizations();
    if (std::find(names.begin(), names.end(), normalization) == names.end()) {
        throw py::value_error("unknown normalization '" + normalization + "' for this " +
                              backend.format() + " file; available: " + join(names));
    }
}

void check_resolution_listed(const Backend& backend, std::int64_t resolution) {
    const std::vector<std::int64_t> list = backend.resolutions();
    if (std::find(list.begin(), list.end(), resolution) == list.end()) {
        std::vector<std::string> names;
        for (const std::int64_t r : list) {
            names.push_back(std::to_string(r));
        }
        throw py::value_error("resolution " + std::to_string(resolution) + " is not in this " +
                              backend.format() + " file; available: " + join(names));
    }
}

[[noreturn]] void refuse_obs_exp(const std::string& format) {
    throw py::value_error(
        "transform 'obs_exp' is only available for .hic files, which store the expected values; "
        "for a " + format + " file observed/expected needs the distance-decay of the whole "
        "chromosome, which a region fetch does not read");
}

// ------------------------------------------------------------------ cool, mcool

class CoolBackend final : public Backend {
  public:
    // uri: a plain path, or "file::/group" for one cooler of a container.
    CoolBackend(const std::string& uri, bool multi) : multi_(multi) {
        if (!multi) {
            add(uri, std::nullopt);
        } else {
            for (const std::string& group : coolercpp::list_coolers(uri)) {
                const std::string prefix = "/resolutions/";
                if (group.rfind(prefix, 0) != 0) {
                    continue;
                }
                const std::string number = group.substr(prefix.size());
                if (number.empty() || number.find_first_not_of("0123456789") != std::string::npos) {
                    continue;
                }
                add(uri + "::" + group, std::stoll(number));
            }
            if (entries_.empty()) {
                throw py::value_error(uri + " holds coolers, but none under /resolutions/<bin size>");
            }
        }
        // The finest resolution describes the chromosomes; zoomify can round
        // the lengths of the coarser ones.
        const coolercpp::ChromSizes& sizes = entries_.begin()->second.cooler->chromsizes();
        for (std::size_t i = 0; i < sizes.size(); ++i) {
            chromosomes_.emplace_back(sizes.names()[i], sizes.lengths()[i]);
        }
        bool weight = false;
        for (const auto& [key, entry] : entries_) {
            for (const std::string& column : entry.columns) {
                if (column == "weight") {
                    weight = true;
                } else if (std::find(extra_.begin(), extra_.end(), column) == extra_.end()) {
                    extra_.push_back(column);
                }
            }
        }
        normalizations_.emplace_back("none");
        if (weight) {
            normalizations_.emplace_back("balanced");
        }
        normalizations_.insert(normalizations_.end(), extra_.begin(), extra_.end());
    }

    [[nodiscard]] std::string format() const override { return multi_ ? "mcool" : "cool"; }
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>> chromosomes() const override {
        return chromosomes_;
    }
    [[nodiscard]] std::vector<std::int64_t> resolutions() const override {
        std::vector<std::int64_t> out;
        for (const auto& [key, entry] : entries_) {
            if (entry.binsize.has_value()) {
                out.push_back(*entry.binsize);
            }
        }
        return out;
    }
    [[nodiscard]] std::vector<std::string> normalizations() const override { return normalizations_; }

    void check(const std::optional<std::int64_t>& resolution, const std::string& normalization,
               Transform transform) const override {
        if (transform == Transform::ObsExp) {
            refuse_obs_exp(format());
        }
        check_normalization(*this, normalization);
        if (multi_ && !resolution.has_value()) {
            throw py::value_error("an mcool file needs a resolution; available: " +
                                  join(resolution_names()));
        }
        if (resolution.has_value()) {
            check_resolution_listed(*this, *resolution);
        }
    }

    Dense fetch(const std::string& region1, const std::string& region2,
                const std::optional<std::int64_t>& resolution, const std::string& normalization,
                Transform) override {
        const Entry& entry = multi_ ? entries_.at(*resolution) : entries_.begin()->second;
        const coolercpp::Cooler& cooler = *entry.cooler;
        parse_region(region1, cooler.chromsizes());
        parse_region(region2, cooler.chromsizes());

        coolercpp::MatrixOptions options;
        if (normalization == "none") {
            options.balance = coolercpp::Balance(false);
        } else {
            const std::string column = normalization == "balanced" ? "weight" : normalization;
            if (std::find(entry.columns.begin(), entry.columns.end(), column) == entry.columns.end()) {
                throw py::value_error("normalization '" + normalization + "' has no bin column '" +
                                      column + "' at resolution " +
                                      std::to_string(entry.binsize.value_or(0)) + "; that resolution has: " +
                                      join(entry.columns));
            }
            options.balance = normalization == "balanced" ? coolercpp::Balance(true)
                                                          : coolercpp::Balance(column);
        }
        coolercpp::MatrixResult result =
            cooler.matrix(options).fetch(coolercpp::Region(region1), coolercpp::Region(region2));
        const coolercpp::DenseMatrix& dense = result.dense();
        Dense out;
        out.rows = dense.rows;
        out.cols = dense.cols;
        out.values = dense.values.as<double>();
        return out;
    }

  private:
    struct Entry {
        std::shared_ptr<coolercpp::Cooler> cooler;
        std::optional<std::int64_t> binsize;
        // The numeric bin columns besides chrom, start and end.
        std::vector<std::string> columns;
    };

    void add(const std::string& uri, std::optional<std::int64_t> key) {
        Entry entry;
        entry.cooler = std::make_shared<coolercpp::Cooler>(uri);
        entry.binsize = entry.cooler->binsize();
        for (const auto& [name, dtype] : entry.cooler->bins().dtypes()) {
            const bool numeric = dtype != coolercpp::DType::String &&
                                 dtype != coolercpp::DType::Categorical &&
                                 dtype != coolercpp::DType::Bool;
            if (numeric && name != "chrom" && name != "start" && name != "end") {
                entry.columns.push_back(name);
            }
        }
        entries_.emplace(key.value_or(0), std::move(entry));
    }

    [[nodiscard]] std::vector<std::string> resolution_names() const {
        std::vector<std::string> out;
        for (const std::int64_t r : resolutions()) {
            out.push_back(std::to_string(r));
        }
        return out;
    }

    bool multi_;
    std::map<std::int64_t, Entry> entries_;
    std::vector<std::pair<std::string, std::int64_t>> chromosomes_;
    std::vector<std::string> extra_;
    std::vector<std::string> normalizations_;
};

// ------------------------------------------------------------------ h5

class H5Backend final : public Backend {
  public:
    explicit H5Backend(const std::string& path) : rows_(path) {
        const h5::File file(path);
        const std::vector<std::string> chroms = file.read_strings("/intervals/chr_list");
        starts_ = file.read_int64("/intervals/start_list");
        ends_ = file.read_int64("/intervals/end_list");
        if (chroms.size() != starts_.size() || chroms.size() != ends_.size() ||
            static_cast<std::int64_t>(chroms.size()) != rows_.nbins()) {
            throw py::value_error(path + ": the bin intervals do not match the matrix size");
        }
        std::vector<std::string> names;
        std::vector<std::int64_t> lengths;
        for (std::size_t i = 0; i < chroms.size(); ++i) {
            if (i == 0 || chroms[i] != chroms[i - 1]) {
                if (std::find(names.begin(), names.end(), chroms[i]) != names.end()) {
                    throw py::value_error(path + ": the bins of chromosome '" + chroms[i] +
                                          "' are not contiguous");
                }
                names.push_back(chroms[i]);
                lengths.push_back(0);
                first_bin_.push_back(static_cast<std::int64_t>(i));
            }
            lengths.back() = std::max(lengths.back(), ends_[i]);
        }
        for (std::size_t c = 0; c < names.size(); ++c) {
            chromosomes_.emplace_back(names[c], lengths[c]);
            chrom_index_[names[c]] = c;
        }
        first_bin_.push_back(static_cast<std::int64_t>(chroms.size()));
        sizes_ = coolercpp::ChromSizes(names, lengths);

        // The bin size is the width every bin but the last of its chromosome
        // shares, as cooler.util.get_binsize decides; none for fragment bins.
        std::optional<std::int64_t> width;
        bool uniform = true;
        for (std::size_t c = 0; c + 1 < first_bin_.size(); ++c) {
            for (std::int64_t i = first_bin_[c]; i + 1 < first_bin_[c + 1]; ++i) {
                const std::int64_t w = ends_[static_cast<std::size_t>(i)] - starts_[static_cast<std::size_t>(i)];
                if (!width.has_value()) {
                    width = w;
                } else if (*width != w) {
                    uniform = false;
                }
            }
        }
        if (uniform && width.has_value()) {
            binsize_ = width;
        }
    }

    [[nodiscard]] std::string format() const override { return "h5"; }
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>> chromosomes() const override {
        return chromosomes_;
    }
    [[nodiscard]] std::vector<std::int64_t> resolutions() const override {
        return binsize_.has_value() ? std::vector<std::int64_t>{*binsize_} : std::vector<std::int64_t>{};
    }
    [[nodiscard]] std::vector<std::string> normalizations() const override { return {"none"}; }

    void check(const std::optional<std::int64_t>& resolution, const std::string& normalization,
               Transform transform) const override {
        if (transform == Transform::ObsExp) {
            refuse_obs_exp(format());
        }
        if (normalization != "none") {
            throw py::value_error("normalization '" + normalization +
                                  "' is not available for an h5 file, which stores its values as "
                                  "they are; available: none");
        }
        if (resolution.has_value()) {
            check_resolution_listed(*this, *resolution);
        }
    }

    Dense fetch(const std::string& region1, const std::string& region2,
                const std::optional<std::int64_t>&, const std::string&, Transform) override {
        const auto [a, b] = bin_range(parse_region(region1, sizes_));
        const auto [c, d] = bin_range(parse_region(region2, sizes_));
        Dense out;
        out.rows = b - a;
        out.cols = d - c;
        out.values.assign(static_cast<std::size_t>(out.rows * out.cols), 0.0);
        const auto at = [&](std::int64_t i, std::int64_t j) -> double& {
            return out.values[static_cast<std::size_t>(i * out.cols + j)];
        };
        // The stored entries S[i, j] of the region's rows.
        rows_.for_each_entry(a, b, [&](std::int64_t row, std::int64_t col, double value) {
            if (col >= c && col < d) {
                at(row - a, col - c) += value;
            }
        });
        // hicmatrix's matrix + triu(matrix, 1).T adds S[j, i] at (i, j) for
        // j < i, which only exists when some column row lies before some row.
        if (out.rows > 0 && out.cols > 0 && c < b - 1 && rows_.fills_lower_triangle()) {
            rows_.for_each_entry(c, d, [&](std::int64_t row, std::int64_t col, double value) {
                if (col >= a && col < b && col > row) {
                    at(col - a, row - c) += value;
                }
            });
        }
        return out;
    }

  private:
    // The bins overlapping [start, end) of a chromosome, as hicmatrix's
    // getRegionBinRange(chrom, start, end - 1) finds them.
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> bin_range(const Region& region) const {
        const std::size_t c = chrom_index_.at(region.chrom);
        const auto first = ends_.begin() + first_bin_[c];
        const auto last = ends_.begin() + first_bin_[c + 1];
        const std::int64_t lo = std::upper_bound(first, last, region.start) - ends_.begin();
        const auto sfirst = starts_.begin() + first_bin_[c];
        const auto slast = starts_.begin() + first_bin_[c + 1];
        const std::int64_t hi = std::lower_bound(sfirst, slast, region.end) - starts_.begin();
        return {lo, std::max(lo, hi)};
    }

    H5Rows rows_;
    std::vector<std::int64_t> starts_;
    std::vector<std::int64_t> ends_;
    std::vector<std::int64_t> first_bin_;
    std::map<std::string, std::size_t> chrom_index_;
    std::vector<std::pair<std::string, std::int64_t>> chromosomes_;
    coolercpp::ChromSizes sizes_;
    std::optional<std::int64_t> binsize_;
};

// ------------------------------------------------------------------ .hic

class HicBackend final : public Backend {
  public:
    explicit HicBackend(const std::string& path) : hic_(path) {
        std::vector<std::string> names;
        std::vector<std::int64_t> lengths;
        for (const hicfilecpp::Chromosome& chrom : hic_.getChromosomes()) {
            std::string upper = chrom.name;
            std::transform(upper.begin(), upper.end(), upper.begin(),
                           [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
            if (upper == "ALL") {
                continue;
            }
            names.push_back(chrom.name);
            lengths.push_back(chrom.length);
            index_[chrom.name] = chrom.index;
            chromosomes_.emplace_back(chrom.name, chrom.length);
        }
        sizes_ = coolercpp::ChromSizes(names, lengths);
        for (const std::int32_t r : hic_.getResolutions()) {
            resolutions_.push_back(r);
        }
        std::sort(resolutions_.begin(), resolutions_.end());
        normalizations_.emplace_back("none");
        for (const std::string& norm : hic_.getNormalizationTypes()) {
            if (norm != "NONE") {
                normalizations_.push_back(norm);
            }
        }
    }

    [[nodiscard]] std::string format() const override { return "hic"; }
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>> chromosomes() const override {
        return chromosomes_;
    }
    [[nodiscard]] std::vector<std::int64_t> resolutions() const override { return resolutions_; }
    [[nodiscard]] std::vector<std::string> normalizations() const override { return normalizations_; }

    void check(const std::optional<std::int64_t>& resolution, const std::string& normalization,
               Transform) const override {
        check_normalization(*this, normalization);
        if (!resolution.has_value()) {
            std::vector<std::string> names;
            for (const std::int64_t r : resolutions_) {
                names.push_back(std::to_string(r));
            }
            throw py::value_error("a .hic file needs a resolution; available: " + join(names));
        }
        check_resolution_listed(*this, *resolution);
    }

    Dense fetch(const std::string& region1, const std::string& region2,
                const std::optional<std::int64_t>& resolution, const std::string& normalization,
                Transform transform) override {
        const std::int64_t res = *resolution;
        const Region r1 = parse_region(region1, sizes_);
        const Region r2 = parse_region(region2, sizes_);
        // Bins overlapping [start, end), as cooler counts them.
        const auto bins = [res](const Region& r) {
            const std::int64_t lo = r.start / res;
            const std::int64_t hi = std::max(lo, (r.end + res - 1) / res);
            return std::pair<std::int64_t, std::int64_t>{lo, hi};
        };
        const auto [a, b] = bins(r1);
        const auto [c, d] = bins(r2);
        Dense out;
        out.rows = b - a;
        out.cols = d - c;
        out.values.assign(static_cast<std::size_t>(out.rows * out.cols), 0.0);
        if (out.rows == 0 || out.cols == 0) {
            return out;
        }
        const std::string norm = normalization == "none" ? "NONE" : normalization;
        const std::string type = transform == Transform::ObsExp ? "oe" : "observed";
        hicfilecpp::MatrixZoomData zoom = [&]() {
            try {
                return hic_.getMatrixZoomData(r1.chrom, r2.chrom, type, norm, "BP",
                                              static_cast<std::int32_t>(res));
            } catch (const hicfilecpp::HicError& e) {
                throw py::value_error(e.what());
            }
        }();
        if (!zoom.found()) {
            // A chromosome pair without a matrix has no contacts; missing
            // expected values make an "oe" query impossible.
            if (zoom.message().rfind("File doesn't have the given chr_chr map", 0) == 0) {
                return out;
            }
            throw py::value_error(zoom.message());
        }
        // hicstraw's x axis is the chromosome with the lower index.
        const bool swap = index_.at(r1.chrom) > index_.at(r2.chrom);
        const std::int64_t x0 = swap ? c : a;
        const std::int64_t x1 = swap ? d : b;
        const std::int64_t y0 = swap ? a : c;
        const std::int64_t y1 = swap ? b : d;
        // Aligned to bins, getRecordsAsMatrix has x1 - x0 rows and y1 - y0
        // columns; it returns a single 0 when no record is inside.
        const hicfilecpp::FloatMatrix matrix =
            zoom.getRecordsAsMatrix(x0 * res, x1 * res - 1, y0 * res, y1 * res - 1);
        if (matrix.rows != x1 - x0 || matrix.cols != y1 - y0) {
            if (matrix.rows == 1 && matrix.cols == 1) {
                return out;
            }
            throw py::value_error("unexpected matrix shape from the .hic reader");
        }
        for (std::int64_t i = 0; i < out.rows; ++i) {
            for (std::int64_t j = 0; j < out.cols; ++j) {
                out.values[static_cast<std::size_t>(i * out.cols + j)] =
                    static_cast<double>(swap ? matrix.at(j, i) : matrix.at(i, j));
            }
        }
        return out;
    }

  private:
    hicfilecpp::HiCFile hic_;
    std::map<std::string, std::int32_t> index_;
    std::vector<std::pair<std::string, std::int64_t>> chromosomes_;
    coolercpp::ChromSizes sizes_;
    std::vector<std::int64_t> resolutions_;
    std::vector<std::string> normalizations_;
};

// ------------------------------------------------------------------ module

bool starts_with_hic_magic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char magic[4] = {0, 0, 0, 0};
    in.read(magic, 4);
    return in.gcount() == 4 && magic[0] == 'H' && magic[1] == 'I' && magic[2] == 'C' && magic[3] == '\0';
}

class File {
  public:
    explicit File(std::string path) : path_(std::move(path)) {
        const std::string::size_type separator = path_.find("::");
        const std::string file = path_.substr(0, separator);
        if (!std::ifstream(file).good()) {
            PyErr_SetString(PyExc_FileNotFoundError, ("no such file: " + file).c_str());
            throw py::error_already_set();
        }
        const std::lock_guard<std::mutex> lock(io_mutex());
        if (separator != std::string::npos) {
            backend_ = std::make_unique<CoolBackend>(path_, false);
        } else if (starts_with_hic_magic(path_)) {
            backend_ = std::make_unique<HicBackend>(path_);
        } else if (!h5::is_hdf5(path_)) {
            throw py::value_error(path_ + " is neither a cool, mcool, h5 nor .hic file");
        } else if (coolercpp::is_cooler(path_)) {
            backend_ = std::make_unique<CoolBackend>(path_, false);
        } else if (!coolercpp::list_coolers(path_).empty()) {
            backend_ = std::make_unique<CoolBackend>(path_, true);
        } else if (is_hicexplorer_h5(path_)) {
            backend_ = std::make_unique<H5Backend>(path_);
        } else {
            throw py::value_error(path_ + " is an HDF5 file but neither a cooler nor a HiCExplorer h5 matrix");
        }
    }

    [[nodiscard]] const std::string& path() const noexcept { return path_; }
    [[nodiscard]] Backend& backend() const noexcept { return *backend_; }

  private:
    std::string path_;
    std::unique_ptr<Backend> backend_;
};

Transform parse_transform(const std::string& name) {
    if (name == "none") {
        return Transform::None;
    }
    if (name == "log1p") {
        return Transform::Log1p;
    }
    if (name == "log") {
        return Transform::Log;
    }
    if (name == "obs_exp") {
        return Transform::ObsExp;
    }
    throw py::value_error("unknown transform '" + name + "'; available: none, log1p, log, obs_exp");
}

py::array_t<double> to_numpy(Dense&& dense) {
    const std::vector<py::ssize_t> shape{static_cast<py::ssize_t>(dense.rows),
                                         static_cast<py::ssize_t>(dense.cols)};
    if (dense.values.empty()) {
        return py::array_t<double>(shape);
    }
    auto* owned = new std::vector<double>(std::move(dense.values));
    const py::capsule base(owned, [](void* p) { delete static_cast<std::vector<double>*>(p); });
    return py::array_t<double>(shape, owned->data(), base);
}

py::array_t<double> fetch(File& file, const std::string& region1, const std::optional<std::string>& region2,
                          const std::optional<std::int64_t>& resolution, const std::string& normalization,
                          const std::string& transform_name) {
    const Transform transform = parse_transform(transform_name);
    Backend& backend = file.backend();
    backend.check(resolution, normalization, transform);
    Dense dense;
    {
        const py::gil_scoped_release release;
        const std::lock_guard<std::mutex> lock(io_mutex());
        dense = backend.fetch(region1, region2.value_or(region1), resolution, normalization, transform);
        if (transform == Transform::Log1p) {
            std::transform(dense.values.begin(), dense.values.end(), dense.values.begin(),
                           [](double v) { return std::log1p(v); });
        } else if (transform == Transform::Log) {
            std::transform(dense.values.begin(), dense.values.end(), dense.values.begin(),
                           [](double v) { return std::log(v); });
        }
    }
    return to_numpy(std::move(dense));
}

}  // namespace

}  // namespace hicx::python

PYBIND11_MODULE(hicx_matrix, m) {
    using hicx::python::File;
    m.doc() = "Region reads of cool, mcool, HiCExplorer h5 and .hic matrices (HiCExplorer v4).";

    py::register_exception_translator([](std::exception_ptr pointer) {
        try {
            if (pointer) {
                std::rethrow_exception(pointer);
            }
        } catch (const coolercpp::KeyError& e) {
            PyErr_SetString(PyExc_KeyError, e.what());
        } catch (const coolercpp::OSError& e) {
            PyErr_SetString(PyExc_OSError, e.what());
        } catch (const coolercpp::Error& e) {
            // ValueError, IndexError, TypeError and the rest: the arguments
            // do not describe something the file holds.
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const hicfilecpp::HicError& e) {
            PyErr_SetString(PyExc_ValueError, e.what());
        } catch (const hicx::h5::Error& e) {
            PyErr_SetString(PyExc_OSError, e.what());
        }
    });

    py::class_<File>(m, "File")
        .def_property_readonly("path", &File::path)
        .def_property_readonly("format", [](const File& f) { return f.backend().format(); },
                               "\"cool\", \"mcool\", \"h5\" or \"hic\".")
        .def("chromosomes", [](const File& f) { return f.backend().chromosomes(); },
             "[(name, length_bp), ...] in file order.")
        .def("resolutions", [](const File& f) { return f.backend().resolutions(); },
             "The bin sizes, ascending. A cool or h5 file has its single bin size (none for "
             "variable bins).")
        .def("normalizations", [](const File& f) { return f.backend().normalizations(); },
             "Names accepted as normalization=, starting with \"none\".")
        .def("fetch", &hicx::python::fetch, py::arg("region1"), py::arg("region2") = py::none(),
             py::arg("resolution") = py::none(), py::arg("normalization") = "none",
             py::arg("transform") = "none",
             "The region1 x region2 block as a float64 array: rows are the bins overlapping "
             "region1, columns the bins overlapping region2.\n\n"
             "normalization: \"none\" (raw), \"balanced\" (cool/mcool weight column), a bin column "
             "of a cool/mcool file, or a .hic normalization such as \"KR\".\n"
             "transform: \"none\", \"log1p\", \"log\" or \"obs_exp\" (.hic only).")
        .def("__repr__", [](const File& f) {
            return "<hicx_matrix.File " + f.backend().format() + " '" + f.path() + "'>";
        });

    m.def(
        "open",
        [](const py::object& path) {
            return std::make_unique<File>(py::str(py::module_::import("os").attr("fspath")(path)));
        },
        py::arg("path"),
        "Opens a .cool, .mcool, \"file.mcool::/resolutions/N\", HiCExplorer .h5 or .hic file.");
}
