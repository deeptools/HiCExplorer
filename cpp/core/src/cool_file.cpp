#include "hicx/cool_file.hpp"

#include <algorithm>
#include <stdexcept>

namespace hicx {

namespace {

std::pair<std::string, std::string> split_uri(const std::string& uri) {
    const std::size_t marker = uri.find("::");
    if (marker == std::string::npos) {
        return {uri, "/"};
    }
    std::string group = uri.substr(marker + 2);
    if (group.empty()) {
        group = "/";
    }
    if (group.front() != '/') {
        group.insert(group.begin(), '/');
    }
    return {uri.substr(0, marker), group};
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

bool is_cooler(const std::string& path) {
    const auto [filename, root] = split_uri(path);
    if (!h5::is_hdf5(filename)) {
        return false;
    }
    try {
        const h5::File file(filename);
        if (!file.exists(root)) {
            return false;
        }
        const std::string prefix = (root == "/") ? "/" : root + "/";
        for (const char* group : {"chroms", "bins", "pixels", "indexes"}) {
            if (!file.exists(prefix + group)) {
                return false;
            }
        }
        return file.exists(prefix + "pixels/bin1_id") &&
               file.exists(prefix + "pixels/bin2_id") &&
               file.exists(prefix + "bins/chrom");
    } catch (const h5::Error&) {
        return false;
    }
}

bool check_cooler(const std::string& path) {
    if (ends_with(path, ".cool")) {
        return true;
    }
    if (is_cooler(path)) {
        return true;
    }
    return path.find(".mcool") != std::string::npos;
}

CoolFile::CoolFile(const std::string& uri) {
    auto [filename, root] = split_uri(uri);
    filename_ = std::move(filename);
    root_ = std::move(root);
    file_ = std::make_shared<h5::File>(filename_);
    if (!file_->exists(root_)) {
        throw h5::Error("no cooler found at " + uri);
    }

    // cooler.api.info: json decode every string attribute, keep the rest.
    for (const auto& [key, value] : file_->attributes(root_)) {
        if (std::holds_alternative<std::string>(value)) {
            const std::string& text = std::get<std::string>(value);
            std::optional<json::Value> decoded = json::parse(text);
            info_.emplace(key, decoded.has_value() ? *decoded : json::Value::string(text));
        } else if (std::holds_alternative<std::int64_t>(value)) {
            info_.emplace(key, json::Value::integer(std::get<std::int64_t>(value)));
        } else {
            info_.emplace(key, json::Value::number(std::get<double>(value)));
        }
    }

    chrom_names_ = file_->read_strings(path_of("chroms/name"));
    chrom_lengths_ = file_->read_int64(path_of("chroms/length"));
}

std::string CoolFile::path_of(const std::string& relative) const {
    return root_ == "/" ? "/" + relative : root_ + "/" + relative;
}

const json::Value* CoolFile::info_value(const std::string& key) const {
    const auto it = info_.find(key);
    return it == info_.end() ? nullptr : &it->second;
}

std::vector<std::string> CoolFile::bin_columns() const {
    std::vector<std::string> columns{"chrom", "start", "end"};
    for (const std::string& name : file_->children(path_of("bins"))) {
        if (std::find(columns.begin(), columns.end(), name) == columns.end()) {
            columns.push_back(name);
        }
    }
    return columns;
}

std::int64_t CoolFile::nbins() const {
    const json::Value* value = info_value("nbins");
    if (value != nullptr && value->is_number()) {
        return value->as_int();
    }
    return static_cast<std::int64_t>(file_->dataset_length(path_of("bins/start")));
}

std::int64_t CoolFile::nnz() const {
    const json::Value* value = info_value("nnz");
    if (value != nullptr && value->is_number()) {
        return value->as_int();
    }
    return static_cast<std::int64_t>(file_->dataset_length(path_of("pixels/bin1_id")));
}

std::vector<CutInterval> CoolFile::read_bins() const {
    const std::vector<std::int32_t> chrom_ids = file_->read_int32(path_of("bins/chrom"));
    const std::vector<std::int64_t> starts = file_->read_int64(path_of("bins/start"));
    const std::vector<std::int64_t> ends = file_->read_int64(path_of("bins/end"));
    std::vector<CutInterval> intervals;
    intervals.reserve(chrom_ids.size());
    for (std::size_t i = 0; i < chrom_ids.size(); ++i) {
        const std::size_t chrom_id = static_cast<std::size_t>(chrom_ids[i]);
        if (chrom_id >= chrom_names_.size()) {
            throw h5::Error("bin " + std::to_string(i) + " refers to an unknown chromosome");
        }
        intervals.push_back(CutInterval{chrom_names_[chrom_id], starts[i], ends[i], 1.0, ""});
    }
    return intervals;
}

bool CoolFile::has_column(const std::string& column) const {
    return file_->exists(path_of("bins/" + column));
}

std::vector<double> CoolFile::read_column(const std::string& column) const {
    return file_->read_doubles(path_of("bins/" + column));
}

CsrMatrix CoolFile::read_matrix() const {
    const std::int64_t bins = nbins();
    const std::string count_path = path_of("pixels/count");
    const std::string dtype = file_->dataset_dtype(count_path);
    const std::int64_t pixels = nnz();

    // indexes/bin1_offset is the CSR row offset array of the pixel table. Using
    // it directly means the bin1_id column never has to be read, which for the
    // 61.8 M pixel gm12878 matrix avoids staging a 247 MB array whose only
    // purpose would be to be counted and thrown away.
    const std::string offset_path = path_of("indexes/bin1_offset");
    if (file_->exists(offset_path) &&
        static_cast<std::int64_t>(file_->dataset_length(offset_path)) == bins + 1) {
        std::vector<std::int64_t> indptr = file_->read_int64(offset_path);
        if (indptr.back() == pixels) {
            std::vector<std::int32_t> bin2 = file_->read_int32(path_of("pixels/bin2_id"));
            std::vector<double> counts = file_->read_doubles(count_path);
            return CsrMatrix(bins, bins, std::move(indptr), std::move(bin2),
                             std::move(counts), dtype);
        }
    }

    // Fall back to the coordinate columns when the index is absent or stale.
    std::vector<std::int32_t> bin1 = file_->read_int32(path_of("pixels/bin1_id"));
    std::vector<std::int32_t> bin2 = file_->read_int32(path_of("pixels/bin2_id"));
    std::vector<double> counts = file_->read_doubles(count_path);

    // cool pixel tables are stored sorted by (bin1_id, bin2_id) without
    // duplicates. Detect that and build the row offsets in place, otherwise
    // fall back to the general coordinate constructor.
    bool sorted = true;
    for (std::size_t k = 1; k < bin1.size(); ++k) {
        if (bin1[k] < bin1[k - 1] || (bin1[k] == bin1[k - 1] && bin2[k] <= bin2[k - 1])) {
            sorted = false;
            break;
        }
    }
    if (!sorted) {
        return CsrMatrix::from_coo(bins, bins, bin1, bin2, std::move(counts), dtype);
    }

    std::vector<std::int64_t> indptr(static_cast<std::size_t>(bins) + 1, 0);
    for (const std::int32_t row : bin1) {
        ++indptr[static_cast<std::size_t>(row) + 1];
    }
    for (std::size_t i = 1; i < indptr.size(); ++i) {
        indptr[i] += indptr[i - 1];
    }
    bin1.clear();
    bin1.shrink_to_fit();
    return CsrMatrix(bins, bins, std::move(indptr), std::move(bin2), std::move(counts),
                     dtype);
}

}  // namespace hicx
