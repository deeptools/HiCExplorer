#include "hicx/cool_file.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "hicx/numpy_compat.hpp"

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

CoolLoadResult read_cool(const std::string& uri, const CoolLoadOptions& options) {
    const CoolFile cool(uri);
    CoolLoadResult result;
    result.data.matrix = cool.read_matrix();
    result.data.cut_intervals = cool.read_bins();
    for (const auto& [key, value] : cool.info()) {
        result.metadata.emplace(key, value.to_python_string());
    }

    std::optional<std::vector<double>> weights;
    if (options.apply_correction && cool.has_column(options.correction_factor_table)) {
        weights = cool.read_column(options.correction_factor_table);
    }

    CsrMatrix& matrix = result.data.matrix;
    if (weights.has_value()) {
        matrix.eliminate_zeros();
        if (matrix.nnz() > 1) {
            const std::vector<double>& factors = *weights;
            const bool all_nan = std::all_of(factors.begin(), factors.end(),
                                             [](double v) { return std::isnan(v); });
            if (!all_nan) {
                // 'weight' is multiplicative, the hic2cool tables KR, VC and
                // SQRT_VC are divisive.
                const bool divide = options.correction_factor_table == "KR" ||
                                    options.correction_factor_table == "VC" ||
                                    options.correction_factor_table == "SQRT_VC";
                result.correction_operator = divide ? '/' : '*';
                const json::Value* generated = cool.info_value("generated-by");
                if (generated != nullptr && generated->is_string()) {
                    const std::string& text = generated->as_string();
                    const std::size_t dash = text.find('-');
                    if (text.find("hic2cool") != std::string::npos &&
                        dash != std::string::npos) {
                        result.hic2cool_version = text.substr(dash + 1);
                    } else if (text.find("hicmatrix") != std::string::npos &&
                               dash != std::string::npos) {
                        result.hicmatrix_version = text.substr(dash + 1);
                    }
                }
                std::vector<double>& values = matrix.mutable_data();
                const std::vector<std::int64_t>& indptr = matrix.indptr();
                const std::vector<std::int32_t>& indices = matrix.indices();
                for (std::int64_t row = 0; row < matrix.rows(); ++row) {
                    const std::size_t begin =
                        static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                    const std::size_t end =
                        static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
                    for (std::size_t k = begin; k < end; ++k) {
                        const double factor =
                            factors[static_cast<std::size_t>(row)] *
                            factors[static_cast<std::size_t>(indices[k])];
                        if (divide) {
                            values[k] /= factor;
                        } else {
                            values[k] *= factor;
                        }
                    }
                }
                matrix.set_dtype("float64");
            }
            result.data.correction_factors = weights;
        }
    }

    // Replace the NaN values introduced by the weights and derive the bins
    // that hold no interaction at all (cool.py:239-255).
    for (double& value : matrix.mutable_data()) {
        if (std::isnan(value)) {
            value = 0.0;
        }
    }
    matrix.eliminate_zeros();

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
            result.data.nan_bins.push_back(bin);
        }
    }
    return result;
}

// --------------------------------------------------------------------------
// Writing

namespace {

// One block of the streamed pixel table. At 262,144 pixels a block costs
// 4 MB across the three columns, which is what keeps the writer's own
// footprint independent of the matrix size.
constexpr std::size_t kPixelBlock = 262144;

// cooler splits a pixel table longer than this into 10,000 parts before
// handing it to create_cooler (hicmatrix/lib/cool.py:366-368), and the 'sum'
// attribute is the sequential total of the per part sums. The boundaries are
// therefore part of the number that ends up in the file.
constexpr std::size_t kSplitThreshold = 10000000;
constexpr std::size_t kSplitFactor = 10000;

// numpy's ufunc reduction buffer. A per part sum is a pairwise sum over
// windows of this size accumulated in order, so a bounded buffer reproduces it
// without ever holding the column (see cpp/PLAN.md 5.2).
constexpr std::size_t kReduceBuffer = 8192;

// The total that lands in the cool 'sum' attribute.
//
// It is numpy's sum of the count column, which the Python computes over the
// values as they sit in the data frame, that is *before* the cast to the
// column dtype. Three things fix the accumulation order and precision:
//
//  * the split into parts above,
//  * numpy's own pairwise blocking inside a part, in the dtype of the column,
//    so a float32 matrix has its parts summed in single precision,
//  * the running total, which is float64 whatever the column dtype is. The
//    Python starts it at the integer literal 0, and `0 + np.float32(x)`
//    promotes to float64 under numpy 1.26's value based casting.
//
// All three are reproduced here on a bounded buffer.
class CountSum {
  public:
    CountSum(std::size_t nnz, bool single_precision)
        : single_precision_(single_precision) {
        if (nnz > kSplitThreshold) {
            // np.array_split(df, n): the first nnz % n parts are one longer.
            const std::size_t base = nnz / kSplitFactor;
            const std::size_t remainder = nnz % kSplitFactor;
            parts_.reserve(kSplitFactor);
            for (std::size_t i = 0; i < kSplitFactor; ++i) {
                parts_.push_back(base + (i < remainder ? 1 : 0));
            }
        } else {
            parts_.push_back(nnz);
        }
        buffer_.reserve(kReduceBuffer);
    }

    void add(double value) {
        buffer_.push_back(value);
        ++in_part_;
        const std::size_t part_length = current_part_length();
        if (buffer_.size() == kReduceBuffer || in_part_ == part_length) {
            flush_window();
        }
        if (in_part_ == part_length) {
            close_part();
        }
    }

    [[nodiscard]] double total() const { return total_; }

  private:
    [[nodiscard]] std::size_t current_part_length() const {
        return part_ < parts_.size() ? parts_[part_] : 0;
    }

    void flush_window() {
        if (buffer_.empty()) {
            return;
        }
        if (single_precision_) {
            std::vector<float> narrowed(buffer_.begin(), buffer_.end());
            part_total_float_ += npy::pairwise_sum(narrowed.data(), narrowed.size());
        } else {
            part_total_ += npy::pairwise_sum(buffer_.data(), buffer_.size());
        }
        buffer_.clear();
    }

    void close_part() {
        total_ += single_precision_ ? static_cast<double>(part_total_float_)
                                    : part_total_;
        part_total_ = 0.0;
        part_total_float_ = 0.0F;
        in_part_ = 0;
        ++part_;
        while (part_ < parts_.size() && parts_[part_] == 0) {
            ++part_;
        }
    }

    bool single_precision_;
    std::vector<std::size_t> parts_;
    std::size_t part_ = 0;
    std::size_t in_part_ = 0;
    std::vector<double> buffer_;
    double part_total_ = 0.0;
    float part_total_float_ = 0.0F;
    double total_ = 0.0;
};

// The chromosome table cooler derives from the bin table:
// bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
// (cooler/util.py:413). The order is the order of the *last* occurrence of
// each name, which is the order of first appearance for a well formed bin
// table and is reproduced literally for one that is not.
struct ChromTable {
    std::vector<std::string> names;
    std::vector<std::int64_t> lengths;
    std::unordered_map<std::string, std::int32_t> ids;
};

ChromTable chrom_table(const std::vector<CutInterval>& bins) {
    std::vector<std::pair<std::size_t, std::string>> last_seen;
    std::unordered_map<std::string, std::size_t> position;
    std::unordered_map<std::string, std::int64_t> length;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        const std::string& name = bins[i].chrom;
        const auto found = position.find(name);
        if (found == position.end()) {
            position.emplace(name, last_seen.size());
            last_seen.emplace_back(i, name);
        } else {
            last_seen[found->second].first = i;
        }
        length[name] = bins[i].end;
    }
    std::stable_sort(last_seen.begin(), last_seen.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });

    ChromTable table;
    table.names.reserve(last_seen.size());
    table.lengths.reserve(last_seen.size());
    for (const auto& [row, name] : last_seen) {
        (void)row;
        table.ids.emplace(name, static_cast<std::int32_t>(table.names.size()));
        table.names.push_back(name);
        table.lengths.push_back(length[name]);
    }
    return table;
}

// cooler.util.get_binsize: the single distinct width among all but the last
// bin of every chromosome, or nothing when the bins are not uniform.
std::optional<std::int64_t> uniform_bin_size(const std::vector<CutInterval>& bins) {
    // The chromosomes are visited as pandas groupby visits them, but the
    // result is a set, so only the set of widths matters.
    std::unordered_map<std::string, std::size_t> last_row;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        last_row[bins[i].chrom] = i;
    }
    std::set<std::int64_t> widths;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        if (last_row[bins[i].chrom] == i) {
            continue;  // .iloc[:-1] drops the last bin of the group
        }
        widths.insert(bins[i].end - bins[i].start);
        if (widths.size() > 1) {
            return std::nullopt;
        }
    }
    if (widths.size() == 1) {
        return *widths.begin();
    }
    return std::nullopt;
}

std::string iso_now() {
    // datetime.now().isoformat(), local time with microseconds.
    const auto now = std::chrono::system_clock::now();
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);
    const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
                            now.time_since_epoch()) %
                        std::chrono::seconds(1);
    std::tm parts{};
    localtime_r(&seconds, &parts);
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%04d-%02d-%02dT%02d:%02d:%02d.%06lld",
                  parts.tm_year + 1900, parts.tm_mon + 1, parts.tm_mday, parts.tm_hour,
                  parts.tm_min, parts.tm_sec,
                  static_cast<long long>(micros.count()));
    return buffer;
}

// Writes a numeric column produced element by element, in blocks.
template <typename T, class Produce>
void write_streamed(h5::FileWriter& file, const std::string& path, hid_t file_type,
                    hid_t mem_type, std::size_t length, h5::Filter filter,
                    Produce&& produce) {
    const h5::Handle dataset =
        file.create_dataset(path, file_type, length, length, filter);
    std::vector<T> buffer;
    buffer.reserve(std::min(length, kPixelBlock));
    std::size_t offset = 0;
    const auto flush = [&]() {
        h5::FileWriter::write_block(dataset.get(), mem_type, offset, buffer.size(),
                                    buffer.data());
        offset += buffer.size();
        buffer.clear();
    };
    produce([&](T value) {
        buffer.push_back(value);
        if (buffer.size() == kPixelBlock) {
            flush();
        }
    });
    flush();
}

}  // namespace

void write_cool(const std::string& path, MatrixData& data,
                const CoolSaveOptions& options) {
    if (path.find("::") != std::string::npos) {
        // cooler.create_cooler accepts a URI and writes into the named group,
        // which is how hicConvertFormat produces mcool and scool. That path is
        // not implemented yet and must fail loudly rather than write a plain
        // cooler to a file whose name happens to contain the separator.
        throw h5::Error("writing into a cooler group ('" + path +
                        "') is not implemented yet; mcool and scool output is "
                        "still open, see cpp/STATUS.md");
    }
    CsrMatrix& matrix = data.matrix;
    const std::int64_t nbins = matrix.rows();
    if (data.cut_intervals.size() != static_cast<std::size_t>(nbins)) {
        throw h5::Error("the bin table has " + std::to_string(data.cut_intervals.size()) +
                        " entries but the matrix has " + std::to_string(nbins) +
                        " rows");
    }

    // ---- create_cooler_input, in the order hicmatrix does it ----
    matrix.eliminate_zeros();

    if (!data.nan_bins.empty() && options.file_was_h5) {
        // cool.py:267-279. An entry is dropped when *both* of its bins are NaN
        // bins, which is what the logical_not(logical_or(...)) computes.
        std::vector<double>& values = matrix.mutable_data();
        std::vector<char> is_nan_bin(static_cast<std::size_t>(nbins), 0);
        for (const std::int64_t bin : data.nan_bins) {
            if (bin >= 0 && bin < nbins) {
                is_nan_bin[static_cast<std::size_t>(bin)] = 1;
            }
        }
        const std::vector<std::int64_t>& indptr = matrix.indptr();
        const std::vector<std::int32_t>& indices = matrix.indices();
        for (std::int64_t row = 0; row < nbins; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                if (is_nan_bin[static_cast<std::size_t>(row)] != 0 &&
                    is_nan_bin[static_cast<std::size_t>(indices[k])] != 0) {
                    values[k] = 0.0;
                }
            }
        }
        matrix.eliminate_zeros();
    }

    for (double& value : matrix.mutable_data()) {
        if (std::isnan(value)) {
            value = 0.0;
        }
    }
    matrix.eliminate_zeros();

    // ---- correction factors ----
    const bool have_factors =
        data.correction_factors.has_value() && !data.correction_factors->empty();
    bool write_weight = false;
    if (have_factors && options.apply_correction) {
        std::vector<double>& factors = *data.correction_factors;
        const bool invert = (options.hic2cool_version.has_value() &&
                             *options.hic2cool_version >= "0.5") ||
                            options.file_was_h5 ||
                            (options.correction_operator.has_value() &&
                             *options.correction_operator == '/');
        char operation = options.correction_operator.value_or('\0');
        if (invert) {
            for (double& factor : factors) {
                factor = 1.0 / factor;
                if (std::isnan(factor) || std::isinf(factor)) {
                    factor = 0.0;
                }
            }
            operation = '*';
        }
        write_weight = true;
        // Revert the correction so that the stored counts are raw again
        // (cool.py:322-341). Only the entries that will be written are
        // touched, because the Python has replaced its matrix with the upper
        // triangle by this point.
        if (operation == '*' || operation == '\0') {
            const std::vector<std::int64_t>& indptr = matrix.indptr();
            const std::vector<std::int32_t>& indices = matrix.indices();
            std::vector<double>& counts = matrix.mutable_data();
            for (std::int64_t row = 0; row < nbins; ++row) {
                const std::size_t begin =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                const std::size_t end =
                    static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
                for (std::size_t k = begin; k < end; ++k) {
                    if (options.symmetric && indices[k] < row) {
                        continue;
                    }
                    counts[k] /= factors[static_cast<std::size_t>(row)] *
                                 factors[static_cast<std::size_t>(indices[k])];
                }
            }
        }
        matrix.set_dtype("float64");
        matrix.eliminate_zeros();
    } else if (have_factors) {
        write_weight = true;
    }

    // ---- the shape of what is about to be written ----
    const std::vector<std::int64_t> bin1_offset =
        options.symmetric ? matrix.upper_triangle_indptr()
                          : matrix.stored_indptr_without_zeros();
    const std::size_t nnz = static_cast<std::size_t>(bin1_offset.back());

    const DType kind = matrix.dtype_kind();
    // cool.py:355-364: an integer matrix keeps cooler's int32 count column, a
    // float matrix carries its own dtype into the file.
    hid_t count_file_type = H5T_STD_I32LE;
    if (!options.enforce_integer) {
        if (kind == DType::Float32) {
            count_file_type = H5T_IEEE_F32LE;
        } else if (kind == DType::Float64) {
            count_file_type = H5T_IEEE_F64LE;
        }
    }
    // The sum is taken over the values as the data frame holds them, which is
    // the matrix dtype, or float64 after np.rint.
    const bool integer_sum = kind == DType::Integer && !options.enforce_integer;
    const bool float32_sum = kind == DType::Float32 && !options.enforce_integer;

    const ChromTable chroms = chrom_table(data.cut_intervals);
    const std::optional<std::int64_t> bin_size = uniform_bin_size(data.cut_intervals);

    // ---- the file ----
    std::remove(path.c_str());
    h5::FileWriter file(path);

    file.create_group("/chroms");
    {
        std::size_t width = 0;
        for (const std::string& name : chroms.names) {
            width = std::max(width, name.size());
        }
        const h5::Handle type = h5::fixed_string_type(width);
        const std::size_t item = std::max<std::size_t>(width, 1);
        const h5::Handle dataset =
            file.create_dataset("/chroms/name", type.get(), chroms.names.size(),
                                chroms.names.size(), h5::Filter::CoolerDefault);
        std::vector<char> buffer(chroms.names.size() * item, '\0');
        for (std::size_t i = 0; i < chroms.names.size(); ++i) {
            std::copy_n(chroms.names[i].data(), chroms.names[i].size(),
                        buffer.begin() + static_cast<std::ptrdiff_t>(i * item));
        }
        h5::FileWriter::write_block(dataset.get(), type.get(), 0, chroms.names.size(),
                                    buffer.data());
    }
    write_streamed<std::int32_t>(file, "/chroms/length", H5T_STD_I32LE, H5T_NATIVE_INT32,
                                 chroms.lengths.size(), h5::Filter::CoolerDefault,
                                 [&](auto emit) {
                                     for (const std::int64_t length : chroms.lengths) {
                                         emit(static_cast<std::int32_t>(length));
                                     }
                                 });

    file.create_group("/bins");
    {
        const h5::Handle type = h5::enum_type(chroms.names);
        const h5::Handle memory_type = h5::enum_type(chroms.names, H5T_NATIVE_INT32);
        const h5::Handle dataset = file.create_dataset(
            "/bins/chrom", type.get(), data.cut_intervals.size(),
            data.cut_intervals.size(), h5::Filter::CoolerDefault);
        std::vector<std::int32_t> buffer;
        buffer.reserve(std::min(data.cut_intervals.size(), kPixelBlock));
        std::size_t offset = 0;
        const auto flush = [&]() {
            h5::FileWriter::write_block(dataset.get(), memory_type.get(), offset,
                                        buffer.size(), buffer.data());
            offset += buffer.size();
            buffer.clear();
        };
        for (const CutInterval& bin : data.cut_intervals) {
            buffer.push_back(chroms.ids.at(bin.chrom));
            if (buffer.size() == kPixelBlock) {
                flush();
            }
        }
        flush();
    }
    write_streamed<std::int32_t>(file, "/bins/start", H5T_STD_I32LE, H5T_NATIVE_INT32,
                                 data.cut_intervals.size(), h5::Filter::CoolerDefault,
                                 [&](auto emit) {
                                     for (const CutInterval& bin : data.cut_intervals) {
                                         emit(static_cast<std::int32_t>(bin.start));
                                     }
                                 });
    write_streamed<std::int32_t>(file, "/bins/end", H5T_STD_I32LE, H5T_NATIVE_INT32,
                                 data.cut_intervals.size(), h5::Filter::CoolerDefault,
                                 [&](auto emit) {
                                     for (const CutInterval& bin : data.cut_intervals) {
                                         emit(static_cast<std::int32_t>(bin.end));
                                     }
                                 });
    if (write_weight) {
        // An extra bin column goes through cooler.core.put, which uses gzip
        // without shuffle and an unlimited maximum shape. convertNansToOnes
        // turns a NaN weight into 1.0.
        const std::vector<double>& factors = *data.correction_factors;
        const h5::Handle dataset = file.create_dataset(
            "/bins/weight", H5T_IEEE_F64LE, factors.size(),
            h5::FileWriter::kUnlimited, h5::Filter::CoolerColumn);
        std::vector<double> buffer;
        buffer.reserve(std::min(factors.size(), kPixelBlock));
        std::size_t offset = 0;
        const auto flush = [&]() {
            h5::FileWriter::write_block(dataset.get(), H5T_NATIVE_DOUBLE, offset,
                                        buffer.size(), buffer.data());
            offset += buffer.size();
            buffer.clear();
        };
        for (const double factor : factors) {
            buffer.push_back(std::isnan(factor) ? 1.0 : factor);
            if (buffer.size() == kPixelBlock) {
                flush();
            }
        }
        flush();
    }

    // ---- pixels ----
    file.create_group("/pixels");
    const std::size_t max_size =
        static_cast<std::size_t>(nbins) * static_cast<std::size_t>(nbins - 1) / 2 +
        static_cast<std::size_t>(nbins);
    const std::size_t init_size =
        std::min(static_cast<std::size_t>(5 * nbins), max_size);
    const h5::Handle bin1 =
        file.create_dataset("/pixels/bin1_id", H5T_STD_I32LE, init_size, max_size,
                            h5::Filter::CoolerDefault);
    const h5::Handle bin2 =
        file.create_dataset("/pixels/bin2_id", H5T_STD_I32LE, init_size, max_size,
                            h5::Filter::CoolerDefault);
    const h5::Handle count =
        file.create_dataset("/pixels/count", count_file_type, init_size, max_size,
                            h5::Filter::CoolerDefault);
    h5::FileWriter::resize(bin1.get(), nnz);
    h5::FileWriter::resize(bin2.get(), nnz);
    h5::FileWriter::resize(count.get(), nnz);

    CountSum sum(nnz, float32_sum);
    std::int64_t integer_total = 0;
    {
        std::vector<std::int32_t> rows;
        std::vector<std::int32_t> columns;
        std::vector<double> counts;
        rows.reserve(kPixelBlock);
        columns.reserve(kPixelBlock);
        counts.reserve(kPixelBlock);
        std::size_t offset = 0;
        const auto flush = [&]() {
            if (rows.empty()) {
                return;
            }
            h5::FileWriter::write_block(bin1.get(), H5T_NATIVE_INT32, offset, rows.size(),
                                        rows.data());
            h5::FileWriter::write_block(bin2.get(), H5T_NATIVE_INT32, offset,
                                        columns.size(), columns.data());
            h5::FileWriter::write_block(count.get(), H5T_NATIVE_DOUBLE, offset,
                                        counts.size(), counts.data());
            offset += rows.size();
            rows.clear();
            columns.clear();
            counts.clear();
        };
        const auto visit = [&](std::int64_t row, std::int64_t column, double value) {
            if (options.enforce_integer) {
                // np.rint is round half to even, which is what nearbyint does
                // under the default rounding mode.
                value = std::nearbyint(value);
            }
            rows.push_back(static_cast<std::int32_t>(row));
            columns.push_back(static_cast<std::int32_t>(column));
            counts.push_back(value);
            if (integer_sum) {
                integer_total += static_cast<std::int64_t>(value);
            } else {
                sum.add(value);
            }
            if (rows.size() == kPixelBlock) {
                flush();
            }
        };
        if (options.symmetric) {
            matrix.for_each_upper(visit);
        } else {
            matrix.for_each_stored(visit);
        }
        flush();
    }

    // ---- indexes ----
    file.create_group("/indexes");
    {
        std::vector<std::int64_t> chrom_offset(chroms.names.size() + 1, 0);
        std::size_t current = 0;
        for (std::size_t i = 0; i < data.cut_intervals.size();) {
            const std::int32_t value = chroms.ids.at(data.cut_intervals[i].chrom);
            const std::size_t start = i;
            while (i < data.cut_intervals.size() &&
                   chroms.ids.at(data.cut_intervals[i].chrom) == value) {
                ++i;
            }
            for (std::size_t k = current;
                 k <= static_cast<std::size_t>(value) && k < chrom_offset.size(); ++k) {
                chrom_offset[k] = static_cast<std::int64_t>(start);
            }
            current = static_cast<std::size_t>(value) + 1;
        }
        for (std::size_t k = current; k < chrom_offset.size(); ++k) {
            chrom_offset[k] = static_cast<std::int64_t>(data.cut_intervals.size());
        }
        write_streamed<std::int64_t>(file, "/indexes/chrom_offset", H5T_STD_I64LE,
                                     H5T_NATIVE_INT64, chrom_offset.size(),
                                     h5::Filter::CoolerDefault, [&](auto emit) {
                                         for (const std::int64_t value : chrom_offset) {
                                             emit(value);
                                         }
                                     });
    }
    write_streamed<std::int64_t>(file, "/indexes/bin1_offset", H5T_STD_I64LE,
                                 H5T_NATIVE_INT64, bin1_offset.size(),
                                 h5::Filter::CoolerDefault, [&](auto emit) {
                                     for (const std::int64_t value : bin1_offset) {
                                         emit(value);
                                     }
                                 });

    // ---- attributes ----
    // The info dictionary hicmatrix builds (cool.py:375-402) in its insertion
    // order, which is also the key order of the 'metadata' JSON string.
    std::vector<std::pair<std::string, std::string>> info{
        {"format", "HDF5::Cooler"},
        {"format-url", options.format_url},
        {"generated-by", options.generated_by},
        {"generated-by-cooler-lib", options.generated_by_cooler_lib},
        {"tool-url", options.tool_url},
    };
    const auto metadata_field = [&](const std::string& key) -> const std::string* {
        const auto found = options.hic_metadata.find(key);
        return found == options.hic_metadata.end() ? nullptr : &found->second;
    };
    if (options.has_hic_metadata) {
        for (const char* key :
             {"matrix-generated-by", "matrix-generated-by-url", "genome-assembly"}) {
            if (const std::string* value = metadata_field(key)) {
                info.emplace_back(key, *value);
            }
        }
    }

    file.set_attribute("/", "bin-type", std::string(bin_size.has_value() ? "fixed"
                                                                        : "variable"));
    if (bin_size.has_value()) {
        file.set_attribute("/", "bin-size", *bin_size);
    } else {
        // cooler writes the literal string "null" here, which is what makes a
        // variable bin size read back as None.
        file.set_attribute("/", "bin-size", std::string("null"));
    }
    file.set_attribute("/", "storage-mode",
                       std::string(options.symmetric ? "symmetric-upper" : "square"));
    file.set_attribute("/", "nchroms", static_cast<std::int64_t>(chroms.names.size()));
    file.set_attribute("/", "nbins", nbins);
    if (integer_sum) {
        file.set_attribute("/", "sum", integer_total);
    } else {
        file.set_attribute("/", "sum", sum.total());
    }
    file.set_attribute("/", "nnz", static_cast<std::int64_t>(nnz));
    {
        const std::string* assembly =
            options.has_hic_metadata ? metadata_field("genome-assembly") : nullptr;
        file.set_attribute("/", "genome-assembly",
                           assembly != nullptr ? *assembly : std::string("unknown"));
    }
    file.set_attribute("/", "metadata", json::dump_object(info));
    file.set_attribute("/", "creation-date",
                       options.creation_date.empty() ? iso_now()
                                                     : options.creation_date);
    file.set_attribute("/", "format-version", static_cast<std::int64_t>(3));
    // cooler writes 'format', 'format-url' and 'generated-by' first and
    // hicmatrix then overwrites them; only the final values are stored here.
    file.set_attribute("/", "format", std::string("HDF5::Cooler"));
    file.set_attribute("/", "format-url", options.format_url);
    file.set_attribute("/", "generated-by", options.generated_by);
    file.set_attribute("/", "generated-by-cooler-lib", options.generated_by_cooler_lib);
    file.set_attribute("/", "tool-url", options.tool_url);
    for (const auto& [key, value] : info) {
        if (key == "matrix-generated-by" || key == "matrix-generated-by-url") {
            file.set_attribute("/", key, value);
        }
    }
}

}  // namespace hicx
