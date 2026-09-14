#include "hicx/cool_adapter.hpp"

#include <coolercpp/coolercpp.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <utility>

#include "hicx/hdf5_util.hpp"
#include "hicx/numpy_compat.hpp"

namespace hicx {

namespace {

// Stored pixel rows per chunk when a pixel table is streamed. A chunk is held
// twice while it is converted (about 60 bytes per pixel), so this bounds the
// transient memory of a read at roughly 12 MB.
constexpr std::int64_t kChunkPixels = 200000;
// Rows per column read of a whole-table read. HDF5 sizes its conversion
// buffers by the request, so one read of a 4 M row slice raised the peak of a
// 4.2 M pixel load by 20 MB; at this size the transient stays small while the
// destination arrays are still sized once.
constexpr std::int64_t kReadRows = kChunkPixels;

// Sorts the columns of a row that is stored out of order, as select_bins
// returns them (cpp/tests/test_cool_chromosome_cut.cpp).
void sort_row_columns(CsrMatrix::Arrays& arrays) {
    std::vector<std::pair<std::int32_t, double>> row_entries;
    for (std::int64_t row = 0; row < arrays.rows; ++row) {
        const auto begin = static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row)]);
        const auto end = static_cast<std::size_t>(arrays.indptr[static_cast<std::size_t>(row) + 1]);
        bool sorted = true;
        for (std::size_t k = begin + 1; k < end && sorted; ++k) {
            sorted = arrays.indices[k - 1] <= arrays.indices[k];
        }
        if (sorted) {
            continue;
        }
        row_entries.clear();
        for (std::size_t k = begin; k < end; ++k) {
            row_entries.emplace_back(arrays.indices[k], arrays.data[k]);
        }
        std::stable_sort(row_entries.begin(), row_entries.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        for (std::size_t k = begin; k < end; ++k) {
            arrays.indices[k] = row_entries[k - begin].first;
            arrays.data[k] = row_entries[k - begin].second;
        }
    }
}

// cooler splits a pixel table longer than this into 10,000 parts before
// handing it to create_cooler (hicmatrix/lib/cool.py:366-368), and the 'sum'
// attribute is the sequential total of the per part sums, so the boundaries
// are part of the number that ends up in the file.
constexpr std::size_t kSplitThreshold = 10000000;
constexpr std::size_t kSplitFactor = 10000;

// Pixels per chunk handed to coolercpp when a table of at most 10^7 pixels is
// streamed instead of passed as one data frame.
constexpr std::size_t kWriteBlock = 262144;

// numpy's ufunc reduction buffer.
constexpr std::size_t kReduceBuffer = 8192;

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// coolercpp raises Python-style exceptions; the tools report h5::Error.
template <typename F>
auto guarded(F&& body) -> decltype(body()) {
    try {
        return body();
    } catch (const coolercpp::Error& error) {
        throw h5::Error(error.what());
    }
}

// The info values in the form the former v4 reader produced them: json_lite
// values, with fixed length string attributes JSON decoded like the variable
// length ones.
json::Value to_json_lite(const coolercpp::json::Value& value) {
    using coolercpp::json::Type;
    switch (value.type()) {
        case Type::Null: return json::Value::null();
        case Type::Bool: return json::Value::boolean(value.as_bool());
        case Type::Int: return json::Value::integer(value.as_int());
        case Type::UInt: return json::Value::integer(static_cast<std::int64_t>(value.as_uint()));
        case Type::Double: return json::Value::number(value.as_double());
        case Type::String: return json::Value::string(value.as_string());
        case Type::Bytes: {
            std::optional<json::Value> decoded = json::parse(value.as_string());
            return decoded.has_value() ? *decoded : json::Value::string(value.as_string());
        }
        case Type::Array: {
            json::Array items;
            for (const coolercpp::json::Value& item : value.as_array()) {
                items.push_back(to_json_lite(item));
            }
            return json::Value::array(std::move(items));
        }
        case Type::Object: {
            json::Object members;
            for (const auto& [key, item] : value.as_object()) {
                members[key] = to_json_lite(item);
            }
            return json::Value::object(std::move(members));
        }
    }
    return json::Value::null();
}

}  // namespace

bool is_cooler(const std::string& path) {
    try {
        return coolercpp::is_cooler(path);
    } catch (const coolercpp::Error&) {
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

// --------------------------------------------------------------------------
// Reading

CoolFile::CoolFile(const std::string& uri) {
    guarded([&] {
        auto cooler = std::make_shared<const coolercpp::Cooler>(uri);
        filename_ = cooler->filename();
        root_ = cooler->root();
        const coolercpp::json::Value info = cooler->info();
        for (const auto& [key, value] : info.as_object()) {
            info_[key] = to_json_lite(value);
        }
        chrom_names_ = cooler->chromnames();
        chrom_lengths_ = cooler->chromsizes().lengths();
        cooler_ = std::move(cooler);
    });
}

const json::Value* CoolFile::info_value(const std::string& key) const {
    const auto it = info_.find(key);
    return it == info_.end() ? nullptr : &it->second;
}

std::vector<std::string> CoolFile::bin_columns() const {
    return guarded([&] { return cooler_->bins().columns(); });
}

std::int64_t CoolFile::nbins() const {
    const json::Value* value = info_value("nbins");
    if (value != nullptr && value->is_number()) {
        return value->as_int();
    }
    return guarded([&] { return static_cast<std::int64_t>(cooler_->bins().all().num_rows()); });
}

std::int64_t CoolFile::nnz() const {
    const json::Value* value = info_value("nnz");
    if (value != nullptr && value->is_number()) {
        return value->as_int();
    }
    return guarded(
        [&] { return static_cast<std::int64_t>(cooler_->pixels()["bin1_id"].all().num_rows()); });
}

std::pair<std::int64_t, std::int64_t> CoolFile::extent(const std::string& region) const {
    return guarded([&] { return cooler_->extent(coolercpp::Region(region)); });
}

std::vector<CutInterval> CoolFile::read_bins() const {
    return guarded([&] {
        const coolercpp::Table table =
            cooler_->bins()[coolercpp::Fields({"chrom", "start", "end"})].all();
        const coolercpp::Column& chrom = table["chrom"];
        const std::vector<std::int64_t> starts = table["start"].as<std::int64_t>();
        const std::vector<std::int64_t> ends = table["end"].as<std::int64_t>();
        std::vector<CutInterval> intervals;
        intervals.reserve(starts.size());
        for (std::size_t i = 0; i < starts.size(); ++i) {
            intervals.push_back(CutInterval{chrom.label(i), starts[i], ends[i], 1.0, ""});
        }
        return intervals;
    });
}

bool CoolFile::has_column(const std::string& column) const {
    const std::vector<std::string> columns = bin_columns();
    return std::find(columns.begin(), columns.end(), column) != columns.end();
}

std::vector<double> CoolFile::read_column(const std::string& column) const {
    return guarded([&] { return cooler_->bins()[column].all()[column].as<double>(); });
}

std::string CoolFile::count_dtype() const {
    return guarded([&] {
        for (const auto& [name, dtype] : cooler_->pixels()[coolercpp::Fields({"count"})].dtypes()) {
            if (name == "count") {
                return std::string(coolercpp::dtype_name(dtype));
            }
        }
        throw h5::Error("the pixel table of " + filename_ + " has no count column");
    });
}

void CoolFile::for_each_pixel_chunk(std::int64_t row_first, std::int64_t row_last,
                                    std::int64_t col_first, std::int64_t col_last,
                                    const std::function<void(const PixelChunk&)>& visit) const {
    guarded([&] {
        const coolercpp::RangeQuery2D query(*cooler_, coolercpp::RangeQuery2D::Kind::Direct,
                                            "count", {row_first, row_last, col_first, col_last},
                                            kChunkPixels);
        PixelChunk chunk;
        for (std::size_t i = 0; i < query.n_chunks(); ++i) {
            const coolercpp::Table table = query.get_chunk(i);
            chunk.bin1 = table["bin1_id"].as<std::int64_t>();
            chunk.bin2 = table["bin2_id"].as<std::int64_t>();
            chunk.count = table["count"].as<double>();
            visit(chunk);
        }
    });
}

std::optional<CsrMatrix> CoolFile::read_from_index(std::int64_t first, std::int64_t last,
                                                  const std::string& dtype) const {
    const std::int64_t bins = nbins();
    return guarded([&]() -> std::optional<CsrMatrix> {
        // indexes/bin1_offset is the CSR row offset array of the pixel table:
        // with it the pixels of rows [first, last) are two column reads into
        // arrays sized once, and bin1_id is never read. A missing, short or
        // inconsistent index leaves the read to the range query, which counts
        // rows from bin1_id.
        if (bins < 0 || bins > std::numeric_limits<std::int32_t>::max() || first < 0 ||
            last > bins || first >= last) {
            return std::nullopt;
        }
        std::optional<coolercpp::DatasetReader> offsets;
        try {
            offsets.emplace(*cooler_, "indexes/bin1_offset");
        } catch (const coolercpp::KeyError&) {
            return std::nullopt;
        }
        const coolercpp::DatasetReader bin2(*cooler_, "pixels/bin2_id");
        const coolercpp::DatasetReader count(*cooler_, "pixels/count");
        const std::int64_t pixels = bin2.size();
        if (offsets->size() != bins + 1 || count.size() != pixels) {
            return std::nullopt;
        }
        std::vector<std::int64_t> index(static_cast<std::size_t>(bins) + 1);
        offsets->read_into(0, bins + 1, std::span<std::int64_t>(index));
        if (index.front() != 0 || index.back() != pixels ||
            !std::is_sorted(index.begin(), index.end())) {
            return std::nullopt;
        }
        const std::int64_t size = last - first;
        const std::int64_t span_first = index[static_cast<std::size_t>(first)];
        const std::int64_t span_last = index[static_cast<std::size_t>(last)];
        CsrMatrix::Arrays arrays;
        arrays.rows = size;
        arrays.cols = size;
        arrays.dtype = dtype;
        if (first == 0 && last == bins) {
            // The whole table: every column is inside and the rows stay as
            // stored, so the columns are read straight into their final
            // arrays.
            arrays.indices.resize(static_cast<std::size_t>(pixels));
            arrays.data.resize(static_cast<std::size_t>(pixels));
            for (std::int64_t lo = 0; lo < pixels; lo += kReadRows) {
                const std::int64_t hi = std::min(pixels, lo + kReadRows);
                const auto at = static_cast<std::size_t>(lo);
                const auto n = static_cast<std::size_t>(hi - lo);
                bin2.read_into(lo, hi, std::span<std::int32_t>(arrays.indices.data() + at, n));
                count.read_into(lo, hi, std::span<double>(arrays.data.data() + at, n));
            }
            arrays.indptr = std::move(index);
            return CsrMatrix::adopt(std::move(arrays));
        }
        // A block: the row span also holds the pixels whose column lies past
        // the block (the other chromosomes of a multi-chromosome file), so it
        // is read in range query sized slices and only the columns inside
        // [first, last) are kept, renumbered. The capacity reserved for the
        // whole span is only resident where it is written, so the peak is the
        // kept pixels plus one slice. Rows stored out of order are sorted
        // afterwards.
        const auto span = static_cast<std::size_t>(span_last - span_first);
        arrays.indices.reserve(span);
        arrays.data.reserve(span);
        arrays.indptr.assign(static_cast<std::size_t>(size) + 1, 0);
        const auto slice = static_cast<std::size_t>(std::min<std::int64_t>(
            kChunkPixels, std::max<std::int64_t>(span_last - span_first, 1)));
        std::vector<std::int32_t> slice_columns(slice);
        std::vector<double> slice_counts(slice);
        std::int64_t row = first;
        for (std::int64_t lo = span_first; lo < span_last; lo += kChunkPixels) {
            const std::int64_t hi = std::min(span_last, lo + kChunkPixels);
            const auto n = static_cast<std::size_t>(hi - lo);
            bin2.read_into(lo, hi, std::span<std::int32_t>(slice_columns.data(), n));
            count.read_into(lo, hi, std::span<double>(slice_counts.data(), n));
            for (std::size_t i = 0; i < n; ++i) {
                const std::int64_t position = lo + static_cast<std::int64_t>(i);
                while (index[static_cast<std::size_t>(row) + 1] <= position) {
                    arrays.indptr[static_cast<std::size_t>(row - first) + 1] =
                        static_cast<std::int64_t>(arrays.indices.size());
                    ++row;
                }
                const std::int32_t column = slice_columns[i];
                if (column >= first && column < last) {
                    arrays.indices.push_back(static_cast<std::int32_t>(column - first));
                    arrays.data.push_back(slice_counts[i]);
                }
            }
        }
        for (; row < last; ++row) {
            arrays.indptr[static_cast<std::size_t>(row - first) + 1] =
                static_cast<std::int64_t>(arrays.indices.size());
        }
        sort_row_columns(arrays);
        return CsrMatrix::adopt(std::move(arrays));
    });
}

CsrMatrix CoolFile::read_matrix() const {
    const std::int64_t bins = nbins();
    const std::string dtype = count_dtype();
    if (std::optional<CsrMatrix> matrix = read_from_index(0, bins, dtype)) {
        return std::move(*matrix);
    }
    // Without a usable index: the rows of the pixel table in storage order
    // through the range query.
    CsrMatrix::Arrays arrays;
    arrays.rows = bins;
    arrays.cols = bins;
    arrays.dtype = dtype;
    arrays.indptr.assign(static_cast<std::size_t>(bins) + 1, 0);
    const auto expected = static_cast<std::size_t>(std::max<std::int64_t>(nnz(), 0));
    arrays.indices.reserve(expected);
    arrays.data.reserve(expected);
    for_each_pixel_chunk(0, bins, 0, bins, [&](const PixelChunk& chunk) {
        for (std::size_t k = 0; k < chunk.bin1.size(); ++k) {
            ++arrays.indptr[static_cast<std::size_t>(chunk.bin1[k]) + 1];
            arrays.indices.push_back(static_cast<std::int32_t>(chunk.bin2[k]));
            arrays.data.push_back(chunk.count[k]);
        }
    });
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }
    return CsrMatrix::adopt(std::move(arrays));
}

CsrMatrix CoolFile::read_block(std::int64_t first, std::int64_t last) const {
    const std::int64_t size = std::max<std::int64_t>(last - first, 0);
    CsrMatrix::Arrays arrays;
    arrays.rows = size;
    arrays.cols = size;
    arrays.dtype = count_dtype();
    if (size > 0) {
        if (std::optional<CsrMatrix> block = read_from_index(first, last, arrays.dtype)) {
            return std::move(*block);
        }
    }
    arrays.indptr.assign(static_cast<std::size_t>(size) + 1, 0);
    if (size > 0) {
        for_each_pixel_chunk(first, last, first, last, [&](const PixelChunk& chunk) {
            for (std::size_t k = 0; k < chunk.bin1.size(); ++k) {
                ++arrays.indptr[static_cast<std::size_t>(chunk.bin1[k] - first) + 1];
                arrays.indices.push_back(static_cast<std::int32_t>(chunk.bin2[k] - first));
                arrays.data.push_back(chunk.count[k]);
            }
        });
    }
    for (std::size_t i = 1; i < arrays.indptr.size(); ++i) {
        arrays.indptr[i] += arrays.indptr[i - 1];
    }
    sort_row_columns(arrays);
    return CsrMatrix::adopt(std::move(arrays));
}

CoolLoadResult read_cool(const std::string& uri, const CoolLoadOptions& options) {
    const CoolFile cool(uri);
    CoolLoadResult result;
    for (const auto& [key, value] : cool.info()) {
        result.metadata.emplace(key, value.to_python_string());
    }
    result.data.cut_intervals = cool.read_bins();

    // The bin range of the requested chromosome. The bins of a chromosome are
    // contiguous, so only that block of the pixel table is read. Every step
    // below acts on single entries or single rows, so doing it on the block
    // gives the entries a whole file load followed by the cut would give.
    std::int64_t first = 0;
    std::int64_t last = static_cast<std::int64_t>(result.data.cut_intervals.size());
    if (options.chrom_name.has_value()) {
        if (options.chrom_name->find(':') != std::string::npos) {
            // hicmatrix hands the name to cooler's fetch, which also takes a
            // region string (hicPlotMatrix --region on a cool): the bins that
            // overlap the region, by cooler's region_to_extent.
            std::tie(first, last) = cool.extent(*options.chrom_name);
        } else {
            first = -1;
            last = -1;
            for (std::size_t bin = 0; bin < result.data.cut_intervals.size(); ++bin) {
                if (result.data.cut_intervals[bin].chrom != *options.chrom_name) {
                    continue;
                }
                if (first < 0) {
                    first = static_cast<std::int64_t>(bin);
                }
                last = static_cast<std::int64_t>(bin) + 1;
            }
            if (first < 0) {
                first = 0;
                last = 0;
            }
        }
        result.data.matrix = cool.read_block(first, last);
    } else {
        result.data.matrix = cool.read_matrix();
    }

    std::optional<std::vector<double>> weights;
    if (options.apply_correction && cool.has_column(options.correction_factor_table)) {
        weights = cool.read_column(options.correction_factor_table);
    }

    CsrMatrix& matrix = result.data.matrix;
    if (weights.has_value()) {
        matrix.eliminate_zeros();
        // hicmatrix tests the whole matrix; the block path takes the stored
        // pixel count of the file, which is the same number for a cooler
        // without explicit zeros.
        const std::size_t whole_nnz = options.chrom_name.has_value()
                                          ? static_cast<std::size_t>(cool.nnz())
                                          : matrix.nnz();
        if (whole_nnz > 1) {
            const std::vector<double>& factors = *weights;
            const bool all_nan = std::all_of(factors.begin(), factors.end(),
                                             [](double v) { return std::isnan(v); });
            if (!all_nan) {
                // The whole derivation is skipped when the caller has already
                // fixed the operator (cool.py:195): the version strings are
                // then never read either.
                if (options.correction_operator.has_value()) {
                    result.correction_operator = options.correction_operator;
                } else {
                    // 'weight' is multiplicative, the hic2cool tables KR, VC
                    // and SQRT_VC are divisive.
                    result.correction_operator =
                        (options.correction_factor_table == "KR" ||
                         options.correction_factor_table == "VC" ||
                         options.correction_factor_table == "SQRT_VC")
                            ? '/'
                            : '*';
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
                }
                const bool divide = *result.correction_operator == '/';
                std::vector<double>& values = matrix.mutable_data();
                const std::vector<std::int64_t>& indptr = matrix.indptr();
                const std::vector<std::int32_t>& indices = matrix.indices();
                const auto offset = static_cast<std::size_t>(first);
                for (std::int64_t row = 0; row < matrix.rows(); ++row) {
                    const auto begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                    const auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
                    for (std::size_t k = begin; k < end; ++k) {
                        const double factor =
                            factors[offset + static_cast<std::size_t>(row)] *
                            factors[offset + static_cast<std::size_t>(indices[k])];
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

    if (options.chrom_name.has_value()) {
        result.data.cut_intervals.assign(
            result.data.cut_intervals.begin() + static_cast<std::ptrdiff_t>(first),
            result.data.cut_intervals.begin() + static_cast<std::ptrdiff_t>(last));
        if (result.data.correction_factors.has_value()) {
            const std::vector<double>& factors = *result.data.correction_factors;
            result.data.correction_factors = std::vector<double>(
                factors.begin() + static_cast<std::ptrdiff_t>(first),
                factors.begin() + static_cast<std::ptrdiff_t>(last));
        }
    }

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

// The total that lands in the cool 'sum' attribute when the count column is a
// float column: numpy's sum of the column cooler is handed, which is the
// sequential total of the per part sums (one part up to 10^7 pixels, 10,000
// parts above), each part a pairwise sum over 8192 element buffers in the
// dtype of the column, accumulated in float64. Reproduced on a bounded
// buffer, as the former v4 writer did.
class CountSum {
  public:
    CountSum(std::size_t nnz, bool single_precision) : single_precision_(single_precision) {
        if (nnz > kSplitThreshold) {
            const std::size_t base = nnz / kSplitFactor;
            const std::size_t remainder = nnz % kSplitFactor;
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
        const std::size_t part_length = part_ < parts_.size() ? parts_[part_] : 0;
        if (buffer_.size() == kReduceBuffer || in_part_ == part_length) {
            flush_window();
        }
        if (in_part_ == part_length) {
            total_ += single_precision_ ? static_cast<double>(part_total_float_) : part_total_;
            part_total_ = 0.0;
            part_total_float_ = 0.0F;
            in_part_ = 0;
            ++part_;
            while (part_ < parts_.size() && parts_[part_] == 0) {
                ++part_;
            }
        }
    }

    [[nodiscard]] double total() const { return total_; }

  private:
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

    bool single_precision_;
    std::vector<std::size_t> parts_;
    std::size_t part_ = 0;
    std::size_t in_part_ = 0;
    std::vector<double> buffer_;
    double part_total_ = 0.0;
    float part_total_float_ = 0.0F;
    double total_ = 0.0;
};

// Walks the entries hicmatrix writes, in the order it writes them: every
// stored nonzero, restricted to the upper triangle for pSymmetric, row by row.
class PixelCursor {
  public:
    PixelCursor(const CsrMatrix& matrix, bool upper_only, bool enforce_integer,
                coolercpp::DType count_dtype)
        : matrix_(matrix), upper_only_(upper_only), enforce_integer_(enforce_integer),
          count_dtype_(count_dtype) {}

    // The next `count` entries as a pixel table with int32 bin IDs and the
    // count column in the dtype the Python data frame holds.
    coolercpp::Table next(std::size_t count) {
        std::vector<std::int32_t> bin1;
        std::vector<std::int32_t> bin2;
        std::vector<double> values;
        bin1.reserve(count);
        bin2.reserve(count);
        values.reserve(count);
        const std::vector<std::int64_t>& indptr = matrix_.indptr();
        const std::vector<std::int32_t>& indices = matrix_.indices();
        const std::vector<double>& data = matrix_.data();
        while (bin1.size() < count && row_ < matrix_.rows()) {
            const auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row_) + 1]);
            if (k_ < static_cast<std::size_t>(indptr[static_cast<std::size_t>(row_)])) {
                k_ = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row_)]);
            }
            while (k_ < end && bin1.size() < count) {
                const std::size_t k = k_++;
                if (data[k] == 0.0) {
                    continue;  // eliminate_zeros
                }
                if (upper_only_ && indices[k] < row_) {
                    continue;  // triu(k=0)
                }
                bin1.push_back(static_cast<std::int32_t>(row_));
                bin2.push_back(indices[k]);
                // np.rint is round half to even, which is nearbyint under the
                // default rounding mode.
                values.push_back(enforce_integer_ ? std::nearbyint(data[k]) : data[k]);
            }
            if (k_ >= end) {
                ++row_;
            }
        }
        coolercpp::Table table;
        table.set("bin1_id", coolercpp::Column(std::move(bin1)));
        table.set("bin2_id", coolercpp::Column(std::move(bin2)));
        switch (count_dtype_) {
            case coolercpp::DType::Int32: {
                std::vector<std::int32_t> ints(values.size());
                for (std::size_t i = 0; i < values.size(); ++i) {
                    ints[i] = static_cast<std::int32_t>(values[i]);
                }
                table.set("count", coolercpp::Column(std::move(ints)));
                break;
            }
            case coolercpp::DType::Int64: {
                std::vector<std::int64_t> ints(values.size());
                for (std::size_t i = 0; i < values.size(); ++i) {
                    ints[i] = static_cast<std::int64_t>(values[i]);
                }
                table.set("count", coolercpp::Column(std::move(ints)));
                break;
            }
            case coolercpp::DType::Float32: {
                std::vector<float> floats(values.size());
                for (std::size_t i = 0; i < values.size(); ++i) {
                    floats[i] = static_cast<float>(values[i]);
                }
                table.set("count", coolercpp::Column(std::move(floats)));
                break;
            }
            default: table.set("count", coolercpp::Column(std::move(values))); break;
        }
        return table;
    }

  private:
    const CsrMatrix& matrix_;
    bool upper_only_;
    bool enforce_integer_;
    coolercpp::DType count_dtype_;
    std::int64_t row_ = 0;
    std::size_t k_ = 0;
};

}  // namespace

void write_cool(const std::string& uri, MatrixData& data, const CoolSaveOptions& options) {
    CsrMatrix& matrix = data.matrix;
    const std::int64_t nbins = matrix.rows();
    if (data.cut_intervals.size() != static_cast<std::size_t>(nbins)) {
        throw h5::Error("the bin table has " + std::to_string(data.cut_intervals.size()) +
                        " entries but the matrix has " + std::to_string(nbins) + " rows");
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
            const auto begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
            const auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
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
                const auto begin = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                const auto end = static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
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

    // ---- the tables handed to cooler ----
    const std::size_t nnz =
        options.symmetric ? matrix.upper_triangle_nnz() : matrix.nonzero_stored_nnz();
    const DType kind = matrix.dtype_kind();
    // cool.py:355-364: an integer matrix keeps cooler's int32 count column, a
    // float matrix carries its own dtype into the file.
    coolercpp::DType count_file = coolercpp::DType::Int32;
    if (!options.enforce_integer) {
        if (kind == DType::Float32) {
            count_file = coolercpp::DType::Float32;
        } else if (kind == DType::Float64) {
            count_file = coolercpp::DType::Float64;
        }
    }
    // The data frame's count column: the matrix dtype, or float64 after
    // np.rint. cooler sums it before the cast to the file dtype.
    // An integer column is held as int32 when every value fits: cooler's sum
    // of an int32 and of an int64 column is the same int64 total.
    coolercpp::DType count_frame = coolercpp::DType::Float64;
    if (!options.enforce_integer) {
        if (kind == DType::Integer) {
            count_frame = coolercpp::DType::Int32;
            for (const double value : matrix.data()) {
                if (value < -2147483648.0 || value > 2147483647.0) {
                    count_frame = coolercpp::DType::Int64;
                    break;
                }
            }
        } else if (kind == DType::Float32) {
            count_frame = coolercpp::DType::Float32;
        }
    }

    coolercpp::Table bins;
    {
        // The chrom column as a categorical in order of first appearance, so
        // that a bin table of many contigs holds each name once.
        std::vector<std::string> names;
        std::map<std::string, std::int32_t> code_of;
        std::vector<std::int32_t> codes;
        std::vector<std::int64_t> starts;
        std::vector<std::int64_t> ends;
        codes.reserve(data.cut_intervals.size());
        starts.reserve(data.cut_intervals.size());
        ends.reserve(data.cut_intervals.size());
        for (const CutInterval& bin : data.cut_intervals) {
            const auto [it, inserted] =
                code_of.emplace(bin.chrom, static_cast<std::int32_t>(names.size()));
            if (inserted) {
                names.push_back(bin.chrom);
            }
            codes.push_back(it->second);
            starts.push_back(bin.start);
            ends.push_back(bin.end);
        }
        bins.set("chrom", coolercpp::Column::categorical(std::move(codes), std::move(names)));
        bins.set("start", coolercpp::Column(std::move(starts)));
        bins.set("end", coolercpp::Column(std::move(ends)));
        if (write_weight) {
            // convertNansToOnes
            std::vector<double> weight = *data.correction_factors;
            for (double& factor : weight) {
                if (std::isnan(factor)) {
                    factor = 1.0;
                }
            }
            bins.set("weight", coolercpp::Column(std::move(weight)));
        }
    }

    // The info dictionary hicmatrix builds (cool.py:375-402) in its insertion
    // order, which is also the key order of the 'metadata' JSON string.
    std::vector<std::pair<std::string, std::string>> info{
        {"format", "HDF5::Cooler"},
        {"format-url", options.format_url},
        {"generated-by", options.generated_by},
        {"generated-by-cooler-lib", options.generated_by_cooler_lib},
        {"tool-url", options.tool_url},
    };
    if (options.has_hic_metadata) {
        for (const char* key :
             {"matrix-generated-by", "matrix-generated-by-url", "genome-assembly"}) {
            const auto found = options.hic_metadata.find(key);
            if (found != options.hic_metadata.end()) {
                info.emplace_back(key, found->second);
            }
        }
    }

    coolercpp::CreateOptions create;
    create.dtypes = {{"bin1_id", coolercpp::DType::Int32},
                     {"bin2_id", coolercpp::DType::Int32},
                     {"count", count_file}};
    coolercpp::json::Value metadata = coolercpp::json::Value::object();
    for (const auto& [key, value] : info) {
        metadata[key] = value;
    }
    create.metadata = metadata;
    create.ordered = true;
    create.symmetric_upper = options.symmetric;
    create.mode = options.append ? "a" : "w";
    if (!options.creation_date.empty()) {
        create.creation_date = options.creation_date;
    }
    // cooler's own identity on the cooler group, as cooler 0.10.2 writes it;
    // hicmatrix replaces it on the file root below.
    create.generated_by = options.generated_by_cooler_lib;

    // The written entries in order: every row's columns strictly increasing,
    // which is the order create_cooler would sort a data frame into, and the
    // float total of the count column hicmatrix's data frame(s) produce.
    const bool integer_sum = kind == DType::Integer && !options.enforce_integer;
    const bool float32_sum = kind == DType::Float32 && !options.enforce_integer;
    bool rows_sorted = true;
    bool ids_in_bounds = true;
    std::optional<double> float_total;
    {
        std::optional<CountSum> sum;
        if (!integer_sum) {
            sum.emplace(nnz, float32_sum);
        }
        std::int64_t last_row = -1;
        std::int64_t last_col = -1;
        const auto visit = [&](std::int64_t row, std::int64_t col, double value) {
            if (row == last_row && col <= last_col) {
                rows_sorted = false;
            }
            if (col < 0 || col >= nbins) {
                ids_in_bounds = false;
            }
            last_row = row;
            last_col = col;
            if (sum.has_value()) {
                sum->add(options.enforce_integer ? std::nearbyint(value) : value);
            }
        };
        if (options.symmetric) {
            matrix.for_each_upper(visit);
        } else {
            matrix.for_each_stored(visit);
        }
        if (sum.has_value()) {
            float_total = sum->total();
        }
    }

    // create_cooler's boundscheck, dupcheck and triucheck only ever raise;
    // they never change what is written. The pass above has already
    // established all three for the entries the cursor yields: rows come in
    // order with strictly increasing columns (no duplicate pixel), every
    // column lies inside the bin table (rows do by construction), and the
    // symmetric case yields the upper triangle only (triucheck applies to
    // symmetric-upper files alone). When that holds, the checks are skipped,
    // which is 0.7 s of CPU on the 61.8 M pixel gm12878 matrix; otherwise they
    // run, and bad input raises as it does in Python.
    if (rows_sorted && ids_in_bounds) {
        create.boundscheck = false;
        create.dupcheck = false;
        create.triucheck = false;
    }

    PixelCursor cursor(matrix, options.symmetric, options.enforce_integer, count_frame);
    guarded([&] {
        if (nnz > kSplitThreshold) {
            // np.array_split(df, 10000): the first nnz % 10000 parts are one
            // longer; each part is one create_cooler chunk.
            const std::size_t base = nnz / kSplitFactor;
            const std::size_t remainder = nnz % kSplitFactor;
            std::size_t part = 0;
            coolercpp::PixelChunks chunks = [&]() -> std::optional<coolercpp::Table> {
                if (part >= kSplitFactor) {
                    return std::nullopt;
                }
                const std::size_t length = base + (part < remainder ? 1 : 0);
                ++part;
                return cursor.next(length);
            };
            coolercpp::create_cooler(uri, bins, chunks, create);
        } else if (rows_sorted && nnz > 0) {
            // hicmatrix hands cooler one data frame, which create_cooler sorts
            // and writes as a single chunk. The entries are already in that
            // order, so they are streamed in blocks, which keeps a copy of the
            // pixel table out of memory; the float 'sum' of a single part is
            // put back below.
            std::size_t written = 0;
            coolercpp::PixelChunks chunks = [&]() -> std::optional<coolercpp::Table> {
                if (written >= nnz) {
                    return std::nullopt;
                }
                const std::size_t length = std::min(kWriteBlock, nnz - written);
                written += length;
                return cursor.next(length);
            };
            coolercpp::create_cooler(uri, bins, chunks, create);
        } else {
            const coolercpp::Table pixels = cursor.next(nnz);
            coolercpp::create_cooler(uri, bins, pixels, create);
        }
    });

    // Cool.save then overwrites the provenance on the *file root*, and only in
    // mode 'w' (cool.py:422-426). For a plain .cool the root is the cooler
    // group; for an mcool it is a different object, so the resolution groups
    // keep cooler's own provenance.
    const auto [path, group] = coolercpp::parse_cooler_uri(uri);
    if (float_total.has_value() || !options.append) {
        h5::FileWriter file(path, h5::WriteMode::Append);
        if (float_total.has_value()) {
            file.set_attribute(group, "sum", *float_total);
        }
        if (!options.append) {
            for (const auto& [key, value] : info) {
                file.set_attribute("/", key, value);
            }
        }
    }
}

}  // namespace hicx
