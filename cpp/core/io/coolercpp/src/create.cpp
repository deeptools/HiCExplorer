// Port of cooler 0.10.2 (BSD-3-Clause):
//   cooler/create/_create.py   create, create_cooler, create_from_unordered,
//                              write_chroms, write_bins, prepare_pixels,
//                              write_pixels, index_bins, index_pixels,
//                              write_indexes, write_info, _set_h5opts
//   cooler/create/_ingest.py   validate_pixels
//   cooler/core/_tableops.py   put
//   cooler/reduce.py           merge_breakpoints, CoolerMerger
//   cooler/util.py             rlencode
// together with the h5py rules they rely on (h5py/_hl/dataset.py
// make_new_dset, h5py/_hl/filters.py fill_dcpl).

#include "coolercpp/create.hpp"

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "coolercpp/api.hpp"
#include "coolercpp/errors.hpp"
#include "coolercpp/util.hpp"
#include "coolercpp/version.hpp"
#include "h5.hpp"
#include "numpy_compat.hpp"

namespace coolercpp {

PixelChunks iterate_chunks(std::vector<Table> chunks) {
    auto state = std::make_shared<std::pair<std::vector<Table>, std::size_t>>(std::move(chunks), 0);
    return [state]() -> std::optional<Table> {
        if (state->second >= state->first.size()) {
            return std::nullopt;
        }
        return std::move(state->first[state->second++]);
    };
}

namespace {

std::string join_path(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

// ---------------------------------------------------------------------------
// h5opts and h5py's dataset creation

class Kwargs {
  public:
    Kwargs() = default;
    Kwargs(std::initializer_list<std::pair<std::string, H5OptValue>> items) {
        for (const auto& [k, v] : items) {
            set(k, v);
        }
    }
    [[nodiscard]] const H5OptValue* get(const std::string& key) const {
        for (const auto& [k, v] : items_) {
            if (k == key) {
                return &v;
            }
        }
        return nullptr;
    }
    [[nodiscard]] bool has(const std::string& key) const { return get(key) != nullptr; }
    void set(const std::string& key, H5OptValue value) {
        for (auto& [k, v] : items_) {
            if (k == key) {
                v = std::move(value);
                return;
            }
        }
        items_.emplace_back(key, std::move(value));
    }
    [[nodiscard]] const std::vector<std::pair<std::string, H5OptValue>>& items() const {
        return items_;
    }

  private:
    std::vector<std::pair<std::string, H5OptValue>> items_;
};

bool is_none(const H5OptValue* v) { return v == nullptr || std::holds_alternative<std::monostate>(*v); }

bool truthy(const H5OptValue* v) {
    if (v == nullptr) {
        return false;
    }
    return std::visit(
        [](const auto& x) -> bool {
            using X = std::decay_t<decltype(x)>;
            if constexpr (std::is_same_v<X, std::monostate>) {
                return false;
            } else if constexpr (std::is_same_v<X, bool>) {
                return x;
            } else if constexpr (std::is_same_v<X, std::int64_t> || std::is_same_v<X, double>) {
                return x != 0;
            } else {
                return !x.empty();
            }
        },
        *v);
}

// A Python int (bool included) out of an option value.
std::optional<std::int64_t> as_int(const H5OptValue* v) {
    if (v == nullptr) {
        return std::nullopt;
    }
    if (const auto* b = std::get_if<bool>(v)) {
        return *b ? 1 : 0;
    }
    if (const auto* i = std::get_if<std::int64_t>(v)) {
        return *i;
    }
    return std::nullopt;
}

std::string repr(const H5OptValue* v) {
    if (v == nullptr) {
        return "None";
    }
    return std::visit(
        [](const auto& x) -> std::string {
            using X = std::decay_t<decltype(x)>;
            if constexpr (std::is_same_v<X, std::monostate>) {
                return "None";
            } else if constexpr (std::is_same_v<X, bool>) {
                return x ? "True" : "False";
            } else if constexpr (std::is_same_v<X, std::int64_t>) {
                return std::to_string(x);
            } else if constexpr (std::is_same_v<X, double>) {
                return npy::float_repr(x);
            } else if constexpr (std::is_same_v<X, std::string>) {
                return "'" + x + "'";
            } else {
                std::string out = "(";
                for (std::size_t i = 0; i < x.size(); ++i) {
                    out += (x[i] < 0 ? std::string("None") : std::to_string(x[i]));
                    out += (x.size() == 1 || i + 1 < x.size()) ? "," : "";
                    if (i + 1 < x.size()) {
                        out += " ";
                    }
                }
                return out + ")";
            }
        },
        *v);
}

// cooler.create._create._set_h5opts
Kwargs set_h5opts(const std::optional<H5Opts>& user) {
    Kwargs result;
    if (user.has_value()) {
        for (const auto& [k, v] : *user) {
            result.set(k, v);
        }
    }
    static const std::set<std::string> kAvailable{"chunks",      "maxshape",   "compression",
                                                  "compression_opts", "scaleoffset", "shuffle",
                                                  "fletcher32", "fillvalue",  "track_times"};
    for (const auto& [k, v] : result.items()) {
        if (kAvailable.count(k) == 0) {
            throw ValueError("Unknown storage option '" + k + "'.");
        }
    }
    if (!result.has("compression")) {
        result.set("compression", std::string("gzip"));
    }
    const H5OptValue* compression = result.get("compression");
    if (std::holds_alternative<std::string>(*compression) &&
        std::get<std::string>(*compression) == "gzip" && !result.has("compression_opts")) {
        result.set("compression_opts", std::int64_t{6});
    }
    if (!result.has("shuffle")) {
        result.set("shuffle", true);
    }
    return result;
}

struct TypeSpec {
    hid_t type = -1;
    std::size_t itemsize = 0;
    // numpy dtype kind: 'i', 'u', 'f', 'b' (bool), 'S', 'e' (enum)
    char kind = 'i';
};

TypeSpec numeric_spec(const h5::Handle& type, DType dtype) {
    TypeSpec spec;
    spec.type = type.get();
    spec.itemsize = itemsize(dtype);
    spec.kind = dtype == DType::Bool       ? 'b'
                : is_float(dtype)          ? 'f'
                : is_unsigned_integer(dtype) ? 'u'
                                           : 'i';
    return spec;
}

// "maxshape=" as a caller of create_dataset passes it: absent, or a length
// with nullopt meaning None (unlimited).
using MaxShapeArg = std::optional<std::optional<std::size_t>>;

// h5py Group.create_dataset(path, shape=(length,), dtype, maxshape=...,
// fillvalue=..., **kwargs), following make_new_dset and fill_dcpl.
h5::Dataset create_dataset(h5::File& file, const std::string& path, const TypeSpec& spec,
                           std::size_t length, MaxShapeArg explicit_maxshape,
                           std::optional<json::Value> explicit_fill, const Kwargs& kwargs) {
    if (explicit_maxshape.has_value() && kwargs.has("maxshape")) {
        throw TypeError("create_dataset() got multiple values for keyword argument 'maxshape'");
    }
    if (explicit_fill.has_value() && kwargs.has("fillvalue")) {
        throw TypeError("create_dataset() got multiple values for keyword argument 'fillvalue'");
    }
    bool has_maxshape = false;
    std::optional<std::size_t> maxlen;
    if (explicit_maxshape.has_value()) {
        has_maxshape = true;
        maxlen = *explicit_maxshape;
    } else if (const H5OptValue* m = kwargs.get("maxshape"); !is_none(m)) {
        has_maxshape = true;
        if (const auto opt = as_int(m)) {
            maxlen = static_cast<std::size_t>(*opt);  // an int n means (n,)
        } else if (const auto* tuple = std::get_if<std::vector<std::int64_t>>(m)) {
            if (tuple->size() != 1) {
                throw ValueError("\"maxshape\" must have same rank as dataset shape");
            }
            if ((*tuple)[0] >= 0) {
                maxlen = static_cast<std::size_t>((*tuple)[0]);
            }
        } else {
            throw TypeError("\"maxshape\" argument must be None or a sequence object");
        }
    }

    // chunks
    const H5OptValue* chunks = kwargs.get("chunks");
    bool chunks_true = false;
    bool chunks_false = false;
    std::optional<std::size_t> user_chunk;
    if (!is_none(chunks)) {
        if (const auto* b = std::get_if<bool>(chunks)) {
            (*b ? chunks_true : chunks_false) = true;
        } else if (const auto* i = std::get_if<std::int64_t>(chunks)) {
            user_chunk = static_cast<std::size_t>(*i);
        } else if (const auto* tuple = std::get_if<std::vector<std::int64_t>>(chunks)) {
            if (tuple->size() != 1) {
                throw ValueError("\"chunks\" must have same rank as dataset shape");
            }
            user_chunk = static_cast<std::size_t>((*tuple)[0]);
        } else {
            throw TypeError("\"chunks\" argument must be None or a sequence object");
        }
    }
    if (user_chunk.has_value()) {
        const std::optional<std::size_t> dim = has_maxshape ? maxlen : std::optional(length);
        if (dim.has_value() && *user_chunk > *dim) {
            throw ValueError(
                "Chunk shape must not be greater than data shape in any dimension. (" +
                std::to_string(*user_chunk) + ",) is not compatible with (" +
                std::to_string(length) + ",)");
        }
    }

    H5OptValue compression = kwargs.get("compression") ? *kwargs.get("compression")
                                                       : H5OptValue{};
    H5OptValue compression_opts = kwargs.get("compression_opts")
                                      ? *kwargs.get("compression_opts")
                                      : H5OptValue{};
    const H5OptValue* shuffle = kwargs.get("shuffle");
    const H5OptValue* fletcher32 = kwargs.get("fletcher32");
    const H5OptValue* scaleoffset = kwargs.get("scaleoffset");

    // Legacy checks of make_new_dset.
    if ((truthy(&compression) || truthy(shuffle) || truthy(fletcher32) || has_maxshape ||
         truthy(scaleoffset)) &&
        chunks_false) {
        throw ValueError("Chunked format required for given storage options");
    }
    if (const auto* b = std::get_if<bool>(&compression); b != nullptr && *b) {
        if (is_none(&compression_opts)) {
            compression_opts = std::int64_t{4};
        }
        compression = std::string("gzip");
    }
    if (const auto level = as_int(&compression); level.has_value() && *level >= 0 && *level < 10) {
        if (!is_none(&compression_opts)) {
            throw TypeError("Conflict in compression options");
        }
        compression_opts = *level;
        compression = std::string("gzip");
    }
    if (chunks_false) {
        throw TypeError("\"chunks\" argument must be None or a sequence object");
    }

    h5::DatasetCreate dc;
    dc.length = length;
    if (has_maxshape) {
        dc.maxlength = maxlen.has_value() ? *maxlen : h5::DatasetCreate::kUnlimited;
    }
    if (!is_none(&compression)) {
        if (const auto* name = std::get_if<std::string>(&compression)) {
            if (*name == "gzip") {
                std::int64_t level = 4;
                if (!is_none(&compression_opts)) {
                    const auto opt = as_int(&compression_opts);
                    if (!opt.has_value() || *opt < 0 || *opt > 9) {
                        throw ValueError("GZIP setting must be an integer from 0-9, not " +
                                         repr(&compression_opts));
                    }
                    level = *opt;
                }
                dc.compression = h5::DatasetCreate::Compression::Gzip;
                dc.gzip_level = static_cast<unsigned>(level);
            } else if (*name == "lzf") {
                if (!is_none(&compression_opts)) {
                    throw ValueError("LZF compression filter accepts no options");
                }
                // Deviation: h5py ships its own LZF filter; libhdf5 does not.
                throw ValueError("coolercpp does not provide the LZF filter (h5py filter 32000)");
            } else if (*name == "szip") {
                if (H5Zfilter_avail(H5Z_FILTER_SZIP) <= 0) {
                    throw ValueError("Compression filter \"szip\" is unavailable");
                }
                if (!is_none(&compression_opts)) {
                    throw TypeError("coolercpp accepts only the default SZIP options");
                }
                dc.compression = h5::DatasetCreate::Compression::Szip;
                dc.szip_mask = H5_SZIP_NN_OPTION_MASK;
                dc.szip_pixels = 8;
            } else {
                throw ValueError("Compression filter \"" + *name + "\" is unavailable");
            }
        } else if (const auto id = as_int(&compression)) {
            if (H5Zfilter_avail(static_cast<H5Z_filter_t>(*id)) <= 0) {
                throw ValueError("Unknown compression filter number: " + std::to_string(*id));
            }
            dc.compression = h5::DatasetCreate::Compression::FilterId;
            dc.filter_id = static_cast<int>(*id);
            if (const auto* opts = std::get_if<std::vector<std::int64_t>>(&compression_opts)) {
                for (const std::int64_t o : *opts) {
                    dc.filter_options.push_back(static_cast<unsigned>(o));
                }
            } else if (const auto single = as_int(&compression_opts)) {
                dc.filter_options.push_back(static_cast<unsigned>(*single));
            }
        } else {
            throw ValueError("Compression filter \"" + repr(&compression) + "\" is unavailable");
        }
    } else if (!is_none(&compression_opts)) {
        throw TypeError("Compression method must be specified");
    }

    std::optional<std::int64_t> so;
    if (!is_none(scaleoffset)) {
        const auto value = as_int(scaleoffset);
        if (!value.has_value()) {
            throw TypeError("scale/offset filter requires an integer");
        }
        const bool is_true = std::holds_alternative<bool>(*scaleoffset) && *value == 1;
        if (*value < 0) {
            throw ValueError("scale factor must be >= 0");
        }
        if (spec.kind == 'f') {
            if (is_true) {
                throw ValueError(
                    "integer scaleoffset must be provided for floating point types");
            }
            so = *value;
        } else if (spec.kind == 'i' || spec.kind == 'u') {
            so = is_true ? 0 : *value;
        } else {
            throw TypeError(
                "scale/offset filter only supported for integer and floating-point types");
        }
        if (truthy(fletcher32)) {
            throw ValueError("fletcher32 cannot be used with potentially lossy scale/offset filter");
        }
    }

    if (chunks_true || (!user_chunk.has_value() &&
                        (truthy(shuffle) || truthy(fletcher32) || !is_none(&compression) ||
                         has_maxshape || so.has_value()))) {
        dc.chunk = h5::guess_chunk(length, spec.itemsize);
    } else if (user_chunk.has_value()) {
        dc.chunk = user_chunk;
    }
    dc.scaleoffset = so;
    dc.shuffle = truthy(shuffle);
    dc.fletcher32 = truthy(fletcher32);

    if (explicit_fill.has_value()) {
        dc.fillvalue = explicit_fill;
    } else if (const H5OptValue* fill = kwargs.get("fillvalue"); !is_none(fill)) {
        if (const auto i = as_int(fill)) {
            dc.fillvalue = json::Value(*i);
        } else if (const auto* d = std::get_if<double>(fill)) {
            dc.fillvalue = json::Value(*d);
        } else {
            throw TypeError("coolercpp accepts numeric fill values only");
        }
    }
    if (const H5OptValue* tt = kwargs.get("track_times"); !is_none(tt)) {
        const auto value = as_int(tt);
        if (!value.has_value() || (*value != 0 && *value != 1)) {
            throw TypeError("track_times must be either True or False");
        }
        dc.track_times = *value == 1;
    }
    return file.create_dataset(path, spec.type, dc);
}

// cooler.core.put for a fresh column of a table group (default h5opts: gzip
// level 6 without shuffle, maxshape (None,)).
void put_column(h5::File& file, const std::string& group, const std::string& name,
                const Column& column) {
    const Kwargs opts{{"compression", std::string("gzip")}, {"compression_opts", std::int64_t{6}}};
    const std::string path = join_path(group, name);
    const std::size_t n = column.size();
    const MaxShapeArg unlimited = std::optional<std::size_t>{};
    switch (column.dtype()) {
        case DType::Categorical: {
            const CategoricalData& data = column.categorical();
            const std::vector<std::string>& categories = *data.categories;
            const std::size_t count = categories.size();
            // pandas coerce_indexer_dtype
            const DType base = count < 127 ? DType::Int8
                               : count < 32767 ? DType::Int16
                               : count < 2147483647U ? DType::Int32
                                                     : DType::Int64;
            const h5::Handle type = h5::enum_type(categories, base);
            TypeSpec spec{type.get(), itemsize(base), 'e'};
            h5::Dataset dataset =
                create_dataset(file, path, spec, n, unlimited, json::Value(-1), opts);
            const h5::Handle native(H5Tget_native_type(type.get(), H5T_DIR_ASCEND),
                                    h5::Handle::Kind::DataType);
            std::vector<unsigned char> buffer(n * spec.itemsize);
            for (std::size_t i = 0; i < n; ++i) {
                const std::int64_t code = data.codes[i];
                switch (spec.itemsize) {
                    case 1: { const auto v = static_cast<std::int8_t>(code); std::memcpy(&buffer[i], &v, 1); break; }
                    case 2: { const auto v = static_cast<std::int16_t>(code); std::memcpy(&buffer[i * 2], &v, 2); break; }
                    case 4: { const auto v = static_cast<std::int32_t>(code); std::memcpy(&buffer[i * 4], &v, 4); break; }
                    default: std::memcpy(&buffer[i * 8], &code, 8); break;
                }
            }
            dataset.write_raw(0, n, native.get(), buffer.data());
            return;
        }
        case DType::String: {
            std::size_t width = 1;
            for (const std::string& s : column.values<std::string>()) {
                width = std::max(width, s.size());
            }
            const h5::Handle type = h5::fixed_string_type(width);
            h5::Dataset dataset = create_dataset(file, path, TypeSpec{type.get(), width, 'S'}, n,
                                                 unlimited, std::nullopt, opts);
            dataset.write_column(0, column);
            return;
        }
        default: {
            const h5::Handle type = h5::file_type(column.dtype());
            h5::Dataset dataset = create_dataset(file, path, numeric_spec(type, column.dtype()),
                                                 n, unlimited, std::nullopt, opts);
            dataset.write_column(0, column);
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// Pixel validation (cooler.create._ingest._validate_pixels)

// Bin ID column values for comparisons, without copying int32 or int64 input.
class Ids {
  public:
    explicit Ids(const Column& column) {
        const DType dtype = column.dtype();
        if (dtype == DType::Int64) {
            ints_ = &column.values<std::int64_t>();
            size_ = ints_->size();
        } else if (dtype == DType::Int32) {
            ints32_ = &column.values<std::int32_t>();
            size_ = ints32_->size();
        } else if (is_float(dtype)) {
            floating_ = true;
            doubles_ = column.as<double>();
            size_ = doubles_.size();
        } else if (is_numeric(dtype)) {
            owned_ = column.as<std::int64_t>();
            ints_ = &owned_;
            size_ = owned_.size();
        } else {
            throw TypeError("'<' not supported between instances of 'str' and 'int'");
        }
    }
    [[nodiscard]] bool floating() const noexcept { return floating_; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] double d(std::size_t k) const {
        return floating_ ? doubles_[k] : static_cast<double>(i(k));
    }
    [[nodiscard]] std::int64_t i(std::size_t k) const {
        return ints32_ != nullptr ? static_cast<std::int64_t>((*ints32_)[k]) : (*ints_)[k];
    }

  private:
    bool floating_ = false;
    std::size_t size_ = 0;
    const std::vector<std::int64_t>* ints_ = nullptr;
    const std::vector<std::int32_t>* ints32_ = nullptr;
    std::vector<std::int64_t> owned_;
    std::vector<double> doubles_;
};

// Lexicographic (bin1, bin2) comparison of rows a and b.
struct PairLess {
    const Ids& b1;
    const Ids& b2;
    bool floating;
    bool operator()(std::size_t a, std::size_t b) const {
        if (floating) {
            if (b1.d(a) != b1.d(b)) {
                return b1.d(a) < b1.d(b);
            }
            return b2.d(a) < b2.d(b);
        }
        if (b1.i(a) != b1.i(b)) {
            return b1.i(a) < b1.i(b);
        }
        return b2.i(a) < b2.i(b);
    }
};

bool sorted_pairs(const Ids& b1, const Ids& b2) {
    const PairLess less{b1, b2, b1.floating() || b2.floating()};
    for (std::size_t k = 1; k < b1.size(); ++k) {
        if (less(k, k - 1)) {
            return false;
        }
    }
    return true;
}

// DataFrame.sort_values(["bin1_id", "bin2_id"]): a stable lexicographic
// sort. Returns nothing when the table is already in order.
std::optional<Table> sort_by_pixel_ids(const Table& table) {
    const Ids b1(table["bin1_id"]);
    const Ids b2(table["bin2_id"]);
    if (sorted_pairs(b1, b2)) {
        return std::nullopt;
    }
    std::vector<std::size_t> order(b1.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), PairLess{b1, b2, b1.floating() || b2.floating()});
    const std::vector<std::int64_t> positions(order.begin(), order.end());
    return table.take(positions);
}

std::string csv_field(std::string text) {
    if (text.find_first_of("\t\"\n\r") != std::string::npos) {
        std::string quoted = "\"";
        for (const char c : text) {
            if (c == '"') {
                quoted += "\"\"";
            } else {
                quoted.push_back(c);
            }
        }
        return quoted + "\"";
    }
    return text;
}

std::string csv_value(const Column& column, std::size_t k) {
    switch (column.dtype()) {
        case DType::Bool: return column.values<Bool8>()[k] ? "True" : "False";
        case DType::Float64: {
            const double v = column.values<double>()[k];
            return std::isnan(v) ? std::string() : npy::float_repr(v);
        }
        case DType::Float32: {
            const float v = column.values<float>()[k];
            return std::isnan(v) ? std::string() : npy::float32_repr(v);
        }
        case DType::String:
        case DType::Categorical: return csv_field(column.label(k));
        case DType::UInt64: return std::to_string(column.values<std::uint64_t>()[k]);
        default: return std::to_string(column.as_int64(k));
    }
}

// err.head().to_csv(sep="\t")
std::string duplicates_csv(const Table& chunk, const std::vector<std::size_t>& rows) {
    std::string out;
    for (const std::string& name : chunk.columns()) {
        out += "\t" + csv_field(name);
    }
    out += "\n";
    for (std::size_t r = 0; r < rows.size() && r < 5; ++r) {
        out += std::to_string(chunk.index()[rows[r]]);
        for (std::size_t c = 0; c < chunk.num_columns(); ++c) {
            out += "\t" + csv_value(chunk.column(c), rows[r]);
        }
        out += "\n";
    }
    return out;
}

struct PairHash {
    std::size_t operator()(const std::pair<std::int64_t, std::int64_t>& p) const noexcept {
        return std::hash<std::int64_t>()(p.first) * 1000003U ^ std::hash<std::int64_t>()(p.second);
    }
};

// Returns a sorted copy when ensure_sorted had to reorder the chunk.
std::optional<Table> validate_pixels(const Table& chunk, std::int64_t n_bins, bool boundscheck,
                                     bool triucheck, bool dupcheck, bool ensure_sorted) {
    const auto need = [&chunk](const char* name) -> const Column& {
        if (!chunk.contains(name)) {
            throw KeyError(name);
        }
        return chunk[name];
    };
    if (boundscheck) {
        const Ids b1(need("bin1_id"));
        const Ids b2(need("bin2_id"));
        const auto nb = static_cast<double>(n_bins);
        for (std::size_t k = 0; k < b1.size(); ++k) {
            if (b1.d(k) < 0 || b2.d(k) < 0) {
                throw BadInputError("Found bin ID < 0");
            }
        }
        for (std::size_t k = 0; k < b1.size(); ++k) {
            const bool excess = (b1.floating() || b2.floating())
                                    ? (b1.d(k) >= nb || b2.d(k) >= nb)
                                    : (b1.i(k) >= n_bins || b2.i(k) >= n_bins);
            if (excess) {
                throw BadInputError(
                    "Found a bin ID that exceeds the declared number of bins. "
                    "Check whether your bin table is correct.");
            }
        }
    }
    if (triucheck) {
        const Ids b1(need("bin1_id"));
        const Ids b2(need("bin2_id"));
        const bool floating = b1.floating() || b2.floating();
        for (std::size_t k = 0; k < b1.size(); ++k) {
            if (floating ? b1.d(k) > b2.d(k) : b1.i(k) > b2.i(k)) {
                throw BadInputError("Found bin1_id greater than bin2_id");
            }
        }
    }
    if (dupcheck) {
        const Ids b1(need("bin1_id"));
        const Ids b2(need("bin2_id"));
        const bool floating = b1.floating() || b2.floating();
        std::vector<std::size_t> duplicates;
        if (sorted_pairs(b1, b2)) {
            for (std::size_t k = 1; k < b1.size(); ++k) {
                const bool same = floating ? (b1.d(k) == b1.d(k - 1) && b2.d(k) == b2.d(k - 1))
                                           : (b1.i(k) == b1.i(k - 1) && b2.i(k) == b2.i(k - 1));
                if (same) {
                    duplicates.push_back(k);
                }
            }
        } else {
            std::unordered_set<std::pair<std::int64_t, std::int64_t>, PairHash> seen;
            seen.reserve(b1.size());
            for (std::size_t k = 0; k < b1.size(); ++k) {
                std::pair<std::int64_t, std::int64_t> key;
                if (floating) {
                    const double x = b1.d(k);
                    const double y = b2.d(k);
                    std::memcpy(&key.first, &x, 8);
                    std::memcpy(&key.second, &y, 8);
                } else {
                    key = {b1.i(k), b2.i(k)};
                }
                if (!seen.insert(key).second) {
                    duplicates.push_back(k);
                }
            }
        }
        if (!duplicates.empty()) {
            throw BadInputError("Found duplicate pixels:\n" + duplicates_csv(chunk, duplicates));
        }
    }
    if (ensure_sorted) {
        need("bin1_id");
        need("bin2_id");
        return sort_by_pixel_ids(chunk);
    }
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// The running 'sum' attribute: total = 0; total += chunk["count"].sum()

class PyTotal {
  public:
    void add(const Column& count) {
        switch (count.dtype()) {
            case DType::Float64: {
                const auto& v = count.values<double>();
                add_float(npy::pairwise_sum(v.data(), v.size()));
                return;
            }
            case DType::Float32: {
                const auto& v = count.values<float>();
                add_float(static_cast<double>(npy::pairwise_sum(v.data(), v.size())));
                return;
            }
            case DType::Bool:
            case DType::Int8:
            case DType::Int16:
            case DType::Int32:
            case DType::Int64: {
                std::uint64_t s = 0;
                for (std::size_t k = 0; k < count.size(); ++k) {
                    s += static_cast<std::uint64_t>(count.as_int64(k));
                }
                add_signed(static_cast<std::int64_t>(s));
                return;
            }
            case DType::UInt8:
            case DType::UInt16:
            case DType::UInt32:
            case DType::UInt64: {
                const std::vector<std::uint64_t> v = count.as<std::uint64_t>();
                std::uint64_t s = 0;
                for (const std::uint64_t x : v) {
                    s += x;
                }
                add_unsigned(s);
                return;
            }
            default: throw TypeError("cannot sum a non-numeric count column");
        }
    }

    [[nodiscard]] json::Value value() const {
        switch (kind_) {
            case Kind::PyInt: return json::Value::integer(0, DType::Int64);
            case Kind::Int64: return json::Value::integer(i_, DType::Int64);
            case Kind::UInt64: return json::Value::unsigned_integer(u_, DType::UInt64);
            case Kind::Float64: return json::Value::number(d_, DType::Float64);
        }
        return json::Value();
    }

  private:
    enum class Kind { PyInt, Int64, UInt64, Float64 };

    void add_float(double s) {
        switch (kind_) {
            case Kind::PyInt: d_ = 0.0 + s; break;
            case Kind::Int64: d_ = static_cast<double>(i_) + s; break;
            case Kind::UInt64: d_ = static_cast<double>(u_) + s; break;
            case Kind::Float64: d_ = d_ + s; break;
        }
        kind_ = Kind::Float64;
    }
    void add_signed(std::int64_t s) {
        switch (kind_) {
            case Kind::PyInt: i_ = s; kind_ = Kind::Int64; break;
            case Kind::Int64:
                i_ = static_cast<std::int64_t>(static_cast<std::uint64_t>(i_) + static_cast<std::uint64_t>(s));
                break;
            case Kind::UInt64: d_ = static_cast<double>(u_) + static_cast<double>(s); kind_ = Kind::Float64; break;
            case Kind::Float64: d_ = d_ + static_cast<double>(s); break;
        }
    }
    void add_unsigned(std::uint64_t s) {
        switch (kind_) {
            case Kind::PyInt: u_ = s; kind_ = Kind::UInt64; break;
            case Kind::Int64: d_ = static_cast<double>(i_) + static_cast<double>(s); kind_ = Kind::Float64; break;
            case Kind::UInt64: u_ += s; break;
            case Kind::Float64: d_ = d_ + static_cast<double>(s); break;
        }
    }

    Kind kind_ = Kind::PyInt;
    std::int64_t i_ = 0;
    std::uint64_t u_ = 0;
    double d_ = 0.0;
};

// ---------------------------------------------------------------------------
// index_bins / index_pixels: rlencode followed by the offset fill loop

class OffsetIndex {
  public:
    explicit OffsetIndex(std::size_t n_bins) : offsets_(n_bins + 1, 0) {}

    // Feeds the next values of the column in order.
    void feed(const std::int64_t* values, std::size_t n) {
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t v = values[k];
            if (!has_last_ || v != last_) {
                assign(position_ + static_cast<std::int64_t>(k), v);
            }
            has_last_ = true;
            last_ = v;
        }
        position_ += static_cast<std::int64_t>(n);
    }

    std::vector<std::int64_t> finish(std::int64_t total) {
        const npy::SliceRange tail = npy::slice_indices(
            current_, std::nullopt, static_cast<std::int64_t>(offsets_.size()));
        for (std::int64_t k = tail.start; k < tail.stop; ++k) {
            offsets_[static_cast<std::size_t>(k)] = total;
        }
        return std::move(offsets_);
    }

  private:
    void assign(std::int64_t start, std::int64_t value) {
        const npy::SliceRange range = npy::slice_indices(
            current_, value + 1, static_cast<std::int64_t>(offsets_.size()));
        for (std::int64_t k = range.start; k < range.stop; ++k) {
            offsets_[static_cast<std::size_t>(k)] = start;
        }
        current_ = value + 1;
    }

    std::vector<std::int64_t> offsets_;
    bool has_last_ = false;
    std::int64_t last_ = 0;
    std::int64_t position_ = 0;
    std::int64_t current_ = 0;
};

// ---------------------------------------------------------------------------
// create()

struct CreateCall {
    std::string cool_uri;
    std::optional<std::vector<std::string>> columns;
    std::vector<std::pair<std::string, DType>> dtypes;
    std::optional<json::Value> metadata;
    std::optional<std::string> assembly;
    bool symmetric_upper = true;
    std::string mode = "w";
    std::optional<H5Opts> h5opts;
    bool boundscheck = true;
    bool triucheck = true;
    bool dupcheck = true;
    bool ensure_sorted = false;
    std::optional<std::string> creation_date;
    std::string generated_by;
};

CreateCall make_call(const std::string& uri, const CreateOptions& o) {
    CreateCall call;
    call.cool_uri = uri;
    call.columns = o.columns;
    call.dtypes = o.dtypes;
    call.metadata = o.metadata;
    call.assembly = o.assembly;
    call.symmetric_upper = o.symmetric_upper;
    call.mode = o.mode;
    call.h5opts = o.h5opts;
    call.boundscheck = o.boundscheck;
    call.triucheck = o.triucheck;
    call.dupcheck = o.dupcheck;
    call.ensure_sorted = o.ensure_sorted;
    call.creation_date = o.creation_date;
    call.generated_by = o.generated_by;
    return call;
}

// The pixel input of create(): one borrowed table, or a generator.
struct Source {
    const Table* single = nullptr;
    PixelChunks chunks;
    std::optional<std::vector<std::string>> input_columns;
};

bool contains(const std::vector<std::string>& v, const std::string& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
}

void create_impl(const CreateCall& call, const Table& bins, Source source) {
    const auto [file_path, group_path] = parse_cooler_uri(call.cool_uri);
    const Kwargs h5opts = set_h5opts(call.h5opts);

    for (const char* col : {"chrom", "start", "end"}) {
        if (!bins.contains(col)) {
            throw ValueError(std::string("Missing column from bin table: '") + col + "'.");
        }
    }
    std::vector<std::string> columns;
    if (!call.columns.has_value()) {
        columns = {"bin1_id", "bin2_id", "count"};
    } else {
        columns = *call.columns;
        for (const char* col : {"bin1_id", "bin2_id"}) {
            if (!contains(columns, col)) {
                columns.insert(columns.begin(), col);
            }
        }
    }
    std::vector<std::pair<std::string, DType>> dtypes{
        {"bin1_id", DType::Int64}, {"bin2_id", DType::Int64}, {"count", DType::Int32}};
    for (const auto& [name, dtype] : call.dtypes) {
        bool replaced = false;
        for (auto& entry : dtypes) {
            if (entry.first == name) {
                entry.second = dtype;
                replaced = true;
            }
        }
        if (!replaced) {
            dtypes.emplace_back(name, dtype);
        }
    }
    const auto meta_dtype = [&dtypes](const std::string& name) {
        for (const auto& [n, d] : dtypes) {
            if (n == name) {
                return d;
            }
        }
        return DType::Float64;
    };
    if (source.input_columns.has_value()) {
        for (const std::string& col : columns) {
            if (!contains(*source.input_columns, col)) {
                const bool standard = col == "bin1_id" || col == "bin2_id" || col == "count";
                throw ValueError(std::string(standard ? "Standard" : "User") +
                                 " column not found in input: '" + col + "'");
            }
        }
    }

    const ChromSizes chromsizes = get_chromsizes(bins);
    if (chromsizes.size() == 0) {
        throw ValueError("not enough values to unpack (expected 2, got 0)");
    }
    const json::Value binsize = get_binsize(bins);
    const std::size_t n_chroms = chromsizes.size();
    const std::size_t n_bins = bins.num_rows();

    bool triucheck = call.triucheck;
    if (!call.symmetric_upper && triucheck) {
        warn("Creating a non-symmetric matrix, but `triucheck` was set to True. Changing to False.");
        triucheck = false;
    }
    const bool validate = call.boundscheck || triucheck || call.dupcheck || call.ensure_sorted;

    h5::File file(file_path, h5::parse_mode(call.mode));
    if (group_path == "/") {
        for (const char* name : {"chroms", "bins", "pixels", "indexes"}) {
            if (file.exists(std::string("/") + name)) {
                file.remove(std::string("/") + name);
            }
        }
    } else {
        if (file.exists(group_path)) {
            file.remove(group_path);
        }
        file.create_group(group_path);
    }

    // ---- chroms ----
    const std::string chroms_group = join_path(group_path, "chroms");
    file.create_group(chroms_group);
    {
        const std::vector<std::string>& names = chromsizes.names();
        std::size_t width = 1;
        for (const std::string& name : names) {
            width = std::max(width, name.size());
        }
        const h5::Handle type = h5::fixed_string_type(width);
        h5::Dataset dataset = create_dataset(file, join_path(chroms_group, "name"),
                                             TypeSpec{type.get(), width, 'S'}, n_chroms,
                                             std::nullopt, std::nullopt, h5opts);
        dataset.write_column(0, Column(names));
        const h5::Handle length_type = h5::file_type(DType::Int32);
        h5::Dataset lengths = create_dataset(file, join_path(chroms_group, "length"),
                                             numeric_spec(length_type, DType::Int32), n_chroms,
                                             std::nullopt, std::nullopt, h5opts);
        lengths.write_column(0, Column(chromsizes.lengths()).astype(DType::Int32));
    }

    // ---- bins ----
    const std::string bins_group = join_path(group_path, "bins");
    file.create_group(bins_group);
    std::vector<std::int32_t> chrom_ids(n_bins);
    {
        std::unordered_map<std::string, std::int32_t> idmap;
        for (std::size_t i = 0; i < n_chroms; ++i) {
            idmap[chromsizes.names()[i]] = static_cast<std::int32_t>(i);
        }
        const Column& chrom = bins["chrom"];
        if (chrom.dtype() == DType::Categorical) {
            const CategoricalData& data = chrom.categorical();
            std::vector<std::int32_t> by_code(data.categories->size(), -1);
            for (std::size_t c = 0; c < by_code.size(); ++c) {
                const auto found = idmap.find((*data.categories)[c]);
                if (found != idmap.end()) {
                    by_code[c] = found->second;
                }
            }
            for (std::size_t i = 0; i < n_bins; ++i) {
                const std::int32_t code = data.codes[i];
                if (code < 0) {
                    throw KeyError("nan");
                }
                chrom_ids[i] = by_code[static_cast<std::size_t>(code)];
            }
        } else {
            for (std::size_t i = 0; i < n_bins; ++i) {
                const std::string label = chrom.dtype() == DType::String
                                              ? chrom.values<std::string>()[i]
                                              : std::to_string(chrom.as_int64(i));
                chrom_ids[i] = idmap.at(label);
            }
        }
        const std::string chrom_path = join_path(bins_group, "chrom");
        bool as_enum = true;
        try {
            const h5::Handle type = h5::enum_type(chromsizes.names(), DType::Int32);
            h5::Dataset dataset = create_dataset(file, chrom_path, TypeSpec{type.get(), 4, 'e'},
                                                 n_bins, std::nullopt, std::nullopt, h5opts);
            const h5::Handle native(H5Tget_native_type(type.get(), H5T_DIR_ASCEND),
                                    h5::Handle::Kind::DataType);
            dataset.write_raw(0, n_bins, native.get(), chrom_ids.data());
        } catch (const ValueError&) {
            // Too many scaffolds for an HDF5 enum header: plain int32 IDs.
            as_enum = false;
            const h5::Handle type = h5::file_type(DType::Int32);
            h5::Dataset dataset = create_dataset(file, chrom_path, numeric_spec(type, DType::Int32),
                                                 n_bins, std::nullopt, std::nullopt, h5opts);
            dataset.write_raw(0, n_bins, H5T_NATIVE_INT32, chrom_ids.data());
        }
        if (!as_enum) {
            file.set_attribute(chrom_path, "enum_path", json::Value("/chroms/name"));
        }
        const h5::Handle coord_type = h5::file_type(DType::Int32);
        for (const char* name : {"start", "end"}) {
            h5::Dataset dataset = create_dataset(file, join_path(bins_group, name),
                                                 numeric_spec(coord_type, DType::Int32), n_bins,
                                                 std::nullopt, std::nullopt, h5opts);
            dataset.write_column(0, bins[name].astype(DType::Int32));
        }
        for (std::size_t c = 0; c < bins.num_columns(); ++c) {
            const std::string& name = bins.columns()[c];
            if (name == "chrom" || name == "start" || name == "end") {
                continue;
            }
            put_column(file, bins_group, name, bins.column(c));
        }
    }

    // ---- pixels ----
    const std::string pixels_group = join_path(group_path, "pixels");
    file.create_group(pixels_group);
    const auto nb = static_cast<std::uint64_t>(n_bins);
    const std::uint64_t max_size = call.symmetric_upper ? nb * (nb - (nb > 0 ? 1 : 0)) / 2 + nb
                                                        : nb * nb;
    const std::uint64_t init_size = std::min<std::uint64_t>(5 * nb, max_size);
    std::vector<std::string> prepared;
    const auto prepare = [&](const std::string& col) {
        const DType dtype = meta_dtype(col);
        const h5::Handle type = h5::file_type(dtype);
        create_dataset(file, join_path(pixels_group, col), numeric_spec(type, dtype),
                       static_cast<std::size_t>(init_size),
                       std::optional<std::size_t>(static_cast<std::size_t>(max_size)),
                       std::nullopt, h5opts);
        prepared.push_back(col);
    };
    prepare("bin1_id");
    prepare("bin2_id");
    if (contains(columns, "count")) {
        prepare("count");
    }
    {
        std::vector<std::string> rest = columns;
        for (const char* col : {"bin1_id", "bin2_id", "count"}) {
            const auto it = std::find(rest.begin(), rest.end(), col);
            if (it != rest.end()) {
                rest.erase(it);
            }
        }
        for (const std::string& col : rest) {
            prepare(col);
        }
    }

    std::vector<h5::Dataset> datasets;
    datasets.reserve(columns.size());
    for (const std::string& col : columns) {
        datasets.push_back(file.open_dataset(join_path(pixels_group, col)));
    }
    const DType bin1_file_dtype = meta_dtype("bin1_id");
    OffsetIndex bin1_index(n_bins);
    bool bin1_in_memory = true;

    std::size_t nnz = 0;
    PyTotal total;
    const auto write_chunk = [&](const Table& input) {
        std::optional<Table> reordered;
        if (validate) {
            reordered = validate_pixels(input, static_cast<std::int64_t>(n_bins), call.boundscheck,
                                        triucheck, call.dupcheck, call.ensure_sorted);
        }
        const Table& chunk = reordered.has_value() ? *reordered : input;
        if (!chunk.contains(columns[0])) {
            throw KeyError(columns[0]);
        }
        const std::size_t n = chunk[columns[0]].size();
        for (std::size_t k = 0; k < columns.size(); ++k) {
            datasets[k].resize(nnz + n);
            if (!chunk.contains(columns[k])) {
                throw KeyError(columns[k]);
            }
            const Column& values = chunk[columns[k]];
            datasets[k].write_column(nnz, values);
            if (columns[k] == "bin1_id" && bin1_in_memory) {
                if (values.dtype() == bin1_file_dtype && is_integer(values.dtype())) {
                    if (values.dtype() == DType::Int64) {
                        const auto& v = values.values<std::int64_t>();
                        bin1_index.feed(v.data(), v.size());
                    } else {
                        const std::vector<std::int64_t> v = values.as<std::int64_t>();
                        bin1_index.feed(v.data(), v.size());
                    }
                } else {
                    bin1_in_memory = false;
                }
            }
        }
        nnz += n;
        if (chunk.contains("count")) {
            total.add(chunk["count"]);
        }
        // cooler reopens and flushes the file for every chunk so that other
        // processes may read it in between; within one process a flush would
        // only evict the partially filled chunks from HDF5's chunk cache and
        // make every following append decompress and recompress them again.
    };
    if (source.single != nullptr) {
        write_chunk(*source.single);
    } else if (source.chunks) {
        while (std::optional<Table> chunk = source.chunks()) {
            write_chunk(*chunk);
        }
    }

    // ---- indexes ----
    const std::string indexes_group = join_path(group_path, "indexes");
    file.create_group(indexes_group);
    OffsetIndex chrom_index(n_chroms);
    {
        std::vector<std::int64_t> ids(chrom_ids.begin(), chrom_ids.end());
        chrom_index.feed(ids.data(), ids.size());
    }
    const std::vector<std::int64_t> chrom_offset = chrom_index.finish(static_cast<std::int64_t>(n_bins));
    std::vector<std::int64_t> bin1_offset;
    if (bin1_in_memory) {
        bin1_offset = bin1_index.finish(static_cast<std::int64_t>(nnz));
    } else {
        OffsetIndex readback(n_bins);
        const h5::Dataset bin1 = file.open_dataset(join_path(pixels_group, "bin1_id"));
        constexpr std::size_t kChunk = 1000000;
        for (std::size_t lo = 0; lo < nnz; lo += kChunk) {
            const std::vector<std::int64_t> v =
                bin1.read<std::int64_t>(lo, std::min(nnz, lo + kChunk));
            readback.feed(v.data(), v.size());
        }
        bin1_offset = readback.finish(static_cast<std::int64_t>(nnz));
    }
    const h5::Handle offset_type = h5::file_type(DType::Int64);
    {
        h5::Dataset dataset = create_dataset(file, join_path(indexes_group, "chrom_offset"),
                                             numeric_spec(offset_type, DType::Int64),
                                             chrom_offset.size(), std::nullopt, std::nullopt, h5opts);
        dataset.write_column(0, Column(chrom_offset));
    }
    {
        h5::Dataset dataset = create_dataset(file, join_path(indexes_group, "bin1_offset"),
                                             numeric_spec(offset_type, DType::Int64),
                                             bin1_offset.size(), std::nullopt, std::nullopt, h5opts);
        dataset.write_column(0, Column(bin1_offset));
    }

    // ---- info ----
    json::Value info = json::Value::object();
    info["bin-type"] = binsize.is_null() ? "variable" : "fixed";
    info["bin-size"] = binsize.is_null() ? json::Value("null") : binsize;
    info["storage-mode"] = call.symmetric_upper ? "symmetric-upper" : "square";
    info["nchroms"] = json::Value::integer(static_cast<std::int64_t>(n_chroms));
    info["nbins"] = json::Value::integer(static_cast<std::int64_t>(n_bins));
    info["sum"] = total.value();
    info["nnz"] = json::Value::integer(static_cast<std::int64_t>(nnz));
    if (call.assembly.has_value()) {
        info["genome-assembly"] = *call.assembly;
    }
    if (call.metadata.has_value()) {
        info["metadata"] = *call.metadata;
    }
    // write_info
    if (info.find("genome-assembly") == nullptr) {
        info["genome-assembly"] = "unknown";
    }
    const json::Value* metadata = info.find("metadata");
    info["metadata"] = json::dumps(metadata != nullptr ? *metadata : json::Value::object());
    info["creation-date"] = call.creation_date.has_value() ? *call.creation_date : npy::iso_now();
    info["generated-by"] = call.generated_by;
    info["format"] = kMagic;
    info["format-version"] = json::Value::integer(kFormatVersion);
    info["format-url"] = kFormatUrl;
    for (const auto& [key, value] : info.as_object()) {
        file.set_attribute(group_path, key, value);
    }
}

// ---------------------------------------------------------------------------
// merge_breakpoints and CoolerMerger (cooler/reduce.py)

std::vector<std::int64_t> merge_breakpoints(const std::vector<std::vector<std::int64_t>>& indexes,
                                            std::int64_t bufsize) {
    std::vector<double> combined(indexes.front().size(), 0.0);
    for (const auto& index : indexes) {
        for (std::size_t k = 0; k < combined.size(); ++k) {
            combined[k] += static_cast<double>(index[k]);
        }
    }
    double combined_start = 0;
    const double combined_nnz = combined.back();
    std::vector<std::int64_t> partition{0};
    std::vector<double> cum{0.0};
    std::size_t lo = 0;
    while (true) {
        const double target = std::min(combined_start + static_cast<double>(bufsize), combined_nnz);
        std::size_t hi =
            static_cast<std::size_t>(std::upper_bound(combined.begin() + static_cast<std::ptrdiff_t>(lo),
                                                      combined.end(), target) -
                                     combined.begin()) - 1;
        if (hi == lo) {
            hi += 1;
        }
        if (hi >= combined.size()) {
            throw IndexError("index " + std::to_string(hi) + " is out of bounds for axis 0 with size " +
                             std::to_string(combined.size()));
        }
        partition.push_back(static_cast<std::int64_t>(hi));
        cum.push_back(combined[hi]);
        if (combined[hi] == combined_nnz) {
            break;
        }
        lo = hi;
        combined_start = combined[hi];
    }
    std::int64_t n_over = 0;
    double largest = 0;
    for (std::size_t k = 1; k < cum.size(); ++k) {
        const double epoch = cum[k] - cum[k - 1];
        if (epoch > static_cast<double>(bufsize)) {
            ++n_over;
        }
        largest = k == 1 ? epoch : std::max(largest, epoch);
    }
    if (n_over > 0) {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%g", largest);
        warn(std::to_string(n_over) + " merge epochs will require buffering more than " +
             std::to_string(bufsize) + " pixel records, with as many as " + buffer + ".");
    }
    return partition;
}

// numpy.concatenate of two columns with dtype promotion.
void append_promoted(Column& into, const Column& more) {
    if (into.dtype() != more.dtype() && is_numeric(into.dtype()) && is_numeric(more.dtype())) {
        const DType common = result_type(into.dtype(), more.dtype());
        into = into.astype(common);
        into.append(more.astype(common));
        return;
    }
    into.append(more);
}

// GroupBy(["bin1_id", "bin2_id"], sort=True).sum() for one value column.
template <typename T>
Column group_sum_float(const std::vector<T>& values, const std::vector<std::size_t>& order,
                       const std::vector<std::size_t>& starts) {
    std::vector<T> out(starts.size() - 1);
    for (std::size_t g = 0; g + 1 < starts.size(); ++g) {
        T sum = 0;
        T compensation = 0;
        for (std::size_t p = starts[g]; p < starts[g + 1]; ++p) {
            const T val = values[order[p]];
            if (std::isnan(val)) {
                continue;
            }
            // pandas group_sum's Kahan summation.
            const T y = val - compensation;
            const T t = sum + y;
            compensation = t - sum - y;
            if (compensation != compensation) {
                compensation = 0;
            }
            sum = t;
        }
        out[g] = sum;
    }
    return Column(std::move(out));
}

Column group_sum(const Column& values, const std::vector<std::size_t>& order,
                 const std::vector<std::size_t>& starts) {
    const DType dtype = values.dtype();
    if (dtype == DType::Float64) {
        return group_sum_float(values.values<double>(), order, starts);
    }
    if (dtype == DType::Float32) {
        return group_sum_float(values.values<float>(), order, starts);
    }
    if (is_integer(dtype) || dtype == DType::Bool) {
        const bool is_unsigned = is_unsigned_integer(dtype);
        std::vector<std::uint64_t> sums(starts.size() - 1, 0);
        for (std::size_t g = 0; g + 1 < starts.size(); ++g) {
            for (std::size_t p = starts[g]; p < starts[g + 1]; ++p) {
                sums[g] += is_unsigned ? values.as<std::uint64_t>().at(order[p])
                                       : static_cast<std::uint64_t>(values.as_int64(order[p]));
            }
        }
        const Column sum_column(sums);
        return dtype == DType::Bool ? sum_column.astype(DType::Int64) : sum_column.astype(dtype);
    }
    throw TypeError("cannot sum a non-numeric column");
}

class CoolerMerger {
  public:
    CoolerMerger(std::vector<Cooler> coolers, std::int64_t mergebuf,
                 const std::optional<std::vector<std::string>>& columns)
        : coolers_(std::move(coolers)), mergebuf_(mergebuf),
          columns_(columns.has_value() ? *columns : std::vector<std::string>{"count"}) {
        if (coolers_.empty()) {
            throw IndexError("list index out of range");
        }
        const std::optional<std::int64_t> binsize = coolers_[0].binsize();
        if (binsize.has_value()) {
            std::set<std::optional<std::int64_t>> sizes;
            for (const Cooler& c : coolers_) {
                sizes.insert(c.binsize());
            }
            if (sizes.size() > 1) {
                throw ValueError("Coolers must have the same resolution");
            }
            const ChromSizes& reference = coolers_[0].chromsizes();
            for (std::size_t i = 1; i < coolers_.size(); ++i) {
                const ChromSizes& other = coolers_[i].chromsizes();
                if (other.names() != reference.names()) {
                    throw ValueError("Can only compare identically-labeled Series objects");
                }
                if (other.lengths() != reference.lengths()) {
                    throw ValueError("Coolers must have the same chromosomes");
                }
            }
        } else {
            const Table reference = coolers_[0].bins()[Fields({"chrom", "start", "end"})].all();
            for (std::size_t i = 1; i < coolers_.size(); ++i) {
                const Table other = coolers_[i].bins()[Fields({"chrom", "start", "end"})].all();
                if (other.num_rows() != reference.num_rows() ||
                    other["chrom"].labels() != reference["chrom"].labels() ||
                    other["start"].as<std::int64_t>() != reference["start"].as<std::int64_t>() ||
                    other["end"].as<std::int64_t>() != reference["end"].as<std::int64_t>()) {
                    throw ValueError("Coolers must have same bin structure");
                }
            }
        }
    }

    PixelChunks chunks() {
        auto state = std::make_shared<State>();
        for (const Cooler& c : coolers_) {
            const h5::File file(c.filename(), h5::Mode::Read);
            const h5::Dataset ds = file.open_dataset(join_path(c.root(), "indexes/bin1_offset"));
            state->indexes.push_back(ds.read<std::int64_t>(0, ds.length()));
        }
        state->partition = merge_breakpoints(state->indexes, mergebuf_);
        for (const Cooler& c : coolers_) {
            (void)c.pixels().size();
        }
        state->starts.assign(coolers_.size(), 0);
        state->epoch = 1;
        auto coolers = coolers_;
        auto columns = columns_;
        return [state, coolers, columns]() -> std::optional<Table> {
            if (state->epoch >= state->partition.size()) {
                return std::nullopt;
            }
            const auto bin1_id = static_cast<std::size_t>(state->partition[state->epoch]);
            std::vector<std::int64_t> stops;
            for (const auto& index : state->indexes) {
                stops.push_back(index.at(bin1_id));
            }
            std::optional<Table> combined;
            for (std::size_t i = 0; i < coolers.size(); ++i) {
                if (stops[i] - state->starts[i] <= 0) {
                    continue;
                }
                Table part = coolers[i].pixels().slice(state->starts[i], stops[i]);
                if (!combined.has_value()) {
                    combined = std::move(part);
                } else {
                    Table merged;
                    for (std::size_t c = 0; c < combined->num_columns(); ++c) {
                        const std::string& name = combined->columns()[c];
                        Column column = combined->column(c);
                        append_promoted(column, part[name]);
                        merged.set(name, std::move(column));
                    }
                    combined = std::move(merged);
                }
            }
            if (!combined.has_value()) {
                throw ValueError("No objects to concatenate");
            }
            const Ids b1((*combined)["bin1_id"]);
            const Ids b2((*combined)["bin2_id"]);
            const bool floating = b1.floating() || b2.floating();
            std::vector<std::size_t> order(b1.size());
            std::iota(order.begin(), order.end(), std::size_t{0});
            std::stable_sort(order.begin(), order.end(), PairLess{b1, b2, floating});
            std::vector<std::size_t> starts;
            std::vector<std::int64_t> first_rows;
            for (std::size_t p = 0; p < order.size(); ++p) {
                const bool new_group =
                    p == 0 || PairLess{b1, b2, floating}(order[p - 1], order[p]);
                if (new_group) {
                    starts.push_back(p);
                    first_rows.push_back(static_cast<std::int64_t>(order[p]));
                }
            }
            starts.push_back(order.size());
            Table out;
            out.set("bin1_id", (*combined)["bin1_id"].take(first_rows));
            out.set("bin2_id", (*combined)["bin2_id"].take(first_rows));
            for (const std::string& name : columns) {
                if (!combined->contains(name)) {
                    throw KeyError("Column(s) ['" + name + "'] do not exist");
                }
                out.set(name, group_sum((*combined)[name], order, starts));
            }
            out.reset_index();
            state->starts = std::move(stops);
            ++state->epoch;
            return out;
        };
    }

  private:
    struct State {
        std::vector<std::vector<std::int64_t>> indexes;
        std::vector<std::int64_t> partition;
        std::vector<std::int64_t> starts;
        std::size_t epoch = 1;
    };

    std::vector<Cooler> coolers_;
    std::int64_t mergebuf_;
    std::vector<std::string> columns_;
};

// tempfile.NamedTemporaryFile(suffix=".multi.cool", dir=temp_dir)
std::string make_temp_file(const std::string& directory) {
    std::string pattern = (directory.empty() ? std::string() : directory + "/") + "tmpXXXXXX.multi.cool";
    std::vector<char> buffer(pattern.begin(), pattern.end());
    buffer.push_back('\0');
    const int fd = mkstemps(buffer.data(), 11);
    if (fd < 0) {
        throw OSError("[Errno " + std::to_string(errno) + "] " + std::strerror(errno) + ": '" +
                      pattern + "'");
    }
    close(fd);
    return std::string(buffer.data());
}

struct TempFiles {
    std::vector<std::string> paths;
    bool remove = true;
    TempFiles() = default;
    TempFiles(const TempFiles&) = delete;
    TempFiles& operator=(const TempFiles&) = delete;
    ~TempFiles() {
        if (remove) {
            for (const std::string& path : paths) {
                ::unlink(path.c_str());
            }
        }
    }
};

void create_from_unordered(const std::string& cool_uri, const Table& bins, PixelChunks chunks,
                           const CreateOptions& options) {
    std::optional<std::vector<std::string>> columns;
    if (options.columns.has_value()) {
        columns.emplace();
        for (const std::string& col : *options.columns) {
            if (col != "bin1_id" && col != "bin2_id") {
                columns->push_back(col);
            }
        }
    }
    std::string directory;
    if (!options.temp_dir.has_value()) {
        const std::string file = parse_cooler_uri(cool_uri).first;
        const std::size_t slash = file.rfind('/');
        directory = slash == std::string::npos ? std::string() : (slash == 0 ? "/" : file.substr(0, slash));
    } else if (*options.temp_dir == "-") {
        const char* env = nullptr;
        for (const char* name : {"TMPDIR", "TEMP", "TMP"}) {
            env = std::getenv(name);
            if (env != nullptr && *env != '\0') {
                break;
            }
        }
        directory = (env != nullptr && *env != '\0') ? env : "/tmp";
    } else {
        directory = *options.temp_dir;
    }

    TempFiles temps;
    temps.remove = options.delete_temp;
    CreateCall sub = make_call(cool_uri, options);
    sub.columns = columns;
    sub.mode = "a";

    const std::string temp1 = make_temp_file(directory);
    temps.paths.push_back(temp1);
    std::vector<std::string> uris;
    for (std::size_t i = 0;; ++i) {
        std::optional<Table> chunk = chunks();
        if (!chunk.has_value()) {
            break;
        }
        const std::string uri = temp1 + "::" + std::to_string(i);
        uris.push_back(uri);
        sub.cool_uri = uri;
        Source source;
        source.single = &*chunk;
        source.input_columns = chunk->columns();
        create_impl(sub, bins, source);
    }

    const auto open_all = [](const std::vector<std::string>& list) {
        std::vector<Cooler> coolers;
        coolers.reserve(list.size());
        for (const std::string& uri : list) {
            coolers.emplace_back(uri);
        }
        return coolers;
    };

    std::vector<std::string> final_uris = uris;
    const auto n = static_cast<std::int64_t>(uris.size());
    if (n > options.max_merge && options.max_merge > 0) {
        const auto k = static_cast<std::int64_t>(std::sqrt(static_cast<double>(n)));
        const std::vector<std::int64_t> edges = npy::linspace_int(0, n, k);
        const std::string temp2 = make_temp_file(directory);
        temps.paths.push_back(temp2);
        std::vector<std::string> uris2;
        for (std::size_t e = 1; e < edges.size(); ++e) {
            const std::int64_t lo = edges[e - 1];
            const std::int64_t hi = edges[e];
            const npy::SliceRange range = npy::slice_indices(lo, hi, n);
            CoolerMerger merger(
                open_all(std::vector<std::string>(uris.begin() + range.start, uris.begin() + range.stop)),
                options.mergebuf, columns);
            const std::string uri = temp2 + "::" + std::to_string(lo) + "-" + std::to_string(hi);
            uris2.push_back(uri);
            sub.cool_uri = uri;
            Source source;
            source.chunks = merger.chunks();
            create_impl(sub, bins, source);
        }
        final_uris = uris2;
    }

    CoolerMerger merger(open_all(final_uris), options.mergebuf, columns);
    CreateCall final_call = make_call(cool_uri, options);
    final_call.columns = columns;
    Source source;
    source.chunks = merger.chunks();
    create_impl(final_call, bins, source);
}

}  // namespace

void create_cooler(const std::string& cool_uri, const Table& bins, const Table& pixels,
                   const CreateOptions& options) {
    for (const char* key : {"bin1_id", "bin2_id"}) {
        if (!pixels.contains(key)) {
            throw KeyError(key);
        }
    }
    const std::optional<Table> sorted = sort_by_pixel_ids(pixels);
    Source source;
    source.single = sorted.has_value() ? &*sorted : &pixels;
    source.input_columns = pixels.columns();
    create_impl(make_call(cool_uri, options), bins, source);
}

void create_cooler(const std::string& cool_uri, const Table& bins, PixelChunks pixels,
                   const CreateOptions& options) {
    if (options.ordered) {
        Source source;
        source.chunks = std::move(pixels);
        create_impl(make_call(cool_uri, options), bins, source);
        return;
    }
    create_from_unordered(cool_uri, bins, std::move(pixels), options);
}

}  // namespace coolercpp
