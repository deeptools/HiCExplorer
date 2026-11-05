// Port of cooler/api.py, cooler/core/_selectors.py, cooler/core/_tableops.py
// (get) and cooler/fileops.py (list_coolers) of cooler 0.10.2
// (BSD-3-Clause).

#include "coolercpp/api.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>

#include "coolercpp/errors.hpp"
#include "coolercpp/util.hpp"
#include "h5.hpp"
#include "numpy_compat.hpp"
#include "rangequery.hpp"

namespace coolercpp {

namespace {

std::string join_path(const std::string& group, const std::string& name) {
    return group == "/" ? "/" + name : group + "/" + name;
}

// repr() of a Python str (ASCII content).
std::string py_repr(const std::string& s) {
    const bool has_single = s.find('\'') != std::string::npos;
    const bool has_double = s.find('"') != std::string::npos;
    const char quote = (has_single && !has_double) ? '"' : '\'';
    std::string out(1, quote);
    for (const char ch : s) {
        const auto c = static_cast<unsigned char>(ch);
        if (ch == quote || ch == '\\') {
            out.push_back('\\');
            out.push_back(ch);
        } else if (ch == '\n') {
            out += "\\n";
        } else if (ch == '\r') {
            out += "\\r";
        } else if (ch == '\t') {
            out += "\\t";
        } else if (c < 0x20 || c == 0x7F) {
            char buffer[8];
            std::snprintf(buffer, sizeof(buffer), "\\x%02x", c);
            out += buffer;
        } else {
            out.push_back(ch);
        }
    }
    out.push_back(quote);
    return out;
}

std::string py_list_repr(const std::vector<std::string>& items) {
    std::string out = "[";
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += py_repr(items[i]);
    }
    return out + "]";
}

// cooler.api.info
json::Value read_info(const h5::File& file, const std::string& path) {
    json::Value out = json::Value::object();
    for (auto& [key, value] : file.attributes(path)) {
        if (value.is_string()) {
            std::optional<json::Value> parsed = json::parse(value.as_string());
            out[key] = parsed.has_value() ? std::move(*parsed) : value;
        } else {
            out[key] = value;
        }
    }
    return out;
}

// _IndexingMixin._process_slice
std::pair<std::int64_t, std::int64_t> process_slice(const Slice& s, std::int64_t nmax) {
    if (s.step.has_value() && *s.step != 1) {
        throw ValueError("slicing with step != 1 not supported");
    }
    std::int64_t i0 = 0;
    if (s.start.has_value()) {
        i0 = *s.start < 0 ? nmax + *s.start : *s.start;
    }
    std::int64_t i1 = nmax;
    if (s.stop.has_value()) {
        i1 = *s.stop < 0 ? nmax + *s.stop : *s.stop;
    }
    return {i0, i1};
}

std::pair<std::int64_t, std::int64_t> process_index(std::int64_t s, std::int64_t nmax) {
    if (s < 0) {
        s += nmax;
    }
    if (s >= nmax) {
        throw IndexError("index is out of bounds");
    }
    return {s, s + 1};
}

std::pair<std::int64_t, std::int64_t> process_key(const AxisKey& key, std::int64_t nmax) {
    if (const auto* slice = std::get_if<Slice>(&key.value())) {
        return process_slice(*slice, nmax);
    }
    return process_index(std::get<std::int64_t>(key.value()), nmax);
}

Column categorical_from_codes(const Column& codes, const std::vector<std::string>& categories) {
    const std::set<std::string> unique(categories.begin(), categories.end());
    if (unique.size() != categories.size()) {
        throw ValueError("Categorical categories must be unique");
    }
    if (!codes.empty() && !is_integer(codes.dtype()) && codes.dtype() != DType::Categorical) {
        throw ValueError("codes need to be array-like integers");
    }
    const auto n = static_cast<std::int64_t>(categories.size());
    std::vector<std::int32_t> out(codes.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::int64_t code = codes.as_int64(i);
        if (code < -1 || code >= n) {
            throw ValueError("codes need to be between -1 and len(categories)-1");
        }
        out[i] = static_cast<std::int32_t>(code);
    }
    return Column::categorical(std::move(out), categories, true);
}

// One column of cooler.core.get.
Column read_field(const h5::File& file, const std::string& path, std::optional<std::int64_t> lo,
                  std::optional<std::int64_t> hi, bool convert_enum, const std::string& name) {
    if (!file.is_dataset(path)) {
        throw KeyError("Unable to synchronously open object (object '" + name +
                       "' doesn't exist)");
    }
    const h5::Dataset dataset = file.open_dataset(path);
    const h5::TypeDesc type = dataset.type();
    const npy::SliceRange range =
        npy::slice_indices(lo, hi, static_cast<std::int64_t>(dataset.length()));
    if (convert_enum && type.kind == h5::TypeDesc::Kind::Enum) {
        const Column codes = dataset.read_column(static_cast<std::size_t>(range.start),
                                                 static_cast<std::size_t>(range.stop));
        auto members = type.enum_members;
        std::stable_sort(members.begin(), members.end(),
                         [](const auto& a, const auto& b) { return a.second < b.second; });
        std::vector<std::string> categories;
        categories.reserve(members.size());
        for (const auto& member : members) {
            categories.push_back(member.first);
        }
        return categorical_from_codes(codes, categories);
    }
    return dataset.read_column(static_cast<std::size_t>(range.start),
                               static_cast<std::size_t>(range.stop));
}

// cooler.core.get(grp, lo, hi, fields, convert_enum)
Table get_table(const h5::File& file, const std::string& group, std::optional<std::int64_t> lo,
                std::optional<std::int64_t> hi, const std::vector<std::string>& fields,
                bool convert_enum) {
    if (!file.is_group(group)) {
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    Table table;
    for (const std::string& field : fields) {
        table.set(field, read_field(file, join_path(group, field), lo, hi, convert_enum, field));
    }
    if (!fields.empty() && lo.has_value()) {
        table.set_index(Index::range(*lo, table.num_rows()));
    }
    return table;
}

std::vector<std::string> default_fields(const h5::File& file, const std::string& group,
                                        std::initializer_list<const char*> leading) {
    if (!file.is_group(group)) {
        throw KeyError("Unable to synchronously open object (component not found)");
    }
    std::vector<std::string> fields(leading.begin(), leading.end());
    for (const std::string& key : file.keys(group)) {
        if (std::find(fields.begin(), fields.end(), key) == fields.end()) {
            fields.push_back(key);
        }
    }
    return fields;
}

std::vector<std::string> field_names(const Fields& fields, const h5::File& file,
                                     const std::string& group,
                                     std::initializer_list<const char*> leading) {
    if (fields.is_none()) {
        return default_fields(file, group, leading);
    }
    if (fields.is_single()) {
        return {fields.name()};
    }
    return fields.names();
}

// ---- element-wise arithmetic with numpy's promotion ----

template <typename T>
T wrapping_multiply(T a, T b) {
    if constexpr (std::is_integral_v<T>) {
        using U = std::make_unsigned_t<T>;
        return static_cast<T>(static_cast<U>(a) * static_cast<U>(b));
    } else {
        return a * b;
    }
}

template <typename T>
Column multiply_as(const Column& a, const Column& b) {
    const std::vector<T> x = a.as<T>();
    const std::vector<T> y = b.as<T>();
    std::vector<T> out(x.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = wrapping_multiply(x[i], y[i]);
    }
    return Column(std::move(out));
}

Column multiply(const Column& a, const Column& b) {
    if (a.size() != b.size()) {
        throw ValueError("operands could not be broadcast together with shapes (" +
                         std::to_string(a.size()) + ",) (" + std::to_string(b.size()) + ",) ");
    }
    switch (result_type(a.dtype(), b.dtype())) {
        case DType::Bool: {
            const std::vector<Bool8> x = a.as<Bool8>();
            const std::vector<Bool8> y = b.as<Bool8>();
            std::vector<Bool8> out(x.size());
            for (std::size_t i = 0; i < out.size(); ++i) {
                out[i] = Bool8(x[i] && y[i]);
            }
            return Column(std::move(out));
        }
        case DType::Int8: return multiply_as<std::int8_t>(a, b);
        case DType::Int16: return multiply_as<std::int16_t>(a, b);
        case DType::Int32: return multiply_as<std::int32_t>(a, b);
        case DType::Int64: return multiply_as<std::int64_t>(a, b);
        case DType::UInt8: return multiply_as<std::uint8_t>(a, b);
        case DType::UInt16: return multiply_as<std::uint16_t>(a, b);
        case DType::UInt32: return multiply_as<std::uint32_t>(a, b);
        case DType::UInt64: return multiply_as<std::uint64_t>(a, b);
        case DType::Float32: return multiply_as<float>(a, b);
        case DType::Float64: return multiply_as<double>(a, b);
        default: break;
    }
    throw TypeError("unsupported operand types for *");
}

// 1 / array
Column reciprocal(const Column& w) {
    if (w.dtype() == DType::Float32) {
        const auto& v = w.values<float>();
        std::vector<float> out(v.size());
        for (std::size_t i = 0; i < v.size(); ++i) {
            out[i] = 1.0F / v[i];
        }
        return Column(std::move(out));
    }
    const std::vector<double> v = w.as<double>();
    std::vector<double> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        out[i] = 1.0 / v[i];
    }
    return Column(std::move(out));
}

// array[indices] with numpy's wrap and bounds rules.
Column take_wrapped(const Column& values, const Column& indices) {
    const auto n = static_cast<std::int64_t>(values.size());
    std::vector<std::int64_t> positions(indices.size());
    for (std::size_t k = 0; k < positions.size(); ++k) {
        positions[k] = npy::wrap_index(indices.as_int64(k), n);
    }
    return values.take(positions);
}

Column read_bias(const h5::Dataset& weights, std::int64_t lo, std::int64_t hi) {
    const npy::SliceRange range =
        npy::slice_indices(lo, hi, static_cast<std::int64_t>(weights.length()));
    return weights.read_column(static_cast<std::size_t>(range.start),
                               static_cast<std::size_t>(range.stop));
}

std::string shape_text(std::int64_t a, std::int64_t b) {
    return "(" + std::to_string(a) + "," + std::to_string(b) + ")";
}

// arr * np.outer(bias1, bias2)
DenseMatrix apply_outer(const DenseMatrix& arr, const Column& bias1, const Column& bias2) {
    const auto b1 = static_cast<std::int64_t>(bias1.size());
    const auto b2 = static_cast<std::int64_t>(bias2.size());
    const auto broadcast = [&](std::int64_t a, std::int64_t b) -> std::int64_t {
        if (a == b || b == 1) {
            return a;
        }
        if (a == 1) {
            return b;
        }
        throw ValueError("operands could not be broadcast together with shapes " +
                         shape_text(arr.rows, arr.cols) + " " + shape_text(b1, b2) + " ");
    };
    const std::int64_t rows = broadcast(arr.rows, b1);
    const std::int64_t cols = broadcast(arr.cols, b2);
    const DType outer_dtype = result_type(bias1.dtype(), bias2.dtype());
    const DType out_dtype = result_type(arr.dtype(), outer_dtype);
    DenseMatrix out;
    out.rows = rows;
    out.cols = cols;
    const auto cell = [&](std::int64_t i, std::int64_t j) -> std::size_t {
        const std::int64_t ai = arr.rows == 1 ? 0 : i;
        const std::int64_t aj = arr.cols == 1 ? 0 : j;
        return static_cast<std::size_t>(ai * arr.cols + aj);
    };
    const auto total = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (out_dtype == DType::Float64 || out_dtype == DType::Float32) {
        const std::vector<double> x1 = bias1.as<double>();
        const std::vector<double> x2 = bias2.as<double>();
        if (out_dtype == DType::Float64) {
            const std::vector<double> a = arr.values.as<double>();
            std::vector<double> values(total);
            for (std::int64_t i = 0; i < rows; ++i) {
                for (std::int64_t j = 0; j < cols; ++j) {
                    const std::size_t i1 = static_cast<std::size_t>(b1 == 1 ? 0 : i);
                    const std::size_t j2 = static_cast<std::size_t>(b2 == 1 ? 0 : j);
                    double outer = 0.0;
                    if (outer_dtype == DType::Float32) {
                        outer = static_cast<double>(static_cast<float>(x1[i1]) *
                                                    static_cast<float>(x2[j2]));
                    } else {
                        outer = x1[i1] * x2[j2];
                    }
                    values[static_cast<std::size_t>(i * cols + j)] = a[cell(i, j)] * outer;
                }
            }
            out.values = Column(std::move(values));
        } else {
            const std::vector<float> a = arr.values.as<float>();
            std::vector<float> values(total);
            for (std::int64_t i = 0; i < rows; ++i) {
                for (std::int64_t j = 0; j < cols; ++j) {
                    const std::size_t i1 = static_cast<std::size_t>(b1 == 1 ? 0 : i);
                    const std::size_t j2 = static_cast<std::size_t>(b2 == 1 ? 0 : j);
                    const float outer = static_cast<float>(x1[i1]) * static_cast<float>(x2[j2]);
                    values[static_cast<std::size_t>(i * cols + j)] = a[cell(i, j)] * outer;
                }
            }
            out.values = Column(std::move(values));
        }
        return out;
    }
    throw TypeError("integer balancing weights are not supported");
}

// cooler.api.annotate for either kind of bin table.
Table annotate_impl(const Table& pixels, std::size_t bins_len,
                    const std::function<Table(std::int64_t, std::optional<std::int64_t>)>& loc_slice,
                    bool replace) {
    std::vector<Table> annotations;
    for (const auto& [column, suffix] :
         {std::pair<const char*, const char*>{"bin1_id", "1"}, {"bin2_id", "2"}}) {
        if (!pixels.contains(column)) {
            continue;
        }
        const Column& ids_column = pixels[column];
        const DType dtype = ids_column.dtype();
        const bool safe = dtype == DType::Bool || is_signed_integer(dtype) ||
                          dtype == DType::UInt8 || dtype == DType::UInt16 ||
                          dtype == DType::UInt32;
        if (!safe) {
            throw TypeError("Cannot cast array data from dtype('" +
                            std::string(dtype_name(dtype)) +
                            "') to dtype('int64') according to the rule 'safe'");
        }
        const std::vector<std::int64_t> ids = ids_column.as<std::int64_t>();
        std::int64_t bmin = 0;
        std::optional<std::int64_t> bmax = 0;
        if (!ids.empty()) {
            if (bins_len > pixels.num_rows()) {
                const auto [mn, mx] = std::minmax_element(ids.begin(), ids.end());
                bmin = *mn;
                bmax = *mx;
            } else {
                bmin = 0;
                bmax = std::nullopt;
            }
        }
        const Table ann = loc_slice(bmin, bmax);
        if (ann.num_rows() == 0) {
            throw IndexError("index 0 is out of bounds for axis 0 with size 0");
        }
        const std::int64_t first = ann.index()[0];
        const auto n = static_cast<std::int64_t>(ann.num_rows());
        std::vector<std::int64_t> positions(ids.size());
        for (std::size_t k = 0; k < ids.size(); ++k) {
            std::int64_t p = ids[k] - first;
            if (p < 0) {
                p += n;
            }
            if (p < 0 || p >= n) {
                throw IndexError("positional indexers are out-of-bounds");
            }
            positions[k] = p;
        }
        const Table taken = ann.take(positions);
        Table renamed;
        for (std::size_t c = 0; c < taken.num_columns(); ++c) {
            renamed.set(taken.columns()[c] + suffix, taken.column(c));
        }
        renamed.reset_index();
        annotations.push_back(std::move(renamed));
    }
    Table out;
    for (const Table& ann : annotations) {
        for (std::size_t c = 0; c < ann.num_columns(); ++c) {
            out.set(ann.columns()[c], ann.column(c));
        }
    }
    for (std::size_t c = 0; c < pixels.num_columns(); ++c) {
        const std::string& name = pixels.columns()[c];
        if (replace && (name == "bin1_id" || name == "bin2_id")) {
            continue;
        }
        out.set(name, pixels.column(c));
    }
    out.set_index(pixels.index());
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Selectors

RangeSelector1D::RangeSelector1D(Fields fields, Slicer slicer, Fetcher fetcher, std::int64_t nmax)
    : fields_(std::move(fields)), slicer_(std::move(slicer)), fetcher_(std::move(fetcher)),
      nmax_(nmax) {}

std::vector<std::string> RangeSelector1D::columns() const {
    const Table empty = slicer_(fields_, 0, 0);
    if (fields_.is_single()) {
        throw AttributeError("'Series' object has no attribute 'columns'");
    }
    return empty.columns();
}

std::vector<std::pair<std::string, DType>> RangeSelector1D::dtypes() const {
    return slicer_(fields_, 0, 0).dtypes();
}

std::vector<std::string> RangeSelector1D::keys() const { return columns(); }

bool RangeSelector1D::contains(const std::string& key) const {
    const std::vector<std::string> names = columns();
    return std::find(names.begin(), names.end(), key) != names.end();
}

RangeSelector1D RangeSelector1D::operator[](Fields fields) const {
    return RangeSelector1D(std::move(fields), slicer_, fetcher_, nmax_);
}

Table RangeSelector1D::operator[](const Slice& rows) const {
    const auto [lo, hi] = process_slice(rows, nmax_);
    return slicer_(fields_, lo, hi);
}

Table RangeSelector1D::row_at(std::int64_t row) const {
    const auto [lo, hi] = process_index(row, nmax_);
    return slicer_(fields_, lo, hi);
}

Table RangeSelector1D::fetch(const Region& region) const {
    if (!fetcher_) {
        throw NotImplementedError("");
    }
    const auto [lo, hi] = fetcher_(region);
    return slicer_(fields_, lo, hi);
}

RangeSelector2D::RangeSelector2D(std::optional<std::string> field, Slicer slicer, Fetcher fetcher,
                                 std::array<std::int64_t, 2> shape)
    : field_(std::move(field)), slicer_(std::move(slicer)), fetcher_(std::move(fetcher)),
      shape_(shape) {}

MatrixResult RangeSelector2D::operator()(const AxisKey& rows, const AxisKey& cols) const {
    const auto [i0, i1] = process_key(rows, shape_[0]);
    const auto [j0, j1] = process_key(cols, shape_[1]);
    return slicer_(field_, i0, i1, j0, j1);
}

MatrixResult RangeSelector2D::fetch(const Region& region,
                                    const std::optional<Region>& region2) const {
    if (!fetcher_) {
        throw NotImplementedError("");
    }
    const auto box = fetcher_(region, region2);
    return slicer_(field_, box[0], box[1], box[2], box[3]);
}

const DenseMatrix& MatrixResult::dense() const {
    if (const auto* v = std::get_if<DenseMatrix>(&value_)) {
        return *v;
    }
    throw TypeError("the matrix query did not return a dense array");
}
const SparseMatrix& MatrixResult::sparse() const {
    if (const auto* v = std::get_if<SparseMatrix>(&value_)) {
        return *v;
    }
    throw TypeError("the matrix query did not return a sparse matrix");
}
const Table& MatrixResult::pixels() const {
    if (const auto* v = std::get_if<Table>(&value_)) {
        return *v;
    }
    throw TypeError("the matrix query did not return pixels");
}
DenseMatrix& MatrixResult::dense() {
    return const_cast<DenseMatrix&>(static_cast<const MatrixResult&>(*this).dense());
}
SparseMatrix& MatrixResult::sparse() {
    return const_cast<SparseMatrix&>(static_cast<const MatrixResult&>(*this).sparse());
}
Table& MatrixResult::pixels() {
    return const_cast<Table&>(static_cast<const MatrixResult&>(*this).pixels());
}

// ---------------------------------------------------------------------------
// Cooler

Cooler::Cooler(const std::string& store) {
    auto [file, group] = parse_cooler_uri(store);
    filename_ = std::move(file);
    root_ = std::move(group);
    uri_ = filename_ + "::" + root_;
    store_ = filename_;
    refresh();
}

Cooler::Cooler(const std::string& store, const std::string& root) {
    if (!h5::is_hdf5(store)) {
        throw ValueError("Not a valid path to a Cooler file");
    }
    filename_ = store;
    root_ = root;
    uri_ = filename_ + "::" + root_;
    store_ = filename_;
    refresh();
}

void Cooler::refresh() {
    const h5::File file(store_, h5::Mode::Read);
    try {
        if (!file.exists(root_)) {
            throw KeyError(root_);
        }
        const std::string group = join_path(root_, "chroms");
        const Table table = get_table(file, group, 0, std::nullopt,
                                      default_fields(file, group, {"name", "length"}), true);
        std::vector<std::string> names = table["name"].labels();
        const Column& length = table["length"];
        std::vector<std::int64_t> lengths = length.as<std::int64_t>();
        chromids_.clear();
        for (std::size_t i = 0; i < names.size(); ++i) {
            chromids_[names[i]] = static_cast<std::int64_t>(i);
        }
        chromsizes_ = ChromSizes(std::move(names), std::move(lengths), length.dtype());
        info_ = read_info(file, root_);
        const json::Value* mode = info_.find("storage-mode");
        is_symm_upper_ =
            mode == nullptr || (mode->is_string() && mode->as_string() == "symmetric-upper");
    } catch (const KeyError&) {
        std::string message = "No cooler found at: " + store_ + ".";
        const std::vector<std::string> listing = list_coolers(store_);
        if (!listing.empty()) {
            message += " Coolers found in " + py_list_repr(listing) +
                       ". Use '::' to specify a group path";
        }
        throw KeyError(message);
    }
}

const json::Value& Cooler::info_item(const std::string& key) const {
    const json::Value* value = info_.find(key);
    if (value == nullptr) {
        throw KeyError(key);
    }
    return *value;
}

std::string Cooler::storage_mode() const {
    const json::Value* mode = info_.find("storage-mode");
    if (mode == nullptr) {
        return "symmetric-upper";
    }
    return mode->is_string() ? mode->as_string() : json::dumps(*mode);
}

std::optional<std::int64_t> Cooler::binsize() const {
    const json::Value& value = info_item("bin-size");
    if (value.is_null()) {
        return std::nullopt;
    }
    if (value.type() == json::Type::Double) {
        return static_cast<std::int64_t>(value.as_double());
    }
    return value.as_int();
}

json::Value Cooler::info() const {
    const h5::File file(store_, h5::Mode::Read);
    return read_info(file, root_);
}

std::array<std::int64_t, 2> Cooler::shape() const {
    const std::int64_t n = info_item("nbins").as_int();
    return {n, n};
}

namespace {

std::pair<std::int64_t, std::int64_t> extent_of(
    const h5::File& file, const std::string& root,
    const std::unordered_map<std::string, std::int64_t>& chromids, const RegionTuple& region,
    std::optional<std::int64_t> binsize) {
    const auto found = chromids.find(region.chrom);
    if (found == chromids.end()) {
        throw KeyError(region.chrom);
    }
    const std::int64_t cid = found->second;
    const h5::Dataset chrom_offset = file.open_dataset(join_path(root, "indexes/chrom_offset"));
    const auto n = static_cast<std::int64_t>(chrom_offset.length());
    const auto scalar = [&](std::int64_t i) {
        const std::int64_t k = npy::wrap_index(i, n);
        return chrom_offset.read<std::int64_t>(static_cast<std::size_t>(k),
                                               static_cast<std::size_t>(k + 1))[0];
    };
    if (binsize.has_value()) {
        if (*binsize == 0) {
            throw Error("division by zero");
        }
        const std::int64_t offset = scalar(cid);
        const double bs = static_cast<double>(*binsize);
        return {offset + static_cast<std::int64_t>(std::floor(static_cast<double>(region.start) / bs)),
                offset + static_cast<std::int64_t>(std::ceil(static_cast<double>(region.end) / bs))};
    }
    const std::int64_t lo = scalar(cid);
    const std::int64_t hi = scalar(cid + 1);
    const h5::Dataset starts_ds = file.open_dataset(join_path(root, "bins/start"));
    const npy::SliceRange range =
        npy::slice_indices(lo, hi, static_cast<std::int64_t>(starts_ds.length()));
    const std::vector<std::int64_t> starts = starts_ds.read<std::int64_t>(
        static_cast<std::size_t>(range.start), static_cast<std::size_t>(range.stop));
    const auto first = static_cast<std::int64_t>(
        std::upper_bound(starts.begin(), starts.end(), region.start) - starts.begin());
    const auto second = static_cast<std::int64_t>(
        std::lower_bound(starts.begin(), starts.end(), region.end) - starts.begin());
    return {lo + first - 1, lo + second};
}

}  // namespace

std::pair<std::int64_t, std::int64_t> Cooler::region_extent(const Region& region) const {
    const RegionTuple parsed = parse_region(region, &chromsizes_);
    const std::optional<std::int64_t> bs = binsize();
    const h5::File file(store_, h5::Mode::Read);
    return extent_of(file, root_, chromids_, parsed, bs);
}

std::int64_t Cooler::offset(const Region& region) const { return region_extent(region).first; }

std::pair<std::int64_t, std::int64_t> Cooler::extent(const Region& region) const {
    return region_extent(region);
}

RangeSelector1D Cooler::chroms(bool convert_enum) const {
    const std::int64_t nmax = info_item("nchroms").as_int();
    const Cooler self = *this;
    auto slicer = [self, convert_enum](const Fields& fields, std::int64_t lo, std::int64_t hi) {
        const h5::File file(self.store_, h5::Mode::Read);
        const std::string group = join_path(self.root_, "chroms");
        return get_table(file, group, lo, hi, field_names(fields, file, group, {"name", "length"}),
                         convert_enum);
    };
    return RangeSelector1D(Fields(), slicer, nullptr, nmax);
}

RangeSelector1D Cooler::bins(bool convert_enum) const {
    const std::int64_t nmax = info_item("nbins").as_int();
    const Cooler self = *this;
    auto slicer = [self, convert_enum](const Fields& fields, std::int64_t lo, std::int64_t hi) {
        const h5::File file(self.store_, h5::Mode::Read);
        const std::string group = join_path(self.root_, "bins");
        const std::vector<std::string> names =
            field_names(fields, file, group, {"chrom", "start", "end"});
        Table out = get_table(file, group, lo, hi, names, convert_enum);
        // `"chrom" in fields` is a substring test when fields is a str.
        const bool has_chrom =
            fields.is_single()
                ? fields.name().find("chrom") != std::string::npos
                : std::find(names.begin(), names.end(), "chrom") != names.end();
        if (has_chrom) {
            const std::string column = fields.is_single() ? fields.name() : "chrom";
            const Column& chrom = out[column];
            if (is_integer(chrom.dtype()) && convert_enum) {
                const std::string chroms_group = join_path(self.root_, "chroms");
                const std::vector<std::string> chromnames =
                    get_table(file, chroms_group, 0, std::nullopt, {"name"}, true)["name"].labels();
                out.set(column, categorical_from_codes(chrom, chromnames));
            }
        }
        return out;
    };
    auto fetcher = [self](const Region& region) { return self.region_extent(region); };
    return RangeSelector1D(Fields(), slicer, fetcher, nmax);
}

RangeSelector1D Cooler::pixels(bool join, bool convert_enum) const {
    const std::int64_t nmax = info_item("nnz").as_int();
    const Cooler self = *this;
    auto slicer = [self, join, convert_enum](const Fields& fields, std::int64_t lo,
                                             std::int64_t hi) {
        const h5::File file(self.store_, h5::Mode::Read);
        const std::string group = join_path(self.root_, "pixels");
        Table df = get_table(file, group, lo, hi,
                             field_names(fields, file, group, {"bin1_id", "bin2_id"}), convert_enum);
        if (join) {
            if (fields.is_single()) {
                throw AttributeError("'Series' object has no attribute 'columns'");
            }
            const Table bins = get_table(file, join_path(self.root_, "bins"), 0, std::nullopt,
                                         {"chrom", "start", "end"}, convert_enum);
            df = annotate(df, bins, true);
        }
        return df;
    };
    auto fetcher = [self](const Region& region) {
        const auto [i0, i1] = self.region_extent(region);
        const h5::File file(self.store_, h5::Mode::Read);
        const h5::Dataset offsets = file.open_dataset(join_path(self.root_, "indexes/bin1_offset"));
        const auto n = static_cast<std::int64_t>(offsets.length());
        const auto scalar = [&](std::int64_t i) {
            const std::int64_t k = npy::wrap_index(i, n);
            return offsets.read<std::int64_t>(static_cast<std::size_t>(k),
                                              static_cast<std::size_t>(k + 1))[0];
        };
        const std::int64_t lo = scalar(i0);
        const std::int64_t hi = scalar(i1);
        return std::pair<std::int64_t, std::int64_t>{lo, hi};
    };
    return RangeSelector1D(Fields(), slicer, fetcher, nmax);
}

RangeSelector2D Cooler::matrix(const MatrixOptions& options) const {
    MatrixOptions opts = options;
    if (opts.balance.is_column_name() && !opts.divisive_weights.has_value()) {
        const std::string name = opts.balance.column();
        // cooler's _4DN_DIVISIVE_WEIGHTS is {"KR", "VC", "VC_SQRT"}.
        if (name == "KR" || name == "VC" || name == "VC_SQRT") {
            opts.divisive_weights = true;
        }
    }
    const Cooler self = *this;
    auto slicer = [self, opts](const std::optional<std::string>& field, std::int64_t i0,
                               std::int64_t i1, std::int64_t j0, std::int64_t j1) {
        return self.query_matrix(field.value_or("count"), i0, i1, j0, j1, opts);
    };
    auto fetcher = [self](const Region& region, const std::optional<Region>& region2) {
        const RegionTuple r1 = parse_region(region, &self.chromsizes_);
        const RegionTuple r2 = parse_region(region2.has_value() ? *region2 : region,
                                            &self.chromsizes_);
        const h5::File file(self.store_, h5::Mode::Read);
        const std::optional<std::int64_t> bs = self.binsize();
        const auto [i0, i1] = extent_of(file, self.root_, self.chromids_, r1, bs);
        const auto [j0, j1] = extent_of(file, self.root_, self.chromids_, r2, bs);
        return std::array<std::int64_t, 4>{i0, i1, j0, j1};
    };
    const std::int64_t n = info_item("nbins").as_int();
    return RangeSelector2D(opts.field, slicer, fetcher, {n, n});
}

MatrixResult Cooler::query_matrix(const std::string& field, std::int64_t i0, std::int64_t i1,
                                  std::int64_t j0, std::int64_t j1,
                                  const MatrixOptions& options) const {
    const h5::File file(store_, h5::Mode::Read);
    const std::string bins_group = join_path(root_, "bins");
    const bool balance = options.balance.enabled();
    const std::string name = options.balance.column();
    if (balance && !file.exists(join_path(bins_group, name))) {
        throw ValueError("No column 'bins/" + name +
                         "'found. Use ``cooler.balance_cooler`` to calculate balancing weights "
                         "or set balance=False.");
    }
    const h5::Dataset offsets_ds = file.open_dataset(join_path(root_, "indexes/bin1_offset"));
    const detail::CSRReader reader(file, join_path(root_, "pixels"),
                                   offsets_ds.read<std::int64_t>(0, offsets_ds.length()));
    const detail::Bbox bbox{i0, i1, j0, j1};
    const bool divisive = options.divisive_weights.value_or(false);

    if (options.as_pixels) {
        Table df = detail::frame_from_dict(
            detail::direct_query(reader, field, bbox, options.chunksize, !options.ignore_index),
            field);
        if (balance) {
            Cooler fresh = *this;
            fresh.refresh();
            const Table df2 =
                annotate(df, fresh.bins()[Fields(std::vector<std::string>{name})], false);
            Column w1 = df2[name + "1"];
            Column w2 = df2[name + "2"];
            if (divisive) {
                w1 = reciprocal(w1);
                w2 = reciprocal(w2);
            }
            df.set("balanced", multiply(multiply(w1, w2), df2[field]));
        }
        if (options.join) {
            Cooler fresh = *this;
            fresh.refresh();
            df = annotate(df, fresh.bins()[Fields({"chrom", "start", "end"})], true);
        }
        return df;
    }

    detail::PixelDict dict =
        is_symm_upper_ ? detail::fill_lower_query(reader, field, bbox, options.chunksize, false)
                       : detail::direct_query(reader, field, bbox, options.chunksize, false);

    const auto biases = [&]() {
        const h5::Dataset weights = file.open_dataset(join_path(bins_group, name));
        Column bias1 = read_bias(weights, i0, i1);
        Column bias2 = (i0 == j0 && i1 == j1) ? bias1 : read_bias(weights, j0, j1);
        if (divisive) {
            bias1 = reciprocal(bias1);
            bias2 = reciprocal(bias2);
        }
        return std::pair<Column, Column>{std::move(bias1), std::move(bias2)};
    };

    if (options.sparse) {
        SparseMatrix mat = detail::sparse_from_dict(std::move(dict), bbox);
        if (balance) {
            const auto [bias1, bias2] = biases();
            mat.data = multiply(multiply(take_wrapped(bias1, mat.row), take_wrapped(bias2, mat.col)),
                                mat.data);
        }
        return mat;
    }
    DenseMatrix arr = detail::dense_from_dict(dict, bbox);
    if (balance) {
        const auto [bias1, bias2] = biases();
        arr = apply_outer(arr, bias1, bias2);
    }
    return arr;
}

// ---------------------------------------------------------------------------
// annotate, list_coolers

Table annotate(const Table& pixels, const Table& bins, bool replace) {
    const auto loc_slice = [&bins](std::int64_t beg, std::optional<std::int64_t> end) {
        const Index& index = bins.index();
        const std::size_t n = index.size();
        std::vector<std::int64_t> labels;
        if (!index.is_range()) {
            labels = index.to_vector();
            if (!std::is_sorted(labels.begin(), labels.end())) {
                throw KeyError("label slicing requires a monotonic index");
            }
        }
        const auto position = [&](std::int64_t label, bool upper) -> std::size_t {
            if (index.is_range()) {
                std::int64_t p = label - index.start() + (upper ? 1 : 0);
                p = std::clamp<std::int64_t>(p, 0, static_cast<std::int64_t>(n));
                return static_cast<std::size_t>(p);
            }
            const auto it = upper ? std::upper_bound(labels.begin(), labels.end(), label)
                                  : std::lower_bound(labels.begin(), labels.end(), label);
            return static_cast<std::size_t>(it - labels.begin());
        };
        const std::size_t lo = position(beg, false);
        std::size_t hi = end.has_value() ? position(*end, true) : n;
        hi = std::max(hi, lo);
        return bins.slice(lo, hi);
    };
    return annotate_impl(pixels, bins.num_rows(), loc_slice, replace);
}

Table annotate(const Table& pixels, const RangeSelector1D& bins, bool replace) {
    const auto loc_slice = [&bins](std::int64_t beg, std::optional<std::int64_t> end) {
        std::optional<std::int64_t> stop;
        if (end.has_value()) {
            stop = *end + 1;
        }
        return bins[Slice{beg, stop, std::nullopt}];
    };
    return annotate_impl(pixels, bins.size(), loc_slice, replace);
}

std::vector<std::string> list_coolers(const std::string& filepath) {
    if (!h5::is_hdf5(filepath)) {
        throw OSError("'" + filepath + "' is not an HDF5 file.");
    }
    const h5::File file(filepath, h5::Mode::Read);
    std::vector<std::string> listing;
    const auto check = [&](const std::string& path) {
        for (const auto& [key, value] : file.attributes(path)) {
            if (key == "format" && value.is_string() && value.as_string() == "HDF5::Cooler") {
                bool complete = file.is_group(path);
                for (const char* name : {"chroms", "bins", "pixels", "indexes"}) {
                    complete = complete && file.exists(join_path(path, name));
                }
                if (!complete) {
                    warn("Cooler path " + path + " appears to be corrupt");
                }
                listing.push_back(path);
            }
        }
    };
    std::function<void(const std::string&)> visit = [&](const std::string& group) {
        for (const std::string& name : file.keys(group)) {
            const std::string child = join_path(group, name);
            check(child);
            if (file.is_group(child)) {
                visit(child);
            }
        }
    };
    check("/");
    visit("/");
    return natsorted(std::move(listing));
}

}  // namespace coolercpp
