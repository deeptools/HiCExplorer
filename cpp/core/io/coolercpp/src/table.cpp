#include "coolercpp/table.hpp"

#include <algorithm>

namespace coolercpp {

Column Column::empty(DType dtype, std::size_t size) {
    switch (dtype) {
        case DType::Bool: return Column(std::vector<Bool8>(size));
        case DType::Int8: return Column(std::vector<std::int8_t>(size));
        case DType::Int16: return Column(std::vector<std::int16_t>(size));
        case DType::Int32: return Column(std::vector<std::int32_t>(size));
        case DType::Int64: return Column(std::vector<std::int64_t>(size));
        case DType::UInt8: return Column(std::vector<std::uint8_t>(size));
        case DType::UInt16: return Column(std::vector<std::uint16_t>(size));
        case DType::UInt32: return Column(std::vector<std::uint32_t>(size));
        case DType::UInt64: return Column(std::vector<std::uint64_t>(size));
        case DType::Float32: return Column(std::vector<float>(size));
        case DType::Float64: return Column(std::vector<double>(size));
        case DType::String: return Column(std::vector<std::string>(size));
        case DType::Categorical: {
            CategoricalData data;
            data.codes.assign(size, 0);
            data.categories = std::make_shared<const std::vector<std::string>>();
            return Column(std::move(data));
        }
    }
    throw TypeError("unknown dtype");
}

Column Column::categorical(std::vector<std::int32_t> codes,
                           std::vector<std::string> categories, bool ordered) {
    return categorical(std::move(codes),
                       std::make_shared<const std::vector<std::string>>(std::move(categories)),
                       ordered);
}

Column Column::categorical(std::vector<std::int32_t> codes,
                           std::shared_ptr<const std::vector<std::string>> categories,
                           bool ordered) {
    CategoricalData data;
    data.codes = std::move(codes);
    data.categories = std::move(categories);
    data.ordered = ordered;
    return Column(std::move(data));
}

std::size_t Column::size() const noexcept {
    return std::visit(
        [](const auto& v) -> std::size_t {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, CategoricalData>) {
                return v.codes.size();
            } else {
                return v.size();
            }
        },
        data_);
}

const CategoricalData& Column::categorical() const {
    if (const auto* c = std::get_if<CategoricalData>(&data_)) {
        return *c;
    }
    throw TypeError(type_mismatch(DType::Categorical));
}

std::string Column::type_mismatch(DType requested) const {
    return "column has dtype " + std::string(dtype_name(dtype())) + ", not " +
           std::string(dtype_name(requested));
}

Column Column::astype(DType target) const {
    if (target == dtype()) {
        return *this;
    }
    switch (target) {
        case DType::Bool: return Column(as<Bool8>());
        case DType::Int8: return Column(as<std::int8_t>());
        case DType::Int16: return Column(as<std::int16_t>());
        case DType::Int32: return Column(as<std::int32_t>());
        case DType::Int64: return Column(as<std::int64_t>());
        case DType::UInt8: return Column(as<std::uint8_t>());
        case DType::UInt16: return Column(as<std::uint16_t>());
        case DType::UInt32: return Column(as<std::uint32_t>());
        case DType::UInt64: return Column(as<std::uint64_t>());
        case DType::Float32: return Column(as<float>());
        case DType::Float64: return Column(as<double>());
        case DType::String:
            if (dtype() == DType::Categorical) {
                return Column(labels());
            }
            break;
        case DType::Categorical: break;
    }
    throw TypeError("cannot convert a " + std::string(dtype_name(dtype())) + " column to " +
                    std::string(dtype_name(target)));
}

double Column::as_double(std::size_t i) const {
    return std::visit(
        [i](const auto& v) -> double {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, std::vector<std::string>>) {
                throw TypeError("a string column has no numeric value");
            } else if constexpr (std::is_same_v<V, CategoricalData>) {
                return static_cast<double>(v.codes.at(i));
            } else if constexpr (std::is_same_v<V, std::vector<Bool8>>) {
                return v.at(i).value;
            } else {
                return static_cast<double>(v.at(i));
            }
        },
        data_);
}

std::int64_t Column::as_int64(std::size_t i) const {
    return std::visit(
        [i](const auto& v) -> std::int64_t {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, std::vector<std::string>>) {
                throw TypeError("a string column has no numeric value");
            } else if constexpr (std::is_same_v<V, CategoricalData>) {
                return v.codes.at(i);
            } else if constexpr (std::is_same_v<V, std::vector<Bool8>>) {
                return v.at(i).value;
            } else {
                return static_cast<std::int64_t>(v.at(i));
            }
        },
        data_);
}

std::string Column::label(std::size_t i) const {
    if (const auto* s = std::get_if<std::vector<std::string>>(&data_)) {
        return s->at(i);
    }
    if (const auto* c = std::get_if<CategoricalData>(&data_)) {
        const std::int32_t code = c->codes.at(i);
        if (code < 0 || c->categories == nullptr ||
            static_cast<std::size_t>(code) >= c->categories->size()) {
            return {};
        }
        return (*c->categories)[static_cast<std::size_t>(code)];
    }
    throw TypeError("column of dtype " + std::string(dtype_name(dtype())) + " has no labels");
}

std::vector<std::string> Column::labels() const {
    if (const auto* s = std::get_if<std::vector<std::string>>(&data_)) {
        return *s;
    }
    std::vector<std::string> out(size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = label(i);
    }
    return out;
}

Column Column::slice(std::size_t lo, std::size_t hi) const {
    if (lo > hi || hi > size()) {
        throw IndexError("column slice [" + std::to_string(lo) + ", " + std::to_string(hi) +
                         ") out of range for size " + std::to_string(size()));
    }
    return std::visit(
        [lo, hi](const auto& v) -> Column {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, CategoricalData>) {
                CategoricalData out;
                out.codes.assign(v.codes.begin() + static_cast<std::ptrdiff_t>(lo),
                                 v.codes.begin() + static_cast<std::ptrdiff_t>(hi));
                out.categories = v.categories;
                out.ordered = v.ordered;
                return Column(std::move(out));
            } else {
                return Column(V(v.begin() + static_cast<std::ptrdiff_t>(lo),
                                v.begin() + static_cast<std::ptrdiff_t>(hi)));
            }
        },
        data_);
}

Column Column::take(std::span<const std::int64_t> positions) const {
    const std::size_t n = size();
    for (const std::int64_t p : positions) {
        if (p < 0 || static_cast<std::size_t>(p) >= n) {
            throw IndexError("positional indexer " + std::to_string(p) + " is out-of-bounds");
        }
    }
    return std::visit(
        [positions](const auto& v) -> Column {
            using V = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<V, CategoricalData>) {
                CategoricalData out;
                out.codes.reserve(positions.size());
                for (const std::int64_t p : positions) {
                    out.codes.push_back(v.codes[static_cast<std::size_t>(p)]);
                }
                out.categories = v.categories;
                out.ordered = v.ordered;
                return Column(std::move(out));
            } else {
                V out;
                out.reserve(positions.size());
                for (const std::int64_t p : positions) {
                    out.push_back(v[static_cast<std::size_t>(p)]);
                }
                return Column(std::move(out));
            }
        },
        data_);
}

void Column::append(const Column& other) {
    if (other.dtype() != dtype()) {
        throw TypeError("cannot append a " + std::string(dtype_name(other.dtype())) +
                        " column to a " + std::string(dtype_name(dtype())) + " column");
    }
    std::visit(
        [&other](auto& v) {
            using V = std::decay_t<decltype(v)>;
            const auto& o = std::get<V>(other.data_);
            if constexpr (std::is_same_v<V, CategoricalData>) {
                if (v.categories != o.categories &&
                    (v.categories == nullptr || o.categories == nullptr ||
                     *v.categories != *o.categories)) {
                    if (v.codes.empty() && v.categories != nullptr && v.categories->empty()) {
                        v.categories = o.categories;
                        v.ordered = o.ordered;
                    } else {
                        throw ValueError("cannot concatenate categoricals with different "
                                         "categories");
                    }
                }
                v.codes.insert(v.codes.end(), o.codes.begin(), o.codes.end());
            } else /* NOLINT */ {
                v.insert(v.end(), o.begin(), o.end());
            }
        },
        data_);
}

bool operator==(const CategoricalData& x, const CategoricalData& y) {
    const bool same_categories =
        x.categories == y.categories ||
        (x.categories != nullptr && y.categories != nullptr && *x.categories == *y.categories);
    return same_categories && x.ordered == y.ordered && x.codes == y.codes;
}

bool operator==(const Column& a, const Column& b) {
    return a.data_ == b.data_;
}

Index Index::range(std::int64_t start, std::size_t size) {
    Index out;
    out.start_ = start;
    out.size_ = size;
    return out;
}

Index Index::values(std::vector<std::int64_t> values) {
    Index out;
    out.size_ = values.size();
    out.values_ = std::move(values);
    return out;
}

std::vector<std::int64_t> Index::to_vector() const {
    if (values_.has_value()) {
        return *values_;
    }
    std::vector<std::int64_t> out(size_);
    for (std::size_t i = 0; i < size_; ++i) {
        out[i] = start_ + static_cast<std::int64_t>(i);
    }
    return out;
}

bool operator==(const Index& a, const Index& b) {
    if (a.size() != b.size()) {
        return false;
    }
    if (a.is_range() && b.is_range()) {
        return a.size() == 0 || a.start() == b.start();
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

Table::Table(std::initializer_list<std::pair<std::string, Column>> columns) {
    for (const auto& [name, column] : columns) {
        set(name, column);
    }
}

bool Table::contains(std::string_view name) const noexcept {
    return std::find(names_.begin(), names_.end(), name) != names_.end();
}

const Column& Table::operator[](std::string_view name) const {
    const auto it = std::find(names_.begin(), names_.end(), name);
    if (it == names_.end()) {
        throw KeyError(std::string(name));
    }
    return columns_[static_cast<std::size_t>(it - names_.begin())];
}

Column& Table::operator[](std::string_view name) {
    const auto it = std::find(names_.begin(), names_.end(), name);
    if (it == names_.end()) {
        throw KeyError(std::string(name));
    }
    return columns_[static_cast<std::size_t>(it - names_.begin())];
}

Table& Table::set(std::string name, Column column) {
    if (names_.empty() && index_.size() == 0) {
        index_ = Index::range(0, column.size());
    } else if (column.size() != num_rows()) {
        throw ValueError("Length of values (" + std::to_string(column.size()) +
                         ") does not match length of index (" +
                         std::to_string(num_rows()) + ")");
    }
    const auto it = std::find(names_.begin(), names_.end(), name);
    if (it != names_.end()) {
        columns_[static_cast<std::size_t>(it - names_.begin())] = std::move(column);
    } else {
        names_.push_back(std::move(name));
        columns_.push_back(std::move(column));
    }
    return *this;
}

Table& Table::insert(std::size_t position, std::string name, Column column) {
    if (contains(name)) {
        throw ValueError("cannot insert " + name + ", already exists");
    }
    if (names_.empty() && index_.size() == 0) {
        index_ = Index::range(0, column.size());
    } else if (column.size() != num_rows()) {
        throw ValueError("Length of values (" + std::to_string(column.size()) +
                         ") does not match length of index (" +
                         std::to_string(num_rows()) + ")");
    }
    position = std::min(position, names_.size());
    names_.insert(names_.begin() + static_cast<std::ptrdiff_t>(position), std::move(name));
    columns_.insert(columns_.begin() + static_cast<std::ptrdiff_t>(position), std::move(column));
    return *this;
}

Table& Table::drop(std::string_view name) {
    const auto it = std::find(names_.begin(), names_.end(), name);
    if (it == names_.end()) {
        throw KeyError("\"['" + std::string(name) + "'] not found in axis\"");
    }
    const auto position = it - names_.begin();
    names_.erase(it);
    columns_.erase(columns_.begin() + position);
    return *this;
}

Table Table::select(const std::vector<std::string>& names) const {
    Table out;
    out.index_ = index_;
    for (const std::string& name : names) {
        const auto it = std::find(names_.begin(), names_.end(), name);
        if (it == names_.end()) {
            throw KeyError("\"['" + name + "'] not in index\"");
        }
        out.names_.push_back(name);
        out.columns_.push_back(columns_[static_cast<std::size_t>(it - names_.begin())]);
    }
    return out;
}

Table Table::slice(std::size_t lo, std::size_t hi) const {
    if (lo > hi || hi > num_rows()) {
        throw IndexError("row slice out of range");
    }
    Table out;
    out.names_ = names_;
    out.columns_.reserve(columns_.size());
    for (const Column& column : columns_) {
        out.columns_.push_back(column.slice(lo, hi));
    }
    if (index_.is_range()) {
        out.index_ = Index::range(index_.start() + static_cast<std::int64_t>(lo), hi - lo);
    } else {
        const std::vector<std::int64_t> all = index_.to_vector();
        out.index_ = Index::values(
            std::vector<std::int64_t>(all.begin() + static_cast<std::ptrdiff_t>(lo),
                                      all.begin() + static_cast<std::ptrdiff_t>(hi)));
    }
    return out;
}

Table Table::take(std::span<const std::int64_t> positions) const {
    Table out;
    out.names_ = names_;
    out.columns_.reserve(columns_.size());
    for (const Column& column : columns_) {
        out.columns_.push_back(column.take(positions));
    }
    std::vector<std::int64_t> index(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const std::int64_t p = positions[i];
        if (p < 0 || static_cast<std::size_t>(p) >= index_.size()) {
            throw IndexError("positional indexer is out-of-bounds");
        }
        index[i] = index_[static_cast<std::size_t>(p)];
    }
    out.index_ = Index::values(std::move(index));
    return out;
}

void Table::set_index(Index index) {
    if (!names_.empty() && index.size() != num_rows()) {
        throw ValueError("Length mismatch: Expected axis has " + std::to_string(num_rows()) +
                         " elements, new values have " + std::to_string(index.size()) +
                         " elements");
    }
    index_ = std::move(index);
}

std::vector<std::pair<std::string, DType>> Table::dtypes() const {
    std::vector<std::pair<std::string, DType>> out;
    out.reserve(names_.size());
    for (std::size_t i = 0; i < names_.size(); ++i) {
        out.emplace_back(names_[i], columns_[i].dtype());
    }
    return out;
}

bool operator==(const Table& a, const Table& b) {
    return a.names_ == b.names_ && a.columns_ == b.columns_ && a.index_ == b.index_;
}

}  // namespace coolercpp
