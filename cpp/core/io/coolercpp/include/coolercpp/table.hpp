// A typed columnar table: the C++ counterpart of the pandas DataFrames and
// Series that cooler's Python API returns and accepts.
//
// A Column is one contiguous typed buffer (numpy semantics) or a list of
// strings or a Categorical. A Table is an ordered list of named columns plus
// an index. Tables carry arbitrary extra columns, which is how cooler's bin
// table holds "weight", "KR", "VC" or any user column, and how a pixel table
// holds value columns besides "count".
//
// Mapping to pandas:
//   DataFrame            -> Table
//   Series (one column)  -> Table with a single column (the Series name)
//   df.index             -> Table::index() (RangeIndex or explicit int64)
//   df["x"]              -> table["x"] (a Column)
//   df.dtypes            -> Table::dtypes()
//   pd.Categorical       -> Column with DType::Categorical

#ifndef COOLERCPP_TABLE_HPP
#define COOLERCPP_TABLE_HPP

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "coolercpp/dtype.hpp"
#include "coolercpp/errors.hpp"

namespace coolercpp {

struct CategoricalData {
    std::vector<std::int32_t> codes;
    // Shared, because every slice of a bin table refers to the same labels.
    std::shared_ptr<const std::vector<std::string>> categories;
    bool ordered = true;
};

// Same codes, same labels (by content) and the same ordered flag.
bool operator==(const CategoricalData& a, const CategoricalData& b);

// The element types a Column can hold.
template <typename T>
concept ColumnElement =
    std::is_same_v<T, Bool8> || std::is_same_v<T, std::int8_t> ||
    std::is_same_v<T, std::int16_t> || std::is_same_v<T, std::int32_t> ||
    std::is_same_v<T, std::int64_t> || std::is_same_v<T, std::uint8_t> ||
    std::is_same_v<T, std::uint16_t> || std::is_same_v<T, std::uint32_t> ||
    std::is_same_v<T, std::uint64_t> || std::is_same_v<T, float> ||
    std::is_same_v<T, double> || std::is_same_v<T, std::string>;

// Alternative i holds the element type of DType i.
using ColumnData =
    std::variant<std::vector<Bool8>, std::vector<std::int8_t>, std::vector<std::int16_t>,
                 std::vector<std::int32_t>, std::vector<std::int64_t>,
                 std::vector<std::uint8_t>, std::vector<std::uint16_t>,
                 std::vector<std::uint32_t>, std::vector<std::uint64_t>, std::vector<float>,
                 std::vector<double>, std::vector<std::string>, CategoricalData>;

class Column {
  public:
    // An empty float64 column.
    Column() : data_(std::vector<double>{}) {}

    template <ColumnElement T>
    explicit Column(std::vector<T> values) : data_(std::move(values)) {}

    template <ColumnElement T>
    Column(std::initializer_list<T> values) : data_(std::vector<T>(values)) {}

    explicit Column(CategoricalData categorical) : data_(std::move(categorical)) {}

    explicit Column(ColumnData data) : data_(std::move(data)) {}

    // An empty column of the given dtype (Categorical gets no labels).
    [[nodiscard]] static Column empty(DType dtype, std::size_t size = 0);

    [[nodiscard]] static Column categorical(std::vector<std::int32_t> codes,
                                            std::vector<std::string> categories,
                                            bool ordered = true);
    [[nodiscard]] static Column categorical(
        std::vector<std::int32_t> codes,
        std::shared_ptr<const std::vector<std::string>> categories, bool ordered = true);

    [[nodiscard]] DType dtype() const noexcept { return static_cast<DType>(data_.index()); }
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }

    [[nodiscard]] const ColumnData& data() const noexcept { return data_; }
    [[nodiscard]] ColumnData& data() noexcept { return data_; }

    // The underlying buffer. Throws TypeError when T is not the column's type.
    template <typename T>
    [[nodiscard]] const std::vector<T>& values() const {
        if (const auto* v = std::get_if<std::vector<T>>(&data_)) {
            return *v;
        }
        throw TypeError(type_mismatch(dtype_v<T>));
    }
    template <typename T>
    [[nodiscard]] std::vector<T>& values() {
        if (auto* v = std::get_if<std::vector<T>>(&data_)) {
            return *v;
        }
        throw TypeError(type_mismatch(dtype_v<T>));
    }

    [[nodiscard]] const CategoricalData& categorical() const;

    // numpy astype for numeric columns (C casts, float to int truncates);
    // a Categorical converts its codes. Throws TypeError for strings.
    template <typename T>
    [[nodiscard]] std::vector<T> as() const {
        std::vector<T> out(size());
        std::visit(
            [&](const auto& v) {
                using V = std::decay_t<decltype(v)>;
                if constexpr (std::is_same_v<V, std::vector<std::string>>) {
                    throw TypeError("cannot convert a string column to a number");
                } else if constexpr (std::is_same_v<V, CategoricalData>) {
                    for (std::size_t i = 0; i < v.codes.size(); ++i) {
                        out[i] = static_cast<T>(v.codes[i]);
                    }
                } else {
                    for (std::size_t i = 0; i < v.size(); ++i) {
                        out[i] = cast_element<T>(v[i]);
                    }
                }
            },
            data_);
        return out;
    }

    // A converted copy as a Column of the target dtype. String and
    // Categorical targets are only reachable from the same kind.
    [[nodiscard]] Column astype(DType target) const;

    [[nodiscard]] double as_double(std::size_t i) const;
    [[nodiscard]] std::int64_t as_int64(std::size_t i) const;
    // The label of a String or Categorical element (a NaN code, -1, gives "").
    [[nodiscard]] std::string label(std::size_t i) const;
    // Every element as a label (String or Categorical columns only).
    [[nodiscard]] std::vector<std::string> labels() const;

    [[nodiscard]] Column slice(std::size_t lo, std::size_t hi) const;
    // Elements at the given positions (numpy fancy indexing, no bounds wrap).
    [[nodiscard]] Column take(std::span<const std::int64_t> positions) const;
    // Appends another column of the same dtype (numpy.concatenate).
    void append(const Column& other);

    friend bool operator==(const Column& a, const Column& b);

  private:
    template <typename T, typename U>
    static T cast_element(U value) {
        if constexpr (std::is_same_v<U, Bool8>) {
            return static_cast<T>(value.value);
        } else if constexpr (std::is_same_v<T, Bool8>) {
            return Bool8(value != U{});
        } else {
            return static_cast<T>(value);
        }
    }
    [[nodiscard]] std::string type_mismatch(DType requested) const;

    ColumnData data_;
};

// pandas RangeIndex(start, start + size) or an explicit int64 index.
class Index {
  public:
    Index() = default;
    [[nodiscard]] static Index range(std::int64_t start, std::size_t size);
    [[nodiscard]] static Index values(std::vector<std::int64_t> values);

    [[nodiscard]] bool is_range() const noexcept { return !values_.has_value(); }
    [[nodiscard]] std::int64_t start() const noexcept { return start_; }
    [[nodiscard]] std::size_t size() const noexcept {
        return values_.has_value() ? values_->size() : size_;
    }
    [[nodiscard]] std::int64_t operator[](std::size_t i) const {
        return values_.has_value() ? (*values_)[i] : start_ + static_cast<std::int64_t>(i);
    }
    [[nodiscard]] std::vector<std::int64_t> to_vector() const;

    friend bool operator==(const Index& a, const Index& b);

  private:
    std::int64_t start_ = 0;
    std::size_t size_ = 0;
    std::optional<std::vector<std::int64_t>> values_;
};

class Table {
  public:
    Table() = default;
    Table(std::initializer_list<std::pair<std::string, Column>> columns);

    // len(df): the length of the index.
    [[nodiscard]] std::size_t num_rows() const noexcept { return index_.size(); }
    [[nodiscard]] std::size_t size() const noexcept { return num_rows(); }
    [[nodiscard]] std::size_t num_columns() const noexcept { return names_.size(); }
    [[nodiscard]] const std::vector<std::string>& columns() const noexcept { return names_; }
    [[nodiscard]] bool contains(std::string_view name) const noexcept;

    // df[name]; throws KeyError when absent.
    [[nodiscard]] const Column& operator[](std::string_view name) const;
    [[nodiscard]] Column& operator[](std::string_view name);
    [[nodiscard]] const Column& column(std::size_t position) const { return columns_.at(position); }
    [[nodiscard]] Column& column(std::size_t position) { return columns_.at(position); }

    // df[name] = column: replaces in place, or appends a new last column. The
    // first column of an empty table defines the row count with a RangeIndex.
    Table& set(std::string name, Column column);
    // Inserts at a position, like DataFrame.insert.
    Table& insert(std::size_t position, std::string name, Column column);
    Table& drop(std::string_view name);

    // df[[names]] (keeps the index).
    [[nodiscard]] Table select(const std::vector<std::string>& names) const;
    // Positional rows [lo, hi) with the index sliced alongside (df.iloc[lo:hi]).
    [[nodiscard]] Table slice(std::size_t lo, std::size_t hi) const;
    // df.iloc[positions].
    [[nodiscard]] Table take(std::span<const std::int64_t> positions) const;

    [[nodiscard]] const Index& index() const noexcept { return index_; }
    void set_index(Index index);
    // df.reset_index(drop=True)
    void reset_index() { index_ = Index::range(0, num_rows()); }

    [[nodiscard]] std::vector<std::pair<std::string, DType>> dtypes() const;

    friend bool operator==(const Table& a, const Table& b);

  private:
    std::vector<std::string> names_;
    std::vector<Column> columns_;
    Index index_;
};

}  // namespace coolercpp

#endif  // COOLERCPP_TABLE_HPP
