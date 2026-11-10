// The Cooler class and its table and matrix selectors (cooler/api.py and
// cooler/core/_selectors.py of cooler 0.10.2).

#ifndef COOLERCPP_API_HPP
#define COOLERCPP_API_HPP

#include <array>
#include <concepts>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "coolercpp/json.hpp"
#include "coolercpp/region.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp {

// A Python slice. Only step 1 (or None) is accepted, as in cooler.
struct Slice {
    std::optional<std::int64_t> start;
    std::optional<std::int64_t> stop;
    std::optional<std::int64_t> step;
};

// The column selection of a table selector: None (every column), a single
// name (the Python str form, which yields a Series: a one column Table) or a
// list of names.
class Fields {
  public:
    Fields() = default;
    Fields(const char* name) : kind_(Kind::Single), name_(name) {}           // NOLINT
    Fields(std::string name) : kind_(Kind::Single), name_(std::move(name)) {}  // NOLINT
    Fields(std::vector<std::string> names) : kind_(Kind::List), names_(std::move(names)) {}  // NOLINT
    Fields(std::initializer_list<std::string> names) : kind_(Kind::List), names_(names) {}

    [[nodiscard]] bool is_none() const noexcept { return kind_ == Kind::None; }
    [[nodiscard]] bool is_single() const noexcept { return kind_ == Kind::Single; }
    [[nodiscard]] bool is_list() const noexcept { return kind_ == Kind::List; }
    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] const std::vector<std::string>& names() const noexcept { return names_; }

  private:
    enum class Kind { None, Single, List };
    Kind kind_ = Kind::None;
    std::string name_;
    std::vector<std::string> names_;
};

// RangeSelector1D: DataFrame-like column selection and list-like row access
// over an out-of-core table.
class RangeSelector1D {
  public:
    using Slicer = std::function<Table(const Fields& fields, std::int64_t lo, std::int64_t hi)>;
    using Fetcher = std::function<std::pair<std::int64_t, std::int64_t>(const Region& region)>;

    RangeSelector1D(Fields fields, Slicer slicer, Fetcher fetcher, std::int64_t nmax);

    [[nodiscard]] const Fields& fields() const noexcept { return fields_; }
    [[nodiscard]] std::array<std::int64_t, 1> shape() const noexcept { return {nmax_}; }
    // len(selector)
    [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(nmax_); }

    // .columns, .dtypes, .keys() and `key in selector`, evaluated on an empty
    // slice as cooler does. A single-name selector raises AttributeError for
    // columns, keys and contains, as a pandas Series does.
    [[nodiscard]] std::vector<std::string> columns() const;
    [[nodiscard]] std::vector<std::pair<std::string, DType>> dtypes() const;
    [[nodiscard]] std::vector<std::string> keys() const;
    [[nodiscard]] bool contains(const std::string& key) const;

    // sel["A"], sel[{"A", "B"}]: a selector over a column subset.
    [[nodiscard]] RangeSelector1D operator[](Fields fields) const;
    [[nodiscard]] RangeSelector1D operator[](const char* field) const {
        return (*this)[Fields(field)];
    }
    // sel[start:stop]
    [[nodiscard]] Table operator[](const Slice& rows) const;
    // sel[i]: a one row table (cooler slices [i, i + 1)).
    template <std::integral I>
    [[nodiscard]] Table operator[](I row) const {
        return row_at(static_cast<std::int64_t>(row));
    }
    [[nodiscard]] Table slice(std::optional<std::int64_t> start,
                              std::optional<std::int64_t> stop) const {
        return (*this)[Slice{start, stop, std::nullopt}];
    }
    // sel[:]
    [[nodiscard]] Table all() const { return (*this)[Slice{}]; }
    // sel.fetch(region); NotImplementedError for the chroms selector.
    [[nodiscard]] Table fetch(const Region& region) const;

  private:
    [[nodiscard]] Table row_at(std::int64_t row) const;

    Fields fields_;
    Slicer slicer_;
    Fetcher fetcher_;
    std::int64_t nmax_ = 0;
};

// A 2D numpy.ndarray in row major order.
struct DenseMatrix {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    Column values;
    [[nodiscard]] DType dtype() const noexcept { return values.dtype(); }
    [[nodiscard]] double value(std::int64_t i, std::int64_t j) const {
        return values.as_double(static_cast<std::size_t>(i * cols + j));
    }
};

// scipy.sparse.coo_matrix: entries in the order cooler assembles them,
// duplicates not summed. row and col are int32 unless the shape or the
// coordinates need int64, as scipy decides.
struct SparseMatrix {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    Column row;
    Column col;
    Column data;
    [[nodiscard]] std::size_t nnz() const noexcept { return data.size(); }
};

// What a matrix selector returns: a dense array, a sparse matrix, or a pixel
// table (as_pixels=True).
class MatrixResult {
  public:
    MatrixResult(DenseMatrix dense) : value_(std::move(dense)) {}   // NOLINT
    MatrixResult(SparseMatrix sparse) : value_(std::move(sparse)) {}  // NOLINT
    MatrixResult(Table pixels) : value_(std::move(pixels)) {}       // NOLINT

    [[nodiscard]] bool is_dense() const noexcept { return value_.index() == 0; }
    [[nodiscard]] bool is_sparse() const noexcept { return value_.index() == 1; }
    [[nodiscard]] bool is_pixels() const noexcept { return value_.index() == 2; }
    [[nodiscard]] const DenseMatrix& dense() const;
    [[nodiscard]] const SparseMatrix& sparse() const;
    [[nodiscard]] const Table& pixels() const;
    [[nodiscard]] DenseMatrix& dense();
    [[nodiscard]] SparseMatrix& sparse();
    [[nodiscard]] Table& pixels();

  private:
    std::variant<DenseMatrix, SparseMatrix, Table> value_;
};

// One axis of a 2D subscript: a slice or a scalar.
class AxisKey {
  public:
    AxisKey(Slice slice) : value_(slice) {}  // NOLINT
    template <std::integral I>
    AxisKey(I index) : value_(static_cast<std::int64_t>(index)) {}  // NOLINT
    [[nodiscard]] const std::variant<Slice, std::int64_t>& value() const noexcept {
        return value_;
    }

  private:
    std::variant<Slice, std::int64_t> value_;
};

class RangeSelector2D {
  public:
    using Slicer = std::function<MatrixResult(const std::optional<std::string>& field,
                                              std::int64_t i0, std::int64_t i1,
                                              std::int64_t j0, std::int64_t j1)>;
    using Fetcher = std::function<std::array<std::int64_t, 4>(
        const Region& region, const std::optional<Region>& region2)>;

    RangeSelector2D(std::optional<std::string> field, Slicer slicer, Fetcher fetcher,
                    std::array<std::int64_t, 2> shape);

    [[nodiscard]] std::array<std::int64_t, 2> shape() const noexcept { return shape_; }
    [[nodiscard]] std::size_t size() const noexcept { return static_cast<std::size_t>(shape_[0]); }

    // sel[rows, cols]; sel[rows] is sel[rows, :].
    [[nodiscard]] MatrixResult operator()(const AxisKey& rows,
                                          const AxisKey& cols = Slice{}) const;
    [[nodiscard]] MatrixResult operator[](const AxisKey& rows) const { return (*this)(rows); }
    [[nodiscard]] MatrixResult all() const { return (*this)(Slice{}, Slice{}); }
    [[nodiscard]] MatrixResult fetch(const Region& region,
                                     const std::optional<Region>& region2 = std::nullopt) const;

  private:
    std::optional<std::string> field_;
    Slicer slicer_;
    Fetcher fetcher_;
    std::array<std::int64_t, 2> shape_{};
};

// The balance argument of Cooler.matrix: a bool, or the name of a bin table
// column holding the weights.
class Balance {
  public:
    Balance(bool enabled = true) : enabled_(enabled) {}              // NOLINT
    Balance(const char* column) : column_(std::string(column)) {}    // NOLINT
    Balance(std::string column) : column_(std::move(column)) {}      // NOLINT

    [[nodiscard]] bool is_column_name() const noexcept { return column_.has_value(); }
    // Python truthiness: an empty column name disables balancing.
    [[nodiscard]] bool enabled() const noexcept {
        return column_.has_value() ? !column_->empty() : enabled_;
    }
    [[nodiscard]] std::string column() const { return column_.value_or("weight"); }

  private:
    bool enabled_ = true;
    std::optional<std::string> column_;
};

// Cooler.matrix keyword arguments.
struct MatrixOptions {
    std::optional<std::string> field;
    Balance balance = true;
    bool sparse = false;
    bool as_pixels = false;
    bool join = false;
    bool ignore_index = true;
    std::optional<bool> divisive_weights;
    // cooler 0.10.2 calls this `chunksize` (earlier releases `max_chunk`).
    std::int64_t chunksize = 10000000;
};

class Cooler {
  public:
    // A path or a URI "file.mcool::/resolutions/10000".
    explicit Cooler(const std::string& store);
    // The deprecated root= form: a file path plus a group path.
    Cooler(const std::string& store, const std::string& root);

    [[nodiscard]] const std::string& filename() const noexcept { return filename_; }
    [[nodiscard]] const std::string& root() const noexcept { return root_; }
    [[nodiscard]] const std::string& uri() const noexcept { return uri_; }
    [[nodiscard]] const std::string& store() const noexcept { return store_; }

    // Re-reads the chromosome table and the metadata (Cooler._refresh).
    void refresh();

    [[nodiscard]] std::string storage_mode() const;
    // Cooler.binsize: nullopt for None; KeyError when "bin-size" is absent.
    [[nodiscard]] std::optional<std::int64_t> binsize() const;
    [[nodiscard]] const ChromSizes& chromsizes() const noexcept { return chromsizes_; }
    [[nodiscard]] std::vector<std::string> chromnames() const { return chromsizes_.names(); }
    [[nodiscard]] std::int64_t offset(const Region& region) const;
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> extent(const Region& region) const;
    // Cooler.info, read from the file on every call.
    [[nodiscard]] json::Value info() const;
    [[nodiscard]] std::array<std::int64_t, 2> shape() const;

    // **kwargs of the selectors reach cooler.core.get; convert_enum is the one
    // that changes results.
    [[nodiscard]] RangeSelector1D chroms(bool convert_enum = true) const;
    [[nodiscard]] RangeSelector1D bins(bool convert_enum = true) const;
    [[nodiscard]] RangeSelector1D pixels(bool join = false, bool convert_enum = true) const;
    [[nodiscard]] RangeSelector2D matrix(const MatrixOptions& options = {}) const;

    // The block of absolute bin indices [i0, i1) x [j0, j1), as cooler's
    // clr.matrix()[i0:i1, j0:j1] returns it (symmetric coolers are filled on
    // both sides). Added for downstream tile servers that address the genome
    // by absolute bin rather than by named region.
    [[nodiscard]] MatrixResult matrix_block(const std::string& field, std::int64_t i0,
                                            std::int64_t i1, std::int64_t j0,
                                            std::int64_t j1,
                                            const MatrixOptions& options = {}) const {
        return query_matrix(field, i0, i1, j0, j1, options);
    }

  private:
    [[nodiscard]] std::pair<std::int64_t, std::int64_t> region_extent(const Region& region) const;
    [[nodiscard]] MatrixResult query_matrix(const std::string& field, std::int64_t i0,
                                            std::int64_t i1, std::int64_t j0, std::int64_t j1,
                                            const MatrixOptions& options) const;
    [[nodiscard]] const json::Value& info_item(const std::string& key) const;

    std::string filename_;
    std::string root_ = "/";
    std::string uri_;
    std::string store_;
    ChromSizes chromsizes_;
    std::unordered_map<std::string, std::int64_t> chromids_;
    json::Value info_;
    bool is_symm_upper_ = true;
};

// cooler.annotate: joins bin annotations onto the bin1_id/bin2_id columns of
// a pixel table, suffixing them with 1 and 2.
[[nodiscard]] Table annotate(const Table& pixels, const Table& bins, bool replace = false);
[[nodiscard]] Table annotate(const Table& pixels, const RangeSelector1D& bins,
                             bool replace = false);

// cooler.fileops.list_coolers (needed by Cooler's error message).
[[nodiscard]] std::vector<std::string> list_coolers(const std::string& filepath);

}  // namespace coolercpp

#endif  // COOLERCPP_API_HPP
