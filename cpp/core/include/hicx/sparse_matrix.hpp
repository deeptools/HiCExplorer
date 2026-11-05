// Sparse contact matrix with scipy.sparse.csr_matrix semantics.
//
// HiCExplorer stores a Hi-C contact map as a square scipy CSR matrix. Only the
// operations the tools actually rely on are implemented, but they follow scipy
// exactly, down to the order in which floating point values are accumulated:
//
//   * matrix.sum() is not a plain loop over data. scipy multiplies by a dense
//     vector of ones, which sums every row sequentially, and then reduces the
//     resulting dense vector with numpy's pairwise summation.
//   * adding two CSR matrices drops entries whose sum is exactly zero.
//   * matrix.nnz counts stored entries, including explicit zeros.
//
// Memory layout, and why it is what it is:
//
//   * Column indices are int32, one per stored entry. A cooler cannot have
//     more than 2^31 bins in practice, so int64 indices would double the
//     largest array in the process for no gain.
//   * Row offsets are int64, one per row plus one. That array is O(rows), not
//     O(nnz): 200 kB for the 24,926 bin gm12878 matrix. Narrowing it to int32
//     would save nothing measurable, while widening it is mandatory as soon as
//     nnz exceeds 2^31, so int64 is the memory-optimal choice at every size.
//   * Symmetry is a property of the matrix. A Hi-C map is symmetric and both
//     the cool and the h5 format store only the upper triangle. Marking the
//     matrix Symmetry::UpperTriangle keeps it that way and answers every query
//     as if the lower triangle were present, which halves the resident set
//     compared with materialising it. materialize_full() is the only operation
//     that allocates a second full matrix, and it says so in its name.
//   * Values are held as double regardless of the stored dtype. This costs a
//     factor of two on int32 and float32 matrices and is the one known
//     remaining inefficiency of this type; see cpp/STATUS.md.

#ifndef HICX_SPARSE_MATRIX_HPP
#define HICX_SPARSE_MATRIX_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace hicx {

// The numpy dtype of the stored counts. It decides in which precision scipy
// accumulates, and therefore what hicInfo prints.
enum class DType { Integer, Float32, Float64 };

[[nodiscard]] DType dtype_from_name(const std::string& name);

// A numeric value that remembers the dtype it came from.
struct Scalar {
    DType kind = DType::Float64;
    std::int64_t integer_value = 0;
    double float_value = 0.0;

    [[nodiscard]] bool is_integer() const { return kind == DType::Integer; }
    [[nodiscard]] double as_double() const {
        return kind == DType::Integer ? static_cast<double>(integer_value) : float_value;
    }
    static Scalar from_int(std::int64_t value) {
        return Scalar{DType::Integer, value, 0.0};
    }
    static Scalar from_double(double value) {
        return Scalar{DType::Float64, 0, value};
    }
    static Scalar from_float(float value) {
        return Scalar{DType::Float32, 0, static_cast<double>(value)};
    }
};

// Whether the stored entries are the whole matrix or only its upper triangle.
enum class Symmetry {
    // Every stored entry is an entry of the matrix.
    Full,
    // Only entries with column >= row are stored. The represented matrix is
    // U + triu(U, 1).T with entries that are exactly zero dropped, which is
    // what hicmatrix.HiCMatrix.fillLowerTriangle produces.
    UpperTriangle,
};

class CsrMatrix {
  public:
    CsrMatrix() = default;
    CsrMatrix(std::int64_t rows, std::int64_t cols, std::vector<std::int64_t> indptr,
              std::vector<std::int32_t> indices, std::vector<double> data,
              std::string dtype);

    // Builds a CSR matrix from coordinate triplets, like
    // csr_matrix((data, (row, col)), shape=(n, n)): duplicate entries are
    // summed and the column indices of every row end up sorted.
    static CsrMatrix from_coo(std::int64_t rows, std::int64_t cols,
                              const std::vector<std::int32_t>& row,
                              const std::vector<std::int32_t>& col,
                              std::vector<double> data, std::string dtype);

    [[nodiscard]] std::int64_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::int64_t cols() const noexcept { return cols_; }

    [[nodiscard]] Symmetry symmetry() const noexcept { return symmetry_; }
    // Declares what the stored entries mean. Only an operation that knows it
    // preserves the triangle may use this: a selection whose index list is
    // strictly increasing is order preserving, so it maps triu onto triu
    // (hicx::select_bins). Everything else must go through materialize_full.
    void set_symmetry(Symmetry symmetry) noexcept { symmetry_ = symmetry; }

    // Entries physically held in memory.
    [[nodiscard]] std::size_t stored_nnz() const noexcept { return data_.size(); }
    // matrix.nnz of the represented matrix. For Symmetry::UpperTriangle the
    // off diagonal entries count twice and exact zeros do not count, because
    // that is what the CSR addition in fillLowerTriangle produces.
    [[nodiscard]] std::size_t nnz() const;
    [[nodiscard]] const std::string& dtype() const noexcept { return dtype_; }
    void set_dtype(std::string dtype) { dtype_ = std::move(dtype); }
    [[nodiscard]] DType dtype_kind() const { return dtype_from_name(dtype_); }
    [[nodiscard]] bool integral_dtype() const {
        return dtype_kind() == DType::Integer;
    }

    [[nodiscard]] const std::vector<std::int64_t>& indptr() const noexcept {
        return indptr_;
    }
    [[nodiscard]] const std::vector<std::int32_t>& indices() const noexcept {
        return indices_;
    }
    [[nodiscard]] const std::vector<double>& data() const noexcept { return data_; }
    [[nodiscard]] std::vector<double>& mutable_data() noexcept { return data_; }

    // scipy's eliminate_zeros: drops stored entries that are exactly zero.
    void eliminate_zeros();

    // True when the strict lower triangle sums to zero, the condition
    // hicmatrix.HiCMatrix.fillLowerTriangle tests.
    [[nodiscard]] bool lower_triangle_is_zero() const;
    // True when at least one entry with column < row is stored.
    [[nodiscard]] bool has_lower_entries() const;

    // hicmatrix.HiCMatrix.fillLowerTriangle without the allocation: when
    // nothing is stored below the diagonal the matrix is simply relabelled
    // Symmetry::UpperTriangle. Falls back to materialize_full() in the corner
    // case where entries below the diagonal exist but cancel to zero, because
    // the Python code adds the transpose there as well.
    void symmetrize_in_place();

    // Expands to the explicit symmetric matrix: self + triu(self, 1).T.
    // Allocates a second full matrix, so call it only when a tool needs
    // explicit lower triangle rows.
    void materialize_full();

    // matrix.sum() with scipy's evaluation order.
    [[nodiscard]] Scalar sum() const;

    // matrix.diagonal(), a dense vector of length min(rows, cols).
    [[nodiscard]] std::vector<double> diagonal() const;
    // matrix.diagonal().sum() with numpy's evaluation order.
    [[nodiscard]] Scalar diagonal_sum() const;

    [[nodiscard]] Scalar data_min() const;
    [[nodiscard]] Scalar data_max() const;

    // Dense element access, for tests and small matrices only.
    [[nodiscard]] double at(std::int64_t row, std::int64_t col) const;

    // triu(self, k=0) after eliminate_zeros, which is what both file writers
    // store, expressed without building it.
    //
    // The row offsets are the CSR indptr of that triangle, and they are also
    // the cool /indexes/bin1_offset array. The array is O(rows), so producing
    // it costs 8 bytes per bin and no copy of the matrix.
    [[nodiscard]] std::vector<std::int64_t> upper_triangle_indptr() const;
    // Number of entries that triu(self, k=0) would store.
    [[nodiscard]] std::size_t upper_triangle_nnz() const;

    // Visits every entry of triu(self, k=0) in row major order, exact zeros
    // skipped, calling visit(row, column, value). Nothing is allocated, which
    // is what lets the writers stream.
    template <class F>
    void for_each_upper(F&& visit) const {
        for_each_selected(true, std::forward<F>(visit));
    }
    // The same over every stored entry, for pSymmetric=False.
    template <class F>
    void for_each_stored(F&& visit) const {
        for_each_selected(false, std::forward<F>(visit));
    }
    [[nodiscard]] std::vector<std::int64_t> stored_indptr_without_zeros() const;
    [[nodiscard]] std::size_t nonzero_stored_nnz() const;

    // The three CSR arrays and the metadata that goes with them, detached from
    // the matrix.
    //
    // release() leaves the matrix empty and hands the caller the buffers
    // themselves, and adopt() puts them back. An operation that only shrinks
    // the matrix or renumbers its rows can then rewrite the arrays in place and
    // reuse the same allocation, instead of building a second matrix beside the
    // first. That matters for one case in particular: masking the zero and
    // outlier bins out of the 61.8 M pixel gm12878 matrix would otherwise cost
    // a transient 1.5 GB against a 954 MB budget (cpp/PLAN.md 4.5), because the
    // value array alone is 494 MB.
    struct Arrays {
        std::int64_t rows = 0;
        std::int64_t cols = 0;
        std::vector<std::int64_t> indptr;
        std::vector<std::int32_t> indices;
        std::vector<double> data;
        std::string dtype = "float64";
        Symmetry symmetry = Symmetry::Full;
    };
    [[nodiscard]] Arrays release();
    [[nodiscard]] static CsrMatrix adopt(Arrays arrays);

  private:
    template <class F>
    void for_each_selected(bool upper_only, F&& visit) const {
        for (std::int64_t row = 0; row < rows_; ++row) {
            const std::size_t begin =
                static_cast<std::size_t>(indptr_[static_cast<std::size_t>(row)]);
            const std::size_t end =
                static_cast<std::size_t>(indptr_[static_cast<std::size_t>(row) + 1]);
            for (std::size_t k = begin; k < end; ++k) {
                if (data_[k] == 0.0) {
                    continue;  // eliminate_zeros
                }
                if (upper_only && indices_[k] < row) {
                    continue;  // triu(k=0)
                }
                visit(row, static_cast<std::int64_t>(indices_[k]), data_[k]);
            }
        }
    }

    std::int64_t rows_ = 0;
    std::int64_t cols_ = 0;
    std::vector<std::int64_t> indptr_{0};
    std::vector<std::int32_t> indices_;
    std::vector<double> data_;
    std::string dtype_ = "float64";
    Symmetry symmetry_ = Symmetry::Full;
};

}  // namespace hicx

#endif  // HICX_SPARSE_MATRIX_HPP
