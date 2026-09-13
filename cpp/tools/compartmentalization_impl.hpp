// Computational core of hicCompartmentalization, in its own translation unit so
// that cpp/tests/test_compartmentalization.cpp links against exactly the code
// the tool runs. The design notes and the reproduced quirks are in
// hicCompartmentalization.cpp.

#ifndef HICX_TOOLS_COMPARTMENTALIZATION_IMPL_HPP
#define HICX_TOOLS_COMPARTMENTALIZATION_IMPL_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hicx/sparse_matrix.hpp"

namespace hicx::compartments {

// One line of the --pca bedgraph as pandas reads it with
// dtype={0: object, 1: Int64, 2: Int64, 3: float32}.
struct PcaRow {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    float pc1 = 0.0F;
};

// pd.read_table(path, header=None, sep="\t", dtype=...). Blank lines are
// skipped, as pandas' skip_blank_lines does; the fourth column is parsed as a
// float64 and narrowed to float32, which is what pandas' C parser does for a
// float32 dtype, and pandas' default NA spellings become NaN. Throws
// std::runtime_error for anything pandas would reject or turn into an Int64
// NA, because the Python then fails in getRegionBinRange.
[[nodiscard]] std::vector<PcaRow> read_pca_bedgraph(const std::string& path);

// np.nanquantile(values, quantiles) with the default method 'linear', numpy
// 1.26 arithmetic: virtual index (n - 1) * q, the two neighbours, and _lerp's
// two sided formula. Throws std::invalid_argument with numpy's message when a
// quantile is outside [0, 1] or NaN.
[[nodiscard]] std::vector<double> nanquantile_linear(std::vector<double> values,
                                                     const std::vector<double>& quantiles);

// np.linspace(start, stop, num) for float64 scalars, numpy 1.26 arithmetic.
// Throws std::invalid_argument for a negative num, as numpy does.
[[nodiscard]] std::vector<double> linspace(double start, double stop, std::int64_t num);

// np.searchsorted(sorted, keys, side='right'), including numpy's binary search
// that seeds each search from the previous key and its ordering in which NaN
// sorts after everything. The seeding only changes results for an array that
// is not sorted, which quantile boundaries can in principle be by an ulp.
[[nodiscard]] std::vector<std::int64_t> searchsorted_right(const std::vector<double>& sorted,
                                                           const std::vector<double>& keys);

// hicCompartmentalization.main:186-198: the quantile boundaries.
[[nodiscard]] std::vector<double> quantile_boundaries(const std::vector<PcaRow>& rows,
                                                      std::int64_t quantiles, double outliers);

// Full rows of a square matrix held as scipy's fillLowerTriangle represents it.
// For Symmetry::UpperTriangle the lower part of row r is column r of the upper
// triangle, found through a column index of 4 bytes per off diagonal entry
// plus 8 bytes per bin, which is a third of the working set instead of the
// whole of it that materialize_full() would cost.
class SymmetricRows {
  public:
    explicit SymmetricRows(const CsrMatrix& matrix);

    [[nodiscard]] std::int64_t size() const noexcept { return n_; }

    // Visits every stored entry of the represented row r as (column, value).
    // Columns arrive lower part first, each part in ascending order.
    template <class F>
    void for_each_in_row(std::int64_t r, F&& visit) const {
        const auto& indptr = matrix_->indptr();
        const auto& indices = matrix_->indices();
        const auto& data = matrix_->data();
        if (upper_) {
            const auto first = column_indptr_[static_cast<std::size_t>(r)];
            const auto last = column_indptr_[static_cast<std::size_t>(r) + 1];
            for (auto p = first; p < last; ++p) {
                const std::int32_t k = column_rows_[static_cast<std::size_t>(p)];
                visit(static_cast<std::int64_t>(k), data[position_in_row(k, r)]);
            }
        }
        const auto begin = indptr[static_cast<std::size_t>(r)];
        const auto end = indptr[static_cast<std::size_t>(r) + 1];
        for (auto p = begin; p < end; ++p) {
            visit(static_cast<std::int64_t>(indices[static_cast<std::size_t>(p)]),
                  data[static_cast<std::size_t>(p)]);
        }
    }

  private:
    [[nodiscard]] std::size_t position_in_row(std::int32_t row, std::int64_t column) const;

    const CsrMatrix* matrix_ = nullptr;
    std::int64_t n_ = 0;
    bool upper_ = false;
    std::vector<std::int64_t> column_indptr_;
    std::vector<std::int32_t> column_rows_;
};

// count_interactions (hicCompartmentalization.py:93-129) followed by the
// np.nan_to_num of :212, as a quantiles-by-quantiles row major array.
//
//   bin_ids[k]          the bins of pca row k, get_indices of :87
//   quantile_of_row[k]  the searchsorted quantile of pca row k
//   offsets             --offset, every value already checked >= 0
//   chromosome_count    len(pc1["chr"].unique()); see quirk 3 in the tool
[[nodiscard]] std::vector<double> normalised_sum_per_quantile(
    const CsrMatrix& matrix, const std::vector<std::vector<std::int64_t>>& bin_ids,
    const std::vector<std::int64_t>& quantile_of_row, std::int64_t quantiles,
    const std::vector<std::int64_t>& offsets, std::int64_t chromosome_count);

// np.nan_to_num with its defaults: NaN to 0, +-inf to +-DBL_MAX.
void nan_to_num_in_place(std::vector<double>& values);

// within_vs_between_compartments (hicCompartmentalization.py:132-153).
[[nodiscard]] std::vector<double> within_vs_between(const std::vector<double>& normalised,
                                                    std::int64_t quantiles);

// One line of np.savetxt(fname, X) with its default fmt '%.18e', delimiter ' '
// and newline '\n'. Python's '%e' spells every NaN 'nan', never '-nan'.
[[nodiscard]] std::string savetxt_line(const std::vector<double>& row);

}  // namespace hicx::compartments

#endif  // HICX_TOOLS_COMPARTMENTALIZATION_IMPL_HPP
