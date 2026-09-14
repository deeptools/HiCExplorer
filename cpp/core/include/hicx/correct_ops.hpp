// The matrix surgery hicCorrectMatrix performs around the two balancing
// algorithms: the modified z-score outlier filter, the bin masking, the
// diagonal removal and the diagonal completion Knight-Ruiz needs.
//
// The masking deserves a note, because it is where the memory goes. The Python
// calls maskBins twice, once for the zero coverage bins and once for the MAD
// outliers, and each call slices the matrix twice
// (self.matrix[rows, :][:, cols]) while the previous matrix is still live. The
// second call first restores the earlier mask, so the net effect of the pair is
// a single removal of the union of the two sets, and hiCMatrix.save then puts
// the removed bins back as empty rows.
//
// This port performs that net effect in one in place compaction. The surviving
// entries are a subsequence of the stored ones in the same order, so they are
// rewritten from the front of the same arrays, and the arrays are then shrunk
// without releasing their capacity. Nothing is allocated that scales with the
// pixel count, and the restore at the end only renumbers columns and rebuilds
// the row offsets, both of which are also in place.

#ifndef HICX_CORRECT_OPS_HPP
#define HICX_CORRECT_OPS_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx::correct {

// hicCorrectMatrix.MAD. `points` is the per bin coverage vector.
//
//   median            = np.median(points[points > 0])
//   diff              = points - median
//   med_abs_deviation = np.median(np.abs(diff))
//   modified_z_score  = 0.6745 * diff / med_abs_deviation
//
// The median is numpy's: the middle element for an odd count, the mean of the
// two middle elements for an even one, and NaN with a warning for an empty
// selection.
class Mad {
  public:
    explicit Mad(const std::vector<double>& points);

    [[nodiscard]] const std::vector<double>& modified_z_scores() const noexcept {
        return z_;
    }
    [[nodiscard]] double median() const noexcept { return median_; }
    [[nodiscard]] double median_absolute_deviation() const noexcept { return mad_; }

  private:
    std::vector<double> z_;
    double median_ = 0.0;
    double mad_ = 0.0;
};

// numpy.median over a copy of `values`.
[[nodiscard]] double numpy_median(std::vector<double> values);

// hicCorrectMatrix.filter_by_zscore. The per bin coverage is the row sum of the
// symmetric matrix minus its diagonal, computed per chromosome when `perchr` is
// set and over the whole matrix otherwise, and a bin is an outlier when its
// modified z-score falls outside the two thresholds. Returns the bin ids in
// ascending order, as the Python's sorted() does.
[[nodiscard]] std::vector<std::int64_t> filter_by_zscore(
    CsrMatrix& matrix, const std::vector<std::pair<std::string, BinRange>>& boundaries,
    double lower_threshold, double upper_threshold, bool perchr, int threads,
    std::vector<std::string>* empty_chromosome_warnings);

// np.flatnonzero(matrix.sum(axis=1) == 0), the zero coverage bins ICE masks
// before it starts.
[[nodiscard]] std::vector<std::int64_t> zero_coverage_bins(CsrMatrix& matrix, int threads);

// Row sums of the symmetric matrix, matrix.sum(axis=1).
[[nodiscard]] std::vector<double> row_sums(CsrMatrix& matrix, int threads);

// The coverage per bin without self contacts of rows and columns [first, last):
// block.sum(axis=1) - block.diagonal(), the vector filter_by_zscore and
// hicCorrectMatrix's diagnostic plot compute the MAD of.
[[nodiscard]] std::vector<double> coverage_without_diagonal(CsrMatrix& matrix, std::int64_t first,
                                                            std::int64_t last, int threads);

// hiCMatrix.diagflat(value=0): drops every stored entry on the main diagonal.
// The Python subtracts the diagonal and adds a zero one, and the CSR addition
// drops the exactly zero results, so the entries disappear rather than becoming
// explicit zeros.
void remove_diagonal(CsrMatrix& matrix);

// The union of maskBins(zero bins) and maskBins(outliers) as a single in place
// compaction. `masked` marks the bins to remove and must have one entry per
// row. The bin table, the NaN bins and the correction factors follow, and the
// original bin ids are recorded in `removed` so that restore_masked_bins can
// undo it.
struct MaskState {
    // Original id of every surviving bin, ascending.
    std::vector<std::int64_t> kept;
    // Original id of every removed bin, ascending. This is what becomes
    // nan_bins after the restore.
    std::vector<std::int64_t> removed;
    // The bin table before the mask, needed to put it back.
    std::vector<CutInterval> original_intervals;
};

[[nodiscard]] MaskState mask_bins_in_place(MatrixData& data, const std::vector<char>& masked);

// hiCMatrix.restoreMaskedBins, which hiCMatrix.save performs before writing:
// the removed bins come back as empty rows and columns, nan_bins becomes the
// removed set, and the correction factors are padded with NaN at the removed
// positions. Also in place; only the row offsets are rebuilt.
void restore_masked_bins(MatrixData& data, const MaskState& state);

// Puts an explicit entry on every position of the main diagonal that does not
// have one, with value `value`. Knight-Ruiz balances A = M + 1e-5 * I and
// returns triu(A), so its output carries an entry on every diagonal position
// even where the input had none; reproducing the output pattern means creating
// them.
//
// Returns the number of entries inserted. When there is nothing to insert this
// is free; otherwise the two value arrays are rebuilt, which is the one place
// in this tool that allocates a second copy of the matrix. It is reached only
// when the corrected matrix is what gets written, that is when --outFileName
// ends in '.h5' (hicCorrectMatrix.py:756,759).
std::size_t add_missing_diagonal(CsrMatrix& matrix, double value);

// The order CPython iterates `set(bin_ids)` in.
//
// hicCorrectMatrix.py:668 writes the --filteredBed file with
// `for outlier_region in set(outlier_regions)`, so the order of the lines is
// the slot order of a CPython set, not the ascending order of the bin ids. The
// checked-in reference hicCorrectMatrix/filtered.bed has that order and the
// existing Python test compares it line by line, so it is observable output and
// is reproduced here rather than sorted.
//
// The simulation is of Objects/setobject.c: an open addressed table of 8 slots
// that grows to `used * 4` rounded up to a power of two whenever fill * 5
// reaches mask * 3, a probe sequence of up to ten consecutive slots followed by
// i = (i * 5 + 1 + perturb) & mask with perturb shifted right by five each
// round, and a rehash in old slot order on every resize. The hash of a Python
// int below 2^61 is the int itself, which is what makes this reproducible at
// all.
[[nodiscard]] std::vector<std::int64_t> cpython_set_order(
    const std::vector<std::int64_t>& keys);

// Drops every stored entry outside the diagonal chromosome blocks. --perchr
// assembles its result into a lil_matrix that is only ever written inside those
// blocks (hicCorrectMatrix.py:711,728), so every inter chromosomal contact is
// absent from the output.
void keep_only_diagonal_blocks(CsrMatrix& matrix,
                               const std::vector<std::pair<std::string, BinRange>>& boundaries);

}  // namespace hicx::correct

#endif  // HICX_CORRECT_OPS_HPP
