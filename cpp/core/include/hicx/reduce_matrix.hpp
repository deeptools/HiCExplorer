// Bin merging: the reduction that turns a matrix at one resolution into a
// matrix at a coarser one.
//
// Two Python functions live here, both of which more than one tool needs:
//
//   * hicexplorer/reduceMatrix.py reduce_matrix, which sums the rows and
//     columns of each group of bins. hicMergeMatrixBins, hicMergeTADbins and
//     hicBuildMatrix all call it.
//   * hicexplorer/hicMergeMatrixBins.py merge_bins, which builds the groups
//     from a bin count and then calls reduce_matrix. hicConvertFormat calls it
//     directly for every resolution of an mcool
//     (hicConvertFormat.py:314-315), and the hicMergeMatrixBins tool is a thin
//     command line shell over it.
//
// What is deliberately *not* here is hicMergeMatrixBins.remove_nans_if_needed,
// which drops NaN bins through hiCMatrix.maskBins before merging. That is a
// hiCMatrix operation, it belongs with the masking semantics, and
// hicConvertFormat never reaches it: the mcool path rebuilds the matrix with
// setMatrix, which leaves nan_bins unset, so the guard at
// hicMergeMatrixBins.py:67 is false. The caller runs it, merge_bins does not.

#ifndef HICX_REDUCE_MATRIX_HPP
#define HICX_REDUCE_MATRIX_HPP

#include <cstdint>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// Port of reduceMatrix.reduce_matrix.
//
// `bins_to_merge[k]` lists the row and column indices of the input that become
// index k of the output. An index that appears in no group is dropped, which
// is how the Python removes the bins of a chromosome that is too short to
// merge. Behaviours that are reproduced rather than cleaned up:
//
//   * when the number of groups equals the number of rows the input is
//     returned unchanged, without the triangle selection and without the
//     symmetrisation (reduceMatrix.py:154-155)
//   * the summed values are cast back to the dtype of the input, so an int32
//     matrix truncates the float64 sums numpy's bincount produced
//   * with pDiagonal false the main diagonal ends up exactly zero, because the
//     Python subtracts twice the diagonal from the symmetrised result
CsrMatrix reduce_matrix(const CsrMatrix& matrix,
                        const std::vector<std::vector<std::int64_t>>& bins_to_merge,
                        bool use_triu = true, bool diagonal = false);

// The groups and the new bin table that merge_bins derives from a bin count.
struct BinMergePlan {
    std::vector<CutInterval> intervals;
    std::vector<std::vector<std::int64_t>> bins_to_merge;
};

// The first half of hicMergeMatrixBins.merge_bins: walk the bin table, cut a
// new bin every `num_bins` bins and at every chromosome change, and drop a
// trailing group of fewer than num_bins/2 bins. The last group is kept
// whatever its size, which the Python does by appending it outside the loop.
// The coverage of a merged bin is the mean of the coverages it absorbs.
[[nodiscard]] BinMergePlan plan_bin_merge(const std::vector<CutInterval>& intervals,
                                          std::int64_t num_bins);

// hicMergeMatrixBins.merge_bins from the point where remove_nans_if_needed has
// already run. Returns the merged matrix, the merged bin table and the nan
// bins the Python recomputes as the all zero columns of the result. Correction
// factors and distance counts are carried through untouched, which is what the
// Python does: merge_bins never permutes them, so they only stay meaningful
// because the caller cleared them.
[[nodiscard]] MatrixData merge_bins(const MatrixData& input, std::int64_t num_bins);

}  // namespace hicx

#endif  // HICX_REDUCE_MATRIX_HPP
