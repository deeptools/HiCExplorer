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

// The grouping hicMergeTADbins.merge_tad_bins builds: cut a new bin at every
// domain boundary bin id and at every chromosome change, and keep every group
// however short it is. That last point is the difference from plan_bin_merge,
// which drops a trailing group of fewer than num_bins/2 bins; a TAD is a TAD
// whatever its size. `boundaries` is the sorted set of bin ids
// get_boundary_bin_id returns, and the count that gates the cut is reset at
// every group, so a boundary at the first bin of a group is ignored.
[[nodiscard]] BinMergePlan plan_tad_merge(const std::vector<CutInterval>& intervals,
                                          const std::vector<std::int64_t>& boundaries);

// hicMergeMatrixBins.running_window_merge from the point where
// remove_nans_if_needed has already run and the trivial num_bins == 1 case has
// been handled by the caller.
//
// Every entry of the upper triangle is added into every cell of the
// num_bins x num_bins window centred on it, positions outside the matrix are
// dropped, the accumulated result is folded back to the upper triangle and
// mirrored as R + R.T - diag(R). The resolution and the bin table do not
// change, only the values.
//
// Two properties are reproduced rather than fixed:
//   * the window is applied to raw bin indices and therefore ignores
//     chromosome borders, so the last bins of one chromosome pick up counts
//     from the first bins of the next
//   * num_bins must be odd, which the Python asserts
[[nodiscard]] CsrMatrix running_window(const CsrMatrix& matrix, std::int64_t num_bins);

// hicMergeMatrixBins.merge_bins from the point where remove_nans_if_needed has
// already run. Returns the merged matrix, the merged bin table and the nan
// bins the Python recomputes as the all zero columns of the result. Correction
// factors and distance counts are carried through untouched, which is what the
// Python does: merge_bins never permutes them, so they only stay meaningful
// because the caller cleared them.
[[nodiscard]] MatrixData merge_bins(const MatrixData& input, std::int64_t num_bins);

// Chromosome name -> length pairs, in the shape hicx::read_chromosome_sizes
// returns.
using ChromosomeLengths = std::vector<std::pair<std::string, std::int64_t>>;

// v4-only deviation from hicMergeMatrixBins.py, used by the hicMergeMatrixBins
// tool only (see the deviation note at the top of hicMergeMatrixBins.cpp).
// plan_bin_merge builds its groups by walking `intervals`, the bin table the
// *input* actually has, so a bin missing from one input (for example because
// no reads were observed there) shrinks the group count and shifts every
// group after it. plan_bin_merge_genome instead builds the groups by walking
// the full set of bins the genome should have: every chromosome of
// `chromosome_lengths` tiled from position 0 in steps of `resolution`,
// exactly as read_two_dimensional_text does. A bin of `intervals` is mapped
// into that layout by (chromosome, start) and contributes to whichever group
// it falls in; a bin the layout expects but `intervals` does not have simply
// contributes nothing to its group, exactly as a present but all-zero bin
// already would. Two inputs covering the same genome at the same resolution
// therefore always produce the same group layout with the same
// `--chromosomeSizes` and the same num_bins, whatever bins either of them is
// missing.
//
// This only makes sense for a matrix with a single, fixed bin size: a
// restriction-fragment matrix has no such thing as "the bin a chromosome
// position belongs to" without the restriction cut positions, which
// chromosome_lengths does not carry, so the caller does not call this
// function for one (BinTable::bin_size_homogeneous() is false) and uses
// plan_bin_merge instead.
//
// Throws std::invalid_argument if num_bins, resolution or chromosome_lengths
// is unusable, the same as plan_bin_merge.
[[nodiscard]] BinMergePlan plan_bin_merge_genome(const std::vector<CutInterval>& intervals,
                                                  std::int64_t num_bins,
                                                  const ChromosomeLengths& chromosome_lengths,
                                                  std::int64_t resolution);

// merge_bins, built on plan_bin_merge_genome instead of plan_bin_merge. Same
// contract otherwise: nan_bins is recomputed from the result, correction
// factors and distance counts pass through untouched.
[[nodiscard]] MatrixData merge_bins_genome(const MatrixData& input, std::int64_t num_bins,
                                           const ChromosomeLengths& chromosome_lengths,
                                           std::int64_t resolution);

}  // namespace hicx

#endif  // HICX_REDUCE_MATRIX_HPP
