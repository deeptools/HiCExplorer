// The hiCMatrix mutators that select, reorder and delete bins.
//
// matrix_ops.hpp already carries the one combination every tool that masks
// before saving performs, mask_and_restore_bins. These are the other three,
// and they are what hicAdjustMatrix and hicMergeMatrixBins are built out of:
//
//   * hiCMatrix.reorderBins (HiCMatrix.py:723), matrix[order, :][:, order]
//     with the bin table and the NaN bins carried along. It deliberately does
//     not permute the correction factors, which is quirk 3 of cpp/PLAN.md 2.7,
//     and that is reproduced here rather than fixed.
//   * hiCMatrix.maskBins (HiCMatrix.py:761) *without* the restoreMaskedBins
//     that hiCMatrix.save would perform. Two callers clear orig_bin_ids right
//     after masking, so for them the bins really disappear:
//     hicAdjustMatrix.py:155-159 (--action remove) and
//     hicMergeMatrixBins.remove_nans_if_needed.
//   * the block zeroing of hicAdjustMatrix --interIntraHandling.
//
// All three work on the represented matrix, so a matrix held as an upper
// triangle stays one whenever the operation allows it: a selection whose index
// list is strictly increasing is order preserving and therefore maps the
// triangle onto a triangle. Only a genuine reordering needs both triangles,
// which is the exemption hicAdjustMatrix takes in cpp/PLAN.md 4.4 rule 2.

#ifndef HICX_ADJUST_OPS_HPP
#define HICX_ADJUST_OPS_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// matrix[order, :][:, order] with scipy's fancy indexing semantics: an index
// that appears twice duplicates its row and its column, and an index that
// appears in no position is dropped. The column indices of every output row
// come out sorted, as scipy's do.
[[nodiscard]] CsrMatrix select_bins(const CsrMatrix& matrix,
                                    const std::vector<std::int64_t>& order);

// hiCMatrix.reorderBins. The matrix, the bin table and the NaN bins follow the
// new order; the correction factors are left exactly as they were, which is
// the Python behaviour and not an oversight of the port.
void reorder_bins(MatrixData& data, const std::vector<std::int64_t>& order);

// hiCMatrix.maskBins with no restore afterwards, so the bins are gone for
// good. The bins that are dropped are the union of `bin_ids` and the NaN bins
// the matrix carries, matching the join at HiCMatrix.py:794-799, and nan_bins
// ends up empty. An empty `bin_ids` is a no operation, because maskBins
// returns before touching anything.
void delete_bins(MatrixData& data, const std::vector<std::int64_t>& bin_ids);

// np.flatnonzero(matrix.sum(0).A == 0): the columns of the represented matrix
// whose sum is exactly zero. hicMergeMatrixBins and hicMergeTADbins both use
// it to rebuild nan_bins after a merge.
[[nodiscard]] std::vector<std::int64_t> empty_column_bins(const CsrMatrix& matrix);

enum class InterIntra {
    // hic_matrix.matrix[start:end, end:] = 0 for every chromosome, which zeroes
    // every entry whose row is on an earlier chromosome than its column.
    Inter,
    // hic_matrix.matrix[start:end, start:end] = 0, the diagonal blocks.
    Intra,
};

// hicAdjustMatrix.py:170-189. The Python writes explicit zeros into the CSR
// and both writers call eliminate_zeros before storing, so dropping the
// entries outright produces the same file.
//
// Note that the Inter branch zeroes only the blocks *above* the diagonal: for
// a chromosome occupying [start, end) it zeroes the columns from `end` on and
// never the columns before `start`. The mirrored entries survive in memory and
// are then discarded by the triu(k=0) both writers apply, so the file is the
// one the tool's help text promises. Reproduced literally, on the upper
// triangle, which is where those entries are.
void zero_inter_or_intra(CsrMatrix& matrix,
                         const std::vector<std::pair<std::string, BinRange>>& boundaries,
                         InterIntra mode);

}  // namespace hicx

#endif  // HICX_ADJUST_OPS_HPP
