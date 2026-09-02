// Observed over expected and z-score transforms of a contact matrix.
//
// Port of hicmatrix.HiCMatrix.convert_to_obs_exp_matrix (HiCMatrix.py:349-556),
// which convert_to_zscore_matrix is a one line wrapper around. It is a core
// component rather than part of a tool because hicFindTADs, hicTransform,
// hicPCA and hicDetectLoops all go through it.
//
// What the Python does, and what this reproduces exactly:
//
//  1. Truncate to the band 0 <= column - row < int(maxdepth * 1.5 / binsize)
//     of the upper triangle. Entries outside the band cancel to zero and are
//     dropped, *unless* they are NaN, because NaN - NaN is NaN and scipy's
//     sparse subtraction stores every result that is not exactly zero.
//  2. For the z-score, add a sparse band of ones over the same width and
//     subtract one again after converting to COO. The net effect is that every
//     position of the band becomes an explicitly stored entry, zero where the
//     matrix had nothing, so that the per diagonal mean and standard deviation
//     count the zeros. This port never materialises that band: it walks the
//     band positions directly, which is where most of the memory saving over
//     the Python comes from.
//  3. Per diagonal, sum the values with np.bincount, which accumulates
//     sequentially in COO order, and divide by the diagonal length. The
//     diagonal length is the number of positions the diagonal would have in a
//     dense matrix, raised to the observed count when that is larger.
//  4. For the z-score, the standard deviation counts the missing zeros too:
//     sum((v - mu)^2) over the stored values, plus (length - stored) * mu^2,
//     divided by the length. The stored part is summed in numpy's pairwise
//     order, not sequentially.
//  5. Divide, or standardise, every value in the band; a value further from
//     the diagonal than int(maxdepth * 1.5 / binsize) + 1 becomes zero.
//  6. Write the result back through a LIL matrix, which silently drops every
//     exact zero. NaN survives, so a diagonal with zero standard deviation
//     leaves NaN entries behind.
//
// Only the branch the tools actually take is implemented: maxdepth is
// required. Without it the Python densifies the whole matrix, which no ported
// tool asks for, and pretending to support it would mean shipping an untested
// path.

#ifndef HICX_OBSEXP_OPS_HPP
#define HICX_OBSEXP_OPS_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"

namespace hicx {

struct ObsExpOptions {
    // The `maxdepth` argument in base pairs. hicFindTADs passes
    // max_depth * 2.5.
    double max_depth_bp = 0.0;
    // zscore=True standardises, zscore=False divides by the expected value.
    bool zscore = true;
    // perchr=True processes one chromosome block at a time, which is what
    // every caller in HiCExplorer does.
    bool perchr = true;
};

// hiCMatrix.convert_to_obs_exp_matrix. `bin_size` is hiCMatrix.getBinSize(),
// which the caller has to pass because hiCMatrix caches it from the bin table
// as it was when the object was constructed, not as it is now.
//
// On return `data.matrix` holds the transformed upper triangle with
// Symmetry::Full, because the transform destroys the symmetry the loader
// established: only the upper triangle carries values.
// This runs on one thread. Threading it over the chromosome blocks was
// implemented and measured, and reverted because it cost memory and bought no
// time; the numbers are in the source.
void convert_to_obs_exp_matrix(MatrixData& data, std::int64_t bin_size,
                               const ObsExpOptions& options);

// hiCMatrix.fit_cut_intervals: when more than one percent of the bins deviate
// from the median bin length, every start and end is snapped to the nearest
// multiple of that median. Exposed because the distance list the transform
// builds is derived from the snapped starts, and a caller that wants to
// predict a distance has to snap the same way.
[[nodiscard]] std::vector<CutInterval> fit_cut_intervals(
    const std::vector<CutInterval>& intervals);

// hicexplorer.utilities.enlarge_bins: closes the gaps a masked bin leaves by
// moving the boundary between two neighbouring bins to the middle of the gap,
// and starts every chromosome at zero. Mutates and returns, as the Python does.
void enlarge_bins(std::vector<CutInterval>& intervals);

}  // namespace hicx

#endif  // HICX_OBSEXP_OPS_HPP
