// The short versus long range sums of hicPlotSVL
// (compute_relation_short_long_range of hicexplorer/hicPlotSVL.py), shared by
// the tool and the Python API (hicexplorer._hicx.hicPlotSVL).
//
// For a chromosome block matrix[first:last, first:last] of a symmetric
// matrix: the sum of the entries at most max_distance bins from the diagonal,
// the sum of the others, and their ratio, with numpy's arithmetic:
//  * an integer matrix sums in int64 (numpy's accumulator for int32 and
//    int64, wrapping), and the ratio is the float64 division;
//  * a float32 matrix sums in float32 with numpy's pairwise sum and divides
//    in float32; a float64 matrix the same in float64;
//  * the entries are summed in scipy's row major order of the block, which
//    for a matrix held as its upper triangle (Symmetry::UpperTriangle) is the
//    order of the filled symmetric matrix.
// hic_matrix.nonzero() skips stored zeros but hic_matrix.data does not, so a
// full block that stores explicit zeros makes the Python index the data with
// a mask of a different length; that is ExplicitZeros here.

#ifndef HICX_SVL_HPP
#define HICX_SVL_HPP

#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include "hicx/sparse_matrix.hpp"

namespace hicx::svl {

struct BlockSums {
    DType kind = DType::Float64;
    // Integer kind
    std::int64_t smaller_int = 0;
    std::int64_t greater_int = 0;
    // Float kinds, a float32 widened
    double smaller = 0.0;
    double greater = 0.0;
    // smaller / greater as the Python divides (float32 widened for Float32)
    double ratio = 0.0;
    // the ratio is neither inf nor nan, so the Python keeps the chromosome
    bool kept = false;
};

class ExplicitZeros : public std::runtime_error {
  public:
    ExplicitZeros(std::size_t stored, std::size_t nonzero)
        : std::runtime_error("the block stores explicit zeros"), stored_(stored), nonzero_(nonzero) {}
    // len(hic_matrix.data) and len(hic_matrix.nonzero()[0])
    [[nodiscard]] std::size_t stored() const noexcept { return stored_; }
    [[nodiscard]] std::size_t nonzero() const noexcept { return nonzero_; }

  private:
    std::size_t stored_;
    std::size_t nonzero_;
};

[[nodiscard]] BlockSums block_sums(const CsrMatrix& matrix, std::int64_t first, std::int64_t last,
                                   double max_distance);

}  // namespace hicx::svl

#endif  // HICX_SVL_HPP
