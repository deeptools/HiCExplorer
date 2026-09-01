// Matrix arithmetic and bin masking, at the level hicmatrix.HiCMatrix works at.
//
// Everything here follows scipy and numpy semantics rather than the obvious
// C++ one, because the results are compared against the Python reference:
//
//   * a binary operation on two CSR matrices walks the union of the two
//     sparsity patterns and stores a result only when it is not exactly zero,
//     which is what scipy's csr_binop_csr does. An operation that produces NaN
//     from a missing operand (0 * inf) therefore does end up stored.
//   * the result dtype is numpy's promotion of the two operand dtypes, and the
//     arithmetic happens in that dtype, so a float32 pair is added in single
//     precision.
//   * matrix.data.sum() is not matrix.sum(). It is numpy's pairwise reduction
//     over the data array in storage order, and for a matrix held as an upper
//     triangle that order is the order the symmetric matrix would have, not
//     the order of the stored triangle. data_sum() reproduces it without
//     materialising the symmetric matrix.
//   * hiCMatrix.maskBins followed by hiCMatrix.restoreMaskedBins, which is what
//     every tool that masks before saving performs, is a single operation on
//     the represented matrix: drop every entry in a masked row or column, keep
//     the shape, and record the masked bins as the new NaN bins. It also
//     upcasts the matrix to float64, because restoreMaskedBins pads with an
//     empty float64 block.

#ifndef HICX_MATRIX_OPS_HPP
#define HICX_MATRIX_OPS_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// numpy.promote_types over the dtype names that occur in the corpus.
[[nodiscard]] std::string promote_dtype(const std::string& a, const std::string& b);

// scipy's csr_matrix + csr_matrix, - and .multiply(). Both operands must carry
// the same Symmetry, and the result carries it too.
[[nodiscard]] CsrMatrix add(const CsrMatrix& a, const CsrMatrix& b);
[[nodiscard]] CsrMatrix subtract(const CsrMatrix& a, const CsrMatrix& b);
[[nodiscard]] CsrMatrix multiply_elementwise(const CsrMatrix& a, const CsrMatrix& b);

// numpy's ndarray.sum() over matrix.data, in the dtype of the matrix and in
// scipy's storage order. Integer matrices reduce in int64 and are therefore
// order independent; float matrices go through the pairwise reduction of
// numpy_compat, fed in the exact order of the represented matrix.
[[nodiscard]] Scalar data_sum(const CsrMatrix& matrix);

// matrix.data = matrix.data.astype(float) / divisor. Always float64, because
// astype(float) is float64 whatever the divisor is.
void divide_data_in_place(CsrMatrix& matrix, const Scalar& divisor);

// matrix.data = float(1) / matrix.data. A Python float is a weak scalar, so a
// float32 array stays float32 and everything else becomes float64.
void reciprocal_data_in_place(CsrMatrix& matrix);

// matrix.data = np.log2(matrix.data). numpy has its own float64 log2 which
// differs from libm's by up to an ulp on about a fifth of the values measured
// on the corpus; that is far inside the ED gate but it is why this tool is not
// claimed byte identical on the log2ratio path (cpp/PLAN.md 5.0).
void log2_data_in_place(CsrMatrix& matrix);

// chrBinBoundaries: hicmatrix.HiCMatrix.intervalListToIntervalTree's second
// return value, an ordered chromosome -> [first bin, last bin) mapping. The
// tools compare two of these for equality, so both the order and the ranges
// are observable. A chromosome that reappears after a different one overwrites
// its entry and keeps its original position, which is what an OrderedDict does.
[[nodiscard]] std::vector<std::pair<std::string, BinRange>> chrom_bin_boundaries(
    const std::vector<CutInterval>& cut_intervals);

// hiCMatrix.maskBins(bin_ids) followed by the hiCMatrix.restoreMaskedBins that
// hiCMatrix.save performs before writing.
//
// Net effect on the file that is written:
//   * every stored entry whose row or column is masked is discarded. On the
//     GSM2644945 plus GSM2644947 pair that silently loses 5,142 counts; it is
//     reproduced, not fixed (cpp/PLAN.md, hicSumMatrices.py:72).
//   * the shape and the bin table are unchanged.
//   * nan_bins becomes the masked set, which is the union of the requested
//     bins and the NaN bins the matrix already carried.
//   * the matrix dtype becomes float64, because restoreMaskedBins vstacks an
//     empty float64 block onto it.
//   * masked positions of the correction factors become NaN.
// An empty bin_ids is a no operation, including the dtype change, because
// maskBins returns before touching anything.
void mask_and_restore_bins(MatrixData& data, const std::vector<std::int64_t>& bin_ids);

}  // namespace hicx

#endif  // HICX_MATRIX_OPS_HPP
