// The five values every hicmatrix reader returns and every writer consumes.
//
// hicmatrix.lib.matrixFile.MatrixFile keeps exactly these as instance state
// and MatrixFileHandler.set_matrix_variables fills them in before a save. The
// C++ writers take the same aggregate so that a port of a tool can hand over
// what it loaded without an intermediate representation.
//
// Note the load time field swap of cpp/PLAN.md 2.7 quirk 1: the loaders return
// (matrix, cut_intervals, nan_bins, distance_counts, correction_factors) while
// hiCMatrix.__init__ unpacks the last two the other way round. This struct
// names the fields after what the *loader* means, which is also what
// hicConvertFormat passes back into set_matrix_variables. A caller that goes
// through hiCMatrix semantics has to swap them itself, and that swap is a
// property of hiCMatrix, not of the file layer.

#ifndef HICX_MATRIX_DATA_HPP
#define HICX_MATRIX_DATA_HPP

#include <cstdint>
#include <optional>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

struct MatrixData {
    CsrMatrix matrix;
    std::vector<CutInterval> cut_intervals;
    std::vector<std::int64_t> nan_bins;
    std::optional<std::vector<double>> correction_factors;
    std::optional<std::vector<double>> distance_counts;

    // hicCorrectMatrix --perchr with KR stores the correction factors as an
    // n-by-1 column rather than a flat vector, because krbalancing hands back
    // an Eigen column vector and PyTables writes the array's shape verbatim.
    // One matrix in the corpus is like that
    // (hicCorrectMatrix/small_test_matrix_KRcorrected_chrUextra_chr3LHet.h5),
    // and an h5 round trip is only value exact if the shape survives it.
    bool correction_factors_are_column = false;
};

}  // namespace hicx

#endif  // HICX_MATRIX_DATA_HPP
