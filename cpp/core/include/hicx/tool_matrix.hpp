// hicmatrix.HiCMatrix.hiCMatrix as the tools use it: construct from a file,
// operate on the matrix, save.
//
// This is deliberately not hicx::HiCMatrix (hic_matrix.hpp). That type is the
// read only view hicInfo needs. ToolMatrix is the mutable one, and it carries
// the two pieces of hiCMatrix behaviour that only matter once something is
// written back:
//
//  1. The load time field swap (cpp/PLAN.md 2.7 quirk 1). Every loader returns
//     (matrix, cut_intervals, nan_bins, distance_counts, correction_factors)
//     but hiCMatrix.__init__ unpacks the last two the other way round, and
//     hiCMatrix.save hands its members straight to set_matrix_variables. The
//     net effect is that an h5 to h5 round trip moves /correction_factors into
//     /distance_counts and a cool to cool round trip drops the weight column.
//     `data` holds the fields with hiCMatrix's meaning, which is also the
//     meaning the writers consume, so the swap happens once, at load.
//  2. hiCMatrix.save reuses the file handler built during the load, so the
//     output format is the format of the *input*, whatever the output file
//     name says. `hicSumMatrices -m a.h5 -o out.cool` writes an h5 file, named
//     out.cool.h5 because the h5 writer appends its own suffix. A name ending
//     in neither "cool" nor "h5" is written nowhere at all, silently, which is
//     what the --outFileName help text warns about.

#ifndef HICX_TOOL_MATRIX_HPP
#define HICX_TOOL_MATRIX_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/cool_adapter.hpp"
#include "hicx/matrix_data.hpp"

namespace hicx {

class ToolMatrix {
  public:
    // hm.hiCMatrix(path). `chromosome`, when given, is hiCMatrix's
    // pChrnameList with a single entry and reaches the cool loader and, v4's
    // own addition, the .hic loader (hic_adapter.hpp read_hic): both read only
    // that chromosome's, or "chrom:start-end" region's, block rather than the
    // whole file. An h5 path ignores it exactly as hicmatrix does.
    static ToolMatrix load(const std::string& path,
                           const std::optional<std::string>& chromosome = std::nullopt);

    // hiCMatrix.save(pMatrixName), whose defaults are pSymmetric=True and
    // pApplyCorrection=False. Returns false when the file name ends in neither
    // "cool" nor "h5" and nothing was written.
    bool save(const std::string& path);

    [[nodiscard]] MatrixData& data() noexcept { return data_; }
    [[nodiscard]] const MatrixData& data() const noexcept { return data_; }
    [[nodiscard]] CsrMatrix& matrix() noexcept { return data_.matrix; }
    [[nodiscard]] const CsrMatrix& matrix() const noexcept { return data_.matrix; }
    [[nodiscard]] const std::vector<std::int64_t>& nan_bins() const noexcept {
        return data_.nan_bins;
    }
    [[nodiscard]] const std::vector<CutInterval>& cut_intervals() const noexcept {
        return data_.cut_intervals;
    }
    // Which writer save() will use. hiCMatrix reuses the file handler built
    // during the load, so this is a property of the input, not of the output
    // name. A caller that needs to write the result some other way, as
    // hicTransform does for its streamed dense output, has to ask.
    [[nodiscard]] bool input_is_h5() const noexcept { return input_is_h5_; }

    // chrBinBoundaries, the ordered mapping the tools compare between inputs.
    [[nodiscard]] const std::vector<std::pair<std::string, BinRange>>& boundaries()
        const noexcept {
        return boundaries_;
    }
    // Every hiCMatrix mutator ends with a rebuild of the interval trees and of
    // chrBinBoundaries. A tool that changes the bin table calls this instead.
    void refresh_boundaries();

  private:
    MatrixData data_;
    std::vector<std::pair<std::string, BinRange>> boundaries_;
    bool input_is_h5_ = false;
    CoolSaveOptions cool_options_;
};

}  // namespace hicx

#endif  // HICX_TOOL_MATRIX_HPP
