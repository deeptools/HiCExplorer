// Port of hicmatrix.HiCMatrix.hiCMatrix, the class every HiCExplorer tool
// loads its input through.
//
// Loading a matrix means: pick the format from the file name suffix, read the
// contact matrix and the bin table, apply the balancing weights of a cool file,
// mirror the upper triangle into the lower one, and build the per chromosome
// interval lookup.

#ifndef HICX_HIC_MATRIX_HPP
#define HICX_HIC_MATRIX_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

class HiCMatrix {
  public:
    struct Options {
        // pUpperTriangleOnly=True skips the symmetric completion.
        bool fill_lower_triangle = true;
        // Expand the symmetric matrix into explicit lower triangle rows. The
        // default keeps the upper triangle only, which halves the resident
        // set; every query on CsrMatrix answers for the symmetric matrix
        // either way. Set this only for a tool that needs to walk explicit
        // rows of the lower half.
        bool materialize_full_matrix = false;
        // Cool.applyCorrectionLoad, multiplies the counts with the weights.
        bool apply_correction = true;
        // Cool.correctionFactorTable, the bin column holding the weights.
        std::string correction_factor_table = "weight";
    };

    HiCMatrix() = default;
    static HiCMatrix load(const std::string& path, const Options& options);
    static HiCMatrix load(const std::string& path) { return load(path, Options()); }

    [[nodiscard]] const CsrMatrix& matrix() const noexcept { return matrix_; }
    [[nodiscard]] CsrMatrix& matrix() noexcept { return matrix_; }
    [[nodiscard]] const BinTable& bins() const noexcept { return bins_; }
    [[nodiscard]] const std::vector<std::int64_t>& nan_bins() const noexcept {
        return nan_bins_;
    }
    [[nodiscard]] const std::optional<std::vector<double>>& correction_factors()
        const noexcept {
        return correction_factors_;
    }

    [[nodiscard]] std::int64_t bin_size() const { return bins_.bin_size(); }
    [[nodiscard]] std::vector<std::pair<std::string, std::int64_t>>
    chromosome_sizes() const {
        return bins_.chromosome_sizes();
    }

  private:
    CsrMatrix matrix_;
    BinTable bins_;
    std::vector<std::int64_t> nan_bins_;
    std::optional<std::vector<double>> correction_factors_;
};

}  // namespace hicx

#endif  // HICX_HIC_MATRIX_HPP
