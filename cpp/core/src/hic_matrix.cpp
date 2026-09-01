#include "hicx/hic_matrix.hpp"

#include <algorithm>
#include <cmath>

#include "hicx/cool_file.hpp"
#include "hicx/h5_file.hpp"

namespace hicx {

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// hicmatrix.lib.Cool.load, restricted to the full matrix case that the tools
// use when no chromosome is preselected.
void load_cool(const std::string& path, const HiCMatrix::Options& options,
               CsrMatrix& matrix, std::vector<CutInterval>& cut_intervals,
               std::vector<std::int64_t>& nan_bins,
               std::optional<std::vector<double>>& correction_factors) {
    const CoolFile cool(path);
    matrix = cool.read_matrix();
    cut_intervals = cool.read_bins();

    std::optional<std::vector<double>> weights;
    if (options.apply_correction && cool.has_column(options.correction_factor_table)) {
        weights = cool.read_column(options.correction_factor_table);
    }

    if (weights.has_value()) {
        matrix.eliminate_zeros();
        if (matrix.nnz() > 1) {
            const std::vector<double>& factors = *weights;
            const bool all_nan = std::all_of(factors.begin(), factors.end(),
                                             [](double v) { return std::isnan(v); });
            if (!all_nan) {
                // 'weight' is multiplicative, the hic2cool tables KR, VC and
                // SQRT_VC are divisive.
                const bool divide = options.correction_factor_table == "KR" ||
                                    options.correction_factor_table == "VC" ||
                                    options.correction_factor_table == "SQRT_VC";
                std::vector<double>& data = matrix.mutable_data();
                const std::vector<std::int64_t>& indptr = matrix.indptr();
                const std::vector<std::int32_t>& indices = matrix.indices();
                for (std::int64_t row = 0; row < matrix.rows(); ++row) {
                    const std::size_t begin =
                        static_cast<std::size_t>(indptr[static_cast<std::size_t>(row)]);
                    const std::size_t end =
                        static_cast<std::size_t>(indptr[static_cast<std::size_t>(row) + 1]);
                    for (std::size_t k = begin; k < end; ++k) {
                        const double factor =
                            factors[static_cast<std::size_t>(row)] *
                            factors[static_cast<std::size_t>(indices[k])];
                        if (divide) {
                            data[k] /= factor;
                        } else {
                            data[k] *= factor;
                        }
                    }
                }
                matrix.set_dtype("float64");
            }
            correction_factors = weights;
        }
    }

    // Replace the NaN values introduced by the weights and derive the bins
    // that hold no interaction at all.
    std::vector<double>& data = matrix.mutable_data();
    for (double& value : data) {
        if (std::isnan(value)) {
            value = 0.0;
        }
    }
    matrix.eliminate_zeros();

    const std::int64_t shape = std::min(matrix.rows(), matrix.cols());
    std::vector<char> used_as_column(static_cast<std::size_t>(shape), 0);
    for (const std::int32_t column : matrix.indices()) {
        if (column >= 0 && column < shape) {
            used_as_column[static_cast<std::size_t>(column)] = 1;
        }
    }
    nan_bins.clear();
    for (std::int64_t bin = 0; bin < shape; ++bin) {
        if (used_as_column[static_cast<std::size_t>(bin)] != 0) {
            continue;
        }
        const bool empty_row = matrix.indptr()[static_cast<std::size_t>(bin)] ==
                               matrix.indptr()[static_cast<std::size_t>(bin) + 1];
        if (empty_row) {
            nan_bins.push_back(bin);
        }
    }
}

}  // namespace

HiCMatrix HiCMatrix::load(const std::string& path, const Options& options) {
    HiCMatrix result;
    std::vector<CutInterval> cut_intervals;

    if (ends_with(path, ".h5")) {
        H5MatrixData loaded = read_hicexplorer_h5(path);
        result.matrix_ = std::move(loaded.matrix);
        cut_intervals = std::move(loaded.cut_intervals);
        result.nan_bins_ = std::move(loaded.nan_bins);
        result.correction_factors_ = std::move(loaded.correction_factors);
    } else {
        load_cool(path, options, result.matrix_, cut_intervals, result.nan_bins_,
                  result.correction_factors_);
    }

    if (options.fill_lower_triangle) {
        result.matrix_.symmetrize_in_place();
        if (options.materialize_full_matrix) {
            result.matrix_.materialize_full();
        }
    }
    result.bins_ = BinTable(std::move(cut_intervals));
    return result;
}

}  // namespace hicx
