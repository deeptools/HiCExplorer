#include "hicx/hic_matrix.hpp"

#include <algorithm>
#include <cmath>

#include "hicx/cool_adapter.hpp"
#include "hicx/h5_file.hpp"
#include "hicx/hic_adapter.hpp"

namespace hicx {

namespace {

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

HiCMatrix HiCMatrix::load(const std::string& path, const Options& options) {
    HiCMatrix result;
    std::vector<CutInterval> cut_intervals;

    // The same three-way branch as ToolMatrix::load (tool_matrix.cpp): a .hic
    // source, detected the same way (hic_adapter's is_hic_path), goes through
    // hic_adapter's read_hic with its own resolution/normalisation selector
    // in `path` rather than through the generic 'weight' column Options
    // above, which no hic2cool cool file has (see hic_adapter.hpp).
    if (is_hic_path(path)) {
        HicLoadResult loaded = read_hic(path);
        result.matrix_ = std::move(loaded.data.matrix);
        cut_intervals = std::move(loaded.data.cut_intervals);
        result.nan_bins_ = std::move(loaded.data.nan_bins);
        result.correction_factors_ = std::move(loaded.data.correction_factors);
    } else if (ends_with(path, ".h5")) {
        H5MatrixData loaded = read_hicexplorer_h5(path);
        result.matrix_ = std::move(loaded.matrix);
        cut_intervals = std::move(loaded.cut_intervals);
        result.nan_bins_ = std::move(loaded.nan_bins);
        result.correction_factors_ = std::move(loaded.correction_factors);
    } else {
        CoolLoadOptions cool_options;
        cool_options.apply_correction = options.apply_correction;
        cool_options.correction_factor_table = options.correction_factor_table;
        CoolLoadResult loaded = read_cool(path, cool_options);
        result.matrix_ = std::move(loaded.data.matrix);
        cut_intervals = std::move(loaded.data.cut_intervals);
        result.nan_bins_ = std::move(loaded.data.nan_bins);
        result.correction_factors_ = std::move(loaded.data.correction_factors);
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
