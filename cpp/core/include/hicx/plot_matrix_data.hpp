// The data computations of the public functions of
// hicexplorer/hicPlotMatrix.py (Python API, hicexplorer._hicx.hicPlotMatrix):
// tick labels, region strings and the --log/--log1p clearing of a dense
// chromosome block.

#ifndef HICX_PLOT_MATRIX_DATA_HPP
#define HICX_PLOT_MATRIX_DATA_HPP

#include <string>
#include <string_view>
#include <vector>

namespace hicx::plot_matrix {

// relabel_ticks: '{:.2f} ' of x / 1e6 with ' Mbp' on the second last label
// when the last tick is above 1.5e6, '{:.0f}' of x / 1e3 with ' Kbp' above
// 1500, else '{:.2f} ' of x with ' bp'. Throws std::out_of_range for fewer than
// two ticks (labels[-2]).
[[nodiscard]] std::vector<std::string> relabel_ticks(const std::vector<double>& ticks);

// translate_region's cleanup: ',' ';' '!' removed, '-' becomes ':', then the
// fields of split(':').
[[nodiscard]] std::vector<std::string> region_fields(std::string_view region);

// plotPerChr's clearing of a dense float64 block for --log and --log1p
// (hicPlotMatrix.py:500-517): the masks of the zeros, NaN and infinite values
// are taken first, then each is replaced in turn by np.nanmin of the values
// outside it; an empty selection (numpy's ValueError) stops the replacements.
struct LogClearing {
    bool any_nan_before = false;
    bool any_inf_before = false;
    bool failed = false;
    bool any_nan_after = false;
    bool any_inf_after = false;
};
[[nodiscard]] LogClearing clear_for_log(std::vector<double>& matrix);

// np.clip(values, low, high) as numpy 1.26 computes it, minimum(maximum(v,
// low), high), NaN propagating; a missing bound is not applied.
void clip(std::vector<double>& values, const double* low, const double* high);

}  // namespace hicx::plot_matrix

#endif  // HICX_PLOT_MATRIX_DATA_HPP
