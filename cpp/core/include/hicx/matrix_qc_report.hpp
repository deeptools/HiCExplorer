// The outputs of hicMatrixQC: a statistics table (tab separated, one column
// per sample) and an HTML report.
//
// The report follows hicQC.html after commit 103a2658: it needs no
// JavaScript, SVG or external file. Charts are HTML tables whose bars are
// cells with a percentage width and a bgcolor attribute, so they survive
// Galaxy's HTML sanitizer, which removes <script>, <style> and <svg>; the
// style sheet (inside an HTML comment) only sets the type and the spacing.
// Tables with more than three samples are split into blocks of three
// sample columns, which keeps every table readable at 400 px width.

#ifndef HICX_MATRIX_QC_REPORT_HPP
#define HICX_MATRIX_QC_REPORT_HPP

#include <string>
#include <utility>
#include <vector>

#include "hicx/matrix_qc.hpp"

namespace hicx::matrixqc {

// The statistics of one sample as (name, value) rows, in the order of the
// table. Values are text: integral counts without decimals, other numbers
// with 6 significant digits, "NA" for a statistic that does not exist.
[[nodiscard]] std::vector<std::pair<std::string, std::string>> flatten_stats(
    const MatrixStats& stats, const MapResolution* map_resolution);

// "statistic<TAB>sample 1<TAB>sample 2 ..." and one row per statistic; a
// sample without a statistic has an empty cell.
[[nodiscard]] std::string render_table(const std::vector<FileStats>& files);

[[nodiscard]] std::string render_html(const std::vector<FileStats>& files, const Options& options,
                                      const std::string& version);

}  // namespace hicx::matrixqc

#endif  // HICX_MATRIX_QC_REPORT_HPP
