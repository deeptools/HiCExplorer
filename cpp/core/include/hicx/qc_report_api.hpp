// The arithmetic of hicexplorer.hicPrepareQCreport's make_figure_* functions
// (HiCExplorer 3.7.6, upstream c2ad8630) on the columns of a QC table: the
// counts divided by 1e6, by a denominator column, or by the sum of several
// columns. The binding reads the DataFrame, calls these, and builds the
// DataFrames the Python functions return.
//
// pandas divides integers and floats as float64 numbers; a text cell makes
// numpy raise TypeError naming the operand types, which is reproduced.

#ifndef HICX_QC_REPORT_API_HPP
#define HICX_QC_REPORT_API_HPP

#include <string>
#include <vector>

#include "hicx/qc_html.hpp"

namespace hicx::qcapi {

using Column = std::vector<QcStylerCell>;

// column / scalar for every cell, in the column's order. Throws TypeError
// ("unsupported operand type(s) for /: 'str' and 'float'") for a text cell.
[[nodiscard]] std::vector<double> divide_by_scalar(const Column& column, double scalar);

// numerators / denominators cell by cell (the two columns have as many cells).
[[nodiscard]] std::vector<double> divide_columns(const Column& numerators,
                                                 const Column& denominators);

// The per row sum of the columns: the integer sum when every column holds
// integers (so that the division is exact), the float sum otherwise. The
// returned cells are integers or floats.
[[nodiscard]] Column sum_columns(const std::vector<Column>& columns);

}  // namespace hicx::qcapi

#endif  // HICX_QC_REPORT_API_HPP
