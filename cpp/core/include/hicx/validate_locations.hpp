// The computations of hicexplorer/hicValidateLocations.py: pandas.read_csv
// tables (hicx::TextTable), pybedtools' sort, merge and intersect -c
// (hicx::bedtools), the binning of applyBinning and the cool lookup of
// correlateCool. Used by the hicValidateLocations tool and by the public
// Python functions of the module (cpp/python/API_BINDINGS.md).
//
// Every BedTool.from_dataframe(...).<op>().to_dataframe() in the Python
// writes a temporary file and reads it back, so a column's dtype is
// re-inferred from the text at each step; TextTable::reparse does the same in
// memory.

#ifndef HICX_VALIDATE_LOCATIONS_HPP
#define HICX_VALIDATE_LOCATIONS_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hicx/bedtools_ops.hpp"
#include "hicx/text_table.hpp"

namespace hicx::validate_locations {

// BedTool.from_dataframe(df).sort().to_dataframe(...)
[[nodiscard]] TextTable sorted_table(const TextTable& table);

// The intervals of three columns of a table (chrom as written, start, end).
[[nodiscard]] bedtools::Intervals intervals_of(const TextTable& table, std::size_t chrom_col,
                                               std::size_t start_col, std::size_t end_col);

// A table of chrom (object), start and end (int64).
[[nodiscard]] TextTable table_of(const bedtools::Intervals& intervals);

// applyBinning step by step, for the debug records of the Python:
//   binned   (df[1] / bin).astype(int) * bin, ((df[2] / bin).astype(int) + 1) * bin
//   kept     the rows binned.drop_duplicates() keeps
//   deduped  binned restricted to kept
//   merged   BedTool.from_dataframe(deduped).merge().to_dataframe()   (merge only)
//   sorted   ...merge().sort().to_dataframe()                          (merge only)
struct BinningSteps {
    TextTable binned;
    std::vector<std::size_t> kept;
    TextTable deduped;
    TextTable merged;
    TextTable sorted;
};
[[nodiscard]] BinningSteps binning_steps(const TextTable& table, std::int64_t bin_size, bool merge);

// applyBinning's result: sorted for a merge, else deduped.
[[nodiscard]] TextTable apply_binning(const TextTable& table, std::int64_t bin_size, bool merge);

// 'add': 'chr' + column.astype(str); 'remove': column.str.lstrip('chr');
// anything else: unchanged.
void apply_chr_prefix(TextTable& table, std::size_t column, const std::string& mode);

// pd.read_csv(file, sep='\t', header=None)[[0, 1, 2]] with the prefix applied.
[[nodiscard]] TextTable read_three_column(const std::string& path, const std::string& prefix);

// correlateCool: cooler.Cooler(...).matrix(balance=False, sparse=True).fetch()
// for every loop, reduced to the only question the caller asks, whether the
// fetched block holds a stored pixel. Any exception inside the Python loop is
// swallowed and recorded as no match (hicValidateLocations.py:145-146), which
// covers a chromosome the cooler does not have and a coordinate past the end
// of a chromosome.
struct CoolPeaks {
    std::int64_t number_of_peaks = 0;
    std::vector<bool> matched;
};
[[nodiscard]] CoolPeaks correlate_cool(const std::string& path, const TextTable& loops);

// bedtools intersect -a A -b B -c of two tables (the first three columns of
// each), the count column appended to A, read back as to_dataframe() reads it.
[[nodiscard]] TextTable intersect_count_table(const TextTable& a, const TextTable& b);

}  // namespace hicx::validate_locations

#endif  // HICX_VALIDATE_LOCATIONS_HPP
