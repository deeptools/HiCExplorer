// The reduction of compute_distance_mean (hicexplorer/hicPlotDistVsCounts.py
// :117-271), shared by the tool hicPlotDistVsCounts and the Python API
// (hicexplorer._hicx.hicPlotDistVsCounts).
//
// For the stored upper triangle entries of one unit (the whole matrix or one
// chromosome) in scipy's row major order: the bin distance of every entry from
// the bin starts as hicmatrix.HiCMatrix.fit_cut_intervals snaps them, the
// np.bincount sums and counts per distance, and the per distance loop of
// :219-261 with its diagonal lengths, its log records and the early stop after
// more than ten adjacent empty diagonals.

#ifndef HICX_DIST_VS_COUNTS_HPP
#define HICX_DIST_VS_COUNTS_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include "hicx/bins.hpp"

namespace hicx::dist_counts {

struct Entry {
    std::int32_t row;
    std::int32_t col;
    double value;
};

// HiCMatrix.fit_cut_intervals: the bin starts, snapped to the nearest
// multiple of the median bin size when more than 1 % of the bins deviate from
// it. `snapped` and `median` report whether that happened (hicmatrix logs it).
// Throws hicx::PythonError("ValueError", "need at least one array to
// concatenate") when no unit holds two bins.
[[nodiscard]] std::vector<std::int64_t> fitted_starts(const std::vector<CutInterval>& intervals,
                                                      bool* snapped = nullptr, std::int64_t* median = nullptr);

// The sizes of intervalListToIntervalTree(intervals)[1] whose names do not
// start with _ignore_: a name seen again later keeps its first position and
// takes the later range.
[[nodiscard]] std::vector<std::int64_t> unit_sizes(const std::vector<CutInterval>& intervals);

// One iteration of the per distance loop, k = bin distance + 1.
struct Distance {
    std::size_t k = 0;
    double sum = 0.0;
    // maxdepth is set and k == 0: mu and nb are NaN, nothing else happens
    bool skipped = false;
    // diagonal_length after the max() with the count of stored values. It is
    // a float (total_intra / 2) for k == 0 unless the count is larger.
    double diagonal_length = 0.0;
    bool diagonal_is_float = false;
    // the mean, NaN when diagonal_length is 0
    double mean = 0.0;
    bool defined = false;
    bool zero_value = false;  // logs "zero value for ..."
    bool stop = false;        // logs "skipping rest of ..." and ends the loop
    bool nan_sum = false;     // logs "nan value found for distance ..."
};

// The loop over entries[begin, end), whose rows and columns are relative to
// `offset`, of a unit of `unit_bins` bins with the given (possibly custom)
// cut intervals. Throws hicx::PythonError for a negative bin distance, as
// np.bincount refuses it.
[[nodiscard]] std::vector<Distance> scan(const std::vector<Entry>& entries, std::size_t begin, std::size_t end,
                                         std::int32_t offset, std::int64_t unit_bins,
                                         const std::vector<CutInterval>& intervals, std::int64_t binsize,
                                         bool maxdepth_set);

}  // namespace hicx::dist_counts

#endif  // HICX_DIST_VS_COUNTS_HPP
