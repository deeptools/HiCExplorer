// The arithmetic of the public functions of hicexplorer/hicAggregateContacts.py
// (Python API, hicexplorer._hicx.hicAggregateContacts): the bin a large region
// maps to, the oriented submatrix and its sum in numpy's order, the averages
// over submatrices, and the outlier rows of get_outlier_indices. The tool
// hicAggregateContacts computes the same in cpp/tools/hicAggregateContacts.cpp.

#ifndef HICX_AGGREGATE_CONTACTS_HPP
#define HICX_AGGREGATE_CONTACTS_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace hicx::aggregate {

enum class LargeRegions { First, Last, Center };

// bin_id[0], bin_id[-1] or int(np.floor(np.mean(bin_id))) of a
// getRegionBinRange pair.
[[nodiscard]] std::int64_t pick_bin(std::int64_t first, std::int64_t last, LargeRegions operation);

// count_contacts' orientation of a submatrix from the two strands.
enum class Orientation { None, FlipLR, FlipUD, Transpose, Unchanged };
[[nodiscard]] Orientation orientation_of(const std::optional<std::string>& strand1,
                                         const std::optional<std::string>& strand2);

// The oriented w by w window (row major) and np.matrix.sum() of it: numpy
// sums a fliplr or flipud view in its logical order and the unflipped and the
// transposed matrix in their memory order, the window's.
struct Oriented {
    std::vector<double> values;
    double sum = 0.0;
};
[[nodiscard]] Oriented orient(const std::vector<double>& window, std::size_t w, Orientation orientation);

// np.median(values) of a flat copy: NaN when any value is NaN, else the
// middle value or the mean of the two middle values.
[[nodiscard]] double median(std::vector<double> values);

// compute_avg: np.median, np.mean or np.sum along axis 0 of `count` arrays of
// `cell` values each (row major, array after array). np.mean and np.sum add
// the arrays in order.
enum class Average { Median, Mean, Sum };
[[nodiscard]] std::vector<double> axis0(const std::vector<double>& values, std::size_t count, std::size_t cell,
                                        Average average);

// get_outlier_indices of a rows by cols array: nullopt when the scaled median
// absolute value is 0 (the Python returns None), else per row whether any
// value deviates from the median by more than max_deviation scaled medians.
[[nodiscard]] std::optional<std::vector<bool>> outlier_rows(const std::vector<double>& values, std::size_t rows,
                                                            std::size_t cols, double max_deviation);

}  // namespace hicx::aggregate

#endif  // HICX_AGGREGATE_CONTACTS_HPP
