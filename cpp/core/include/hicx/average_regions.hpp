// The viewpoint windows of hicAverageRegions (calculateViewpointRange and
// calculateViewpointRangeBins of hicexplorer/hicAverageRegions.py), shared by
// the tool and the Python API (hicexplorer._hicx.hicAverageRegions).

#ifndef HICX_AVERAGE_REGIONS_HPP
#define HICX_AVERAGE_REGIONS_HPP

#include <cstdint>

namespace hicx::average_regions {

// int(float(start) + ((float(end) - float(start)) / 2)): the centre of a
// viewpoint, truncated towards zero.
[[nodiscard]] std::int64_t center_of(double start, double end);

struct Window {
    std::int64_t start = 0;
    std::int64_t end = 0;
    bool start_out_of_range = false;
    bool end_out_of_range = false;
};

// calculateViewpointRange: the genomic window anchor - up .. anchor + down,
// clamped to 0 and, when it reaches past max_length (the end of the last bin
// of the chromosome), to max_length - 1.
[[nodiscard]] Window viewpoint_range(std::int64_t anchor, std::int64_t up, std::int64_t down,
                                     std::int64_t max_length);

// calculateViewpointRangeBins: the bin window index - up .. index + down,
// clamped to the chromosome's bins [first_bin, last_bin].
[[nodiscard]] Window viewpoint_range_bins(std::int64_t index, std::int64_t up, std::int64_t down,
                                          std::int64_t first_bin, std::int64_t last_bin);

}  // namespace hicx::average_regions

#endif  // HICX_AVERAGE_REGIONS_HPP
