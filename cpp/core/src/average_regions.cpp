#include "hicx/average_regions.hpp"

namespace hicx::average_regions {

std::int64_t center_of(double start, double end) {
    return static_cast<std::int64_t>(start + ((end - start) / 2.0));
}

Window viewpoint_range(std::int64_t anchor, std::int64_t up, std::int64_t down,
                       std::int64_t max_length) {
    Window window{anchor - up, anchor + down, false, false};
    if (window.start < 0) {
        window.start = 0;
        window.start_out_of_range = true;
    }
    if (window.end > max_length) {
        // -1 keeps getRegionBinRange inside the last bin
        window.end = max_length - 1;
        window.end_out_of_range = true;
    }
    return window;
}

Window viewpoint_range_bins(std::int64_t index, std::int64_t up, std::int64_t down,
                            std::int64_t first_bin, std::int64_t last_bin) {
    Window window{index - up, index + down, false, false};
    if (window.start < first_bin) {
        window.start = first_bin;
        window.start_out_of_range = true;
    }
    if (window.end > last_bin) {
        window.end = last_bin;
        window.end_out_of_range = true;
    }
    return window;
}

}  // namespace hicx::average_regions
