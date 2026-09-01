#include "hicx/h5_file.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "hicx/hdf5_util.hpp"

namespace hicx {

bool is_hicexplorer_h5(const std::string& path) {
    if (!h5::is_hdf5(path)) {
        return false;
    }
    try {
        const h5::File file(path);
        return file.exists("/matrix/data") && file.exists("/matrix/indices") &&
               file.exists("/matrix/indptr") && file.exists("/intervals/chr_list");
    } catch (const h5::Error&) {
        return false;
    }
}

H5MatrixData read_hicexplorer_h5(const std::string& path) {
    const h5::File file(path);

    const std::vector<std::int64_t> shape = file.read_int64("/matrix/shape");
    if (shape.size() != 2) {
        throw h5::Error("unexpected /matrix/shape in " + path);
    }
    const std::string dtype = file.dataset_dtype("/matrix/data");
    std::vector<double> data = file.read_doubles("/matrix/data");
    std::vector<std::int32_t> indices = file.read_int32("/matrix/indices");
    std::vector<std::int64_t> indptr = file.read_int64("/matrix/indptr");

    H5MatrixData result;
    result.matrix = CsrMatrix(shape[0], shape[1], std::move(indptr), std::move(indices),
                              std::move(data), dtype);

    const std::vector<std::string> chroms = file.read_strings("/intervals/chr_list");
    const std::vector<std::int64_t> starts = file.read_int64("/intervals/start_list");
    const std::vector<std::int64_t> ends = file.read_int64("/intervals/end_list");
    // extra_list is normally float64, but matrices written by hicFindTADs
    // store text there (for example the z-score matrices under
    // test_data/find_TADs). Both have to be readable.
    const bool extra_is_text = file.dataset_dtype("/intervals/extra_list") == "string";
    std::vector<double> extra;
    std::vector<std::string> extra_text;
    if (extra_is_text) {
        extra_text = file.read_strings("/intervals/extra_list");
    } else {
        extra = file.read_doubles("/intervals/extra_list");
    }
    const std::size_t extra_size = extra_is_text ? extra_text.size() : extra.size();
    if (chroms.size() != starts.size() || chroms.size() != ends.size() ||
        chroms.size() != extra_size) {
        throw h5::Error("the interval lists of " + path + " have different lengths");
    }
    result.cut_intervals.reserve(chroms.size());
    for (std::size_t i = 0; i < chroms.size(); ++i) {
        CutInterval interval{chroms[i], starts[i], ends[i], 0.0, {}};
        if (extra_is_text) {
            interval.extra = std::numeric_limits<double>::quiet_NaN();
            interval.extra_text = extra_text[i];
        } else {
            interval.extra = extra[i];
        }
        result.cut_intervals.push_back(std::move(interval));
    }
    if (static_cast<std::int64_t>(result.cut_intervals.size()) != shape[0]) {
        throw h5::Error("Error loading matrix. Length of bin intervals (" +
                        std::to_string(result.cut_intervals.size()) +
                        ") is different than the size of the matrix (" +
                        std::to_string(shape[0]) + ")");
    }

    if (file.exists("/nan_bins")) {
        result.nan_bins = file.read_int64("/nan_bins");
    }

    if (file.exists("/correction_factors")) {
        std::vector<double> factors = file.read_doubles("/correction_factors");
        if (static_cast<std::int64_t>(factors.size()) != shape[0]) {
            throw h5::Error(
                "Error loading matrix. Length of correction factors does not"
                "match size of matrix");
        }
        for (double& value : factors) {
            if (std::isnan(value) || std::isinf(value)) {
                value = 0.0;
            }
        }
        result.correction_factors = std::move(factors);
    }

    // The Python reader has a copy and paste defect here: when
    // /distance_counts exists it reads /correction_factors instead. It is
    // reproduced so that the two implementations agree bit for bit.
    if (file.exists("/distance_counts") && file.exists("/correction_factors")) {
        result.distance_counts = file.read_doubles("/correction_factors");
    }

    return result;
}

}  // namespace hicx
