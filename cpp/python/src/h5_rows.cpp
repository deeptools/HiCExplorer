#include "h5_rows.hpp"

#include <algorithm>

namespace hicx::python {

namespace {

// Entries read per hyperslab: 8 MB of indices plus 8 MB of values.
constexpr std::int64_t kChunkEntries = std::int64_t{1} << 20;

}  // namespace

H5Rows::Dataset::Dataset(hid_t file, const std::string& path)
    : dataset_(H5Dopen2(file, path.c_str(), H5P_DEFAULT), h5::Handle::Kind::Dataset), path_(path) {
    if (!dataset_.valid()) {
        throw h5::Error("cannot open dataset " + path);
    }
}

void H5Rows::Dataset::read(std::int64_t lo, std::int64_t hi, hid_t mem_type, void* out) const {
    if (hi <= lo) {
        return;
    }
    const h5::Handle space(H5Dget_space(dataset_.get()), h5::Handle::Kind::DataSpace);
    if (!space.valid()) {
        throw h5::Error("cannot get the dataspace of " + path_);
    }
    const hsize_t start = static_cast<hsize_t>(lo);
    const hsize_t count = static_cast<hsize_t>(hi - lo);
    if (H5Sselect_hyperslab(space.get(), H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0) {
        throw h5::Error("cannot select entries of " + path_);
    }
    const h5::Handle memory(H5Screate_simple(1, &count, nullptr), h5::Handle::Kind::DataSpace);
    if (H5Dread(dataset_.get(), mem_type, memory.get(), space.get(), H5P_DEFAULT, out) < 0) {
        throw h5::Error("cannot read entries of " + path_);
    }
}

H5Rows::H5Rows(const std::string& path) : file_(path) {
    indptr_ = file_.read_int64("/matrix/indptr");
    if (indptr_.empty()) {
        throw h5::Error(path + ": /matrix/indptr is empty");
    }
    indices_ = Dataset(file_.id(), "/matrix/indices");
    data_ = Dataset(file_.id(), "/matrix/data");
}

void H5Rows::for_each_entry(std::int64_t row_first, std::int64_t row_last,
                            const std::function<void(std::int64_t, std::int64_t, double)>& visit) const {
    if (row_last <= row_first) {
        return;
    }
    const std::int64_t end = indptr_[static_cast<std::size_t>(row_last)];
    std::int64_t row = row_first;
    std::vector<std::int64_t> columns;
    std::vector<double> values;
    for (std::int64_t lo = indptr_[static_cast<std::size_t>(row_first)]; lo < end; lo += kChunkEntries) {
        const std::int64_t hi = std::min(end, lo + kChunkEntries);
        columns.resize(static_cast<std::size_t>(hi - lo));
        values.resize(static_cast<std::size_t>(hi - lo));
        indices_.read(lo, hi, H5T_NATIVE_INT64, columns.data());
        data_.read(lo, hi, H5T_NATIVE_DOUBLE, values.data());
        for (std::int64_t k = lo; k < hi; ++k) {
            while (indptr_[static_cast<std::size_t>(row + 1)] <= k) {
                ++row;
            }
            const auto i = static_cast<std::size_t>(k - lo);
            visit(row, columns[i], values[i]);
        }
    }
}

bool H5Rows::fills_lower_triangle() const {
    if (fills_lower_.has_value()) {
        return *fills_lower_;
    }
    const std::int64_t end = indptr_.back();
    std::int64_t row = 0;
    double lower_sum = 0.0;
    std::vector<std::int64_t> columns;
    std::vector<double> values;
    for (std::int64_t lo = indptr_.front(); lo < end; lo += kChunkEntries) {
        const std::int64_t hi = std::min(end, lo + kChunkEntries);
        columns.resize(static_cast<std::size_t>(hi - lo));
        indices_.read(lo, hi, H5T_NATIVE_INT64, columns.data());
        bool values_read = false;
        for (std::int64_t k = lo; k < hi; ++k) {
            while (indptr_[static_cast<std::size_t>(row + 1)] <= k) {
                ++row;
            }
            const auto i = static_cast<std::size_t>(k - lo);
            if (columns[i] < row) {
                if (!values_read) {
                    values.resize(static_cast<std::size_t>(hi - lo));
                    data_.read(lo, hi, H5T_NATIVE_DOUBLE, values.data());
                    values_read = true;
                }
                lower_sum += values[i];
            }
        }
    }
    fills_lower_ = lower_sum == 0.0;
    return *fills_lower_;
}

}  // namespace hicx::python
