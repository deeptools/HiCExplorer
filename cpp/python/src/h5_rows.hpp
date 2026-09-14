// Row range reads of a HiCExplorer h5 matrix for the hicx_matrix module.
//
// hicx::h5::File reads whole datasets only. Showing one region of a matrix
// needs the CSR rows of that region and nothing else, so this reader selects
// hyperslabs of /matrix/indices and /matrix/data, a bounded number of entries
// at a time. The indptr array (one integer per bin) and the bin intervals are
// read whole: they are the index, proportional to the bins and not to the
// stored pixels.

#ifndef HICX_PYTHON_H5_ROWS_HPP
#define HICX_PYTHON_H5_ROWS_HPP

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "hicx/hdf5_util.hpp"

namespace hicx::python {

class H5Rows {
  public:
    explicit H5Rows(const std::string& path);

    [[nodiscard]] std::int64_t nbins() const noexcept {
        return static_cast<std::int64_t>(indptr_.size()) - 1;
    }

    // Calls visit(row, column, value) for every stored entry of the rows
    // [row_first, row_last), in storage order.
    void for_each_entry(std::int64_t row_first, std::int64_t row_last,
                        const std::function<void(std::int64_t, std::int64_t, double)>& visit) const;

    // hicmatrix fills the lower triangle on load when tril(matrix, k=-1).sum()
    // is 0 (hiCMatrix.fillLowerTriangle), which is the case for every matrix
    // hicmatrix saves with pSymmetric. The answer depends on all stored
    // entries, so the first call scans the indices once, a bounded number of
    // entries at a time, and reads data only where an entry lies below the
    // diagonal. The result is cached.
    [[nodiscard]] bool fills_lower_triangle() const;

  private:
    class Dataset {
      public:
        Dataset() = default;
        Dataset(hid_t file, const std::string& path);
        void read(std::int64_t lo, std::int64_t hi, hid_t mem_type, void* out) const;

      private:
        h5::Handle dataset_;
        std::string path_;
    };

    h5::File file_;
    std::vector<std::int64_t> indptr_;
    Dataset indices_;
    Dataset data_;
    mutable std::optional<bool> fills_lower_;
};

}  // namespace hicx::python

#endif  // HICX_PYTHON_H5_ROWS_HPP
