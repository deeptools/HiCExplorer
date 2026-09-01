// Reader for the native HiCExplorer h5 format.
//
// The file is written by PyTables and holds
//   /matrix/{data,indices,indptr,shape}   the CSR arrays of the upper triangle
//   /intervals/{chr_list,start_list,end_list,extra_list}
//   /nan_bins                             optional
//   /correction_factors                   optional
//   /distance_counts                      optional
// All datasets are blosc compressed, see hdf5_util.hpp.

#ifndef HICX_H5_FILE_HPP
#define HICX_H5_FILE_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

struct H5MatrixData {
    CsrMatrix matrix;
    std::vector<CutInterval> cut_intervals;
    std::vector<std::int64_t> nan_bins;
    std::optional<std::vector<double>> correction_factors;
    std::optional<std::vector<double>> distance_counts;
};

[[nodiscard]] bool is_hicexplorer_h5(const std::string& path);

[[nodiscard]] H5MatrixData read_hicexplorer_h5(const std::string& path);

}  // namespace hicx

#endif  // HICX_H5_FILE_HPP
