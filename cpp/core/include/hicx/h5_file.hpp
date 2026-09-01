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
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// The reader used to declare its own payload type; it is the shared one now.
using H5MatrixData = MatrixData;

[[nodiscard]] bool is_hicexplorer_h5(const std::string& path);

[[nodiscard]] H5MatrixData read_hicexplorer_h5(const std::string& path);

struct H5SaveOptions {
    // hicmatrix.lib.H5.save pSymmetric: store triu(matrix, k=0) instead of the
    // stored entries as they are.
    bool symmetric = true;
};

// Port of hicmatrix.lib.H5.save (hicmatrix/lib/h5.py:91).
//
// The file is written with the blosc filter at complevel 5 with shuffle, which
// is the PyTables pipeline, but the chunk shape is h5py's guess rather than
// PyTables' undocumented heuristic. cpp/PLAN.md 2.5 declares h5 output as
// value identical (L3) for exactly this reason.
//
// The matrix is streamed to disk one block at a time and, when it holds the
// full symmetric form, the upper triangle is selected on the fly, so no second
// copy of the matrix and no staged triangle is ever allocated. `data` is read
// but not modified.
void write_hicexplorer_h5(const std::string& path, const MatrixData& data,
                          const H5SaveOptions& options = H5SaveOptions());

}  // namespace hicx

#endif  // HICX_H5_FILE_HPP
