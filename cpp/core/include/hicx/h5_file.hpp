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

// A matrix whose rows are produced on demand instead of being held.
//
// The overload of write_hicexplorer_h5 below exists for one case: a result
// that is dense, so that holding it as a CsrMatrix costs more than the dense
// block it came from. hicTransform --method pearson on Li_et_al_2015.h5 is
// that case. Its output has 123 M stored entries over 11,104 bins, 1.5 GB as
// CSR against a 1,221 MB budget (cpp/PLAN.md 4.5), while the whole computation
// needs nothing but the 40 MB sparse input and a row buffer. This is
// cpp/PLAN.md 4.4 rule 7, "a dense result is emitted row by row instead of
// being assembled in memory first", made available to a caller.
class DenseRowSource {
  public:
    virtual ~DenseRowSource() = default;
    [[nodiscard]] virtual std::int64_t rows() const = 0;
    // The dtype the values are written in, as in CsrMatrix::dtype().
    [[nodiscard]] virtual const std::string& dtype() const = 0;
    // Fills out[0 .. rows()-1] with row `i`. Called concurrently for different
    // rows, so an implementation must be free of shared mutable state, and the
    // values it produces must not depend on how the rows were distributed.
    virtual void fill_row(std::int64_t i, double* out) const = 0;
};

// write_hicexplorer_h5 over a DenseRowSource. `metadata.matrix` is ignored;
// everything else in it is written as by the other overload.
//
// The source is walked twice: once to count the stored entries, because a
// PyTables CArray has a fixed length that has to be known before the dataset
// is created, and once to write them. Both passes split the rows over
// `threads` workers by a fixed contiguous range, and the second pass emits in
// row order regardless of which worker produced a row, so the file is
// byte-identical for any thread count.
//
// The stored entries are the entries that are not exactly zero, which is what
// scipy's csr_matrix(dense) followed by hicmatrix's eliminate_zeros() leaves.
// A NaN is not zero, so it is stored, exactly as in the Python.
void write_hicexplorer_h5(const std::string& path, const MatrixData& metadata,
                          const DenseRowSource& source, int threads,
                          const H5SaveOptions& options = H5SaveOptions());

// The same result as a CsrMatrix, for a caller whose writer cannot stream.
// Only for a source whose result is small; the whole point of DenseRowSource
// is that this is what must not happen on a large one.
[[nodiscard]] CsrMatrix materialize_dense_row_source(const DenseRowSource& source,
                                                     int threads);

}  // namespace hicx

#endif  // HICX_H5_FILE_HPP
