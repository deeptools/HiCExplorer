// Chunked 2D range queries over a cooler's pixel table: cooler.core's
// DirectRangeQuery2D and FillLowerRangeQuery2D (cooler/core/_rangequery.py of
// cooler 0.10.2), which read a bounding box of the matrix one row span at a
// time.

#ifndef COOLERCPP_RANGEQUERY_HPP
#define COOLERCPP_RANGEQUERY_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>

#include "coolercpp/api.hpp"
#include "coolercpp/table.hpp"

namespace coolercpp {

// In Python the engines take a CSRReader built from an open pixel group and
// the bin1_offset index:
//
//   reader = CSRReader(h5["pixels"], h5["indexes/bin1_offset"][:])
//   engine = DirectRangeQuery2D(reader, field, bbox, chunksize, return_index)
//
// coolercpp builds the reader from a Cooler; the object keeps the file open
// until it is destroyed.
class RangeQuery2D {
  public:
    enum class Kind {
        // DirectRangeQuery2D: the pixels exactly as stored.
        Direct,
        // FillLowerRangeQuery2D: for a symmetric-upper cooler, the implicit
        // lower triangle inside the bounding box is generated.
        FillLower,
    };

    // bbox is (row_start, row_stop, col_start, col_stop), half open.
    // chunksize is the rough number of stored pixel rows read per chunk.
    RangeQuery2D(const Cooler& clr, Kind kind, std::string field,
                 std::array<std::int64_t, 4> bbox, std::int64_t chunksize,
                 bool return_index = false);

    // engine.n_chunks
    [[nodiscard]] std::size_t n_chunks() const noexcept;
    // frame_slice_from_dict(engine.get_chunk(i), field): the pixels of one
    // task, with the pixel table index when return_index is set. Chunks come
    // in the order cooler concatenates them. IndexError past the last chunk.
    [[nodiscard]] Table get_chunk(std::size_t i) const;
    // engine.to_frame(), engine.to_sparse_matrix(), engine.to_array()
    [[nodiscard]] Table to_frame() const;
    [[nodiscard]] SparseMatrix to_sparse_matrix() const;
    [[nodiscard]] DenseMatrix to_array() const;

  private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};

// One dataset of a cooler group read slice by slice into buffers the caller
// owns, for whole-table reads that need no range query. In Python:
//
//   with clr.open("r") as grp:
//       dset = grp[path]                 # e.g. "pixels/bin2_id"
//       out[:] = dset[lo:hi]             # out preallocated with dtype T
//
// The object keeps the file open until it is destroyed. HDF5 converts the
// stored type to T; values outside the range of T saturate, where numpy's
// assignment would wrap them. KeyError when the dataset does not exist,
// TypeError when it is not numeric, IndexError when [lo, hi) is not inside
// the dataset (Python slicing would clip), ValueError when out.size() is not
// hi - lo.
class DatasetReader {
  public:
    DatasetReader(const Cooler& clr, const std::string& path);

    // len(dset)
    [[nodiscard]] std::int64_t size() const noexcept;
    void read_into(std::int64_t lo, std::int64_t hi, std::span<std::int32_t> out) const;
    void read_into(std::int64_t lo, std::int64_t hi, std::span<std::int64_t> out) const;
    void read_into(std::int64_t lo, std::int64_t hi, std::span<double> out) const;

  private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};

// cooler.fileops.is_cooler: whether a URI names a cooler group (its "format"
// attribute is "HDF5::Cooler"). False when the file is missing or not HDF5;
// KeyError when the file is HDF5 but the group does not exist.
[[nodiscard]] bool is_cooler(const std::string& uri);

}  // namespace coolercpp

#endif  // COOLERCPP_RANGEQUERY_HPP
