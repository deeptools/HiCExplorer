// Minimal .npy and .npz reader and writer (cpp/PLAN.md 2.4).
//
// npz is not a matrix format in the hicmatrix sense. It is what
// scipy.sparse.save_npz writes and scipy.sparse.load_npz reads, and in
// HiCExplorer it is the interface between hicAverageRegions, which produces
// it, and hicPlotAverageRegions, which consumes it. So the writer's contract
// is that scipy reads the file back, not that it looks like any particular
// matrix on disk.
//
// The container is a ZIP archive holding one `.npy` entry per array. `.npy`
// version 1.0 is a 6 byte magic, a 2 byte version, a 2 byte little endian
// header length, and then a Python dict literal padded with spaces so that the
// whole header is a multiple of 64 bytes and ends in a newline:
//
//     \x93NUMPY\x01\x00\x76\x00{'descr': '<f8', 'fortran_order': False,
//     'shape': (191,), }<spaces>\n
//
// followed by the raw little endian element bytes.
//
// Byte identity with numpy is close but not reached, and the reason is worth
// recording. numpy's `_savez` gives every ZIP entry the fixed timestamp
// 1980-01-01, which is zipfile.ZipInfo's default and not the current time, so
// the file is at least reproducible. What differs is the deflate stream: the
// reference compresses with the CPython interpreter's stock zlib and this
// writer with the zlib in the conda prefix, and two zlib builds do not
// generally emit the same bytes for the same input. The npz comparator
// therefore loads both files and compares the arrays, which is the same
// decision cpp/PLAN.md 2.5 takes for h5.

#ifndef HICX_NPZ_FILE_HPP
#define HICX_NPZ_FILE_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace hicx::npz {

// One array of an npz archive. `dtype` is a numpy dtype string as it appears
// in the .npy header, for example "<f8", "<i4" or "|S3"; `data` holds the
// element bytes in that dtype, little endian. An empty `shape` is a zero
// dimensional array, which is how scipy stores the format tag.
struct Array {
    std::string name;  // without the .npy suffix
    std::string dtype;
    std::vector<std::int64_t> shape;
    std::string data;
};

// The .npy header bytes for one array, exposed so that a unit test can pin it
// against numpy's without writing a file.
[[nodiscard]] std::string npy_header(const std::string& dtype,
                                     const std::vector<std::int64_t>& shape);

// Writes the arrays as a ZIP archive, in the given order, deflated when
// `compressed` (which is scipy.sparse.save_npz's default) and stored
// otherwise.
void write_npz(const std::string& path, const std::vector<Array>& arrays,
               bool compressed = true);

// Reads an npz archive written by numpy or by write_npz. Only the .npy
// features numpy itself writes are supported: version 1.0 and 2.0 headers,
// little endian numeric dtypes and fixed width byte strings, C order. A
// pickled object array is rejected rather than mis-parsed.
[[nodiscard]] std::vector<Array> read_npz(const std::string& path);

// scipy.sparse.save_npz for a CSR matrix, including its array order
// (indices, indptr, format, shape, data) and its dtypes: int32 indices and
// row offsets, an int64 shape pair, the 3 byte tag b'csr', and float64 values.
void save_csr_npz(const std::string& path, std::int64_t rows, std::int64_t cols,
                  const std::vector<std::int32_t>& indptr,
                  const std::vector<std::int32_t>& indices,
                  const std::vector<double>& data);

}  // namespace hicx::npz

#endif  // HICX_NPZ_FILE_HPP
