// The four textual contact matrix formats hicConvertFormat handles.
//
//   homer          a dense, usually gzipped, tab separated table with one
//                  header row of bin names (hicmatrix/lib/homer.py)
//   ginteractions  one line per stored upper triangle entry, seven columns,
//                  write only because the Python load is a stub
//                  (hicmatrix/lib/ginteractions.py:14-15)
//   hicpro         a triplet file with one based bin ids plus a separate bed
//                  file holding the bin table (hicmatrix/lib/hicpro.py)
//   2D-text        seven columns of two genomic intervals and a value, read
//                  inline by hicConvertFormat.py:155-223 and by nothing else
//
// Numbers are written the way Python writes them, which is str() of a numpy
// scalar: an integer matrix prints '0', a float64 matrix prints '0.0' and
// '1.3e-05'. The dtype of the matrix therefore decides the text, and
// value_repr is the one place that knows it.

#ifndef HICX_TEXT_FORMATS_HPP
#define HICX_TEXT_FORMATS_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx {

// str() of one matrix value, in the dtype the matrix carries.
[[nodiscard]] std::string value_repr(double value, const std::string& dtype);

// triu(m, k=0).maximum(triu(m, k=0).T), which is how hicConvertFormat
// symmetrises before writing homer and ginteractions
// (hicConvertFormat.py:258-260). It is not the same as adding the transpose:
// maximum leaves the diagonal alone where a sum would double it, and it
// clamps a negative off diagonal value to zero rather than mirroring it.
[[nodiscard]] CsrMatrix maximum_with_transpose(const CsrMatrix& matrix);

// triu(maximum_with_transpose(m), k=0), built without the mirror.
//
// Ginteractions.save takes the upper triangle of the symmetric matrix
// hicConvertFormat hands it, so the mirror is computed and thrown away again.
// The triangle is derivable from the input directly: the cell (i, j) with
// i < j is max(m[i][j], 0) and the diagonal is unchanged, so this returns the
// same entries at half the peak memory. On Li_et_al_2015.h5, 1.66 M stored
// pixels, that is the difference between a 157 MB and a 41 MB peak.
[[nodiscard]] CsrMatrix upper_triangle_after_maximum(const CsrMatrix& matrix);

// --------------------------------------------------------------------------
// homer

// Reads a homer table, gzipped or not; the magic bytes decide, as
// hicmatrix.utilities.opener does. The bin size is taken from the difference
// between the first two column names, so a table with fewer than two bins
// cannot be read, and a chromosome name containing a '-' makes the name split
// ambiguous and is rejected, both exactly as in the Python.
[[nodiscard]] MatrixData read_homer(const std::string& path);

// Writes the gzipped dense table. Every row of the matrix is materialised in
// turn and released again, so the memory cost is one row, not one matrix.
void write_homer(const std::string& path, const MatrixData& data);

// --------------------------------------------------------------------------
// ginteractions

// Writes `path + ".tsv"`, leaving `path` itself alone. That is what the Python
// does, so --outFileName names a file that is never written.
void write_ginteractions(const std::string& path, const MatrixData& data);

// --------------------------------------------------------------------------
// hicpro

[[nodiscard]] MatrixData read_hicpro(const std::string& matrix_path,
                                     const std::string& bed_path);

// The triplet file and the bed file. Unlike every other writer this one emits
// all stored entries, not the upper triangle, because Hicpro.save ignores
// pSymmetric.
void write_hicpro(const std::string& matrix_path, const std::string& bed_path,
                  const MatrixData& data);

// --------------------------------------------------------------------------
// 2D text

// A two column "name<TAB>size" file, in file order. A repeated name keeps its
// first position and takes its last size, which is what assigning into an
// OrderedDict does.
[[nodiscard]] std::vector<std::pair<std::string, std::int64_t>> read_chromosome_sizes(
    const std::string& path);

// hicConvertFormat.py:169-223. The bin table is generated from the chromosome
// sizes at the requested resolution and the file only supplies values.
//
// Two behaviours of the Python are reproduced rather than fixed. A line whose
// interval reaches into the next bin sets *two* cells, because
// getRegionBinRange returns a (start bin, end bin) pair and the pair is used
// as a numpy fancy index, so `matrix[(a, b), (c, d)] = v` assigns to (a, c)
// and (b, d). And a value is assigned, not accumulated, so a repeated pair
// keeps the last value in the file and a value of zero removes the cell.
[[nodiscard]] MatrixData read_two_dimensional_text(
    const std::string& path,
    const std::vector<std::pair<std::string, std::int64_t>>& chromosome_sizes,
    std::int64_t resolution);

}  // namespace hicx

#endif  // HICX_TEXT_FORMATS_HPP
