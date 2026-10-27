// What the bindings of the TAD modules (hicFindTADs, hicMergeTADbins,
// hicDifferentialTAD, hicInterIntraTAD) share: the hiCMatrix methods the
// 3.7.6 functions called on a hicmatrix.HiCMatrix.hiCMatrix object, done in
// C++ on the object's attributes. Defined in hicFindTADs.cpp.
//
// The object keeps the attributes hicmatrix gives it (matrix, cut_intervals,
// nan_bins, correction_factors, distance_counts, orig_bin_ids,
// orig_cut_intervals, interval_trees, chrBinBoundaries, bin_size); every
// function here reads and writes those attributes exactly where the
// hicmatrix method does, so code that keeps using the object afterwards sees
// the state 3.7.6 left. The positions are looked up in the bin table of
// cut_intervals, which is what interval_trees is built from.

#ifndef HICX_PYTHON_API_HICFINDTADS_HPP
#define HICX_PYTHON_API_HICFINDTADS_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>

#include "hicx/bins.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx::python::hic {

namespace py = pybind11;

// hiCMatrix.intervalListToIntervalTree(intervals): (interval_trees,
// chrBinBoundaries), the dict of intervaltree.IntervalTree per chromosome and
// the OrderedDict chromosome -> (first bin, end bin).
[[nodiscard]] std::pair<py::object, py::object> interval_trees(py::handle intervals);

// hic.cut_intervals = intervals and the two structures rebuilt from them.
void set_intervals(py::handle hic, py::handle intervals);

// The bin table of hic.cut_intervals.
[[nodiscard]] hicx::BinTable bins(py::handle hic);

// hiCMatrix.getRegionBinRange(chrname, startpos, endpos): the (start bin, end
// bin) pair, or nullopt where the Python returns None (logging "Index
// error"). Raises the ValueError of an unknown chromosome.
[[nodiscard]] std::optional<std::pair<std::int64_t, std::int64_t>> region_bin_range(
    py::handle hic, const hicx::BinTable& bins, py::handle chrname, py::handle startpos,
    py::handle endpos);

// hiCMatrix.getBinSize(): the cached hic.bin_size, computed and cached when
// it is None (numpy.int64 for a single bin, int otherwise, as hicmatrix).
[[nodiscard]] py::object bin_size(py::handle hic);

// hiCMatrix.getBinPos(index): cut_intervals[index], ValueError past the end.
[[nodiscard]] py::object bin_pos(py::handle hic, py::handle index);

// hiCMatrix.restoreMaskedBins().
void restore_masked_bins(py::handle hic);

// hic.matrix as a csr matrix with sorted column indices.
[[nodiscard]] hicx::CsrMatrix matrix(py::handle hic);

// A scipy csr_matrix of the represented matrix (an upper triangle matrix is
// written out in full, as hicmatrix holds it).
[[nodiscard]] py::object csr_object(hicx::CsrMatrix matrix);

// open(path, 'w') and write(text), with open()'s errors.
void write_text_file(py::handle path, const std::string& text);

// The `exit(code)` of the 3.7.6 code (site's exit, None for no argument).
[[noreturn]] void python_exit(py::object code);

}  // namespace hicx::python::hic

#endif  // HICX_PYTHON_API_HICFINDTADS_HPP
