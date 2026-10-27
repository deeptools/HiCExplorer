// The C++ behind the public Python names of hicexplorer.lib.buildMatrixMethods
// and hicexplorer.hicBuildMatrixMicroC of HiCExplorer 3.7.6 that the
// hicBuildMatrix machinery (hicx/build_matrix.hpp, hicx/build_matrix_run.hpp)
// does not already provide in the form the Python functions have
// (cpp/python/API_BINDINGS.md). No Python headers.
//
// Every function reproduces its Python counterpart, including its defects;
// those are marked QUIRK with the line of buildMatrixMethods.py.

#ifndef HICX_BUILD_MATRIX_API_HPP
#define HICX_BUILD_MATRIX_API_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hicx::build_matrix_api {

// --- paths -------------------------------------------------------------------

// posixpath.normpath, abspath (relative to `cwd`), dirname and join.
[[nodiscard]] std::string normpath(const std::string& path);
[[nodiscard]] std::string abspath(const std::string& path, const std::string& cwd);
[[nodiscard]] std::string dirname(const std::string& path);
[[nodiscard]] std::string join(const std::string& a, const std::string& b);

// Tester.__init__: root is $HIC_TEST_DATA_DIR when that is set and not empty,
// else the directory `levels` dirname()s above abspath(module_file) plus
// "/test/test_data/"; bam_file_1 is join(root, "hic.bam").
struct TesterPaths {
    std::string root;
    std::string bam_file_1;
};
[[nodiscard]] TesterPaths tester_paths(const std::optional<std::string>& test_data_dir,
                                       const std::string& module_file, const std::string& cwd,
                                       int levels);

// --- ReadPositionMatrix ------------------------------------------------------

// ReadPositionMatrix.is_duplicated's key: "{}-{}" of the two chromosome
// names in the order `chrom1_first` says (chrom1 < chrom2), then "-{}-{}" of
// the two starts in the order `start1_first` says (start1 < start2). The
// arguments are the str() of each value.
[[nodiscard]] std::string read_position_key(const std::string& chrom1, const std::string& chrom2,
                                            bool chrom1_first, const std::string& start1,
                                            const std::string& start2, bool start1_first);

// --- get_bins ----------------------------------------------------------------

// range(start, stop, step) of Python: the values, or nullopt for step 0
// (ValueError "range() arg 3 must not be zero").
[[nodiscard]] std::optional<std::vector<std::int64_t>> python_range(std::int64_t start, std::int64_t stop,
                                                                    std::int64_t step);

// --- get_chrom_sizes -----------------------------------------------------------

// list(OrderedDict(zip(names, values)).items()) as index pairs: for every
// distinct name in the order of its first appearance, (index of that first
// appearance, index of its last value).
[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> ordered_dict_items(
    const std::vector<std::string>& names);

// --- check_dangling_end ----------------------------------------------------------

// str.upper() of ASCII text.
[[nodiscard]] std::string ascii_upper(std::string text);
[[nodiscard]] bool starts_with(const std::string& text, const std::string& prefix);
[[nodiscard]] bool ends_with(const std::string& text, const std::string& suffix);

// --- get_supplementary_alignment ---------------------------------------------------

// len(tag.split(";")[0:-1]): the number of ';' in the SA tag.
[[nodiscard]] std::size_t supplementary_count(const std::string& sa_tag);

// --- get_correct_map ----------------------------------------------------------------

// One CIGAR operation, (op, length) of pysam's cigartuples.
using CigarOperation = std::pair<std::int64_t, std::int64_t>;

// The number of query bases before the first M (pysam.CMATCH, 0) of the
// operations as given (the caller reverses them for a reverse read).
[[nodiscard]] std::int64_t bases_before_match(const std::vector<CigarOperation>& cigar);

// first_mapped.index(min(first_mapped)): the first index of the minimum.
[[nodiscard]] std::size_t first_index_of_minimum(const std::vector<std::int64_t>& values);

// --- enlarge_bins ----------------------------------------------------------------------

// buildMatrixMethods.enlarge_bins (a different function from
// utilities.enlarge_bins) on the coordinates of the bins, as the sequence of
// tuple assignments it makes, so that a caller can rebuild the Python objects
// and stop where the Python raises the KeyError of a missing chromosome size.
//
//   chrom_start = True
//   for idx in range(len(bins) - 1):
//       chrom, start, end = bins[idx]; chrom_next, start_next, end_next = bins[idx + 1]
//       if chr_start: start = 0; chr_start = False        # QUIRK: only once
//       if chrom == chrom_next and end != start_next:
//           middle = start_next - int((start_next - end) / 2)
//           bins[idx] = (chrom, start, middle); bins[idx + 1] = (chrom, middle, end_next)
//       if chrom != chrom_next:
//           bins[idx] = (chrom, start, sizes[chrom]); bins[idx + 1] = (chrom_next, 0, end_next)
//   bins[-1] = (chrom, start, sizes[chrom]) of the last bin
//
// QUIRK (buildMatrixMethods.py:438-442): chr_start is never set back to True,
// so only the first bin of the first chromosome is extended to 0 by that
// branch, and only when one of the two assignments writes it back; the first
// bin of every other chromosome starts at 0 through the chromosome change.
//
// A field of a rebuilt tuple: the original start or end of a bin, the int 0,
// a middle (whose Python type is the type of the original start of
// `type_from`, start_next of the Python), or the size of the chromosome of
// the tuple.
struct EnlargeField {
    enum class Kind : std::uint8_t { OriginalStart, OriginalEnd, Zero, Middle, Size } kind = Kind::OriginalStart;
    std::int64_t value = 0;
    std::size_t from = 0;  // the bin whose original field it is, or type_from of a Middle
};
struct EnlargeAssignment {
    std::size_t index = 0;
    // The chromosome object of the tuple is the original chromosome of bin
    // chrom_from: the same-chromosome branch hands bin idx's object on to
    // bin idx + 1.
    std::size_t chrom_from = 0;
    EnlargeField start;
    EnlargeField end;
};
// `same_chrom[i]` is chrom_i == chrom_{i+1} as Python compares the objects.
[[nodiscard]] std::vector<EnlargeAssignment> enlarge_bins(const std::vector<std::int64_t>& starts,
                                                          const std::vector<std::int64_t>& ends,
                                                          const std::vector<bool>& same_chrom);

}  // namespace hicx::build_matrix_api

#endif  // HICX_BUILD_MATRIX_API_HPP
