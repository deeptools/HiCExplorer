// The C++ implementation behind the public Python module hicexplorer.utilities.
//
// Every function here reproduces one function of hicexplorer/utilities.py of
// HiCExplorer 3.7.6 (upstream c2ad8630), including its defects, because the
// Python names stay public and their results must not change
// (cpp/python/API_BINDINGS.md). The Python wrapper only converts arguments;
// the arithmetic, the control flow and the error conditions live here.
//
// numpy semantics that the results depend on, reproduced as numpy 1.26 has
// them:
//  * np.sum is the pairwise summation of hicx::npy::pairwise_sum, in the dtype
//    of the data (float32 data sums in float32, int32 data in int64);
//    np.add.reduceat, which scipy's csr.sum(axis=1) uses, is the first value
//    plus the pairwise sum of the rest of the segment, in the data's own
//    dtype (int32 wraps in int32);
//  * a float64 to integer cast of NaN, of an infinity or of a value out of
//    range gives the most negative integer (x86 cvttsd2si), cast_from_double;
//  * scalar arithmetic between numpy scalars and Python numbers promotes by
//    type (np.float32 * 1.5 is float64), array-with-scalar arithmetic by value
//    (a float32 array divided by a float64 scalar stays float32).
//
// Defects of the Python that are reproduced, not fixed:
//  * scipy's nonzero() drops explicitly stored zeros, but the Python then
//    indexes the unfiltered data array with positions in the filtered list.
//    With explicit zeros the values and the distances are misaligned
//    (expected_interactions_in_distance, expected_interactions_non_zero,
//    obs_exp_matrix_non_zero) or the lengths disagree and numpy raises
//    (expected_interactions, compute_zscore, obs_exp_matrix with pDistance).
//  * the expected arrays are accumulated at abs(i - j) and read at
//    ceil(abs(i - j) / 2) by obs_exp_matrix and obs_exp_matrix_lieberman.
//  * compute_zscore subtracts the expected value from zero, not from the
//    observed value.

#ifndef HICX_UTILITIES_HPP
#define HICX_UTILITIES_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hicx::utilities {

// --------------------------------------------------------------------------
// numpy casts and Python arithmetic

// numpy's unsafe cast of a float64 to T (setitem, astype).
template <class T>
[[nodiscard]] T cast_from_double(double value) noexcept;

// Python's int % int (floor modulo, result has the sign of the divisor).
[[nodiscard]] std::int64_t python_mod(std::int64_t a, std::int64_t b);
// Python's float % float.
[[nodiscard]] double python_float_mod(double a, double b) noexcept;
// Python's int > float, exact for every int64.
[[nodiscard]] bool python_greater(std::int64_t a, double b) noexcept;

// The numpy messages the kernels raise.
[[nodiscard]] std::string index_out_of_bounds_message(std::int64_t index, std::int64_t size);
[[nodiscard]] std::string boolean_index_mismatch_message(std::size_t dimension,
                                                         std::size_t boolean_dimension);
[[nodiscard]] std::string broadcast_message(std::size_t a, std::size_t b);

// --------------------------------------------------------------------------
// text and files

// readBed(pBedFile): (chrom, start, end) of every line not starting with '#',
// taken from the tab separated fields or, with fewer than three of those, from
// the whitespace separated ones. `lines` are the file's readlines().
// Throws PythonError ValueError for a line with fewer than three fields.
[[nodiscard]] std::vector<std::vector<std::string>> read_bed(const std::vector<std::string>& lines);

// readTargetBed(pBedFile). `genes` lists every gene in the order it first
// appears, which is the key order of all three returned dicts, and
// `positions[k]` belongs to genes[k]. A line with exactly three tab separated
// fields reuses the gene of the previous line, as the Python's unbound
// variable does; on the first line it is an UnboundLocalError.
struct TargetBed {
    struct Positions {
        std::string chromosome;
        std::vector<std::string> starts;
        std::vector<std::string> ends;
    };
    std::vector<std::string> genes;
    std::vector<Positions> positions;
};
[[nodiscard]] TargetBed read_target_bed(const std::vector<std::string>& lines);

// open(path, 'w').close(): 0, or the errno of the failure.
[[nodiscard]] int truncate_for_writing(const std::string& path) noexcept;

// genomicRegion(string): nullopt for a string that is only whitespace.
// Throws PythonError argparse.ArgumentTypeError when nothing remains after
// removing the separators.
[[nodiscard]] std::optional<std::string> genomic_region(std::string_view text);

// change_chrom_names for a str.
[[nodiscard]] std::string change_chrom_name(std::string_view chrom);

// in_units(pBasePosition) after float().
[[nodiscard]] std::string in_units(double base_position);

// remove_non_ascii, that is unidecode 1.3.8's unidecode(text) with
// errors='ignore'. The surrogate code points of the text, for which unidecode
// warns, are appended to `surrogates`.
[[nodiscard]] std::string unidecode(std::string_view text,
                                    std::vector<char32_t>* surrogates = nullptr);

// The first two bytes of a gzip file.
[[nodiscard]] bool gzip_magic(std::string_view first_bytes) noexcept;

// True when unidecode returns its argument unchanged (all ASCII).
[[nodiscard]] bool is_ascii(std::string_view text) noexcept;

// getUserRegion, split in the two steps between which the Python parses the
// tile size. `start` and `end` are the parsed region fields, nullopt when the
// field is absent. Throws PythonError NameError for a start after the end or a
// negative start.
struct UserRegion {
    std::int64_t start = 0;
    std::int64_t end = 0;
    // The end is the chromosome size: the Python returns that object.
    bool end_is_chrom_size = false;
    // The end was moved to a tile boundary, which makes it a new number of
    // the type of the arithmetic.
    bool end_tiled = false;
    std::int64_t chunk_size = 0;
};
[[nodiscard]] UserRegion user_region_bounds(std::int64_t chrom_size,
                                            std::optional<std::int64_t> start,
                                            std::optional<std::int64_t> end,
                                            std::string_view region_string);
// The limit max_chunk_size is an int or a float in Python.
struct ChunkLimit {
    bool is_float = true;
    std::int64_t integer = 0;
    double floating = 1e6;
};
// Applies the tile size (0 or nullopt: none) and the chunk limit.
void user_region_finish(UserRegion& region, std::optional<std::int64_t> tile_size,
                        const ChunkLimit& limit);

// enlarge_bins on the coordinates of the intervals. `same_chrom[i]` is
// chrom_i == chrom_{i+1} as Python compares the two objects.
enum class Coordinate : std::uint8_t { Original, Zero, Middle };
struct EnlargedBin {
    std::int64_t start = 0;
    std::int64_t end = 0;
    // A Middle start takes the type of this bin's original start, a Middle end
    // the type of the next bin's original start (start_next of the Python).
    Coordinate start_kind = Coordinate::Original;
    Coordinate end_kind = Coordinate::Original;
    // The chromosome object of the rebuilt tuple comes from this bin.
    std::size_t chrom_from = 0;
    // The Python assigned a new tuple to this position.
    bool rebuilt = false;
};
[[nodiscard]] std::vector<EnlargedBin> enlarge_bins(const std::vector<std::int64_t>& starts,
                                                    const std::vector<std::int64_t>& ends,
                                                    const std::vector<bool>& same_chrom);

// --------------------------------------------------------------------------
// statistics of flat arrays (C order)

// np.median of the values; the result type is float for float input, double
// otherwise. NaN anywhere gives NaN; an empty array gives NaN.
template <class T>
[[nodiscard]] double median(const std::vector<T>& values);

// np.mean, same result types.
template <class T>
[[nodiscard]] double mean(const std::vector<T>& values);

// remove_outliers: the positions (in C order) of the values it keeps, or
// nullopt when the MAD is not positive and the data is returned unchanged.
template <class T>
[[nodiscard]] std::optional<std::vector<std::size_t>> remove_outliers_keep(
    const std::vector<T>& values, double max_deviation);

// --------------------------------------------------------------------------
// sparse matrices

// The stored entries of a scipy.sparse matrix in tocoo() order.
struct SparseEntries {
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
};

// csr/csc/coo .nonzero(): the coordinates of the stored values != 0.
struct Nonzero {
    std::vector<std::int64_t> row;
    std::vector<std::int64_t> col;
    [[nodiscard]] std::vector<std::int64_t> distances() const;
    [[nodiscard]] std::size_t size() const noexcept { return row.size(); }
};
template <class T>
[[nodiscard]] Nonzero nonzero(const SparseEntries& entries, const T* data);

// expected_interactions_in_distance; `length` and `count` are already the
// np.int32 values. Throws PythonError IndexError for a distance >= n.
template <class T>
[[nodiscard]] std::vector<double> expected_interactions_in_distance(
    std::int64_t n, const std::vector<std::int64_t>& distance, const T* data,
    std::int32_t length, std::int32_t count);

// expected_interactions_non_zero.
template <class T>
[[nodiscard]] std::vector<double> expected_interactions_non_zero(
    std::int64_t n, const std::vector<std::int64_t>& distance, const T* data);

// expected_interactions_thread without the queue: sums[i] = np.sum(data[distance == i])
// for i in [low, high). Throws PythonError IndexError as numpy does.
template <class T>
[[nodiscard]] std::vector<double> distance_sums(std::int64_t size,
                                                const std::vector<std::int64_t>& distance,
                                                const T* data, std::size_t data_size,
                                                std::int64_t low, std::int64_t high);

// expected_interactions with pThreads=None; nullopt where the Python returns
// None (no non-zero entry).
template <class T>
[[nodiscard]] std::optional<std::vector<double>> expected_interactions(
    std::int64_t n, const std::vector<std::int64_t>& distance, const T* data,
    std::size_t data_size);

// compute_zscore on the data array in place. `depth` is pDepth (nullopt for
// None); `depth_is_float` reproduces range() rejecting a float. `expected`
// is what expected_interactions returned.
template <class T>
void compute_zscore(std::int64_t n, const std::vector<std::int64_t>& distance, T* data,
                    std::size_t data_size, std::optional<double> depth, bool depth_is_float,
                    const std::optional<std::vector<double>>& expected);

// The new data of obs_exp_matrix and obs_exp_matrix_lieberman: the data cast
// to float32, divided in float64 by expected[ceil(distance / 2)] (broadcast as
// numpy does), NaN and infinity replaced by `non_finite` and cast back to T.
template <class T>
[[nodiscard]] std::vector<T> obs_exp_divide(const T* data, std::size_t data_size,
                                            const std::vector<double>& expected,
                                            const std::vector<std::int64_t>& distance,
                                            double non_finite);

// pSubmatrix.data[distance >= limit] = 0 (obs_exp_matrix's pDistance).
template <class T>
void zero_from_distance(const std::vector<std::int64_t>& distance, T* data,
                        std::size_t data_size, double limit);

// The float32 data obs_exp_matrix_non_zero leaves before eliminate_zeros.
// `indptr` is the compressed structure for csr (by row) or csc (by column);
// `csr` selects which, because scipy sums rows differently for the two.
template <class T>
[[nodiscard]] std::vector<float> obs_exp_non_zero_values(
    const SparseEntries& entries, const std::vector<std::int64_t>& indptr, bool csr,
    const T* data, std::size_t data_size, const Nonzero& nonzero,
    const std::vector<double>& expected, bool ligation_factor, bool to_epsilon);

// eliminate_zeros of a compressed matrix: the kept positions, and indptr
// rewritten in place.
template <class T>
[[nodiscard]] std::vector<std::size_t> eliminate_zeros(std::vector<std::int64_t>& indptr,
                                                       const T* data);

}  // namespace hicx::utilities

#endif  // HICX_UTILITIES_HPP
