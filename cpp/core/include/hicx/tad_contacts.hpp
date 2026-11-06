// The computational core shared by hicDifferentialTAD and hicInterIntraTAD.
//
// Both Python tools were written from one template (hicDifferentialTAD.py:78-
// 278 and hicInterIntraTAD.py:73-268 are the same loop with a different body),
// so they share every quirk of how a TAD and its two neighbourhoods are cut
// out of the matrix. That geometry lives here once, and so do the dense block
// extraction and the Python number semantics the text output depends on.
// It is part of libhicx4 because the public Python API of both modules
// (hicexplorer._hicx, cpp/python/src/api) runs the same code as the tools.
//
// For every TAD i of a chromosome's list L (n entries) the Python builds:
//
//   intra   the TAD itself, rows and columns.
//   left    rows: the previous TAD, columns: this TAD.
//   right   rows: this TAD, columns: the next TAD.
//
// HiCExplorer 4 cuts all three from the bin range of each TAD (tad_bins in
// tad_contacts.cpp), the same range for the intra block of a TAD and for the
// inter blocks it takes part in, so the blocks of neighbouring TADs tile the
// matrix:
//
//  1. **The bin range of a TAD.** For cool the region cooler fetches for it,
//     bins floor(start / binsize) up to ceil(end / binsize); for h5
//     getRegionBinRange, whose end index is the bin *containing* `end`, used as
//     an exclusive bound. For the bin-aligned TADs of hicFindTADs the two agree.
//
//  2-4. The Python cut the inter blocks wrongly in three ways, fixed in
//     HiCExplorer 4 (project owner, 2026-11-06): the cool right block dropped
//     the last bin of the next TAD (the outer index was the literal -1,
//     hicDifferentialTAD.py:158-159); the last TAD of a chromosome took its
//     left block from the previous iteration's right index
//     (hicDifferentialTAD.py:162-167), an empty block on h5 (hicDifferentialTAD
//     reported NaN, hicInterIntraTAD divided by zero) and the wrong block on
//     cool; and `i - 1 > 0` (hicDifferentialTAD.py:180) gave the second TAD of a
//     two-TAD chromosome no left block. Every TAD but the first now has its left
//     block and every TAD but the last its right block, the same on cool and h5.
//
//  5. **The balancing weights are applied per region, not per matrix.**
//     hicmatrix's cool loader skips the correction for a region holding at
//     most one stored pixel, or whose weights are all NaN (cool.py:215-217),
//     and it zeroes NaN afterwards. None of the corpus matrices has a weight
//     column, so this branch is reproduced from the source, not from data.
//
// What is **not** reproduced, deliberately: the Python's output depends on
// --threads in three situations, all defects of its process partitioning
// (hicDifferentialTAD.py:352-451), pinned by the characterization tests:
// a chromosome with fewer TADs than processes makes the next chromosome lose
// its last TAD; the same situation after a larger chromosome writes stale
// result slots a second time; and exactly one TAD per process skips the last
// TAD's left test. cpp/OPTIMIZATION.md 3 requires the C++ output to be
// independent of the thread count, so this port always produces what the
// Python produces at --threads 1. On the full corpus TAD files the Python's
// output is identical at 1, 4, 11 and 16 processes, so the two agree there.

#ifndef HICX_TAD_CONTACTS_HPP
#define HICX_TAD_CONTACTS_HPP

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx::tads {

// --------------------------------------------------------------------------
// The TAD list

// One row of pd.read_csv(domains, sep='\t', header=None)[[0, 1, 2, 3, 4, 5]].
struct Domain {
    std::string chrom;
    std::int64_t start = 0;
    std::int64_t end = 0;
    // str() of every one of the six cells, which is how both tools print the
    // row back: `list(map(str, row))` of the Python objects `.values.tolist()`
    // returns, or the unicode array np.array(rows) makes of them, which calls
    // the same str().
    std::array<std::string, 6> text;
};

// Throws when the file cannot be read as pandas would read it or has fewer
// than six columns (the column selection raises KeyError in the Python).
[[nodiscard]] std::vector<Domain> read_domains(const std::string& path);

// Consecutive runs of the same chromosome name, in file order. A chromosome
// that reappears after another one starts a new run, exactly like the loop at
// hicDifferentialTAD.py:283-297.
[[nodiscard]] std::vector<std::vector<Domain>> group_by_chromosome(
    const std::vector<Domain>& domains);

// --------------------------------------------------------------------------
// The matrix

enum class Format { Cool, H5 };

struct ContactMatrix {
    Format format = Format::H5;
    // cool: the raw pixel table, upper triangle; h5: the stored matrix after
    // fillLowerTriangle. Column indices are sorted within every row.
    CsrMatrix matrix;
    BinTable bins;
    // cool only: chroms/name and chroms/length, for cooler's region checks.
    std::vector<std::string> chrom_names;
    std::vector<std::int64_t> chrom_lengths;
    // cool only: the 'weight' bin column when the file has one.
    std::optional<std::vector<double>> weights;
    // cool only: the 'bin-size' attribute, 0 for a variable bin size.
    std::int64_t bin_size = 0;
    // h5 only: the nan_bins hicmatrix loads with the matrix.
    std::vector<std::int64_t> nan_bins;
};

// The bins hicmatrix reports as nan_bins after loading the whole matrix, the
// "invalid" bins of hicDifferentialTAD --sharedMask (C++ only). h5: the
// file's nan_bins. cool: the bins that no nonzero pixel touches, in row or
// column, after the weights are applied the way a whole-file load applies
// them (cool.py:236-253). Sorted, without duplicates.
[[nodiscard]] std::vector<std::int64_t> invalid_bins(const ContactMatrix& matrix);

// Whether both matrices have the same bin table (chromosome, start, end).
[[nodiscard]] bool same_bins(const ContactMatrix& a, const ContactMatrix& b);

// Removes every stored pixel in a row or column of `bins`, so that those bins
// read as zero in every block cut from the matrix.
void mask_bins(ContactMatrix& matrix, const std::vector<std::int64_t>& bins);

// hm.hiCMatrix(path) for h5, and the state the per-region cool loads need.
// `is_cooler` is hicexplorer.utilities.check_cooler(path).
[[nodiscard]] ContactMatrix load_contact_matrix(const std::string& path, bool is_cooler);

// Sorts the column indices of every row, a no-op when they already are.
void sort_row_indices(CsrMatrix& matrix);

// --------------------------------------------------------------------------
// Geometry

// (start, stop) of Python's slice(start, stop).indices(length) for step 1,
// with stop raised to start, which is what scipy's _get_submatrix does with an
// inverted slice.
[[nodiscard]] std::pair<std::int64_t, std::int64_t> normalise_slice(std::int64_t start,
                                                                   std::int64_t stop,
                                                                   std::int64_t length);

// A rectangle of the whole matrix in global bin indices, half open.
struct Block {
    std::int64_t row_begin = 0;
    std::int64_t row_end = 0;
    std::int64_t col_begin = 0;
    std::int64_t col_end = 0;
    // cool only: whether the loaded region this block was sliced from had its
    // balancing weights applied (point 5 of the file comment).
    bool apply_correction = false;

    [[nodiscard]] std::int64_t rows() const noexcept { return row_end - row_begin; }
    [[nodiscard]] std::int64_t cols() const noexcept { return col_end - col_begin; }
    friend bool operator==(const Block&, const Block&) = default;
};

struct TadGeometry {
    Block intra;
    std::optional<Block> left;
    std::optional<Block> right;
};

// The three blocks of TAD `index` of `chromosome`. Throws std::runtime_error
// wherever the Python raises: an unknown chromosome, a region cooler rejects,
// or a position no bin covers (getRegionBinRange returns None and the caller
// indexes it).
[[nodiscard]] TadGeometry tad_geometry(const ContactMatrix& matrix,
                                       const std::vector<Domain>& chromosome,
                                       std::size_t index);

// --------------------------------------------------------------------------
// Values

// The dtype scipy carries for a block. Integer sums print without a decimal
// point, float32 sums print float32 reprs.
[[nodiscard]] DType block_dtype(const ContactMatrix& matrix, const Block& block);

struct DenseBlock {
    // Row major, exactly `.toarray().flatten()`.
    std::vector<double> values;
    std::int64_t rows = 0;
    std::int64_t cols = 0;
    // The block's `.nnz`: stored entries of the scipy slice.
    std::int64_t nnz = 0;
};

[[nodiscard]] DenseBlock extract_block(const ContactMatrix& matrix, const Block& block);

// A Python number as the two tools hold it: a plain int where they assign the
// literal 0, a numpy scalar where scipy computed it. The kind decides both the
// arithmetic and the text.
struct PyNumber {
    enum class Kind { PyInt, Int64, Float32, Float64 };
    Kind kind = Kind::PyInt;
    std::int64_t integer = 0;
    double real = 0.0;

    static PyNumber py_int(std::int64_t value) { return {Kind::PyInt, value, 0.0}; }
    static PyNumber int64(std::int64_t value) { return {Kind::Int64, value, 0.0}; }
    static PyNumber float32(double value) { return {Kind::Float32, 0, value}; }
    static PyNumber float64(double value) { return {Kind::Float64, 0, value}; }

    [[nodiscard]] bool is_integer() const noexcept {
        return kind == Kind::PyInt || kind == Kind::Int64;
    }
    [[nodiscard]] double as_double() const noexcept {
        return is_integer() ? static_cast<double>(integer) : real;
    }
    // str() of the value.
    [[nodiscard]] std::string str() const;
};

// a + b and a / b under numpy 1.26's scalar promotion rules. Division by a
// zero numpy scalar gives inf or NaN, as numpy does; the tools never divide
// two plain ints here.
[[nodiscard]] PyNumber py_add(const PyNumber& a, const PyNumber& b);
[[nodiscard]] PyNumber py_divide(const PyNumber& a, const PyNumber& b);

// block.sum() with scipy's evaluation order: every row summed sequentially in
// column order (csr_matvec against a vector of ones), then numpy's pairwise sum
// over the row sums. Integer blocks sum exactly into an int64.
[[nodiscard]] PyNumber block_sum(const DenseBlock& block, DType dtype);

// scipy.stats.ranksums(x, y) including the NaN propagation its
// _axis_nan_policy decorator adds: any NaN in either sample gives a NaN
// statistic and p-value. Two empty samples also give NaN.
struct RankTest {
    double statistic = 0.0;
    double pvalue = 0.0;
};
[[nodiscard]] RankTest rank_sum_test(std::span<const double> x, std::span<const double> y);

}  // namespace hicx::tads

#endif  // HICX_TAD_CONTACTS_HPP
