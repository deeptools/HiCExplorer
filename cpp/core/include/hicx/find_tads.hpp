// The computations of hicexplorer/hicFindTADs.py: the TAD-separation score,
// its local minima, their delta and p-value, the hierarchical clustering of
// the boundaries and the text files the module writes.
//
// The hicFindTADs tool (cpp/tools/hicFindTADs.cpp) and the public Python API
// of hicexplorer.hicFindTADs (cpp/python/src/api/hicFindTADs.cpp) run this
// code. Two kinds of functions live here:
//
//  * the kernels the tool has always used (incremental_step_size,
//    drop_diagonal, mask_bins, truncate_to_distance, the cut weights,
//    peakdetect, delta_wrt_window); they moved here unchanged, so the tool's
//    output did not change;
//  * the general forms the Python API needs where the tool only ever calls
//    a special case: peakdetect with any lookahead and delta, the step sizes
//    over floats, hiCMatrix.diagflat(value=0) with scipy's dtype promotion,
//    and the functions the tool never calls (hierarchical_clustering).
//
// Behaviour of the Python that is reproduced deliberately (also listed in the
// tool) is documented at each function.

#ifndef HICX_FIND_TADS_HPP
#define HICX_FIND_TADS_HPP

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hicx/bins.hpp"
#include "hicx/matrix_data.hpp"
#include "hicx/sparse_matrix.hpp"

namespace hicx::find_tads {

// --------------------------------------------------------------------------
// Window sizes

// get_incremental_step_size for integers: min + int(step * k ** 1.5) for
// k = 0, 1, ..., skipping a repeat after k > 1 and stopping past max.
[[nodiscard]] std::vector<std::int64_t> incremental_step_size(std::int64_t min_window,
                                                              std::int64_t max_window,
                                                              std::int64_t start_step);

// The same over doubles, for float arguments: each entry is
// min + int(step * k ** 1.5) computed as Python does. The Python loops forever
// when step <= 0 and min <= max; this throws PythonError(ValueError) instead.
[[nodiscard]] std::vector<double> incremental_step_size_real(double min_window,
                                                             double max_window,
                                                             double start_step);

// --------------------------------------------------------------------------
// The matrix operations of hiCMatrix that hicFindTADs performs

// The tool's diagflat(value=0): the diagonal entries disappear.
void drop_diagonal(CsrMatrix& matrix);

// hiCMatrix.diagflat(value=0) as scipy computes it: matrix - diag(matrix) +
// dia_matrix(zeros). Every result that is exactly zero is dropped, so stored
// zeros anywhere vanish too; a NaN or infinite diagonal entry becomes NaN
// (x - x), and the result dtype is numpy's promotion of the matrix dtype with
// int64 (the dtype of np.repeat(0, M)).
[[nodiscard]] CsrMatrix diagflat_zero(const CsrMatrix& matrix);

// The state hiCMatrix keeps for masked bins: orig_bin_ids (the kept rows
// followed by the masked ones) and orig_cut_intervals in the same order.
struct MaskState {
    std::vector<std::int64_t> orig_bin_ids;
    std::vector<CutInterval> orig_cut_intervals;
};

// The tool's maskBins on a matrix that has never been masked.
MaskState mask_bins(MatrixData& data, const std::vector<std::int64_t>& bin_ids);
// The tool's restoreMaskedBins.
void restore_masked_bins(MatrixData& data, MaskState& state);

// triu(m, 0) - triu(m, limit), then eliminate_zeros: entries at or beyond
// `limit` cancel, except NaN, which cannot cancel and stays stored.
void truncate_to_distance(CsrMatrix& matrix, std::int64_t limit);

// --------------------------------------------------------------------------
// The TAD-separation score

// A dense block of the matrix, row major, exactly what
// matrix[r0:r1, c0:c1].todense() produces. `stored` is the .nnz of the slice.
struct DenseBlock {
    std::vector<double> values;
    std::size_t stored = 0;
    std::int64_t rows = 0;
    std::int64_t columns = 0;
};

// The block of rows [r0, r1) and columns [c0, c1); an empty range gives an
// empty block. The bounds must lie within the matrix.
void extract_block(const CsrMatrix& matrix, std::int64_t r0, std::int64_t r1,
                   std::int64_t c0, std::int64_t c1, DenseBlock& block);

// The state one TAD-score computation needs: the matrix, the bin table and the
// chromosome end positions (hiCMatrix.get_chromosome_sizes).
struct ScoreContext {
    const CsrMatrix* matrix = nullptr;
    const std::vector<CutInterval>* cut_intervals = nullptr;
    const BinTable* bins = nullptr;
    std::unordered_map<std::string, std::int64_t> chromosome_sizes;
};

// get_idx_of_bins_at_given_distance. Returns false where the Python's
// getRegionBinRange returns None and the caller ends up with a TypeError.
bool bins_at_distance(const ScoreContext& context, std::int64_t bin,
                      std::int64_t window_length, std::int64_t& left, std::int64_t& right);

// get_cut_weight(..., return_mean=True). `defined` is false where the Python
// returns None; a block without stored entries gives 0 (the Python's integer
// 0); otherwise the numpy mean of the dense block.
double cut_weight_mean(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch, bool& defined);

// get_cut_weight(..., return_mean=False): the flattened dense block. False
// where the Python returns None.
bool cut_weight_values(const ScoreContext& context, std::int64_t cut,
                       std::int64_t window_length, DenseBlock& scratch);

// --------------------------------------------------------------------------
// Local minima

struct Peak {
    std::int64_t position = 0;  // an index into y
    double value = 0.0;
};

// HicFindTads.peakdetect with delta 0 and x_axis None, as find_boundaries
// calls it. Throws on an infinite or NaN score (the Python's assert).
void peakdetect(const std::vector<double>& y, const std::vector<std::string>& chrom,
                std::int64_t lookahead, std::vector<Peak>& maxima, std::vector<Peak>& minima);

// HicFindTads.peakdetect in general. The Python iterates over
// y_axis[:-lookahead] (Python slice semantics, so lookahead 0 iterates over
// nothing and a negative lookahead over a prefix) and looks ahead through
// y_axis[index:index + lookahead].max(); positions are indices into y, the
// caller maps them to x_axis. `chrom` is compared by value; nullptr stands for
// chrom=None, which the Python subscripts and so raises TypeError as soon as
// the loop runs. `y_is_array` is false for a y_axis without a max() method
// (a list), which raises AttributeError at the first look ahead. Errors are
// PythonError of the type the Python raises.
// A peak as the Python appends it, [max_pos, max_y]: the position is the index
// of x_axis (-1 for None) and the value the index of y whose value max_y
// holds. The two can differ, because the Python resets max_y to the current
// value after a peak without resetting max_pos.
struct PeakAt {
    std::int64_t position = -1;
    std::int64_t value_index = -1;
};
struct PeakResult {
    std::vector<PeakAt> maxima;
    std::vector<PeakAt> minima;
};
[[nodiscard]] PeakResult peakdetect_general(const std::vector<double>& y,
                                            const std::vector<std::string>* chrom,
                                            std::int64_t lookahead, double delta,
                                            bool y_is_array, bool assertions,
                                            const std::string& y_type_name);

// delta_wrt_window (hicFindTADs.py:574-616): for every minimum not within
// `window_length` of a chromosome range border, the mean of
// matrix_avg[i - w : i + 3] together with matrix_avg[i + 4 : i + w] minus the
// score at i. That window includes the minimum itself and the two bins after
// it and skips i + 3, although the docstring says the minimum is excluded;
// reproduced. NaN close to a border. The ranges are built from the first
// index of every chromosome name in sorted name order (np.unique) plus
// len(chrom) - 1.
std::map<std::int64_t, double> delta_wrt_window(const std::vector<std::int64_t>& minima,
                                                const std::vector<double>& scores,
                                                const std::vector<std::string>& chrom,
                                                std::int64_t window_length = 10);

// The same, one value per entry of `minima` in its order (the Python's dict
// keeps the first insertion position and the last value of a repeated key,
// which for equal keys is the same value).
[[nodiscard]] std::vector<double> delta_wrt_window_values(const std::vector<std::int64_t>& minima,
                                                          const std::vector<double>& scores,
                                                          const std::vector<std::string>& chrom,
                                                          std::int64_t window_length);

// matrix.mean(axis=1) of a row major float64 matrix.
[[nodiscard]] std::vector<double> row_means(const std::vector<double>& values, std::size_t rows,
                                            std::size_t columns);

// --------------------------------------------------------------------------
// Statistics of find_boundaries and min_pvalue

// The Benjamini-Hochberg cutoff of min_pvalue: the largest p with
// p <= q * (i + 1) / n over the sorted p-values, nullopt when none qualifies
// (the Python then keeps its integer 0).
[[nodiscard]] std::optional<double> fdr_cutoff(std::vector<double> pvalues, double q);

// np.median and np.percentile(values, q) with linear interpolation.
[[nodiscard]] double numpy_median(std::vector<double> values);
[[nodiscard]] double numpy_percentile(std::vector<double> values, double q);

// --------------------------------------------------------------------------
// hierarchical_clustering

// The input of one chromosome as np.split hands it to the loop: the boundary
// starts and their values.
// An integer column (int64 in numpy) is sorted by numpy's int64 argsort, a
// float column by its float64 argsort; the two tie orders differ.
struct ClusterInput {
    std::vector<double> start;
    std::vector<double> value;
    bool start_integer = false;
    bool value_integer = false;
};

// One merge: z_value's (left, right, values[order_idx], size, x_left,
// x_right), with the value given as the index into the chromosome's value
// list (the caller keeps the original numbers).
struct ClusterMerge {
    std::int64_t left = 0;
    std::int64_t right = 0;
    std::size_t value_index = 0;  // into ClusterInput::value (values[i] is value[i + 1])
    std::int64_t size = 0;
    std::int64_t x_left = 0;
    std::int64_t x_right = 0;
};

// One clusters_cutoff level of one chromosome: the clusters as the indices
// of the domains whose start and end bound them, ordered by start.
struct ClusterRegion {
    std::size_t start_domain = 0;  // domains[min(set)][0] is start[start_domain]
    std::size_t end_domain = 0;    // domains[max(set)][1] is start[end_domain + 1]
};

struct ClusterResult {
    std::vector<ClusterMerge> merges;
    // One entry per cutoff that was crossed, in crossing order.
    std::vector<std::vector<ClusterRegion>> levels;
};

// The loop body of HicFindTads.hierarchical_clustering for one chromosome.
// `cutoffs` are clusters_cutoff after np.sort(np.unique(...))[::-1].
// Throws PythonError(AssertionError, "error") where the Python's assert fails
// (an empty chromosome piece).
[[nodiscard]] ClusterResult cluster_chromosome(const ClusterInput& input,
                                               const std::vector<double>& cutoffs);

// The pieces np.split(array, indices) produces: (begin, end) of each, Python
// slice semantics, so unsorted split points give empty or overlapping pieces.
[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> numpy_split_bounds(
    std::size_t length, const std::vector<std::int64_t>& indices);

// --------------------------------------------------------------------------
// Text

// '{:.<decimals>f}'.format(value) for a float: nan, inf and -inf as Python
// writes them (never "-nan").
[[nodiscard]] std::string format_fixed(double value, int decimals);

}  // namespace hicx::find_tads

#endif  // HICX_FIND_TADS_HPP
