// ENT3C's entropy-based similarity of two contact matrices. Reference:
// Lainscsek X, Taher L. ENT3C: an entropy-based similarity measure for Hi-C
// and micro-C derived contact matrices. NAR Genom Bioinform. 2024;6(3):lqae076,
// doi: 10.1093/nargab/lqae076.
//
// Written from the method description of the paper. The ENT3C program (2.2.2,
// PyPI) was run as a black box on test matrices to settle what the paper
// leaves open; this implementation reproduces its per-window entropies,
// missing ones included, and its similarity scores to within 1e-14.
//
// Per chromosome:
//  1. Bins whose row is empty in any of the compared matrices are removed
//     from all of them; N is the number of bins left.
//  2. The submatrix size is n = floor(N / c + 0.5) for the chromosome split
//     c, unless n is given. Windows of n x n start every phi bins along the
//     diagonal, Phi = 1 + floor((N - n) / phi) of them; while Phi exceeds
//     Phi_max, phi is increased by one.
//  3. In each window, every entry is replaced by its natural logarithm, a
//     zero by a missing value. A bin whose column of logarithms sums to 0,
//     missing values left out (a column of zeros and ones), is removed from
//     the window, row and column, leaving k bins; the remaining missing values
//     become the smallest remaining value. The k x k Pearson correlation
//     matrix P of the columns, with 0 for a column that is still constant and
//     1 on the diagonal, divided by k, is a density matrix rho, and the
//     window's entropy is the von Neumann entropy -sum(lambda ln lambda) over
//     the eigenvalues lambda > 0 of rho. A window with fewer than two bins
//     left, or without any contact, has no entropy (NaN).
//  4. The entropy values of the windows, in order, form the matrix's entropy
//     signal S, and the similarity Q of two matrices is the Pearson
//     correlation of their signals over the windows with an entropy in both.
//
// When every window of a chromosome has two bins, each entropy is 0 up to
// rounding and Q correlates rounding noise; it then differs between this
// implementation and ENT3C, and means nothing in either. With the default
// split of 7 this happens for chromosomes of 11 to 17 non-empty bins.
//
// The paper averages Q over chromosomes for a pair of samples; hicCorrelate
// --method ent3c reports that mean, over the chromosomes with a finite Q.

#ifndef HICX_ENT3C_HPP
#define HICX_ENT3C_HPP

#include <cstdint>
#include <vector>

#include "hicx/hicspector.hpp"

namespace hicx::ent3c {

struct Options {
    // n; 0 derives it from chromosome_split.
    std::int64_t submatrix_size = 0;
    // c, ENT3C's CHRSPLIT; the paper's default is 7.
    std::int64_t chromosome_split = 7;
    // phi, the shift between consecutive windows in bins.
    std::int64_t window_shift = 1;
    // Phi_max, ENT3C's PHI_MAX; 0 means no limit.
    std::int64_t max_windows = 1000;
    // Windows are independent; each is computed by one thread with one BLAS
    // thread, so the result does not depend on the number of threads.
    int threads = 1;
};

// The upper triangle of one chromosome, diagonal included, local to
// [0, nbins). A position stored twice is not expected.
using hicx::hicspector::ChromosomeEntries;

struct Windows {
    std::int64_t size = 0;   // n
    std::int64_t count = 0;  // Phi
    std::int64_t shift = 0;  // phi after the Phi_max adjustment
};

// Step 2 for a chromosome of `bins` bins left after step 1. count is 0 when
// n exceeds the number of bins. Throws std::invalid_argument on options that
// define no window (a split or shift below 1, a negative size or limit).
Windows windows(std::int64_t bins, const Options& options);

// Step 1: true for every bin whose row (both triangles) holds no contact.
std::vector<bool> empty_bins(const ChromosomeEntries& entries);

// Step 3 for one window: `counts` is the n x n symmetric window, row major.
double window_entropy(const std::vector<double>& counts, std::int64_t n);

// Steps 3 and 4 for one matrix: the entropy signal over `keep`, the bins
// left after step 1 in ascending order.
std::vector<double> entropy_signal(const ChromosomeEntries& entries,
                                   const std::vector<std::int64_t>& keep, const Windows& windows,
                                   int threads);

// The Pearson correlation of two signals over the windows where neither is
// NaN, as numpy.corrcoef computes it: NaN if the signals differ in length,
// fewer than two windows are left, or one of them is constant there.
double similarity(const std::vector<double>& s1, const std::vector<double>& s2);

}  // namespace hicx::ent3c

#endif  // HICX_ENT3C_HPP
