// The computational core of hicDetectStripes, separated from the tool so that
// the unit tests (cpp/tests/test_detect_stripes.cpp) can reach it directly.
//
// hicDetectStripes is a PLAN.md tier 9 feature (section 9.3): HiCExplorer 3.7
// has no stripe caller, so there is no Python behaviour to reproduce and no
// "quirk" to preserve. The algorithm below is free (PLAN.md 5.0 item 3) and
// is documented on its own terms.
//
// Method, in one paragraph. A stripe is an elongated run of enriched contacts
// that starts at an anchor bin on the diagonal and extends away from it along
// one axis: a horizontal stripe holds its row fixed and extends across
// columns, a vertical stripe holds its column fixed and extends down rows.
// For every anchor and orientation, the tool tests a grid of candidate
// lengths against a *local* background: the same distance range read off the
// neighbouring anchors a few bins away (not the genome-wide expected value),
// which is what PLAN.md 9.3 asks the enrichment to be measured against. A
// cheap one-sample z-test against the local background's mean and standard
// deviation preselects candidates in O(1) per anchor and length using a
// precomputed running sum along the distance axis; the survivors, after
// keeping only the best-scoring length per anchor and suppressing anchors
// that lie within one non-maximum window of a stronger call, get a rigorous
// two-sample Wilcoxon rank-sum p-value (hicx::stats::ranksums, the same
// primitive hicDetectLoops' donut test uses) against the pooled local
// background pixels. Benjamini-Hochberg FDR (hicx::stats::
// benjamini_hochberg_adjusted) is then applied across every surviving
// candidate, genome-wide and both orientations together, and calls below the
// requested q-value are kept.
//
// Why this design and not Stripenn's or a HiCCUPS-style approach: Stripenn
// (Yoon et al. 2022) fits a Gaussian mixture to seed candidate stripes from an
// image-derivative pseudo-image and expands them by a decay-fit boundary
// search; Zebra (Zhang lab, unpublished) walks a 1-D score profile along the
// diagonal. Both are single anchor-scan methods at heart: a per-anchor score
// along one axis, thresholded and then locally maximised, exactly the shape
// used here. This port keeps the scan (which is what both external methods
// reduce to) and swaps their imaging or curve-fitting seed step for a direct
// statistical test against the local background, because that is what
// PLAN.md 9.3's recovery rule and precision gate are defined against, and
// because it reuses already-verified, deterministic primitives
// (hicx::parallel::parallel_for and hicx::stats::ranksums) rather than adding
// an OpenCV or curve-fitting dependency for a single tool.
//
// Threading and determinism (cpp/OPTIMIZATION.md 3). The preselection scan is
// one independent computation per (anchor, orientation) pair, threaded with
// hicx::parallel_for and written into a preallocated slot per pair; nothing
// is accumulated across pairs. The final rank-sum pass is the same, over the
// (far smaller) list of survivors. Both stages are independent of the thread
// count and of --threads.

#ifndef HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
#define HICX_TOOLS_DETECT_STRIPES_IMPL_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "hicx/sparse_matrix.hpp"

namespace hicx::stripes {

// A dense per-chromosome band: for anchor bin `a` in [0, n_bins) and distance
// `d` in [1, max_distance], band.raw[a * max_distance + (d - 1)] is the raw
// count at the pixel `d` bins away from `a`, and band.obs_exp the same pixel
// divided by the genome (chromosome)-wide expected count at that distance.
// Horizontal orientation reads the pixel at (row = a, col = a + d); vertical
// orientation reads (row = a - d, col = a), materialised as its own band so
// that both orientations are scanned by the same code.
struct Band {
    std::int64_t n_bins = 0;
    std::int64_t max_distance = 0;
    std::vector<double> raw;
    std::vector<double> obs_exp;

    [[nodiscard]] double raw_at(std::int64_t anchor, std::int64_t distance) const {
        return raw[static_cast<std::size_t>(anchor * max_distance + (distance - 1))];
    }
    [[nodiscard]] double obs_exp_at(std::int64_t anchor, std::int64_t distance) const {
        return obs_exp[static_cast<std::size_t>(anchor * max_distance + (distance - 1))];
    }
};

// Builds the horizontal band from a CSR matrix holding the chromosome's upper
// triangle (row < col, as hicDetectLoops' loader leaves it), restricted to
// entries with column - row <= max_distance. `expected_out`, when given,
// receives the per-distance expected value (the mean raw count over every
// position at that distance, stored or not).
[[nodiscard]] Band build_horizontal_band(const CsrMatrix& upper_triangle, std::int64_t n_bins,
                                         std::int64_t max_distance,
                                         std::vector<double>* expected_out = nullptr);

// The vertical band derived from the horizontal one: vertical.raw_at(c, d) ==
// horizontal.raw_at(c - d, d) for c - d >= 0, and 0 otherwise (no such pixel).
[[nodiscard]] Band build_vertical_band(const Band& horizontal);

// Running sums of a band along the distance axis: row_cum[a * max_distance +
// (L - 1)] is the sum of obs_exp_at(a, 1..L) (and raw_cum the sum of
// raw_at(a, 1..L)), so the mean pixel value of a length-L candidate at anchor
// a is row_cum[...] / L, in O(1) once this is built.
struct RunningSums {
    std::int64_t n_bins = 0;
    std::int64_t max_distance = 0;
    std::vector<double> obs_exp_cum;
    std::vector<double> raw_cum;

    [[nodiscard]] double obs_exp_mean(std::int64_t anchor, std::int64_t length) const {
        return obs_exp_cum[static_cast<std::size_t>(anchor * max_distance + (length - 1))] /
               static_cast<double>(length);
    }
    [[nodiscard]] double raw_mean(std::int64_t anchor, std::int64_t length) const {
        return raw_cum[static_cast<std::size_t>(anchor * max_distance + (length - 1))] /
               static_cast<double>(length);
    }
};

[[nodiscard]] RunningSums build_running_sums(const Band& band);

struct DetectOptions {
    // Candidate lengths to test, in bins, ascending, each >= 1.
    std::vector<std::int64_t> length_grid_bins;
    // Flanking anchors on each side of the background window.
    std::int64_t background_window_bins = 15;
    // Anchors immediately next to the stripe body that are excluded from the
    // background, so that the background never overlaps the stripe itself or
    // a neighbour close enough to share its signal.
    std::int64_t background_gap_bins = 2;
    // The minimum obs/exp enrichment of the stripe body's mean over the local
    // background's mean for a length to be considered at all.
    double min_obs_exp = 1.5;
    // The minimum one-sample z-score (stripe mean against the background
    // anchors' distribution of per-anchor means) for preselection.
    double preselect_z = 2.0;
    // The minimum mean raw count over the stripe body (peakInteractionsThreshold
    // analogue): filters out stripes with too little data to trust.
    double min_raw_count = 1.0;
    // Non-maximum suppression window, in bins: among same-orientation
    // candidates whose anchors are within this many bins of each other, only
    // the best-scoring one (by preselection z) survives.
    std::int64_t merge_window_bins = 5;
    // Benjamini-Hochberg q-value: candidates at or below this after
    // multiple-testing correction are kept in the final call set.
    double fdr_q = 0.05;
};

struct Candidate {
    std::int64_t anchor = 0;
    std::int64_t length_bins = 0;
    bool vertical = false;
    double enrichment = 0.0;   // stripe mean obs/exp / local background mean obs/exp
    double zscore = 0.0;       // preselection z-score at the chosen length
    double pvalue = 1.0;       // final rank-sum p-value
    double qvalue = 1.0;       // Benjamini-Hochberg adjusted, genome-wide
};

// Stage 1: for every anchor, the best-scoring (anchor, orientation) candidate
// over the length grid that clears min_obs_exp, preselect_z and
// min_raw_count, or nothing when none of the grid's lengths clears all three.
// horizontal and vertical must have the same n_bins and max_distance.
[[nodiscard]] std::vector<Candidate> preselect(const Band& horizontal, const Band& vertical,
                                               const RunningSums& horizontal_sums,
                                               const RunningSums& vertical_sums,
                                               const DetectOptions& options,
                                               unsigned int threads);

// Stage 2: non-maximum suppression within options.merge_window_bins, per
// orientation, keeping the candidate with the largest zscore in each cluster
// of anchors that are pairwise within the window of some kept candidate
// (a sweep in anchor order, matching hicDetectLoops' neighbourhood merge in
// spirit but over one axis instead of a 2-D window).
[[nodiscard]] std::vector<Candidate> suppress_non_maximal(std::vector<Candidate> candidates,
                                                          std::int64_t merge_window_bins);

// Stage 3: the final rank-sum p-value for each surviving candidate, computed
// against the pooled pixel values (not the per-anchor means) of the same
// local background window used in preselection.
void compute_pvalues(std::vector<Candidate>& candidates, const Band& horizontal,
                     const Band& vertical, const DetectOptions& options, unsigned int threads);

// Stage 4: Benjamini-Hochberg FDR across every candidate handed in (expected
// to span every chromosome and both orientations), keeping those at or below
// options.fdr_q. Candidates are left in their input order; the caller sorts
// for output if it wants a particular order.
[[nodiscard]] std::vector<Candidate> apply_fdr(std::vector<Candidate> candidates,
                                               double fdr_q);

}  // namespace hicx::stripes

#endif  // HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
