// The computational core of hicArrowhead, separated from the tool so the
// unit tests (cpp/tests/test_arrowhead.cpp) can reach it directly.
//
// hicArrowhead is a PLAN.md tier 9 feature: HiCExplorer has no Arrowhead
// contact-domain caller (hicFindTADs uses an independent insulation-score
// method), so this is a faithful reimplementation of Juicer/Juicebox's own
// Arrowhead (Rao et al. 2014, Cell, doi: 10.1016/j.cell.2014.11.021), ported
// line by line from the real Java source
// (github.com/aidenlab/Juicebox src/juicebox/tools/utils/juicer/arrowhead/
// MatrixTriangles.java, DynamicProgrammingUtils.java, BlockResults.java,
// BlockBuster.java, CumulativeBlockResults.java, BinnedScore.java,
// HighScore.java, and connectedcomponents/BinaryConnectedComponents.java),
// not from the paper's method section alone.
//
// Deliberate deviations from the real tool:
//
//  1. **Input matrix / normalization.** Juicer's Arrowhead reads a .hic file
//     directly and selects a normalization vector (NONE/VC/VC_SQRT/KR) at
//     run time. hicArrowhead instead works the way every other HiCExplorer
//     v4 TAD/loop tool does: on whatever matrix the input file already
//     represents (a cool file's own weight column, an h5's own correction
//     factors, or nothing). There is no runtime normalization selector.
//  2. **Connected components.** BinaryConnectedComponents.java is a
//     two-pass union-find labelling scheme for 8-connectivity components of
//     the thresholded block-score matrix. Any correct connected-components
//     algorithm partitions the same pixels into the same components, so
//     this is a plain BFS flood fill (connected_components below) instead:
//     mathematically equivalent, not a faithful line-by-line port of that
//     specific algorithm's mechanics.
//  3. **Feature/control list scoring** (ArrowheadScoreList's non-empty-list
//     path, Arrowhead's optional 4th/5th CLI arguments): not ported. It
//     only matters when the user supplies an external domain list to score
//     against this run's block-score matrix, which is a secondary
//     diagnostic feature, not the domain-calling pipeline itself.
//  4. **Output format.** Juicer's own BEDPE carries viewer-oriented
//     "name"/"strand1"/"strand2"/"color" columns with placeholder values.
//     hicArrowhead's output keeps only the fields that carry real
//     information: chrom, x1, x2 (== y1, y2: Arrowhead calls square blocks
//     on the diagonal), score, uVarScore, lVarScore, upSign, loSign.
//  5. **Sliding window sizing at the chromosome's last window.**
//     BlockBuster.callSubBlockbuster computes a window of n = limEnd -
//     adjustedLimStart + 1 bins (one more than --matrixWidth) but then
//     extracts the observed matrix at (limStart, limEnd) -- the
//     *unshifted* start, not adjustedLimStart, even though adjustedLimStart
//     is what the window was supposedly resized around. For every window
//     except the last one on a chromosome the two coincide, so this is
//     invisible; for the last, partial window it means the extracted data
//     starts one bin off from where the window nominally starts and the
//     matrix's trailing rows/columns are left as zero padding. This reads
//     as an upstream edge-case inconsistency rather than an intentional
//     design, so it is not reproduced: hicArrowhead uses a single,
//     consistent half-open window (size exactly --matrixWidth bins,
//     shifted left and clamped so the last window on a chromosome still
//     covers --matrixWidth bins) for both the observed-matrix extraction
//     and the coordinate offset.
//
// Everything else -- the directionality-index-upstream transform, the three
// dynamic-programming sum matrices (right/upper), the up/lo triangle means
// and variances, the block-score combination and thresholding, the
// sliding-window scan with its low-confidence/high-confidence two-pass
// scheme and set-difference merge, the distance-based re-binning, and the
// final descending sort by (uVarScore + lVarScore) -- is a faithful,
// line-by-line port.

#ifndef HICX_TOOLS_ARROWHEAD_IMPL_HPP
#define HICX_TOOLS_ARROWHEAD_IMPL_HPP

#include <cstdint>
#include <functional>
#include <vector>

namespace hicx::arrowhead {

// A dense, row-major, square matrix. Arrowhead's own windows are bounded by
// --matrixWidth (Juicer's own default 2000), so a dense representation,
// matching the Java RealMatrix the port is taken from, is both simpler and
// bounded in memory.
struct Matrix {
    int n = 0;
    std::vector<double> v;

    Matrix() = default;
    explicit Matrix(int size) : n(size), v(static_cast<std::size_t>(size) * static_cast<std::size_t>(size), 0.0) {}

    [[nodiscard]] double& at(int r, int c) { return v[static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + c]; }
    [[nodiscard]] double at(int r, int c) const {
        return v[static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + c];
    }
};

// BlockResults.calculateDirectionalityIndexUpstream: for each row i, the
// preference (A - B) / (A + B) between the `window` counts immediately
// behind position i-gap and the `window` counts immediately ahead of i+gap
// (window is the smaller of the room on each side, capped at n), written
// into row i at columns i+gap .. i+window.
[[nodiscard]] Matrix directionality_index_upstream(const Matrix& observed, int gap = 7);

// One connected component's highest-scoring block.
struct HighScore {
    int i = 0;
    int j = 0;
    double score = 0.0;
    double u_var_score = 0.0;
    double l_var_score = 0.0;
    double up_sign = 0.0;
    double lo_sign = 0.0;

    [[nodiscard]] int width() const { return j > i ? j - i : i - j; }
    // (uVarScore + lVarScore): BlockBuster's own final sort key
    // (HighScore.sortValue/compareTo), descending.
    [[nodiscard]] double sort_value() const { return u_var_score + l_var_score; }
};

// MatrixTriangles + BlockResults.<init>: the directionality-index matrix's
// block-score matrix (generateBlockScoreCalculations), thresholded
// (thresholdScoreValues) and reduced to one HighScore per connected
// component of positive entries (extractConnectedComponents +
// calculateResults), for one window. `var_threshold` may be NaN (the low-
// confidence pass skips the variance-sum threshold entirely, matching
// thresholdScoreValues(Double.NaN, signThreshold)).
[[nodiscard]] std::vector<HighScore> block_results(const Matrix& observed, double var_threshold,
                                                    double sign_threshold);

// BlockBuster's tunable thresholds and sliding-window size, its own static
// fields and their defaults.
struct Options {
    double var_threshold = 0.2;
    double high_sign_threshold = 0.5;
    double max_low_sign_threshold = 0.4;
    double min_low_sign_threshold = 0.0;
    double decrement_low_sign_threshold = 0.1;
    int min_block_size = 60;
    int matrix_width = 2000;
};

// BlockBuster.run for one chromosome's `chrom_length` bins (the
// `observed_window` callback supplies the symmetric, diagonal-included,
// already-normalized observed-count window [lim_start, lim_end) the real
// tool reads from the .hic file; here it comes from whatever matrix the
// caller loaded). Indices in the result are in bins, already offset to the
// chromosome's own coordinate system (HighScore.offsetIndex) but not yet
// scaled by resolution (HighScore.scaleIndicesByResolution is the caller's
// job, once bp-per-bin is known) -- scale and sort by sort_value()
// descending are both already applied, matching BlockBuster.run's own
// binScoresByDistance + Collections.reverseOrder() + sort.
[[nodiscard]] std::vector<HighScore> call_domains(
    int chrom_length_bins, const Options& options,
    const std::function<Matrix(int lim_start, int lim_end)>& observed_window);

}  // namespace hicx::arrowhead

#endif  // HICX_TOOLS_ARROWHEAD_IMPL_HPP
