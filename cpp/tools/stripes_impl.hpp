// The computational core of hicDetectStripes, separated from the tool so
// that the unit tests (cpp/tests/test_detect_stripes.cpp) can reach it
// directly.
//
// hicDetectStripes is a PLAN.md tier 9 feature (section 9.3): HiCExplorer 3.7
// has no stripe caller, so there is no Python behaviour to reproduce. The
// project owner directed a faithful C++ reimplementation of Stripenn's own
// method (Yoon et al. 2022, package stripenn 1.1.65.22), not an invented
// stand-in, on the reasoning that a faithful port converges to Stripenn's
// own calls by construction and inherits its validated recall/precision on
// real data, where an ad hoc statistical test did not (see cpp/STATUS.md and
// the report to the orchestrating session for the numbers that motivated
// this). What follows is a line-by-line port of stripenn/getStripe.py's
// detection pipeline (StripeSearch, verticalLine, block, RemoveRedundant),
// with two classes of deliberate, documented simplification where an exact
// port would need weeks rather than hours:
//
//  1. **Edge detection engine.** Stripenn calls skimage.feature.canny with a
//     sigma and no explicit hysteresis thresholds, so skimage's own default
//     threshold selection applies, which is not simple to reproduce without
//     importing skimage. This port implements a standard Canny pipeline
//     (Gaussian smoothing at the given sigma, Sobel gradients, non-maximum
//     suppression, multi-hop hysteresis) and chooses its thresholds with the
//     classic percentile heuristic (strong = 80th percentile of the nonzero
//     gradient magnitude, weak = 0.4 times that), not skimage's own default.
//     A median-based heuristic (Stripenn's own ImageProcessing.auto_canny,
//     which the codebase carries but does not actually call) was tried
//     first and rejected: on a sharp, near-uniform step edge, 1.5 times the
//     median exceeds every pixel's own magnitude, so no pixel ever reaches
//     the strong threshold and nothing is detected (caught by this port's
//     own unit test). Everything downstream of the edge image
//     (verticalLine, block, the column-pairing that turns edges into stripe
//     boxes) is ported faithfully.
//  2. **Background/null model for the p-value.** Stripenn's nulldist()
//     builds an empirical null by a specific stratified, weighted random
//     sample across every chromosome (proportional to chromosome size and to
//     the number of non-empty columns), tracked per genomic distance up to
//     400 bins, separately for stripes anchored "up" and "down" from the
//     diagonal. This port keeps the same statistical idea -- an empirical
//     null of (center minus flanking background) differences, ranked
//     against, per distance -- built from a simpler uniform random sample
//     over valid anchors of each chromosome being processed (not
//     genome-wide, and not stratified by chromosome size), because
//     reproducing the exact sampling weights is not needed for the method
//     to work and would cost disproportionate implementation time. The
//     per-candidate p-value computation itself (StripeSearch.pvalue: the
//     empirical rank of the candidate's own center-vs-flank difference
//     against the null, taking the worse of the left and right side,
//     median over the stripe's rows) is ported faithfully.
//
// Benjamini-Hochberg FDR across every candidate (hicx::stats::
// benjamini_hochberg_adjusted) is applied on top, genome-wide: Stripenn
// itself only thresholds its raw p-value (--pvalue), it has no
// multiple-testing correction, and PLAN.md 9.3 requires one.
//
// Determinism (cpp/OPTIMIZATION.md 3). The frame scan is one independent
// computation per (chromosome, maxpixel percentile, frame index) triple,
// threaded with hicx::parallel_for and written into a preallocated slot;
// results are combined in index order. The background sample uses a
// deterministic seeded RNG (std::mt19937_64), not Python's random.Random, so
// it is not bit-identical to a Stripenn run with the "same" seed, only to
// itself at a fixed --threads.

#ifndef HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
#define HICX_TOOLS_DETECT_STRIPES_IMPL_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace hicx::stripes {

// A dense image, row major, values usually in [0, 1] once normalised.
struct Image {
    int rows = 0;
    int cols = 0;
    std::vector<double> v;

    Image() = default;
    Image(int r, int c, double fill = 0.0)
        : rows(r), cols(c), v(static_cast<std::size_t>(r) * static_cast<std::size_t>(c), fill) {}

    [[nodiscard]] double& at(int r, int c) { return v[static_cast<std::size_t>(r) * cols + c]; }
    [[nodiscard]] double at(int r, int c) const { return v[static_cast<std::size_t>(r) * cols + c]; }
};

// numpy.quantile(values, q, interpolation='linear') on the given values
// (NaN not expected; empty input returns 0).
[[nodiscard]] double quantile(std::vector<double> values, double q);

// cv.filter2D(img, -1, ones(k,k)/(k*k)): a box (mean) blur, "same" size,
// edge-replicated border (Stripenn/OpenCV's own default is BORDER_REFLECT_101;
// replication is used here as a documented simplification that matters only
// at the few-pixel image border).
[[nodiscard]] Image box_blur(const Image& img, int k);

// A standard Canny edge detector on a grayscale image already scaled to
// [0, 1]: Gaussian smoothing at `sigma`, Sobel gradients, non-maximum
// suppression, hysteresis thresholding with the median-based auto
// thresholds documented at the top of this file. Output is 0/1.
[[nodiscard]] Image canny(const Image& gray01, double sigma);

// stripenn.ImageProcessing.verticalLine: a Sobel-like gradient-direction
// filter over a binary edge image, keeping edge pixels whose gradient
// orientation (atan2(Filtered_X, Filtered_Y), degrees, wrapped to [0, 360))
// falls in (L, H); ported with its column shift (y -= 1).
[[nodiscard]] Image vertical_line(const Image& edges, double low_degrees = 60.0,
                                  double high_degrees = 120.0);

struct BlockResult {
    int length = 0;
    int end = 0;
};

// stripenn.ImageProcessing.block: the longest contiguous run of 1s in
// column c (ORed with its immediate neighbours), tolerating gaps of up to 4,
// and the row index the run ends at.
[[nodiscard]] BlockResult block_scan(const Image& vert, int column);

// One candidate stripe box, in the frame's local (row, column) coordinates:
// column range [x, x + w), row range [y, y + h). w is always the narrow
// (anchor) axis, bounded by --maxWidth; h is the long (extent) axis.
struct FrameCandidate {
    int x = 0;
    int y = 0;
    int w = 0;
    int h = 0;
};

// stripenn.getStripe.StripeSearch, minus the medpixel/Mean/total columns
// (not needed downstream) and minus the zero-row/column compaction Stripenn
// applies before calling it (documented simplification: candidates are
// found in the frame's own coordinates directly, which very sparse
// telomeric frames can dilute slightly; this is bounded by the same
// --minRawCount style implicit filtering the maxpixel/edge thresholds
// already provide). `submat` is the frame's dense observed-count block
// (symmetric, diagonal included), `maxpixel_value` is that percentile's
// count value M for the chromosome.
[[nodiscard]] std::vector<FrameCandidate> stripe_search_frame(const Image& submat,
                                                               double maxpixel_value,
                                                               double canny_sigma, int min_length,
                                                               int max_width, int blur_filter);

// A candidate mapped to genomic coordinates and chromosome, the unit the
// rest of the pipeline (redundancy removal, background sampling, FDR)
// works on.
struct Candidate {
    std::string chrom;
    std::int64_t pos1 = 0;  // narrow axis start (bp, or a bin index before genomic mapping)
    std::int64_t pos2 = 0;  // narrow axis end
    std::int64_t pos3 = 0;  // long axis start
    std::int64_t pos4 = 0;  // long axis end
    int frame_index = 0;    // which 200-bin frame this came from (for cross-frame dedup)
    double mean = 0.0;      // mean raw count over the box, stripenn.py 'Mean'
    // The same box's mean count shifted by the background window in the
    // narrow-axis direction, computed once at detection time from the same
    // dense frame the box was found in -- an apples-to-apples comparison,
    // the box's own shape against itself shifted, rather than a single
    // (anchor, distance) lookup or a small fixed window (see
    // hicDetectStripes.cpp candidate_pvalue call site).
    double left_mean = 0.0;
    double right_mean = 0.0;
    double pvalue = 1.0;
    double qvalue = 1.0;
};

// stripenn.getStripe.RemoveRedundant(by='size'): among candidates that
// overlap by more than 20% on both axes, keep the more elongated one
// (larger h/w). `by_pvalue` selects the 'pvalue' variant (keep the smaller
// p-value) used for the final cross-percentile merge.
[[nodiscard]] std::vector<Candidate> remove_redundant(std::vector<Candidate> candidates,
                                                       bool by_pvalue, unsigned int threads);

// A background sample built from many random (anchor, distance) draws on
// the chromosome's own obs/exp-free raw band: at distance d (0 up to the
// band's max distance, capped at 399 as Stripenn's own tables are), the
// pooled empirical sample of "this pixel's count minus a same-distance
// pixel a few anchors to the left/right" (see the file header, point 2, for
// why this replaces Stripenn's 2-D window-averaged, direction-split
// nulldist() table: it reuses this tool's own already-verified band
// construction and keeps the same statistical idea -- an empirical local
// background contrast, ranked -- at a fraction of the implementation cost).
struct BackgroundModel {
    // background_model.left[d] / .right[d] is the sample at distance d.
    std::vector<std::vector<double>> left;
    std::vector<std::vector<double>> right;
};

// The candidate's p-value: the worse (larger) of the empirical rank of its
// own mean center-minus-left and center-minus-right differences (averaged
// over the candidate's rows, at its own distance from the diagonal) against
// the background model at that distance -- the fraction of the background
// sample at least as large as the observed difference, matching
// stripenn.getStripe.pvalue's per-row p1/p2/max(p1,p2) and its "never
// exactly zero" floor of 1/len(sample).
[[nodiscard]] double candidate_pvalue(const BackgroundModel& model, int distance,
                                      double left_diff, double right_diff);

}  // namespace hicx::stripes

#endif  // HICX_TOOLS_DETECT_STRIPES_IMPL_HPP
