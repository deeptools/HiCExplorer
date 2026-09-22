// CHiCAGO scoring (Cairns et al. 2016, Genome Biology 17:127), added to the
// chic* tools as an alternative to HiCExplorer's own negative-binomial score
// (PLAN.md 9.15).
//
// Scope of this file, as actually implemented and validated so far: the
// closed-form parts of the CHiCAGO pipeline that can be computed from
// already-fitted parameters, independent of the genome-wide binning that
// produces those parameters in the first place:
//
//   * the Brownian/technical-noise p-value (getPvals in the R package,
//     R/Chicago/R/pval.R): a Delaporte (Poisson + Gamma-Poisson) survival
//     function, with R's own fallback to a plain negative-binomial
//     approximation for p-values too small for the direct convolution;
//   * the distance-based p-value weighting and score (getScores,
//     R/Chicago/R/scores.R): a logistic weight curve and its baseline value;
//   * the cubic log-log distance function fit (estimateDistFun,
//     R/Chicago/R/distFun.R): ordinary least squares of log(refBinMean) on
//     log(midpoint) up to a cubic term, with linear extrapolation outside the
//     observed range.
//
// NOT implemented here (see the final report handed to the orchestrating
// session for the exact accounting): estimateTechnicalNoise's genome-wide
// trans-count binning, normaliseBaits/normaliseOtherEnds, and
// estimateBrownianComponent's dispersion sampling. Those require the .npb/
// .nbpb/.poe genome-wide design tables and an iterative estimation procedure
// that was not reproduced and verified against R in the time available. The
// functions below take Bmean, Tmean and the dispersion (alpha) as already
// estimated inputs, exactly as R's getPvals and getScores do (they too take
// cd@params$dispersion and the per-row Bmean/Tmean columns as given).
//
// Validated against real R Chicago 1.38.0 / Delaporte 8.4.3 output on
// PCHiCdata 1.38.0's GM12878 (chr20/21) and mouse ES (chr18/19) chinput
// files: see cpp/tests/test_chicago.cpp.

#ifndef HICX_CHICAGO_HPP
#define HICX_CHICAGO_HPP

#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hicx/matrix_data.hpp"

namespace hicx::chicago {

// getPvals: log of the upper-tail survival probability P(N' >= N), N' being
// the Brownian-plus-technical-noise count.
//
//   * Bmean < DBL_EPSILON: pure technical noise, log ppois(N - 1, Tmean,
//     lower.tail = FALSE).
//   * otherwise: log pdelap(N - 1, alpha, beta = Bmean / alpha, lambda =
//     Tmean, lower.tail = FALSE), alpha being the fitted dispersion; when
//     that convolution underflows to 0 or is not finite, R substitutes
//     pnbinom(N - 1, size = min(alpha * (1 + Tmean / Bmean)^2, 1e10),
//     mu = Bmean + Tmean, lower.tail = FALSE), reproduced identically here.
[[nodiscard]] double log_pvalue(double N, double alpha, double Bmean, double Tmean);

// R's weightAlpha/weightBeta/weightGamma/weightDelta defaults
// (34.1157346557331, -2.58688050486759, -17.1347845819659, -7.07609217973722
// as read from a real chicagoData object's settings; not literal defaults in
// this header, the caller reads them from wherever it reads --scoring
// parameters, this struct just carries them to log_weight/log_weight_min).
struct WeightSettings {
    double alpha = 34.1157346557331;
    double beta = -2.58688050486759;
    double gamma = -17.1347845819659;
    double delta = -7.07609217973722;
};

// .getWeights for a single (finite, non-negative) distance, given eta.bar
// (the genome-wide average weight, computed once per experiment by
// eta_bar_from_design below).
[[nodiscard]] double log_weight(double abs_dist, const WeightSettings& w, double eta_bar);

// score = max(-log_weight(0, w, eta_bar) - (log_p - log_weight(dist, w,
// eta_bar)), 0), i.e. getScores' method = "weightedRelative" (the CHiCAGO
// default and the only mode PLAN.md 9.15 asks for).
[[nodiscard]] double score_from_pvalue(double log_p, double abs_dist, const WeightSettings& w,
                                        double eta_bar);

// A restriction fragment, as one line of a .rmap file: chromosome, start, end
// (1-based, inclusive, as CHiCAGO's own files are), fragment ID.
struct RmapFragment {
    std::string chrom;
    long start = 0;
    long end = 0;
    long id = 0;
};

// A baited fragment, as one line of a .baitmap file: the same four fields
// plus the bait's name.
struct BaitmapFragment {
    std::string chrom;
    long start = 0;
    long end = 0;
    long id = 0;
    std::string name;
};

// One interaction record from a .chinput file (baitID, otherEndID, N,
// otherEndLen, distSign), as produced by chicagoTools/bam2chicago or an
// equivalent upstream pipeline.
struct ChinputRecord {
    long bait_id = 0;
    long other_end_id = 0;
    double N = 0.0;
    long other_end_len = 0;
    // distSign is NA (represented here as has_dist_sign = false) for trans
    // interactions.
    bool has_dist_sign = false;
    long dist_sign = 0;
};

[[nodiscard]] std::vector<RmapFragment> read_rmap(const std::string& path);
[[nodiscard]] std::vector<BaitmapFragment> read_baitmap(const std::string& path);
[[nodiscard]] std::vector<ChinputRecord> read_chinput(const std::string& path, int threads = 1);

// The reverse direction of chinput_from_matrices: builds a real Hi-C matrix
// (one bin per .rmap restriction fragment, in the .rmap's own order) with
// one entry per .chinput row, so a .chinput file (however it was produced -
// bam2chicago, a prior chinput_from_matrices run, or any other CHiCAGO
// pipeline) can be written out as a real .cool/.h5/.hic file through
// hicConvertFormat --inputFormat chinput, instead of staying in a format
// only R's Chicago package reads.
//
// A .rmap fragment's 1-based, inclusive [start, end] becomes the bin table's
// usual 0-based, half-open [start - 1, end), the same convention
// chinput_from_matrices' own header comment documents (and the inverse of
// it: round-tripping a matrix through chinput_from_matrices and back through
// this function reproduces the same N values at the same fragment pairs).
// Every row is one observed direction of what a Hi-C matrix represents
// symmetrically, so entries are canonicalised to (min bin, max bin) and the
// result is marked Symmetry::UpperTriangle: a query at either (bait, other
// end) or (other end, bait) then finds the same value, matching how every
// other Hi-C matrix this project reads or writes is stored.
//
// A .chinput row whose baitID or otherEndID is not one of the .rmap's own
// fragment ids fails loudly (a real design/data mismatch, not something to
// silently drop rows over).
[[nodiscard]] hicx::MatrixData matrix_from_chinput(const std::string& chinput_path,
                                                     const std::string& rmap_path);

// eta.bar (.getEtaBar in R): the genome-wide average of exp(log_weight) over
// every possible bait-fragment pair up to the chromosome ends, weighted the
// way R enumerates them (one arithmetic sequence of distances per bait,
// spaced by the average fragment length). avg_frag_len is .getAvgFragLength,
// the mean fragment length over the whole .rmap excluding chrMT/chrM.
[[nodiscard]] double avg_frag_length(const std::vector<RmapFragment>& rmap);

// .getNoOfHypotheses(cd, includeTrans = TRUE): sum(nBaits) * (2 * nrow(rmap)
// - sum(nBaits) - 1) / 2, the number of possible bait/other-end pairs.
[[nodiscard]] std::size_t n_hypotheses(const std::vector<RmapFragment>& rmap,
                                        const std::vector<BaitmapFragment>& baitmap);

[[nodiscard]] double eta_bar_from_design(const std::vector<RmapFragment>& rmap,
                                          const std::vector<BaitmapFragment>& baitmap,
                                          const WeightSettings& w, double avg_frag_len,
                                          bool include_trans, std::size_t n_hypotheses);

// estimateDistFun's cubic fit: log(refBinMean) ~ log(midpoint) +
// log(midpoint)^2 + log(midpoint)^3 over the bins with a defined refBinMean,
// with linear extrapolation of the fitted curve below the smallest and above
// the largest observed midpoint (the head/tail coefficients R also stores).
struct DistFunFit {
    // cubicFit: intercept, linear, quadratic, cubic coefficients of the OLS
    // fit, in that order (matches R's `fit` vector).
    double cubic[4] = {0, 0, 0, 0};
    double obs_min_log = 0;  // log(min(midpoint))
    double obs_max_log = 0;  // log(max(midpoint))
    double head_coef[2] = {0, 0};  // alpha, beta of the linear extrapolation below obs_min
    double tail_coef[2] = {0, 0};  // alpha, beta of the linear extrapolation above obs_max
};

// midpoints and refBinMean must be the same length, one entry per distance
// bin, refBinMean already the observed mean strictly positive (rows with
// NA/undefined refBinMean excluded before calling, as R's is.na(refBinMean)
// == FALSE filter does).
[[nodiscard]] DistFunFit fit_distance_function(const std::vector<double>& midpoints,
                                                const std::vector<double>& ref_bin_mean);

// Evaluates the fitted distance function at a distance (predicts
// log(refBinMean); the caller exponentiates for a mean).
[[nodiscard]] double eval_distance_function(const DistFunFit& fit, double distance);

// ---------------------------------------------------------------------------
// Genome-wide parameter estimation: readSample's per-file ingestion filter,
// estimateTechnicalNoise's trans-count binning, normaliseBaits/
// normaliseOtherEnds' non-shrunken (default) scaling factors, and
// estimateBrownianComponent's dispersion.
//
// Scope note: R's readAndMerge merges several replicate .chinput files
// (mergeSamples) before readSample's own filtering. That multi-file merge is
// not reproduced here; callers pass one already-merged .chinput (summed
// per bait/otherEnd pair across replicates by the caller, or a single
// replicate file, exactly what read_chinput already reads). readSample's own
// per-file filtering (fragment length, self-ligation, minNPerBait,
// removeAdjacent, bait2bait-only baits) is reproduced in full below.
//
// Also not reproduced: normaliseBaits/normaliseOtherEnds' shrink=TRUE path
// (gamma MLE + loess smoothing). Chicago's own default is shrink=FALSE for
// normaliseBaits and normaliseOtherEnds always calls it with shrink=FALSE,
// so this covers Chicago's actual default pipeline, not every option it
// exposes.
// ---------------------------------------------------------------------------

// Settings readSample/estimateTechnicalNoise/estimateBrownianComponent need,
// all Chicago default values (Chicago::defaultSettings()).
struct FilterSettings {
    long min_frag_len = 150;
    long max_frag_len = 40000;
    long min_n_per_bait = 250;
    bool remove_adjacent = true;
    long max_l_brown_est = 1500000;
    long binsize = 20000;
    bool adj_bait2bait = true;
};

// One row of x after readSample's filtering: a bait/other-end interaction
// with its trans flag (has_dist_sign == false) preserved, and isBait2bait
// precomputed (wb2b: otherEndID appears in the baitmap).
struct ChiInteraction {
    long bait_id = 0;
    long other_end_id = 0;
    double N = 0.0;
    long other_end_len = 0;
    bool has_dist_sign = false;
    long dist_sign = 0;
    bool is_bait2bait = false;
};

// readSample (R/Chicago/R/readData.R), the ingestion path for ONE chinput
// file already read by read_chinput: otherEndLen range filter, self-ligation
// removal, minNPerBait filter, optional removeAdjacent, and dropping baits
// whose only proximal (within max_l_brown_est) interactions are all
// bait2bait.
//
// Takes raw by value (move it in with std::move at the call site) for the
// same reason fit_chicago_background does: raw is only read in this
// function's first stage, and freed right after, before the three remaining
// filter stages run without it. A const& parameter would keep it alive for
// the whole call (the caller's storage does not go away just because this
// function stops reading it), which is most of read_sample's own cost.
[[nodiscard]] std::vector<ChiInteraction> read_sample(std::vector<ChinputRecord> raw,
                                                        const std::vector<BaitmapFragment>& baitmap,
                                                        const FilterSettings& fs, int threads = 1);

// Derives ChinputRecord rows directly from a Hi-C contact matrix (cool, h5 or
// .hic, through hicx::ToolMatrix::load's genuinely partial, chromosome-scoped
// path: one load per chromosome the .baitmap actually uses, never a
// whole-genome load) instead of a .chinput file: for every bait fragment and
// every .rmap fragment on the same chromosome within fs.max_l_brown_est of
// it, the matrix's contact value between the two fragments' genomic spans
// becomes N. Several matrix paths are summed per (baitID, otherEndID) pair,
// the matrix-input equivalent of pre-summing several replicate .chinput
// files before this pipeline sees them (see the note on chicChicagoBackground
// Model's own --chinput cardinality).
//
// Coordinate convention: a .rmap fragment is CHiCAGO's own 1-based,
// inclusive-of-both-ends span (contiguous fragments: the next fragment's
// start is the previous fragment's end + 1). A Hi-C matrix's bin table is
// 0-based, half-open, the convention every other coordinate this project
// feeds to BinTable::region_bin_range already uses (reference point BED
// files, --region strings). A fragment [start, end] is therefore queried as
// the half-open region [start - 1, end), exactly the span
// hicBuildMatrix --restrictionCutFile produces from a BED digest whose
// 1-based .rmap is that same BED shifted by one. A matrix at a fixed bin
// resolution (not restriction-fragment resolution) is handled the same way:
// every bin overlapping [start - 1, end) contributes, and their values are
// summed.
//
// Scope: only cis pairs within fs.max_l_brown_est are enumerated, matching
// the distance bound the tool's own Brownian estimation already uses
// (FilterSettings.max_l_brown_est; see read_sample's own proximal-only use
// of it). Trans pairs and cis pairs beyond that bound are not derived from
// the matrix: reproducing them would mean walking the whole genome for every
// bait, which this tool's proximal design does not need. A run that also
// needs trans-based technical-noise estimation from real reads still needs a
// .chinput file for that (--chinput remains the only way to feed trans
// counts in).
[[nodiscard]] std::vector<ChinputRecord> chinput_from_matrices(
    const std::vector<std::string>& matrix_paths, const std::vector<RmapFragment>& rmap,
    const std::vector<BaitmapFragment>& baitmap, const FilterSettings& fs);

// One (bait, other-end) fragment pool assignment, from Hmisc::cut2's default
// quantile-binning algorithm (cuts missing, onlycuts = FALSE): group index
// (1-based) per input element, replicating cut2's own y-vector semantics
// (the native grouping technicalNoise's tblb binning uses directly).
[[nodiscard]] std::vector<int> cut2_native_groups(const std::vector<double>& x, long m);

// cut2(x, m = m, onlycuts = TRUE): the g roughly-equal-count breakpoints
// (unique(c(low, max(x)))), used as input to a plain R cut() elsewhere
// (addTLB calls cut2 this way, then cut() on a different, larger vector).
[[nodiscard]] std::vector<double> cut2_cuts(const std::vector<double>& x, long m);

// R's cut(x, breaks, right = TRUE, include.lowest = TRUE): 1-based bin index,
// (breaks[i-1], breaks[i]] except the first bin which also includes
// breaks[0] itself; 0 if x falls outside [breaks.front(), breaks.back()].
[[nodiscard]] int cut_with_breaks(double x, const std::vector<double>& breaks);

// .addTLB: bins other ends by their (bait2bait-adjusted) trans-interaction
// count into pools, returned as a 1-based pool id per otherEndID (the same
// tlb column normaliseOtherEnds and estimateTechnicalNoise both consume).
// Non-bait2bait and bait2bait other ends are binned separately and never
// share a pool id (bait2bait ids are offset past the last non-B2B id, the
// same separation R's "B2B" suffix gives its factor levels).
struct TlbResult {
    std::unordered_map<long, int> pool_of_other_end;  // otherEndID -> tlb pool id
    int n_non_b2b_pools = 0;
    int n_b2b_pools = 0;
};
[[nodiscard]] TlbResult add_tlb(const std::vector<ChiInteraction>& x, const FilterSettings& fs,
                                 double tlb_filter_top_percent, long tlb_min_prox_oe_per_bin,
                                 long tlb_min_prox_b2b_per_bin, int threads = 1);

// estimateTechnicalNoise: bins baits by their observed trans-interaction
// count (tblb, separate from the other-end tlb pools above) and computes the
// Poisson mean trans-count Tmean for every observed (tlb, tblb) pool,
// dividing the observed trans counts in that pool by the total number of
// possible bait/other-end pairs the pool could have produced.
struct TechnicalNoiseResult {
    std::unordered_map<long, int> tblb_of_bait;             // baitID -> tblb pool id
    std::map<std::pair<int, int>, double> tmean_by_pool;     // (tlb, tblb) -> Tmean
};
[[nodiscard]] TechnicalNoiseResult estimate_technical_noise(
    const std::vector<ChiInteraction>& x, const TlbResult& tlb,
    const std::vector<RmapFragment>& rmap, const std::vector<BaitmapFragment>& baitmap,
    long min_baits_per_bin, int threads = 1);

// normaliseFragmentSets' non-shrunken path (shrink = FALSE), specialised to
// the bait side: for every baitID, s_j = median over distance bins of
// (binwise sum of N / total possible other ends in that bin) / (geometric
// mean of that ratio over all baits in the bin); also returns refBinMean,
// the geometric mean per distance bin (estimateDistFun's own input).
struct BaitFactors {
    std::unordered_map<long, double> s_j;              // baitID -> s_j
    std::vector<double> ref_bin_mean_by_distbin;        // 1-based distbin -> refBinMean (NaN if undefined)
};
[[nodiscard]] BaitFactors normalise_baits(const std::vector<ChiInteraction>& x,
                                           const std::vector<RmapFragment>& rmap,
                                           const std::vector<BaitmapFragment>& baitmap,
                                           const std::string& npb_path, const FilterSettings& fs,
                                           int threads = 1);

// normaliseFragmentSets' non-shrunken path, other-end side: for every tlb
// pool, s_i = median over distance bins of (binwise sum of N / total
// possible baits in that bin) / (geometric mean of that ratio over all tlb
// pools in the bin).
// bait_s_j is normalise_baits' own s_j output: normaliseOtherEnds sums NNb =
// pmax(1, round(N / s_j)), the bait-normalised count, not the raw N (R:
// normaliseOtherEnds(cd, Ncol = "NNb", ...), NNb set by normaliseBaits).
[[nodiscard]] std::unordered_map<int, double> normalise_other_ends(
    const std::vector<ChiInteraction>& x, const TlbResult& tlb,
    const std::unordered_map<long, double>& bait_s_j, const std::string& nbpb_path,
    const FilterSettings& fs, int threads = 1);

// A proximal (baitID, otherEndID) pair from the precomputed .poe design
// file: every pair within max_l_brown_est of each other, after removeb2b
// and removeAdjacent, with its exact distance.
struct ProxOePair {
    long bait_id = 0;
    long other_end_id = 0;
    double dist = 0.0;
};
[[nodiscard]] std::vector<ProxOePair> read_poe(const std::string& path);

// estimateBrownianComponent's dispersion (MASS::glm.nb's theta, fit by
// MASS::theta.ml's Newton iteration on the negative-binomial profile
// likelihood, exactly as R computes it): N ~ NB(mean = Bmean, dispersion =
// alpha) with Bmean fixed by the offset (s_j * s_i * distance function), no
// free regression coefficients, so glm.nb's IRLS loop never changes Bmean
// and this reduces to a single theta.ml(N, Bmean) call.
[[nodiscard]] double estimate_dispersion_theta_ml(const std::vector<double>& N,
                                                    const std::vector<double>& Bmean, int threads = 1);

// s_j * s_i * eval_distance_function(fit, |dist|), i.e. estimateBMean; NaN
// distSign (trans) maps to Bmean = 0, matching R's x[is.na(distSign), Bmean := 0].
[[nodiscard]] double estimate_bmean(double s_j, double s_i, double abs_dist_or_nan,
                                     const DistFunFit& fit);

// The full parameter set chicChicagoBackgroundModel estimates and
// chicChicagoScores consumes: everything above, wired into one pipeline in
// R's own order (normaliseBaits -> normaliseOtherEnds -> estimateTechnicalNoise
// -> estimateDistFun -> estimateBrownianComponent).
struct BackgroundModel {
    FilterSettings filters;
    TlbResult tlb;
    BaitFactors bait_factors;
    std::unordered_map<int, double> s_i_by_tlb_pool;
    TechnicalNoiseResult tech_noise;
    DistFunFit dist_fun;
    double dispersion = 0.0;
    // Diagnostics: how many (bait, other-end) pairs from the .poe design fed
    // the dispersion fit, and whether R's own brownianNoise.subset (1000
    // baits by default) would have triggered sub-sampling on this input
    // (report only; this implementation always fits on the full design,
    // see the header comment on estimate_dispersion_theta_ml's caller).
    std::size_t dispersion_n_pairs = 0;
    bool subset_would_trigger_in_r = false;
};

// Takes raw by value (move it in with std::move at the call site to avoid a
// copy): every interaction the .chinput file has, cis and trans alike, is
// only needed to build x (read_sample's filtered subset), and freed
// immediately after, before the rest of this function's work (addTLB,
// normaliseBaits/normaliseOtherEnds, estimateTechnicalNoise, the dispersion
// fit) runs without it. A const& parameter could not do that: the caller's
// storage would stay alive, unused, for the whole call.
[[nodiscard]] BackgroundModel fit_chicago_background(
    std::vector<ChinputRecord> raw, const std::vector<RmapFragment>& rmap,
    const std::vector<BaitmapFragment>& baitmap, const std::string& npb_path,
    const std::string& nbpb_path, const std::string& poe_path, const FilterSettings& fs,
    long tlb_min_baits_per_bin = 1000, double tlb_filter_top_percent = 0.01,
    long tlb_min_prox_oe_per_bin = 50000, long tlb_min_prox_b2b_per_bin = 2500,
    long brownian_noise_subset = 1000, int threads = 1);

}  // namespace hicx::chicago

#endif  // HICX_CHICAGO_HPP
