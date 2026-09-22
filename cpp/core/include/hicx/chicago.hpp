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
#include <string>
#include <vector>

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
[[nodiscard]] std::vector<ChinputRecord> read_chinput(const std::string& path);

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

}  // namespace hicx::chicago

#endif  // HICX_CHICAGO_HPP
