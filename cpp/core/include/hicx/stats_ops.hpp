// Statistics shared by the TAD, loop and cHi-C tools.
//
// This is a general core component, not a part of any one tool. Everything in
// it is reusable and is written against the Python function it replaces:
//
//   scipy.special.erf / erfc / ndtr        the cephes routines, vendored so
//                                          that a p-value printed with twelve
//                                          decimals matches the reference
//   scipy.stats.rankdata(method='average') tie corrected ranks
//   scipy.stats.ranksums                   Wilcoxon rank sum, two sided
//   the Benjamini-Hochberg cutoff and the Bonferroni scaling as
//   hicFindTADs.min_pvalue applies them
//   fit_nbinom.fit                         negative binomial MLE by L-BFGS-B
//   L-BFGS-B itself                        scipy.optimize.fmin_l_bfgs_b
//
// Who uses what, so that the next tool does not reimplement any of it:
//
//   hicFindTADs           ranksums, benjamini_hochberg_cutoff,
//                         bonferroni_in_place
//   hicDetectLoops        ranksums (its Wilcoxon preselection), fit_nbinom and
//                         nbinom_sf for the p-value of a candidate loop
//   chicViewpointBackgroundModel, chicViewpoint, chicSignificantInteractions
//                         fit_nbinom over the per distance distributions
//   hicDifferentialTAD    ranksums
//
// Determinism: every reduction here is sequential and in a fixed order, and
// nothing in this file is threaded. A caller that parallelises does so over
// independent problems (one p-value per boundary, one fit per distance) and
// combines the results in index order, which is what cpp/OPTIMIZATION.md 3
// requires. Where a sum must match numpy's, hicx::npy::pairwise_sum is used
// rather than a plain loop.

#ifndef HICX_STATS_OPS_HPP
#define HICX_STATS_OPS_HPP

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <vector>

namespace hicx::stats {

// --------------------------------------------------------------------------
// Cephes special functions, translated from
// scipy/special/special/cephes/ndtr.h (Cephes Math Library, Moshier;
// redistributed by SciPy under the BSD licence).

[[nodiscard]] double erf(double x);
[[nodiscard]] double erfc(double x);
// The normal cumulative distribution function, scipy.special.ndtr.
[[nodiscard]] double ndtr(double x);
// scipy.stats.norm.sf, which is ndtr(-x). This is the exact route scipy takes:
// _SimpleNormal.sf calls special.ndtr(-x).
[[nodiscard]] double normal_sf(double x);

// log(Gamma(x)) for x > 0, scipy.special.gammaln.
[[nodiscard]] double gammaln(double x);
// The digamma function, scipy.special.psi.
[[nodiscard]] double digamma(double x);

// --------------------------------------------------------------------------
// Ranks and the rank sum test

// scipy.stats.rankdata(values, method='average'): ties share the average of
// the ranks they span. Ranks are one based.
[[nodiscard]] std::vector<double> rankdata_average(std::span<const double> values);

struct RanksumsResult {
    double statistic = 0.0;
    double pvalue = 0.0;
    // scipy raises when both samples are empty; hicFindTADs guards against
    // that itself, but a caller that does not can test this.
    bool defined = true;
};

// scipy.stats.ranksums(x, y) with the default two sided alternative.
//
//   ranked   = rankdata(concatenate(x, y))
//   s        = sum(ranked[:len(x)])            numpy pairwise order
//   expected = n1 * (n1 + n2 + 1) / 2
//   z        = (s - expected) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
//   p        = 2 * norm.sf(|z|)
[[nodiscard]] RanksumsResult ranksums(std::span<const double> x,
                                      std::span<const double> y);

// --------------------------------------------------------------------------
// Multiple testing

// The Benjamini-Hochberg step up cutoff exactly as hicFindTADs.py:1232-1241
// computes it: the largest p that satisfies p <= q * rank / n over the sorted
// p-values, or 0.0 when none does. NaN must already have been replaced by 1,
// which is what the caller does, so this function does not do it again.
//
// Note that this is a *cutoff*, not a set of adjusted p-values. The caller
// then keeps every boundary whose raw p-value is at most the cutoff.
[[nodiscard]] double benjamini_hochberg_cutoff(std::vector<double> pvalues, double q);

// pvalues *= n, then clamp values above 1 to exactly 1. NaN is left alone,
// because `e > 1 if ~np.isnan(e) else False` selects nothing for NaN.
void bonferroni_in_place(std::vector<double>& pvalues);

// --------------------------------------------------------------------------
// L-BFGS-B and the negative binomial fit

// The outcome of an L-BFGS-B run. `status` mirrors the task string of
// scipy.optimize.fmin_l_bfgs_b: 0 converged, 1 iteration or evaluation limit
// reached, 2 abnormal termination in the line search.
struct LbfgsbResult {
    std::vector<double> x;
    double f = 0.0;
    int iterations = 0;
    int function_evaluations = 0;
    int status = 0;
};

struct LbfgsbOptions {
    // scipy's defaults for fmin_l_bfgs_b.
    int memory = 10;              // m
    double factr = 1e7;           // stop when (f_k - f_{k+1}) / max(...) <= factr * eps
    double pgtol = 1e-5;          // stop when max |projected gradient| <= pgtol
    int max_iterations = 15000;
    int max_function_evaluations = 15000;
    // approx_grad=1 uses a forward difference gradient with this step.
    double epsilon = 1e-8;
};

// A box constrained variable. An absent bound is +/- infinity.
struct Bound {
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();
};

// Minimise `objective` subject to the bounds, starting from x0, with a forward
// difference gradient of step `epsilon`, which is what approx_grad=1 does.
//
// This is a **projected limited-memory BFGS**, not a translation of the
// Fortran L-BFGS-B that scipy calls, and the difference is stated here rather
// than glossed over. It keeps the same limited-memory history size, the same
// two convergence tests (projected gradient against pgtol, relative function
// decrease against factr * eps) and the same forward-difference gradient, and
// it builds its search direction with the two loop recursion restricted to the
// variables that are free at the current point, followed by a projected
// backtracking line search with the Armijo condition. What it does not do is
// compute a generalised Cauchy point and minimise over the resulting subspace,
// which is how Byrd, Lu, Nocedal and Zhu let several bounds become active in
// one iteration.
//
// The consequence is that the iterate the two implementations stop at is not
// the same one. It is the same minimiser, reached along a different path and
// left at a different point inside the stopping tolerance, so the difference
// is of the order of the tolerance rather than of the answer. It is measured
// on real data rather than assumed: see the note on fit_nbinom below and the
// figure recorded for hicFindTADs in cpp/STATUS.md.
[[nodiscard]] LbfgsbResult minimise_lbfgsb(
    const std::function<double(std::span<const double>)>& objective,
    std::span<const double> x0, std::span<const Bound> bounds,
    const LbfgsbOptions& options = LbfgsbOptions());

struct NBinomFit {
    double size = 0.0;   // r
    double prob = 0.0;   // p
    int status = 0;
    int iterations = 0;
};

// fit_nbinom.fit(X): maximum likelihood estimate of the negative binomial
// parameters, with the initial values of R's fitdistr and the same bounds,
// objective and optimiser settings the Python package uses.
[[nodiscard]] NBinomFit fit_nbinom(std::span<const double> data);

// Still missing for hicDetectLoops and the cHi-C background model, and left
// out deliberately rather than stubbed: the negative binomial survival
// function needs the regularised incomplete beta function (cephes incbet),
// which nothing ported so far calls. It belongs next to gammaln and digamma in
// this file when that tool lands.

}  // namespace hicx::stats

#endif  // HICX_STATS_OPS_HPP
