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
//   hicDetectLoops        ranksums (its donut tests), fit_nbinom, betainc and
//                         nbinom_sf for the per distance preselection p-value
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

// The precision the objective is evaluated in, which is **not** a detail.
//
// fit_nbinom.fit is handed the `.data` array of a scipy matrix, so its dtype
// is whatever the caller's matrix holds. Under numpy's value based promotion
// (1.26, which the reference environment pins) a float32 array combined with a
// float64 scalar stays float32, so for a float32 matrix two of the five terms
// of the log likelihood -- sum(gammaln(X + r)) and sum(X * log(1 - p)) -- are
// computed and reduced in **single precision**, while the other three stay in
// float64.
//
// The consequence is not a rounding difference, it is a different answer. At
// an objective value around 100 the float32 terms resolve to about 8e-06,
// while approx_grad probes the objective with a step of 1e-08, so the forward
// difference gradient is noise of order 1e+03 rather than the true gradient of
// order 1. scipy's L-BFGS-B then fails its line search and returns after three
// iterations, a few parts in 1e-06 away from its starting point, with
// warnflag 2 (ABNORMAL_TERMINATION_IN_LNSRCH). Measured on the 19 distance
// distributions of GSE63525_GM12878_insitu_primary_2_5mb.cool chromosome 1:
// every one of them, and the fitted `size` stays at the initial 10 while the
// true maximum likelihood value is between 120 and 220.
//
// An implementation that evaluates the same objective in float64 converges
// instead, finds an objective up to 4 percent lower, and reports a `size` 20
// times larger. It is the better fit and it is not the reference. Choosing
// Float32 here reproduces the reference, including its failure.
enum class NBinomPrecision {
    // X is an integer or float64 array; everything is float64.
    Float64,
    // X is a float32 array; the two X dependent sums, and the moment
    // estimator that seeds the optimiser, are evaluated in float32.
    Float32,
};

[[nodiscard]] NBinomFit fit_nbinom(std::span<const double> data,
                                   NBinomPrecision precision);

// --------------------------------------------------------------------------
// The negative binomial tail
//
// Added for hicDetectLoops; chicViewpointBackgroundModel and the rest of the
// cHi-C suite need the same two functions.

// The regularised incomplete beta function I_x(a, b), scipy.special.betainc.
// Translated from scipy/special/special/cephes/incbet.h: the same two
// continued fractions, the same power series for the small-b corner, the same
// switch between them and the same symmetry reflection. Only the gamma
// functions differ, being libm's rather than a second vendored table, which
// costs at most an ulp; the measured agreement against scipy over the corpus
// is recorded in tests/test_stats_ops.cpp.
//
// Returns NaN for a <= 0, b <= 0 or x outside [0, 1], as cephes does.
[[nodiscard]] double betainc(double a, double b, double x);

// The survival function of the *continuous* generalisation of the negative
// binomial that hicexplorer/lib/cnb.py defines:
//
//     cnb.cdf(x, r, p) = betainc(r, x + 1, p)
//     nbinom_sf(x, r, p) = 1 - cnb.cdf(x, r, p)
//
// hicDetectLoops.py:163 writes exactly `1 - cnb.cdf(...)`, so the subtraction
// is performed here in the same place and in float64, rather than folded into
// the incomplete beta as a betaincc would. A p-value that lands within an ulp
// of the preselection threshold decides whether a pixel becomes a loop
// candidate, so where the rounding happens is observable.
[[nodiscard]] double nbinom_sf(double x, double r, double p);

}  // namespace hicx::stats

#endif  // HICX_STATS_OPS_HPP
