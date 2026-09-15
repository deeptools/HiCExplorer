// The count model of hicDifferentialAnalysis (cpp/PLAN.md 9.7, work item 2).
//
// Not a port: HiCExplorer has no such tool. This file holds the statistics
// only, the same for every tested unit type (TADs, TAD boundaries, loops and
// compartment bins). The tool turns contact matrices into a *family* of units,
// each unit one count per sample plus one log offset per sample, and this
// engine tests every unit for a difference between two conditions.
//
// The model, per unit, is a negative binomial generalised linear model
//
//     log mu_s = offset_s + intercept + block effects + beta * [s in B]
//
// with the edgeR quasi-likelihood pipeline (Lund et al. 2012, Chen, Lun and
// Smyth 2016) around it:
//
//  1. NB dispersion. A per-unit moment estimate (the Pearson statistic of the
//     Poisson fit set equal to its residual degrees of freedom), smoothed into
//     a trend over the unit's abundance covariate by binned 10 % trimmed means.
//     The trend, not the per-unit value, is used from here on.
//  2. Quasi-likelihood dispersion. s2 = deviance / residual df of the NB fit
//     with the trended dispersion, moderated by limma's empirical Bayes
//     (squeezeVar: a scaled F distribution fitted to s2 with a trended prior
//     over the covariate). This yields prior df d0 and posterior s2.
//  3. Test. The likelihood ratio of the condition coefficient, divided by the
//     posterior s2, is F(1, d0 + df) distributed; its signed root is a t
//     statistic. The minimum fold change is tested, not filtered: the root is
//     taken against the shifted nulls beta = +tau and beta = -tau (constrained
//     fits with tau folded into the offsets), and the p-value is limma's TREAT
//     (McCarthy and Smyth 2009) with those roots in place of Wald statistics,
//     p = P(T > r_near) + P(T > r_far). tau = 0 is the ordinary two sided test.
//
// Why quasi-likelihood rather than plug-in NB dispersions with a likelihood
// ratio test: with two or three replicates per condition the dispersion of a
// single unit is barely estimated, and a plug-in test ignores that
// uncertainty, which makes it liberal. The QL F test carries the uncertainty
// in its denominator degrees of freedom. Measured on GSE234292 (cpp/scripts/
// diff_calibration.py), the replicate nulls reach 3 to 7 % p <= 0.05 this way.
//
// Exploratory mode (one sample in a condition). The residual of a model with
// the condition term has no degrees of freedom, so steps 1 and 2 use the model
// without it: every difference between the samples, including a real one, is
// read as noise, which is why the tool labels such output exploratory. The
// unit's own s2 then *contains* the tested contrast, so the test uses the EB
// prior alone: LR / s0^2 is F(1, d0) distributed under the fitted prior.
//
// Abundance covariate. For a unit whose offset is itself a count (a local
// background, the flanks of a boundary, the contacts with the other
// compartment), the offset's Poisson noise is part of the unit's variance, so
// the tool passes -log(1 / (mean count + 0.5) + 1 / (mean offset count + 0.5))
// as the covariate. Without it the dispersion looked heterogeneous (prior df
// about 3 instead of 15 to 80 on the compartment units) and the test lost its
// power.
//
// Determinism: every per-unit step writes its own slot through
// hicx::parallel_for and every family-level reduction runs sequentially in
// unit order, so the result is byte-identical at any thread count.

#ifndef HICX_TOOLS_DIFF_ENGINE_IMPL_HPP
#define HICX_TOOLS_DIFF_ENGINE_IMPL_HPP

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace hicx::diff {

// --------------------------------------------------------------------------
// Special functions and distributions

// psi'(x) and psi''(x) for x > 0: recurrence to x >= 6, then the asymptotic
// series.
[[nodiscard]] double trigamma(double x);
[[nodiscard]] double tetragamma(double x);
// The x with trigamma(x) = y, limma's trigammaInverse (Newton, from the same
// starting point and with the same limits).
[[nodiscard]] double trigamma_inverse(double y);
// P(T > t) for Student's t with df degrees of freedom; df = infinity is the
// normal distribution.
[[nodiscard]] double student_t_sf(double t, double df);

// --------------------------------------------------------------------------
// Design

struct Design {
    // 0 for condition A, 1 for condition B, one entry per sample.
    std::vector<int> condition;
    // Empty, or one block label per sample (an additive factor).
    std::vector<std::string> block;
};

struct DesignMatrices {
    std::size_t samples = 0;
    // Intercept plus one indicator per block level after the first.
    std::size_t null_columns = 0;
    // Row-major, samples x null_columns and samples x (null_columns + 1); the
    // last column of `full` is the condition indicator.
    std::vector<double> null;
    std::vector<double> full;
    std::vector<double> condition;
    // Residual degrees of freedom of the full and of the null model.
    int full_df = 0;
    int null_df = 0;
    // True when a condition has a single sample: dispersions come from the
    // null model (see the file comment).
    bool exploratory = false;
};

// Throws std::runtime_error when the design cannot be estimated: a condition
// without samples, a block confounded with the condition, or no residual
// degrees of freedom for the model the dispersion is estimated from.
[[nodiscard]] DesignMatrices build_design(const Design& design);

// --------------------------------------------------------------------------
// The NB GLM of one unit

struct GlmFit {
    std::vector<double> beta;
    double deviance = 0.0;
    int iterations = 0;
    bool converged = false;
};

// Unit deviance of the NB distribution with dispersion phi (Poisson at 0).
[[nodiscard]] double nb_deviance(std::span<const double> y, std::span<const double> mu,
                                 double phi);

// Fisher scoring with step halving on the deviance. `x` is row-major, one row
// per sample, `columns` wide. Coefficients are kept within +-50.
[[nodiscard]] GlmFit fit_nb_glm(std::span<const double> y, std::span<const double> log_offset,
                                std::span<const double> x, std::size_t columns, double phi);

// --------------------------------------------------------------------------
// Families

struct Family {
    std::size_t samples = 0;
    // units x samples, row-major.
    std::vector<double> counts;
    std::vector<double> log_offsets;
    // One value per unit, or empty for log(mean count + 0.5).
    std::vector<double> covariate;

    [[nodiscard]] std::size_t units() const noexcept {
        return samples == 0 ? 0 : counts.size() / samples;
    }
};

struct FamilyOptions {
    // Minimum absolute log fold change tau (natural log), 0 for none.
    double min_log_fold = 0.0;
    unsigned threads = 1;
};

struct FamilyResult {
    // Per unit; NaN for a unit that was not testable (a non-finite offset, a
    // negative count, or no counts at all).
    std::vector<double> log_fold;
    std::vector<double> pvalue;
    std::vector<double> nb_dispersion;
    std::vector<double> ql_dispersion;
    double prior_df = 0.0;
    double test_df = 0.0;
    std::size_t tested = 0;
};

// Units with fewer testable members than kMinimumFamilyUnits get NaN.
inline constexpr std::size_t kMinimumFamilyUnits = 10;

[[nodiscard]] FamilyResult test_family(const Family& family, const DesignMatrices& design,
                                       const FamilyOptions& options);

// --------------------------------------------------------------------------
// Multiple testing

// Simes' combination of the non-NaN values; NaN when there are none.
[[nodiscard]] double simes(std::span<const double> pvalues);

// A smooth of values over covariate: the units sorted by covariate cut into
// min(20, n / 50) (at least 1) equal-count bins, each summarised by its median
// covariate and its mean value (10 % trimmed when `trimmed`), and linear
// interpolation between the bin centres, constant beyond them.
[[nodiscard]] std::vector<double> binned_trend(std::span<const double> covariate,
                                               std::span<const double> values, bool trimmed);

}  // namespace hicx::diff

#endif  // HICX_TOOLS_DIFF_ENGINE_IMPL_HPP
