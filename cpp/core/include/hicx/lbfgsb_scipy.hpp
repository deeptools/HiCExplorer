// L-BFGS-B as scipy 1.14 runs it, and the negative binomial fit on top of it.
//
// This is a reusable core component, separate from the projected L-BFGS in
// stats_ops.hpp, and it exists because of a measurement rather than a
// preference.
//
// hicx::stats::minimise_lbfgsb is a projected limited-memory BFGS: the same
// objective, bounds, gradient and stopping tests as scipy, but a different
// search direction and line search. For hicDetectLoops that is good enough.
// For chicViewpointBackgroundModel it is not. That tool fits a negative
// binomial to about a thousand small distributions of smoothed counts, and on
// that data the likelihood has a long flat ridge along `size`: refitting the
// same values in a different order, which changes nothing but the last bits of
// the pairwise sums, moves scipy's `size` by a median of 55 percent (measured
// over 11 orderings of the 1,001 distributions of the cHi-C test data). The
// projected BFGS lands outside that whole envelope on 12 percent of the
// distributions, because it walks the ridge along a different path: it stops at
// the starting point where scipy moves, or travels much further. No tolerance
// can absorb that, so the optimiser itself is reproduced.
//
// What is translated, and from where:
//
//   setulb, mainlb, active, bmv, cauchy, cmprlb, formk, formt, freev, hpsolb,
//   lnsrlb, matupd, projgr, subsm, dcsrch, dcstep
//       scipy/optimize/lbfgsb_src/lbfgsb.f at scipy v1.14.1, which is
//       L-BFGS-B 3.0 (Zhu, Byrd, Lu, Nocedal; the subspace minimisation of
//       Morales and Nocedal, 2011), distributed under the BSD licence
//   dpofa, dtrsl
//       scipy/optimize/lbfgsb_src/linpack.f at the same tag
//   the driver loop, including the maxiter and maxfun tests and the
//   warnflag, and the forward difference gradient
//       scipy/optimize/_lbfgsb_py.py:_minimize_lbfgsb,
//       scipy/optimize/_differentiable_functions.py:ScalarFunction (whose
//       memoisation decides how often the objective is evaluated) and
//       scipy/optimize/_numdiff.py:approx_derivative with method '2-point',
//       an absolute step and the bound aware step flip
//
// The printing and timing code is left out; nothing else is. The arithmetic
// is performed in the order the Fortran performs it, with the reference BLAS
// loops for ddot, daxpy, dscal and dcopy, which reduce sequentially. scipy
// links OpenBLAS, whose SIMD dot kernels group their sums differently, so the
// last bits of an iterate can still differ; on this ridge that is exactly the
// perturbation the reference exhibits between its own runs, which is why the
// fitted parameters are compared at class EN (cpp/PLAN.md 5.7).
//
// Determinism: no threads, no global state, no dependence on anything but the
// arguments. Independent problems may be solved concurrently.

#ifndef HICX_LBFGSB_SCIPY_HPP
#define HICX_LBFGSB_SCIPY_HPP

#include <functional>
#include <span>

#include "hicx/stats_ops.hpp"

namespace hicx::stats {

// scipy.optimize.fmin_l_bfgs_b(func, x0, approx_grad=1, bounds=bounds) with
// the options in `options`, maxls = 20 as in scipy. `status` is scipy's
// warnflag: 0 converged, 1 iteration or evaluation limit, 2 anything else
// (ABNORMAL_TERMINATION_IN_LNSRCH in practice). `function_evaluations` is
// ScalarFunction.nfev, which counts the gradient probes as well.
[[nodiscard]] LbfgsbResult minimise_lbfgsb_scipy(
    const std::function<double(std::span<const double>)>& objective,
    std::span<const double> x0, std::span<const Bound> bounds,
    const LbfgsbOptions& options = LbfgsbOptions(), int max_line_search_steps = 20);

// fit_nbinom.fit(X) for a float64 X, driven by minimise_lbfgsb_scipy. The
// objective, the moment estimator start and the bounds are those of
// fit_nbinom/fit_nbinom.py; the objective's sums use numpy's pairwise order.
//
// Two inputs have no optimisation to run, and both are returned the way the
// Python returns them rather than improved:
//   * an empty X gives size 10 and prob NaN (np.mean of an empty array is NaN,
//     so the start is [10, 10 / (10 + NaN)]), with status 2;
//   * a value above 170 makes scipy's factorial overflow, the objective is
//     +inf for every parameter pair, and scipy returns its starting point.
[[nodiscard]] NBinomFit fit_nbinom_scipy(std::span<const double> data);

}  // namespace hicx::stats

#endif  // HICX_LBFGSB_SCIPY_HPP
