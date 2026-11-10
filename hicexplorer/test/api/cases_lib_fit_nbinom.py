"""Equivalence cases for hicexplorer.lib.fit_nbinom (apieq.py).

The reference is the external fit_nbinom 1.2 package of the 3.7.6
environment, which HiCExplorer 4 now carries as hicexplorer.lib.fit_nbinom.
Its only public function is fit(X, initial_params=None), returning
{'size': np.float64, 'prob': np.float64}.

Inputs are the shapes the callers use: hicDetectLoops hands it the obs/exp
values of one genomic distance (float64, or float32 for a float32 matrix),
chicViewpointBackgroundModel the counts of one relative distance. The count
vectors are drawn from negative binomials with fixed seeds, the real-data
ones are the per distance rows of the cHi-C test data.

Tolerance: the fixed cases compare size and prob at a per case rtol (RTOL
below). The float64 and integer fits are the exact maximum likelihood solution
(hicx/nbinom_mle.hpp, the root of the digamma profile score), while the
reference stops L-BFGS-B inside its stopping tolerance, so the two differ by
the optimiser's residual error: 1e-7 to 8e-4 relative on these cases, measured
against the reference of the 3.7.6 environment (scipy 1.14.1). The rtol of
each case is its measured maximum rounded up to one significant digit. The
float32 cases are fitted by the projected optimiser of hicx/stats_ops.hpp,
which stops at a different point of the same tolerance (F24 of cpp/STATUS.md);
they differ by 2e-8 to 2e-4. The 'constant' case has no finite maximum
likelihood estimate (variance not above the mean): the reference stops at a
size of a few hundred on the likelihood ridge towards the Poisson limit, the
exact solution reports the limit itself (size 1e6), so that case compares the
negative log-likelihood, like the likelihood/* cases. In general, on 600 random
fits (150 seeds, float64, float32, int64 and gamma distributed obs/exp
values, 5 to 2000 values) the median relative difference of size and prob was
3e-6, the 95th percentile 5e-4 and the maximum 5 percent for the old L-BFGS-B
translation (float32 10 percent, obs/exp float32 66 percent on a flat
likelihood). The tool equivalence of hicDetectLoops and
chicViewpointBackgroundModel judges the fits by likelihood. The likelihood/*
cases do the same: the driver returns the negative log-likelihood of the fitted
parameters, evaluated in float64, and the cases compare it at rtol 2e-3.

Deviations: int8, int16, uint8, uint16, uint32, uint64, float16 and complex
arrays raise TypeError naming HiCExplorer 4 (3.7.6 fitted them with
single precision arithmetic); a list or other non-ndarray X raises the
AttributeError for the missing .size without first evaluating np.mean on it.
"""
import numpy as np

from apieq import Case

# The drivers pick fit of either side and ignore the target, but apieq resolves
# it before calling the driver, on both sides. The target must therefore name
# something importable in the reference tree and here: with the former
# 'hicexplorer' it did not resolve (a dotted name needs a module and an
# attribute), every case raised ImportError on both sides, and the comparison
# passed without fitting anything.
F = 'hicexplorer.utilities.toString'
TOL = {}


def run(target, *args, **kwargs):
    """fit of hicexplorer.lib.fit_nbinom here, of the external fit_nbinom 1.2
    package on the reference side (the 3.7.6 tree has no such module)."""
    try:
        from hicexplorer.lib.fit_nbinom import fit
    except ImportError:
        from fit_nbinom import fit
    return fit(*args, **kwargs)


def run_nll(target, X, initial_params=None):
    """The negative log-likelihood, in float64, at the parameters fit returns."""
    from scipy.special import gammaln
    d = run(target, X, initial_params)
    x = np.asarray(X, dtype=np.float64)
    r, p, n = float(d['size']), float(d['prob']), x.size
    return float(-(gammaln(x + r).sum() - gammaln(x + 1).sum() - n * gammaln(r)
                   + n * r * np.log(p) + (x * np.log(1 - p)).sum()))


def case(id, target, **kw):
    kw.setdefault('driver', 'run')
    return Case(id, target, **kw)


def nb_counts(seed, size, prob, n, dtype=np.float64):
    rng = np.random.default_rng(seed)
    return rng.negative_binomial(size, prob, n).astype(dtype)


def obs_exp(seed, n, dtype=np.float64):
    rng = np.random.default_rng(seed)
    return np.round(rng.gamma(2.0, 3.0, n), 3).astype(dtype)


CASES = [
    case('counts/float64', F, args=lambda: (nb_counts(1, 4.0, 0.3, 400),), **TOL),
    case('counts/overdispersed', F, args=lambda: (nb_counts(2, 1.5, 0.1, 800),), **TOL),
    case('counts/small_sample', F, args=lambda: (np.array([1., 2, 3, 5, 8, 2, 1, 0, 4]),), **TOL),
    case('counts/int64', F, args=lambda: (nb_counts(3, 3.0, 0.4, 300, np.int64),), **TOL),
    case('counts/int32', F, args=lambda: (nb_counts(4, 3.0, 0.4, 300, np.int32),), **TOL),
    case('obs_exp/float64', F, args=lambda: (obs_exp(5, 500),), **TOL),
    case('obs_exp/float32', F, args=lambda: (obs_exp(6, 500, np.float32),), **TOL),
    case('obs_exp/float32_small', F, args=lambda: (np.array([1, 2, 3], dtype=np.float32),), **TOL),
    case('two_dimensional', F, args=lambda: (np.array([[1., 2], [3, 5], [0, 4]]),), **TOL),
    case('non_contiguous', F, args=lambda: (nb_counts(7, 4.0, 0.3, 400)[::2],), **TOL),
    case('constant', F, args=lambda: (np.array([3., 3, 3]),), driver='run_nll', rtol=4e-3,
         note='no finite maximum likelihood estimate (variance not above the mean): compared by likelihood; '
              'the exact solution (size 1e6) has the lower value, 4.4878 against 4.5044 of the reference, 3.7e-3 relative'),
    case('large_values', F, args=lambda: (np.array([300., 2, 3, 6]),), **TOL),
    case('negative_value', F, args=lambda: (np.array([-1., 2, 3]),), **TOL),
    case('nan_value', F, args=lambda: (np.array([1., 2, np.nan]),)),
    case('empty', F, args=lambda: (np.array([]),)),
    case('factorial_overflow', F, args=lambda: (np.array([1., 2, 175.]),)),
    case('initial/tuple', F, args=lambda: (nb_counts(8, 4.0, 0.3, 300), (2, 0.5)), **TOL),
    case('initial/array', F, args=lambda: (nb_counts(8, 4.0, 0.3, 300), np.array([2.0, 0.5])), **TOL),
    case('initial/keyword', F, args=lambda: (nb_counts(8, 4.0, 0.3, 300),), kwargs={'initial_params': [5.0, 0.4]}, **TOL),
    case('initial/column', F, args=lambda: (nb_counts(8, 4.0, 0.3, 300), np.array([[2.0, 0.5]])), **TOL),
    case('initial/float32_data', F, args=lambda: (obs_exp(9, 300, np.float32), (2.0, 0.5)), **TOL),
    case('initial/out_of_bounds', F, args=lambda: (np.array([1., 2, 3, 6]), (-1, 5)), **TOL),
    case('initial/none', F, args=lambda: (np.array([1., 2, 3, 5, 8]), None), **TOL),
    case('initial/too_long', F, args=lambda: (np.array([1., 2, 3]), [1, 2, 3])),
    case('initial/too_short', F, args=lambda: (np.array([1., 2, 3, 6]), [1])),
    case('initial/empty_data', F, args=lambda: (np.array([]), (2.0, 0.5))),
    case('input/list', F, args=([1, 2, 3],)),
    case('input/scalar', F, args=(5.0,)),
    case('input/none', F, args=(None,)),
    case('input/bool_array', F, args=lambda: (np.array([1, 2, 3, 7], dtype=bool),)),
]

# Measured maximum relative difference of size and prob against the reference,
# rounded up to one significant digit (see the docstring).
RTOL = {
    'counts/float64': 7e-6, 'counts/overdispersed': 1e-6, 'counts/small_sample': 2e-7,
    'counts/int64': 3e-6, 'counts/int32': 8e-7, 'obs_exp/float64': 2e-6,
    'obs_exp/float32': 2e-4, 'obs_exp/float32_small': 2e-4, 'two_dimensional': 9e-5,
    'non_contiguous': 1e-6, 'initial/tuple': 4e-6, 'initial/array': 4e-6,
    'initial/keyword': 4e-6, 'initial/column': 4e-6, 'initial/float32_data': 3e-8,
    'initial/out_of_bounds': 8e-4, 'initial/none': 4e-5,
}
for _case in CASES:
    if _case.id in RTOL:
        _case.rtol = RTOL[_case.id]
        _case.note = ('rtol is the measured maximum difference of the fit from the reference, rounded up: '
                      'the reference optimiser stops inside its tolerance, the float64 and integer fits are '
                      'the exact maximum likelihood solution')

NLL = dict(rtol=2e-3, note='both optimisers stop inside their tolerance at different points; the objective '
                           'is the comparable quantity (F24 of cpp/STATUS.md)')
for seed in range(12):
    _rng = np.random.default_rng(100 + seed)
    _n = int(_rng.integers(5, 2000))
    _size, _prob = _rng.uniform(.5, 30), _rng.uniform(.02, .9)
    for _name, _dtype in (('float64', np.float64), ('float32', np.float32), ('int64', np.int64)):
        CASES.append(Case('likelihood/{}_{}'.format(_name, seed), F, driver='run_nll',
                          args=lambda s=_size, p=_prob, n=_n, seed=seed, d=_dtype: (nb_counts(100 + seed, s, p, n, d),),
                          **NLL))
    CASES.append(Case('likelihood/obs_exp_float32_{}'.format(seed), F, driver='run_nll',
                      args=lambda n=_n, seed=seed: (obs_exp(200 + seed, n, np.float32),), **NLL))
