// The continuous generalisation of the negative binomial distribution of
// hicexplorer/lib/cnb.py, element by element in float64:
//
//     pdf(x, r, p) = exp(gammaln(r + x) - gammaln(x + 1) - gammaln(r)
//                        + r * log(p) + xlog1py(x, -p))
//     cdf(x, r, p) = betainc(r, x + 1, p)
//
// The special functions are the ones scipy 1.14 evaluates: gammaln is cephes
// lgam (from scipy's own headers), betainc
// is Boost's ibeta (hicx::scipy::betainc). xlog1py is libm's log1p, as
// scipy's own xlog1py is. The sums are formed in the order of
// the Python expression. np.exp and np.log are libm's here; numpy 1.26 has
// its own AVX512F float64 exp and log, which differ from libm by one ulp on
// about a quarter of their arguments, so pdf can differ from 3.7.6 in the
// last bits. cdf has no such step and is exact.

#ifndef HICX_CNB_HPP
#define HICX_CNB_HPP

namespace hicx::cnb {

// scipy.special.gammaln(x) for a float64 argument (cephes lgam).
[[nodiscard]] double gammaln(double x);
// scipy.special.xlog1py(x, y) for float64 arguments: 0 for x == 0 unless y
// is NaN, x * log1p(y) otherwise.
[[nodiscard]] double xlog1py(double x, double y);

[[nodiscard]] double pdf(double x, double r, double p);
[[nodiscard]] double cdf(double x, double r, double p);

}  // namespace hicx::cnb

#endif  // HICX_CNB_HPP
