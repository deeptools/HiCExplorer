// scipy.special functions exactly as scipy 1.14.1 computes them.
//
// This is a reusable core component for the places where the last bit of a
// special function is observable in a tool's output, and where scipy does not
// use the cephes routine hicx/stats_ops.hpp translates.
//
// betainc. Since scipy 1.12, scipy.special.betainc is not cephes incbet but
// Boost.Math's ibeta (scipy/special/functions.json maps betainc to
// ibeta_double in boost_special_functions.h), evaluated with scipy's
// SpecialPolicy: no promotion of float or double to long double, and at most
// 400 root iterations. The two differ in the upper tail. For
// betainc(58.254577634279, 10.8, 0.998604709282) cephes returns
// 1 - MACHEP, 0.9999999999999999, where Boost returns 1.0, and chicViewpoint
// stores 1 - betainc as its p-value, so the difference is a stored 0.0
// against 1.1e-16, which the acceptance gate rejects.
//
// The implementation calls Boost.Math itself, header only and pinned to the
// commit scipy v1.14.1 vendors as scipy/_lib/boost_math
// (a53b013c735caa98179532a32ad24d34569b9710), with the same policy and the
// same exception mapping as scipy's ibeta_wrap. Measured against
// scipy.special.betainc on every (size, x + 1, prob) triple the cHi-C
// interaction files of the test data evaluate, in float64 and float32 form,
// plus 20,000 random triples: bit identical.
//
// hicx::stats::betainc (cephes) is left as it is, because hicDetectLoops
// depends on it and was validated with it.

#ifndef HICX_SCIPY_SPECIAL_HPP
#define HICX_SCIPY_SPECIAL_HPP

namespace hicx::scipy {

// scipy.special.betainc(a, b, x) for float64 arguments: NaN for a NaN
// argument, for a <= 0, b <= 0 or x outside [0, 1]; otherwise Boost's ibeta,
// with a domain error mapped to NaN, an overflow to infinity and an underflow
// to 0, as scipy's ibeta_wrap maps them.
[[nodiscard]] double betainc(double a, double b, double x);

}  // namespace hicx::scipy

#endif  // HICX_SCIPY_SPECIAL_HPP
