#!/usr/bin/env python3
"""Reference values of the special functions in hicx/stats_ops.hpp.

Prints the C++ table that tests/test_stats_ops.cpp compares against. The values
are evaluated by mpmath at 50 digits, independently of Boost.Math, and rounded
to the nearest float64, so a test can bound the error of the Boost evaluation
in units of the last place.

    python3 cpp/scripts/special_function_reference.py
"""

import mpmath as mp

mp.mp.dps = 50


def d(value):
    return repr(float(value))


def ndtr(x):
    return mp.erfc(-mp.mpf(x) / mp.sqrt(2)) / 2


def gamma_p_inverse(a, p):
    a, p = mp.mpf(a), mp.mpf(p)
    low, high = mp.mpf(0), mp.mpf(1)
    while mp.gammainc(a, 0, high, regularized=True) < p:
        high *= 2
    for _ in range(700):
        middle = (low + high) / 2
        if mp.gammainc(a, 0, middle, regularized=True) < p:
            low = middle
        else:
            high = middle
    return (low + high) / 2


TABLES = {
    "ndtr": ([-37.0, -20.0, -6.0, -2.5, -1.0, -0.5, 0.0, 0.3, 0.7071, 1.0, 2.5, 6.0, 9.0], ndtr),
    "erf": ([0.1, 0.5, 1.0, 2.0, 5.0, -1.5], mp.erf),
    "erfc": ([0.1, 1.0, 2.0, 5.0, 10.0, -1.5], mp.erfc),
    "digamma": ([1.0, 1.5, 2.0, 0.25, 7.5, 100.0, -0.5, -2.5, 1e-5, 1e6], mp.digamma),
    "gammaln": ([0.5, 1.0, 2.0, 7.25, 100.5, 1e5, -1.5, 1e-10],
                lambda x: mp.log(abs(mp.gamma(mp.mpf(x))))),
}
TABLES2 = {
    "gamma_q": ([(0.5, 0.5), (0.5, 1.9207294103470620), (1.0, 3.0), (2.5, 0.1), (10.0, 5.0),
                 (10.0, 50.0), (0.5, 200.0), (100.0, 100.0)],
                lambda a, x: mp.gammainc(mp.mpf(a), mp.mpf(x), mp.inf, regularized=True)),
    "gamma_p_inverse": ([(0.5, 0.5), (0.5, 0.95), (0.5, 0.999), (2.5, 1e-10), (10.0, 0.5),
                         (100.0, 0.01), (1.0, 1e-20), (0.5, 1 - 2.0 ** -53)],
                        gamma_p_inverse),
    "betainc": ([(0.5, 0.5, 0.3), (1.0, 2.0, 0.25), (2.5, 3.5, 0.5), (10.0, 1.0, 0.9),
                 (0.1, 0.2, 0.05), (200.0, 3.0, 0.995), (0.001, 0.001, 0.5), (5.0, 5.0, 0.999),
                 (0.15883824665748106, 1.0, 0.00033612374545108738), (100.0, 0.5, 0.2),
                 (3.0, 700.0, 0.01), (700.0, 3.0, 0.99), (2.5, 4.0, 0.3)],
                lambda a, b, x: mp.betainc(mp.mpf(a), mp.mpf(b), 0, mp.mpf(x), regularized=True)),
}

for name, (points, function) in TABLES.items():
    print(f"// {name}")
    for x in points:
        print(f"    {{{d(x)}, {d(function(x))}}},")
for name, (points, function) in TABLES2.items():
    print(f"// {name}")
    for args in points:
        print("    {" + ", ".join(d(v) for v in args) + ", " + d(function(*args)) + "},")
