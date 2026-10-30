#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/float_format.txt, the oracle of the
fixed-point and exponent float formatting helpers (find_tads / chic_viewpoint
format_fixed for '{:.nf}', compartmentalization savetxt_line for '%.18e').

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_pyfloatformat_fixture.py

Reference: Python 3.12.7. Doubles are IEEE bits in hex, text in hex.
  f <bits> <decimals> <hex of format(x, '.<decimals>f')>
  e <bits> <hex of '%.18e' % x>
"""
import math
import random
import struct
import sys


def bits(x):
    return '%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]


def rand_double(rng):
    r = rng.random()
    if r < 0.12:
        return rng.choice([0.0, -0.0, math.inf, -math.inf, math.nan, -math.nan, 0.5, 1.5, 2.5, -0.5, 0.05, 0.005,
                           1e-7, 1e15, 1e16, 1e22, 1e300, -1e300, 5e-324, 1.7976931348623157e308,
                           0.125, 0.375, 2.675, 1.005])
    if r < 0.4:
        return struct.unpack('<d', struct.pack('<Q', rng.getrandbits(64)))[0]
    if r < 0.7:
        return rng.uniform(-100, 100)
    return rng.uniform(-1, 1) * 10.0 ** rng.randint(-20, 25)


def main():
    rng = random.Random(1034)
    out = ['version ' + sys.version.split()[0]]
    for _ in range(6000):
        x = rand_double(rng)
        d = rng.choice([0, 1, 2, 3, 4, 6, 8, 12, 15, 20])
        out.append('f %s %d %s' % (bits(x), d, format(x, '.%df' % d).encode().hex()))
        out.append('e %s %s' % (bits(x), ('%.18e' % x).encode().hex()))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
