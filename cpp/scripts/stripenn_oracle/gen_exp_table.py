#!/usr/bin/env python
"""Regenerates kExpSigma2 of tools/stripes_numpy.cpp: numpy.exp(-0.5 / sigma^2 * x^2)
for sigma = 2 and x = 0..8, the exponentials of the Gaussian kernel of Stripenn's
default Canny sigma (radius 8), as hexadecimal floating point literals.

numpy 1.26 evaluates the float64 exponential with a vectorised routine on CPUs
with AVX512, which differs from the C library's exp in the last bit for some
arguments. The table in the C++ source holds the values of such a CPU (the one
the comparison with Stripenn was made on). On a CPU without AVX512 this script
reports the entries that differ; that is expected and no reason to change the
table.

Usage (with the Stripenn environment's interpreter, numpy 1.26.4):
  ~/miniconda3/envs/stripenn-bench/bin/python gen_exp_table.py [path/to/stripes_numpy.cpp]
"""
import math
import re
import sys

import numpy as np

x = np.arange(0, 9, dtype=np.float64)
table = np.exp(-0.5 / 4 * x ** 2)
print("numpy", np.__version__)
print("constexpr double kExpSigma2[9] = {")
print(",\n".join("    %s" % float(v).hex() for v in table))
print("};")

# the value of the C library's exp, to show which entries depend on the CPU's routine
libm = np.array([math.exp(-0.5 / 4 * float(i) ** 2) for i in range(9)])
print("entries where numpy.exp differs from the C library's exp:", [int(i) for i in np.nonzero(table != libm)[0]])

if len(sys.argv) > 1:
    text = open(sys.argv[1]).read()
    block = text[text.index("kExpSigma2[9]"):]
    block = block[:block.index("};")]
    literals = re.findall(r"0x[0-9a-f.]+p[+-]?\d+", block)
    committed = np.array([float.fromhex(h) for h in literals])
    same = np.array_equal(committed, table)
    print("committed table equals numpy.exp on this CPU:", same)
    if not same:
        print("differing entries:", [int(i) for i in np.nonzero(committed != table)[0]])
