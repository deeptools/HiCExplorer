#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/pyformat.txt, the oracle of
core/io/hic2cool/src/pyformat.cpp: repr(float), repr of lists and dicts, and
the pandas DataFrame preview that extract-norms prints.

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_pyformat_fixture.py

Reference: Python 3.12.7, pandas 2.2.3 (first line). Doubles are given as the
hex of their IEEE bits, text in hex ('-' for empty).
Records:
  version <python> <pandas>
  f <bits> <repr hex>
  l <n> <int>*n <repr hex>            repr of a list of ints
  ls <n> <str hex>*n <repr hex>       repr of a list of str
  d <n> (<int> <str hex>)*n <repr hex>  repr of a dict int -> str
  df <norm hex> <n> (<chrom hex> <start> <end> <bits>)*n <str(df.head()) hex>
"""
import math
import random
import struct
import sys

import pandas as pd


def hx(s):
    b = s.encode('utf-8') if isinstance(s, str) else s
    return b.hex() if b else '-'


def bits(x):
    return '%016x' % struct.unpack('<Q', struct.pack('<d', x))[0]


def rand_double(rng):
    r = rng.random()
    if r < 0.1:
        return rng.choice([0.0, -0.0, math.inf, -math.inf, math.nan, 1.0, 0.1, 5e-324, 1.7976931348623157e308,
                           1e16, 1e15, 9999999999999998.0, 1e-4, 1e-5, 0.0001, 123456789012345678.0, 2.5, 100.0])
    if r < 0.4:
        return struct.unpack('<d', struct.pack('<Q', rng.getrandbits(64)))[0]
    if r < 0.7:
        return rng.uniform(-10, 10)
    if r < 0.85:
        return float(rng.randint(-10**6, 10**6)) / 10 ** rng.randint(0, 8)
    return rng.uniform(0, 1) * 10.0 ** rng.randint(-25, 25)


def main():
    rng = random.Random(1031)
    out = ['version %s %s' % (sys.version.split()[0], pd.__version__)]
    for _ in range(3000):
        x = rand_double(rng)
        out.append('f %s %s' % (bits(x), hx(repr(x))))
    alphabet = list('abXY0_-.:') + ['\'', '"', '\\', ' ', '\t', '\n', '\x7f', 'é', '中', '\x1b']
    for _ in range(300):
        ints = [rng.randint(-10**12, 10**12) for _ in range(rng.randint(0, 6))]
        out.append('l %d %s %s' % (len(ints), ' '.join(map(str, ints)), hx(repr(ints))))
        strs = [''.join(rng.choice(alphabet) for _ in range(rng.randint(0, 6))) for _ in range(rng.randint(0, 5))]
        out.append('ls %d %s %s' % (len(strs), ' '.join(hx(s) for s in strs), hx(repr(strs))))
        keys = rng.sample(range(0, 40), rng.randint(0, 5))
        dd = {k: ''.join(rng.choice(alphabet) for _ in range(rng.randint(0, 5))) for k in keys}
        out.append('d %d %s %s' % (len(dd), ' '.join('%d %s' % (k, hx(v)) for k, v in dd.items()), hx(repr(dd))))
    chroms = ['chr1', 'chr2', 'chrX', 'chrUn_gl000220', '1', '22', 'X', 'chrM', 'scaffold_1234']
    norms = ['KR', 'VC', 'VC_SQRT', 'SCALE', 'GW_KR', 'INTER_VC']
    for _ in range(1500):
        n = rng.randint(1, 5)
        mode = rng.random()
        values = []
        for _ in range(n):
            if mode < 0.35:
                values.append(rng.choice([math.nan, rng.uniform(0.3, 3.5)]) if rng.random() < 0.3 else rng.uniform(0.3, 3.5))
            elif mode < 0.55:
                values.append(math.nan if rng.random() < 0.5 else 1.0)
            elif mode < 0.7:
                values.append(round(rng.uniform(0, 5), rng.randint(0, 6)))
            else:
                values.append(rand_double(rng))
        rows = [(rng.choice(chroms), rng.randint(0, 3 * 10**8), rng.randint(0, 3 * 10**8)) for _ in range(n)]
        norm = rng.choice(norms)
        df = pd.DataFrame({'chrom': [r[0] for r in rows], 'start': [r[1] for r in rows],
                           'end': [r[2] for r in rows], norm: values})
        out.append('df %s %d %s %s' % (hx(norm), n,
                                       ' '.join('%s %d %d %s' % (hx(r[0]), r[1], r[2], bits(v))
                                                for r, v in zip(rows, values)),
                                       hx(str(df.head()))))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
