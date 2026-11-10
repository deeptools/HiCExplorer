#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/py_number.txt, the oracle of
hicx::cli::python_int and python_float (core/src/argparse.cpp): which strings
int() and float() accept and the value they give.

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_pynumber_fixture.py

Reference: Python 3.12.7. Records (text in hex):
  i <hex> <value>     int(text) == value, or  <hex> B  (valid but beyond int64),
                      or <hex> X (ValueError)
  g <hex> <bits>      float(text) as IEEE bits, or <hex> X
"""
import random
import struct
import sys


def hx(s):
    return s.encode('utf-8', 'surrogatepass').hex() or '-'


def main():
    rng = random.Random(1032)
    out = ['version ' + sys.version.split()[0]]
    atoms = (['0', '1', '2', '7', '9', '10', '00', '123', '1_0', '1__0', '_1', '1_', '+', '-', '+-', '.', '..',
              'e', 'E', 'e5', 'e-5', 'e+5', '1e5', '1E-3', '.5', '5.', '1_0.5', '1.5_5', '1e1_0', '0x10', '0b1',
              'inf', 'Inf', 'INF', 'infinity', 'Infinity', 'nan', 'NaN', 'nan_', 'in', 'infi', 'x', 'a', ' ',
              '  ', '\t', '\n', '\x0b', '\x0c', '\r', '\x1c', '\x1f', '\xa0', ' ', '　', '​',
              '٣', '१', '１', '０', '²', '①', '০', '\U0001d7ce', '٠',
              '9223372036854775807', '9223372036854775808', '-9223372036854775808', '-9223372036854775809',
              '99999999999999999999', '1e308', '1e309', '1e-400', '4.9e-324', '2.5e-324', '0.1', '1e400',
              '00012', '+7', '-7', '١٢', '１２'])
    for _ in range(40000):
        text = ''.join(rng.choice(atoms) for _ in range(rng.choice([1, 1, 2, 2, 3, 4])))
        try:
            v = int(text)
            rec = str(v) if -(1 << 63) <= v < (1 << 63) else 'B'
        except ValueError:
            rec = 'X'
        out.append('i %s %s' % (hx(text), rec))
        try:
            f = float(text)
            rec = '%016x' % struct.unpack('<Q', struct.pack('<d', f))[0]
        except ValueError:
            rec = 'X'
        out.append('g %s %s' % (hx(text), rec))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
