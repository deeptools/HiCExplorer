#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/intervaltree.txt, the oracle of
hicx/py_intervaltree.hpp.

Each case builds an intervaltree.IntervalTree from half-open integer intervals
in the way hicMergeLoops does (Interval(begin, end, running index) added one by
one), removes and queries interleaved, and records for every query the running
indices of the returned set in the set's iteration order. Run with the
interpreter whose behaviour is reproduced:

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_intervaltree_fixture.py

Reference: Python 3.12.7, intervaltree 3.1.0 (reported on the first line).

Format:
  version <python> <intervaltree>
  case <name>
  a <begin> <end>        tree.add(Interval(begin, end, next index))
  r <index>              tree.remove(the interval with that index)
  q <begin> <end> <n> <index>*n   tree.overlap(begin, end) in iteration order
"""
import importlib.metadata
import random
import sys

from intervaltree import Interval, IntervalTree


def make_case(rng, name, nops, span, width, remove_rate, query_rate, resolution):
    lines = ['case ' + name]
    tree = IntervalTree()
    alive = {}
    count = 0
    for _ in range(nops):
        r = rng.random()
        if r < remove_rate and alive:
            idx = rng.choice(sorted(alive))
            tree.remove(alive.pop(idx))
            lines.append('r %d' % idx)
        elif r < remove_rate + query_rate:
            b = rng.randint(0, span) // resolution * resolution
            e = b + rng.randint(1, width) * resolution
            hits = [iv.data for iv in tree.overlap(b, e)]
            lines.append('q %d %d %d %s' % (b, e, len(hits), ' '.join(map(str, hits))))
        else:
            b = rng.randint(0, span) // resolution * resolution
            e = b + rng.randint(1, width) * resolution
            iv = Interval(b, e, count)
            tree.add(iv)
            alive[count] = iv
            lines.append('a %d %d' % (b, e))
            count += 1
    return lines


def main():
    rng = random.Random(7062026)
    out = ['version %s %s' % (sys.version.split()[0], importlib.metadata.version('intervaltree'))]
    n = 0
    # (ops, span, width in resolution units, remove rate, query rate, resolution)
    shapes = [(30, 50, 5, 0.0, 0.2, 1), (80, 200, 20, 0.15, 0.2, 1),
              (200, 100, 10, 0.2, 0.2, 1), (300, 5000, 400, 0.25, 0.15, 1),
              (400, 3000, 6, 0.3, 0.2, 5000), (600, 60, 3, 0.35, 0.2, 5000),
              (800, 2000000, 100, 0.3, 0.2, 1000), (1500, 300, 30, 0.4, 0.15, 25000),
              (1200, 400, 8, 0.1, 0.3, 10000)]
    for shape in shapes:
        for rep in range(8 if shape[0] < 500 else 4):
            n += 1
            out += make_case(rng, 'rand%d' % n, *shape)
    print('\n'.join(out))


if __name__ == '__main__':
    main()
