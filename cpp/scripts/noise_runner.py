#!/usr/bin/env python3
"""Run a Python reference tool with a declared source of oracle noise switched on.

    noise_runner.py TOOL_SCRIPT [ARGUMENTS...]

Used by equiv.py for class EN (cpp/PLAN.md 5.7). The environment selects the
noise source and the run:

    HICX_NOISE_SHIM   the name of a shim below
    HICX_NOISE_RUN    the run number; run 0 is the unmodified tool

A shim never changes what the tool computes, only a property of the
computation that the tool itself leaves to chance, and it must say which
property that is and where the tool leaves it to chance.

fit_nbinom_order
    Presents the values of every distribution to fit_nbinom.fit in a
    different, seeded order. This is the nondeterminism the tool already has:
    chicViewpointBackgroundModel.py:173-203 extends each relative position's
    value list in the order its worker processes happen to report, so at
    --threads > 1 the order varies from run to run (measured: --threads 8
    writes different digits from --threads 1 on 987 of 1,001 lines of the
    cHi-C test data). The order changes nothing but the last bits of numpy's
    pairwise sums inside the likelihood, and L-BFGS-B on this flat likelihood
    amplifies those bits. Reordering explicitly makes the noise reproducible
    and independent of how busy the machine is.
"""
from __future__ import annotations

import os
import runpy
import sys


def _fit_nbinom_order(run):
    import numpy as np
    import fit_nbinom

    original = fit_nbinom.fit
    calls = [0]

    def fit(X, initial_params=None):
        calls[0] += 1
        values = np.asarray(X)
        if run > 0 and values.size > 1:
            order = np.random.default_rng([run, calls[0]]).permutation(values.size)
            values = values[order]
        return original(values, initial_params)

    fit_nbinom.fit = fit


SHIMS = {"fit_nbinom_order": _fit_nbinom_order}


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    shim = os.environ.get("HICX_NOISE_SHIM", "")
    run = int(os.environ.get("HICX_NOISE_RUN", "0"))
    if shim:
        if shim not in SHIMS:
            print(f"noise_runner.py: unknown shim {shim!r}; known: {sorted(SHIMS)}",
                  file=sys.stderr)
            return 2
        SHIMS[shim](run)
    tool = sys.argv[1]
    sys.argv = [tool] + sys.argv[2:]
    runpy.run_path(tool, run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
