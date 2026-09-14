#!/usr/bin/env python3
"""Python reference for hicDifferentialTAD's C++ only options (PLAN.md 9.7
step 1, dual mode as in section 5.8).

    py_hicDifferentialTAD_calibrated.py <hicDifferentialTAD arguments>
        [--sharedMask] [--correctForMultipleTesting {none,fdr,bonferroni}]

Without the two options this is hicexplorer.hicDifferentialTAD itself. With
them it builds what the C++ port must produce out of hicmatrix and the
unmodified Python tool:

--sharedMask
    Both matrices are loaded whole through hicmatrix. Their nan_bins (for a
    cool file the bins whose row is all zero after loading, for h5 the file's
    nan_bins) are united, every pixel in a row or column of that union is set
    to zero in both, and the two masked matrices are saved through hicmatrix
    into a temporary directory. The Python tool then runs on those files. Two
    matrices whose bin tables differ are refused (exit 1): one mask cannot
    describe both. A cool input with a weight column is refused as well,
    because saving the loaded matrix writes the balanced values and the tool
    would balance them again; no case uses one.

--correctForMultipleTesting fdr|bonferroni
    The p-values the Python tool wrote are adjusted across all TADs, separately
    for the left inter-TAD, right inter-TAD and intra-TAD test, over the values
    that are not NaN (hicx::stats::benjamini_hochberg_adjusted and
    bonferroni_adjusted, computed operation for operation as here). A test
    rejects when its adjusted p-value is at most --pValue, and --mode and
    --modeReject combine the tests as in the tool. Both files gain the three
    adjusted p-values as the last columns.

Either option adds a comment line to both headers. Run with --threads 1: the
rows are put back into the order of the domain file, which is the order the
Python writes at one process.
"""
from __future__ import annotations

import math
import os
import shutil
import sys
import tempfile

import numpy as np

CORRECTIONS = ("none", "fdr", "bonferroni")
TESTS = ("left-inter-TAD", "right-inter-TAD", "intra-TAD")


def split_options(argv):
    shared = False
    correction = "none"
    rest = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == "--sharedMask":
            shared = True
        elif token == "--correctForMultipleTesting" or token.startswith("--correctForMultipleTesting="):
            if "=" in token:
                value = token.split("=", 1)[1]
            elif i + 1 < len(argv):
                i += 1
                value = argv[i]
            else:
                fail("argument --correctForMultipleTesting: expected one argument")
            if value not in CORRECTIONS:
                fail("argument --correctForMultipleTesting: invalid choice: '{}' (choose from "
                     "'none', 'fdr', 'bonferroni')".format(value))
            correction = value
        else:
            rest.append(token)
        i += 1
    return shared, correction, rest


def fail(message):
    print("hicDifferentialTAD: error: " + message, file=sys.stderr)
    sys.exit(2)


def benjamini_hochberg_adjusted(pvalues):
    order = [i for i, p in enumerate(pvalues) if not math.isnan(p)]
    order.sort(key=lambda i: pvalues[i])
    m = float(len(order))
    adjusted = [pvalues[i] * m / float(k + 1) for k, i in enumerate(order)]
    for k in range(len(adjusted) - 1, 0, -1):
        adjusted[k - 1] = min(adjusted[k - 1], adjusted[k])
    out = list(pvalues)
    for k, i in enumerate(order):
        out[i] = min(adjusted[k], 1.0)
    return out


def bonferroni_adjusted(pvalues):
    m = float(sum(1 for p in pvalues if not math.isnan(p)))
    return [p if math.isnan(p) else min(p * m, 1.0) for p in pvalues]


def masked_copies(target, control, workdir):
    """The two matrices with the union of their invalid bins zeroed, saved
    through hicmatrix; returns the new paths and the number of masked bins."""
    import cooler
    from scipy.sparse import diags
    from hicmatrix import HiCMatrix as hm
    from hicexplorer.utilities import check_cooler

    for path in (target, control):
        if check_cooler(path) and "weight" in cooler.Cooler(path).bins().columns:
            print("hicDifferentialTAD reference: a cool input with a weight column is not "
                  "supported: " + path, file=sys.stderr)
            sys.exit(3)
    matrices = [hm.hiCMatrix(target), hm.hiCMatrix(control)]
    bins = [[(str(c), int(s), int(e)) for c, s, e, *_ in ma.cut_intervals] for ma in matrices]
    if bins[0] != bins[1]:
        print("ERROR:hicexplorer.hicDifferentialTAD:--sharedMask needs the target and the "
              "control matrix on the same bins", file=sys.stderr)
        sys.exit(1)
    invalid = set()
    for ma in matrices:
        if ma.nan_bins is not None:
            invalid.update(int(b) for b in np.asarray(ma.nan_bins).ravel())
    keep = np.ones(len(bins[0]))
    keep[sorted(invalid)] = 0.0
    scale = diags(keep)
    paths = []
    for name, path, ma in (("target", target, matrices[0]), ("control", control, matrices[1])):
        ma.matrix = (scale @ ma.matrix @ scale).tocsr()
        ma.matrix.eliminate_zeros()
        suffix = ".cool" if check_cooler(path) else ".h5"
        out = os.path.join(workdir, name + suffix)
        ma.save(out)
        paths.append(out)
    return paths, len(invalid)


def main(argv):
    shared, correction, rest = split_options(argv)
    from hicexplorer import hicDifferentialTAD as tool

    if not shared and correction == "none":
        tool.main(rest)
        return 0

    args = tool.parse_arguments().parse_args(rest)
    workdir = tempfile.mkdtemp(prefix="difftad-reference-")
    try:
        target, control = args.targetMatrix, args.controlMatrix
        masked = None
        if shared and target is not None and control is not None:
            from hicexplorer.utilities import check_cooler
            if check_cooler(target) == check_cooler(control):
                (target, control), masked = masked_copies(target, control, workdir)
        prefix = os.path.join(workdir, "diff")
        call = []
        for option, value in (("-tm", target), ("-cm", control), ("-td", args.tadDomains)):
            if value is not None:
                call += [option, value]
        call += ["-o", prefix, "-p", repr(args.pValue), "-m", args.mode, "-mr", args.modeReject,
                 "-t", str(args.threads)]
        tool.main(call)
        write_outputs(tool, args, prefix, masked, correction)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


def write_outputs(tool, args, prefix, masked, correction):
    headers = {}
    rows = []
    for kind in ("accepted", "rejected"):
        with open(prefix + "_" + kind + ".diff_tad") as handle:
            lines = handle.read().split("\n")
        headers[kind] = [line for line in lines if line.startswith("#")]
        rows += [line.split("\t") for line in lines if line and not line.startswith("#")]

    # Back into the order of the domain file.
    domains = tool.readDomainBoundaries(args.tadDomains).values.tolist()
    positions = {}
    for index, domain in enumerate(domains):
        positions.setdefault(tuple(map(str, domain)), []).append(index)
    rows.sort(key=lambda fields: positions[tuple(fields[:6])].pop(0))

    pvalues = [[float(fields[6 + t]) for fields in rows] for t in range(3)]
    if correction == "fdr":
        adjusted = [benjamini_hochberg_adjusted(column) for column in pvalues]
    elif correction == "bonferroni":
        adjusted = [bonferroni_adjusted(column) for column in pvalues]
    else:
        adjusted = pvalues
    reject = [[not math.isnan(p) and p <= args.pValue for p in column] for column in adjusted]

    extra = []
    if masked is not None:
        extra.append("# Shared bin mask: {} bins invalid in the target or the control matrix are "
                     "masked in both\n".format(masked))
    if correction != "none":
        extra.append("# Multiple testing correction: {}, across all TADs and separately for each "
                     "test; the p-value threshold applies to the adjusted p-values\n"
                     .format(correction))
    outputs = {"accepted": "", "rejected": ""}
    for kind in outputs:
        header = headers[kind]
        text = header[0] + "\n" + header[1] + "\n" + header[2] + "\n" + "".join(extra)
        columns = header[3]
        if correction != "none":
            columns += "".join("\tadjusted p-value " + test for test in TESTS)
        outputs[kind] = text + columns + "\n"
    reject_all = args.modeReject == "all"
    for k, fields in enumerate(rows):
        left, right, intra = reject[0][k], reject[1][k], reject[2][k]
        if args.mode == "intra-TAD":
            mask = intra
        elif args.mode == "left-inter-TAD":
            mask = (left and intra) if reject_all else (left or intra)
        elif args.mode == "right-inter-TAD":
            mask = (intra and right) if reject_all else (intra or right)
        else:
            mask = (left and right and intra) if reject_all else (left or right or intra)
        line = "\t".join(fields[:12])
        if correction != "none":
            line += "".join("\t" + str(np.float64(adjusted[t][k])) for t in range(3))
        outputs["rejected" if mask else "accepted"] += line + "\n"
    for kind, text in outputs.items():
        with open(args.outFileNamePrefix + "_" + kind + ".diff_tad", "w") as handle:
            handle.write(text)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
