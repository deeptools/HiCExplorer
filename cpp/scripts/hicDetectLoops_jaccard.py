#!/usr/bin/env python3
"""Class E5 measurement for hicDetectLoops: Jaccard of the called loops.

cpp/PLAN.md 5.1 defines E5 as "Jaccard index of the called intervals >= 0.99,
and every disagreeing call has a score within 1 % of its threshold". The
harness has no interval comparator yet (cpp/STATUS.md, "known open work"), and
adding one is not this workstream's to do, so the measurement lives here and
its result is quoted in the report.

Why it is needed even though every harness case for this tool passes at ED or
E0 with the bedgraph compared line for line: the Python's line *order* is not
reproducible. main() starts one process per chromosome and appends each result
as that process happens to finish, and for a cool input without --chromosomes
it also permutes the chromosome list by size in a way that depends on
--threads. The harness cases therefore all pass --threads 1 and a fixed
--chromosomes, where the reference's order is defined. This script covers the
case they cannot: the whole file, where only the set is comparable.

    cpp/scripts/hicDetectLoops_jaccard.py PY_BEDGRAPH CPP_BEDGRAPH [--pvalue P]

A loop is identified by its six interval fields; the seventh field is the
Wilcoxon rank sum p-value, which is the score --pValue thresholds. A loop
present in one output and not the other is reported together with
|p - pValue| / pValue, the side condition E5 puts on a disagreement.
"""
from __future__ import annotations

import argparse
import sys


def read_loops(path):
    """{(chrom, start, end, chrom, start, end): pvalue} for one bedgraph.

    A missing file means the tool found no loop and wrote nothing, which is a
    valid outcome and is treated as the empty set rather than as an error.
    """
    loops = {}
    try:
        handle = open(path)
    except FileNotFoundError:
        return loops
    with handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 7:
                raise SystemExit(f"{path}:{number}: expected 7 fields, got "
                                 f"{len(fields)}")
            key = (fields[0], int(fields[1]), int(fields[2]),
                   fields[3], int(fields[4]), int(fields[5]))
            loops[key] = float(fields[6])
    return loops


def compare(python_path, cpp_path, pvalue):
    a = read_loops(python_path)
    b = read_loops(cpp_path)
    union = set(a) | set(b)
    shared = set(a) & set(b)
    jaccard = 1.0 if not union else len(shared) / len(union)

    result = {
        "python_loops": len(a),
        "cpp_loops": len(b),
        "shared": len(shared),
        "union": len(union),
        "jaccard": jaccard,
        "only_python": [],
        "only_cpp": [],
        "score_differences": [],
    }
    for key in sorted(set(a) - shared):
        result["only_python"].append((key, a[key], abs(a[key] - pvalue) / pvalue))
    for key in sorted(set(b) - shared):
        result["only_cpp"].append((key, b[key], abs(b[key] - pvalue) / pvalue))
    for key in sorted(shared):
        if a[key] == b[key]:
            continue
        denominator = abs(a[key]) if a[key] != 0.0 else 1.0
        result["score_differences"].append((key, a[key], b[key],
                                            abs(a[key] - b[key]) / denominator))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("python_bedgraph")
    parser.add_argument("cpp_bedgraph")
    parser.add_argument("--pvalue", type=float, default=0.025,
                        help="the --pValue the run used, which is the "
                             "threshold a disagreeing call is scored against "
                             "(default: %(default)s)")
    parser.add_argument("--jaccard-floor", type=float, default=0.99)
    parser.add_argument("--margin", type=float, default=0.01,
                        help="the fraction of the threshold a disagreeing "
                             "call must sit within (default: %(default)s)")
    options = parser.parse_args(argv)

    result = compare(options.python_bedgraph, options.cpp_bedgraph,
                     options.pvalue)
    print(f"python {result['python_loops']} loops, cpp {result['cpp_loops']} "
          f"loops, {result['shared']} shared, Jaccard {result['jaccard']:.6f}")
    if result["score_differences"]:
        worst = max(d for _, _, _, d in result["score_differences"])
        print(f"  {len(result['score_differences'])} shared loops differ in "
              f"p-value, worst relative {worst:.3e}")
    else:
        print("  every shared loop has a bit identical p-value")

    failed = result["jaccard"] < options.jaccard_floor
    for label, rows in (("only in python", result["only_python"]),
                        ("only in cpp", result["only_cpp"])):
        for key, score, margin in rows:
            inside = margin <= options.margin
            failed = failed or not inside
            print(f"  {label}: {key} p={score!r} "
                  f"|p - threshold| / threshold = {margin:.3e} "
                  f"{'within' if inside else 'OUTSIDE'} the E5 margin")

    print("E5:", "fail" if failed else "pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
