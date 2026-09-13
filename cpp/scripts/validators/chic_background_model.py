#!/usr/bin/env python3
"""Acceptance of a chicViewpointBackgroundModel output (coordinator decision).

Why the fitted size and prob are not compared as numbers
--------------------------------------------------------
fit_nbinom.fit runs scipy's L-BFGS-B on a negative binomial likelihood that is
nearly flat along `size` for these distributions. The fitted parameters are not
a reproducible quantity of the reference: presenting the same values to the fit
in another order, which changes nothing but the last bits of numpy's pairwise
sums (and which the tool does itself at --threads > 1), moves the fitted size
by a median of 55 percent, and a held-out run of the Python reference lies
outside the elementwise class EN envelope of the other runs as often as the
C++ output does. Class EN is therefore withdrawn for these two columns, and
the model is judged by what it is used for. The position, max value and mean
value columns stay exact; the comparator format chic_background_model checks
them. The two validators here judge the fit.

chic_background_likelihood
--------------------------
The data each distribution was fitted to is rebuilt with the Python reference
itself (chicViewpointBackgroundModel.compute_background on the case's own
arguments, in --threads 1 order, with --truncateZeros applied). For every
distribution i and every reference run r, NLL[r, i] is fit_nbinom's objective,
evaluated in numpy exactly as fit_nbinom.fit evaluates it, at the size and prob
that run printed (the printed 12 decimals are what every consumer reads). Then,
over every distribution with a finite objective (see below):

    spread[i] = max_r NLL[r, i] - min_r NLL[r, i]
    T         = max_i spread[i]
    worst[i]  = max_r NLL[r, i]
    passes    = NLL[cpp, i] - worst[i] <= T      for every such i

T is the measured run-to-run likelihood spread of the reference on this case,
recomputed on every run of the harness from its N reference runs, never a
chosen constant. A C++ fit that reaches a better likelihood than every
reference run passes. The per distribution spread is not used as the
tolerance, because where all N runs happen to stop at the same bits it is 0,
and a fit a hundredth of a nat better than the others elsewhere would be held
to bit identity there; how many distributions would fail that stricter rule is
reported.

Degenerate distributions, and how each is handled:

  empty           --truncateZeros removed every value. There is nothing to
                  fit; fit_nbinom returns size 10 and prob nan. The C++ fields
                  must be the same text as every reference run's.
  unbounded       a value above 170 makes scipy's factorial overflow, so the
                  objective is +inf at every parameter pair and no likelihood
                  comparison exists. fit_nbinom returns its moment estimator
                  start. The C++ size and prob must lie within the closed range
                  of the reference runs' printed values, which is exact
                  equality when the runs agree.
  two or fewer distinct values
                  kept in the rule and in T, not skipped. They are where the
                  reference is least reproducible: on the one matrix edge case
                  a distribution of six values, two distinct, has three
                  reference runs at size 7.2e14 (NLL -4.748) and two at size 10
                  (NLL 3.229), a spread of 7.98 nats, which then is T for that
                  case. Because one such distribution can dominate T, the
                  metrics also give T recomputed without them and how many C++
                  fits would fail against that stricter T; that figure is
                  informational and does not decide the verdict.

Also reported, informational: the spread and the C++ excess separately for
distributions with np.var(X) > np.mean(X) (fit_nbinom's own test for its moment
start) and for the rest. Measured, the run-to-run spread lives almost entirely
in the var <= mean group, where the likelihood is flat along size.

A reference run whose own fitted parameters give a non finite likelihood makes
worst[i] infinite; the C++ passes there only if it is finite or equally non
finite, and the count is reported.

chic_background_downstream
--------------------------
Class E5. The Python chicViewpoint and the Python chicSignificantInteractions
are run on the case's matrices and reference points, once with the Python
model and once with the C++ model, and additionally with the second reference
run's model as a control. Parameters: --range min(fixateRange, 200000) on both
sides, the case's --fixateRange and --averageContactBin, --pValue 0.2
--xFoldBackground 1.5 (the values test_chicSignificantInteractions.py uses),
--combinationMode dual for two or more matrices and single for one, --threads 1.
chicSignificantInteractions looks up each merged peak's fitted distribution by
its genomic relative position, so every fitted position a peak lands on
decides a call.

A called interaction is (group path, start, end) of the significant file,
excluding the genes links. The case passes when the Jaccard index of the C++
model's calls against the Python model's is at least 0.99. A Python call set
that is empty makes the Jaccard index vacuous, and fails.

The module is also a script, used by the validators to rebuild the fitted
distributions inside the reference's environment:

    chic_background_model.py distributions OUT_NPZ -- TOOL_ARGUMENTS...
"""
from __future__ import annotations

import concurrent.futures
import math
import os
import subprocess
import sys
from pathlib import Path

JACCARD_FLOOR = 0.99
DOWNSTREAM_PVALUE = "0.2"
DOWNSTREAM_XFOLD = "1.5"
DOWNSTREAM_MAX_RANGE = 200000


# ---------------------------------------------------------------------------
# shared helpers


def _result(passed, metrics, diffs):
    from comparators.base import Result  # pylint: disable=C0415

    return Result(passed, None if not passed else ("NLL" if "tolerance_T" in metrics else "E5"),
                  metrics, diffs[:20])


def _read_model(path):
    """[(position, size text, prob text)], header excluded."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        rows.append((int(fields[0]), fields[1], fields[2]))
    return rows


def _tool_options(tool_args):
    """The chicViewpointBackgroundModel arguments the validators need.

    The options and defaults are those of chicViewpointBackgroundModel.py's
    parser; argparse's prefix matching accepts the same abbreviations.
    """
    import argparse  # pylint: disable=C0415

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--matrices", "-m", nargs="+", required=True)
    parser.add_argument("--referencePoints", "-rp", required=True)
    parser.add_argument("--averageContactBin", type=int, default=5)
    parser.add_argument("--truncateZeros", "-tz", action="store_true")
    parser.add_argument("--outFileName", "-o", default="background_model.txt")
    parser.add_argument("--threads", "-t", type=int, default=4)
    parser.add_argument("--fixateRange", "-fs", type=int, default=500000)
    return parser.parse_args(tool_args)


# ---------------------------------------------------------------------------
# the fitted distributions, rebuilt by the reference (runs in a subprocess)


def _dump_distributions(out_path, tool_args):
    import numpy as np  # pylint: disable=C0415
    from hicmatrix import HiCMatrix as hm  # pylint: disable=C0415

    from hicexplorer import chicViewpointBackgroundModel as tool  # pylint: disable=C0415
    from hicexplorer.lib import Viewpoint  # pylint: disable=C0415

    class Collect:
        item = None

        def put(self, item):
            self.item = item

    args = tool.parse_arguments().parse_args(tool_args)
    viewpoint = Viewpoint()
    reference_points, _ = viewpoint.readReferencePointFile(args.referencePoints)
    collected = {}
    bin_size = 0
    for matrix in args.matrices:
        viewpoint.hicMatrix = hm.hiCMatrix(matrix)
        bin_size = viewpoint.hicMatrix.getBinSize()
        queue = Collect()
        tool.compute_background(reference_points, viewpoint, args, queue)
        if isinstance(queue.item, str):
            raise SystemExit(queue.item)
        for position, values in queue.item[0].items():
            collected.setdefault(position, []).extend(values)
    keys = sorted(collected)
    arrays = {}
    for position in keys:
        values = np.array(collected[position], dtype=np.float64)
        if args.truncateZeros:
            values = values[values > 0.0]
        arrays[f"d{position}"] = values
    np.savez(out_path, positions=np.array(keys, dtype=np.int64) * bin_size,
             bin_size=np.int64(bin_size), **arrays)


def _distributions(context):
    import numpy as np  # pylint: disable=C0415

    target = Path(context["workdir"]) / "fitted_distributions.npz"
    if not target.exists():
        command = [str(context["py_python"]), str(Path(__file__).resolve()), "distributions",
                   str(target), "--"] + list(context["args_py"])
        completed = subprocess.run(command, cwd=str(context["workdir"]), env=context["env"],
                                   capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError("rebuilding the fitted distributions failed: "
                               + completed.stderr[-2000:])
    payload = np.load(target)
    positions = [int(value) for value in payload["positions"]]
    bin_size = int(payload["bin_size"])
    return {position: payload[f"d{position // bin_size}"] for position in positions}


def negative_log_likelihood(values, size, prob):
    """fit_nbinom.fit's objective, evaluated as fit_nbinom evaluates it."""
    import numpy as np  # pylint: disable=C0415
    from scipy.special import factorial, gammaln  # pylint: disable=C0415

    infinitesimal = np.finfo(np.float64).eps
    count = values.size
    with np.errstate(all="ignore"):
        result = np.sum(gammaln(values + size)) \
            - np.sum(np.log(factorial(values))) \
            - count * gammaln(size) \
            + count * size * np.log(prob) \
            + np.sum(values * np.log(1 - (prob if prob < 1 else 1 - infinitesimal)))
    return float(-result)


# ---------------------------------------------------------------------------
# validator 1: fit likelihood


def chic_background_likelihood(context):
    import numpy as np  # pylint: disable=C0415
    from scipy.special import factorial  # pylint: disable=C0415

    name = context["output"]
    reference_paths = [Path(context["out_py"]) / name] + [
        Path(directory) / name for directory in context["noise_dirs"]]
    candidate_path = Path(context["out_cpp"]) / name
    diffs = []
    if len(reference_paths) < 2:
        return _result(False, {}, ["the likelihood rule needs at least two reference runs; "
                                   "declare a noise source in the case"])
    for path in reference_paths + [candidate_path]:
        if not path.exists():
            return _result(False, {}, [f"missing {path}"])

    references = [_read_model(path) for path in reference_paths]
    candidate = _read_model(candidate_path)
    positions = [row[0] for row in references[0]]
    for rows in references[1:] + [candidate]:
        if [row[0] for row in rows] != positions:
            return _result(False, {}, ["the files do not list the same positions"])
    distributions = _distributions(context)
    if sorted(distributions) != positions:
        return _result(False, {}, ["the rebuilt distributions do not match the positions "
                                   "of the model file"])

    counts = {"finite_objective": 0, "empty": 0, "unbounded": 0,
              "two_or_fewer_distinct": 0, "var_above_mean": 0, "var_at_most_mean": 0,
              "reference_non_finite": 0}
    rows = []   # dict per distribution with a finite objective
    for index, position in enumerate(positions):
        values = distributions[position]
        reference_text = [(run[index][1], run[index][2]) for run in references]
        candidate_text = (candidate[index][1], candidate[index][2])
        if values.size == 0:
            counts["empty"] += 1
            if any(text != candidate_text for text in reference_text):
                diffs.append(f"position {position}: empty distribution, cpp {candidate_text} "
                             f"against the reference {sorted(set(reference_text))}")
            continue
        with np.errstate(all="ignore"):
            unbounded = not np.isfinite(np.sum(np.log(factorial(values))))
        if unbounded:
            counts["unbounded"] += 1
            for column, label in ((0, "size"), (1, "prob")):
                reference_values = [float(text[column]) for text in reference_text]
                value = float(candidate_text[column])
                if not min(reference_values) <= value <= max(reference_values):
                    diffs.append(f"position {position}: unbounded objective, cpp {label} "
                                 f"{candidate_text[column]} outside the reference range "
                                 f"[{min(reference_values)!r}, {max(reference_values)!r}]")
            continue
        counts["finite_objective"] += 1
        few = int(np.unique(values).size) <= 2
        over = bool(np.var(values) > np.mean(values))
        counts["two_or_fewer_distinct"] += few
        counts["var_above_mean" if over else "var_at_most_mean"] += 1
        nll = [negative_log_likelihood(values, float(size), float(prob))
               for size, prob in reference_text]
        nll_cpp = negative_log_likelihood(values, float(candidate_text[0]),
                                          float(candidate_text[1]))
        finite = all(math.isfinite(value) for value in nll)
        if not finite:
            counts["reference_non_finite"] += 1
        worst = max(nll)
        if not math.isfinite(worst):
            excess = -math.inf if (math.isfinite(nll_cpp) or nll_cpp == worst) else math.inf
        elif not math.isfinite(nll_cpp):
            excess = math.inf
        else:
            excess = nll_cpp - worst
        rows.append({"position": position, "few": few, "over": over, "finite": finite,
                     "spread": (worst - min(nll)) if finite else None, "best": min(nll),
                     "worst": worst, "cpp": nll_cpp, "excess": excess})

    def tolerance(selected):
        spreads = [row["spread"] for row in selected if row["spread"] is not None]
        return max(spreads) if spreads else 0.0

    def excess_max(selected):
        values = [row["excess"] for row in selected if math.isfinite(row["excess"])]
        return max(values) if values else None

    t_all = tolerance(rows)
    failing = [row for row in rows if row["excess"] > t_all]
    for row in sorted(failing, key=lambda item: -item["excess"])[:20]:
        diffs.append(f"position {row['position']}: cpp NLL {row['cpp']:.6f}, worst reference "
                     f"{row['worst']:.6f}, best {row['best']:.6f}: {row['excess']:.4g} nats "
                     f"worse than the worst, above the measured spread T = {t_all:.4g}")

    without_few = [row for row in rows if not row["few"]]
    t_without_few = tolerance(without_few)
    over_rows = [row for row in rows if row["over"]]
    under_rows = [row for row in rows if not row["over"]]
    metrics = {
        "reference_runs": len(references),
        "distributions": len(positions),
        **{f"distributions_{key}": value for key, value in counts.items()},
        "tolerance_T": t_all,
        "cpp_minus_worst_max": excess_max(rows),
        "cpp_worse_than_worst_by_more_than_T": len(failing),
        "cpp_better_than_best_reference": sum(1 for row in rows if row["cpp"] < row["best"]),
        "cpp_worse_than_worst_reference": sum(1 for row in rows if row["excess"] > 0.0),
        "informational_T_without_two_or_fewer_distinct": t_without_few,
        "informational_cpp_above_that_T": sum(1 for row in rows
                                              if row["excess"] > t_without_few),
        "informational_T_var_above_mean": tolerance(over_rows),
        "informational_T_var_at_most_mean": tolerance(under_rows),
        "informational_cpp_minus_worst_max_var_above_mean": excess_max(over_rows),
        "informational_cpp_minus_worst_max_var_at_most_mean": excess_max(under_rows),
        "informational_worse_than_worst_by_more_than_own_spread": sum(
            1 for row in rows if row["spread"] is not None and math.isfinite(row["excess"])
            and row["excess"] > row["spread"]),
    }
    return _result(not diffs, metrics, diffs)


# ---------------------------------------------------------------------------
# validator 2: downstream significance calls


def _significant_calls(path):
    import h5py  # pylint: disable=C0415

    calls = set()
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Group) and "start_list" in obj and "end_list" in obj:
                if "genes" in name.split("/"):
                    return
                for start, end in zip(obj["start_list"][()], obj["end_list"][()]):
                    calls.add((name, int(start), int(end)))
        handle.visititems(visit)
    return calls


def _interaction_pvalues(path):
    import h5py  # pylint: disable=C0415

    values = {}
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and name.endswith("/pvalue") \
                    and "genes" not in name.split("/"):
                values[name] = obj[()]
        handle.visititems(visit)
    return values


def _run_downstream(context, options, model, label):
    workdir = Path(context["workdir"]) / "downstream" / label
    workdir.mkdir(parents=True, exist_ok=True)
    python = str(context["py_python"])
    bin_dir = Path(context["repo_root"]) / "bin"
    distance = str(min(options.fixateRange, DOWNSTREAM_MAX_RANGE))
    interactions = workdir / "interactions.hdf5"
    significant = workdir / "significant.hdf5"
    commands = [
        [python, str(bin_dir / "chicViewpoint"), "--matrices", *options.matrices,
         "--referencePoints", options.referencePoints, "--backgroundModelFile", str(model),
         "--range", distance, distance, "--fixateRange", str(options.fixateRange),
         "--averageContactBin", str(options.averageContactBin),
         "-o", str(interactions), "-t", "1"],
        [python, str(bin_dir / "chicSignificantInteractions"),
         "--interactionFile", str(interactions), "--backgroundModelFile", str(model),
         "--range", distance, distance, "--fixateRange", str(options.fixateRange),
         "--pValue", DOWNSTREAM_PVALUE, "--xFoldBackground", DOWNSTREAM_XFOLD,
         "--combinationMode", "dual" if len(options.matrices) > 1 else "single",
         "--outFileNameSignificant", str(significant),
         "--outFileNameTarget", str(workdir / "target.hdf5"), "-t", "1"],
    ]
    for command in commands:
        completed = subprocess.run(command, cwd=str(workdir), env=context["env"],
                                   capture_output=True, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"{Path(command[1]).name} on the {label} model exited "
                               f"{completed.returncode}: {completed.stderr[-1500:]}")
    return _significant_calls(significant), _interaction_pvalues(interactions)


def _jaccard(left, right):
    union = left | right
    return 1.0 if not union else len(left & right) / len(union)


def chic_background_downstream(context):
    import numpy as np  # pylint: disable=C0415

    name = context["output"]
    options = _tool_options(context["args_py"])
    models = {"python": Path(context["out_py"]) / name,
              "cpp": Path(context["out_cpp"]) / name}
    if context["noise_dirs"]:
        models["python_run2"] = Path(context["noise_dirs"][0]) / name
    for label, path in models.items():
        if not path.exists():
            return _result(False, {}, [f"missing the {label} model {path}"])
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(models)) as pool:
            futures = {label: pool.submit(_run_downstream, context, options, path, label)
                       for label, path in models.items()}
            outcome = {label: future.result() for label, future in futures.items()}
    except RuntimeError as error:
        return _result(False, {}, [str(error)])

    python_calls, python_pvalues = outcome["python"]
    metrics = {"pvalue": DOWNSTREAM_PVALUE, "xfold": DOWNSTREAM_XFOLD,
               "range": min(options.fixateRange, DOWNSTREAM_MAX_RANGE),
               "combination_mode": "dual" if len(options.matrices) > 1 else "single",
               "calls_python": len(python_calls)}
    for label in [label for label in models if label != "python"]:
        calls, pvalues = outcome[label]
        differing = sum(int(np.count_nonzero(~((pvalues[key] == python_pvalues[key])
                                                | (np.isnan(pvalues[key])
                                                   & np.isnan(python_pvalues[key])))))
                        for key in python_pvalues)
        total = sum(values.size for values in python_pvalues.values())
        metrics[f"calls_{label}"] = len(calls)
        metrics[f"jaccard_{label}"] = _jaccard(python_calls, calls)
        metrics[f"only_python_vs_{label}"] = len(python_calls - calls)
        metrics[f"only_{label}"] = len(calls - python_calls)
        metrics[f"interaction_pvalues_differing_{label}"] = f"{differing} of {total}"
    diffs = []
    if not python_calls:
        diffs.append("the Python model yields no significant interactions, so the Jaccard "
                     "index would be vacuous")
    if metrics["jaccard_cpp"] < JACCARD_FLOOR:
        cpp_calls = outcome["cpp"][0]
        diffs.append(f"Jaccard {metrics['jaccard_cpp']:.6f} < {JACCARD_FLOOR}")
        diffs += [f"only with the Python model: {call}" for call in
                  sorted(python_calls - cpp_calls)[:8]]
        diffs += [f"only with the C++ model: {call}" for call in
                  sorted(cpp_calls - python_calls)[:8]]
    return _result(not diffs, metrics, diffs)


def main(argv):
    if len(argv) >= 3 and argv[0] == "distributions" and "--" in argv:
        split = argv.index("--")
        _dump_distributions(argv[1], argv[split + 1:])
        return 0
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
