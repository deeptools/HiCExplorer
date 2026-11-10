#!/usr/bin/env python3
"""False-positive calibration of hicDifferentialAnalysis (cpp/PLAN.md 9.7, gate EX).

    diff_calibration.py --cpp-bin BUILD/tools --out DIR
        [--tads WT1 WT2 KO1 KO2 --domains BED]
        [--loops WT1 WT2 KO1 KO2 --loop-calls FILE [FILE ...] [--peak-width BINS]]
        [--compartments WT1 WT2 KO1 KO2 --gc-track BEDGRAPH [--plant-run-bins N]]
        [--chromosomes CHROM ...] [--threads N] [--seed N]
        [--plant-fraction F] [--min-fold-changes 1.0 1.1]

WT1 WT2 are two biological replicates of one condition and KO1 KO2 two of
another, as cool files of raw counts on one bin table. The large data stay out
of the repository; the script takes their paths. Run it with the reference
venv interpreter (cpp/AGENTS_CONTRACT.md), which provides h5py.

Designs, for every unit type and every minimum fold change:

  null wt        WT1 against WT2, --exploratory (one sample per condition)
  null ko        KO1 against KO2, --exploratory
  null swap      {WT1, KO2} against {WT2, KO1}, unpaired
  null swap/blk  the same with the genotype as block (--blocks wt ko wt ko)
  plant          WT1 and WT2 each split into two binomial halves
                 (--splitReplicates): the first halves are condition A, the
                 second halves condition B, so both conditions hold both
                 replicates at half depth and differ only by sampling. A random
                 --plant-fraction of the units is planted (--plantRegions): the
                 contacts of a unit are thinned to 1 / fold in the condition
                 drawn for it, which is an exact count operation. Unpaired, and
                 with the replicate as block (--blocks r1 r2 r1 r2).
  wt vs ko       {WT1, WT2} against {KO1, KO2}; WT1 against KO1 (--exploratory)

The gate, fixed in PLAN.md 9.7 before the implementation:

  every null     at most 1 % of the tested units called at FDR 0.05, and the
                 fraction of p <= 0.05 between 3 % and 7 %
  plants         recall >= 0.8 at 2-fold and observed FDR <= 0.10

The p-value criterion is evaluated at minimum fold change 1, the ordinary
test. A TREAT p-value is by construction conservative under a point null (it
tests |log fold change| <= tau), so at a minimum fold change above 1 the
fraction is reported and marked n/a.

Unit types (families):
  tads           TADs (Simes over the total and the distance strata) and
                 TAD boundaries; a boundary counts as planted when a planted TAD
                 boundaries. TADs are planted as the square of a TAD's contacts;
                 boundaries separately, as the rectangle of contacts crossing
                 the boundary within the window the tool reports.
  loops          united loop calls, each against its local background; planted
                 as the peak square (+- --peak-width bins) of a tested loop.
  compartments   every bin's contacts with A against B bins (GC-oriented
                 consensus compartments); planted as runs of --plant-run-bins
                 consecutive tested bins whose contacts with A bins are thinned.

Output: DIR/report.md (every table), DIR/results.json, and the tool's output
tables and logs under DIR/runs.
"""

import argparse
import json
import math
import os
import random
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_CHROMOSOMES = [f"chr{i}" for i in range(1, 20)] + ["chrX"]
FDR = 0.05


def read_table(path):
    columns, rows = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#chrom"):
                columns = line[1:].split("\t")
            elif line and not line.startswith("#"):
                rows.append(dict(zip(columns, line.split("\t"))))
    return rows


def bin_size(path):
    import h5py  # noqa: PLC0415
    with h5py.File(path, "r") as handle:
        return int(handle.attrs["bin-size"])


class Runner:
    def __init__(self, options):
        self.tool = str(Path(options.cpp_bin) / "hicDifferentialAnalysis")
        self.runs = Path(options.out) / "runs"
        self.runs.mkdir(parents=True, exist_ok=True)
        self.threads = options.threads
        self.chromosomes = options.chromosomes
        self.resources = []

    def run(self, name, command, arguments):
        prefix = self.runs / name
        argv = [self.tool, command] + [str(a) for a in arguments] + [
            "-o", str(prefix), "--threads", str(self.threads), "--chromosomes"] + self.chromosomes
        with open(str(prefix) + ".log", "w") as log:
            log.write(" ".join(argv) + "\n")
            log.flush()
            started = time.monotonic()
            process = subprocess.Popen(argv, stdout=log, stderr=log)
            _, status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
        record = {"run": name, "wall_s": time.monotonic() - started,
                  "cpu_s": usage.ru_utime + usage.ru_stime,
                  "peak_rss_mb": usage.ru_maxrss * 1024 / 1e6}
        self.resources.append(record)
        if process.returncode != 0:
            raise RuntimeError(f"{name} exited {process.returncode}, see {prefix}.log")
        return prefix


def null_stats(rows):
    tested = [r for r in rows if r["pvalue"] != "nan"]
    n = len(tested)
    low = sum(float(r["pvalue"]) <= 0.05 for r in tested)
    called = sum(r["differential"] == "1" for r in tested)
    return {"tested": n, "p_at_most_0.05": low / n if n else None, "called": called,
            "called_fraction": called / n if n else None}


def plant_stats(rows, is_planted):
    tested = [r for r in rows if r["pvalue"] != "nan"]
    planted = [r for r in tested if is_planted(r)]
    called = [r for r in tested if r["differential"] == "1"]
    true_calls = sum(1 for r in called if is_planted(r))
    return {"tested": len(tested), "planted": len(planted), "called": len(called),
            "true_calls": true_calls,
            "recall": true_calls / len(planted) if planted else None,
            "observed_fdr": (len(called) - true_calls) / len(called) if called else 0.0}


def percent(value):
    return "n/a" if value is None else f"{100 * value:.2f} %"


def verdict(ok):
    return "PASS" if ok else "FAIL"


# --------------------------------------------------------------------------
# unit types


def tads_unit(options):
    wt1, wt2, ko1, ko2 = options.tads
    resolution = bin_size(wt1)
    domains = []
    with open(options.domains) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith("#") or fields[0] not in options.chromosomes:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end - start >= 2 * resolution:
                domains.append((fields[0], start, end))

    def plant_tads(null_prefix, path, seed):
        # The square of a TAD's contacts with itself.
        rng = random.Random(seed)
        chosen = sorted(rng.sample(domains, max(1, round(options.plant_fraction * len(domains)))))
        with open(path, "w") as handle:
            for chrom, start, end in chosen:
                handle.write(f"{chrom}\t{start}\t{end}\t{rng.choice('AB')}\n")
        keys = set(chosen)
        return lambda r: (r["chrom"], int(r["start"]), int(r["end"])) in keys

    def plant_boundaries(null_prefix, path, seed):
        # The rectangle of contacts crossing a boundary within its window, as
        # the tool reports it in the null run: a change of insulation. Thinning
        # a TAD instead would change a bordering boundary by far less than the
        # planted fold, because only one flank changes.
        rows = [r for r in read_table(str(null_prefix) + "_boundaries.tsv") if r["pvalue"] != "nan"]
        rng = random.Random(seed + 7)
        chosen = rng.sample(rows, max(1, round(options.plant_fraction * len(rows))))
        chosen.sort(key=lambda r: (r["chrom"], int(r["position"])))
        with open(path, "w") as handle:
            for r in chosen:
                handle.write(f"{r['chrom']}\t{r['windowStart']}\t{r['position']}\t{r['chrom']}\t"
                             f"{r['position']}\t{r['windowEnd']}\t{rng.choice('AB')}\n")
        keys = {(r["chrom"], r["position"]) for r in chosen}
        return lambda r: (r["chrom"], r["position"]) in keys

    return {"name": "tads", "command": "tads", "matrices": (wt1, wt2, ko1, ko2),
            "extra": ["--domains", options.domains],
            "families": {"tads": "_tads.tsv", "boundaries": "_boundaries.tsv"},
            "plants": {"tads": plant_tads, "boundaries": plant_boundaries},
            "resolution": resolution}


def loops_unit(options):
    wt1, wt2, ko1, ko2 = options.loops
    resolution = bin_size(wt1)
    width = options.peak_width

    def plant_loops(null_prefix, path, seed):
        # The peak square of a united loop, as the tool reports it in the null
        # run: the loop's enrichment over its (unchanged) background.
        rows = [r for r in read_table(str(null_prefix) + "_loops.tsv") if r["pvalue"] != "nan"]
        rng = random.Random(seed + 11)
        chosen = rng.sample(rows, max(1, round(options.plant_fraction * len(rows))))
        chosen.sort(key=lambda r: (r["chrom1"], int(r["start1"]), int(r["start2"])))
        with open(path, "w") as handle:
            for r in chosen:
                s1, s2 = int(r["start1"]), int(r["start2"])
                handle.write(f"{r['chrom1']}\t{s1 - width * resolution}\t{s1 + (width + 1) * resolution}\t"
                             f"{r['chrom2']}\t{s2 - width * resolution}\t{s2 + (width + 1) * resolution}\t"
                             f"{rng.choice('AB')}\n")
        keys = {(r["chrom1"], r["start1"], r["start2"]) for r in chosen}
        return lambda r: (r["chrom1"], r["start1"], r["start2"]) in keys

    return {"name": "loops", "command": "loops", "matrices": (wt1, wt2, ko1, ko2),
            "extra": ["--loops"] + options.loop_calls + ["--peakWidth", width],
            "families": {"loops": "_loops.tsv"}, "plants": {"loops": plant_loops},
            "resolution": resolution}


def compartments_unit(options):
    wt1, wt2, ko1, ko2 = options.compartments
    resolution = bin_size(wt1)
    run = options.plant_run_bins

    def plant_bins(null_prefix, path, seed):
        # Runs of consecutive tested bins (no masked bin in between); a BED
        # region thins the region's contacts with the A compartment.
        #
        # Contacts are symmetric, so this also thins the contacts of the
        # region's A partners (a real change of the neighbours, counted as
        # false calls here). Three design rules keep every planted bin at the
        # planted fold, checked with identical samples in both conditions:
        # a run lies within one compartment (the tool leaves contacts between
        # region bins of opposite compartments unthinned, which would
        # otherwise dilute the change of transition bins), one condition per
        # chromosome (a contact between two runs planted in opposite
        # conditions would be thinned in both and cancel), and at least
        # --plant-run-gap bins between runs. A first design without these
        # rules realised as little as 0.74 of the planted log fold change for
        # a fifth of the planted A bins, which then failed the recall gate.
        rows = [r for r in read_table(str(null_prefix) + "_compartments.tsv") if r["pvalue"] != "nan"]
        gap = options.plant_run_gap
        starts = []
        for k in range(len(rows) - run + 1):
            first, last = rows[k], rows[k + run - 1]
            if first["chrom"] == last["chrom"] and \
                    int(last["start"]) - int(first["start"]) == (run - 1) * resolution and \
                    len({rows[k + i]["compartment"] for i in range(run)}) == 1:
                starts.append(k)
        rng = random.Random(seed + 13)
        condition = {chrom: rng.choice("AB") for chrom in sorted({r["chrom"] for r in rows})}
        wanted = max(1, round(options.plant_fraction * len(rows) / run))
        chosen, used = [], set()
        for k in rng.sample(starts, len(starts)):
            if len(chosen) == wanted:
                break
            if any(k + i in used for i in range(-gap, run + gap)):
                continue
            chosen.append(k)
            used.update(range(k, k + run))
        chosen.sort()
        keys = set()
        with open(path, "w") as handle:
            for k in chosen:
                first, last = rows[k], rows[k + run - 1]
                handle.write(f"{first['chrom']}\t{first['start']}\t{last['end']}\t"
                             f"{condition[first['chrom']]}\n")
                keys.update((rows[k + i]["chrom"], rows[k + i]["start"]) for i in range(run))
        return lambda r: (r["chrom"], r["start"]) in keys

    return {"name": "compartments", "command": "compartments",
            "matrices": (wt1, wt2, ko1, ko2), "extra": ["--gcTrack", options.gc_track],
            "families": {"compartments": "_compartments.tsv"},
            "plants": {"compartments": plant_bins}, "resolution": resolution}


# --------------------------------------------------------------------------


def calibrate(runner, unit, options):
    wt1, wt2, ko1, ko2 = unit["matrices"]
    results = {"nulls": [], "plants": [], "conditions": []}
    for fold_change in options.min_fold_changes:
        tag = f"{unit['name']}_mfc{fold_change:g}"
        common = unit["extra"] + ["--minFoldChange", fold_change, "--fdr", FDR]
        nulls = [
            ("null wt (exploratory)", ["-a", wt1, "-b", wt2, "--exploratory"]),
            ("null ko (exploratory)", ["-a", ko1, "-b", ko2, "--exploratory"]),
            ("null swap, unpaired", ["-a", wt1, ko2, "-b", wt2, ko1]),
            ("null swap, genotype block", ["-a", wt1, ko2, "-b", wt2, ko1,
                                           "--blocks", "wt", "ko", "wt", "ko"]),
        ]
        null_prefix = None
        for index, (design, arguments) in enumerate(nulls):
            prefix = runner.run(f"{tag}_null{index}", unit["command"], arguments + common)
            null_prefix = null_prefix or prefix
            for family, suffix in unit["families"].items():
                stats = null_stats(read_table(str(prefix) + suffix))
                stats.update(unit=unit["name"], min_fold_change=fold_change, design=design,
                             family=family)
                stats["gate_calls"] = verdict(stats["called_fraction"] <= 0.01)
                stats["gate_p"] = (verdict(0.03 <= stats["p_at_most_0.05"] <= 0.07)
                                   if fold_change == 1.0 else "n/a")
                results["nulls"].append(stats)
        for family, builder in unit["plants"].items():
            suffix = unit["families"][family]
            for fold, paired in [(f, p) for f in options.plant_folds for p in (False, True)]:
                name = f"{tag}_plant_{family}{fold:g}_{'block' if paired else 'unpaired'}"
                plant_path = runner.runs / f"{name}_regions.bed"
                truth = builder(null_prefix, plant_path, options.seed)
                arguments = ["-a", wt1, wt2, "-b", wt1, wt2,
                             "--splitReplicates", options.seed, "--plantRegions", plant_path,
                             "--plantFold", fold, "--plantSeed", options.seed + 1]
                if paired:
                    arguments += ["--blocks", "r1", "r2", "r1", "r2"]
                prefix = runner.run(name, unit["command"], arguments + common)
                stats = plant_stats(read_table(str(prefix) + suffix), truth)
                stats.update(unit=unit["name"], min_fold_change=fold_change, fold=fold,
                             design="replicate block" if paired else "unpaired", family=family)
                if fold == 2.0:
                    stats["gate"] = verdict(stats["recall"] is not None and stats["recall"] >= 0.8
                                            and stats["observed_fdr"] <= 0.10)
                else:
                    stats["gate"] = "n/a"
                results["plants"].append(stats)
        conditions = [
            ("wt vs ko, 2 against 2", ["-a", wt1, wt2, "-b", ko1, ko2]),
            ("wt1 vs ko1 (exploratory)", ["-a", wt1, "-b", ko1, "--exploratory"]),
        ]
        for index, (design, arguments) in enumerate(conditions):
            prefix = runner.run(f"{tag}_condition{index}", unit["command"], arguments + common)
            for family, suffix in unit["families"].items():
                stats = null_stats(read_table(str(prefix) + suffix))
                stats.update(unit=unit["name"], min_fold_change=fold_change, design=design,
                             family=family)
                results["conditions"].append(stats)
    return results


def render(all_results, resources, options):
    lines = ["# hicDifferentialAnalysis calibration (cpp/PLAN.md 9.7)", ""]
    lines.append(f"- chromosomes: {' '.join(options.chromosomes)}")
    lines.append(f"- seed {options.seed}, planted fraction {options.plant_fraction}, FDR {FDR}")
    lines.append("- gate: nulls call at most 1 % of tested units, p <= 0.05 in 3 % to 7 % "
                 "(at minimum fold change 1); plants at 2-fold recall >= 0.8 and observed "
                 "FDR <= 0.10")
    lines.append("")
    for unit, results in all_results.items():
        lines.append(f"## {unit}: nulls")
        lines.append("")
        lines.append("| min fold change | design | family | tested | p <= 0.05 | called | called % "
                     "| gate: calls <= 1 % | gate: p in [3, 7] % |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for s in results["nulls"]:
            lines.append(f"| {s['min_fold_change']:g} | {s['design']} | {s['family']} | {s['tested']} | "
                         f"{percent(s['p_at_most_0.05'])} | {s['called']} | "
                         f"{percent(s['called_fraction'])} | {s['gate_calls']} | {s['gate_p']} |")
        lines.append("")
        lines.append(f"## {unit}: planted differences")
        lines.append("")
        lines.append("| min fold change | fold | design | family | tested | planted | called | "
                     "true calls | recall | observed FDR | gate (2-fold) |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for s in results["plants"]:
            recall = "n/a" if s["recall"] is None else f"{s['recall']:.3f}"
            lines.append(f"| {s['min_fold_change']:g} | {s['fold']:g} | {s['design']} | {s['family']} | "
                         f"{s['tested']} | {s['planted']} | {s['called']} | {s['true_calls']} | "
                         f"{recall} | {s['observed_fdr']:.3f} | {s['gate']} |")
        lines.append("")
        lines.append(f"## {unit}: wild type against knockout")
        lines.append("")
        lines.append("| min fold change | design | family | tested | p <= 0.05 | called | called % |")
        lines.append("|---|---|---|---|---|---|---|")
        for s in results["conditions"]:
            lines.append(f"| {s['min_fold_change']:g} | {s['design']} | {s['family']} | {s['tested']} | "
                         f"{percent(s['p_at_most_0.05'])} | {s['called']} | "
                         f"{percent(s['called_fraction'])} |")
        lines.append("")
    lines.append("## Resources per run")
    lines.append("")
    lines.append(f"`--threads {options.threads}`; peak RSS and CPU time of the process (wait4).")
    lines.append("")
    lines.append("| run | wall s | CPU s | peak RSS MB |")
    lines.append("|---|---|---|---|")
    for r in resources:
        lines.append(f"| {r['run']} | {r['wall_s']:.1f} | {r['cpu_s']:.1f} | {r['peak_rss_mb']:.1f} |")
    lines.append("")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cpp-bin", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--tads", nargs=4, metavar=("WT1", "WT2", "KO1", "KO2"))
    parser.add_argument("--domains")
    parser.add_argument("--loops", nargs=4, metavar=("WT1", "WT2", "KO1", "KO2"))
    parser.add_argument("--loop-calls", nargs="+", default=[])
    parser.add_argument("--peak-width", type=int, default=1)
    parser.add_argument("--compartments", nargs=4, metavar=("WT1", "WT2", "KO1", "KO2"))
    parser.add_argument("--gc-track")
    parser.add_argument("--plant-run-bins", type=int, default=10)
    parser.add_argument("--plant-run-gap", type=int, default=30)
    parser.add_argument("--chromosomes", nargs="+", default=DEFAULT_CHROMOSOMES)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--plant-fraction", type=float, default=0.05)
    parser.add_argument("--plant-folds", type=float, nargs="+", default=[1.5, 2.0])
    parser.add_argument("--min-fold-changes", type=float, nargs="+", default=[1.0, 1.1])
    options = parser.parse_args(argv)
    units = []
    if options.tads:
        if not options.domains:
            parser.error("--tads needs --domains")
        units.append(tads_unit(options))
    if options.loops:
        if not options.loop_calls:
            parser.error("--loops needs --loop-calls")
        units.append(loops_unit(options))
    if options.compartments:
        if not options.gc_track:
            parser.error("--compartments needs --gc-track")
        units.append(compartments_unit(options))
    if not units:
        parser.error("give at least one unit type")
    Path(options.out).mkdir(parents=True, exist_ok=True)
    runner = Runner(options)
    all_results = {unit["name"]: calibrate(runner, unit, options) for unit in units}
    report = render(all_results, runner.resources, options)
    (Path(options.out) / "report.md").write_text(report)
    (Path(options.out) / "results.json").write_text(
        json.dumps({"results": all_results, "resources": runner.resources}, indent=1))
    print(report)
    failed = any(s.get("gate_calls") == "FAIL" or s.get("gate_p") == "FAIL"
                 for r in all_results.values() for s in r["nulls"])
    failed |= any(s.get("gate") == "FAIL" for r in all_results.values() for s in r["plants"])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
