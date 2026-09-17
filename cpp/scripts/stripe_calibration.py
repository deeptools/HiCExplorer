#!/usr/bin/env python3
"""hicDetectStripes calibration, PLAN.md 9.3.

This script performs the whole gate fixed in PLAN.md 9.3 before the method was
chosen: it plants 60 stripes (vertical and horizontal, at 300 kb, 1 Mb and
2 Mb, at 1.5-fold and 2-fold enrichment over the local background) into a
real GM12878 10 kb matrix with a fixed seed, runs the C++ hicDetectStripes on
the planted and the unplanted matrix, and reports recall, precision, the
GSE234292 replicate Jaccard, and (when a Stripenn interpreter is given) the
Jaccard against Stripenn together with memory and CPU.

It takes every data path as an argument; nothing here is a hardcoded machine
path, and it writes no large file into the repository (cpp/AGENTS_CONTRACT.md
rule 10, and the "committed test inputs must be small" note of PLAN.md 9.3).

Usage:

    stripe_calibration.py \\
        --gm12878 GM12878_10kb.mcool::resolutions/10000 \\
        --cpp-bin BUILD/tools \\
        --out OUTDIR \\
        [--gse234292-rep1 rep1_10000.cool --gse234292-rep2 rep2_10000.cool] \\
        [--stripenn-python STRIPENN_VENV/bin/python] \\
        [--chromosomes 1 2 3 ...] [--seed 20260915] [--threads 6]

GM12878 chromosome names follow hic2cool's convention (bare numbers, no
"chr" prefix); GSE234292 uses UCSC names ("chr1"). Both are handled by name,
so pass whatever the input file actually uses via --chromosomes /
--gse234292-chromosomes.
"""

import argparse
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time

import cooler
import numpy as np
import pandas as pd

DEFAULT_GM12878_CHROMOSOMES = [str(i) for i in range(1, 23)]
DEFAULT_GSE234292_CHROMOSOMES = ["chr" + str(i) for i in range(1, 20)]

# The plant grid: PLAN.md 9.3 fixes lengths {300kb, 1Mb, 2Mb} and enrichment
# {1.5, 2}. Five replicates of each (length, enrichment, orientation) gives
# 3 * 2 * 2 * 5 = 60 total, 30 per orientation, matching "60 stripes, half
# vertical and half horizontal".
LENGTHS_BP = [300_000, 1_000_000, 2_000_000]
ENRICHMENTS = [1.5, 2.0]
REPLICATES = 5

# The default hicDetectStripes parameters used throughout this calibration,
# recorded here once so that every run (unplanted, planted, GSE234292,
# Stripenn timing) uses the same settings.
TOOL_ARGS = [
    "--minStripeLength", "100000",
    "--maxStripeLength", "3000000",
    "--stripeLengthStep", "100000",
    "--backgroundWindow", "15",
    "--backgroundGap", "2",
    "--obsExpThreshold", "2.0",
    "--zScoreThreshold", "4.0",
    "--minRawCount", "2.0",
    "--mergeWindow", "5",
    "--fdr", "0.05",
]

# The background window used to plant stripes, matching the tool's own
# defaults above so that "enrichment over the local background" means the
# same thing at planting time and at measurement time.
BG_WINDOW = 15
BG_GAP = 2


def log(message):
    print(f"[stripe_calibration] {message}", file=sys.stderr, flush=True)


def measured_run(cmd, out_path=None, stdin_text="y\n"):
    """Runs cmd under /usr/bin/time -v, returns (returncode, wall_seconds,
    cpu_seconds, peak_rss_kb, stderr).

    Not resource.getrusage(RUSAGE_CHILDREN): its ru_maxrss is a high-water
    mark for the whole calling process's lifetime, not resettable between
    calls, so a "before/after" delta around one subprocess call is only
    correct for the first such call ever made -- every later call silently
    returns the earlier call's peak once it is not exceeded again. Measured
    here: two unrelated subprocesses (an EX call and, minutes later,
    Stripenn) reported the identical peak_rss_kb, which is what exposed it.
    /usr/bin/time -v isolates exactly the one child process.

    `stdin_text` answers a yes/no prompt some external tools print (Stripenn
    asks to overwrite an existing output directory); it is otherwise unread
    and harmless."""
    with tempfile.NamedTemporaryFile(mode="r", prefix="stripe_calib_time_",
                                     suffix=".log", delete=False) as handle:
        time_log = handle.name
    start = time.time()
    result = subprocess.run(["/usr/bin/time", "-v", "-o", time_log] + list(cmd),
                            capture_output=True, text=True, input=stdin_text)
    wall = time.time() - start
    cpu = 0.0
    peak_rss_kb = 0
    try:
        with open(time_log) as handle:
            time_text = handle.read()
        user = re.search(r"User time \(seconds\): ([\d.]+)", time_text)
        sys_ = re.search(r"System time \(seconds\): ([\d.]+)", time_text)
        rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", time_text)
        cpu = (float(user.group(1)) if user else 0.0) + (float(sys_.group(1)) if sys_ else 0.0)
        peak_rss_kb = int(rss.group(1)) if rss else 0
    except OSError:
        pass
    finally:
        try:
            os.remove(time_log)
        except OSError:
            pass
    if result.returncode != 0:
        log(f"command failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr[-4000:]}")
    return result.returncode, wall, cpu, peak_rss_kb, result.stderr


def read_calls(path):
    """Parses hicDetectStripes' tab separated output into a list of dicts."""
    calls = []
    if not os.path.exists(path):
        return calls
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, a_start, a_end, orientation, e_start, e_end, enrichment, pvalue, qvalue = fields
            calls.append({
                "chrom": chrom, "anchor_start": int(a_start), "anchor_end": int(a_end),
                "orientation": orientation, "extent_start": int(e_start), "extent_end": int(e_end),
                "enrichment": float(enrichment), "pvalue": float(pvalue), "qvalue": float(qvalue),
            })
    return calls


def run_tool(cpp_bin, matrix_uri, chromosomes, out_path, threads, extra_args=None):
    cmd = [os.path.join(cpp_bin, "hicDetectStripes"), "-m", matrix_uri, "-o", out_path,
          "--chromosomes"] + list(chromosomes) + TOOL_ARGS + \
        ["--threads", str(threads)] + list(extra_args or [])
    code, wall, cpu, peak_rss_kb, stderr = measured_run(cmd)
    calls = read_calls(out_path) if code == 0 else []
    return {"returncode": code, "wall_seconds": wall, "cpu_seconds": cpu,
           "peak_rss_kb": peak_rss_kb, "calls": calls, "stderr_tail": stderr[-2000:]}


# --------------------------------------------------------------------------
# Matrix extraction and planting


def load_chromosome_pixels(clr, chrom, bin_size, band_bins):
    """The chromosome's bins (local table) and its pixels restricted to the
    band, as a dense (n_bins, band_bins) numpy array of raw counts, row a,
    column index d-1 holding the count at distance d = column - row bins.

    Uses cooler's region-restricted pixel selector (fetch), not
    Cooler.pixels()[:], which would materialise the whole genome's pixel
    table (over a billion rows for GM12878) on every call."""
    bins = clr.bins().fetch(chrom)
    n_bins = len(bins)
    offset = clr.offset(chrom)
    # Region-restricted pixel selector (fetch), not Cooler.pixels()[:], which
    # would materialise the whole genome's pixel table (over a billion rows
    # for GM12878) on every call.
    block = clr.pixels(join=False).fetch(chrom).copy()
    # fetch(chrom) restricts bin1 to the chromosome but not bin2, which can
    # spill into the next chromosome in the genome-wide bin table for an
    # anchor near the chromosome's end; restrict bin2 too, so this matches
    # the C++ loader (cpp/tools/hicDetectStripes.cpp load_cool_chromosome),
    # which keeps only pixels with both ends inside [first, last).
    block = block.loc[block["bin2_id"] < offset + n_bins]
    block["bin1_id"] -= offset
    block["bin2_id"] -= offset
    distance = block["bin2_id"] - block["bin1_id"]
    block = block.loc[(distance >= 1) & (distance <= band_bins)]
    band = np.zeros((n_bins, band_bins), dtype=np.float64)
    band[block["bin1_id"].to_numpy(), (block["bin2_id"] - block["bin1_id"] - 1).to_numpy()] = \
        block["count"].to_numpy()
    return bins.reset_index(drop=True), n_bins, band


def local_background(band, anchor, length, vertical, n_bins, window=BG_WINDOW, gap=BG_GAP):
    """The mean raw count per distance, over the flanking anchors used as the
    background at planting time, for distances 1..length. Returns an array of
    length `length`, or None when there is not enough context."""
    neighbours = []
    for offset in range(gap + 1, gap + window + 1):
        for candidate in (anchor - offset, anchor + offset):
            if 0 <= candidate < n_bins:
                if vertical and candidate - length < 0:
                    continue
                if not vertical and candidate + length > n_bins:
                    continue
                neighbours.append(candidate)
    if len(neighbours) < 4:
        return None
    rows = band[neighbours, :length]
    return rows.mean(axis=0)


def plant_one(band, n_bins, anchor, length_bins, enrichment, vertical, rng):
    """Adds Poisson((enrichment - 1) * local_background) counts to the
    pixels of one candidate stripe, in place. Returns True on success, False
    when there was not enough local context to compute a background (the
    caller should pick a different anchor)."""
    bg = local_background(band, anchor, length_bins, vertical, n_bins)
    if bg is None:
        return False
    extra_mean = np.clip((enrichment - 1.0) * bg, 0.0, None)
    extra = rng.poisson(extra_mean)
    for d in range(1, length_bins + 1):
        if vertical:
            row = anchor - d
        else:
            row = anchor
        col_index = d - 1
        band[row, col_index] += extra[d - 1]
    return True


def calls_overlap_region(calls, chrom, anchor_bin, length_bins, bin_size, vertical, tolerance_bins=2):
    """Whether any call of the same orientation on `chrom` sits within
    `tolerance_bins` of `anchor_bin`, used to keep planted anchors away from
    what the tool already calls on the unplanted matrix (PLAN.md 9.3)."""
    orientation = "vertical" if vertical else "horizontal"
    anchor_bp = anchor_bin * bin_size
    for call in calls:
        if call["chrom"] != chrom or call["orientation"] != orientation:
            continue
        if abs(call["anchor_start"] - anchor_bp) <= tolerance_bins * bin_size:
            return True
    return False


# --------------------------------------------------------------------------
# Recovery rule and set matching (PLAN.md 9.3)


def recovered(plant, calls, bin_size, tolerance_bins=2):
    """A planted stripe counts as recovered when a same-orientation call
    overlaps its anchor within 2 bins and covers at least half its length."""
    for call in calls:
        if call["chrom"] != plant["chrom"] or call["orientation"] != plant["orientation"]:
            continue
        if abs(call["anchor_start"] - plant["anchor_bp"]) > tolerance_bins * bin_size:
            continue
        call_len = call["extent_end"] - call["extent_start"]
        overlap_start = max(call["extent_start"], plant["extent_start"])
        overlap_end = min(call["extent_end"], plant["extent_end"])
        overlap = max(0, overlap_end - overlap_start)
        if overlap >= 0.5 * plant["length_bp"]:
            return call
    return None


def calls_key(call, bin_size, tolerance_bins=2):
    """A coarse key so two calls from different runs can be compared as "the
    same call": chrom, orientation, and the anchor bin rounded down to the
    tolerance window."""
    return (call["chrom"], call["orientation"],
           call["anchor_start"] // (tolerance_bins * bin_size + 1))


def jaccard(calls_a, calls_b, bin_size, tolerance_bins=2):
    """Reciprocal Jaccard between two call sets: two calls match when they
    share chrom, orientation, and their anchors are within `tolerance_bins`
    of each other. Each call matches at most one call of the other set."""
    used_b = [False] * len(calls_b)
    matched = 0
    for call_a in calls_a:
        best = None
        for index, call_b in enumerate(calls_b):
            if used_b[index]:
                continue
            if call_a["chrom"] != call_b["chrom"] or call_a["orientation"] != call_b["orientation"]:
                continue
            if abs(call_a["anchor_start"] - call_b["anchor_start"]) <= tolerance_bins * bin_size:
                best = index
                break
        if best is not None:
            used_b[best] = True
            matched += 1
    union = len(calls_a) + len(calls_b) - matched
    return matched / union if union > 0 else 1.0


# --------------------------------------------------------------------------


def build_plant_plan(chromosomes, chrom_lengths, bin_size, seed):
    """The 60 (chrom, anchor placeholder, length, enrichment, orientation)
    slots, without anchors yet: those are assigned during planting once each
    chromosome's band and the unplanted call set are known. Chromosomes are
    drawn with probability proportional to their length."""
    rng = random.Random(seed)
    weights = [chrom_lengths[c] for c in chromosomes]
    total = sum(weights)
    slots = []
    for length_bp in LENGTHS_BP:
        for enrichment in ENRICHMENTS:
            for orientation in ("horizontal", "vertical"):
                for _ in range(REPLICATES):
                    r = rng.uniform(0, total)
                    acc = 0.0
                    chrom = chromosomes[-1]
                    for c, w in zip(chromosomes, weights):
                        acc += w
                        if r <= acc:
                            chrom = c
                            break
                    slots.append({"chrom": chrom, "length_bp": length_bp,
                                 "enrichment": enrichment,
                                 "vertical": orientation == "vertical"})
    return slots


def plant_and_write(clr, chromosomes, out_uri, plan, unplanted_calls, seed, bin_size):
    """Plants every slot in `plan` (mutating chromosome bands in memory) and
    writes the resulting matrix, restricted to `chromosomes`, to `out_uri`.
    Returns the finished plant list, each with its genomic anchor and extent,
    and the count of slots that could not be placed after retrying."""
    rng = np.random.default_rng(seed)
    py_rng = random.Random(seed + 1)

    bands = {}
    bin_tables = {}
    n_bins_by_chrom = {}
    for chrom in chromosomes:
        length_bp_max = max(LENGTHS_BP)
        band_bins = length_bp_max // bin_size + BG_GAP + BG_WINDOW + 1
        bins, n_bins, band = load_chromosome_pixels(clr, chrom, bin_size, band_bins)
        bands[chrom] = band
        bin_tables[chrom] = bins
        n_bins_by_chrom[chrom] = n_bins

    used_regions = {chrom: [] for chrom in chromosomes}
    plants = []
    failures = 0
    for slot in plan:
        chrom = slot["chrom"]
        n_bins = n_bins_by_chrom[chrom]
        length_bins = slot["length_bp"] // bin_size
        vertical = slot["vertical"]
        margin = BG_GAP + BG_WINDOW + 1
        placed = False
        for _ in range(200):
            if vertical:
                anchor = py_rng.randint(length_bins + margin, n_bins - margin - 1)
            else:
                anchor = py_rng.randint(margin, n_bins - length_bins - margin - 1)
            if calls_overlap_region(unplanted_calls, chrom, anchor, length_bins, bin_size,
                                    vertical, tolerance_bins=2):
                continue
            span = (anchor - length_bins - margin, anchor + length_bins + margin) if not vertical \
                else (anchor - length_bins - margin, anchor + margin)
            collides = any(not (span[1] < lo or span[0] > hi) for lo, hi in used_regions[chrom])
            if collides:
                continue
            if not plant_one(bands[chrom], n_bins, anchor, length_bins, slot["enrichment"],
                             vertical, rng):
                continue
            used_regions[chrom].append(span)
            bins = bin_tables[chrom]
            anchor_start = int(bins.iloc[anchor]["start"])
            anchor_end = int(bins.iloc[anchor]["end"])
            if vertical:
                far = max(anchor - length_bins, 0)
                extent_start = int(bins.iloc[far]["start"])
                extent_end = anchor_end
            else:
                far = min(anchor + length_bins, n_bins - 1)
                extent_start = anchor_start
                extent_end = int(bins.iloc[far]["end"])
            plants.append({
                "chrom": chrom, "anchor_bin": anchor, "anchor_bp": anchor_start,
                "orientation": "vertical" if vertical else "horizontal",
                "length_bp": slot["length_bp"], "enrichment": slot["enrichment"],
                "extent_start": extent_start, "extent_end": extent_end,
            })
            placed = True
            break
        if not placed:
            failures += 1
            log(f"could not place a slot on {chrom} ({'vertical' if vertical else 'horizontal'}, "
                f"{slot['length_bp']} bp, {slot['enrichment']}x) after 200 tries")

    write_bands(clr, chromosomes, bands, bin_tables, n_bins_by_chrom, out_uri, bin_size)
    return plants, failures


def write_bands(clr, chromosomes, bands, bin_tables, n_bins_by_chrom, out_uri, bin_size):
    """Writes the (possibly planted) bands for `chromosomes` as a new cool
    file, autosomes only, bins renumbered from 0."""
    all_bins = []
    offset = 0
    offsets = {}
    for chrom in chromosomes:
        bins = bin_tables[chrom][["chrom", "start", "end"]].copy()
        offsets[chrom] = offset
        offset += len(bins)
        all_bins.append(bins)
    bins_df = pd.concat(all_bins, ignore_index=True)

    pixel_frames = []
    for chrom in chromosomes:
        band = bands[chrom]
        n_bins = n_bins_by_chrom[chrom]
        rows, cols = np.nonzero(band)
        if len(rows) == 0:
            continue
        counts = band[rows, cols]
        bin1 = rows + offsets[chrom]
        bin2 = rows + (cols + 1) + offsets[chrom]
        valid = (cols + 1 + rows) < n_bins
        pixel_frames.append(pd.DataFrame({
            "bin1_id": bin1[valid], "bin2_id": bin2[valid], "count": counts[valid],
        }))
    pixels_df = pd.concat(pixel_frames, ignore_index=True) if pixel_frames else \
        pd.DataFrame({"bin1_id": [], "bin2_id": [], "count": []})
    pixels_df = pixels_df.sort_values(["bin1_id", "bin2_id"]).reset_index(drop=True)
    pixels_df["count"] = pixels_df["count"].round().astype(np.int32)
    pixels_df = pixels_df.loc[pixels_df["count"] > 0]

    cooler.create_cooler(out_uri, bins_df, pixels_df, ordered=True, assembly=clr.info.get("genome-assembly"))


# --------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gm12878", required=True,
                        help="GM12878 10kb cool/mcool URI, e.g. file.mcool::resolutions/10000")
    parser.add_argument("--cpp-bin", required=True, help="directory holding hicDetectStripes")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--chromosomes", nargs="+", default=DEFAULT_GM12878_CHROMOSOMES)
    parser.add_argument("--gse234292-rep1", default=None)
    parser.add_argument("--gse234292-rep2", default=None)
    parser.add_argument("--gse234292-chromosomes", nargs="+", default=DEFAULT_GSE234292_CHROMOSOMES)
    parser.add_argument("--stripenn-python", default=None,
                        help="python interpreter with stripenn installed (a separate venv, "
                             "never the reference oracle venv)")
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    report = {"seed": args.seed, "chromosomes": args.chromosomes, "tool_args": TOOL_ARGS}

    clr = cooler.Cooler(args.gm12878)
    bin_size = clr.binsize
    report["bin_size"] = bin_size
    chrom_lengths = dict(zip(clr.chromnames, clr.chromsizes))

    # --- extract the autosome-only matrix (unplanted) ---
    log("extracting the unplanted autosome matrix")
    unplanted_uri = os.path.join(args.out, "unplanted.cool")
    bands = {}
    bin_tables = {}
    n_bins_by_chrom = {}
    for chrom in args.chromosomes:
        band_bins = max(LENGTHS_BP) // bin_size + BG_GAP + BG_WINDOW + 1
        bins, n_bins, band = load_chromosome_pixels(clr, chrom, bin_size, band_bins)
        bands[chrom] = band
        bin_tables[chrom] = bins
        n_bins_by_chrom[chrom] = n_bins
    write_bands(clr, args.chromosomes, bands, bin_tables, n_bins_by_chrom, unplanted_uri, bin_size)
    del bands  # freed before planting re-reads from the extracted file

    log("running hicDetectStripes on the unplanted matrix")
    unplanted_run = run_tool(args.cpp_bin, unplanted_uri, args.chromosomes,
                             os.path.join(args.out, "unplanted_calls.tsv"), args.threads)
    report["unplanted"] = {k: v for k, v in unplanted_run.items() if k != "calls"}
    report["unplanted"]["n_calls"] = len(unplanted_run["calls"])
    log(f"unplanted: {len(unplanted_run['calls'])} calls, "
        f"{unplanted_run['wall_seconds']:.1f}s wall, {unplanted_run['peak_rss_kb'] / 1024:.0f} MB")

    # --- plant 60 stripes and write the planted matrix ---
    log("building the plant plan and planting 60 stripes")
    plan = build_plant_plan(args.chromosomes, chrom_lengths, bin_size, args.seed)
    planted_uri = os.path.join(args.out, "planted.cool")
    plants, failures = plant_and_write(cooler.Cooler(unplanted_uri), args.chromosomes, planted_uri,
                                       plan, unplanted_run["calls"], args.seed, bin_size)
    report["plants_placed"] = len(plants)
    report["plants_failed_to_place"] = failures

    log("running hicDetectStripes on the planted matrix")
    planted_run = run_tool(args.cpp_bin, planted_uri, args.chromosomes,
                           os.path.join(args.out, "planted_calls.tsv"), args.threads)
    report["planted"] = {k: v for k, v in planted_run.items() if k != "calls"}
    report["planted"]["n_calls"] = len(planted_run["calls"])
    log(f"planted: {len(planted_run['calls'])} calls, "
        f"{planted_run['wall_seconds']:.1f}s wall, {planted_run['peak_rss_kb'] / 1024:.0f} MB")

    # --- recall and precision (PLAN.md 9.3) ---
    for plant in plants:
        call = recovered(plant, planted_run["calls"], bin_size)
        plant["recovered"] = call is not None

    def recall(length_min_bp, enrichment):
        subset = [p for p in plants if p["length_bp"] >= length_min_bp and
                 p["enrichment"] == enrichment]
        if not subset:
            return None, 0
        hits = sum(1 for p in subset if p["recovered"])
        return hits / len(subset), len(subset)

    report["recall"] = {}
    for length_bp in [300_000] + LENGTHS_BP:
        for enrichment in ENRICHMENTS:
            label = f"len>={length_bp}_x{enrichment}"
            value, n = recall(length_bp, enrichment)
            report["recall"][label] = {"recall": value, "n": n}
    gated_recall, gated_n = recall(1_000_000, 2.0)
    report["recall_gate"] = {"criterion": ">=0.8 for 2-fold stripes of 1 Mb or longer",
                             "value": gated_recall, "n": gated_n,
                             "passed": (gated_recall is not None and gated_recall >= 0.8)}

    unplanted_keys = {calls_key(c, bin_size) for c in unplanted_run["calls"]}
    new_calls = [c for c in planted_run["calls"] if calls_key(c, bin_size) not in unplanted_keys]
    plant_lookup_hits = set()
    for plant in plants:
        call = recovered(plant, new_calls, bin_size)
        if call is not None:
            plant_lookup_hits.add(id(call))
    precision_hits = sum(1 for c in new_calls if id(c) in plant_lookup_hits)
    precision = precision_hits / len(new_calls) if new_calls else None
    report["precision_gate"] = {
        "criterion": ">=0.9 of calls present on the planted matrix but not the unplanted one "
                    "are recovered plants",
        "n_new_calls": len(new_calls), "n_matched_plants": precision_hits,
        "value": precision, "passed": (precision is not None and precision >= 0.9),
    }

    report["plants"] = plants

    # --- GSE234292 reproducibility (reported, no gate) ---
    if args.gse234292_rep1 and args.gse234292_rep2:
        log("running hicDetectStripes on the GSE234292 replicates")
        rep1_run = run_tool(args.cpp_bin, args.gse234292_rep1, args.gse234292_chromosomes,
                            os.path.join(args.out, "gse234292_rep1_calls.tsv"), args.threads)
        rep2_run = run_tool(args.cpp_bin, args.gse234292_rep2, args.gse234292_chromosomes,
                            os.path.join(args.out, "gse234292_rep2_calls.tsv"), args.threads)
        rep_bin_size = cooler.Cooler(args.gse234292_rep1).binsize
        rep_jaccard = jaccard(rep1_run["calls"], rep2_run["calls"], rep_bin_size)
        report["gse234292_reproducibility"] = {
            "rep1_n_calls": len(rep1_run["calls"]), "rep2_n_calls": len(rep2_run["calls"]),
            "jaccard": rep_jaccard,
        }
        log(f"GSE234292 rep1 vs rep2 Jaccard: {rep_jaccard:.3f} "
            f"({len(rep1_run['calls'])} vs {len(rep2_run['calls'])} calls)")
    else:
        report["gse234292_reproducibility"] = None
        log("GSE234292 replicates not given, skipping reproducibility")

    # --- Stripenn agreement, memory and CPU (reported, no gate) ---
    #
    # Stripenn 1.1.65.22's `compute` writes result_filtered.tsv with columns
    # chr, pos1, pos2, chr2, pos3, pos4, length, width, Mean, maxpixel,
    # pvalue, Stripiness: a bounding box (pos1..pos2) x (pos3..pos4) on one
    # chromosome (chr == chr2 always, intra-chromosomal), with no explicit
    # orientation label. The narrower side of the box (closest to its `width`
    # column) is the anchor axis and the wider side (closest to `length`) is
    # the extent axis; whichever axis is the anchor decides horizontal
    # (pos1/pos2 narrow) versus vertical (pos3/pos4 narrow). This mapping is
    # a documented heuristic, not part of Stripenn's own interface.
    if args.stripenn_python:
        log("running Stripenn on the unplanted matrix's chromosomes")
        stripenn_out = os.path.join(args.out, "stripenn") + os.sep
        # Stripenn prompts to overwrite an existing output directory; start
        # clean so there is nothing to confirm (measured_run also answers
        # "y" defensively).
        if os.path.isdir(stripenn_out):
            import shutil as _shutil
            _shutil.rmtree(stripenn_out)
        os.makedirs(stripenn_out, exist_ok=True)
        cmd = [os.path.join(os.path.dirname(args.stripenn_python), "stripenn"), "compute",
              "--cool", unplanted_uri, "--out", stripenn_out,
              "--chrom", ",".join(args.chromosomes), "--numcores", str(args.threads)]
        code, wall, cpu, peak_rss_kb, stderr = measured_run(cmd)
        report["stripenn"] = {"returncode": code, "wall_seconds": wall, "cpu_seconds": cpu,
                              "peak_rss_kb": peak_rss_kb, "stderr_tail": stderr[-2000:],
                              "command": cmd}
        result_path = os.path.join(stripenn_out, "result_filtered.tsv")
        if code == 0 and os.path.exists(result_path):
            stripenn_table = pd.read_csv(result_path, sep="\t")
            stripenn_calls = []
            for _, row in stripenn_table.iterrows():
                width_12 = abs(int(row["pos2"]) - int(row["pos1"]))
                width_34 = abs(int(row["pos4"]) - int(row["pos3"]))
                if width_12 <= width_34:
                    vertical = False
                    anchor_start = int(row["pos1"])
                else:
                    vertical = True
                    anchor_start = int(row["pos3"])
                stripenn_calls.append({
                    "chrom": str(row["chr"]), "anchor_start": anchor_start,
                    "orientation": "vertical" if vertical else "horizontal",
                })
            report["stripenn"]["n_calls"] = len(stripenn_calls)
            report["stripenn"]["jaccard_vs_hicDetectStripes"] = jaccard(
                unplanted_run["calls"], stripenn_calls, bin_size)
            log(f"Stripenn: {len(stripenn_calls)} calls, {wall:.1f}s wall, "
                f"{peak_rss_kb / 1024:.0f} MB, Jaccard vs hicDetectStripes "
                f"(both on the unplanted matrix) "
                f"{report['stripenn']['jaccard_vs_hicDetectStripes']:.3f}")
        else:
            report["stripenn"]["n_calls"] = None
            report["stripenn"]["jaccard_vs_hicDetectStripes"] = None
            log("Stripenn did not produce result_filtered.tsv; see stderr_tail in the report")
    else:
        report["stripenn"] = None
        log("no --stripenn-python given, skipping the Stripenn comparison")

    with open(os.path.join(args.out, "report.json"), "w") as handle:
        json.dump(report, handle, indent=2, default=str)
    log(f"wrote {os.path.join(args.out, 'report.json')}")
    log(f"recall gate: {report['recall_gate']}")
    log(f"precision gate: {report['precision_gate']}")


if __name__ == "__main__":
    main()
