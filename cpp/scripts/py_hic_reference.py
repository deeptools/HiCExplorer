#!/usr/bin/env python3
"""The Python reference of the hicConvertFormat directions the Python tool
does not have (cpp/PLAN.md tier 9, item 9.1, class EX).

    py_hic_reference.py <hicConvertFormat arguments>

Run it as equiv.py runs every py_script: with the reference interpreter and
the upstream reference tree first on PYTHONPATH (reference_tree.py).

- --inputFormat hic --outputFormat cool: hicConvertFormat itself.
- --inputFormat hic --outputFormat mcool without --resolutions: hicConvertFormat
  with --outputFormat cool, whose hic2cool call writes every resolution into the
  .mcool name it is given.
- --inputFormat hic with h5, homer, ginteractions or hicpro: hic2cool writes the
  one resolution of --resolutions into a temporary cool file, and
  hicConvertFormat converts that file with --inputFormat cool. That is the route
  a Python user takes in two steps.
- --outputFormat hic: Python has no .hic writer. The script loads every input
  matrix the way hicConvertFormat loads it (MatrixFileHandler,
  hicConvertFormat.py:229-242; every resolution group of a multi resolution
  cool file, or those in --resolutions) and writes, at the output path, an
  .npz file with what the .hic file must hold: the version, the
  normalizations, the chromosomes, the resolutions and, per resolution, the
  non-zero finite pixels of the upper triangle. A cool input is loaded with
  raw counts: the C++ tool never applies a balancing weight to the pixels of
  a .hic file. Coarser --resolutions of a single matrix are summed from its
  pixels in row order, in float32. comparators/hic.py reads the .hic file the
  C++ tool writes and compares it with this file.

  --hicNormalizations auto (the default) carries the normalization columns a
  cool input stores (cool2hic's rule: Juicer names, or every column of a file
  hic2cool wrote, never 'weight') and computes VC, VC_SQRT, KR and SCALE only
  when there are none. Carried vectors exist only at the stored resolutions,
  and comparators/hic.py checks every expected name at every resolution, so
  the expected list holds the carried names only when no resolution is added;
  otherwise it is empty and the carried vectors go unchecked here.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile

import h5py
import numpy as np
from hic2cool import hic2cool_convert
from hicmatrix.lib import MatrixFileHandler
from scipy.sparse import triu

from hicexplorer import hicConvertFormat

NORMALIZATION_ORDER = ["VC", "VC_SQRT", "KR", "SCALE"]


def parse(argv):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--matrices", "-m", nargs="+")
    parser.add_argument("--outFileName", "-o", nargs="+")
    parser.add_argument("--inputFormat")
    parser.add_argument("--outputFormat", default="cool")
    parser.add_argument("--resolutions", "-r", nargs="+")
    parser.add_argument("--correction_name", default="weight")
    parser.add_argument("--correction_division", action="store_true")
    parser.add_argument("--store_applied_correction", action="store_true")
    parser.add_argument("--enforce_integer", action="store_true")
    parser.add_argument("--load_raw_values", action="store_true")
    parser.add_argument("--chromosome")
    parser.add_argument("--bedFileHicpro", "-bf", nargs="+")
    parser.add_argument("--hicVersion", default="8")
    parser.add_argument("--hicNormalizations", nargs="+", default=["auto"])
    parser.add_argument("--threads", default="1")
    return parser.parse_args(argv)


def passthrough(args):
    options = ["--correction_name", args.correction_name]
    if args.chromosome:
        options += ["--chromosome", args.chromosome]
    for flag in ("correction_division", "store_applied_correction", "enforce_integer", "load_raw_values"):
        if getattr(args, flag):
            options.append("--" + flag)
    return options


def hic_input(args, argv):
    if args.outputFormat == "cool":
        return hicConvertFormat.main(argv)
    if args.outputFormat == "mcool":
        if args.resolutions:
            sys.exit("the reference covers --outputFormat mcool without --resolutions only")
        for matrix, out in zip(args.matrices, args.outFileName):
            hicConvertFormat.main(["--matrices", matrix, "--outFileName", out,
                                   "--inputFormat", "hic", "--outputFormat", "cool"])
        return 0
    if not args.resolutions or len(args.resolutions) != 1:
        sys.exit(1)
    for index, (matrix, out) in enumerate(zip(args.matrices, args.outFileName)):
        directory = tempfile.mkdtemp()
        try:
            cool = os.path.join(directory, "matrix.cool")
            hic2cool_convert(matrix, cool, int(args.resolutions[0]), silent=True)
            call = ["--matrices", cool, "--outFileName", out, "--inputFormat", "cool",
                    "--outputFormat", args.outputFormat] + passthrough(args)
            if args.outputFormat == "hicpro":
                call += ["--bedFileHicpro", args.bedFileHicpro[index]]
            hicConvertFormat.main(call)
        finally:
            shutil.rmtree(directory, ignore_errors=True)
    return 0


def load(args, uri):
    handler = MatrixFileHandler(pFileType=args.inputFormat, pMatrixFile=uri,
                                pCorrectionFactorTable=args.correction_name,
                                pCorrectionOperator="/" if args.correction_division else None,
                                pChrnameList=[args.chromosome] if args.chromosome else None,
                                pEnforceInteger=args.enforce_integer,
                                pApplyCorrectionCoolerLoad=False)
    matrix, cut_intervals, _nan_bins, _distance_counts, _correction_factors = handler.load()
    return matrix, cut_intervals


def layout(cut_intervals):
    names, first, lengths = [], [], []
    bin_size = max(int(end) - int(start) for _, start, end, *_ in cut_intervals)
    for index, (chrom, _start, end, *_rest) in enumerate(cut_intervals):
        if not names or names[-1] != chrom:
            names.append(chrom)
            first.append(index)
            lengths.append(0)
        lengths[-1] = int(end)
    chrom_of_bin = np.zeros(len(cut_intervals), dtype=np.int64)
    local = np.zeros(len(cut_intervals), dtype=np.int64)
    bounds = first + [len(cut_intervals)]
    for c in range(len(names)):
        chrom_of_bin[bounds[c]:bounds[c + 1]] = c
        local[bounds[c]:bounds[c + 1]] = np.arange(bounds[c + 1] - bounds[c])
    return names, lengths, bin_size, chrom_of_bin, local


def upper_pixels(matrix, chrom_of_bin, local):
    upper = triu(matrix, k=0, format="coo")
    order = np.lexsort((upper.col, upper.row))
    row, col, data = upper.row[order], upper.col[order], upper.data[order]
    keep = np.isfinite(data) & (data != 0)
    row, col, data = row[keep], col[keep], data[keep]
    return (chrom_of_bin[row], chrom_of_bin[col], local[row], local[col], data.astype(np.float32))


def coarsen(pixels, factor):
    chr1, chr2, bin1, bin2, count = pixels
    keys = np.stack([chr1, chr2, bin1 // factor, bin2 // factor], axis=1)
    unique, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    sums = np.zeros(len(unique), dtype=np.float32)
    np.add.at(sums, inverse.ravel(), count)
    return unique[:, 0], unique[:, 1], unique[:, 2], unique[:, 3], sums


JUICER_NORMS = {"VC", "VC_SQRT", "KR", "SCALE", "GW_KR", "GW_VC", "GW_SCALE", "INTER_KR",
                "INTER_VC", "INTER_SCALE"}


def carried_norms(uris):
    """The bins columns cool2hic carries out of these coolers, in its order."""
    carried = []
    for uri in uris:
        path, _, group = uri.partition("::")
        with h5py.File(path, "r") as handle:
            root = handle[group or "/"]
            generated = root.attrs.get("generated-by", "")
            if isinstance(generated, bytes):
                generated = generated.decode()
            n_bins = int(root["indexes/chrom_offset"][-1])
            for column in root["bins"].keys():
                if column in ("chrom", "start", "end", "weight"):
                    continue
                if not str(generated).startswith("hic2cool") and column not in JUICER_NORMS:
                    continue
                item = root["bins"][column]
                if not isinstance(item, h5py.Dataset) or len(item) != n_bins:
                    continue
                if column not in carried:
                    carried.append(column)
    return carried


def hic_output(args):
    computed = [] if "none" in args.hicNormalizations else \
        [n for n in NORMALIZATION_ORDER if n in args.hicNormalizations]
    auto = "auto" in args.hicNormalizations and "none" not in args.hicNormalizations
    extra = [int(r) for r in (args.resolutions or [])]
    for index, (matrix_path, out) in enumerate(zip(args.matrices, args.outFileName)):
        levels = []
        groups = []
        if args.inputFormat == "cool" and "::" not in matrix_path:
            with h5py.File(matrix_path, "r") as handle:
                if "resolutions" in handle:
                    groups = sorted(handle["resolutions"].keys(), key=int)
        uris = []
        if groups:
            for group in groups:
                if extra and int(group) not in extra:
                    continue
                uris.append(matrix_path + "::/resolutions/" + group)
        else:
            uris.append(matrix_path)
        for uri in uris:
            levels.append(load(args, uri))
        wanted = computed
        if auto:
            # The C++ tool streams a cool input through cool2hic, except a
            # --chromosome cut and a subset of several mcool groups, which it
            # loads; a loaded matrix stores no vectors.
            streamed = args.inputFormat == "cool" and not args.chromosome and \
                (len(uris) == 1 or len(uris) == len(groups))
            carried = carried_norms(uris) if streamed else []
            if not carried:
                wanted = list(NORMALIZATION_ORDER)
            elif groups or not [r for r in extra if r not in {layout(c)[2] for _, c in levels}]:
                wanted = carried
            else:
                wanted = []
        arrays = {"version": np.int64(int(args.hicVersion)),
                  "normalizations": np.array(wanted, dtype="U16")}
        resolutions = {}
        names = lengths = None
        for matrix, cut_intervals in levels:
            level_names, level_lengths, bin_size, chrom_of_bin, local = layout(cut_intervals)
            if names is None:
                names, lengths = level_names, level_lengths
            else:
                # The largest end of a chromosome's last bin over the resolutions.
                lengths = [max(a, b) for a, b in zip(lengths, level_lengths)]
            resolutions[bin_size] = upper_pixels(matrix, chrom_of_bin, local)
        if not groups:
            base = next(iter(resolutions))
            for resolution in extra:
                if resolution != base:
                    resolutions[resolution] = coarsen(resolutions[base], resolution // base)
        arrays["chrom_names"] = np.array(names, dtype="U64")
        arrays["chrom_lengths"] = np.array(lengths, dtype=np.int64)
        arrays["resolutions"] = np.array(sorted(resolutions, reverse=True), dtype=np.int64)
        for resolution, (chr1, chr2, bin1, bin2, count) in resolutions.items():
            arrays[f"r{resolution}_chr1"] = np.asarray(chr1, dtype=np.int64)
            arrays[f"r{resolution}_chr2"] = np.asarray(chr2, dtype=np.int64)
            arrays[f"r{resolution}_bin1"] = np.asarray(bin1, dtype=np.int64)
            arrays[f"r{resolution}_bin2"] = np.asarray(bin2, dtype=np.int64)
            arrays[f"r{resolution}_count"] = np.asarray(count, dtype=np.float32)
        with open(out, "wb") as handle:
            np.savez(handle, **arrays)
    return 0


def main(argv):
    args = parse(argv)
    if args.inputFormat == "hic":
        return hic_input(args, argv)
    if args.outputFormat == "hic":
        return hic_output(args)
    return hicConvertFormat.main(argv)


if __name__ == "__main__":
    import reference_tree  # the reference tree, not the package (reference_tree.py)
    reference_tree.check_imported()
    sys.exit(main(sys.argv[1:]))
