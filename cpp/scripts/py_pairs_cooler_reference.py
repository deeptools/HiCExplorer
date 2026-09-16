#!/usr/bin/env python3
"""cooler cload pairs as the reference of hicBuildMatrix --pairsFile
(cpp/PLAN.md tier 9, item 9.2, class EX).

    py_pairs_cooler_reference.py <hicBuildMatrix --pairsFile arguments>

Run it as equiv.py runs every py_script: with the reference venv interpreter
(cooler 0.10.2) and PYTHONPATH at the repository root.

cooler bins every pair of a file. hicBuildMatrix --pairsFile first applies the
filters that its pairs route keeps (build_matrix_pairs_impl.hpp, P1 to P4), so
this script applies them first, with its own implementation, and hands the
pairs that remain to cooler:

- a pair with chromosome '!' on a side, or a pair_type with a letter other
  than U and R, is dropped (unmapped, not unique); so is pair_type DD;
- with columns mapq1 and mapq2, a pair with either below --minMappingQuality
  (default 15) is dropped;
- unless --skipDuplicationCheck is given, a pair whose two (chromosome,
  position) ends repeat those of an earlier pair, in either order, is dropped.

When the file has neither pair_type nor mapq columns and --skipDuplicationCheck
is given, nothing needs filtering and cooler reads the file itself. Otherwise
the pairs that remain are written with the file's header to a temporary file.
cooler then gets the column numbers of the #columns line (or of the reserved
layout), one-based positions, and the #chromsize lines (or
--chromosomeSizes) in their order with the first --binSize as the bins. The
output is the .cool file --outFileName names; the QC folder is created empty,
since cooler writes no QC.

Supported: one or more --binSize values of which the first is used, a .cool
output name, no --region and no restriction fragment bins. Anything else exits
2 with the reason.
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RESERVED = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]


def parse(argv):
    parser = argparse.ArgumentParser(prog="py_pairs_cooler_reference.py", add_help=False)
    parser.add_argument("--pairsFile", required=True)
    parser.add_argument("--outFileName", "-o", required=True)
    parser.add_argument("--QCfolder", required=True)
    parser.add_argument("--binSize", "-bs", type=int, nargs="+")
    parser.add_argument("--chromosomeSizes", "-cs")
    parser.add_argument("--minMappingQuality", type=int, default=15)
    parser.add_argument("--skipDuplicationCheck", action="store_true")
    parser.add_argument("--genomeAssembly", "-ga")
    parser.add_argument("--threads")
    parser.add_argument("--inputBufferSize")
    return parser.parse_args(argv)


def open_text(path):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def read_header(path):
    header, columns, sizes, assembly = [], None, [], None
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            header.append(line)
            if line.startswith("#columns:"):
                columns = line.split(":", 1)[1].split()
            elif line.startswith("#chromsize:"):
                name, size = line.split(":", 1)[1].split()
                sizes.append((name, int(size)))
            elif line.startswith("#genome_assembly:"):
                assembly = line.split(":", 1)[1].strip()
    columns = [{"chr1": "chrom1", "chr2": "chrom2"}.get(name, name) for name in (columns or RESERVED)]
    return header, columns, sizes, assembly


def main(argv):
    args = parse(argv)
    if not args.binSize:
        sys.exit("the cooler reference covers --binSize only")
    if not args.outFileName.endswith(".cool"):
        sys.exit("the cooler reference writes .cool output only")
    header, columns, sizes, assembly = read_header(args.pairsFile)
    index = {name: k for k, name in enumerate(columns)}
    work = Path(tempfile.mkdtemp(prefix="pairs_cooler_reference_"))
    try:
        chromsizes = work / "chrom.sizes"
        if args.chromosomeSizes:
            shutil.copyfile(args.chromosomeSizes, chromsizes)
        else:
            chromsizes.write_text("".join(f"{name}\t{size}\n" for name, size in sizes))

        has_type = "pair_type" in index
        has_mapq = "mapq1" in index and "mapq2" in index
        pairs = args.pairsFile
        if has_type or has_mapq or not args.skipDuplicationCheck:
            pairs = str(work / "selected.pairs")
            c1, p1, c2, p2 = (index[k] for k in ("chrom1", "pos1", "chrom2", "pos2"))
            seen = set()
            kept = dropped = 0
            with open_text(args.pairsFile) as source, open(pairs, "w") as target:
                target.writelines(header)
                for line in source:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip("\r\n").split("\t")
                    keep = fields[c1] != "!" and fields[c2] != "!"
                    if keep and has_type:
                        kind = fields[index["pair_type"]]
                        keep = kind != "DD" and (kind == "." or all(letter in "UR" for letter in kind))
                    if keep and has_mapq:
                        keep = int(fields[index["mapq1"]]) >= args.minMappingQuality and \
                            int(fields[index["mapq2"]]) >= args.minMappingQuality
                    if keep and not args.skipDuplicationCheck:
                        end1 = (fields[c1], int(fields[p1]))
                        end2 = (fields[c2], int(fields[p2]))
                        key = (end1, end2) if end1 <= end2 else (end2, end1)
                        keep = key not in seen
                        seen.add(key)
                    if keep:
                        target.write(line)
                        kept += 1
                    else:
                        dropped += 1
            print(f"selected {kept} pairs, dropped {dropped}", file=sys.stderr)

        command = [sys.executable, "-m", "cooler", "cload", "pairs",
                   "-c1", str(index["chrom1"] + 1), "-p1", str(index["pos1"] + 1),
                   "-c2", str(index["chrom2"] + 1), "-p2", str(index["pos2"] + 1)]
        if args.genomeAssembly or assembly:
            command += ["--assembly", args.genomeAssembly or assembly]
        command += [f"{chromsizes}:{args.binSize[0]}", pairs, args.outFileName]
        print(" ".join(command), file=sys.stderr)
        if os.path.exists(args.outFileName):
            os.remove(args.outFileName)
        subprocess.run(command, check=True)
        os.makedirs(args.QCfolder, exist_ok=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
