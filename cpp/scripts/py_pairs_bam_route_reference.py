#!/usr/bin/env python3
"""The Python BAM route as the reference of hicBuildMatrix --pairsFile
(cpp/PLAN.md tier 9, item 9.2, class EX).

    py_pairs_bam_route_reference.py <hicBuildMatrix --pairsFile arguments>

Run it as equiv.py runs every py_script: with the reference venv interpreter
and PYTHONPATH at the repository root.

The .pairs file must be one that cpp/scripts/make_bam_route_pairs.py --mode
valid wrote: it carries exactly the read pairs the Python hicBuildMatrix keeps
from two BAM files, and its header line '#hicexplorer_bam_route:' holds the
arguments of that BAM route. The script runs the Python hicBuildMatrix on the
BAM files with those arguments and the case's --outFileName, --QCfolder,
--threads and --genomeAssembly.

Both runs must build the same bins and keep the same pairs, so the script
first checks the case against the recorded arguments and exits 2 on a
mismatch:
- the first --binSize equals the recorded one (further sizes of the case, for
  an mcool, are passed on), or neither has --binSize and --restrictionCutFile,
  --minDistance, --maxLibraryInsertSize and --maxDistance are equal;
- --region, --chromosomeSizes and --minMappingQuality are equal.
--skipDuplicationCheck of the case is not passed on: the file holds the pairs
of the recorded run, whose duplicate check has already removed duplicates.
--inputBufferSize only sizes chunks and is not passed on either.
"""
from __future__ import annotations

import argparse
import gzip
import shlex
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data"


def options_parser():
    parser = argparse.ArgumentParser(prog="py_pairs_bam_route_reference.py", add_help=False)
    parser.add_argument("--binSize", "-bs", type=int, nargs="+")
    parser.add_argument("--restrictionCutFile", "-rs", nargs="+")
    parser.add_argument("--minDistance", type=int, default=300)
    parser.add_argument("--maxLibraryInsertSize", type=int, default=1000)
    parser.add_argument("--maxDistance", type=int)
    parser.add_argument("--region", "-r")
    parser.add_argument("--chromosomeSizes", "-cs")
    parser.add_argument("--minMappingQuality", type=int, default=15)
    return parser


def recorded_route(path):
    opener = gzip.open if open(path, "rb").read(2) == b"\x1f\x8b" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("#hicexplorer_bam_route:"):
                return shlex.split(line.split(":", 1)[1])
    sys.exit(f"{path} has no '#hicexplorer_bam_route:' header line; it was not written by "
             "make_bam_route_pairs.py --mode valid")


def main(argv):
    case_parser = options_parser()
    case_parser.add_argument("--pairsFile", required=True)
    case_parser.add_argument("--outFileName", "-o", required=True)
    case_parser.add_argument("--QCfolder", required=True)
    case_parser.add_argument("--threads", default="4")
    case_parser.add_argument("--genomeAssembly", "-ga")
    case_parser.add_argument("--skipDuplicationCheck", action="store_true")
    case_parser.add_argument("--inputBufferSize")
    case = case_parser.parse_args(argv)

    tokens = [token.replace("{data}", str(DATA)) for token in recorded_route(case.pairsFile)]
    recorded, rest = options_parser().parse_known_args(tokens)

    def mismatch(what, left, right):
        sys.exit(f"the case's {what} {left!r} differs from the recorded BAM route's {right!r}")

    if bool(case.binSize) != bool(recorded.binSize):
        mismatch("--binSize", case.binSize, recorded.binSize)
    if case.binSize and case.binSize[0] != recorded.binSize[0]:
        mismatch("first --binSize", case.binSize[0], recorded.binSize[0])
    if not case.binSize:
        for name in ("restrictionCutFile", "minDistance", "maxLibraryInsertSize", "maxDistance"):
            if getattr(case, name) != getattr(recorded, name):
                mismatch("--" + name, getattr(case, name), getattr(recorded, name))
    for name in ("region", "chromosomeSizes", "minMappingQuality"):
        if getattr(case, name) != getattr(recorded, name):
            mismatch("--" + name, getattr(case, name), getattr(recorded, name))

    run = list(rest)
    if recorded.restrictionCutFile:
        run += ["--restrictionCutFile"] + recorded.restrictionCutFile
    if recorded.binSize:
        run += ["--binSize"] + [str(size) for size in case.binSize]
    run += ["--minDistance", str(recorded.minDistance),
            "--maxLibraryInsertSize", str(recorded.maxLibraryInsertSize),
            "--minMappingQuality", str(recorded.minMappingQuality)]
    for name in ("maxDistance", "region", "chromosomeSizes"):
        if getattr(recorded, name) is not None:
            run += ["--" + name, str(getattr(recorded, name))]
    run += ["--outFileName", case.outFileName, "--QCfolder", case.QCfolder,
            "--threads", case.threads]
    if case.genomeAssembly:
        run += ["--genomeAssembly", case.genomeAssembly]
    print("hicBuildMatrix " + shlex.join(run), file=sys.stderr)

    from hicexplorer import hicBuildMatrix
    hicBuildMatrix.main(run)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
