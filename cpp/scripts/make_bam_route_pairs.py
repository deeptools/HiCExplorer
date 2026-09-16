#!/usr/bin/env python3
"""Writes the .pairs inputs of the hicBuildMatrix --pairsFile cases from the
committed test BAM files (cpp/PLAN.md tier 9, item 9.2, class EX).

    cd <repository root>
    PYTHONPATH=. $VENV/bin/python cpp/scripts/make_bam_route_pairs.py \\
        --mode {valid,all} --order {sorted,input} [--columns {full,reserved}] \\
        [--chrom-names {chrom,chr}] --out hicexplorer/test/test_data/hicBuildMatrix/pairs/NAME \\
        -- <hicBuildMatrix arguments>

The hicBuildMatrix arguments write {data} for hicexplorer/test/test_data, and
leave out --outFileName, --QCfolder, --outBam and --threads. An --out name
ending in .gz is bgzip compressed. $VENV is the reference venv of
cpp/AGENTS_CONTRACT.md. Every read ID is written as '.', which the
specification allows, to keep the files small.

--mode valid
    The read pairs the Python hicBuildMatrix keeps. The script runs it with
    --outBam and --threads 2 (a single worker, so the records keep the input
    order) and converts the pairs of the written BAM file. The arguments go to
    the header line '#hicexplorer_bam_route:', which
    cpp/scripts/py_pairs_bam_route_reference.py reads to run the same BAM route
    as the reference of a case.

    Position. hicBuildMatrix places a read at its middle, pos + int(qlen / 2)
    (zero-based, buildMatrixMethods.py:661), in the bin whose closed interval
    [begin, end] holds it, found by a binary search (:662-688). hicBuildMatrix
    --pairsFile, like cooler, bins position - 1 into the half-open interval
    [start, end). The two agree except where the middle equals the end of one
    bin and the start of the next, where the closed search may return either
    bin. The script finds the bin the Python search returns, with a copy of
    that search over the interval array createMatrix builds (:946-962), and
    writes middle + 1 where the half-open bin is the same and middle, one base
    to the left, where the search returned the left bin. The number of moved
    read ends goes to '#hicexplorer_moved_ends:'. Any other disagreement stops
    the script.

--mode all
    Every primary read pair of the two BAM files before any filter, in
    pairtools' conventions: an unmapped side is pair type letter N with
    chromosome '!', position 0 and strand '-'; a side with mapping quality 0 is
    M with chromosome '!', position 0 and its strand; any other side is U at its
    read middle + 1. mapq1 and mapq2 carry the mapping qualities. Only
    --samFiles and --chromosomeSizes of the arguments are used.

--order sorted
    Each pair is flipped so that its first side has the smaller (chromosome
    index in the chromosome size list, position), '!' first, and the file is
    sorted by chromosome name 1, chromosome name 2, position 1, position 2,
    as pairtools sort writes it; the header says '#sorted: chr1-chr2-pos1-pos2'
    and '#shape: upper triangle'. --order input keeps the BAM order and the
    first BAM file on side 1, with '#sorted: none'.

--columns full writes readID chrom1 pos1 chrom2 pos2 strand1 strand2
    pair_type mapq1 mapq2; --columns reserved writes the seven reserved columns
    and no #columns line, which the specification makes the default layout.
--chrom-names chr names the chromosome columns chr1 and chr2, as ENCODE's
    files do, in #columns and #sorted.
"""
from __future__ import annotations

import argparse
import bisect
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "hicexplorer" / "test" / "test_data"


def expand(argv):
    return [argument.replace("{data}", str(DATA)) for argument in argv]


def chrom_sizes_of(options, bam_path):
    if options.chromosomeSizes:
        sizes = OrderedDict()
        with open(options.chromosomeSizes) as handle:
            for line in handle:
                line = line.strip()
                if line:
                    fields = line.split("\t")
                    sizes[fields[0]] = int(fields[1])
        return list(sizes.items())
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return list(zip(bam.references, bam.lengths))


def python_search_array(options, chrom_sizes):
    """The shared interval array and index of createMatrix, :925-962."""
    from hicexplorer.lib import buildMatrixMethods as methods

    rf_interval = []
    for path in options.restrictionCutFile or []:
        with open(path) as handle:
            rf_interval.extend(methods.bed2interval_list(handle, chrom_sizes, options.region))
    max_insert = options.maxDistance if options.maxDistance else options.maxLibraryInsertSize
    if options.binSize:
        bins = methods.get_bins(options.binSize[0], chrom_sizes, options.region)
    else:
        bins = methods.get_rf_bins(rf_interval, min_distance=options.minDistance,
                                   max_distance=max_insert)
    tree = methods.intervalListToIntervalTree(bins)
    shared, index, end = [], {}, -1
    for seq in tree:
        start = end + 1
        items = sorted((interval.begin, interval.end, interval.data) for interval in tree[seq])
        end = start + len(tree[seq]) - 1
        index[seq] = (start, end)
        # C_Interval holds c_uint fields.
        shared.extend((b & 0xFFFFFFFF, e & 0xFFFFFFFF, d) for b, e, d in items)
    return shared, index


def python_search(shared, index, chrom, read_middle):
    """buildMatrixMethods.py:662-688: the flat index of the closed interval found."""
    if chrom not in index:
        return None
    start, end = index[chrom]
    middle = int((start + end) / 2)
    while not start > end:
        begin, stop, _ = shared[middle]
        if begin <= read_middle <= stop:
            return middle
        if begin > read_middle:
            end = middle - 1
        else:
            start = middle + 1
        middle = int((start + end) / 2)
    return None


def half_open(shared, index, chrom, position):
    start, end = index[chrom]
    begins = [shared[k][0] for k in range(start, end + 1)]
    k = bisect.bisect_right(begins, position) - 1
    if k < 0 or position >= shared[start + k][1]:
        return None
    return start + k


def valid_pairs(argv, workdir):
    from hicexplorer import hicBuildMatrix

    expanded = expand(argv)
    for forbidden in ("--outFileName", "-o", "--QCfolder", "--outBam", "-b", "--threads"):
        if forbidden in expanded:
            sys.exit(f"leave {forbidden} out of the hicBuildMatrix arguments")
    valid_bam = workdir / "valid.bam"
    run = expanded + ["--outFileName", str(workdir / "matrix.h5"), "--QCfolder",
                      str(workdir / "qc"), "--outBam", str(valid_bam), "--threads", "2"]
    hicBuildMatrix.main(run)

    # Parsed again for the bins, with placeholder outputs: argparse's
    # FileType('w') opens what it names, and --outBam would truncate the BAM.
    options = hicBuildMatrix.parse_arguments().parse_args(
        expanded + ["--outFileName", str(workdir / "parse.h5"), "--QCfolder", str(workdir / "qc")])
    options.outFileName.close()
    for handle in list(options.samFiles) + list(options.restrictionCutFile or []):
        handle.close()
    options.samFiles = [handle.name for handle in options.samFiles]
    options.restrictionCutFile = [handle.name for handle in options.restrictionCutFile or []]
    if options.chromosomeSizes:
        options.chromosomeSizes.close()
        options.chromosomeSizes = options.chromosomeSizes.name
    chrom_sizes = chrom_sizes_of(options, options.samFiles[0])
    shared, index = python_search_array(options, chrom_sizes)

    pairs, moved = [], 0
    with pysam.AlignmentFile(str(valid_bam), "rb") as bam:
        records = iter(bam)
        for mate1 in records:
            mate2 = next(records)
            if mate1.query_name != mate2.query_name:
                sys.exit(f"--outBam records out of pairs: {mate1.query_name} {mate2.query_name}")
            sides = []
            for mate in (mate1, mate2):
                chrom = bam.get_reference_name(mate.reference_id)
                middle = mate.pos + int(mate.qlen / 2)
                found = python_search(shared, index, chrom, middle)
                if found is None:
                    sys.exit(f"a kept read of {mate.query_name} finds no bin")
                position = middle
                if half_open(shared, index, chrom, middle) != found:
                    if shared[found][1] == middle and \
                            half_open(shared, index, chrom, middle - 1) == found:
                        position = middle - 1
                        moved += 1
                    else:
                        sys.exit(f"{mate.query_name}: bin {found} of the closed search cannot "
                                 f"be reached half-open from {middle}")
                sides.append((chrom, position + 1, "-" if mate.is_reverse else "+", "U",
                              mate.mapq))
            pairs.append(sides)
    return pairs, chrom_sizes, moved


def all_pairs(argv):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--samFiles", "-s", nargs=2, required=True)
    parser.add_argument("--chromosomeSizes", "-cs")
    options, _ = parser.parse_known_args(expand(argv))
    chrom_sizes = chrom_sizes_of(options, options.samFiles[0])
    pairs = []
    with pysam.AlignmentFile(options.samFiles[0], "rb") as bam1, \
            pysam.AlignmentFile(options.samFiles[1], "rb") as bam2:
        records1, records2 = iter(bam1), iter(bam2)
        for mate1 in records1:
            mate2 = next(records2)
            # readBamFiles, buildMatrixMethods.py:481-492.
            while mate1.flag & 256 == 256:
                mate1 = next(records1)
            while mate2.flag & 256 == 256:
                mate2 = next(records2)
            if mate1.query_name != mate2.query_name:
                sys.exit(f"the BAM files are out of step: {mate1.query_name} {mate2.query_name}")
            sides = []
            for bam, mate in ((bam1, mate1), (bam2, mate2)):
                if mate.has_tag("SA"):
                    sys.exit(f"{mate.query_name} has a supplementary alignment, which this "
                             "conversion does not resolve")
                strand = "-" if mate.is_reverse else "+"
                if mate.is_unmapped:
                    sides.append(("!", 0, "-", "N", 0))
                elif mate.mapq == 0:
                    sides.append(("!", 0, strand, "M", 0))
                else:
                    sides.append((bam.get_reference_name(mate.reference_id),
                                  mate.pos + int(mate.qlen / 2) + 1, strand, "U", mate.mapq))
            pairs.append(sides)
    return pairs, chrom_sizes


def write(pairs, chrom_sizes, options, extra_header):
    rank = {name: k for k, (name, _) in enumerate(chrom_sizes)}
    rank["!"] = -1
    if options.order == "sorted":
        flipped = []
        for side1, side2 in pairs:
            if (rank.get(side2[0], len(rank)), side2[1]) < (rank.get(side1[0], len(rank)), side1[1]):
                side1, side2 = side2, side1
            flipped.append((side1, side2))
        pairs = sorted(flipped, key=lambda p: (p[0][0], p[1][0], p[0][1], p[1][1]))
    c1, c2 = ("chr1", "chr2") if options.chrom_names == "chr" else ("chrom1", "chrom2")
    lines = ["## pairs format v1.0"]
    if options.order == "sorted":
        lines += [f"#sorted: {c1}-{c2}-pos1-pos2", "#shape: upper triangle"]
    else:
        lines += ["#sorted: none"]
    lines += [f"#chromsize: {name} {size}" for name, size in chrom_sizes]
    lines += extra_header
    if options.columns == "full":
        lines.append(f"#columns: readID {c1} pos1 {c2} pos2 strand1 strand2 pair_type mapq1 mapq2")
    for side1, side2 in pairs:
        fields = [".", side1[0], str(side1[1]), side2[0], str(side2[1]), side1[2], side2[2]]
        if options.columns == "full":
            fields += [side1[3] + side2[3], str(side1[4]), str(side2[4])]
        lines.append("\t".join(fields))
    text = "\n".join(lines) + "\n"
    out = Path(options.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".gz":
        plain = out.with_suffix("")
        plain.write_text(text)
        pysam.tabix_compress(str(plain), str(out), force=True)
        plain.unlink()
    else:
        out.write_text(text)
    return len(pairs)


def main():
    if "--" not in sys.argv:
        sys.exit(__doc__)
    split = sys.argv.index("--")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["valid", "all"], required=True)
    parser.add_argument("--order", choices=["sorted", "input"], required=True)
    parser.add_argument("--columns", choices=["full", "reserved"], default="full")
    parser.add_argument("--chrom-names", choices=["chrom", "chr"], default="chrom")
    parser.add_argument("--out", required=True)
    options = parser.parse_args(sys.argv[1:split])
    argv = sys.argv[split + 1:]
    recorded = shlex.join(argv)
    workdir = Path(tempfile.mkdtemp(prefix="bam_route_pairs_"))
    try:
        if options.mode == "valid":
            pairs, chrom_sizes, moved = valid_pairs(argv, workdir)
            extra = [f"#hicexplorer_bam_route: {recorded}",
                     "#hicexplorer_position: the read middle of hicBuildMatrix plus one "
                     "(cpp/scripts/make_bam_route_pairs.py)",
                     f"#hicexplorer_moved_ends: {moved}"]
        else:
            pairs, chrom_sizes = all_pairs(argv)
            moved = 0
            extra = [f"#hicexplorer_bam_source: {recorded}",
                     "#hicexplorer_position: the read middle of hicBuildMatrix plus one "
                     "(cpp/scripts/make_bam_route_pairs.py)"]
        count = write(pairs, chrom_sizes, options, extra)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    print(f"{options.out}: {count} pairs, {moved} read ends moved one base to the left")


if __name__ == "__main__":
    main()
