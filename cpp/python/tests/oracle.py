"""Reference values for the hicx_matrix tests, from the reference readers.

The tests run this script in a separate interpreter (support.oracle_python), so
the reference packages and the module under test never share a process or an
HDF5 library.

    python oracle.py {cooler|hicmatrix|hicstraw} queries.json out.npz

Exit code 3 means the reference package cannot be imported.
"""

import json
import sys

import numpy as np


def region_bins(start, end, resolution):
    """The bins overlapping [start, end), as cooler counts them."""
    return start // resolution, max(start // resolution, -(-end // resolution))


def cooler_query(cooler, query):
    clr = cooler.Cooler(query["uri"])
    if query["what"] == "chromosomes":
        return {"names": np.array(clr.chromnames, dtype=str),
                "lengths": clr.chromsizes.values.astype(np.int64)}
    matrix = clr.matrix(balance=query["balance"]).fetch(query["region1"], query["region2"])
    return {"matrix": np.asarray(matrix, dtype=np.float64)}


def hicmatrix_query(hicmatrix, cache, query):
    path = query["path"]
    if path not in cache:
        cache[path] = hicmatrix.hiCMatrix(path)
    hic = cache[path]
    sizes = hic.get_chromosome_sizes()
    if query["what"] == "chromosomes":
        return {"names": np.array(list(sizes.keys()), dtype=str),
                "lengths": np.array(list(sizes.values()), dtype=np.int64)}
    ranges = []
    for chrom, start, end in (query["region1"], query["region2"]):
        start = 0 if start is None else start
        end = sizes[chrom] if end is None else end
        # getRegionBinRange returns the bins holding the two positions, both
        # included; end - 1 is the last base of the half open region.
        first, last = hic.getRegionBinRange(chrom, start, end - 1)
        ranges.append((first, last + 1))
    (a, b), (c, d) = ranges
    return {"matrix": hic.matrix[a:b, c:d].toarray().astype(np.float64)}


def hicstraw_query(hicstraw, cache, query):
    path = query["path"]
    if path not in cache:
        cache[path] = hicstraw.HiCFile(path)
    hic = cache[path]
    chromosomes = hic.getChromosomes()
    if query["what"] == "chromosomes":
        kept = [c for c in chromosomes if c.name.upper() != "ALL"]
        return {"names": np.array([c.name for c in kept], dtype=str),
                "lengths": np.array([c.length for c in kept], dtype=np.int64),
                "resolutions": np.array(hic.getResolutions(), dtype=np.int64)}
    index = {c.name: i for i, c in enumerate(chromosomes)}
    length = {c.name: c.length for c in chromosomes}
    res = query["resolution"]
    spans = []
    for chrom, start, end in (query["region1"], query["region2"]):
        start = 0 if start is None else start
        end = length[chrom] if end is None else end
        spans.append((chrom,) + region_bins(start, end, res))
    (chrom1, a, b), (chrom2, c, d) = spans
    # hicstraw's x axis is the chromosome with the lower index.
    swap = index[chrom1] > index[chrom2]
    (cx, x0, x1), (cy, y0, y1) = (spans[1], spans[0]) if swap else (spans[0], spans[1])
    zoom = hic.getMatrixZoomData(cx, cy, query["matrix_type"], query["norm"], "BP", res)
    matrix = np.asarray(zoom.getRecordsAsMatrix(x0 * res, x1 * res - 1, y0 * res, y1 * res - 1),
                        dtype=np.float64)
    if matrix.shape == (1, 1) and (x1 - x0, y1 - y0) != (1, 1):
        # getRecordsAsMatrix returns a single 0 when no record is inside.
        matrix = np.zeros((x1 - x0, y1 - y0))
    if swap:
        matrix = matrix.T
    return {"matrix": matrix}


def main():
    kind, query_file, out_file = sys.argv[1:4]
    with open(query_file) as handle:
        queries = json.load(handle)
    try:
        if kind == "cooler":
            import cooler as reader
        elif kind == "hicmatrix":
            from hicmatrix import HiCMatrix as reader
        elif kind == "hicstraw":
            import hicstraw as reader
        else:
            raise SystemExit(f"unknown oracle kind {kind!r}")
    except ImportError as exc:
        print(f"{kind} is not importable: {exc}", file=sys.stderr)
        sys.exit(3)
    cache = {}
    results = {}
    for query in queries:
        if kind == "cooler":
            values = cooler_query(reader, query)
        elif kind == "hicmatrix":
            values = hicmatrix_query(reader, cache, query)
        else:
            values = hicstraw_query(reader, cache, query)
        for key, value in values.items():
            results[f"{query['id']}__{key}"] = value
    np.savez(out_file, **results)


if __name__ == "__main__":
    main()
