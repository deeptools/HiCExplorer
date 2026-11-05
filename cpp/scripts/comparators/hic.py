"""Comparator for the .hic files the C++ hicConvertFormat writes.

No Python program writes the Juicer .hic format, so the Python half of such a
case (cpp/scripts/py_hic_reference.py) writes, at the output path, an .npz
file that holds what the .hic file must contain: the version, the chromosomes,
the normalizations asked for, and the pixels of every resolution. This
comparator reads the .hic file with its own small reader (header, footer
index, the blocks of versions 8 and 9) and checks

- the version, the chromosome names and lengths after "All", and the base
  pair resolutions;
- at every resolution, the pixels (chromosome pair, bins, count) against the
  expected ones: bit-identical counts for E0 to E2, the class's float
  tolerance above that;
- that every normalization asked for has a vector for every chromosome with an
  intra-chromosomal matrix at every resolution, with as many values as Juicer
  tools gives it: the bin count of the matrix's grid axis, block bin count
  times block column count.

The vectors' values are validated against Juicer tools pre and addNorm in
hicfilecpp's own harness, not here.
"""
from __future__ import annotations

import struct
import zlib

import numpy as np

from .base import Result, fail, values_agree


class _Reader:
    def __init__(self, data, offset=0):
        self.data = data
        self.pos = offset

    def get(self, fmt):
        value = struct.unpack_from("<" + fmt, self.data, self.pos)[0]
        self.pos += struct.calcsize("<" + fmt)
        return value

    def cstr(self):
        end = self.data.index(b"\0", self.pos)
        text = self.data[self.pos:end].decode("utf-8")
        self.pos = end + 1
        return text


def read_hic(path):
    with open(path, "rb") as handle:
        data = handle.read()
    r = _Reader(data)
    if r.cstr()[:3] != "HIC":
        raise ValueError("not a .hic file")
    version = r.get("i")
    master = r.get("q")
    genome = r.cstr()
    if version > 8:
        r.get("q")
        r.get("q")
    for _ in range(r.get("i")):
        r.cstr()
        r.cstr()
    chromosomes = []
    for _ in range(r.get("i")):
        name = r.cstr()
        chromosomes.append((name, r.get("q") if version > 8 else r.get("i")))
    resolutions = [r.get("i") for _ in range(r.get("i"))]
    f = _Reader(data, master)
    f.get("q" if version > 8 else "i")
    matrices = {}
    for _ in range(f.get("i")):
        key = f.cstr()
        matrices[key] = f.get("q")
        f.get("i")
    value = 4 if version > 8 else 8
    for normalized in (False, True):
        for _ in range(f.get("i")):
            if normalized:
                f.cstr()
            f.cstr()
            f.get("i")
            n = f.get("q") if version > 8 else f.get("i")
            f.pos += n * value
            factors = f.get("i")
            f.pos += factors * (4 + value)
    norm_vectors = {}
    for _ in range(f.get("i")):
        norm = f.cstr()
        chrom = f.get("i")
        unit = f.cstr()
        resolution = f.get("i")
        position = f.get("q")
        f.get("q" if version > 8 else "i")
        count = struct.unpack_from("<q" if version > 8 else "<i", data, position)[0]
        norm_vectors[(norm, chrom, unit, resolution)] = count
    zooms = {}
    layouts = {}
    for key, position in matrices.items():
        m = _Reader(data, position)
        m.get("i")
        m.get("i")
        for _ in range(m.get("i")):
            unit = m.cstr()
            m.pos += 20
            bin_size = m.get("i")
            block_bin_count = m.get("i")
            block_column_count = m.get("i")
            blocks = [struct.unpack_from("<iqi", data, m.pos + 16 * b) for b in range(m.get("i"))]
            m.pos += 16 * len(blocks)
            zooms[(key, unit, bin_size)] = blocks
            layouts[(key, unit, bin_size)] = block_bin_count * block_column_count
    return {"data": data, "version": version, "genome": genome, "chromosomes": chromosomes,
            "resolutions": resolutions, "zooms": zooms, "layouts": layouts, "norm_vectors": norm_vectors}


def _decode(raw, version):
    data = zlib.decompress(raw)
    _count, x_offset, y_offset = struct.unpack_from("<iii", data, 0)
    pos = 12
    use_short = data[pos] == 0
    pos += 1
    short_x = short_y = True
    if version > 8:
        short_x, short_y = data[pos] == 0, data[pos + 1] == 0
        pos += 2
    kind = data[pos]
    pos += 1
    xs, ys, cs = [], [], []
    if kind == 1:
        row_fmt = "<h" if short_y else "<i"
        col_fmt = "<h" if short_x else "<i"
        row_size, col_size = struct.calcsize(row_fmt), struct.calcsize(col_fmt)
        record = np.dtype([("x", "<i2" if short_x else "<i4"), ("c", "<i2" if use_short else "<f4")])
        rows = struct.unpack_from(row_fmt, data, pos)[0]
        pos += row_size
        for _ in range(rows):
            y = struct.unpack_from(row_fmt, data, pos)[0]
            pos += row_size
            n = struct.unpack_from(col_fmt, data, pos)[0]
            pos += col_size
            values = np.frombuffer(data, dtype=record, count=n, offset=pos)
            pos += n * record.itemsize
            xs.append(values["x"].astype(np.int64) + x_offset)
            ys.append(np.full(n, y + y_offset, dtype=np.int64))
            cs.append(values["c"].astype(np.float32))
    elif kind == 2:
        n, width = struct.unpack_from("<ih", data, pos)
        pos += 6
        values = np.frombuffer(data, dtype="<i2" if use_short else "<f4", count=n, offset=pos)
        keep = values != -32768 if use_short else ~np.isnan(values)
        index = np.flatnonzero(keep)
        xs.append(x_offset + index % width)
        ys.append(y_offset + index // width)
        cs.append(values[keep].astype(np.float32))
    if not xs:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.float32)
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(cs)


def compare(path_a, path_b, cls, opts=None):
    try:
        expected = np.load(path_a, allow_pickle=False)
    except Exception as error:  # pylint: disable=W0718
        return fail(f"cannot read the expected pixels {path_a}: {error}")
    try:
        hic = read_hic(path_b)
    except Exception as error:  # pylint: disable=W0718
        return fail(f"cannot read {path_b} as a .hic file: {error}")
    diffs = []
    metrics = {"version": hic["version"]}
    if hic["version"] != int(expected["version"]):
        diffs.append(f"version {hic['version']} vs expected {int(expected['version'])}")
    names = [str(name) for name in expected["chrom_names"]]
    lengths = [int(length) for length in expected["chrom_lengths"]]
    header = hic["chromosomes"]
    if not header or header[0][0].lower() != "all" or \
            [h[0] for h in header[1:]] != names or [int(h[1]) for h in header[1:]] != lengths:
        diffs.append(f"chromosomes {header[:4]}... vs expected {list(zip(names, lengths))[:3]}...")
        return Result(False, None, metrics, diffs)
    resolutions = [int(r) for r in expected["resolutions"]]
    if hic["resolutions"] != resolutions:
        diffs.append(f"resolutions {hic['resolutions']} vs expected {resolutions}")
        return Result(False, None, metrics, diffs)

    data = hic["data"]
    for resolution in resolutions:
        chr1, chr2, bin1, bin2, count = [], [], [], [], []
        for (key, unit, bin_size), blocks in hic["zooms"].items():
            if unit != "BP" or bin_size != resolution or key == "0_0":
                continue
            c1, c2 = (int(v) - 1 for v in key.split("_"))
            for _number, position, size in blocks:
                x, y, c = _decode(data[position:position + size], hic["version"])
                chr1.append(np.full(len(x), c1))
                chr2.append(np.full(len(x), c2))
                bin1.append(x)
                bin2.append(y)
                count.append(c)
        got = [np.concatenate(v) if v else np.zeros(0) for v in (chr1, chr2, bin1, bin2)]
        got_count = np.concatenate(count) if count else np.zeros(0, np.float32)
        want = [expected[f"r{resolution}_{k}"].astype(np.int64) for k in ("chr1", "chr2", "bin1", "bin2")]
        want_count = expected[f"r{resolution}_count"].astype(np.float32)
        order_got = np.lexsort((got[3], got[2], got[1], got[0]))
        order_want = np.lexsort((want[3], want[2], want[1], want[0]))
        metrics[f"pixels_{resolution}"] = int(len(got_count))
        if len(got_count) != len(want_count) or not all(
                np.array_equal(g.astype(np.int64)[order_got], w[order_want]) for g, w in zip(got, want)):
            diffs.append(f"resolution {resolution}: {len(got_count)} pixels vs expected "
                         f"{len(want_count)}, or different coordinates")
            continue
        a = got_count[order_got]
        b = want_count[order_want]
        if cls in ("E0", "E1", "E2"):
            if a.tobytes() != b.tobytes():
                diffs.append(f"resolution {resolution}: {int(np.count_nonzero(a != b))} counts differ")
        else:
            bad = [i for i in np.flatnonzero(a != b) if not values_agree(float(b[i]), float(a[i]), cls)]
            if bad:
                diffs.append(f"resolution {resolution}: {len(bad)} counts differ beyond {cls}")

    wanted = [str(n) for n in expected["normalizations"]]
    missing = 0
    for resolution in resolutions:
        for index, length in enumerate(lengths, start=1):
            zoom = (f"{index}_{index}", "BP", resolution)
            if zoom not in hic["zooms"]:
                continue
            expected_values = hic["layouts"][zoom]
            for norm in wanted:
                count_values = hic["norm_vectors"].get((norm, index, "BP", resolution))
                if count_values != expected_values:
                    missing += 1
                    if missing <= 5:
                        diffs.append(f"{norm} vector of {names[index - 1]} at {resolution}: "
                                     f"{count_values} values, expected {expected_values}")
    metrics["normalization_vectors"] = len(hic["norm_vectors"])
    if diffs:
        return Result(False, None, metrics, diffs)
    return Result(True, cls, metrics)
