#!/usr/bin/env bash
# The tiled, reduced handoff of hicPlotMatrix (--plotData shows the files the
# drawing layer reads). Checked on small_test_matrix_50kb_res.cool (3383 bins,
# integer counts) and on its h5 twin, which has NaN bins on the whole-matrix
# path:
#
#  1. full resolution in one tile and in tiles of 700 bins: the tiles, put
#     together at their block ranges, are the one matrix (transposed, as
#     pcolormesh wants them), and the cell edges are the same;
#  2. the element type: int32 for counts, also after --log1p (+1 keeps them
#     integers); float32 where NaN bins must be NaN;
#  3. the reduction (--dpi 20 leaves 100 x 96 pixels): the blocks are the sums
#     of the full resolution cells (the total is the same, blocks are cut at
#     chromosome ends and cover the matrix edges exactly), NaN bins do not spread over the blocks that hold valid bins and --noDownsample switches the reduction off;
#  4. the colour limits of several tiles (vmin and vmax of all of them).
#
# Usage: hicPlotMatrix_tiles_cli.sh <hicPlotMatrix binary> <test_data directory>

set -u
tool=$1
data=$2
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
failures=0

expect() {
    local description=$1
    shift
    if "$@"; then
        echo "ok: $description"
    else
        echo "FAILED: $description"
        failures=$((failures + 1))
    fi
}

matrix="$data/small_test_matrix_50kb_res.cool"
matrix_nan="$data/small_test_matrix_50kb_res.h5"
cd "$work" || exit 1

run() {  # run NAME ARGS...: --plotData NAME.json, quiet
    local name=$1
    shift
    "$tool" -m "${source:-$matrix}" "$@" --plotData "$name.json" -out "$name.png" > "$name.out" 2> "$name.err"
}

run full1 --noDownsample --tileBins 4096 --log1p
expect "one tile at full resolution" [ $? -eq 0 ]
run full700 --noDownsample --tileBins 700 --log1p
expect "tiles of 700 bins" [ $? -eq 0 ]
source=$matrix_nan run nolog --noDownsample --tileBins 4096
expect "no log, one tile" [ $? -eq 0 ]
run reduced --dpi 20 --log1p
expect "reduced" [ $? -eq 0 ]
source=$matrix_nan run reduced_nolog --dpi 20
expect "reduced without log" [ $? -eq 0 ]
run reduced_off --dpi 20 --noDownsample --tileBins 4096 --log1p
expect "reduction switched off" [ $? -eq 0 ]
expect "the reduction is announced" grep -q "drawn as" reduced.err
expect "no reduction is announced with --noDownsample" bash -c '! grep -q "drawn as" reduced_off.err'

python3 - "$work" <<'PY'
import json
import sys
import numpy as np

work = sys.argv[1]
failed = []


def check(description, condition):
    print(("ok: " if condition else "FAILED: ") + description)
    if not condition:
        failed.append(description)


def load(name):
    return json.load(open(f"{work}/{name}.json"))


def assembled(name):
    """The matrix of a figure as rows by columns, from its tiles."""
    data = load(name)
    spec = data.get("tile_spec")
    if spec is None:
        return np.load(data["matrix"]), data, None
    tiles = spec["tiles"]
    rows, cols = spec["blocks"]
    first = np.load(tiles[0]["file"])
    out = np.zeros((rows, cols), dtype=first.dtype)
    for tile in tiles:
        array = np.load(tile["file"])
        assert array.shape == (tile["c1"] - tile["c0"], tile["r1"] - tile["r0"]), array.shape
        out[tile["r0"]:tile["r1"], tile["c0"]:tile["c1"]] = array.T
    return out, data, spec


full, full_data, full_spec = assembled("full1")
tiled, tiled_data, tiled_spec = assembled("full700")
check("the one tile at full resolution is the single matrix path", full_spec is None)
check("tiles of 700 bins give 5 x 5 tiles", tiled_spec is not None and len(tiled_spec["tiles"]) == 25)
check("tiles put together are the single matrix", full.shape == tiled.shape and (full == tiled).all())
check("counts after --log1p are int32", full.dtype == np.int32 and tiled.dtype == np.int32)
check("--log1p adds 1 to a count matrix whose zeros took the smallest non-zero value",
      full.min() >= 2)
check("no tile spans more than 700 blocks", all(t["r1"] - t["r0"] <= 700 and t["c1"] - t["c0"] <= 700
                                                  for t in tiled_spec["tiles"]))
x_edges = np.array(tiled_spec["x_edges"])
check("edges are increasing and one more than blocks", len(x_edges) == tiled.shape[0] + 1
      and (np.diff(x_edges) > 0).all())
# edges of the single matrix path: the nearest shading edges of start_pos
start = np.array(full_data["start_pos"], dtype=np.float64)
half = np.diff(start) * 0.5
expected = np.concatenate([[start[0] - half[0]], start[:-1] + half, [start[-1] + half[-1]]])
check("block edges at full resolution are the pcolormesh edges of start_pos",
      np.allclose(x_edges, expected))

check("the colour limits of several tiles are those of all of them",
      tiled_spec["vmin"] == full.min() and tiled_spec["vmax"] == full.max())

nolog, _, _ = assembled("nolog")
check("NaN bins keep a float32 matrix without --log1p", nolog.dtype == np.float32 and np.isnan(nolog).any())
nan_rows = np.isnan(nolog).all(axis=1)
check("NaN bins are whole NaN rows and columns", (np.isnan(nolog).all(axis=0) == nan_rows).all()
      and nan_rows.sum() > 0)
check("--log1p leaves no NaN", not np.isnan(full.astype(np.float64)).any())

reduced, reduced_data, spec = assembled("reduced")
off, _, off_spec = assembled("reduced_off")
check("the reduction produces blocks", spec is not None and spec["aggregate"] == "sum"
      and max(spec["bins_per_block"]) > 1)
check("the reduction is switched off by --noDownsample", off_spec is None and off.shape == (3383, 3383))
check("reduced blocks cut at chromosome ends", spec["blocks_by_chromosome"])
check("reduced counts are int32", reduced.dtype == np.int32)

# Without --log1p the sum of all blocks is the sum of all cells of the matrix
# (NaN bins excluded on both sides).
nolog_reduced, _, spec2 = assembled("reduced_nolog")
check("a block with valid bins is a number, so the NaN bins did not spread", np.isnan(nolog_reduced).sum() < np.isnan(nolog).sum() / 100)
check("the total of the blocks is the total of the cells",
      np.nansum(nolog_reduced.astype(np.float64)) == np.nansum(nolog.astype(np.float64)))
full_edges_end = np.array(tiled_spec["x_edges"])[[0, -1]]
reduced_edges = np.array(spec["x_edges"])
check("the first and last edge are those of the full resolution plot",
      reduced_edges[0] == full_edges_end[0] and reduced_edges[-1] == full_edges_end[1])
check("block edges are increasing", (np.diff(reduced_edges) > 0).all())
check("--log1p applies to the sums", reduced.min() >= 2 and reduced.max() > off.max())

if failed:
    sys.exit(1)
PY
status=$?
expect "the python checks passed" [ $status -eq 0 ]

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
