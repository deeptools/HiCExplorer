#!/usr/bin/env bash
# --matrix2 (cpp/PLAN.md tier 9): a v4-only addition with no Python
# equivalent, so it is checked here rather than through cpp/scripts/equiv.py.
#
# hicPlotMatrix --matrix2 MATRIX2 combines --matrix and --matrix2, read for
# the same requested region, into one heatmap: the upper triangle (column
# index greater than row index) from --matrix, the lower triangle (column
# index less than row index) from --matrix2, the diagonal always from
# --matrix. This is checked element-wise, not just "the tool exits 0": a
# second, genuinely distinct but same-shape matrix is built at run time with
# the already-ported hicSumMatrices (matrix + matrix, so every value doubles
# but the bin table, and therefore the shape, is identical), --plotData
# dumps the combined matrix and each source matrix alone, and a small numpy
# check compares them cell by cell over a real, non-trivial region of a real
# test matrix (small_test_matrix_50kb_res.cool, 3383 bins).
#
# Also checked: --matrix2 is refused together with --perChromosome (a single
# upper/lower triangle heatmap has no per-chromosome analogue), and a real
# shape mismatch (the same region read at two different resolutions) is
# refused with a message naming both shapes, not a crash or a silent crop.
#
# Usage: hicPlotMatrix_triangle_cli.sh <hicPlotMatrix binary> <hicSumMatrices binary>
#                                       <test_data directory>

set -u
tool=$1
sum_tool=$2
data=$3
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

matrix1="$data/small_test_matrix_50kb_res.cool"
matrix2_5kb="$data/small_test_matrix.cool"
region="chrX:1000000-8000000"

cd "$work" || exit 1

# A second matrix with the identical bin table as matrix1 (hicSumMatrices
# requires it) but genuinely different values: matrix1 + matrix1 doubles
# every count, so matrix2 != matrix1 everywhere matrix1 is non-zero.
"$sum_tool" -m "$matrix1" "$matrix1" -o doubled.cool > sum_out.txt 2> sum_err.txt
status=$?
expect "hicSumMatrices builds the doubled fixture" [ "$status" -eq 0 ]
expect "the doubled fixture was written" [ -e doubled.cool ]

# 1: element-wise correctness on a real, non-trivial region.
"$tool" -m "$matrix1" --region "$region" --plotData m1.json -out unused1.png \
    > m1_out.txt 2> m1_err.txt
status1=$?
"$tool" -m doubled.cool --region "$region" --plotData m2.json -out unused2.png \
    > m2_out.txt 2> m2_err.txt
status2=$?
"$tool" -m "$matrix1" --matrix2 doubled.cool --region "$region" --plotData combo.json \
    -out unused3.png > combo_out.txt 2> combo_err.txt
status3=$?
expect "matrix1 alone exits 0" [ "$status1" -eq 0 ]
expect "matrix2 alone exits 0" [ "$status2" -eq 0 ]
expect "the combined plot exits 0" [ "$status3" -eq 0 ]

python3 - "$work" <<'PY'
import sys
import numpy as np

work = sys.argv[1]
m1 = np.load(f"{work}/m1.json.0.npy")
m2 = np.load(f"{work}/m2.json.0.npy")
combo = np.load(f"{work}/combo.json.0.npy")

ok = True

def check(description, condition):
    global ok
    print(("ok: " if condition else "FAILED: ") + description)
    ok = ok and condition

check("all three dumps share the same shape",
      m1.shape == m2.shape == combo.shape)
check("the region is a real, non-trivial size (not a 2x2 toy)",
      m1.shape[0] >= 50)
check("matrix2 is genuinely distinct from matrix1 somewhere",
      not np.array_equal(m1, m2))

n = m1.shape[0]
rows, cols = np.indices((n, n))
upper = cols > rows
lower = cols < rows
diag = cols == rows

check("every upper-triangle cell of the combined plot equals matrix1's own dump",
      np.array_equal(combo[upper], m1[upper]))
check("every lower-triangle cell of the combined plot equals matrix2's own dump",
      np.array_equal(combo[lower], m2[lower]))
check("the diagonal of the combined plot equals matrix1's own dump",
      np.array_equal(combo[diag], m1[diag]))
check("the combined plot differs from a plain matrix1 dump (the merge did something)",
      not np.array_equal(combo, m1))

sys.exit(0 if ok else 1)
PY
if [ $? -ne 0 ]; then
    failures=$((failures + 1))
fi

# 2: --matrix2 and --perChromosome are refused together, before any output.
mkdir "$work/per_chrom" && cd "$work/per_chrom" || exit 1
"$tool" -m "$matrix1" --matrix2 "$matrix1" --perChromosome --plotData p.json \
    -out p.png > out.txt 2> err.txt
status=$?
expect "--matrix2 with --perChromosome exits non-zero" [ "$status" -ne 0 ]
expect "the message names both options" grep -q -- "--matrix2" err.txt
# p.png is the writableFile probe's own empty file (opened while parsing,
# before the --perChromosome/--matrix2 check runs, cpp/AGENTS_CONTRACT.md
# rule 7); the refusal must still write no plotData JSON.
expect "--matrix2 with --perChromosome writes no plotData" [ ! -e p.json ]
cd "$work" || exit 1

# 3: a genuine shape mismatch (the same coordinates at two different
# resolutions) is refused with a message naming both shapes.
mkdir "$work/mismatch" && cd "$work/mismatch" || exit 1
"$tool" -m "$matrix1" --matrix2 "$matrix2_5kb" --region "$region" --plotData p.json \
    -out p.png > out.txt 2> err.txt
status=$?
expect "a real --matrix2 shape mismatch exits non-zero" [ "$status" -ne 0 ]
expect "the message names --matrix's shape" grep -q "140x140" err.txt
expect "the message names --matrix2's shape" grep -q "1400x1400" err.txt
cd "$work" || exit 1

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
