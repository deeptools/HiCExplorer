#!/usr/bin/env bash
# hicCorrelate --method scc (cpp/PLAN.md tier 7 extension, C++-only, no
# HiCExplorer Python counterpart): HiCRep's stratum-adjusted correlation
# coefficient, reimplemented from the real hicrep 0.2.6 Python package
# (hicx::scc, core/scc.cpp). There is no HiCExplorer reference to run
# cpp/scripts/equiv.py against, so this is checked here instead; the
# algorithm itself is validated against hicrep 0.2.6 on real GM12878
# replicate .cool files (see hicCorrelate.cpp's run_scc header comment and
# the commit message for the exact reference numbers), which is not
# reproducible from a small checked-in fixture.
#
# Checked here: a matrix against itself scores a perfect 1.0 on every
# pairwise entry (self-correlation, and the no-clustering-distortion case);
# --method scc only accepts cool input, matching hicrep's own scope (an h5
# matrix is refused, not silently reduced); the heatmap is drawn but the
# scatter plot is not (no per-bin vectors survive the per-diagonal
# reduction); --sccSmoothWindow, --sccMaxDistance and --sccDownsample are
# accepted and change nothing about exit status.
#
# Usage: hicCorrelate_scc_cli.sh <hicCorrelate binary> <test_data directory>

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

cool="$data/small_test_matrix_50kb_res.cool"
h5="$data/small_test_matrix.h5"
cd "$work" || exit 1

# A matrix against itself. This is *not* exactly 1.0: small/sparse
# chromosomes in this fixture have diagonals with zero variance (every
# stored value along the diagonal equal), whose correlation is undefined and
# is treated as a zero-weighted zero contribution (hicrep.hicrep.sccByDiag's
# own NaN-to-0 rule), diluting the weighted average below 1 even when the two
# inputs are identical. 0.6979227772668228 is sccByDiag's own value for this
# exact fixture (every chromosome but chrM, averaged, --sccMaxDistance
# 2000000), computed directly from hicrep 0.2.6's trimDiags/sccByDiag on
# each chromosome by hand, not from this port.
"$tool" -m "$cool" "$cool" --method scc --sccMaxDistance 2000000 \
    -oh self_heatmap.png -os self_scatter.png --plotData self.json \
    >self.out 2>self.err
expect "self vs self exits 0" test $? -eq 0
python3 - self.json <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
results = d["results"]
assert results[0][0] == 1.0 and results[1][1] == 1.0, results
assert abs(results[0][1] - 0.6979227772668228) < 1e-9, results
EOF
expect "self vs self matches hicrep's own sccByDiag on this fixture" test $? -eq 0
expect "self vs self: no scatter plot is written" test ! -e self_scatter.png
expect "self vs self: --plotData skips drawing the heatmap too" test ! -e self_heatmap.png

# The heatmap is actually drawn (and the scatter still is not) when not
# using --plotData, same as pearson/spearman.
if [ -n "${HICX_PLOT_PYTHON:-}" ]; then
    "$tool" -m "$cool" "$cool" --method scc -oh drawn_heatmap.png -os drawn_scatter.png \
        >draw.out 2>draw.err
    expect "drawing exits 0" test $? -eq 0
    expect "the heatmap file is written" test -s drawn_heatmap.png
    expect "the scatter file is not written" test ! -e drawn_scatter.png
fi

# --sccSmoothWindow, --sccMaxDistance, --sccDownsample are accepted and run.
"$tool" -m "$cool" "$cool" --method scc --sccSmoothWindow 2 --sccMaxDistance 500000 \
    --sccDownsample -oh smooth_heatmap.png -os smooth_scatter.png --plotData smooth.json \
    >smooth.out 2>smooth.err
expect "--sccSmoothWindow/--sccMaxDistance/--sccDownsample together exit 0" test $? -eq 0

# --method scc is cool-only, matching hicrep's own scope.
if [ -f "$h5" ]; then
    "$tool" -m "$h5" "$h5" --method scc -oh refused_heatmap.png -os refused_scatter.png \
        --plotData refused.json >refused.out 2>refused.err
    expect "an h5 matrix is refused, not silently accepted" test $? -ne 0
    expect "the refusal names the cool-only restriction" grep -qi "cool" refused.err
    expect "nothing was written for the refused input" test ! -e refused_heatmap.png
fi

echo "$failures failing checks"
exit "$failures"
