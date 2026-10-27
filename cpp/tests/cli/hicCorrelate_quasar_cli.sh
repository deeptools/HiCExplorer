#!/usr/bin/env bash
# hicCorrelate --method quasar (cpp/PLAN.md tier 7 extension, C++-only, no
# HiCExplorer Python counterpart): QuASAR-Rep's replicate reproducibility
# score, reimplemented from the real hifive 1.5.7 package's compiled core
# (hicx::quasar, core/quasar.cpp). hifive is Python 2 only, so unlike
# hicx::scc/hicx::hicspector's references this cannot be cross-checked by
# the tests or build here, but it does install and run (under a plain
# Python 2.7 + numpy/scipy environment) elsewhere, and was called there
# directly -- hifive.libraries._quasar.filter_bins,
# .find_bg_dist_norm and .find_correlations, bypassing hifive's own
# Project/Fend/HiC file-format machinery -- as this port's own validation
# oracle on real GM12878/A549 data (see hicCorrelate.cpp's run_quasar
# header comment and the commit message for the exact reference numbers),
# not reproducible from a small checked-in fixture.
#
# Checked here: a matrix against itself scores exactly 1.0 on every
# pairwise entry (QuASAR-Rep's score, like GenomeDISCO's and HiC-Spector's
# and unlike SCC's, is exactly 1.0 for identical inputs, since the two
# "trans" value sets are then identical and their Pearson correlation is
# exactly 1); --method quasar only accepts cool input, matching the real
# tool's own scope; the heatmap is drawn but the scatter plot is not;
# --quasarWidth is accepted and changes nothing about exit status.
#
# Usage: hicCorrelate_quasar_cli.sh <hicCorrelate binary> <test_data directory>

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

# A matrix against itself: every pairwise score is exactly 1.0.
"$tool" -m "$cool" "$cool" --method quasar \
    -oh self_heatmap.png -os self_scatter.png --plotData self.json \
    >self.out 2>self.err
expect "self vs self exits 0" test $? -eq 0
python3 - self.json <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
results = d["results"]
assert results == [[1.0, 1.0], [1.0, 1.0]], results  # symmetric: the C++ passes the full matrix
EOF
expect "self vs self is exactly 1.0 on every entry" test $? -eq 0
expect "self vs self: no scatter plot is written" test ! -e self_scatter.png
expect "self vs self: --plotData skips drawing the heatmap too" test ! -e self_heatmap.png

# The heatmap is actually drawn (and the scatter still is not) when not
# using --plotData, same as pearson/spearman/scc/genomedisco/hicspector.
if [ -n "${HICX_PLOT_PYTHON:-}" ]; then
    "$tool" -m "$cool" "$cool" --method quasar -oh drawn_heatmap.png -os drawn_scatter.png \
        >draw.out 2>draw.err
    expect "drawing exits 0" test $? -eq 0
    expect "the heatmap file is written" test -s drawn_heatmap.png
    expect "the scatter file is not written" test ! -e drawn_scatter.png
fi

# --quasarWidth is accepted.
"$tool" -m "$cool" "$cool" --method quasar --quasarWidth 20 \
    -oh width_heatmap.png -os width_scatter.png --plotData width.json \
    >width.out 2>width.err
expect "--quasarWidth exits 0" test $? -eq 0

# --method quasar is cool-only, matching the real tool's own scope.
if [ -f "$h5" ]; then
    "$tool" -m "$h5" "$h5" --method quasar -oh refused_heatmap.png -os refused_scatter.png \
        --plotData refused.json >refused.out 2>refused.err
    expect "an h5 matrix is refused, not silently accepted" test $? -ne 0
    expect "the refusal names the cool-only restriction" grep -qi "cool" refused.err
    expect "nothing was written for the refused input" test ! -e refused_heatmap.png
fi

echo "$failures failing checks"
exit "$failures"
