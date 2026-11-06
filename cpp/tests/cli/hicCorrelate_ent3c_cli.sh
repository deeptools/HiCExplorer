#!/usr/bin/env bash
# hicCorrelate --method ent3c (C++-only, no HiCExplorer Python counterpart):
# ENT3C's entropy-based similarity (hicx::ent3c, core/ent3c.cpp).
#
# The reference values are ENT3C 2.2.2 (PyPI) run on the two GM12878 matrices
# of the test data at 2.5 Mb, with CHRSPLIT 7, phi 1, PHI_MAX 1000 and
# NormM 0, the defaults of hicCorrelate:
#   chromosome 1                Q = 0.9993969242632847
#   chromosome Y                Q = 0.6614378277661476 (two windows of each
#                               signal have no entropy and are left out)
#   mean over 1-20, X and Y     0.9826061586237302
# Chromosomes 21 and 22 are left out of the mean: at 2.5 Mb their windows
# have two bins, every entropy is 0 up to rounding, and Q correlates rounding
# noise (see hicx/ent3c.hpp).
#
# Usage: hicCorrelate_ent3c_cli.sh <hicCorrelate binary> <test_data directory>

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

# close <plot.json> <expected> : results[0][1] within 1e-12 of the expected value
close() {
    python3 - "$1" "$2" <<'EOF'
import json, sys
value = json.load(open(sys.argv[1]))["results"][0][1]
expected = float(sys.argv[2])
if abs(value - expected) > 1e-12:
    print(f"got {value!r}, ENT3C gives {expected!r}")
    sys.exit(1)
EOF
}

primary="$data/hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool"
raw="$data/hicCorrectMatrix/gm12878_raw_values.cool"
small="$data/small_test_matrix_50kb_res.cool"
h5="$data/small_test_matrix.h5"
cd "$work" || exit 1

run() {
    local name=$1
    shift
    "$tool" --method ent3c -oh "$name.png" -os "${name}_scatter.png" --plotData "$name.json" "$@" \
        >"$name.out" 2>"$name.err"
}

run chr1 -m "$primary" "$raw" --chromosomes 1
expect "chromosome 1 exits 0" test $? -eq 0
expect "chromosome 1 matches ENT3C" close chr1.json 0.9993969242632847

run chrY -m "$primary" "$raw" --chromosomes Y
expect "chromosome Y exits 0" test $? -eq 0
expect "chromosome Y matches ENT3C, windows without entropy left out" \
    close chrY.json 0.6614378277661476

autosomes=$(seq 1 20)
# shellcheck disable=SC2086
run mean -m "$primary" "$raw" --chromosomes $autosomes X Y --threads 3
expect "the chromosome mean exits 0" test $? -eq 0
expect "the chromosome mean matches ENT3C's mean" close mean.json 0.9826061586237302
expect "no scatter plot is written" test ! -e mean_scatter.png

# The windows are computed concurrently; the result does not depend on it.
# shellcheck disable=SC2086
run mean1 -m "$primary" "$raw" --chromosomes $autosomes X Y --threads 1
expect "one thread gives the same scores as three, bit for bit" \
    python3 -c "import json,sys; r=[json.load(open(f))['results'] for f in sys.argv[1:]]; sys.exit(r[0] != r[1])" \
    mean.json mean1.json

# A matrix against itself.
run self -m "$small" "$small"
expect "self vs self exits 0" test $? -eq 0
expect "self vs self is exactly 1" close self.json 1.0

# Options reach the computation: a fixed submatrix size changes the score.
run fixed -m "$primary" "$raw" --chromosomes 1 --ent3cSubmatrixSize 5 --ent3cMaxWindows 20
expect "--ent3cSubmatrixSize/--ent3cMaxWindows exit 0" test $? -eq 0
expect "a fixed submatrix size changes the score" \
    python3 -c "import json,sys; sys.exit(json.load(open('fixed.json'))['results'][0][1] == json.load(open('chr1.json'))['results'][0][1])"

# Chromosomes without any score end with a message, not a crash.
run none -m "$primary" "$raw" --chromosomes 22 --ent3cSubmatrixSize 1
expect "no finite score exits 1" test $? -eq 1
expect "the message names the pair" grep -q "no finite score for" none.err

run badshift -m "$primary" "$raw" --ent3cWindowShift 0
expect "a window shift of 0 is refused" test $? -eq 1

# --method ent3c is cool-only, as ENT3C is.
if [ -f "$h5" ]; then
    run refused -m "$h5" "$h5"
    expect "an h5 matrix is refused" test $? -ne 0
    expect "the refusal names the cool-only restriction" grep -qi "cool" refused.err
fi

echo "$failures failing checks"
exit "$failures"
