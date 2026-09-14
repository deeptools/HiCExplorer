#!/usr/bin/env bash
# Command line behaviour of the C++ hicPlotSVL that the harness cannot compare
# against the Python, because the Python never ends:
#
#  1. A chromosome that is not in an h5 matrix: the Python raises inside a
#     worker process and then polls for it forever (hicPlotSVL.py:183-196).
#     The port exits non-zero with the cause and writes nothing (a deliberate
#     deviation approved 2026-09-13).
#  2. --threads 0 fails as the Python's ZeroDivisionError does.
#
# Usage: hicPlotSVL_cli.sh <hicPlotSVL binary> <test_data directory>

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

only_logs_in() {
    [ -z "$(ls -A "$1" | grep -v -x -e out.txt -e err.txt)" ]
}

# 1
mkdir "$work/missing" && cd "$work/missing" || exit 1
"$tool" -m "$data/small_test_matrix_50kb_res.h5" --chromosomes chr2L chrNotThere \
    -o p_values.txt -od data.txt > out.txt 2> err.txt
status=$?
expect "a chromosome missing from an h5 matrix exits non-zero" [ "$status" -ne 0 ]
expect "the message names the chromosome" grep -q "chrNotThere" err.txt
expect "the message names the reference behaviour" grep -q "waits for that worker forever" err.txt
expect "a missing chromosome writes nothing" only_logs_in "$work/missing"

"$tool" -m "$data/small_test_matrix_50kb_res.cool" --chromosomes chrNotThere \
    -o p_values.txt -od data.txt > out.txt 2> err.txt
status=$?
expect "a chromosome missing from a cool matrix exits non-zero" [ "$status" -ne 0 ]
expect "a missing cool chromosome writes nothing" only_logs_in "$work/missing"

# 2
mkdir "$work/threads" && cd "$work/threads" || exit 1
"$tool" -m "$data/small_test_matrix_50kb_res.h5" --threads 0 -o p_values.txt -od data.txt \
    > out.txt 2> err.txt
status=$?
expect "--threads 0 exits non-zero" [ "$status" -ne 0 ]
expect "--threads 0 writes nothing" only_logs_in "$work/threads"

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
