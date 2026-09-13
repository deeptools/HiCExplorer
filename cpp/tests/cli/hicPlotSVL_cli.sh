#!/usr/bin/env bash
# Command line behaviour of the C++ hicPlotSVL that the harness cannot compare
# against the Python, because the Python either draws a figure or never ends:
#
#  1. --plotFileName is refused (contract rule 7): non-zero exit, a message
#     naming the cause, and nothing written into the working directory.
#  2. Without --plotFileName the data file is written, the plot is skipped with
#     a note on stderr, and no plot.png appears.
#  3. A chromosome that is not in an h5 matrix: the Python raises inside a
#     worker process and then polls for it forever (hicPlotSVL.py:183-196).
#     The port exits non-zero with the cause and writes nothing (a deliberate
#     deviation approved 2026-09-13).
#  4. --threads 0 fails as the Python's ZeroDivisionError does.
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
mkdir "$work/plot" && cd "$work/plot" || exit 1
"$tool" -m "$data/small_test_matrix_50kb_res.h5" --plotFileName plot.png \
    -o p_values.txt -od data.txt > out.txt 2> err.txt
status=$?
expect "--plotFileName exits non-zero" [ "$status" -ne 0 ]
expect "--plotFileName names the cause" grep -q "plotting is not yet available" err.txt
expect "--plotFileName writes nothing" only_logs_in "$work/plot"

"$tool" -m "$data/small_test_matrix_50kb_res.h5" -pfn other.svg > out.txt 2> err.txt
status=$?
expect "-pfn exits non-zero" [ "$status" -ne 0 ]
expect "-pfn writes nothing, not even the default outputs" only_logs_in "$work/plot"

# 2
mkdir "$work/data" && cd "$work/data" || exit 1
"$tool" -m "$data/small_test_matrix_50kb_res.h5" "$data/small_test_matrix_50kb_res.h5" \
    > out.txt 2> err.txt
status=$?
expect "without --plotFileName exits 0" [ "$status" -eq 0 ]
expect "the default data file is written" [ -s data.txt ]
expect "the default p-value file is written" [ -s p_values.txt ]
expect "no plot is written" [ ! -e plot.png ]
expect "the skipped plot is noted on stderr" grep -q "box plot" err.txt

# 3
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

# 4
mkdir "$work/threads" && cd "$work/threads" || exit 1
"$tool" -m "$data/small_test_matrix_50kb_res.h5" --threads 0 -o p_values.txt -od data.txt \
    > out.txt 2> err.txt
status=$?
expect "--threads 0 exits non-zero" [ "$status" -ne 0 ]
expect "--threads 0 writes nothing" only_logs_in "$work/threads"

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
