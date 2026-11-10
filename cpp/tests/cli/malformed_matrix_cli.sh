#!/usr/bin/env bash
# Matrix files whose sparse structure points outside the matrix are rejected
# with an error message and exit status 1. Before the structure check such a
# file made the tools index past their arrays.
#
#  1. A cool file with a pixel whose bin2_id lies far outside the bin table.
#  2. An h5 file with a column index outside the matrix.
#  3. An h5 file whose last row offset points past the stored entries.
#  4. Well-formed cool and h5 files of the same layout load.
#
# The files are test_data/malformed, written by
# cpp/scripts/make_malformed_fixtures.py.
#
# Usage: malformed_matrix_cli.sh <hicTransform binary> <test_data directory>

set -u
transform=$1
data=$2/malformed
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

# rejected <file> <label>: hicTransform exits 1 (not by a signal) and names the problem.
rejected() {
    local file=$1 label=$2
    "$transform" -m "$file" -o "$work/out_$label.cool" --method obs_exp \
        > /dev/null 2> "$work/err_$label.txt"
    local status=$?
    expect "$label: hicTransform exits 1, not by a signal (got $status)" [ "$status" -eq 1 ]
    expect "$label: the error names the invalid structure" \
        grep -q "not a valid" "$work/err_$label.txt"
}

# accepted <file> <label>: a well-formed file of the same layout loads.
accepted() {
    local file=$1 label=$2
    "$transform" -m "$file" -o "$work/out_$label.${file##*.}" --method obs_exp \
        > /dev/null 2> "$work/err_$label.txt"
    local status=$?
    expect "$label: hicTransform reads the well-formed file (got $status)" [ "$status" -eq 0 ]
}

rejected "$data/bad_column.cool" cool_column
rejected "$data/bad_column.h5" h5_column
rejected "$data/bad_indptr.h5" h5_indptr
accepted "$data/good.cool" cool_good
accepted "$data/good.h5" h5_good

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
