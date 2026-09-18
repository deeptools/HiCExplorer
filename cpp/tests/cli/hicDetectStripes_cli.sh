#!/usr/bin/env bash
# Command line behaviour of hicDetectStripes (cpp/PLAN.md 9.3, a faithful C++
# reimplementation of Stripenn 1.1.65.22's own detection method), as a ctest
# complement to cpp/scripts/cases/hicDetectStripes.json. The statistical
# validation itself runs through cpp/scripts/stripe_calibration.py on real
# data.
#
# Checks:
#  1. --help-json succeeds and looks like a tool spec.
#  2. A cool matrix with a real chromosome runs to completion, writes its
#     output file and reports a count on stderr.
#  3. An unknown --chromosomes name exits 1 and writes nothing.
#  4. Two runs of the same cool input, at --threads 1 and at --threads 4, are
#     byte identical (cpp/OPTIMIZATION.md 3 determinism).
#  5. An h5 matrix runs the same way as a cool matrix (the whole-matrix
#     loader path, cpp/tools/hicDetectStripes.cpp).
#
# Usage: hicDetectStripes_cli.sh <hicDetectStripes binary> <test_data directory>

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

lenient_args=(--minStripeLength 25000 --maxWidth 6 --canny 1.5 --backgroundSamples 2000
             --pValue 1.0 --fdr 1.0)

# 1
mkdir "$work/help" && cd "$work/help" || exit 1
"$tool" --help-json > out.txt 2> err.txt
status=$?
expect "--help-json exits 0" [ "$status" -eq 0 ]
expect "--help-json names the tool" grep -q '"tool": "hicDetectStripes"' out.txt
expect "--help-json names the schema" grep -q '"schema": "hicexplorer-tool-spec"' out.txt

# 2
mkdir "$work/cool" && cd "$work/cool" || exit 1
"$tool" -m "$data/hicNormalize/small_test_matrix.cool" --chromosomes chr2L \
    "${lenient_args[@]}" -t 1 -o stripes.tsv > out.txt 2> err.txt
status=$?
expect "a cool matrix with a real chromosome exits 0" [ "$status" -eq 0 ]
expect "stderr reports a call count" grep -q "Number of detected stripes:" err.txt
expect "the output file exists" [ -f stripes.tsv ]

# 3
mkdir "$work/missing" && cd "$work/missing" || exit 1
"$tool" -m "$data/hicNormalize/small_test_matrix.cool" --chromosomes chrDoesNotExist \
    -o stripes.tsv > out.txt 2> err.txt
status=$?
expect "an unknown chromosome exits non-zero" [ "$status" -ne 0 ]
expect "the message names the chromosome" grep -q "chrDoesNotExist" err.txt
expect "an unknown chromosome writes nothing" only_logs_in "$work/missing"

# 4
mkdir "$work/determinism" && cd "$work/determinism" || exit 1
"$tool" -m "$data/hicNormalize/small_test_matrix.cool" "${lenient_args[@]}" \
    -t 1 -o t1.tsv > out_t1.txt 2> err_t1.txt
"$tool" -m "$data/hicNormalize/small_test_matrix.cool" "${lenient_args[@]}" \
    -t 4 -o t4.tsv > out_t4.txt 2> err_t4.txt
expect "output is independent of --threads" cmp -s t1.tsv t4.tsv

# 5
mkdir "$work/h5" && cd "$work/h5" || exit 1
"$tool" -m "$data/small_test_matrix.h5" --chromosomes chr2L \
    "${lenient_args[@]}" -t 1 -o stripes.tsv > out.txt 2> err.txt
status=$?
expect "an h5 matrix with a real chromosome exits 0" [ "$status" -eq 0 ]
expect "stderr reports a call count for h5 input" grep -q "Number of detected stripes:" err.txt

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
