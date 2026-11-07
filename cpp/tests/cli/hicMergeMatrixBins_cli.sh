#!/usr/bin/env bash
# Command line behaviour of hicMergeMatrixBins without --chromosomeSizes:
#
#  1. A run without --chromosomeSizes succeeds on a cool and on an h5 input and
#     writes a matrix.
#  2. On a cool input the result is the same as the result of a run that passes
#     a sizes file holding the chromosome lengths of that cool file.
#  3. --chromosomeSizes is optional in the help text and a sizes file that is
#     given is still used: lengths larger than the matrix's own give a result
#     with more bins.
#
# Usage: hicMergeMatrixBins_cli.sh <hicMergeMatrixBins binary> <test_data directory>
#                                  <hicInfo binary> <hicConvertFormat binary>

set -u
tool=$1
data=$2
info=$3
convert=$4
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
h5="$data/small_test_matrix_50kb_res.h5"

# 1
"$tool" -m "$cool" -o "$work/cool_default.cool" --numBins 4 > /dev/null 2> "$work/err1.txt"
expect "a cool input without --chromosomeSizes succeeds" [ $? -eq 0 ]
expect "the cool result is written" [ -s "$work/cool_default.cool" ]
"$tool" -m "$h5" -o "$work/h5_default.h5" --numBins 4 > /dev/null 2> "$work/err2.txt"
expect "an h5 input without --chromosomeSizes succeeds" [ $? -eq 0 ]
expect "the h5 result is written" [ -s "$work/h5_default.h5" ]

# 2
"$info" -m "$cool" | sed -n 's/^Chromosomes:length: //p' | tr ';' '\n' \
    | sed 's/^ *//; /^$/d; s/: /\t/; s/ bp$//' > "$work/sizes.txt"
expect "the sizes file of the cool input is not empty" [ -s "$work/sizes.txt" ]
"$tool" -m "$cool" -o "$work/cool_sizes.cool" --numBins 4 --chromosomeSizes "$work/sizes.txt" \
    > /dev/null 2> "$work/err3.txt"
expect "a cool input with its own chromosome sizes succeeds" [ $? -eq 0 ]
"$convert" -m "$work/cool_default.cool" --inputFormat cool --outputFormat ginteractions \
    -o "$work/default.gi" > /dev/null 2>&1
"$convert" -m "$work/cool_sizes.cool" --inputFormat cool --outputFormat ginteractions \
    -o "$work/sizes.gi" > /dev/null 2>&1
expect "both results hold the same contacts" [ -s "$work/default.gi" ]
expect "the result without --chromosomeSizes equals the result with the cool's own sizes" \
    cmp -s "$work/default.gi" "$work/sizes.gi"
expect "both results have the same chromosomes and bins" \
    [ "$("$info" -m "$work/cool_default.cool" | grep -E '^(Size|Chromosomes)')" = \
      "$("$info" -m "$work/cool_sizes.cool" | grep -E '^(Size|Chromosomes)')" ]

# 3
expect "--help lists the option among the optional arguments" \
    bash -c "'$tool' --help | sed -n '/Optional arguments/,\$p' | grep -q -- '--chromosomeSizes'"
awk 'BEGIN{OFS="\t"} {print $1, $2 + 1000000}' "$work/sizes.txt" > "$work/longer.txt"
"$tool" -m "$cool" -o "$work/cool_longer.cool" --numBins 4 --chromosomeSizes "$work/longer.txt" \
    > /dev/null 2> "$work/err4.txt"
expect "a longer sizes file is used" [ $? -eq 0 ]
expect "a longer sizes file gives more bins" \
    [ "$("$info" -m "$work/cool_longer.cool" | sed -n 's/^Size:[[:space:]]*//p' | tr -d ',')" -gt \
      "$("$info" -m "$work/cool_default.cool" | sed -n 's/^Size:[[:space:]]*//p' | tr -d ',')" ]

if [ "$failures" -ne 0 ]; then
    echo "$failures failure(s)"
    exit 1
fi
echo "all hicMergeMatrixBins checks passed"
