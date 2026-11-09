#!/usr/bin/env bash
# hicAdjustMatrix with --regions (issue 934):
#  1. a region whose end equals the chromosome length keeps the last bin, so
#     the result is the whole chromosome, as with an end one base earlier;
#  2. --action keep, remove and mask exit 1 and write nothing when the
#     regions match no bin of the matrix.
#
# Usage: hicAdjustMatrix_cli.sh <hicAdjustMatrix> <hicInfo> <test_data directory>

set -u
tool=$1
info=$2
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

size_of() { "$info" -m "$1" | sed -n 's/^Size:[[:space:]]*//p'; }

printf 'chr4\t0\t1351857\n' > "$work/end.bed"
printf 'chr4\t0\t1351856\n' > "$work/before.bed"
"$tool" -m "$data/small_test_matrix.h5" --action keep --regions "$work/end.bed" -o "$work/end.h5" 2> /dev/null
expect "an end at the chromosome length exits 0" [ $? -eq 0 ]
"$tool" -m "$data/small_test_matrix.h5" --action keep --regions "$work/before.bed" -o "$work/before.h5" 2> /dev/null
expect "the end region keeps all bins of the chromosome" [ "$(size_of "$work/end.h5")" = "271" ]
expect "the same size as an end one base earlier" [ "$(size_of "$work/end.h5")" = "$(size_of "$work/before.h5")" ]

printf 'chr4\t5000000\t6000000\n' > "$work/outside.bed"
for action in keep remove mask; do
    "$tool" -m "$data/small_test_matrix.h5" --action "$action" --regions "$work/outside.bed" \
        -o "$work/out_$action.h5" > /dev/null 2> "$work/err_$action.txt"
    code=$?
    expect "$action with a region that matches no bin exits 1" [ "$code" -eq 1 ]
    expect "$action with a region that matches no bin writes nothing" [ ! -e "$work/out_$action.h5" ]
done

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
