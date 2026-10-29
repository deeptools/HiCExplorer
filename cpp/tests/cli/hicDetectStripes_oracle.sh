#!/usr/bin/env bash
# hicDetectStripes against Stripenn 1.1.65.22 on a small synthetic cooler
# (tests/fixtures/stripenn/tiny_stripes.cool: three chromosomes with planted
# stripes, balancing weights with NaN bins, chromosome sizes that are not a
# multiple of the bin size). The expected tables are what `stripenn compute`
# wrote for it (cpp/scripts/stripenn_oracle/make_fixtures.py), with Stripenn's
# defaults and two cores. Acceptance gate: identical output.
#
# Checks:
#  1. The unfiltered candidate table is byte identical to Stripenn's
#     result_unfiltered.tsv: same candidates, p-values, Mean and Stripiness.
#  2. The called stripes are the rows of Stripenn's result_filtered.tsv, in
#     its order (descending Stripiness).
#  3. The same at --threads 1 and 4 (and with --sequentialRandom the table
#     differs, which shows the option is wired).
#
# Usage: hicDetectStripes_oracle.sh <hicDetectStripes binary> <fixture directory>

set -u
tool=$1
fixtures=$2
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

cd "$work" || exit 1
"$tool" -m "$fixtures/tiny_stripes.cool" -o called_t1.tsv --unfilteredFileName unfiltered_t1.tsv -t 1 \
    > out1.txt 2> err1.txt
expect "the run exits 0" [ -f called_t1.tsv ]
expect "the unfiltered table is byte identical to Stripenn's result_unfiltered.tsv" \
    cmp -s unfiltered_t1.tsv "$fixtures/tiny_unfiltered.tsv"

# called stripes: chrom, anchor_start + 1 = pos1, anchor_end = pos2, in Stripenn's order
awk -F'\t' 'NR > 1 {print $1 "\t" $2 + 1 "\t" $3}' called_t1.tsv > called_keys.txt
awk -F'\t' 'NR > 1 {print $1 "\t" $2 "\t" $3}' "$fixtures/tiny_filtered.tsv" > expected_keys.txt
expect "the called stripes are Stripenn's result_filtered rows in its order" cmp -s called_keys.txt expected_keys.txt
expect "the output has the nine documented columns" [ "$(head -1 called_t1.tsv | awk -F'\t' '{print NF}')" = "9" ]

"$tool" -m "$fixtures/tiny_stripes.cool" -o called_t4.tsv --unfilteredFileName unfiltered_t4.tsv -t 4 \
    > out4.txt 2> err4.txt
expect "output is independent of --threads" cmp -s unfiltered_t1.tsv unfiltered_t4.tsv

"$tool" -m "$fixtures/tiny_stripes.cool" -o called_seq.tsv --unfilteredFileName unfiltered_seq.tsv -t 2 \
    --sequentialRandom > outs.txt 2> errs.txt
expect "--sequentialRandom draws another background" bash -c '! cmp -s unfiltered_seq.tsv unfiltered_t1.tsv'

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
