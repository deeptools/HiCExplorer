#!/usr/bin/env bash
# --nestedMinDepth/--nestedMaxDepth/--nestedStep (C++-only, optional sub-TAD
# detection): within each top-level domain wide enough to plausibly hold
# internal structure, reruns the same TAD-separation-score pipeline at a
# finer depth confined to that domain's own bin range, writing
# <outPrefix>_nested_boundaries.bed, _nested_domains.bed and
# _nested_boundaries.gff. No Python counterpart, so checked here rather than
# through cpp/scripts/equiv.py.
#
# small_test_matrix.h5 (bin size 5000 here) is tuned to produce exactly one
# real boundary and never two boundaries on the same chromosome, so no
# top-level domain ever forms on it and this script cannot exercise genuine
# sub-TAD discovery end to end. That was instead verified manually against
# real GM12878 Hi-C data (gm12878_ENCSR946AKI_r8_25kb.cool, chr1): with
# --nestedMinDepth 75000 --nestedMaxDepth 150000 --nestedStep 25000 on top of
# --minDepth 75000 --maxDepth 250000 --step 25000, 15 sub-boundaries split 11
# top-level domains into 26 sub-domains, each sub-boundary landing on a real
# convergent CTCF pair with p<0.01. What is checked here is the CLI contract:
#  1. --nestedMinDepth/--nestedMaxDepth/--nestedStep must be given together;
#     any two of the three without the third is refused.
#  2. --nestedMinDepth below 3 times the matrix bin size is refused (a
#     Wilcoxon rank-sum test below that has too few values per side to carry
#     any power, the same rule --minDepth itself enforces).
#  3. --nestedMaxDepth below 3 times the matrix bin size is refused likewise.
#  4. Without any --nested* flag, the top-level output is unaffected (same
#     boundary the existing equivalence cases already pin) and no nested
#     files are written.
#  5. With all three given but no top-level domain ever forms (this
#     fixture's case), the run still exits 0, writes the nested output files,
#     and leaves them empty rather than crashing or leaving them missing.
#
# Usage: hicFindTADs_nested_cli.sh <hicFindTADs binary> <test_data directory>

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

matrix="$data/small_test_matrix.h5"
common=(--minDepth 60000 --maxDepth 180000 --step 20000 --minBoundaryDistance 20000
        --correctForMultipleTesting fdr --thresholdComparisons 0.5 --chromosomes chr2L)

cd "$work" || exit 1

"$tool" -m "$matrix" "${common[@]}" --outPrefix base >/dev/null 2>&1
expect "the baseline run (no nested flags) still finds the one real boundary" \
    test "$(wc -l < base_boundaries.bed)" -eq 1
expect "the baseline run writes no nested files" test ! -e base_nested_boundaries.bed

expect "--nestedMinDepth and --nestedMaxDepth without --nestedStep is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 60000 --nestedMaxDepth 120000 \
        --outPrefix partial1 >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e partial1_boundaries.bed

expect "--nestedMinDepth and --nestedStep without --nestedMaxDepth is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 60000 --nestedStep 20000 \
        --outPrefix partial2 >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e partial2_boundaries.bed

expect "--nestedMinDepth below 3x the bin size (5000) is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 10000 --nestedMaxDepth 60000 \
        --nestedStep 10000 --outPrefix toofine_min >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e toofine_min_boundaries.bed

expect "--nestedMaxDepth below 3x the bin size (5000) is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 10000 --nestedMaxDepth 14000 \
        --nestedStep 2000 --outPrefix toofine_max >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e toofine_max_boundaries.bed

"$tool" -m "$matrix" "${common[@]}" --nestedMinDepth 60000 --nestedMaxDepth 120000 \
    --nestedStep 20000 --outPrefix nested_ok >/dev/null 2>&1
expect "a valid nested run (no domain to search) still exits cleanly" \
    test "$(wc -l < nested_ok_boundaries.bed)" -eq 1
expect "it writes an empty nested boundaries file" \
    test -e nested_ok_nested_boundaries.bed -a ! -s nested_ok_nested_boundaries.bed
expect "it writes an empty nested domains file" \
    test -e nested_ok_nested_domains.bed -a ! -s nested_ok_nested_domains.bed
expect "it writes an empty nested boundaries GFF" \
    test -e nested_ok_nested_boundaries.gff -a ! -s nested_ok_nested_boundaries.gff

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
