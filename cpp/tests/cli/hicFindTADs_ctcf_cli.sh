#!/usr/bin/env bash
# --ctcfFile / --ctcfWindow (C++-only, cpp/PLAN.md tier 9): requires a CTCF
# binding site in the correct orientation for a TAD-separation score minimum
# to be accepted as a boundary (the convergent-motif signature of domain
# boundaries). No Python counterpart, so checked here rather than through
# cpp/scripts/equiv.py.
#
# Checks, all on small_test_matrix.h5 chr2L, the single real boundary
# (B04585, breakpoint near 22925000) the default parameters below produce:
#  1. Without --ctcfFile the boundary is unaffected (checked by the existing
#     equivalence case; not repeated here).
#  2. A '+' site at or before the breakpoint, within the window, keeps the
#     boundary and the GFF records which site matched.
#  3. A '-' site at or after the breakpoint, within the window, keeps it too
#     (the same point serves as the right edge of the left domain and the
#     left edge of the right domain).
#  4. The wrong orientation on either side drops the boundary.
#  5. A correctly oriented site outside --ctcfWindow drops it.
#  6. A site with strand '.' (no orientation) drops it.
#  7. An empty --ctcfFile drops every boundary.
#  9. --ctcfBoundarySlack: a site on the wrong side of the exact breakpoint is
#     dropped at slack 0 (the default, the same test as without the flag) and
#     kept once the slack reaches past the site, on both sides, while a site
#     outside --ctcfWindow is still dropped no matter how large the slack is.
# 10. A real convergent pair (both a forward and a reverse site present)
#     still passes under the single-side matching rule, and the GFF records
#     both matched sites.
#  8. A missing --ctcfFile, or one with a line short of the required 6
#     columns, is refused and writes nothing.
#
# Usage: hicFindTADs_ctcf_cli.sh <hicFindTADs binary> <test_data directory>

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
expect "the baseline run finds the one real boundary" test "$(wc -l < base_boundaries.bed)" -eq 1

printf 'chr2L\t22923900\t22924100\tctcf1\t500\t+\n' > left_plus.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile left_plus.bed --ctcfWindow 5000 \
    --outPrefix left_plus >/dev/null 2>&1
expect "a '+' site left of the breakpoint keeps the boundary" \
    test "$(wc -l < left_plus_boundaries.bed)" -eq 1
expect "the GFF records the matched site's position and strand" \
    grep -q 'ctcf_site_forward=22924000' left_plus_boundaries.gff

printf 'chr2L\t22925900\t22926100\tctcf2\t500\t-\n' > right_minus.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile right_minus.bed --ctcfWindow 5000 \
    --outPrefix right_minus >/dev/null 2>&1
expect "a '-' site right of the breakpoint also keeps the boundary" \
    test "$(wc -l < right_minus_boundaries.bed)" -eq 1

printf 'chr2L\t22923900\t22924100\tctcf3\t500\t-\n' > left_minus.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile left_minus.bed --ctcfWindow 5000 \
    --outPrefix left_minus >/dev/null 2>&1
expect "a '-' site left of the breakpoint does not support it" \
    test "$(wc -l < left_minus_boundaries.bed)" -eq 0

printf 'chr2L\t22925900\t22926100\tctcf4\t500\t+\n' > right_plus.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile right_plus.bed --ctcfWindow 5000 \
    --outPrefix right_plus >/dev/null 2>&1
expect "a '+' site right of the breakpoint does not support it" \
    test "$(wc -l < right_plus_boundaries.bed)" -eq 0

printf 'chr2L\t22900000\t22900200\tctcf5\t500\t+\n' > too_far.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile too_far.bed --ctcfWindow 5000 \
    --outPrefix too_far >/dev/null 2>&1
expect "a correctly oriented site outside the window does not support it" \
    test "$(wc -l < too_far_boundaries.bed)" -eq 0

printf 'chr2L\t22923900\t22924100\tctcf6\t500\t.\n' > dot_strand.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile dot_strand.bed --ctcfWindow 5000 \
    --outPrefix dot_strand >/dev/null 2>&1
expect "a '.' (no orientation) site does not support the boundary" \
    test "$(wc -l < dot_strand_boundaries.bed)" -eq 0

: > empty.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile empty.bed --outPrefix empty >/dev/null 2>&1
expect "an empty --ctcfFile drops every boundary" test "$(wc -l < empty_boundaries.bed)" -eq 0

expect "a missing --ctcfFile is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --ctcfFile /nonexistent.bed --outPrefix missing >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e missing_boundaries.bed

printf 'chr2L\t100\t200\n' > short.bed
expect "a short line (fewer than 6 columns) is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --ctcfFile short.bed --outPrefix short_run >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e short_run_boundaries.bed

cat left_plus.bed right_minus.bed > both.bed
"$tool" -m "$matrix" "${common[@]}" --ctcfFile both.bed --ctcfWindow 5000 --outPrefix pair_both >/dev/null 2>&1
expect "a real convergent pair (both sites present) still passes" \
    test "$(wc -l < pair_both_boundaries.bed)" -eq 1
expect "the GFF records both matched sites" bash -c \
    "grep -q 'ctcf_site_forward=22924000' pair_both_boundaries.gff && \
     grep -q 'ctcf_site_reverse=22926000' pair_both_boundaries.gff"

"$tool" -m "$matrix" "${common[@]}" --ctcfFile right_plus.bed --ctcfWindow 5000 \
    --ctcfBoundarySlack 0 --outPrefix slack0 >/dev/null 2>&1
expect "slack 0 rejects a '+' site right of the breakpoint, same as no slack flag" \
    test "$(wc -l < slack0_boundaries.bed)" -eq 0
"$tool" -m "$matrix" "${common[@]}" --ctcfFile right_plus.bed --ctcfWindow 5000 \
    --ctcfBoundarySlack 5000 --outPrefix slack_right >/dev/null 2>&1
expect "slack past the site accepts a '+' site right of the breakpoint" \
    test "$(wc -l < slack_right_boundaries.bed)" -eq 1
"$tool" -m "$matrix" "${common[@]}" --ctcfFile left_minus.bed --ctcfWindow 5000 \
    --ctcfBoundarySlack 5000 --outPrefix slack_left >/dev/null 2>&1
expect "slack past the site accepts a '-' site left of the breakpoint" \
    test "$(wc -l < slack_left_boundaries.bed)" -eq 1
"$tool" -m "$matrix" "${common[@]}" --ctcfFile too_far.bed --ctcfWindow 5000 \
    --ctcfBoundarySlack 50000 --outPrefix slack_far >/dev/null 2>&1
expect "--ctcfWindow still bounds distance no matter how large the slack is" \
    test "$(wc -l < slack_far_boundaries.bed)" -eq 0

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
