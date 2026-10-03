#!/usr/bin/env bash
# --nestedMinDepth/--nestedMaxDepth/--nestedStep/--nestedCoarsenFactor
# (C++-only, optional sub-TAD detection): coarsens the matrix's raw,
# pre-z-score counts by --nestedCoarsenFactor adjacent bins (once, for the
# whole chromosome, not per domain) and re-derives the observed/expected
# z-score fresh on that coarsened matrix, independently per diagonal, the
# same way and with the same chromosome-wide independence the top level's
# own z-score already has. For each top-level domain wide enough to
# plausibly hold internal structure, the TAD-separation-score -> minima ->
# delta -> Wilcoxon -> FDR/Bonferroni -> (optional CTCF) pipeline then reruns
# confined to that domain's own bin range within this coarsened matrix,
# writing <outPrefix>_nested_boundaries.bed, _nested_domains.bed and
# _nested_boundaries.gff. No Python counterpart, so checked here rather than
# through cpp/scripts/equiv.py.
#
# small_test_matrix.h5 (bin size 5000 here) is tuned to produce exactly one
# real boundary and never two boundaries on the same chromosome, so no
# top-level domain ever forms on it and this script cannot exercise genuine
# sub-TAD discovery end to end. That was instead verified manually against
# real GM12878 Hi-C data (gm12878_ENCSR946AKI_r8_25kb.cool, chr1): with
# --nestedMinDepth 300000 --nestedMaxDepth 600000 --nestedStep 100000
# --nestedCoarsenFactor 4 on top of --minDepth 75000 --maxDepth 250000
# --step 25000 and --ctcfRequireConvergentPair, 11 sub-boundaries split 9 of
# 19 qualifying top-level domains into 20 sub-domains, each sub-boundary
# landing on a real convergent CTCF pair with p<0.01 (several below 1e-3).
#
# An earlier version of this feature re-derived the z-score per domain, from
# only that domain's own diamond (bin range plus the --nestedMaxDepth
# margin). That is circular: estimating a diagonal's mean and standard
# deviation from the very region being tested for internal structure lets a
# real sub-boundary partly determine its own "expected" baseline, which gave
# implausibly strong, not-to-be-trusted p-values and is not a sound null
# model. A version before that reused the genome-wide z-scored matrix at a
# finer window instead of recomputing anything; that kept a sound,
# domain-independent baseline but only ever evaluated single raw-resolution
# cells per window, which gave much weaker p-values (~0.01-0.3) and put
# sub-boundaries in areas that did not track real local structure. Coarsening
# first and then z-scoring chromosome-wide, as this version does, keeps the
# sound, domain-independent baseline while still reducing the per-cell noise
# the finer window struggled with.
#
# What is checked here is the CLI contract:
#  1. --nestedMinDepth/--nestedMaxDepth/--nestedStep must be given together;
#     any two of the three without the third is refused.
#  2. --nestedMinDepth below 3 times the coarsened bin size (the matrix bin
#     size times --nestedCoarsenFactor, default 4) is refused: a Wilcoxon
#     rank-sum test below that has too few values per side to carry any
#     power, the same rule --minDepth itself enforces against the matrix's
#     own bin size.
#  3. --nestedMaxDepth below that same floor is refused likewise.
#  4. --nestedCoarsenFactor below 1 is refused.
#  5. --nested* combined with --TAD_sep_score_prefix (reusing a cached
#     z-score matrix/score file instead of computing one fresh) is refused,
#     because that cached file keeps no raw counts for the per-domain
#     coarsened re-analysis to extract.
#  6. Without any --nested* flag, the top-level output is unaffected (same
#     boundary the existing equivalence cases already pin) and no nested
#     files are written.
#  7. With all three depth flags given but no top-level domain ever forms
#     (this fixture's case), the run still exits 0, writes the nested output
#     files, and leaves them empty rather than crashing or leaving them
#     missing.
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

# Default --nestedCoarsenFactor is 4, so the floor is 3 * (5000 * 4) = 60000.
expect "--nestedMinDepth below the coarsened-bin-size floor is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 40000 --nestedMaxDepth 120000 \
        --nestedStep 20000 --outPrefix toofine_min >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e toofine_min_boundaries.bed

expect "--nestedMaxDepth below the coarsened-bin-size floor is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 20000 --nestedMaxDepth 40000 \
        --nestedStep 10000 --outPrefix toofine_max >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e toofine_max_boundaries.bed

# With --nestedCoarsenFactor 1 the floor drops back to 3 * 5000 = 15000.
expect "--nestedCoarsenFactor 1 lowers the floor back to 3x the raw bin size" bash -c \
    "'$tool' -m '$matrix' ${common[*]} --nestedMinDepth 20000 --nestedMaxDepth 60000 \
        --nestedStep 10000 --nestedCoarsenFactor 1 --outPrefix coarsen1 >/dev/null 2>&1"
expect "it still finds the one real top-level boundary" \
    test "$(wc -l < coarsen1_boundaries.bed)" -eq 1

expect "--nestedCoarsenFactor 0 is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --nestedMinDepth 60000 --nestedMaxDepth 120000 \
        --nestedStep 20000 --nestedCoarsenFactor 0 --outPrefix coarsen0 >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e coarsen0_boundaries.bed

"$tool" -m "$matrix" "${common[@]}" --outPrefix cached >/dev/null 2>&1
expect "--nested* combined with --TAD_sep_score_prefix is refused" bash -c \
    "! '$tool' -m '$matrix' ${common[*]} --TAD_sep_score_prefix cached \
        --nestedMinDepth 60000 --nestedMaxDepth 120000 --nestedStep 20000 \
        --outPrefix reuse_nested >/dev/null 2>&1"
expect "the refused run wrote no boundaries file" test ! -e reuse_nested_boundaries.bed

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
