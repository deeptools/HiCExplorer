#!/usr/bin/env bash
# hicArrowhead (cpp/PLAN.md tier 9): a faithful reimplementation of
# Juicer/Juicebox's own Arrowhead contact-domain caller -- HiCExplorer has
# no Arrowhead caller of its own (hicFindTADs uses an independent
# insulation-score method), so there is no Python behaviour to reproduce
# and no cpp/scripts/equiv.py case; checked here instead. There is also no
# local Java/juicer_tools environment to generate reference domain calls
# from the real tool, so this exercises the CLI end to end (exit codes,
# option plumbing, output format) and a basic sanity check that the
# algorithm calls real-looking domains on real data, rather than pinning
# exact values; see tools/arrowhead_impl.hpp for the algorithm's own unit
# tests (hand-computed primitives, synthetic planted-block pipeline check)
# and its documented deviations from the real tool.
#
# Usage: hicArrowhead_cli.sh <hicArrowhead binary> <test_data directory>

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

cool="$data/small_test_matrix_50kb_res.cool"
cd "$work" || exit 1

# Default settings (matrixWidth=2000, minBlockSize=60) assume real TAD-scale
# windows; this fixture is too small for them, but the tool should still
# exit cleanly and write just the header.
"$tool" -m "$cool" -o default.bedpe >default.out 2>default.err
expect "default settings exit 0" test $? -eq 0
expect "the header is always written" test -s default.bedpe
expect "the header names every column" grep -q "uVarScore" default.bedpe

# Scaled down to the fixture's own size, it finds real domains.
"$tool" -m "$cool" -o scaled.bedpe --matrixWidth 40 --minBlockSize 2 \
    >scaled.out 2>scaled.err
expect "scaled-down settings exit 0" test $? -eq 0
domain_lines=$(grep -vc '^#' scaled.bedpe)
expect "scaled-down settings find at least one domain" test "$domain_lines" -gt 0
expect "every domain row has 8 tab-separated fields" bash -c \
    "awk -F'\t' '!/^#/ && NF != 8 {exit 1}' scaled.bedpe"
expect "x2 is always at or past x1 (domains are square, on the diagonal)" bash -c \
    "awk -F'\t' '!/^#/ && \$3 < \$2 {exit 1}' scaled.bedpe"

# --chromosomes restricts the output to the requested chromosome(s).
"$tool" -m "$cool" -o restricted.bedpe --matrixWidth 40 --minBlockSize 2 \
    --chromosomes chr2L >restricted.out 2>restricted.err
expect "--chromosomes exits 0" test $? -eq 0
expect "every domain is on the requested chromosome" bash -c \
    "awk -F'\t' '!/^#/ && \$1 != \"chr2L\" {exit 1}' restricted.bedpe"

# A nonexistent matrix file is refused, not silently accepted.
"$tool" -m "$work/does_not_exist.cool" -o refused.bedpe >refused.out 2>refused.err
expect "a missing matrix file is refused" test $? -ne 0
expect "nothing is written for the refused input" test ! -e refused.bedpe

echo "$failures failing checks"
exit "$failures"
