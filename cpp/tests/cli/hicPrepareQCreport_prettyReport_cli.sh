#!/usr/bin/env bash
# hicPrepareQCreport --prettyReport (C++-only, no Python counterpart): a
# second HTML report at OUTPUTFOLDER/hicQC_report.html, independent of and
# unchanged from the existing hicQC.html the tool already draws through the
# Python plotting layer. It mirrors hicQC.html's own sections (pairs
# sequenced, unmappable/non-unique, pairs discarded with its definitions
# list, contact distance, read orientation, the full QC_table.txt) as native
# CSS bar charts and sortable tables instead of static matplotlib PNGs, adds
# an overview with five metrics per sample judged PASS/WARN/FAIL against
# --min*/--max*Percent cutoffs, and still links to the same underlying .txt
# tables write_qc_report_tables writes for download. Checked here rather
# than through cpp/scripts/equiv.py since it has no reference to match.
#
# Usage: hicPrepareQCreport_prettyReport_cli.sh <hicPrepareQCreport binary> <test_data directory>

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

log="$data/QC/QC.log"
cd "$work" || exit 1

"$tool" --logfiles "$log" --outputFolder base --plotData base/plotdata.json >/dev/null 2>&1
expect "without --prettyReport, no hicQC_report.html is written" test ! -e base/hicQC_report.html

"$tool" --logfiles "$log" --outputFolder pretty --prettyReport \
    --plotData pretty/plotdata.json >/dev/null 2>&1
expect "with --prettyReport, hicQC_report.html is written" test -e pretty/hicQC_report.html
expect "it writes a PASS badge for at least one metric" \
    grep -q 'qc-badge qc-pass">PASS' pretty/hicQC_report.html
expect "hicQC.html itself is unaffected (still delegated, not drawn without a Python environment here)" \
    test ! -e pretty/hicQC.html

"$tool" --logfiles "$log" --outputFolder strict --prettyReport \
    --maxInterChromosomalPercent 0.0001 --plotData strict/plotdata.json >/dev/null 2>&1
expect "an unreachable --maxInterChromosomalPercent cutoff FAILs that metric" \
    grep -q 'Inter-chromosomal contacts</span><span class="qc-badge qc-fail">FAIL' strict/hicQC_report.html

"$tool" --logfiles "$log" --outputFolder lenient --prettyReport \
    --maxInterChromosomalPercent 100 --plotData lenient/plotdata.json >/dev/null 2>&1
expect "a 100% --maxInterChromosomalPercent cutoff PASSes that metric" \
    grep -q 'Inter-chromosomal contacts</span><span class="qc-badge qc-pass">PASS' lenient/hicQC_report.html

for section in "Pairs sequenced" "Unmappable and non-unique pairs" "Pairs discarded" \
    "Contact distance" "Read orientation" "Full table"; do
    expect "the $section section is present" bash -c \
        "grep -q '<h2>$section</h2>' pretty/hicQC_report.html"
done
expect "the discarded-pairs definitions list is present" \
    grep -q 'Dangling ends</dt>' pretty/hicQC_report.html
expect "a download link to the underlying tab-separated table is present" \
    grep -q 'href="QC_table.txt"' pretty/hicQC_report.html
expect "the full table includes a raw QC.log column" \
    grep -q '>Sequenced reads<' pretty/hicQC_report.html
expect "the sortable-table script is present" \
    grep -q 'table.qc-table' pretty/hicQC_report.html

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
