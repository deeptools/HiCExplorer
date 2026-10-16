#!/usr/bin/env bash
# hicPrepareQCreport --prettyReport (C++-only, no Python counterpart): a
# second, card-based HTML report at OUTPUTFOLDER/hicQC_report.html, judging
# five metrics per sample as PASS, WARN or FAIL against --min*/--max*Percent
# cutoffs, independent of and unchanged from the existing hicQC.html the
# tool already draws through the Python plotting layer. Checked here rather
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

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
