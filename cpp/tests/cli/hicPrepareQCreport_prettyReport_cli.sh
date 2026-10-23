#!/usr/bin/env bash
# hicPrepareQCreport --prettyReport (C++-only, no Python counterpart): a
# second HTML report at OUTPUTFOLDER/hicQC_report.html, independent of and
# unchanged from the existing hicQC.html the tool already draws through the
# Python plotting layer. It follows the falcoAmadeus report design: a header
# with a light/dark theme toggle and a plot palette picker, a summary sidebar
# with the worst verdict of every section, and one collapsible card per
# section (overview with five metrics per sample judged PASS/WARN/FAIL against
# --min*/--max*Percent cutoffs; pairs sequenced; unmappable/non-unique; pairs
# discarded with its definitions list; contact distance; read orientation;
# the full QC_table.txt). Each section draws grouped Plotly bars and keeps its
# numbers in plain HTML tables, which link to the same underlying .txt tables
# write_qc_report_tables writes. Checked here rather than through
# cpp/scripts/equiv.py since it has no reference to match.
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

report=pretty/hicQC_report.html
for entry in "overview:Overview" "sequenced:Pairs sequenced" \
    "unmappable:Unmappable and non-unique pairs" "discarded:Pairs discarded" \
    "distance:Contact distance" "orientation:Read orientation" "fulltable:Full table"; do
    id=${entry%%:*}
    section=${entry#*:}
    expect "the $section card is present, titled with its verdict" bash -c \
        "grep -Eq '<h2 class=\"(pass|warn|fail|info)\" id=\"$id\"><span class=\"title\">$section: (pass|warn|fail|not judged)</span></h2>' $report"
    expect "the sidebar lists $section with a status class" bash -c \
        "grep -Eq '<li><a class=\"(pass|warn|fail|info)\" href=\"#$id\"' $report"
done
expect "there is exactly one card per section (7)" \
    test "$(grep -c '<div class="module" id="module-' $report)" -eq 7
expect "the overview card carries the worst verdict (FAIL for the default cutoffs on this log)" \
    grep -q '<h2 class="fail" id="overview"><span class="title">Overview: fail' $report
expect "the lenient run's contact distance card passes" \
    grep -q '<h2 class="pass" id="distance"><span class="title">Contact distance: pass' lenient/hicQC_report.html
expect "the discarded-pairs definitions list is present" \
    grep -q 'Dangling ends</dt>' $report
expect "a download link to the underlying tab-separated table is present" \
    grep -q 'href="QC_table.txt"' $report
for table in unmapable_table.txt discarded_table.txt distance_table.txt read_orientation_table.txt; do
    expect "a download link to $table is present" grep -q "href=\"$table\"" $report
done
expect "the full table includes a raw QC.log column" \
    grep -q '>Sequenced reads<' $report
for table in sequenced-table unmappable-table discarded-table distance-table orientation-table full-table; do
    expect "the $table data table is plain HTML" grep -q "<table class=\"qc-table\" id=\"$table\">" $report
done
for plot in sequenced unmappable discarded distance orientation; do
    expect "the plot-$plot Plotly div and its chart data are present" bash -c \
        "grep -q '<div class=\"qc-plot\" id=\"plot-$plot\"></div>' $report && grep -q '{\"id\":\"plot-$plot\"' $report"
done
expect "Plotly is loaded from a pinned version" \
    grep -q 'src="https://cdn.plot.ly/plotly-2.35.2.min.js"' $report
expect "the theme toggle, its stored key and the ?theme= override are present" bash -c \
    "grep -q 'id=\"theme-toggle\"' $report && grep -q \"'hicx-theme'\" $report && grep -q 'theme=(light|dark)' $report"
expect "the palette picker offers Ocean, Sunset and Ink" bash -c \
    "grep -q 'id=\"palette-select\"' $report && grep -q '>Ocean<' $report && grep -q '>Sunset<' $report && grep -q '>Ink<' $report"
expect "the sidebar collapse state is stored" grep -q "'hicx-nav'" $report
expect "the header names the report and the number of samples" bash -c \
    "grep -q '<div id=\"header_title\">HiCExplorer QC</div>' $report && grep -q '<br>1 sample</div>' $report"
expect "the footer names the HiCExplorer version" grep -q '<div class="footer">HiCExplorer [0-9]' $report
expect "a link to the HiCExplorer documentation is present" \
    grep -q 'href="https://hicexplorer.readthedocs.io/' $report
expect "the sortable-table script is present" \
    grep -q 'table.qc-table' $report
expect "the report contains no em dash" bash -c "! grep -q $'\xe2\x80\x94' $report"

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
