#!/usr/bin/env bash
# The drawing environment check of the C++ plotting tools
# (plot/hicexplorer_plot/environment.py, hicx::plot::preflight):
#
#  1. A wrong matplotlib version, faked by a stub package on PYTHONPATH, is
#     refused with exit status 3 before anything is read or written; the
#     message names the found and the required version, HICX_PLOT_PYTHON and
#     the opt-out.
#  2. A missing matplotlib is refused the same way, also with the opt-out set.
#  3. A pygenometracks that is not the installed pyGenomeTracks 3.9 is refused
#     for hicPlotTADs.
#  4. HICX_PLOT_ALLOW_UNPINNED=1 accepts a wrong version with a warning.
#  5. --plotData draws nothing and therefore checks nothing.
#  5b. For every drawing tool that has FileType('w') outputs or creates an
#     output folder (hicPlotDistVsCounts, hicPlotMatrix, hicAggregateContacts,
#     hicBuildMatrix, hicBuildMatrixMicroC, hicQuickQC), a refusal leaves the
#     output directory
#     byte for byte as it was: an existing output keeps its content and no new
#     file or folder appears.
#  6. With the pinned interpreter the figure is drawn without a warning, and
#     the pins equal those of plot/pyproject.toml. Needs the pinned
#     interpreter (CMake HICX_PLOT_TEST_PYTHON); without it this part is
#     reported and the test ends with the skip status 77.
#
# Usage: plot_environment_cli.sh <tools directory> <test_data directory>
#                                <pinned python or ""> <plot/pyproject.toml>

set -u
tools=$1
data=$2
pinned=$3
pyproject=$4
python=${pinned:-python3}
package=$(cd "$tools/../plot" && pwd)
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

absent() {
    local path
    for path in "$@"; do
        [ ! -e "$path" ] || return 1
    done
}

stub() {
    mkdir -p "$work/stub_$1/$2"
    printf '%s\n' "$3" > "$work/stub_$1/$2/__init__.py"
    echo "$work/stub_$1"
}

wrong=$(stub wrong matplotlib '__version__ = "3.11.0"')
missing=$(stub missing matplotlib 'raise ModuleNotFoundError("No module named '"'"'matplotlib'"'"'", name="matplotlib")')
foreign_pgt=$(stub pgt pygenometracks '')
iit=(-m "$data/hicInterIntraTAD/GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool"
     --tadDomains "$data/hicInterIntraTAD/untreated_R1_domains_chr1_chr2.bed" -t 1)

# 1
mkdir "$work/wrong" && cd "$work/wrong" || exit 1
HICX_PLOT_PYTHON=$python PYTHONPATH=$wrong "$tools/hicInterIntraTAD" "${iit[@]}" \
    -o table.txt -op ratio.png > out.txt 2> err.txt
status=$?
expect "a wrong matplotlib exits 3" [ "$status" -eq 3 ]
expect "the message names both versions" \
    grep -q "matplotlib 3.11.0 is installed, but matplotlib 3.8.4 is required" err.txt
expect "the message names HICX_PLOT_PYTHON" grep -q "HICX_PLOT_PYTHON" err.txt
expect "the message names the opt-out" grep -q "HICX_PLOT_ALLOW_UNPINNED=1" err.txt
expect "the tool says nothing was written" grep -q "nothing was read or written" err.txt
expect "neither the table nor the figure is written" absent table.txt ratio.png

# 2
mkdir "$work/missing" && cd "$work/missing" || exit 1
HICX_PLOT_PYTHON=$python PYTHONPATH=$missing "$tools/hicInterIntraTAD" "${iit[@]}" \
    -o table.txt -op ratio.png > out.txt 2> err.txt
status=$?
expect "a missing matplotlib exits 3" [ "$status" -eq 3 ]
expect "the message says it is not installed" \
    grep -q "matplotlib is not installed .*(required: matplotlib 3.8.4)" err.txt
expect "a missing package writes nothing" absent table.txt ratio.png
HICX_PLOT_ALLOW_UNPINNED=1 HICX_PLOT_PYTHON=$python PYTHONPATH=$missing \
    "$tools/hicInterIntraTAD" "${iit[@]}" -o table.txt -op ratio.png > out.txt 2> err.txt
status=$?
expect "the opt-out does not accept a missing package" [ "$status" -eq 3 ]
expect "and still writes nothing" absent table.txt ratio.png

# 3
mkdir "$work/pgt" && cd "$work/pgt" || exit 1
printf '[x-axis]\n' > tracks.ini
HICX_PLOT_PYTHON=$python PYTHONPATH=$foreign_pgt "$tools/hicPlotTADs" --tracks tracks.ini \
    --region chr1:1000-2000 --outFileName tads.png > out.txt 2> err.txt
status=$?
expect "a foreign pygenometracks exits 3 for hicPlotTADs" [ "$status" -eq 3 ]
expect "the message names pyGenomeTracks 3.9" grep -q "required: pyGenomeTracks 3.9" err.txt
expect "no figure is written" absent tads.png

# 4
cd "$work" || exit 1
HICX_PLOT_ALLOW_UNPINNED=1 PYTHONPATH=$package:$wrong "$python" -m hicexplorer_plot \
    --check hicInterIntraTAD > out.txt 2> err.txt
status=$?
expect "the opt-out accepts a wrong version" [ "$status" -eq 0 ]
expect "and warns" grep -q "warning: hicInterIntraTAD draws with matplotlib 3.11.0" err.txt

# 5
mkdir "$work/plotdata" && cd "$work/plotdata" || exit 1
HICX_PLOT_PYTHON=$python PYTHONPATH=$wrong "$tools/hicInterIntraTAD" "${iit[@]}" \
    -o table.txt -op ratio.png --plotData data.json > out.txt 2> err.txt
status=$?
expect "--plotData runs without the check" [ "$status" -eq 0 ]
expect "--plotData writes the table and the data" [ -s table.txt ] && [ -s data.json ]

# 5b
snapshot() {
    (cd "$1" && find . -printf '%y %p %s\n' | sort && find . -type f -exec sha256sum {} + | sort)
}

unchanged_after_refusal() {
    local name=$1
    shift
    local dir="$work/untouched_$name"
    "$@" > "$work/untouched_$name.out" 2> "$work/untouched_$name.err"
    local status=$?
    expect "$name: a refusal exits 3" [ "$status" -eq 3 ]
    expect "$name: the refusal names the versions" \
        grep -q "matplotlib 3.11.0 is installed, but matplotlib 3.8.4 is required" \
        "$work/untouched_$name.err"
    snapshot "$dir" > "$work/untouched_$name.after"
    expect "$name: the output directory is byte for byte unchanged" \
        cmp -s "$work/untouched_$name.before" "$work/untouched_$name.after"
}

prepare() {
    local dir="$work/untouched_$1"
    shift
    mkdir -p "$dir"
    local file
    for file in "$@"; do
        mkdir -p "$(dirname "$dir/$file")"
        printf 'existing content of %s\n' "$file" > "$dir/$file"
    done
    snapshot "$dir" > "$work/untouched_$(basename "$dir" | sed 's/^untouched_//').before"
    echo "$dir"
}

refused() {
    HICX_PLOT_PYTHON=$python PYTHONPATH=$wrong "$@"
}

d=$(prepare hicPlotDistVsCounts dist_vs_counts.png)
unchanged_after_refusal hicPlotDistVsCounts refused "$tools/hicPlotDistVsCounts" \
    --matrices "$data/small_test_matrix_50kb_res.h5" --plotFile "$d/dist_vs_counts.png" \
    --outFileData "$d/data.txt" --plotsize 8 4

d=$(prepare hicAggregateContacts aggregate.png m_genome.tab)
unchanged_after_refusal hicAggregateContacts refused "$tools/hicAggregateContacts" \
    --matrix "$data/Li_et_al_2015.h5" --BED "$data/hicAggregateContacts/test_regions.bed" \
    --mode intra-chr --range 50000:900000 --numberOfBins 30 --outFileName "$d/aggregate.png" \
    --outFilePrefixMatrix "$d/m" --outFileContactPairs "$d/p" \
    --diagnosticHeatmapFile "$d/heatmap.png"

d=$(prepare hicBuildMatrix matrix.h5 qc/QC.log)
unchanged_after_refusal hicBuildMatrix refused "$tools/hicBuildMatrix" \
    -s "$data/R1_1000.bam" "$data/R2_1000.bam" --outFileName "$d/matrix.h5" -bs 100000 \
    --QCfolder "$d/qc" --restrictionSequence AAGCTT --danglingSequence AGCT \
    -rs "$data/hicFindRestSite/hindIII.bed" --outBam "$d/valid.bam"

d=$(prepare hicBuildMatrixMicroC valid.bam)
unchanged_after_refusal hicBuildMatrixMicroC refused "$tools/hicBuildMatrixMicroC" \
    -s "$data/R1_1000.bam" "$data/R2_1000.bam" --outFileName "$d/matrix.h5" \
    --QCfolder "$d/qc" --binSize 100000 --outBam "$d/valid.bam"

d=$(prepare hicPlotMatrix matrix.png)
unchanged_after_refusal hicPlotMatrix refused "$tools/hicPlotMatrix" \
    --matrix "$data/small_test_matrix_50kb_res.h5" --outFileName "$d/matrix.png" --log1p
d=$(prepare hicPlotMatrix_new)
unchanged_after_refusal hicPlotMatrix_new refused "$tools/hicPlotMatrix" \
    --matrix "$data/small_test_matrix_50kb_res.h5" --outFileName "$d/matrix.png"

d=$(prepare hicQuickQC)
unchanged_after_refusal hicQuickQC refused "$tools/hicQuickQC" \
    -s "$data/R1_1000.bam" "$data/R2_1000.bam" --QCfolder "$d/qc" -seq AAGCTT \
    --danglingSequence AGCT -rs "$data/hicFindRestSite/hindIII.bed" --lines 500

# 6
if [ -z "$pinned" ]; then
    echo "skipped: the pinned interpreter part (HICX_PLOT_TEST_PYTHON is not configured)"
    [ "$failures" -eq 0 ] && exit 77
else
    mkdir "$work/pinned" && cd "$work/pinned" || exit 1
    HICX_PLOT_PYTHON=$pinned "$tools/hicInterIntraTAD" "${iit[@]}" -o table.txt -op ratio.png \
        > out.txt 2> err.txt
    status=$?
    expect "the pinned interpreter draws" [ "$status" -eq 0 ]
    expect "the table and the figure are written" [ -s table.txt ] && [ -s ratio.png ]
    expect "without a warning" bash -c '! grep -q "warning" err.txt'
    PYTHONPATH=$package "$pinned" - "$pyproject" <<'PY' > pins.txt 2>&1
import sys
import tomllib
from hicexplorer_plot import environment
with open(sys.argv[1], "rb") as handle:
    pins = {d.split("==")[0]: d.split("==")[1] for d in tomllib.load(handle)["project"]["dependencies"] if "==" in d}
wanted = list(environment.PINS.values()) + [p for tool in environment.TOOL_PINS.values() for p in tool.values()]
bad = [(name, version, pins.get(name)) for name, version in wanted if pins.get(name) != version]
print(bad)
sys.exit(1 if bad else 0)
PY
    expect "the pins equal plot/pyproject.toml" [ "$?" -eq 0 ]
fi

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
