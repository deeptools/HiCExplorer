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
