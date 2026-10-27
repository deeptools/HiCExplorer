#!/usr/bin/env bash
# The `hicexplorer` command (cpp/tools/hicexplorer.cpp):
#
#  1. without arguments it prints the same help as --help and exits 0;
#  2. --version prints "hicexplorer <version>";
#  3. the list names every executable of the tools directory but itself, each
#     once, and the Python drivers, in sorted order;
#  4. an unknown option is an argparse error, exit status 2.
#
# Usage: hicexplorer_cli.sh <hicexplorer binary>

set -u
tool=$1
tools_dir=$(dirname "$tool")
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

"$tool" > "$work/bare.txt" 2> "$work/bare.err"
bare_status=$?
"$tool" --help > "$work/help.txt" 2> "$work/help.err"
expect "no arguments exits 0" [ "$bare_status" -eq 0 ]
expect "no arguments prints the help" cmp -s "$work/bare.txt" "$work/help.txt"
expect "the help starts with the usage" [ "$(head -n 1 "$work/help.txt")" = "usage: hicexplorer [-h] [--version]" ]

"$tool" --version > "$work/version.txt"
expect "--version names the command" grep -q -E '^hicexplorer [0-9]' "$work/version.txt"

sed -n '/^The following is the list of tools:$/,/^For more information/p' "$work/help.txt" |
    sed '1,2d;$d' | sed '$d' | cut -d' ' -f1 > "$work/listed.txt"
for path in "$tools_dir"/hic* "$tools_dir"/chic*; do
    name=$(basename "$path")
    if [ -f "$path" ] && [ -x "$path" ] && [ "$name" != hicexplorer ]; then
        echo "$name"
    fi
done | sort > "$work/expected_exe.txt"
{ cat "$work/expected_exe.txt"; printf '%s\n' hicHyperoptDetectLoops hicHyperoptDetectLoopsHiCCUPS; } |
    LC_ALL=C sort > "$work/expected.txt"
expect "the list is the executables and the Python drivers, sorted" cmp -s "$work/listed.txt" "$work/expected.txt"
if ! cmp -s "$work/listed.txt" "$work/expected.txt"; then
    diff "$work/expected.txt" "$work/listed.txt"
fi

"$tool" --no-such-option > /dev/null 2> "$work/error.txt"
expect "an unknown option exits 2" [ $? -eq 2 ]
expect "an unknown option is reported as argparse does" \
    grep -q -x 'hicexplorer: error: unrecognized arguments: --no-such-option' "$work/error.txt"

[ "$failures" -eq 0 ]
