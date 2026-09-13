#!/usr/bin/env bash
# Command line behaviour of the C++ hicMergeDomains around graphviz rendering,
# which the harness does not see:
#
#  1. With graphviz's dot on PATH, one <prefix>_<chromosome>.<format> is
#     rendered per chromosome with a valid signature, the DOT source that
#     graphviz's render(cleanup=True) deletes is deleted too, a missing prefix
#     directory is created as graphviz does, and "Saved relation tree of" is
#     printed once per graph.
#  2. Without dot on PATH, or with an unknown --outputTreePlotFormat, the tool
#     exits non-zero before writing anything (contract rule 7: never exit 0
#     without a requested file).
#  3. A dot that fails makes the tool fail, and the source stays on disk, as
#     graphviz raises CalledProcessError before its cleanup.
#  4. One domain file without a protein file exits 1, as the Python does.
#
# Usage: hicMergeDomains_cli.sh <hicMergeDomains binary> <test_data directory>
#                               <conda prefix holding bin/dot>

set -u
tool=$1
data=$2/hicMergeDomains
deps=$3
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

only_logs_in() {
    [ -z "$(ls -A "$1" | grep -v -x -e out.txt -e err.txt)" ]
}

dot_dir=""
if command -v dot > /dev/null 2>&1; then
    dot_dir=$(dirname "$(command -v dot)")
elif [ -x "$deps/bin/dot" ]; then
    dot_dir="$deps/bin"
fi

# 1
if [ -n "$dot_dir" ]; then
    mkdir "$work/render" && cd "$work/render" || exit 1
    PATH="$dot_dir:$PATH" "$tool" -d "$data/10kbtad_domains.bed" "$data/50kbtad_domains.bed" \
        -om merged.bed -or relations.txt -ot trees/tree -of png > out.txt 2> err.txt
    status=$?
    expect "rendering exits 0" [ "$status" -eq 0 ]
    expect "23 trees are rendered" [ "$(ls trees | grep -c '^tree_.*\.png$')" -eq 23 ]
    expect "nothing but the rendered trees is left in the prefix directory" \
        [ "$(ls trees | grep -vc '\.png$')" -eq 0 ]
    expect "the tree of chromosome X is a png" \
        bash -c "head -c 8 trees/tree_X.png | od -An -c | grep -q 'P   N   G'"
    expect "one 'Saved relation tree of' line per graph" \
        [ "$(grep -c '^Saved relation tree of ' out.txt)" -eq 23 ]
    expect "the merged list is written" [ -s merged.bed ]
    expect "the relation list is written" [ -s relations.txt ]
else
    echo "FAILED: no dot executable on PATH or in $deps/bin"
    failures=$((failures + 1))
fi

# 2
mkdir "$work/nodot" && cd "$work/nodot" || exit 1
if PATH=/usr/bin:/bin command -v dot > /dev/null 2>&1; then
    echo "skipped: /usr/bin or /bin holds a dot, so its absence cannot be tested"
else
    PATH=/usr/bin:/bin "$tool" -d "$data/10kbtad_domains.bed" "$data/50kbtad_domains.bed" \
        -om merged.bed -or relations.txt -ot tree -of png > out.txt 2> err.txt
    status=$?
    expect "a missing dot exits non-zero" [ "$status" -ne 0 ]
    expect "a missing dot is named" grep -q "dot executable was not found" err.txt
    expect "a missing dot writes nothing" only_logs_in "$work/nodot"
fi
PATH="${dot_dir:-/usr/bin}:$PATH" "$tool" -d "$data/10kbtad_domains.bed" \
    "$data/50kbtad_domains.bed" -om merged.bed -or relations.txt -ot tree \
    -of jpeg2000 > out.txt 2> err.txt
status=$?
expect "an unknown format exits non-zero" [ "$status" -ne 0 ]
expect "an unknown format writes nothing" only_logs_in "$work/nodot"

# 3
mkdir -p "$work/faildot/bin" && cd "$work/faildot" || exit 1
printf '#!/bin/sh\necho "dot: simulated failure" >&2\nexit 1\n' > bin/dot
chmod +x bin/dot
PATH="$work/faildot/bin:$PATH" "$tool" -d "$data/10kbtad_domains.bed" \
    "$data/50kbtad_domains.bed" -om merged.bed -or relations.txt -ot tree -of png \
    > out.txt 2> err.txt
status=$?
expect "a failing dot makes the tool fail" [ "$status" -ne 0 ]
expect "the source of the failed rendering stays, as graphviz leaves it" [ -s tree_1 ]
expect "no rendered file is claimed" [ ! -e tree_1.png ]

# 4
mkdir "$work/single" && cd "$work/single" || exit 1
"$tool" -d "$data/10kbtad_domains.bed" -om merged.bed > out.txt 2> err.txt
status=$?
expect "one domain file without a protein file exits 1" [ "$status" -eq 1 ]
expect "one domain file without a protein file writes nothing" only_logs_in "$work/single"

echo "$failures failure(s)"
[ "$failures" -eq 0 ]
