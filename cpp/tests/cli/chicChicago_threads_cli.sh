#!/usr/bin/env bash
# chicChicagoBackgroundModel, chicChicagoScores and
# chicChicagoSignificantInteractions give byte-identical output for every
# --threads value on an input file under 1 MiB.
#
# Their chunked readers (read_chinput, read_scores) parse a file under 1 MiB
# with one thread. They used to compute the chunk bounds for the requested
# thread count first, so that one thread read only the first of `threads`
# chunks and the rest of the file was lost (half the rows at --threads 2).
# Every other test input is larger than 1 MiB and did not reach that path.
#
# It also checks that chicChicagoBackgroundModel ends with a message, not a
# crash, when no bait has enough reads to fit a model.
#
# Usage: chicChicago_threads_cli.sh <directory of the built tools> <test_data directory>

set -u
tools=$1
data=$2/chicago
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

cd "$work" || exit 1
# About 600 kB of GM_rep1.chinput: its comment line, the header and complete
# rows only.
head -c 600000 "$data/GM_rep1.chinput" | sed '$d' >small.chinput
design=(--baitmap "$data/h19_chr20and21.baitmap" --rmap "$data/h19_chr20and21.rmap")
rows=$(($(wc -l <small.chinput) - 2))

for t in 1 4; do
    "$tools/chicChicagoBackgroundModel" --chinput small.chinput "${design[@]}" \
        --nperbin "$data/h19_chr20and21.npb" --nbaitsperbin "$data/h19_chr20and21.nbpb" \
        --proxOE "$data/h19_chr20and21.poe" --threads "$t" -o "bg$t.txt" >/dev/null 2>&1
    expect "chicChicagoBackgroundModel --threads $t exits 0" test $? -eq 0
    "$tools/chicChicagoScores" --chinput small.chinput --backgroundModel bg1.txt "${design[@]}" \
        --threads "$t" -o "scores$t.txt" >/dev/null 2>&1
    expect "chicChicagoScores --threads $t exits 0" test $? -eq 0
    "$tools/chicChicagoSignificantInteractions" --scores scores1.txt "${design[@]}" \
        --threads "$t" -o "significant$t.txt" >/dev/null 2>&1
    expect "chicChicagoSignificantInteractions --threads $t exits 0" test $? -eq 0
done

expect "the background model does not depend on --threads" cmp -s bg1.txt bg4.txt
expect "the scores do not depend on --threads" cmp -s scores1.txt scores4.txt
expect "the significant interactions do not depend on --threads" \
    cmp -s significant1.txt significant4.txt
# The filters drop rows, but most of the file's rows are scored: a reader
# that lost a chunk scores far fewer.
scored=$(($(wc -l <scores4.txt) - 1))
expect "most chinput rows are scored ($scored of $rows)" test "$scored" -gt $((rows / 2))

# Too few reads per bait: every bait falls below --minNPerBait, nothing is
# left to fit. This used to crash with a segmentation fault (an empty sample
# reached the quantile of add_tlb); it is an error message and exit 1 now.
head -n 200 small.chinput >tiny.chinput
"$tools/chicChicagoBackgroundModel" --chinput tiny.chinput "${design[@]}" \
    --nperbin "$data/h19_chr20and21.npb" --nbaitsperbin "$data/h19_chr20and21.nbpb" \
    --proxOE "$data/h19_chr20and21.poe" -o tiny_bg.txt >tiny.out 2>tiny.err
expect "an input with too few reads per bait exits 1" test $? -eq 1
expect "the message names --minNPerBait" grep -q "minNPerBait" tiny.err

echo "$failures failing checks"
exit "$failures"
