#!/usr/bin/env bash
# chicChicagoPlotViewpoint (no Python counterpart): the data the figure is
# drawn from is checked through --plotData, so no matplotlib is needed. Scores
# come from the ported chicChicago tools on the chr20/chr21 GM12878 fixture.
#
# Usage: chicChicagoPlotViewpoint_cli.sh <plot tool> <background model tool>
#                                        <scores tool> <test_data directory>

set -u
tool=$1
model_tool=$2
scores_tool=$3
data=$4/chicago
p=$data/h19_chr20and21
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
"$model_tool" --rmap "$p.rmap" --baitmap "$p.baitmap" --chinput "$data/GM_rep1.chinput" \
    --nperbin "$p.npb" --nbaitsperbin "$p.nbpb" --proxOE "$p.poe" -o bg.txt >/dev/null 2>&1
"$scores_tool" --chinput "$data/GM_rep1.chinput" --backgroundModel bg.txt \
    --baitmap "$p.baitmap" --rmap "$p.rmap" -o scores.txt >/dev/null 2>&1
expect "scores were produced" test -s scores.txt

json() { python3 -c "import json,sys; d=json.load(open('$1')); $2"; }

# One bait, scatter, with the Brownian overlay.
expect "scatter run exits 0" "$tool" --scores scores.txt --baitID 417632 --baitmap "$p.baitmap" \
    --backgroundModel bg.txt --plotData scatter.json
expect "scatter holds 168 points for RBM38 with matching x, y, score" json scatter.json \
    "assert d['style']=='scatter' and len(d['x'])==len(d['y'])==len(d['score'])==168 and d['title'].startswith('RBM38')"
expect "the background model adds a mean and an upper band" json scatter.json \
    "assert d['hasBackground'] and len(d['bmeanX'])==len(d['bmeanY'])==len(d['bmeanUpperY'])>0 and all(u>=m for u,m in zip(d['bmeanUpperY'],d['bmeanY']))"
"$tool" --scores scores.txt --baitID 417632 --baitmap "$p.baitmap" --plotData scatter2.json
expect "without --backgroundModel there is no overlay" json scatter2.json "assert not d['hasBackground']"

# One bait as arcs: only significant rows.
"$tool" --scores scores.txt --baitID 417632 --baitmap "$p.baitmap" --style arcs --onlySignificant \
    --plotData arcs1.json
expect "--onlySignificant keeps only scores of at least plevel2" json arcs1.json \
    "assert len(d['score'])>0 and min(d['score'])>=3 and d['anchors']==[0.0]"
"$tool" --scores scores.txt --baitID 417632 --baitmap "$p.baitmap" --style arcs --plevel2 5 \
    --onlySignificant --plotData arcs5.json
expect "raising --plevel2 keeps fewer rows" json arcs5.json \
    "import json as j; n=len(j.load(open('arcs1.json'))['score']); assert min(d['score'])>=5 and len(d['score'])<n"

# Region: all baits, links file.
"$tool" --scores scores.txt --region 20 0 3000000 --baitmap "$p.baitmap" --rmap "$p.rmap" \
    --style arcs --onlySignificant --linksFile chicago.links --plotData region.json
expect "region arcs hold 220 significant interactions" json region.json \
    "assert len(d['score'])==len(d['arcX1'])==len(d['arcX2'])==len(d['arcHeight'])==220 and min(d['score'])>=3"
expect "the links file has 220 rows of 7 fields" bash -c \
    "test \$(wc -l < chicago.links) -eq 220 && test \$(awk -F'\t' 'NF!=7' chicago.links | wc -l) -eq 0"
expect "link ends lie on the region chromosome and start before they end" bash -c \
    "awk -F'\t' '\$1!=\"20\" || \$4!=\"20\" || \$2>=\$3 || \$5>=\$6 {bad=1} END{exit bad}' chicago.links"
"$tool" --scores scores.txt --region 20 0 3000000 --baitmap "$p.baitmap" --style arcs \
    --plotData region_all.json
expect "without --onlySignificant the region holds more rows" json region_all.json \
    "import json as j; assert len(d['score'])>len(j.load(open('region.json'))['score'])"
"$tool" --scores scores.txt --region 20 0 3000000 --baitmap "$p.baitmap" --rmap "$p.rmap" \
    --style arcs --onlySignificant --linksFile again.links --plotData again.json
expect "repeated runs are identical" cmp -s chicago.links again.links

# Refusals leave no output behind.
expect "scatter with --region is refused" bash -c \
    "! '$tool' --scores scores.txt --region 20 0 3000000 --baitmap '$p.baitmap' --plotData x.json 2>/dev/null"
expect "--region without --baitmap is refused" bash -c \
    "! '$tool' --scores scores.txt --region 20 0 3000000 --style arcs --plotData x.json 2>/dev/null"
expect "--linksFile without --baitmap is refused" bash -c \
    "! '$tool' --scores scores.txt --baitID 417632 --linksFile nolinks.txt --plotData x.json 2>/dev/null"
expect "--baitID and --region together are refused" bash -c \
    "! '$tool' --scores scores.txt --baitID 417632 --region 20 0 10 --plotData x.json 2>/dev/null"
expect "the refusals wrote nothing" test ! -e x.json -a ! -e nolinks.txt

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
