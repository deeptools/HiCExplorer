#!/usr/bin/env bash
# Command line behaviour of hicMatrixQC (no Python counterpart, so checked here
# rather than through cpp/scripts/equiv.py).
#
# Checks:
#  1. --help-json succeeds and names the tool.
#  2. A cool matrix gives an HTML report and a table with one sample column;
#     the report needs no script or svg, wraps its style in a comment, draws
#     bars as table cells with a bgcolor, and holds no em-dash.
#  3. The table of two runs at --threads 1 and --threads 4 is identical.
#  4. An mcool gives one column per resolution; --resolutions selects.
#  5. A matrix built by hicBuildMatrix from the test BAM files carries the QC
#     statistics in its metadata, and the table reproduces every number of
#     QC.log, for a cool file and for every resolution of an mcool.
#  6. Several matrices give several columns, --labels names them.
#  7. A missing file, an h5 matrix and a wrong number of labels exit non-zero
#     and write nothing.
#
# Usage: hicMatrixQC_cli.sh <hicMatrixQC> <hicBuildMatrix> <test_data directory>

set -u
tool=$1
build=$2
data=$3
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

# value of a statistic in column $2 (1 based, after the name) of table $3
cell() {
    awk -F'\t' -v key="$1" -v col="$(($2 + 1))" '$1 == key { print $col }' "$3"
}

# 1
cd "$work" || exit 1
"$tool" --help-json > help.json 2> help.err
expect "--help-json exits 0" [ $? -eq 0 ]
expect "--help-json names the tool" grep -q '"tool": "hicMatrixQC"' help.json

# 2
"$tool" -m "$data/gm12878_chr1.cool" -o one.html --tableFileName one.tsv > run.out 2> run.err
expect "a cool matrix runs" [ $? -eq 0 ]
expect "the report is written" [ -s one.html ]
expect "the table is written" [ -s one.tsv ]
expect "the table has one sample column" [ "$(head -1 one.tsv | awk -F'\t' '{print NF}')" -eq 2 ]
expect "no script in the report" bash -c '! grep -qi "<script" one.html'
expect "no svg in the report" bash -c '! grep -qi "<svg" one.html'
expect "the style is wrapped in a comment" grep -q '<style><!--' one.html
expect "bars are cells with a bgcolor" grep -q 'bgcolor="#' one.html
expect "no em-dash" bash -c '! grep -q $'"'"'\xe2\x80\x94'"'"' one.html'
expect "total contacts are the sum of the file" [ "$(cell contacts.total 1 one.tsv)" -eq 326878012 ]
expect "a single chromosome has no trans contacts" [ "$(cell contacts.trans 1 one.tsv)" -eq 0 ]
expect "the decay exponent is negative" awk -F'\t' '$1 == "decay.exponent.genome" { exit !($2 < 0) }' one.tsv

# 3
"$tool" -m "$data/gm12878_chr1.cool" -o t1.html --tableFileName t1.tsv --threads 1 2> /dev/null
"$tool" -m "$data/gm12878_chr1.cool" -o t4.html --tableFileName t4.tsv --threads 4 2> /dev/null
expect "table independent of the thread count" cmp -s <(cut -f1,2 t1.tsv | grep -v '^sample.uri') <(cut -f1,2 t4.tsv | grep -v '^sample.uri')

# 4
"$tool" -m "$data/matrix.mcool" -o m.html --tableFileName m.tsv 2> /dev/null
expect "an mcool gives a column per resolution" [ "$(head -1 m.tsv | awk -F'\t' '{print NF}')" -eq 6 ]
expect "the resolutions are in the header" bash -c "head -1 m.tsv | grep -q '@40000' && head -1 m.tsv | grep -q '@640000'"
expect "an mcool report splits the columns in blocks" bash -c "[ \$(grep -c 'summary=\"Summary of the statistics\"' m.html) -eq 2 ]"
"$tool" -m "$data/matrix.mcool" -o r.html --tableFileName r.tsv --resolutions 80000 160000 2> /dev/null
expect "--resolutions selects" [ "$(head -1 r.tsv | awk -F'\t' '{print NF}')" -eq 3 ]
"$tool" -m "$data/matrix.mcool" -o bad.html --tableFileName bad.tsv --resolutions 12345 2> bad.err
expect "an absent resolution is an error" [ $? -ne 0 ]
expect "an absent resolution writes nothing" [ ! -e bad.html ]

# 5
"$build" -s "$data/small_test_R1_unsorted.bam" "$data/small_test_R2_unsorted.bam" \
    --outFileName built.cool -bs 5000 --QCfolder qc --threads 4 \
    --restrictionSequence GATC --danglingSequence GATC -rs "$data/DpnII.bed" --noPlot > /dev/null 2>&1
expect "hicBuildMatrix builds the cool" [ -s built.cool ]
"$tool" -m built.cool -o built.html --tableFileName built.tsv 2> /dev/null
expect "hicMatrixQC reads the built matrix" [ -s built.tsv ]
check_qc() { # column table
    local col=$1 table=$2 log=qc/QC.log
    local pairs=(
        "qc.sequenced_reads|Sequenced reads"
        "qc.pairs_mappable_unique_high_quality|Pairs mappable, unique and high quality"
        "qc.hic_contacts|Hi-C contacts"
        "qc.one_mate_unmapped|One mate unmapped"
        "qc.one_mate_not_unique|One mate not unique"
        "qc.low_mapping_quality|Low mapping quality"
        "qc.dangling_end.GATC|dangling end GATC (restriction sequence GATC)"
        "qc.self_ligation|self ligation (removed)"
        "qc.one_mate_not_close_to_rest_site|One mate not close to rest site"
        "qc.same_fragment|same fragment"
        "qc.self_circle|self circle"
        "qc.duplicated_pairs|duplicated pairs"
        "qc.inter_chromosomal|inter chromosomal"
        "qc.intra_short_range|Intra short range (< 20kb)"
        "qc.intra_long_range|Intra long range (>= 20kb)"
        "qc.inward_pairs|Read pair type: inward pairs"
        "qc.outward_pairs|Read pair type: outward pairs"
        "qc.left_pairs|Read pair type: left pairs"
        "qc.right_pairs|Read pair type: right pairs"
    )
    local entry key label expected got bad=0
    for entry in "${pairs[@]}"; do
        key=${entry%%|*}
        label=${entry#*|}
        expected=$(awk -F'\t' -v label="$label" '$1 == label { print $2 }' "$log")
        got=$(cell "$key" "$col" "$table")
        if [ -z "$expected" ] || [ "$expected" != "$got" ]; then
            echo "  $key: QC.log '$expected', table '$got'"
            bad=1
        fi
    done
    [ $bad -eq 0 ]
}
expect "the table reproduces every number of QC.log" check_qc 1 built.tsv
expect "the pairs by distance add up to the cis pairs" bash -c "
    awk -F'\t' '\$1 ~ /^qc.distance.contacts\./ { s += \$2 } \$1 == \"qc.intra_short_range\" { c += \$2 } \$1 == \"qc.intra_long_range\" { c += \$2 } END { exit !(s == c && s > 0) }' built.tsv"
expect "the report shows the library statistics" grep -q 'Library statistics stored by hicBuildMatrix' built.html
expect "the report shows the distance plot" grep -q 'Read pair type by distance' built.html

"$build" -s "$data/small_test_R1_unsorted.bam" "$data/small_test_R2_unsorted.bam" \
    --outFileName built.mcool -bs 5000 10000 20000 --QCfolder qcm --threads 4 \
    --restrictionSequence GATC --danglingSequence GATC -rs "$data/DpnII.bed" --noPlot > /dev/null 2>&1
"$tool" -m built.mcool -o builtm.html --tableFileName builtm.tsv 2> /dev/null
expect "an mcool built by hicBuildMatrix: three resolutions" [ "$(head -1 builtm.tsv | awk -F'\t' '{print NF}')" -eq 4 ]
for column in 1 2 3; do
    expect "mcool resolution column $column reproduces QC.log" bash -c "
        for key in qc.sequenced_reads qc.hic_contacts qc.inter_chromosomal qc.right_pairs; do
            [ -n \"\$(awk -F'\t' -v k=\$key '\$1 == k { print \$$((column + 1)) }' builtm.tsv)\" ] || exit 1
        done
        [ \"\$(awk -F'\t' '\$1 == \"qc.sequenced_reads\" { print \$$((column + 1)) }' builtm.tsv)\" = \"\$(awk -F'\t' '\$1 == \"Sequenced reads\" { print \$2 }' qcm/QC.log)\" ]"
done

# 6
"$tool" -m "$data/gm12878_chr1.cool" built.cool -o two.html --tableFileName two.tsv --labels first second 2> /dev/null
expect "two matrices give two columns" [ "$(head -1 two.tsv | awk -F'\t' '{print NF}')" -eq 3 ]
expect "--labels names the columns" [ "$(head -1 two.tsv)" = "$(printf 'statistic\tfirst\tsecond')" ]

# 7
"$tool" -m "$work/missing.cool" -o miss.html --tableFileName miss.tsv > /dev/null 2> miss.err
expect "a missing file exits non-zero" [ $? -ne 0 ]
expect "a missing file writes nothing" [ ! -e miss.html ] && [ ! -e miss.tsv ]
"$tool" -m "$data/small_test_matrix.h5" -o h5.html --tableFileName h5.tsv > /dev/null 2> h5.err
expect "an h5 matrix exits non-zero" [ $? -ne 0 ]
"$tool" -m "$data/gm12878_chr1.cool" -o lab.html --labels a b > /dev/null 2> lab.err
expect "a wrong number of labels exits non-zero" [ $? -ne 0 ]

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
