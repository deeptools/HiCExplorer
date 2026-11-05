#!/usr/bin/env bash
# hicBuildMatrix --pairsFile at the command line (cpp/PLAN.md tier 9, item
# 9.2). The equivalence harness compares its matrices with cooler and with the
# BAM route; what only the command line shows is checked here:
#
#  1. The options a .pairs file cannot serve are refused with exit status 2
#     and argparse's message, and nothing is created: --outBam,
#     --restrictionSequence, --danglingSequence, --keepSelfLigation and
#     --keepSelfCircles; --restrictionCutFile, --minDistance, --maxDistance and
#     --maxLibraryInsertSize together with --binSize; neither --binSize nor
#     --restrictionCutFile; --pairsFile together with --samFiles. A BAM run
#     without --restrictionCutFile, --restrictionSequence and
#     --danglingSequence is refused the same way.
#  2. --minMappingQuality on a file without mapq columns exits 1 and creates
#     no QC folder.
#  3. The duplicate check in file order (a file declared sorted and upper
#     triangle) and with a hash set (the same pairs shuffled, without those
#     header lines) give the same h5 matrix and the same QC counts, at
#     --threads 1 and 4 and with --inputBufferSize 997.
#  4. A file declared sorted whose pairs are not sorted is refused with exit
#     status 1 and no matrix.
#  5. The QC counters of R1_1000_all.pairs at --minMappingQuality 30 equal the
#     counts awk takes from the file.
#
# Parts 3 to 5 draw the QC report and need the pinned drawing interpreter;
# without it they are reported and the test ends with the skip status 77.
#
# Usage: hicBuildMatrix_pairs_cli.sh <hicBuildMatrix> <test_data directory> <pinned python or "">

set -u
tool=$1
data=$2
pinned=${3:-}
pairs=$data/hicBuildMatrix/pairs
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
failures=0

fail() {
    echo "FAIL: $*"
    failures=$((failures + 1))
}

# expect_refusal <status> <message part> <arguments...>: the run exits with
# <status>, stderr holds <message part>, and the output directory stays empty.
expect_refusal() {
    local status=$1 message=$2
    shift 2
    local out=$work/refusal
    rm -rf "$out"
    mkdir -p "$out"
    (cd "$out" && "$tool" "$@" -o "$out/matrix.h5" --QCfolder "$out/qc" > "$work/stdout" 2> "$work/stderr")
    local code=$?
    if [ "$code" != "$status" ]; then
        fail "exit $code instead of $status for: $*"
        sed 's/^/    /' "$work/stderr"
    elif ! grep -qF -- "$message" "$work/stderr"; then
        fail "no '$message' in the message for: $*"
        sed 's/^/    /' "$work/stderr"
    elif [ -n "$(ls -A "$out")" ]; then
        fail "the refused run created $(ls -A "$out") for: $*"
    fi
}

valid=$pairs/small_test_valid_5kb.pairs.gz
reserved=$pairs/small_test_valid_rf.pairs.gz

# 1
refused="not allowed with argument --pairsFile"
expect_refusal 2 "argument --outBam/-b: $refused" --pairsFile "$valid" -bs 5000 --outBam "$work/refusal/x.bam"
expect_refusal 2 "argument --restrictionSequence/-seq: $refused" --pairsFile "$valid" -bs 5000 --restrictionSequence GATC
expect_refusal 2 "argument --danglingSequence: $refused" --pairsFile "$valid" -bs 5000 --danglingSequence GATC
expect_refusal 2 "argument --keepSelfLigation: $refused" --pairsFile "$valid" -bs 5000 --keepSelfLigation
expect_refusal 2 "argument --keepSelfCircles: $refused" --pairsFile "$valid" -bs 5000 --keepSelfCircles
with_bins="not allowed with arguments --pairsFile and --binSize"
expect_refusal 2 "argument --restrictionCutFile/-rs: $with_bins" --pairsFile "$valid" -bs 5000 -rs "$data/DpnII.bed"
expect_refusal 2 "argument --minDistance: $with_bins" --pairsFile "$valid" -bs 5000 --minDistance 100
expect_refusal 2 "argument --maxDistance: $with_bins" --pairsFile "$valid" -bs 5000 --maxDistance 100
expect_refusal 2 "argument --maxLibraryInsertSize: $with_bins" --pairsFile "$valid" -bs 5000 --maxLibraryInsertSize 100
expect_refusal 2 "needs --binSize/-bs for fixed bins or --restrictionCutFile/-rs" --pairsFile "$valid"
expect_refusal 2 "argument --samFiles/-s: not allowed with argument --pairsFile" --pairsFile "$valid" -s "$data/R1_1000.bam" "$data/R2_1000.bam" -bs 5000
expect_refusal 2 "the following arguments are required: --restrictionCutFile/-rs, --restrictionSequence/-seq, --danglingSequence" -s "$data/R1_1000.bam" "$data/R2_1000.bam" -bs 5000
expect_refusal 2 "one of the arguments --samFiles/-s --pairsFile is required" -bs 5000

if [ -z "$pinned" ]; then
    echo "SKIP: parts 2 to 5 need the pinned drawing interpreter (HICX_PLOT_TEST_PYTHON)"
    [ "$failures" = 0 ] && exit 77
    exit 1
fi
export HICX_PLOT_PYTHON=$pinned

# 2
HICX_PLOT_PYTHON=$pinned expect_refusal 1 "--minMappingQuality needs the columns mapq1 and mapq2" --pairsFile "$reserved" -rs "$data/DpnII.bed" --minMappingQuality 20

# 3
all=$pairs/small_test_all.pairs.gz
shuffled=$work/shuffled.pairs
gzip -dc "$all" | grep '^#' | grep -v -e '^#sorted:' -e '^#shape:' > "$shuffled"
gzip -dc "$all" | grep -v '^#' | awk 'BEGIN { srand(7) } { print rand() "\t" $0 }' | sort -k1,1 | cut -f2- >> "$shuffled"

run_matrix() {
    local name=$1
    shift
    mkdir -p "$work/$name"
    if ! "$tool" "$@" -bs 5000 -o "$work/$name/matrix.h5" --QCfolder "$work/$name/qc" > "$work/$name/stdout" 2> "$work/$name/stderr"; then
        fail "$name exited non-zero"
        sed 's/^/    /' "$work/$name/stderr"
        return 1
    fi
    grep -v '^File' "$work/$name/qc/QC.log" > "$work/$name/counts"
}

if run_matrix sorted_t1 --pairsFile "$all" --threads 1 &&
   run_matrix sorted_t4 --pairsFile "$all" --threads 4 --inputBufferSize 997 &&
   run_matrix shuffled_t4 --pairsFile "$shuffled" --threads 4 --inputBufferSize 997; then
    grep -q "in file order" "$work/sorted_t1/stderr" || fail "the sorted file was not checked in file order"
    grep -q "with a hash set" "$work/shuffled_t4/stderr" || fail "the shuffled file was not checked with a hash set"
    for other in sorted_t4 shuffled_t4; do
        cmp -s "$work/sorted_t1/matrix.h5" "$work/$other/matrix.h5" || fail "matrix of $other differs from sorted_t1"
        cmp -s "$work/sorted_t1/counts" "$work/$other/counts" || {
            fail "QC counts of $other differ from sorted_t1"
            diff "$work/sorted_t1/counts" "$work/$other/counts" | sed 's/^/    /'
        }
    done
    grep -q "^duplicated pairs	[1-9]" "$work/sorted_t1/counts" || fail "small_test_all.pairs.gz should hold duplicated pairs"
fi

# 4
unsorted=$work/declared_sorted.pairs
gzip -dc "$valid" | grep '^#' > "$unsorted"
gzip -dc "$valid" | grep -v '^#' | tac >> "$unsorted"
mkdir -p "$work/unsorted"
"$tool" --pairsFile "$unsorted" -bs 5000 -o "$work/unsorted/matrix.h5" --QCfolder "$work/unsorted/qc" > /dev/null 2> "$work/unsorted/stderr"
code=$?
[ "$code" = 1 ] || fail "a file declared sorted but not sorted exited $code instead of 1"
grep -q "although the header declares" "$work/unsorted/stderr" || fail "no order message: $(tail -1 "$work/unsorted/stderr")"
[ -e "$work/unsorted/matrix.h5" ] && fail "a matrix was written for the unsorted file"

# 5
plain=$pairs/R1_1000_all.pairs
mkdir -p "$work/counts"
if "$tool" --pairsFile "$plain" -bs 100000 --minMappingQuality 30 -o "$work/counts/matrix.h5" --QCfolder "$work/counts/qc" > /dev/null 2> "$work/counts/stderr"; then
    read -r unmapped not_unique low kept <<< "$(awk -F'\t' '!/^#/ {
        if ($8 ~ /[^URM]/) u++; else if ($8 ~ /M/) n++; else if ($9 < 30 || $10 < 30) l++; else k++
    } END { print u + 0, n + 0, l + 0, k + 0 }' "$plain")"
    log=$work/counts/qc/QC.log
    for check in "Sequenced reads	983" "One mate unmapped	$unmapped	" "One mate not unique	$not_unique	" \
                 "Low mapping quality	$low	" "Hi-C contacts	$kept	"; do
        grep -qF "$check" "$log" || fail "QC.log of R1_1000_all.pairs lacks '$check'"
    done
else
    fail "R1_1000_all.pairs exited non-zero"
    sed 's/^/    /' "$work/counts/stderr"
fi

if [ "$failures" != 0 ]; then
    echo "$failures failure(s)"
    exit 1
fi
echo "all checks passed"
