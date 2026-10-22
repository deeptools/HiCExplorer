#!/usr/bin/env bash
# Derives the inputs of the external capture Hi-C cases (chic*, hicConvertFormat
# --inputFormat chinput) from public data under the data root:
#   derive_capture.sh DATA_ROOT [CPP_BIN]
# Writes DATA_ROOT/derived/capture/. Needs the HiCExplorer 3.7.6 environment
# (PY_BIN) and the C++ hicConvertFormat (CPP_BIN, for the chinput -> cool step,
# which has no 3.7.6 counterpart). Every Python-side derived file is produced by
# the 3.7.6 reference, so the chain of a case's inputs never depends on the
# C++ tool under test except for the fragment-level cools.
#
# Data:
#   capture/EP_genomewide.chinput   genome-wide promoter capture Hi-C sample on
#                                   the CHiCAGO hg19 HindIII design (22,076 baits)
#   capture/k562_rep{1,2}.chinput   ENCODE K562 promoter capture Hi-C, replicate 1
#                                   (ENCFF427HFJ + ENCFF695CPL) and 2 (ENCFF255JIM +
#                                   ENCFF192DGQ), HiCUP + bam2chicago on the same
#                                   design. Linked here from K562_CHINPUT_DIR.
#   capture/K562_capture_rep{1,2}_10kb.cool
#                                   the same two K562 chinputs binned at 10 kb
#                                   (capturehic-to-hic scripts/chinput_to_cool.py,
#                                   fragment midpoints). Linked from K562_COOL_DIR.
set -euo pipefail
D=${1:?usage: derive_capture.sh DATA_ROOT [CPP_BIN]}
CPP_BIN=${2:-$HOME/src/HiCExplorer-v4/build/tools}
PY_BIN=${PY_BIN:-$HOME/miniconda3/envs/__hicexplorer@3.7.6/bin}
K562_CHINPUT_DIR=${K562_CHINPUT_DIR:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/hicup/k562_output}
K562_COOL_DIR=${K562_COOL_DIR:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/encode_cool}

C=$D/capture
O=$D/derived/capture
mkdir -p "$O"
for r in 1 2; do
    ln -sfn "$K562_CHINPUT_DIR/k562_rep$r/k562_rep$r.chinput" "$C/k562_rep$r.chinput"
    ln -sfn "$K562_COOL_DIR/K562_capture_rep${r}_10kb.cool" "$C/K562_capture_rep${r}_10kb.cool"
done
cd "$O"

# 1. Fragment-level cools (one bin per .rmap fragment) from the chinputs.
#    Bait-to-bait pairs are listed once per bait in a chinput with the same N
#    (100 % of pairs in K562; EP was filtered by R before export: 97.9 % of its
#    pairs in both directions with the same N); hicConvertFormat stores N once
#    (the larger N where the two rows differ). Same matrix for both sides of every case.
"$CPP_BIN/hicConvertFormat" --matrices "$C/EP_genomewide.chinput" --rmap "$C/Digest_Human_HindIII.rmap" \
    --inputFormat chinput --outputFormat cool -o EP_genomewide_cpp.cool
for r in 1 2; do
    "$CPP_BIN/hicConvertFormat" --matrices "$C/k562_rep$r.chinput" --rmap "$C/Digest_Human_HindIII.rmap" \
        --inputFormat chinput --outputFormat cool -o K562_rep${r}_frag.cool
done

# 2. Chromosome subsets with the 3.7.6 hicAdjustMatrix.
"$PY_BIN/hicAdjustMatrix" -m EP_genomewide_cpp.cool --chromosomes 20 21 22 --action keep -o EP_frag_chr20-22.cool
"$PY_BIN/hicAdjustMatrix" -m EP_genomewide_cpp.cool --chromosomes 1 --action keep -o EP_frag_chr1.cool
for r in 1 2; do
    "$PY_BIN/hicAdjustMatrix" -m K562_rep${r}_frag.cool --chromosomes 20 21 22 --action keep -o K562_rep${r}_frag_chr20-22.cool
    "$PY_BIN/hicAdjustMatrix" -m "$C/K562_capture_rep${r}_10kb.cool" --chromosomes 20 --action keep -o K562_rep${r}_10kb_chr20.cool
done

# 3. Reference points: every bait of the baitmap, at the bait fragment's midpoint,
#    named <first gene>_<baitID> so names are unique.
bait_bed() { awk -v c="$1" 'BEGIN{OFS="\t"} c=="all" || $1==c {split($5,g,";"); m=int(($2+$3)/2); print $1, m, m, g[1]"_"$4}' \
                 "$C/Digest_Human_HindIII_baits_e75_ID.baitmap"; }
bait_bed all > baits_genomewide.bed
for c in 1 20 21 22; do bait_bed $c > baits_chr$c.bed; done
cat baits_chr20.bed baits_chr21.bed baits_chr22.bed > baits_chr20-22.bed
awk '{if ($3 > m[$1]) m[$1] = $3} END {for (k in m) print k"\t"m[k]}' "$C/Digest_Human_HindIII.rmap" | sort -k1,1V > hg19_nochr.chrom.sizes
# The 10 kb K562 reference points that pass the 3.7.6 chicQualityControl at sparsity 0.05 (573 of 619).
tmp=$(mktemp -d "$O/qc.XXXX")
"$PY_BIN/chicQualityControl" -m K562_rep1_10kb_chr20.cool K562_rep2_10kb_chr20.cool -rp baits_chr20.bed --sparsity 0.05 \
    -o "$tmp/qc.bed" -oh "$tmp/h.png" -os "$tmp/s.png" -t 4
cp "$tmp/qc.bed" baits_chr20_qc_10kb.bed
rm -rf "$tmp"
# Twelve chr20 baits with significant interactions in EP (step 6 below), for the plotting cases.
for g in DHX35_412750 YWHAB_414406 UBE2C_414672 STAU1_415456 SRC_412345 SNHG17_412613 SCAND1_412068 \
         RPRD1B_412510 PIGU_411799 NCOA6_411808 MYT1_419121 CTCFL_417659; do
    awk -v g=$g '$4==g' baits_chr20.bed
done > baits_plot12.bed

# 4. Background models (3.7.6).
"$PY_BIN/chicViewpointBackgroundModel" -m EP_frag_chr20-22.cool K562_rep1_frag_chr20-22.cool K562_rep2_frag_chr20-22.cool \
    -rp baits_chr20-22.bed -o background_frag_chr20-22.txt -t 4
"$PY_BIN/chicViewpointBackgroundModel" -m K562_rep1_10kb_chr20.cool K562_rep2_10kb_chr20.cool \
    -rp baits_chr20_qc_10kb.bed -o background_10kb_chr20.txt -t 4

# 5-8. Viewpoints, significant interactions, aggregation and differential test (3.7.6).
chain() {  # name matrices... -- refpoints background
    local n=$1; shift; local m=(); while [ "$1" != "--" ]; do m+=("$1"); shift; done; shift
    local rp=$1 bg=$2
    "$PY_BIN/chicViewpoint" -m "${m[@]}" -rp "$rp" -bmf "$bg" --range 500000 500000 -o viewpoints_$n.hdf5 -t 4
    "$PY_BIN/chicSignificantInteractions" --interactionFile viewpoints_$n.hdf5 -bmf "$bg" --range 500000 500000 \
        --outFileNameSignificant significant_$n.hdf5 --outFileNameTarget target_$n.hdf5 \
        --xFoldBackground 2 --pValue 0.01 --combinationMode dual -t 4
    "$PY_BIN/chicAggregateStatistic" --interactionFile viewpoints_$n.hdf5 --targetFile target_$n.hdf5 \
        --outFileName aggregate_$n.hdf5 -t 4
    "$PY_BIN/chicDifferentialTest" --aggregatedFile aggregate_$n.hdf5 --alpha 0.05 --statisticTest fisher \
        --outFileName differential_$n.hdf5 -t 4
}
chain EP_K562r1_frag_chr20-22 EP_frag_chr20-22.cool K562_rep1_frag_chr20-22.cool -- baits_chr20-22.bed background_frag_chr20-22.txt
chain K562r1_r2_10kb_chr20 K562_rep1_10kb_chr20.cool K562_rep2_10kb_chr20.cool -- baits_chr20_qc_10kb.bed background_10kb_chr20.txt
chain EP_K562r1_frag_plot12 EP_frag_chr20-22.cool K562_rep1_frag_chr20-22.cool -- baits_plot12.bed background_frag_chr20-22.txt

# 9. A BED target list for chicAggregateStatistic: the chr20 refGene promoters
#    (TSS +-2 kb) of tracks/chicago_tracks/promoters.bed, three columns, sorted.
cut -f1-3 "$D/tracks/chicago_tracks/promoters.bed" | sort -k1,1 -k2,2n -u > targets_promoters_chr20.bed
rm -f errorLog.txt   # side output of the 3.7.6 chicAggregateStatistic

# 10. Inputs of the C++-only CHiCAGO cases (no 3.7.6 counterpart, so these come
#     from the C++ tools themselves): the background model of the genome-wide EP
#     sample from the chinput and from the fragment cool, and the scores from each.
CHICAGO_DESIGN_ARGS=(--rmap "$C/Digest_Human_HindIII.rmap" --baitmap "$C/Digest_Human_HindIII_baits_e75_ID.baitmap"
                     --nperbin "$C/Human_HindIII_NperBin.npb" --nbaitsperbin "$C/Human_HindIII_NBaitsPerBin.nbpb"
                     --proxOE "$C/Human_HindIII_proxOE.poe" --threads 4)
"$CPP_BIN/chicChicagoBackgroundModel" "${CHICAGO_DESIGN_ARGS[@]}" --chinput "$C/EP_genomewide.chinput" -o chicago_bg_EP_chinput.txt
"$CPP_BIN/chicChicagoBackgroundModel" "${CHICAGO_DESIGN_ARGS[@]}" --matrices EP_genomewide_cpp.cool -o chicago_bg_EP_matrix.txt
SCORE_ARGS=(--baitmap "$C/Digest_Human_HindIII_baits_e75_ID.baitmap" --rmap "$C/Digest_Human_HindIII.rmap" --threads 4)
"$CPP_BIN/chicChicagoScores" "${SCORE_ARGS[@]}" --backgroundModel chicago_bg_EP_chinput.txt --chinput "$C/EP_genomewide.chinput" -o chicago_scores_EP_chinput.txt
"$CPP_BIN/chicChicagoScores" "${SCORE_ARGS[@]}" --backgroundModel chicago_bg_EP_matrix.txt --matrices EP_genomewide_cpp.cool -o chicago_scores_EP_matrix.txt

# 11. ADAMTS1_423843 (chr21) makes the 3.7.6 chicViewpoint fail at --range 200000
#     200000 on the fragment cools ("operands could not be broadcast together with
#     shapes (158,) (157,)"); the C++ port fails the same way. One file with only
#     that bait, and the chr20-22 set without it for the background-model case,
#     whose downstream validator runs chicViewpoint at 200 kb.
awk '$4=="ADAMTS1_423843"' baits_chr20-22.bed > baits_ADAMTS1.bed
awk '$4!="ADAMTS1_423843"' baits_chr20-22.bed > baits_chr20-22_noADAMTS1.bed
