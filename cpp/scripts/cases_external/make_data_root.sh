#!/usr/bin/env bash
# Builds the data root of the external-data cases: symlinks to public datasets
# plus inputs derived from them. Every source path can be overridden.
#   make_data_root.sh DATA_ROOT
# then: equiv.py run --cases-dir cpp/scripts/cases_external --data DATA_ROOT ...
set -euo pipefail
D=${1:?usage: make_data_root.sh DATA_ROOT}

ENCODE_COOL_DIR=${ENCODE_COOL_DIR:-/mnt/processed/wolffjoa/hic-multiome/encode_hic}   # ENCODE GM12878 ENCSR946AKI r8/r12, A549
GSE63525_DIR=${GSE63525_DIR:-/mnt/processed/wolffjoa/restart_of_prediction_project}   # GEO GSE63525 combined_30 .hic
MCOOL_MICROC=${MCOOL_MICROC:-/mnt/processed/wolffjoa/hic-multiome/matrices/hff_microc.mcool}   # 4DN HFF Micro-C
MCOOL_IMR90=${MCOOL_IMR90:-/mnt/processed/wolffjoa/hyperparameter_publication/imr90.mcool}
READS_DIR=${READS_DIR:-/mnt/processed/wolffjoa/hicexplorer_bench/fastq}               # SRA SRR927086, hg19, minibwa --hic
HG19_FASTA=${HG19_FASTA:-/mnt/processed/wolffjoa/vardictcpp_validation_zenodo/reference/hg19.fa}
CHICAGO_DESIGN=${CHICAGO_DESIGN:-/mnt/processed/wolffjoa/capturehic-to-hic-raw/chicago_design/Human_hg19}
CHINPUT=${CHINPUT:-/mnt/processed/wolffjoa/hicexplorer_bench/chicago_genomewide/EP_genomewide.chinput}
CHICAGO_TRACKS=${CHICAGO_TRACKS:-/mnt/processed/wolffjoa/hicexplorer_bench/chicago_genomewide/tracks}
REPLICATES=${REPLICATES:-$HOME/data/Hi-cGAN/biological_replicates/GSE234292}            # GEO GSE234292
DOWNLOADS=${DOWNLOADS:-/mnt/processed/wolffjoa/hicexplorer4_extval/data/downloads}     # GSE63525 loop and domain lists, ENCODE CTCF peaks

mkdir -p "$D"/{hic,cool,mcool,reads,genome,capture,replicates,tracks,derived}
for f in gm12878_ENCSR946AKI_r8_25kb gm12878_ENCSR946AKI_r12_25kb gm12878_ENCSR946AKI_r8_10kb \
         gm12878_ENCSR946AKI_r12_10kb a549_25kb a549_10kb; do
    ln -sfn "$ENCODE_COOL_DIR/$f.cool" "$D/cool/$f.cool"
done
ln -sfn "$GSE63525_DIR/GSE63525_GM12878_insitu_primary+replicate_combined_30.hic" "$D/hic/GSE63525_GM12878_combined_30.hic"
ln -sfn "$GSE63525_DIR/GSE63525_K562_combined_30.hic" "$D/hic/GSE63525_K562_combined_30.hic"
ln -sfn "$MCOOL_MICROC" "$D/mcool/hff_microc.mcool"
ln -sfn "$MCOOL_IMR90" "$D/mcool/imr90.mcool"
for s in SRR927086_1 SRR927086_2; do ln -sfn "$READS_DIR/$s.sam" "$D/reads/$s.sam"; done
ln -sfn "$HG19_FASTA" "$D/genome/hg19.fa"
for f in "$CHICAGO_DESIGN"/*; do ln -sfn "$f" "$D/capture/$(basename "$f")"; done
ln -sfn "$CHINPUT" "$D/capture/EP_genomewide.chinput"
ln -sfn "$REPLICATES" "$D/replicates/GSE234292"
ln -sfn "$CHICAGO_TRACKS" "$D/tracks/chicago_tracks"
ln -sfn "$DOWNLOADS" "$D/downloads"

# Stranded CTCF motifs at the anchors of the GSE63525 GM12878 HiCCUPS loops (hg19).
zcat "$D/downloads/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt.gz" \
  | awk -F'\t' 'NR>1 {
        if ($21 != "na" && $21 != "" && ($24 == "p" || $24 == "n")) print $1"\t"$21"\t"$22"\tCTCF\t0\t"($24 == "p" ? "+" : "-");
        if ($26 != "na" && $26 != "" && ($29 == "p" || $29 == "n")) print $4"\t"$26"\t"$27"\tCTCF\t0\t"($29 == "p" ? "+" : "-") }' \
  | sort -k1,1 -k2,2n -u > "$D/derived/gm12878_ctcf_motifs_hg19_nochr.bed"
awk 'BEGIN { OFS = "\t" } { $1 = "chr"$1; print }' "$D/derived/gm12878_ctcf_motifs_hg19_nochr.bed" \
  > "$D/derived/gm12878_ctcf_motifs_hg19.bed"
