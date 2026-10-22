#!/usr/bin/env bash
# Derived inputs of the matrix-tool external cases (ids ext.hic*), written to
# DATA_ROOT/derived/matrix. Every input is a deterministic subset of a public
# dataset in the data root built by make_data_root.sh.
#   derive_matrix.sh DATA_ROOT
# Needs the HiCExplorer 3.7.6 conda environment (cooler, samtools, the 3.7.6
# tools) at $HICX376 (default ~/miniconda3/envs/__hicexplorer@3.7.6).
set -euo pipefail
D=${1:?usage: derive_matrix.sh DATA_ROOT}
HICX376=${HICX376:-$HOME/miniconda3/envs/__hicexplorer@3.7.6}
PY=$HICX376/bin/python3
O=$D/derived/matrix
mkdir -p "$O"

# 1. Raw-count subsets of the ENCODE ENCSR946AKI GM12878 replicates r8 and r12 and
#    of ENCODE A549 (hg38): chr19 to chr22 at 25 kb, chr21 and chr22 at 10 kb.
#    Bins keep chrom/start/end only (the hic2cool SCALE/VC/weight columns are
#    dropped), pixels are copied unchanged, int32 counts.
subset() {   # in.cool out.cool chrom...
    "$PY" -I - "$@" <<'PYEOF'
import sys, cooler, pandas as pd
src, dst, chroms = sys.argv[1], sys.argv[2], sys.argv[3:]
c = cooler.Cooler(src)
bins = c.bins()[:]
keep = bins.index[bins.chrom.isin(chroms)].values
# order the bins by the requested chromosome order
parts = [bins[bins.chrom == ch] for ch in chroms]
new_bins = pd.concat(parts)[["chrom", "start", "end"]]
old_ids = new_bins.index.values
remap = pd.Series(range(len(old_ids)), index=old_ids)
pix = []
for i, ch in enumerate(chroms):
    for ch2 in chroms[i:]:
        p = c.matrix(balance=False, as_pixels=True).fetch(ch, ch2)
        if c.chromnames.index(ch) > c.chromnames.index(ch2):
            p = p.rename(columns={"bin1_id": "bin2_id", "bin2_id": "bin1_id"})
        pix.append(p[["bin1_id", "bin2_id", "count"]])
pix = pd.concat(pix)
pix["bin1_id"] = remap.loc[pix.bin1_id.values].values
pix["bin2_id"] = remap.loc[pix.bin2_id.values].values
swap = pix.bin1_id > pix.bin2_id
pix.loc[swap, ["bin1_id", "bin2_id"]] = pix.loc[swap, ["bin2_id", "bin1_id"]].values
pix = pix.sort_values(["bin1_id", "bin2_id"]).reset_index(drop=True)
new_bins = new_bins.reset_index(drop=True)
new_bins["chrom"] = new_bins["chrom"].astype(str)
chromsizes = c.chromsizes.loc[chroms]
cooler.create_cooler(dst, new_bins, pix, dtypes={"count": "int32"}, ordered=True,
                     assembly=c.info.get("genome-assembly"))
PYEOF
}
C=$D/cool
for s in gm12878_ENCSR946AKI_r8 gm12878_ENCSR946AKI_r12 a549; do
    subset "$C/${s}_25kb.cool" "$O/${s}_25kb_chr19-22.cool" chr19 chr20 chr21 chr22
    subset "$C/${s}_10kb.cool" "$O/${s}_10kb_chr21-22.cool" chr21 chr22
done
subset "$C/gm12878_ENCSR946AKI_r8_25kb.cool" "$O/gm12878_ENCSR946AKI_r8_25kb_chr22.cool" chr22

# 2. The first N read pairs of SRA SRR927086 (GM12878 HindIII Hi-C, hg19, read 1
#    and read 2 aligned separately with minibwa). Records are kept while the
#    read number in the name (SRR927086.<n>) is <= N, so both files hold the
#    same read names in the same order, supplementary records included; then
#    converted to BAM unchanged (no sorting).
N_PAIRS=${N_PAIRS:-2000000}
for s in 1 2; do
    awk -v n="$N_PAIRS" 'BEGIN { FS = "\t" }
        /^@/ { print; next }
        { split($1, a, "."); if (a[2] + 0 > n) exit; print }' "$D/reads/SRR927086_$s.sam" \
      | "$HICX376/bin/samtools" view -b -o "$O/SRR927086_${s}_first${N_PAIRS}.bam" -
done

# 3. hg19 chr19, chr21 and chr22 (UCSC hg19.fa, soft-masked) for hicFindRestSite.
"$HICX376/bin/samtools" faidx "$D/genome/hg19.fa" chr19 chr21 chr22 > "$O/hg19_chr19_chr21_chr22.fa"

# 4. Genome-wide hg19 HindIII (AAGCTT) cut sites, the restriction site file of the
#    hicBuildMatrix and hicQuickQC cases, written by the HiCExplorer 3.7.6 tool.
"$HICX376/bin/hicFindRestSite" -f "$D/genome/hg19.fa" -p AAGCTT -o "$O/hg19_HindIII.bed"

# 5. h5 copies of the 25 kb subsets (the 3.7.6 hicConvertFormat, cool -> h5) and
#    the hg38 sizes of chr19 to chr22 taken from the subset's own chroms table.
for s in gm12878_ENCSR946AKI_r8 gm12878_ENCSR946AKI_r12 a549; do
    "$HICX376/bin/hicConvertFormat" -m "$O/${s}_25kb_chr19-22.cool" --inputFormat cool \
        --outputFormat h5 -o "$O/${s}_25kb_chr19-22.h5"
done
"$PY" -I -c 'import sys, cooler; c = cooler.Cooler(sys.argv[1]); print(c.chromsizes.to_csv(sep="\t", header=False), end="")' \
    "$O/gm12878_ENCSR946AKI_r8_25kb_chr19-22.cool" > "$O/hg38_chr19-22.chrom.sizes"

# 6. hic2cool (through the 3.7.6 hicConvertFormat) of GSE63525 GM12878
#    combined_30 at 1 Mb and 2.5 Mb: the independent pixels the C++
#    hic -> mcool case is checked against (gm12878_combined30_<res>.cool).
"$HICX376/bin/hicConvertFormat" -m "$D/hic/GSE63525_GM12878_combined_30.hic" --inputFormat hic \
    --outputFormat cool -o "$O/gm12878_combined30.cool" -r 1000000 2500000

# 7. Regions for hicAdjustMatrix --regions (hg38): the chr19 and chr20
#    centromeres with their flanks, and the acrocentric p arm of chr21.
printf 'chr19\t24000000\t28000000\nchr20\t26000000\t30000000\nchr21\t0\t12000000\n' > "$O/hg38_chr19-21_regions.bed"

# 8. QC logs for the hicQC / hicPrepareQCreport cases: the 3.7.6 hicBuildMatrix
#    (HindIII, 100 kb) and hicBuildMatrixMicroC (100 kb) on the reads of step 2.
B1=$O/SRR927086_1_first${N_PAIRS}.bam; B2=$O/SRR927086_2_first${N_PAIRS}.bam
T=$(mktemp -d)
"$HICX376/bin/hicBuildMatrix" -s "$B1" "$B2" -rs "$O/hg19_HindIII.bed" --restrictionSequence AAGCTT \
    --danglingSequence AGCT -bs 100000 --outFileName "$T/hindiii.h5" \
    --QCfolder "$O/qc_SRR927086_${N_PAIRS:0:1}M_HindIII" --threads 4
"$HICX376/bin/hicBuildMatrixMicroC" -s "$B1" "$B2" -bs 100000 --outFileName "$T/microc.h5" \
    --QCfolder "$O/qc_SRR927086_${N_PAIRS:0:1}M_microc" --threads 4
rm -rf "$T"

# 9. Independent coarsening oracle for the C++ hicMergeMatrixBins cases, whose
#    bin layout deliberately deviates from 3.7.6 (the merge groups follow
#    --chromosomeSizes, so no bin is dropped): cooler coarsen of the subsets.
"$HICX376/bin/cooler" coarsen -k 4 -o "$O/gm12878_ENCSR946AKI_r8_25kb_chr19-22.coarsen4.cool" \
    "$O/gm12878_ENCSR946AKI_r8_25kb_chr19-22.cool"
"$HICX376/bin/cooler" coarsen -k 10 -o "$O/a549_25kb_chr19-22.coarsen10.cool" "$O/a549_25kb_chr19-22.cool"

# 10. A second Hi-C sample for the multi-log hicQC cases: read pairs N+1 to 2N of
#     SRR927086 (selected by read number as in step 2) and the QC log of the
#     3.7.6 hicBuildMatrix on them (HindIII, 100 kb).
for s in 1 2; do
    awk -v n="$N_PAIRS" 'BEGIN { FS = "\t" }
        /^@/ { print; next }
        { split($1, a, "."); k = a[2] + 0; if (k > 2 * n) exit; if (k > n) print }' "$D/reads/SRR927086_$s.sam" \
      | "$HICX376/bin/samtools" view -b -o "$O/SRR927086_${s}_second${N_PAIRS}.bam" -
done
T=$(mktemp -d)
"$HICX376/bin/hicBuildMatrix" -s "$O/SRR927086_1_second${N_PAIRS}.bam" "$O/SRR927086_2_second${N_PAIRS}.bam" \
    -rs "$O/hg19_HindIII.bed" --restrictionSequence AAGCTT --danglingSequence AGCT -bs 100000 \
    --outFileName "$T/hindiii_second.h5" --QCfolder "$O/qc_SRR927086_second2M_HindIII" --threads 4
rm -rf "$T"
