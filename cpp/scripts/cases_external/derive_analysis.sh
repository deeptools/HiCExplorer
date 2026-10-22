#!/usr/bin/env bash
# Recreates the inputs of the analysis-tool external cases (hicFindTADs, hicArrowhead,
# hicDetectLoops, hicPCA, ...) under DATA_ROOT/derived/analysis from the data root
# that make_data_root.sh builds.
#   derive_analysis.sh DATA_ROOT
# Every input is real public data: chromosome subsets of the GSE63525, ENCODE and
# GSE234292 matrices, plus a few files computed from them with HiCExplorer 3.7.6.
# Override PY (3.7.6 interpreter), HICX_BIN (built C++ tools), MM10_FASTA.
set -euo pipefail
D=${1:?usage: derive_analysis.sh DATA_ROOT}
PY=${PY:-$HOME/miniconda3/envs/__hicexplorer@3.7.6/bin/python3}
PYBIN=$(dirname "$PY")
HICX_BIN=${HICX_BIN:-$HOME/src/HiCExplorer-v4/build/tools}
MM10_FASTA=${MM10_FASTA:-$HOME/data/facount/mm10.fa}   # UCSC mm10, only for the GC track
O="$D/derived/analysis"
mkdir -p "$O"
TMP=$(mktemp -d "$O/tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# subset IN OUT CHROM...: keep the given chromosomes (adjacent in the file's bin
# order) with every bin column (raw counts plus weight/KR/VC...) and the intra- and
# inter-chromosomal pixels among them. With RENAME_WEIGHT=KR the KR column (hic2cool
# convention, multiplicative) becomes the 'weight' column that cooler and HiCExplorer read.
subset() {
    "$PY" -I - "$@" <<'EOF'
import os, sys, h5py, numpy as np, pandas as pd, cooler
src, dst, chroms = sys.argv[1], sys.argv[2], sys.argv[3:]
c = cooler.Cooler(src)
bins = c.bins()[:]
keep = bins['chrom'].isin(chroms).to_numpy()
idx = np.flatnonzero(keep)
lo, hi = idx[0], idx[-1] + 1
assert len(idx) == hi - lo, 'chromosomes must be adjacent in the bin table'
with h5py.File(src, 'r') as h5:
    off = h5['indexes/bin1_offset']
    a, b = int(off[lo]), int(off[hi])
    b1 = h5['pixels/bin1_id'][a:b]; b2 = h5['pixels/bin2_id'][a:b]; ct = h5['pixels/count'][a:b]
m = b2 < hi
pixels = pd.DataFrame({'bin1_id': b1[m] - lo, 'bin2_id': b2[m] - lo, 'count': ct[m]})
out_bins = bins.iloc[lo:hi].reset_index(drop=True)
out_bins['chrom'] = out_bins['chrom'].astype(str)
rename = os.environ.get('RENAME_WEIGHT')
if rename:
    out_bins = out_bins[['chrom', 'start', 'end', rename]].rename(columns={rename: 'weight'})
chromsizes = c.chromsizes[chroms]
cooler.create_cooler(dst, out_bins, pixels, assembly=c.info.get('genome-assembly'),
                     ordered=True, symmetric_upper=True, dtypes={'count': ct.dtype},
                     metadata={'source': os.path.basename(src), 'chromosomes': chroms})
print(dst, len(out_bins), 'bins', len(pixels), 'pixels')
EOF
}

# --- GSE63525 GM12878 in situ combined_30 (hg19, names without chr), KR as weight
for res in 25000 10000 100000; do
    [ -s "$O/gse63525_gm12878_chr17_19_${res}.cool" ] && continue   # rerun: keep what exists
    "$HICX_BIN/hicConvertFormat" -m "$D/hic/GSE63525_GM12878_combined_30.hic" --inputFormat hic \
        --outputFormat cool -r $res -o "$TMP/gm_gse63525.cool"
    RENAME_WEIGHT=KR subset "$TMP/gm_gse63525_$res.cool" "$O/gse63525_gm12878_chr17_19_${res}.cool" 17 18 19
    rm -f "$TMP/gm_gse63525_$res.cool"
done

# --- ENCODE ENCSR946AKI GM12878 r8, r12 and A549 (hg38, chr prefix), 25 kb, with weight
for s in gm12878_ENCSR946AKI_r8 gm12878_ENCSR946AKI_r12 a549; do
    [ -s "$O/${s}_chr17_19_25kb.cool" ] || subset "$D/cool/${s}_25kb.cool" "$O/${s}_chr17_19_25kb.cool" chr17 chr18 chr19
done

# --- GSE234292 mESC wild type and BPTF knockout, two replicates each (mm10), raw int counts + weight
for s in GSM7461656_HiC_wt_rep1 GSM7461657_HiC_wt_rep2 GSM7461658_HiC_BPTFko_rep1 GSM7461659_HiC_BPTFko_rep2; do
    for res in 25000 100000; do
        [ -s "$O/${s}_chr17_19_${res}.cool" ] || subset "$D/replicates/GSE234292/${s}_${res}.cool" "$O/${s}_chr17_19_${res}.cool" chr17 chr18 chr19
    done
done

# --- copies without the weight column (HiCExplorer then reads the raw counts)
noweight() {
    cp "$1" "$2"
    "$PY" -I -c "import h5py, sys; f = h5py.File(sys.argv[1], 'r+'); del f['bins/weight']; f.close()" "$2"
}
noweight "$O/gse63525_gm12878_chr17_19_25000.cool" "$O/gse63525_gm12878_chr17_19_25000_raw.cool"
noweight "$O/gse63525_gm12878_chr17_19_10000.cool" "$O/gse63525_gm12878_chr17_19_10000_raw.cool"
for s in GSM7461656_HiC_wt_rep1 GSM7461658_HiC_BPTFko_rep1; do
    noweight "$O/${s}_chr17_19_25000.cool" "$O/${s}_chr17_19_25000_raw.cool"
done

# --- annotation tracks
# ENCODE GM12878 CTCF IDR peaks hg19 without the chr prefix (GSE63525 naming), sorted
zcat "$D/downloads/gm12878_ctcf_hg19.narrowPeak.gz" | sed 's/^chr//' | sort -k1,1 -k2,2n > "$O/gm12878_ctcf_hg19_nochr.narrowPeak"
zcat "$D/downloads/gm12878_ctcf_hg19.narrowPeak.gz" | sort -k1,1 -k2,2n > "$O/gm12878_ctcf_hg19.narrowPeak"
zcat "$D/downloads/gm12878_ctcf_GRCh38.narrowPeak.gz" | sort -k1,1 -k2,2n > "$O/gm12878_ctcf_GRCh38.narrowPeak"
# hg19 refGene transcripts on chr17-19 as BED6 without chr prefix (hicPCA --extraTrack)
zcat "$D/tracks/chicago_tracks/refGene.txt.gz" \
  | awk 'BEGIN { OFS = "\t" } $3 == "chr17" || $3 == "chr18" || $3 == "chr19" { sub(/^chr/, "", $3); print $3, $5, $6, $13, 0, $4 }' \
  | sort -k1,1 -k2,2n -u > "$O/hg19_refGene_chr17_19_nochr.bed"
# stranded CTCF motifs at GSE63525 loop anchors on chr17-19 (hicAverageRegions, hicAggregateContacts)
awk '$1 == "17" || $1 == "18" || $1 == "19"' "$D/derived/gm12878_ctcf_motifs_hg19_nochr.bed" > "$O/gm12878_ctcf_motifs_hg19_chr17_19_nochr.bed"
# mm10 GC content of chr17-19 in 100 kb bins (hicDifferentialAnalysis compartments --gcTrack)
"$PY" -I - "$MM10_FASTA" "$O/mm10_gc_chr17_19_100kb.bedgraph" <<'PYEOF'
import sys, pysam
fa = pysam.FastaFile(sys.argv[1])
with open(sys.argv[2], 'w') as out:
    for chrom in ('chr17', 'chr18', 'chr19'):
        n = fa.get_reference_length(chrom)
        for start in range(0, n, 100000):
            s = fa.fetch(chrom, start, min(n, start + 100000)).upper()
            acgt = sum(s.count(b) for b in 'ACGT')
            if acgt:
                out.write(f'{chrom}\t{start}\t{min(n, start + 100000)}\t{(s.count("G") + s.count("C")) / acgt:.6f}\n')
PYEOF

# --- results of HiCExplorer 3.7.6 used as inputs of downstream tools
# TAD domains
"$PYBIN/hicFindTADs" -m "$O/gse63525_gm12878_chr17_19_25000.cool" --correctForMultipleTesting fdr -p 4 --outPrefix "$TMP/g25"
"$PYBIN/hicFindTADs" -m "$O/gse63525_gm12878_chr17_19_10000.cool" --correctForMultipleTesting fdr -p 4 --outPrefix "$TMP/g10"
"$PYBIN/hicFindTADs" -m "$O/gse63525_gm12878_chr17_19_100000.cool" --correctForMultipleTesting fdr -p 4 --outPrefix "$TMP/g100"
cp "$TMP/g25_domains.bed" "$O/gse63525_gm12878_25kb_domains.bed"
cp "$TMP/g25_boundaries.bed" "$O/gse63525_gm12878_25kb_boundaries.bed"
cp "$TMP/g10_domains.bed" "$O/gse63525_gm12878_10kb_domains.bed"
cp "$TMP/g100_domains.bed" "$O/gse63525_gm12878_100kb_domains.bed"
"$PYBIN/hicFindTADs" -m "$O/gm12878_ENCSR946AKI_r8_chr17_19_25kb.cool" --minDepth 75000 --maxDepth 500000 --step 25000 \
    --minBoundaryDistance 75000 --correctForMultipleTesting fdr --thresholdComparisons 0.05 -p 4 --outPrefix "$TMP/r8"
cp "$TMP/r8_domains.bed" "$O/gm12878_ENCSR946AKI_r8_25kb_domains.bed"
cp "$TMP/r8_boundaries.bed" "$O/gm12878_ENCSR946AKI_r8_25kb_boundaries.bed"
"$PYBIN/hicFindTADs" -m "$O/GSM7461656_HiC_wt_rep1_chr17_19_25000.cool" --minDepth 100000 --maxDepth 400000 --step 50000 \
    --correctForMultipleTesting fdr --thresholdComparisons 0.05 -p 4 --outPrefix "$TMP/wt"
cp "$TMP/wt_domains.bed" "$O/GSM7461656_HiC_wt_rep1_25kb_domains.bed"
"$PYBIN/hicFindTADs" -m "$O/GSM7461656_HiC_wt_rep1_chr17_19_100000.cool" --minDepth 300000 --maxDepth 1000000 --step 100000 \
    --correctForMultipleTesting fdr --thresholdComparisons 0.05 -p 4 --outPrefix "$TMP/wt100"
cp "$TMP/wt100_domains.bed" "$O/GSM7461656_HiC_wt_rep1_100kb_domains.bed"
# loops
"$PYBIN/hicDetectLoops" -m "$O/gse63525_gm12878_chr17_19_10000.cool" --chromosomes 17 -t 1 -tpc 4 -o "$O/gse63525_gm12878_10kb_chr17_loops.bedgraph"
"$PYBIN/hicDetectLoops" -m "$O/gse63525_gm12878_chr17_19_25000.cool" --maxLoopDistance 3000000 -pit 5 -p 0.05 -t 3 -tpc 1 -o "$O/gse63525_gm12878_25kb_loops.bedgraph"
for s in GSM7461656_HiC_wt_rep1 GSM7461658_HiC_BPTFko_rep1; do
    "$PYBIN/hicDetectLoops" -m "$O/${s}_chr17_19_25000_raw.cool" -pit 5 -p 0.05 -t 3 -tpc 1 -o "$O/${s}_25kb_loops.bedgraph"
done
# compartments: obs/exp matrix and first eigenvector of GSE63525 100 kb
"$PYBIN/hicPCA" -m "$O/gse63525_gm12878_chr17_19_100000.cool" -o "$O/gse63525_gm12878_100kb_pca1.bedgraph" "$TMP/pca2.bedgraph" \
    -f bedgraph --whichEigenvectors 1 2 --method lieberman --extraTrack "$O/hg19_refGene_chr17_19_nochr.bed" \
    --obsexpMatrix "$O/gse63525_gm12878_100kb_obsexp.cool"
# averaged regions around the CTCF motifs (hicPlotAverageRegions input)
"$PYBIN/hicAverageRegions" -m "$O/gse63525_gm12878_chr17_19_25000.cool" --regions "$O/gm12878_ctcf_motifs_hg19_chr17_19_nochr.bed" \
    --range 500000 500000 -o "$O/gse63525_gm12878_25kb_ctcf_average.npz"
"$PYBIN/hicAverageRegions" -m "$O/gse63525_gm12878_chr17_19_10000.cool" --regions "$O/gm12878_ctcf_motifs_hg19_chr17_19_nochr.bed" \
    --rangeInBins 30 30 -cb center -o "$O/gse63525_gm12878_10kb_ctcf_average_bins.npz"
# a distance-dependent p-value threshold file for hicDetectLoops -pp (10 kb, 0 to 2 Mb)
"$PYBIN/hicCreateThresholdFile" --range 0 2000000 -tv 0.05 -r 10000 -o "$O/threshold_10kb_0_2mb.txt"
