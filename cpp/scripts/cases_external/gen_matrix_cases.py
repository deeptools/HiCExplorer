"""Writes cpp/scripts/cases_external/<tool>.json for the matrix-tool group (ext.hic* cases).

    python3 gen_matrix_cases.py cpp/scripts/cases_external

Inputs under {data}/derived/matrix come from derive_matrix.sh."""
import json
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
DM = "{data}/derived/matrix"
R8 = f"{DM}/gm12878_ENCSR946AKI_r8_25kb_chr19-22.cool"
R12 = f"{DM}/gm12878_ENCSR946AKI_r12_25kb_chr19-22.cool"
A5 = f"{DM}/a549_25kb_chr19-22.cool"
R8H = f"{DM}/gm12878_ENCSR946AKI_r8_25kb_chr19-22.h5"
R12H = f"{DM}/gm12878_ENCSR946AKI_r12_25kb_chr19-22.h5"
A5H = f"{DM}/a549_25kb_chr19-22.h5"
R8_10 = f"{DM}/gm12878_ENCSR946AKI_r8_10kb_chr21-22.cool"
R12_10 = f"{DM}/gm12878_ENCSR946AKI_r12_10kb_chr21-22.cool"
A5_10 = f"{DM}/a549_10kb_chr21-22.cool"
SIZES = f"{DM}/hg38_chr19-22.chrom.sizes"
HIC_GM = "{data}/hic/GSE63525_GM12878_combined_30.hic"
HIC_K562 = "{data}/hic/GSE63525_K562_combined_30.hic"
BAM1 = f"{DM}/SRR927086_1_first2000000.bam"
BAM2 = f"{DM}/SRR927086_2_first2000000.bam"
HINDIII = f"{DM}/hg19_HindIII.bed"

SUB25 = ("ENCODE ENCSR946AKI GM12878 replicate {r} (hg38, 25 kb), raw counts of chr19 to "
         "chr22, derived by derive_matrix.sh step 1")
N_R8 = SUB25.format(r=8)
N_SUBSETS = ("ENCODE ENCSR946AKI GM12878 replicates 8 and 12 and ENCODE A549 (hg38, 25 kb), raw "
             "counts of chr19 to chr22 (8,825 bins), derive_matrix.sh step 1")

cases = {}


def add(tool, cid, args, outputs, notes, **extra):
    case = {"id": f"ext.{tool}.{cid}", "tool": tool, "tier": extra.pop("tier", 1),
            "args": args, "outputs": outputs, "notes": notes}
    case.update(extra)
    cases.setdefault(tool, []).append(case)


def o(path, fmt, cls, **kw):
    entry = {"path": "{out}/" + path, "format": fmt, "class": cls}
    entry.update(kw)
    return entry


# ---------------------------------------------------------------- hicFindRestSite
add("hicFindRestSite", "hg19_chr19_21_22.HindIII",
    ["-f", f"{DM}/hg19_chr19_chr21_chr22.fa", "-p", "AAGCTT", "-o", "{out}/sites.bed"],
    [o("sites.bed", "plain", "E0")],
    "UCSC hg19 chr19, chr21 and chr22 (soft-masked, samtools faidx, derive_matrix.sh step 3), "
    "HindIII AAGCTT.", tier=2)
add("hicFindRestSite", "hg19_genome.HindIII_DpnII",
    ["-f", "{data}/genome/hg19.fa", "-p", "AAGCTT", "GATC", "-o", "{out}/sites.bed"],
    [o("sites.bed", "plain", "E0")],
    "The whole UCSC hg19 FASTA (93 sequences, 3.1 Gb) with HindIII and DpnII; about 8 M sites.",
    tier=2, large=True)

# ---------------------------------------------------------------- hicBuildMatrix
QC_FILES = ["QC.log", "QC_table.txt", "discarded_table.txt", "distance_table.txt",
            "read_orientation_table.txt", "unmapable_table.txt"]
# The repository's hicBuildMatrix budget formula (cases/hicBuildMatrix.json):
# 2 * nnz_out * 12 B + threads * 64 MB + 64 MB, nnz_out about 1.17 M here.
BM_BUDGET = 348
BM_NOTE = ("SRA SRR927086 (GM12878 in situ HindIII Hi-C, hg19), the first 2,000,000 read pairs of "
           "the minibwa-aligned read 1 and read 2 SAM files as BAM (derive_matrix.sh step 2), "
           "against the genome-wide hg19 HindIII sites of the 3.7.6 hicFindRestSite (step 4). ")
add("hicBuildMatrix", "SRR927086_2M.HindIII.bs100kb",
    ["-s", BAM1, BAM2, "-rs", HINDIII, "--restrictionSequence", "AAGCTT",
     "--danglingSequence", "AGCT", "-bs", "100000", "--outFileName", "{out}/matrix.h5",
     "-b", "{out}/valid.bam", "--QCfolder", "{out}/qc", "--threads", "4"],
    [o("matrix.h5", "h5", "E2")],
    BM_NOTE + "100 kb bins, h5 output; -b writes the valid-pair BAM, which is not declared: the "
    "plain comparator does not finish on a 232 MB binary within 200 s (the two BAMs were "
    "byte-identical by cmp, 2026-10-22). The QC files name the output path "
    "and are compared in the qc_log case instead.",
    tier=4, threads_arg="--threads", large=True,
    memory={"budget_mb": BM_BUDGET})
add("hicBuildMatrix", "SRR927086_2M.HindIII.restrictionFragments",
    ["-s", BAM1, BAM2, "-rs", HINDIII, "--restrictionSequence", "AAGCTT",
     "--danglingSequence", "AGCT", "--minDistance", "300", "--maxLibraryInsertSize", "1000",
     "--outFileName", "{out}/matrix.cool", "--QCfolder", "{out}/qc", "--threads", "4"],
    [o("matrix.cool", "cool", "E1")],
    BM_NOTE + "No --binSize: one bin per HindIII fragment (846,132 sites, fragments closer than "
    "300 bp merged), cool output.", tier=4, threads_arg="--threads", large=True,
    memory={"budget_mb": BM_BUDGET})
add("hicBuildMatrix", "SRR927086_2M.HindIII.qc_log",
    ["-s", BAM1, BAM2, "-rs", HINDIII, "--restrictionSequence", "AAGCTT",
     "--danglingSequence", "AGCT", "-bs", "1000000", "--outFileName", "ext_bm_qc_matrix.h5",
     "--QCfolder", "{out}/qc", "--threads", "4"],
    [o("qc/" + f, "plain", "E0") for f in QC_FILES],
    BM_NOTE + "1 Mb bins; --outFileName relative to the shared work directory, so the File "
    "line of QC.log is the same string on both sides and QC.log is compared byte for byte.",
    tier=4, threads_arg="--threads", large=True,
    memory={"budget_mb": BM_BUDGET})

# ---------------------------------------------------------------- hicBuildMatrixMicroC
add("hicBuildMatrixMicroC", "SRR927086_2M.bs100kb",
    ["-s", BAM1, BAM2, "-bs", "100000", "--outFileName", "{out}/matrix.h5",
     "-b", "{out}/valid.bam", "--QCfolder", "{out}/qc", "--threads", "4"],
    [o("matrix.h5", "h5", "E2")],
    "The same 2,000,000 SRR927086 read pairs as the hicBuildMatrix cases (derive_matrix.sh step 2), "
    "Micro-C mode with bin size only (no restriction sites). The data is HindIII Hi-C; the "
    "case checks the Micro-C code path on real alignments, not Micro-C biology. valid.bam is "
    "written but not declared (see the hicBuildMatrix bs100kb case).",
    tier=4, threads_arg="--threads", large=True,
    memory={"budget_mb": BM_BUDGET})
add("hicBuildMatrixMicroC", "SRR927086_2M.qc_log.bs1Mb",
    ["-s", BAM1, BAM2, "-bs", "1000000", "--outFileName", "ext_bmmc_qc_matrix.cool",
     "--QCfolder", "{out}/qc", "--threads", "4"],
    [o("qc/" + f, "plain", "E0") for f in QC_FILES],
    "As above at 1 Mb with the output name relative to the work directory, so QC.log compares "
    "byte for byte.", tier=4, threads_arg="--threads", large=True,
    memory={"budget_mb": BM_BUDGET})

# ---------------------------------------------------------------- hicQuickQC
add("hicQuickQC", "SRR927086_2M.HindIII.lines1000000",
    ["-s", BAM1, BAM2, "--QCfolder", "{out}/qc", "-seq", "AAGCTT", "--danglingSequence", "AGCT",
     "--restrictionCutFile", HINDIII, "--lines", "1000000"],
    [o("qc/" + f, "plain", "E0", options={"normalise": ["quickqc_temporary_matrix"]})
     for f in QC_FILES]
    + [o("qc/hicQC.html", "text", "E0",
         options={"normalise": ["quickqc_temporary_name", "pandas_styler_uuid"]})],
    BM_NOTE + "The first 1,000,000 pairs (--lines), the QC log and tables with the temporary "
    "matrix name normalised, as in the repository case. Budget 192 MB as in the repository "
    "hicQuickQC cases.", tier=4, large=True, memory={"budget_mb": 192})

# ---------------------------------------------------------------- hicQC / hicPrepareQCreport
QC_TABLES = ["QC_table.txt", "unmapable_table.txt", "discarded_table.txt", "distance_table.txt",
             "read_orientation_table.txt"]
QC_PNG = ["pairs_sequenced.png", "unmappable_and_non_unique.png", "pairs_discarded.png",
          "distance.png", "read_orientation.png"]
QCLOG = f"{DM}/qc_SRR927086_2M_HindIII/QC.log"
QCLOG_MC = f"{DM}/qc_SRR927086_2M_microc/QC.log"
QCLOG2 = f"{DM}/qc_SRR927086_second2M_HindIII/QC.log"
QC_NOTE = ("QC.log of the 3.7.6 hicBuildMatrix on the first 2,000,000 SRR927086 pairs against hg19 "
           "HindIII at 100 kb (derive_matrix.sh step 8)")
for tool in ("hicQC", "hicPrepareQCreport"):
    add(tool, "SRR927086_two_HindIII_samples",
        ["--logfiles", QCLOG, QCLOG2, "--labels", "pairs_1-2M", "pairs_2-4M",
         "--outputFolder", "{out}/qc"],
        [o("qc/" + f, "text", "E0") for f in QC_TABLES] + [o("qc/" + f, "png", "E6") for f in QC_PNG]
        + [o("qc/hicQC.html", "text", "E0", options={"normalise": ["pandas_styler_uuid"]})],
        QC_NOTE + ", and the log of read pairs 2,000,001 to 4,000,000 (step 10).",
        tier=7, memory={"alpha": 0.0, "beta": 0.0})
    add(tool, "SRR927086_2M.HindIII_and_microc_refused",
        ["--logfiles", QCLOG, QCLOG_MC, "--labels", "HindIII", "MicroC_mode",
         "--outputFolder", "{out}/qc"], [], QC_NOTE + ", and the 3.7.6 hicBuildMatrixMicroC log "
        "of the same reads (step 8). The Micro-C log lacks the restriction-site lines, and 3.7.6 "
        "stops with pandas' 'All arrays must be of the same length'; the C++ tool must exit 1 "
        "as well.", tier=7, expect_exit=1, memory={"alpha": 0.0, "beta": 0.0})
    add(tool, "SRR927086_2M.prettyReport",
        ["--logfiles", QCLOG, QCLOG2, "--labels", "pairs_1-2M", "pairs_2-4M",
         "--outputFolder", "{out}/qc", "--prettyReport"],
        [o("qc/hicQC_report.html", "html", "EX")] + [o("qc/QC_table.txt", "text", "EX")],
        QC_NOTE + ". C++ only (--prettyReport): the card report is written and every PASS/WARN/"
        "FAIL badge equals the status recomputed from QC_table.txt with the default cutoffs "
        "and the 20 % WARN margin (validator qc_pretty_report). Budget 256 MB: the run also draws "
        "the five charts in the drawing interpreter (about 160 MB measured), and a case without "
        "a reference has no Python peak to pair it with.",
        tier=7, reference="none", memory={"budget_mb": 256},
        validators=[{"name": "qc_pretty_report", "output": "{out}/qc", "options": {"folder": "qc"}}])
    add(tool, "SRR927086_2M.prettyReport.strict_cutoffs",
        ["--logfiles", QCLOG, "--outputFolder", "{out}/qc", "--prettyReport",
         "--minValidPairsPercent", "60", "--maxInterChromosomalPercent", "5",
         "--maxDuplicatedPairsPercent", "1", "--maxDanglingEndPercent", "1",
         "--maxSelfLigationPercent", "1"],
        [o("qc/hicQC_report.html", "html", "EX")],
        QC_NOTE + ". C++ only: cutoffs moved so that the real values land on the WARN and FAIL "
        "sides; the badges must follow the given cutoffs.",
        tier=7, reference="none", memory={"budget_mb": 256},
        validators=[{"name": "qc_pretty_report", "output": "{out}/qc",
                     "options": {"folder": "qc", "thresholds": {
                         "valid": 60.0, "inter": 5.0, "duplicated": 1.0, "dangling": 1.0,
                         "self_ligation": 1.0}}}])

for tool in ("hicQC", "hicPrepareQCreport"):
    add(tool, "SRR927086_2M.prettyReport.warn_band",
        ["--logfiles", QCLOG, "--outputFolder", "{out}/qc", "--prettyReport",
         "--minValidPairsPercent", "70", "--maxInterChromosomalPercent", "60",
         "--maxDanglingEndPercent", "1.7", "--maxSelfLigationPercent", "1.5"],
        [o("qc/hicQC_report.html", "html", "EX")],
        QC_NOTE + ". C++ only: cutoffs placed so that valid pairs (62.3 %), inter-chromosomal "
        "(65.1 %), dangling ends (1.86 %) and self-ligation (1.73 %) fall inside the 20 % WARN "
        "band; every badge must read WARN except duplicated pairs.",
        tier=7, reference="none", memory={"budget_mb": 256},
        validators=[{"name": "qc_pretty_report", "output": "{out}/qc",
                     "options": {"folder": "qc", "thresholds": {
                         "valid": 70.0, "inter": 60.0, "dangling": 1.7, "self_ligation": 1.5}}}])

# ---------------------------------------------------------------- hicInfo
add("hicInfo", "gm12878_r8_25kb_chr19-22.h5", ["-m", R8H],
    [o("stdout.txt", "text", "E0")], N_R8 + ", as h5 (step 5).", memory={"alpha": 1.05, "beta": 0.0})
add("hicInfo", "a549_25kb.whole_genome", ["-m", "{data}/cool/a549_25kb.cool"],
    [o("stdout.txt", "text", "E0")],
    "ENCODE A549 25 kb, hg38, whole genome, 318 M stored pixels, chrM first.",
    memory={"alpha": 0.0, "beta": 0.0})
add("hicInfo", "imr90_mcool.1Mb", ["-m", "{data}/mcool/imr90.mcool::/resolutions/1000000"],
    [o("stdout.txt", "text", "E0")], "IMR90 mcool (hg38), the 1 Mb resolution group by URI.",
    memory={"alpha": 0.0, "beta": 0.0})
add("hicInfo", "hff_microc_mcool.25kb", ["-m", "{data}/mcool/hff_microc.mcool::/resolutions/25000"],
    [o("stdout.txt", "text", "E0")], "4DN HFF Micro-C mcool, 25 kb group.",
    memory={"alpha": 0.0, "beta": 0.0})

# ---------------------------------------------------------------- hicConvertFormat
REF = "cpp/scripts/py_hic_reference.py"
add("hicConvertFormat", "hic_to_cool.gm12878_combined30.1Mb",
    ["--matrices", HIC_GM, "--outFileName", "{out}/matrix.cool", "--inputFormat", "hic",
     "--outputFormat", "cool", "--resolutions", "1000000"],
    [o("matrix_1000000.cool", "cool", "E1")],
    "GEO GSE63525 GM12878 in situ combined_30 .hic (hg19, 40 GB, chromosome names without chr), "
    "1 Mb, against hic2cool through the 3.7.6 tool.", large=True)
add("hicConvertFormat", "hic_to_cool.k562_combined30.500kb_and_2500kb",
    ["--matrices", HIC_K562, "--outFileName", "{out}/matrix.cool", "--inputFormat", "hic",
     "--outputFormat", "cool", "--resolutions", "500000", "2500000"],
    [o("matrix_500000.cool", "cool", "E1"), o("matrix_2500000.cool", "cool", "E1")],
    "GEO GSE63525 K562 combined_30 .hic (hg19), 500 kb and 2.5 Mb.", large=True)
for fmt, outs, extra in (
        ("h5", [o("matrix.h5", "h5", "E2")], []),
        ("ginteractions", [o("matrix.ginteractions.tsv", "text", "E0")], []),
        ("hicpro", [o("matrix.hicpro", "text", "E0"), o("matrix.bed", "text", "E0")],
         ["--bedFileHicpro", "{out}/matrix.bed"])):
    name = "matrix." + fmt
    add("hicConvertFormat", f"hic_to_{fmt}.gm12878_combined30.1Mb",
        ["--matrices", HIC_GM, "--outFileName", "{out}/" + name, "--inputFormat", "hic",
         "--outputFormat", fmt, "--resolutions", "1000000"] + extra,
        outs,
        "GSE63525 GM12878 combined_30 at 1 Mb. 3.7.6 cannot write this direction from .hic; the "
        "reference is the two-step route of py_hic_reference.py (hic2cool, then the 3.7.6 "
        "hicConvertFormat --inputFormat cool).", py_script=REF, large=True)
add("hicConvertFormat", "hic_to_mcool.gm12878_combined30.1Mb_2500kb",
    ["--matrices", HIC_GM, "--outFileName", "{out}/matrix.mcool", "--inputFormat", "hic",
     "--outputFormat", "mcool", "--resolutions", "1000000", "2500000"],
    [o("matrix.mcool", "mcool", "EX")],
    "v4 only: 3.7.6 has no hic -> mcool with --resolutions. Each resolution group must hold the "
    "same bins and bit-identical pixels as hic2cool's own conversion of the same file "
    "(derive_matrix.sh step 6, written by the 3.7.6 hicConvertFormat).",
    reference="none", large=True,
    validators=[{"name": "cool_pixels_match", "output": "{out}/matrix.mcool",
                 "options": {"expected": f"{DM}/gm12878_combined30_1000000.cool",
                             "uri": "::/resolutions/1000000"}},
                {"name": "cool_pixels_match", "output": "{out}/matrix.mcool",
                 "options": {"expected": f"{DM}/gm12878_combined30_2500000.cool",
                             "uri": "::/resolutions/2500000"}}])
add("hicConvertFormat", "cool_to_hic.v8.gm12878_r8_25kb_chr19-22",
    ["--matrices", R8, "--outFileName", "{out}/matrix.hic", "--inputFormat", "cool",
     "--outputFormat", "hic", "--resolutions", "50000", "100000", "250000"],
    [o("matrix.hic", "hic", "E2")],
    N_R8 + ". .hic version 8 at 25, 50, 100 and 250 kb; comparators/hic.py reads the C++ .hic "
    "back and requires the pixels the 3.7.6 loader sees (py_hic_reference.py) bit for bit, plus "
    "full-length normalization vectors: the pixel round trip.", py_script=REF)
add("hicConvertFormat", "cool_to_hic.v9.gm12878_r8_10kb_chr21-22",
    ["--matrices", R8_10, "--outFileName", "{out}/matrix.hic", "--inputFormat", "cool",
     "--outputFormat", "hic", "--hicVersion", "9", "--resolutions", "50000", "--threads", "4"],
    [o("matrix.hic", "hic", "E2")],
    "ENCODE ENCSR946AKI GM12878 r8 10 kb, chr21 and chr22 (step 1), .hic version 9 at 10 and "
    "50 kb, checked as above.", py_script=REF)
add("hicConvertFormat", "cool_to_hic.v9.weighted_whole_genome_r12_25kb",
    ["--matrices", "{data}/cool/gm12878_ENCSR946AKI_r12_25kb.cool", "--outFileName",
     "{out}/matrix.hic", "--inputFormat", "cool", "--outputFormat", "hic", "--hicVersion", "9",
     "--load_raw_values", "--resolutions", "1000000", "--threads", "4"],
    [o("matrix.hic", "hic", "E2")],
    "The whole-genome ENCODE r12 25 kb cool as published (hic2cool output with SCALE/VC/weight "
    "columns, 25 chromosomes), raw counts, version 9 at 25 kb and 1 Mb.", py_script=REF, large=True)
add("hicConvertFormat", "cool_to_h5.gm12878_r8_25kb_chr19-22",
    ["--matrices", R8, "--outFileName", "{out}/matrix.h5", "--inputFormat", "cool",
     "--outputFormat", "h5"], [o("matrix.h5", "h5", "E2")], N_R8 + ".")
add("hicConvertFormat", "h5_to_cool.a549_25kb_chr19-22",
    ["--matrices", A5H, "--outFileName", "{out}/matrix.cool", "--inputFormat", "h5",
     "--outputFormat", "cool"], [o("matrix.cool", "cool", "E1")],
    "ENCODE A549 25 kb chr19 to chr22 as h5 (3.7.6 cool -> h5, step 5), back to cool.")
add("hicConvertFormat", "cool_to_ginteractions.gm12878_r12_25kb_chr19-22",
    ["--matrices", R12, "--outFileName", "{out}/matrix.ginteractions", "--inputFormat", "cool",
     "--outputFormat", "ginteractions"], [o("matrix.ginteractions.tsv", "text", "E0")],
    SUB25.format(r=12) + ".")
add("hicConvertFormat", "cool_to_homer.gm12878_r8_25kb_chr22",
    ["--matrices", R8, "--outFileName", "{out}/matrix.homer", "--inputFormat", "cool",
     "--outputFormat", "homer", "--chromosome", "chr22"], [o("matrix.homer", "text", "E0")],
    N_R8 + ", --chromosome chr22 (2,033 bins, dense table).")
add("hicConvertFormat", "cool_to_homer.gm12878_r8_25kb_chr22_subset",
    ["--matrices", f"{DM}/gm12878_ENCSR946AKI_r8_25kb_chr22.cool", "--outFileName",
     "{out}/matrix.homer", "--inputFormat", "cool", "--outputFormat", "homer"],
    [o("matrix.homer", "text", "E0")],
    "ENCODE GM12878 r8 25 kb, chr22 only (2,033 bins, derive_matrix.sh step 1), dense homer table.")
add("hicConvertFormat", "cool_to_hicpro.a549_25kb_chr19-22",
    ["--matrices", A5, "--outFileName", "{out}/matrix.hicpro", "--inputFormat", "cool",
     "--outputFormat", "hicpro", "--bedFileHicpro", "{out}/matrix.bed"],
    [o("matrix.hicpro", "text", "E0"), o("matrix.bed", "text", "E0")],
    "ENCODE A549 25 kb chr19 to chr22 (step 1), 6.1 M pixels.")
add("hicConvertFormat", "cool_to_mcool.gm12878_r8_25kb_chr19-22",
    ["--matrices", R8, "--outFileName", "{out}/matrix.mcool", "--inputFormat", "cool",
     "--outputFormat", "mcool", "-r", "25000", "50000", "100000", "1000000"],
    [o("matrix.mcool", "mcool", "E1")], N_R8 + ", coarsened to 50 kb, 100 kb and 1 Mb.")
add("hicConvertFormat", "mcool_group_to_cool.imr90_1Mb",
    ["--matrices", "{data}/mcool/imr90.mcool::/resolutions/1000000", "--outFileName",
     "{out}/matrix.cool", "--inputFormat", "cool", "--outputFormat", "cool"],
    [o("matrix.cool", "cool", "E1")],
    "IMR90 mcool (hg38, cooler 0.7.6 file), the balanced 1 Mb group read by URI and written as cool.")
add("hicConvertFormat", "mcool_group_to_cool.imr90_250kb",
    ["--matrices", "{data}/mcool/imr90.mcool::/resolutions/250000", "--outFileName",
     "{out}/matrix.cool", "--inputFormat", "cool", "--outputFormat", "cool"],
    [o("matrix.cool", "cool", "E1")],
    "IMR90 mcool (hg38), the balanced 250 kb group read by URI and written as cool. This group "
    "stores 34,561 duplicate (bin1, bin2) pixels; 3.7.6 sums them (scipy COO to CSR).", large=True)

# ---------------------------------------------------------------- hicMergeMatrixBins
MERGE_NOTE = (" Reference none: the C++ merge groups follow --chromosomeSizes (required by the "
              "port only) and keep every bin, a deliberate deviation from 3.7.6, whose merge on "
              "this real matrix yields {py} bins instead of {cpp} (measured: the 3.7.6 output "
              "drops bins). The C++ output must instead hold the bins and bit-identical pixels "
              "of cooler coarsen of the same input (derive_matrix.sh step 9).")
add("hicMergeMatrixBins", "gm12878_r8_25kb_chr19-22.numBins4",
    ["-m", R8, "-o", "{out}/merged.cool", "--numBins", "4"],
    [o("merged.cool", "cool", "EX")],
    N_R8 + ", merged to 100 kb." + MERGE_NOTE.format(py="1,926", cpp="2,209"),
    cpp_args=["--chromosomeSizes", SIZES], reference="none",
    validators=[{"name": "cool_pixels_match", "output": "{out}/merged.cool",
                 "options": {"expected": f"{DM}/gm12878_ENCSR946AKI_r8_25kb_chr19-22.coarsen4.cool"}}])
add("hicMergeMatrixBins", "a549_25kb_chr19-22.numBins10",
    ["-m", A5, "-o", "{out}/merged.cool", "--numBins", "10"],
    [o("merged.cool", "cool", "EX")],
    "ENCODE A549 25 kb chr19 to chr22, merged to 250 kb."
    + MERGE_NOTE.format(py="754", cpp="884"),
    cpp_args=["--chromosomeSizes", SIZES], reference="none",
    validators=[{"name": "cool_pixels_match", "output": "{out}/merged.cool",
                 "options": {"expected": f"{DM}/a549_25kb_chr19-22.coarsen10.cool"}}])
add("hicMergeMatrixBins", "gm12878_r12_25kb_chr19-22.runningWindow3",
    ["-m", R12, "-o", "{out}/merged.cool", "--numBins", "3", "--runningWindow"],
    [o("merged.cool", "cool", "E1")], SUB25.format(r=12) + ", running window of 3 bins.",
    cpp_args=["--chromosomeSizes", SIZES])

# ---------------------------------------------------------------- hicSumMatrices
add("hicSumMatrices", "gm12878_r8_r12_25kb_chr19-22",
    ["-m", R8, R12, "-o", "{out}/sum.cool"], [o("sum.cool", "cool", "E1")],
    "ENCODE ENCSR946AKI GM12878 replicates 8 + 12, 25 kb, chr19 to chr22 (step 1).")
add("hicSumMatrices", "gm12878_r8_r12_25kb_chr19-22.h5",
    ["-m", R8H, R12H, "-o", "{out}/sum.h5"], [o("sum.h5", "h5", "E2")],
    "The same as h5 (step 5).")
add("hicSumMatrices", "gm12878_r8_r12_10kb.whole_genome",
    ["-m", "{data}/cool/gm12878_ENCSR946AKI_r8_10kb.cool",
     "{data}/cool/gm12878_ENCSR946AKI_r12_10kb.cool", "-o", "{out}/sum.cool"],
    [o("sum.cool", "cool", "E1")],
    "ENCODE ENCSR946AKI GM12878 replicates 8 + 12 as published, 10 kb, whole genome (hg38, "
    "161 M stored pixels in r8; SCALE/VC/weight bin columns).", large=True)

# ---------------------------------------------------------------- hicNormalize
add("hicNormalize", "norm_range.gm12878_r8_r12",
    ["--matrices", R8, R12, "--normalize", "norm_range", "-o", "{out}/r8.cool", "{out}/r12.cool"],
    [o("r8.cool", "cool", "E1"), o("r12.cool", "cool", "E1")], N_SUBSETS + ".")
add("hicNormalize", "smallest.gm12878_r8_r12_a549",
    ["--matrices", R8, R12, A5, "--normalize", "smallest", "-o", "{out}/r8.cool",
     "{out}/r12.cool", "{out}/a549.cool"],
    [o("r8.cool", "cool", "E1"), o("r12.cool", "cool", "E1"), o("a549.cool", "cool", "E1")],
    N_SUBSETS + "; A549 (87 M contacts) scaled down to r8 (10 M).")
add("hicNormalize", "multiplicative.a549.h5",
    ["--matrices", A5H, "--normalize", "multiplicative", "--multiplicativeValue", "0.25",
     "-o", "{out}/a549.h5"], [o("a549.h5", "h5", "E2")],
    "ENCODE A549 25 kb chr19 to chr22 as h5, times 0.25.")
add("hicNormalize", "smallest.setToZeroThreshold.r8_10kb_a549_10kb",
    ["--matrices", R8_10, A5_10, "--normalize", "smallest", "--setToZeroThreshold", "1.0",
     "-o", "{out}/r8.cool", "{out}/a549.cool"],
    [o("r8.cool", "cool", "E1"), o("a549.cool", "cool", "E1")],
    "ENCODE GM12878 r8 and A549, 10 kb, chr21 and chr22 (step 1), smallest with values below 1 "
    "set to zero.")

# ---------------------------------------------------------------- hicCorrectMatrix
add("hicCorrectMatrix", "diagnostic_plot.gm12878_r8_25kb_chr19-22",
    ["diagnostic_plot", "--matrix", R8, "--plotName", "{out}/diagnostic.png"],
    [o("diagnostic.png", "png", "E6")], N_R8 + ".", tier=7)
add("hicCorrectMatrix", "diagnostic_plot.a549.perchr",
    ["diagnostic_plot", "--matrix", A5, "--perchr", "--plotName", "{out}/diagnostic.png"],
    [o("diagnostic.png", "png", "E6")], "ENCODE A549 25 kb chr19 to chr22, per chromosome.", tier=7)
add("hicCorrectMatrix", "ICE.gm12878_r8_25kb_chr19-22",
    ["correct", "-m", R8H, "--correctionMethod", "ICE", "--filterThreshold", "-1.5", "5.0",
     "--iterNum", "500", "--filteredBed", "{out}/filtered.bed", "-o", "{out}/ice.h5"],
    [o("ice.h5", "h5", "E3"), o("filtered.bed", "bed", "E0")],
    N_R8 + ", as h5 (step 5); ICE with the usual -1.5 5 MAD filter.", tier=3,
    threads_arg="--threads")
add("hicCorrectMatrix", "KR.gm12878_r8_25kb_chr19-22",
    ["correct", "-m", R8H, "--correctionMethod", "KR", "-o", "{out}/kr.h5"],
    [o("kr.h5", "h5", "ED")], N_R8 + ", as h5; Knight-Ruiz.", tier=3, threads_arg="--threads")
add("hicCorrectMatrix", "KR.perchr.a549_25kb_chr19-22",
    ["correct", "-m", A5H, "--correctionMethod", "KR", "--perchr", "-o", "{out}/kr.h5"],
    [o("kr.h5", "h5", "ED")], "ENCODE A549 25 kb chr19 to chr22 as h5; KR per chromosome, "
    "default --compatMode v4 (float64 balancing).",
    tier=3, threads_arg="--threads")
add("hicCorrectMatrix", "KR.perchr.a549_25kb_chr19-22.compatMode_v3",
    ["correct", "-m", A5H, "--correctionMethod", "KR", "--perchr", "-o", "{out}/kr.h5"],
    [o("kr.h5", "h5", "ED")], "As the case above with the C++-only --compatMode v3 (cpp_args), "
    "which reproduces krbalancing's float32 input rounding and float32 rescaling accumulators.",
    tier=3, threads_arg="--threads", cpp_args=["--compatMode", "v3"])
add("hicCorrectMatrix", "ICE.a549_25kb_chr19-22.cool_input",
    ["correct", "-m", A5, "--correctionMethod", "ICE", "--filterThreshold", "-2", "4",
     "-o", "{out}/ice.h5"], [o("ice.h5", "cool", "E3")],
    "ENCODE A549 25 kb chr19 to chr22 read from cool, ICE, output named .h5: both tools reuse "
    "the cool handler of the load and write a cool file under the .h5 name.", tier=3,
    threads_arg="--threads")

# ---------------------------------------------------------------- hicAdjustMatrix
REG = f"{DM}/hg38_chr19-21_regions.bed"
add("hicAdjustMatrix", "keep.chromosomes.reordered",
    ["-m", R8, "-o", "{out}/adjusted.cool", "--chromosomes", "chr22", "chr19", "--action", "keep"],
    [o("adjusted.cool", "cool", "E1")], N_R8 + "; keep chr22 and chr19 in that order.")
add("hicAdjustMatrix", "remove.chromosomes",
    ["-m", A5, "-o", "{out}/adjusted.cool", "--chromosomes", "chr20", "--action", "remove"],
    [o("adjusted.cool", "cool", "E1")], "ENCODE A549 25 kb chr19 to chr22; remove chr20.")
add("hicAdjustMatrix", "mask.chromosomes",
    ["-m", R12, "-o", "{out}/adjusted.cool", "--chromosomes", "chr21", "--action", "mask"],
    [o("adjusted.cool", "cool", "E1")], SUB25.format(r=12) + "; mask chr21.")
for action in ("keep", "remove", "mask"):
    add("hicAdjustMatrix", f"{action}.regions.h5",
        ["-m", R8H, "-o", "{out}/adjusted.h5", "--regions", REG, "--action", action],
        [o("adjusted.h5", "h5", "E2")],
        N_R8 + ", as h5; --regions: the chr19 and chr20 centromeres with flanks and the chr21 "
        "p arm (step 7).")
add("hicAdjustMatrix", "keep.chromosome.whole_genome_a549",
    ["-m", "{data}/cool/a549_25kb.cool", "-o", "{out}/adjusted.cool", "--chromosomes", "chr1",
     "--action", "keep"], [o("adjusted.cool", "cool", "E1")],
    "ENCODE A549 25 kb whole genome (chrM first, weight column), keep chr1.", large=True)

# ---------------------------------------------------------------- hicTransform
for method in ("obs_exp", "obs_exp_non_zero", "obs_exp_lieberman"):
    add("hicTransform", f"{method}.gm12878_r8_25kb_chr19-22",
        ["-m", R8, "-o", f"{{out}}/{method}.cool", "--method", method],
        [o(f"{method}.cool", "cool", "ED")], N_R8 + ".")
    add("hicTransform", f"{method}.perChromosome.a549.h5",
        ["-m", A5H, "-o", f"{{out}}/{method}.h5", "--method", method, "--perChromosome"],
        [o(f"{method}.h5", "h5", "ED")], "ENCODE A549 25 kb chr19 to chr22 as h5, per chromosome.")
add("hicTransform", "obs_exp_non_zero.ligation_factor.r12.h5",
    ["-m", R12H, "-o", "{out}/oe.h5", "--method", "obs_exp_non_zero", "--ligation_factor"],
    [o("oe.h5", "h5", "ED")], SUB25.format(r=12) + ", as h5.")
add("hicTransform", "pearson.gm12878_r8_25kb_chr19-22.h5",
    ["-m", R8H, "-o", "{out}/pearson.h5", "--method", "pearson"],
    [o("pearson.h5", "h5", "E3")], N_R8 + ", as h5. E3, not ED: about ten of 29.6 M "
    "correlations are cancellation residues of order 1e-17 whose sign and size depend on the "
    "summation order (measured 1e-17 vs -3e-21), which a pure relative tolerance cannot "
    "accept; the largest absolute difference is 2.8e-14.")
add("hicTransform", "pearson.perChromosome.a549.h5",
    ["-m", A5H, "-o", "{out}/pearson.h5", "--method", "pearson", "--perChromosome"],
    [o("pearson.h5", "h5", "E3")], "ENCODE A549 25 kb chr19 to chr22 as h5.")
add("hicTransform", "covariance.gm12878_r12_25kb_chr19-22.h5",
    ["-m", R12H, "-o", "{out}/covariance.h5", "--method", "covariance"],
    [o("covariance.h5", "h5", "E3")], SUB25.format(r=12) + ", as h5. E3, not ED: as for pearson, "
    "12 of 29.5 M values are cancellation residues near 1e-18 (largest absolute difference "
    "9.5e-12).")
add("hicTransform", "covariance.perChromosome.r8_10kb",
    ["-m", R8_10, "-o", "{out}/covariance.cool", "--method", "covariance", "--perChromosome"],
    [o("covariance.cool", "cool", "E3")],
    "ENCODE GM12878 r8 10 kb, chr21 and chr22 (step 1).")

# ---------------------------------------------------------------- hicCompareMatrices
for op, cls in (("log2ratio", "ED"), ("ratio", "E1"), ("diff", "E1")):
    add("hicCompareMatrices", f"{op}.gm12878_r8_vs_a549",
        ["-m", R8, A5, "--operation", op, "-o", "{out}/result.cool"],
        [o("result.cool", "cool", cls)], "GM12878 r8 against A549, 25 kb, chr19 to chr22 (step 1).")
add("hicCompareMatrices", "log2ratio.gm12878_r8_vs_r12.h5.noNorm",
    ["-m", R8H, R12H, "--operation", "log2ratio", "--noNorm", "-o", "{out}/result.h5"],
    [o("result.h5", "h5", "ED")], "GM12878 r8 against r12 as h5, without depth normalisation.")

# ---------------------------------------------------------------- hicCorrelate
THREE = [R8, R12, A5]
add("hicCorrelate", "pearson.log1p.range.r8_r12_a549",
    ["--matrices"] + THREE + ["--labels", "GM12878_r8", "GM12878_r12", "A549", "--method",
                              "pearson", "--log1p", "--range", "50000:2000000",
                              "--outFileNameHeatmap", "{out}/heatmap.png",
                              "--outFileNameScatter", "{out}/scatter.png", "--plotNumbers"],
    [o("heatmap.png", "png", "E6"), o("scatter.png", "png", "E6")], N_SUBSETS + ".", tier=7)
add("hicCorrelate", "spearman.r8_r12_a549",
    ["--matrices"] + THREE + ["--labels", "GM12878_r8", "GM12878_r12", "A549", "--method",
                              "spearman", "--outFileNameHeatmap", "{out}/heatmap.png",
                              "--outFileNameScatter", "{out}/scatter.png"],
    [o("heatmap.png", "png", "E6"), o("scatter.png", "png", "E6")], N_SUBSETS + ".", tier=7)
add("hicCorrelate", "pearson.chromosomes.10kb",
    ["--matrices", R8_10, R12_10, A5_10, "--method", "pearson", "--log1p", "--chromosomes",
     "chr21", "--range", "10000:1000000", "--outFileNameHeatmap", "{out}/heatmap.png",
     "--outFileNameScatter", "{out}/scatter.png"],
    [o("heatmap.png", "png", "E6"), o("scatter.png", "png", "E6")],
    "GM12878 r8, r12 and A549 at 10 kb, chr21 only.", tier=7)
V4_PARAMS = {
    "scc": ["--sccSmoothWindow", "10", "--sccMaxDistance", "5000000"],
    "genomedisco": [],
    "hicspector": [],
    "quasar": [],
}
for method, params in V4_PARAMS.items():
    add("hicCorrelate", f"{method}.r8_r12_a549",
        ["--matrices"] + THREE + ["--labels", "GM12878_r8", "GM12878_r12", "A549", "--method",
                                  method] + params + ["--outFileNameHeatmap", "{out}/heatmap.png",
                                  "--outFileNameScatter", "{out}/scatter.png",
                                  "--plotData", "{out}/plot.json"],
        [o("plot.json", "json", "EX")],
        N_SUBSETS + f". v4 only (--method {method}, no 3.7.6 counterpart): the GM12878 replicate "
        "pair must score above both GM12878-vs-A549 pairs (validator plotdata_ranking)."
        + (" hicrep's suggested h=10 for 25 kb and dBPMax 5 Mb; the C++ scores match "
           "hicrep 0.2.6 run on the same files (mean over chromosomes, checked by hand)."
           if method == "scc" else ""),
        tier=7, reference="none",
        validators=[{"name": "plotdata_ranking", "output": "{out}/plot.json",
                     "options": {"higher": [[0, 1]], "lower": [[0, 2], [1, 2]]}}])

# ---------------------------------------------------------------- hicPlotMatrix
add("hicPlotMatrix", "region.log1p.gm12878_r8",
    ["--matrix", R8, "--region", "chr21:30000000-36000000", "--log1p",
     "--outFileName", "{out}/matrix.png"], [o("matrix.png", "png", "E6")], N_R8 + ".", tier=7)
add("hicPlotMatrix", "region_region2.a549.h5",
    ["--matrix", A5H, "--region", "chr19:40000000-45000000", "--region2", "chr19:42000000-47000000",
     "--log", "--outFileName", "{out}/matrix.png", "--colorMap", "hot_r"],
    [o("matrix.png", "png", "E6")], "ENCODE A549 25 kb as h5, two regions, log.", tier=7)
add("hicPlotMatrix", "perChr.gm12878_r12",
    ["--matrix", R12, "--perChromosome", "--log1p", "--outFileName", "{out}/matrix.png"],
    [o("matrix.png", "png", "E6")], SUB25.format(r=12) + ", every chromosome.", tier=7)
add("hicPlotMatrix", "region.weighted_whole_genome_10kb",
    ["--matrix", "{data}/cool/gm12878_ENCSR946AKI_r8_10kb.cool", "--region",
     "chr2:100000000-102000000", "--log1p", "--outFileName", "{out}/matrix.png"],
    [o("matrix.png", "png", "E6")],
    "The published whole-genome ENCODE r8 10 kb cool (weight column), a 2 Mb region of chr2.",
    tier=7)
add("hicPlotMatrix", "matrix2.gm12878_r8_vs_a549",
    ["--matrix", R8, "--matrix2", A5, "--region", "chr21:30000000-36000000", "--log1p",
     "--outFileName", "{out}/matrix.png"], [o("matrix.png", "png", "EX")],
    "v4 only (--matrix2): GM12878 r8 in the upper and A549 in the lower triangle of one heatmap; "
    "the figure must be written.", tier=7, reference="none")

# ---------------------------------------------------------------- hicPlotDistVsCounts
add("hicPlotDistVsCounts", "perchr.r8_r12_a549",
    ["--matrices"] + THREE + ["--labels", "GM12878_r8", "GM12878_r12", "A549", "--perchr",
                              "--plotFile", "{out}/plot.png", "--outFileData", "{out}/data.txt",
                              "--maxdepth", "10000000"],
    [o("plot.png", "png", "E6"), o("data.txt", "text", "E0")], N_SUBSETS + ".", tier=7)
add("hicPlotDistVsCounts", "genome.skipDiagonal.10kb",
    ["--matrices", R8_10, A5_10, "--plotFile", "{out}/plot.png", "--skipDiagonal",
     "--outFileData", "{out}/data.txt", "--maxdepth", "3000000", "--chromosomeExclude", "chr22"],
    [o("plot.png", "png", "E6"), o("data.txt", "text", "E0")],
    "GM12878 r8 and A549 at 10 kb, chr21 (chr22 excluded).", tier=7)

for tool, payload in cases.items():
    path = OUT / f"{tool}.json"
    if tool == "hicInfo" and path.exists():
        existing = [c for c in json.loads(path.read_text()) if c["id"] == "ext.hicInfo.gm12878_r8_25kb"]
        payload = existing + payload
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print(tool, len(payload))
