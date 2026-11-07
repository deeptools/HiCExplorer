#!/bin/bash
# Wall-clock time and peak memory of HiCExplorer 4.0 against 3.7.6, the numbers of
# docs/benchmarks.md. Every run is timed with /usr/bin/time -v and its log is kept.
#
# Required environment:
#   OUT        output folder; one <run>.log per run plus summary.tsv
#   CPP        folder with the HiCExplorer 4 tools (e.g. cpp/build/tools)
#   PY376      bin folder of a HiCExplorer 3.7.6 environment
#   SAM1 SAM2  the two mate files for hicBuildMatrix
#   CUTFILE    restriction cut file (DpnII) for hicBuildMatrix
#   CHR1       single-chromosome 10 kb cool file (loops, TADs, PCA)
#   CHINPUT RMAP BAITMAP NPB NBPB POE   CHiCAGO design files and .chinput file
# Optional:
#   STEPS      subset of: build loops tads pca chicago (default: all)
#   REPS       repetitions of the short runs (default 3)
#   MAXLOAD    wait until the 1-minute load average is below this (default 4)
set -euo pipefail
: "${OUT:?}" "${CPP:?}" "${PY376:?}"
STEPS=${STEPS:-"build loops tads pca chicago"}
REPS=${REPS:-3}
MAXLOAD=${MAXLOAD:-4}
mkdir -p "$OUT/work"
W=$OUT/work

wait_quiet() {
    while :; do
        load=$(cut -d' ' -f1 /proc/loadavg)
        awk -v l="$load" -v m="$MAXLOAD" 'BEGIN{exit !(l < m)}' && break
        sleep 60
    done
}

# run <name> <command...>: one timed run, log in $OUT/<name>.log, a row in summary.tsv
run() {
    local name=$1; shift
    wait_quiet
    echo "# start $(date -Is) load $(cut -d' ' -f1-3 /proc/loadavg)" > "$OUT/$name.log"
    /usr/bin/time -v "$@" >> "$OUT/$name.log" 2>&1
    local wall rss
    wall=$(grep 'Elapsed (wall clock)' "$OUT/$name.log" | awk '{print $NF}')
    rss=$(grep 'Maximum resident set size' "$OUT/$name.log" | awk '{print $NF}')
    printf '%s\t%s\t%s\n' "$name" "$wall" "$rss" >> "$OUT/summary.tsv"
}

[ -f "$OUT/summary.tsv" ] || printf 'run\twall_clock\tmax_rss_kib\n' > "$OUT/summary.tsv"
has() { [[ " $STEPS " == *" $1 "* ]]; }

if has build; then
    : "${SAM1:?}" "${SAM2:?}" "${CUTFILE:?}"
    for t in 1 4 8; do
        args=(--samFiles "$SAM1" "$SAM2" --restrictionCutFile "$CUTFILE"
              --restrictionSequence GATC --danglingSequence GATC --binSize 10000
              --threads "$t")
        run "build_v4_t$t" "$CPP/hicBuildMatrix" "${args[@]}" --noPlot \
            --outFileName "$W/build_v4_t$t.cool" --QCfolder "$W/qc_v4_t$t"
        run "build_376_t$t" "$PY376/hicBuildMatrix" "${args[@]}" \
            --outFileName "$W/build_376_t$t.cool" --QCfolder "$W/qc_376_t$t"
    done
fi

# Loops and TADs run on the corrected 10 kb matrix, PCA on the corrected 100 kb matrix.
# Both versions read the same files, prepared once with HiCExplorer 4.
if has loops || has tads || has pca; then
    : "${CHR1:?}"
    [ -f "$W/chr1_10kb_corrected.cool" ] || "$CPP/hicCorrectMatrix" correct \
        --matrix "$CHR1" --outFileName "$W/chr1_10kb_corrected.cool" > "$W/prep.log" 2>&1
    "$PY376/python" -c 'import sys, h5py
f = h5py.File(sys.argv[1], "r")
for n, l in zip(f["chroms/name"][:], f["chroms/length"][:]):
    print(n.decode(), l, sep="\t")' "$CHR1" > "$W/chrom.sizes"
    [ -f "$W/chr1_100kb.cool" ] || "$CPP/hicMergeMatrixBins" --matrix "$CHR1" \
        --numBins 10 --chromosomeSizes "$W/chrom.sizes" \
        --outFileName "$W/chr1_100kb.cool" >> "$W/prep.log" 2>&1
    [ -f "$W/chr1_100kb_corrected.cool" ] || "$CPP/hicCorrectMatrix" correct \
        --matrix "$W/chr1_100kb.cool" --outFileName "$W/chr1_100kb_corrected.cool" \
        >> "$W/prep.log" 2>&1
fi

for r in $(seq 1 "$REPS"); do
    if has loops; then
        for v in v4 376; do
            bin=$CPP; [ $v = 376 ] && bin=$PY376
            run "loops_${v}_r$r" "$bin/hicDetectLoops" --matrix "$W/chr1_10kb_corrected.cool" \
                --outFileName "$W/loops_${v}_r$r.bedgraph" --threads 1
        done
    fi
    if has tads; then
        for v in v4 376; do
            bin=$CPP; [ $v = 376 ] && bin=$PY376
            run "tads_${v}_r$r" "$bin/hicFindTADs" --matrix "$W/chr1_10kb_corrected.cool" \
                --outPrefix "$W/tads_${v}_r$r" --correctForMultipleTesting fdr \
                --numberOfProcessors 1
        done
    fi
    if has pca; then
        for v in v4 376; do
            bin=$CPP; [ $v = 376 ] && bin=$PY376
            run "pca_${v}_r$r" "$bin/hicPCA" --matrix "$W/chr1_100kb_corrected.cool" \
                --outputFileName "$W/pca_${v}_r$r.bedgraph" --whichEigenvectors 1
        done
    fi
done

if has chicago; then
    : "${CHINPUT:?}" "${RMAP:?}" "${BAITMAP:?}" "${NPB:?}" "${NBPB:?}" "${POE:?}"
    design=(--rmap "$RMAP" --baitmap "$BAITMAP")
    model=(--nperbin "$NPB" --nbaitsperbin "$NBPB" --proxOE "$POE")
    run chicago_convert "$CPP/hicConvertFormat" --matrices "$CHINPUT" --rmap "$RMAP" \
        --inputFormat chinput --outputFormat cool --outFileName "$W/chicago.cool"
    for src in chinput matrix; do
        if [ $src = chinput ]; then in=(--chinput "$CHINPUT"); else in=(--matrices "$W/chicago.cool"); fi
        run "chicago_bg_$src" "$CPP/chicChicagoBackgroundModel" "${design[@]}" "${model[@]}" \
            "${in[@]}" -o "$W/chicago_bg_$src.txt" --threads 16
        run "chicago_scores_$src" "$CPP/chicChicagoScores" "${design[@]}" "${in[@]}" \
            --backgroundModel "$W/chicago_bg_$src.txt" -o "$W/chicago_scores_$src.txt" --threads 16
    done
fi
echo "done $(date -Is)" >> "$OUT/summary.tsv.done"
