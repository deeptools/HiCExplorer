"""End to end on committed real data with the C++ tools.

hicBuildMatrix -> hicCorrectMatrix correct -> hicFindTADs / hicDetectLoops /
hicPCA through ``hicexplorer-workflow run``; then every step's logged argv is
executed by hand in a separate directory and each output is compared at class
E0 with ``cpp/scripts/equiv.py``.

Needs ``HICX_CPP_BIN`` (tools with ``--help-json``) and
``HICX_REFERENCE_PYTHON`` (a Python that can run equiv.py). The bin size is
100 kb: at 10 kb the 37,321 contacts of these BAMs are too sparse for ICE,
which stops with "matrix correction produced extremely large values".
"""

import json
import os
import subprocess
import sys

import pytest
import yaml

from conftest import GUI_DIR

REPO = os.path.dirname(GUI_DIR)
DATA = os.path.join(REPO, "hicexplorer", "test", "test_data")
EQUIV = os.path.join(REPO, "cpp", "scripts", "equiv.py")
CPP_BIN = os.environ.get("HICX_CPP_BIN")
REFERENCE_PYTHON = os.environ.get("HICX_REFERENCE_PYTHON")


def _skip_reason():
    if not CPP_BIN:
        return "HICX_CPP_BIN is not set"
    tool = os.path.join(CPP_BIN, "hicBuildMatrix")
    try:
        proc = subprocess.run([tool, "--help-json"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        spec = json.loads(proc.stdout.decode()) if proc.returncode == 0 else None
    except (OSError, ValueError, subprocess.TimeoutExpired):
        spec = None
    if not spec or spec.get("schema") != "hicexplorer-tool-spec":
        return "{} --help-json does not print a tool specification (tier 10.1 not built)".format(tool)
    if not REFERENCE_PYTHON:
        return "HICX_REFERENCE_PYTHON (a Python able to run cpp/scripts/equiv.py) is not set"
    return None


SKIP_REASON = _skip_reason()
pytestmark = pytest.mark.skipif(SKIP_REASON is not None, reason=SKIP_REASON or "")

# (step, output name) -> (equiv.py format, file names to compare; None: the output path itself)
COMPARISONS = {
    ("build", "matrix"): ("h5", None),
    ("build", "bam"): ("plain", None),
    ("build", "qc"): ("text", ["QC_table.txt", "distance_table.txt", "read_orientation_table.txt",
                               "unmapable_table.txt", "discarded_table.txt"]),
    ("correct", "matrix"): ("h5", None),
    ("tads", "t"): (None, [("text", "_tad_score.bm"), ("text", "_score.bedgraph"), ("text", "_boundaries.bed"),
                           ("text", "_boundaries.gff"), ("text", "_domains.bed"), ("h5", "_zscore_matrix.h5")]),
    ("loops", "loops"): ("text", None),
    ("pca", "pca1"): ("bedgraph", None),
    ("pca", "pca2"): ("bedgraph", None),
}


def workflow():
    return {
        "version": 1, "name": "hic-basic-e2e", "threads": 4,
        "inputs": {
            "r1": os.path.join(DATA, "small_test_R1_unsorted.bam"),
            "r2": os.path.join(DATA, "small_test_R2_unsorted.bam"),
            "dpnii": os.path.join(DATA, "DpnII.bed"),
            "genes": os.path.join(DATA, "dm3_genes.bed.gz"),
        },
        "steps": [
            {"id": "build", "tool": "hicBuildMatrix", "threads": 4,
             "args": {"samFiles": ["${inputs.r1}", "${inputs.r2}"], "binSize": [100000],
                      "restrictionCutFile": ["${inputs.dpnii}"], "restrictionSequence": ["GATC"],
                      "danglingSequence": ["GATC"], "outBam": "${outputs.bam}",
                      "QCfolder": "${outputs.qc}", "outFileName": "${outputs.matrix}"},
             "outputs": {"matrix": "build/matrix.h5", "bam": "build/valid.bam", "qc": "build/qc"}},
            {"id": "correct", "tool": "hicCorrectMatrix", "subcommand": "correct",
             "args": {"matrix": "${steps.build.outputs.matrix}", "correctionMethod": "ICE",
                      "filterThreshold": [-1.5, 5], "outFileName": "${outputs.matrix}"},
             "outputs": {"matrix": "correct/matrix.h5"}},
            {"id": "tads", "tool": "hicFindTADs", "threads": 2,
             "args": {"matrix": "${steps.correct.outputs.matrix}", "correctForMultipleTesting": "fdr",
                      "outPrefix": "${outputs.t}"},
             "outputs": {"t": "tads/t"}},
            {"id": "loops", "tool": "hicDetectLoops",
             "args": {"matrix": "${steps.correct.outputs.matrix}", "peakInteractionsThreshold": 1,
                      "pValue": 0.5, "pValuePreselection": 0.5, "threadsPerChromosome": 1,
                      "outFileName": "${outputs.loops}"},
             "outputs": {"loops": "loops/loops.bedgraph"}},
            {"id": "pca", "tool": "hicPCA",
             "args": {"matrix": "${steps.correct.outputs.matrix}",
                      "outputFileName": ["${outputs.pca1}", "${outputs.pca2}"], "format": "bedgraph",
                      "whichEigenvectors": [1, 2], "method": "lieberman", "extraTrack": "${inputs.genes}",
                      "chromosomes": ["chrX"]},
             "outputs": {"pca1": "pca/pca1.bedgraph", "pca2": "pca/pca2.bedgraph"}},
        ],
    }


def equiv(fmt, a, b):
    proc = subprocess.run([REFERENCE_PYTHON, EQUIV, "compare", "--format", fmt, "--class", "E0", a, b],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    return proc.returncode, proc.stdout


def cli(*args):
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    return subprocess.run([sys.executable, "-m", "hicexplorer_gui.workflow.cli"] + list(args),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env)


def test_hic_workflow_outputs_match_logged_command_lines(tmp_path):
    work = tmp_path / "workflow"
    work.mkdir()
    path = str(work / "hic.yaml")
    with open(path, "w") as handle:
        yaml.safe_dump(workflow(), handle, sort_keys=False)

    proc = cli("validate", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0, proc.stdout
    proc = cli("run", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0, proc.stdout
    proc = cli("run", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0 and proc.stdout.count("[skipped]") == 5, proc.stdout

    manual = tmp_path / "manual"
    manual.mkdir()
    wf = workflow()
    failures = []
    for step in wf["steps"]:
        with open(str(work / ".hicexplorer-workflow" / "runs" / step["id"] / "run.json")) as handle:
            record = json.load(handle)
        assert record["exit_code"] == 0 and record["peak_rss_kb"] > 0
        for rel in step["outputs"].values():
            os.makedirs(os.path.dirname(os.path.join(str(manual), rel)), exist_ok=True)
        # the run's outputs are relative to its workdir; upstream outputs are
        # found at the same relative paths in the manual directory
        subprocess.run(record["argv"], cwd=str(manual), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        for name, rel in step["outputs"].items():
            fmt, files = COMPARISONS[(step["id"], name)]
            a, b = os.path.join(str(work), rel), os.path.join(str(manual), rel)
            if files is None:
                pairs = [(fmt, a, b)]
            elif fmt is None:
                pairs = [(f, a + suffix, b + suffix) for f, suffix in files]
            else:
                pairs = [(fmt, os.path.join(a, f), os.path.join(b, f)) for f in files]
            for f, x, y in pairs:
                assert os.path.exists(x), x
                code, out = equiv(f, x, y)
                if code != 0:
                    failures.append("{} {}: {}".format(step["id"], os.path.relpath(x, str(work)), out.strip()))
    assert not failures, "\n".join(failures)
