"""Workflow templates (PLAN 10.7): the picker, and every template end to end.

Each template is instantiated with small committed real data and run with
``hicexplorer-workflow run``. Every step's logged argv is then run by hand in
a separate directory, in step order, and every declared output (a file, a
directory's files, or the files of an output prefix) is compared: byte for
byte, else at class E0 with ``cpp/scripts/equiv.py``. Two named
normalisations are applied before that and nowhere else: hicQC.html names its
tables after a random pandas Styler uuid (``T_`` and five hexadecimal
characters), and hicInfo writes the time it ran on its ``Date:`` line.

Needs HICX_CPP_BIN and HICX_REFERENCE_PYTHON (gui_support also takes the
latter as HICX_PLOT_PYTHON for the drawing tools).
"""

import filecmp
import glob
import json
import os
import re
import subprocess
import sys

import pytest
import yaml

from gui_support import CPP_BIN, DATA, EQUIV, GUI_DIR, REFERENCE_PYTHON, needs_reference, needs_tools
from hicexplorer_gui.templates import TEMPLATE_DIR, Template, load_templates

pytestmark = [needs_tools, needs_reference]

VALUES = {
    "hic": {"r1": os.path.join(DATA, "small_test_R1_unsorted.bam"),
            "r2": os.path.join(DATA, "small_test_R2_unsorted.bam"),
            "restrictionCutFile": os.path.join(DATA, "DpnII.bed"),
            "geneTrack": os.path.join(DATA, "dm3_genes.bed.gz")},
    "differential": {
        "targetReplicates": [os.path.join(DATA, "hicDifferentialTAD", "GSM2644945_Untreated-R1.100000_chr1.cool")],
        "controlReplicates": [os.path.join(DATA, "hicDifferentialTAD", "GSM2644947_Auxin2days-R1.100000_chr1.cool")]},
    "capture_hic": {"matrices": [os.path.join(DATA, "cHi-C", "FL-E13-5_chr1.cool"),
                                 os.path.join(DATA, "cHi-C", "MB-E10-5_chr1.cool")],
                    "referencePoints": os.path.join(DATA, "cHi-C", "referencePoints.bed")},
    "conversion_qc": {"matrixA": os.path.join(DATA, "hicConvertFormat", "GM12878_combined_30.chr21_chr22.v6.hic"),
                      "matrixB": os.path.join(DATA, "hicConvertFormat", "GM12878_combined_30.chr21_chr22.v7.hic")},
}

FORMATS = [("_export.tar.gz", "tar_members"), (".tar.gz", "tar_images"), (".hdf5", "chic_hdf5"), (".h5", "h5"), (".cool", "cool"),
           (".bedgraph", "text"), (".npz", "npz"), (".png", "png"), (".bam", "plain")]

NORMALISATIONS = {
    "pandas_styler_uuid": (re.compile(rb"T_[0-9a-f]{5}"), b"T_xxxxx"),
    "hicinfo_date": (re.compile(rb"^Date:\t[^\n]*$", re.M), b"Date:\t<time of the run>"),
}


def _format(path):
    for suffix, fmt in FORMATS:
        if path.endswith(suffix):
            return fmt
    return "text"


def _normalised(path):
    with open(path, "rb") as handle:
        data = handle.read()
    applied = []
    for name, (pattern, replacement) in NORMALISATIONS.items():
        if (name == "pandas_styler_uuid" and path.endswith(".html")) or \
                (name == "hicinfo_date" and data.startswith(b"# Matrix information file")):
            data, count = pattern.subn(replacement, data)
            if count:
                applied.append(name)
    return data, applied


def compare(a, b):
    """None when equal; otherwise what differs."""
    if filecmp.cmp(a, b, shallow=False):
        return None
    data_a, applied = _normalised(a)
    data_b, _ = _normalised(b)
    if applied:
        return None if data_a == data_b else "differs after {}".format(", ".join(applied))
    proc = subprocess.run([REFERENCE_PYTHON, EQUIV, "compare", "--format", _format(a), "--class", "E0", a, b],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    return None if proc.returncode == 0 else proc.stdout.strip()[-800:]


def output_files(root, rel):
    """The files of a declared output: the file, a directory's files, or the
    files starting with an output prefix."""
    path = os.path.join(root, rel)
    if os.path.isfile(path):
        return [rel]
    if os.path.isdir(path):
        return sorted(os.path.relpath(f, root) for f in glob.glob(os.path.join(path, "**"), recursive=True)
                      if os.path.isfile(f))
    return sorted(os.path.relpath(f, root) for f in glob.glob(path + "*") if os.path.isfile(f))


def cli(*args):
    env = dict(os.environ, PYTHONPATH=GUI_DIR)
    return subprocess.run([sys.executable, "-m", "hicexplorer_gui.workflow.cli"] + list(args),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env)


@pytest.mark.parametrize("name", sorted(VALUES))
def test_template_runs_end_to_end_and_every_output_matches_the_logged_commands(tmp_path, name):
    template = Template(os.path.join(TEMPLATE_DIR, name + ".yaml"))
    work = tmp_path / "workflow"
    work.mkdir()
    path = str(work / (name + ".yaml"))
    workflow = template.instantiate(VALUES[name])
    with open(path, "w") as handle:
        yaml.safe_dump(workflow, handle, sort_keys=False)
    proc = cli("validate", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0, proc.stdout
    proc = cli("run", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0, proc.stdout

    manual = tmp_path / "manual"
    manual.mkdir()
    failures, compared = [], 0
    for step in workflow["steps"]:
        with open(str(work / ".hicexplorer-workflow" / "runs" / step["id"] / "run.json")) as handle:
            record = json.load(handle)
        assert record["exit_code"] == 0, (step["id"], record)
        for rel in step["outputs"].values():
            os.makedirs(os.path.dirname(os.path.join(str(manual), rel)) or str(manual), exist_ok=True)
        subprocess.run(record["argv"], cwd=str(manual), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for rel in step["outputs"].values():
            files = output_files(str(work), rel)
            assert files, "{}: no files for output {}".format(step["id"], rel)
            assert files == output_files(str(manual), rel), (step["id"], rel)
            for file_rel in files:
                problem = compare(os.path.join(str(work), file_rel), os.path.join(str(manual), file_rel))
                compared += 1
                if problem:
                    failures.append("{} {}: {}".format(step["id"], file_rel, problem))
    print("{}: {} steps, {} output files compared".format(name, len(workflow["steps"]), compared))
    assert not failures, "\n".join(failures)


def test_every_template_loads_with_documented_parameters():
    templates = load_templates()
    assert sorted(t.name for t in templates) == sorted(VALUES)
    for template in templates:
        assert template.description
        for parameter in template.parameters:
            assert parameter.get("description"), (template.name, parameter["name"])
        required = {p["name"] for p in template.parameters if "default" not in p}
        assert required == set(VALUES[template.name]), template.name


def test_capture_template_runs_the_whole_chic_suite():
    from hicexplorer_gui.catalog import tool_entries
    entries, _ = tool_entries(CPP_BIN)
    template = Template(os.path.join(TEMPLATE_DIR, "capture_hic.yaml"))
    assert list(template.unavailable_steps(entries)) == []
    assert set(template.tools()) == {"chicQualityControl", "chicViewpointBackgroundModel", "chicViewpoint",
                                     "chicPlotViewpoint", "chicSignificantInteractions", "chicAggregateStatistic",
                                     "chicDifferentialTest", "chicExportData"}
    available = {entry.name for entry in entries if entry.available}
    assert set(template.tools()) <= available


def test_picker_creates_the_workflow_in_the_project(qtbot, tmp_path):
    from PySide6 import QtCore
    from hicexplorer_gui.main_window import MainWindow
    from hicexplorer_gui.settings import Settings
    from hicexplorer_gui.template_picker import TemplatePicker

    qsettings = QtCore.QSettings(str(tmp_path / "settings.ini"), QtCore.QSettings.IniFormat)
    window = MainWindow(Settings(qsettings))
    qtbot.addWidget(window)
    window.set_tools_dir(CPP_BIN)
    window.create_project(str(tmp_path / "project"))
    picker = TemplatePicker(window.entries, window)
    qtbot.addWidget(picker)
    titles = [picker.list.item(i).text() for i in range(picker.list.count())]
    picker.list.setCurrentRow(titles.index("Capture Hi-C"))
    assert picker.unavailable.text() == "All steps are available."
    assert set(picker.fields) == {"matrices", "referencePoints", "sparsity", "range", "fixateRange", "plotGene",
                                  "xFoldBackground", "pValue", "alpha", "statisticTest", "correction"}
    assert picker.fields["sparsity"][1].text() == "0.05"
    picker.create()
    assert picker.workflow is None and "needs a value" in picker.status.text()
    values = VALUES["capture_hic"]
    picker.fields["matrices"][1].setText(yaml.safe_dump(values["matrices"], default_flow_style=True).strip())
    picker.fields["referencePoints"][1].setText(values["referencePoints"])
    picker.name.setText("my capture run")
    picker.create()
    assert picker.workflow is not None
    assert picker.workflow["steps"][0]["args"]["matrices"] == values["matrices"]
    assert picker.workflow["steps"][0]["args"]["sparsity"] == 0.05
    path = window.create_workflow(picker.workflow)
    assert os.path.isfile(path) and path.startswith(window.project.workflows_dir)
    assert window.tabs.currentWidget() is window.editor
    proc = cli("validate", path, "--tools-dir", CPP_BIN)
    assert proc.returncode == 0, proc.stdout
