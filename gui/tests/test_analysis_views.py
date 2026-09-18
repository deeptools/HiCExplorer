"""Analysis views (PLAN 10.6): every view on the data outputs of a real run,
its link to the matrix browser, and figure export against the CLI.

Each test runs the C++ tool into a temporary directory, as the run history
records it, opens the view offscreen, checks what it shows and, where a
genomic position is involved, that ``MainWindow.navigate_browser`` moves the
matrix browser to it. Export equivalence: the exported figure is compared with
the figure the command-line tool writes for the same parameters, byte for
byte; the chicPlotViewpoint archive through ``equiv.py`` at E0 (tar members in
order, each image byte-identical), and hicQC.html after the named
normalisation of its random pandas table ids.

Needs HICX_CPP_BIN, HICX_REFERENCE_PYTHON (the drawing interpreter) and the
hicx_matrix module (the browser).
"""

import filecmp
import glob
import json
import os
import re
import subprocess

import pytest
from PySide6 import QtCore

from gui_support import CPP_BIN, DATA, EQUIV, REFERENCE_PYTHON, needs_reference, needs_tools

pytestmark = [needs_tools, needs_reference]


def run_tool(tmp_path, tool, args, name="cli"):
    cwd = tmp_path / name
    cwd.mkdir(exist_ok=True)
    argv = [os.path.join(CPP_BIN, tool)] + [a.replace("{data}", DATA) for a in args]
    proc = subprocess.run(argv, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)
    assert proc.returncode == 0, proc.stderr[-2000:]
    from hicexplorer_gui.analysis_views import StepRecord
    return StepRecord(tool, argv, str(cwd), name="{} test run".format(tool))


@pytest.fixture
def loader():
    from hicexplorer_gui.workflow.spec import SpecLoader
    return SpecLoader(CPP_BIN)


@pytest.fixture
def window(qtbot, tmp_path):
    """A ProjectTab of a freshly created project (named "window" so the rest
    of this file, written against the pre-redesign single-project MainWindow,
    needs no further changes: a ProjectTab exposes the same tabs/browser/
    open_analysis_view(s) surface the old MainWindow did, now scoped to one
    project instead of the whole window)."""
    import sys
    module_dir = os.environ.get("HICX_PYTHON_MODULE_DIR")
    if module_dir and module_dir not in sys.path:
        sys.path.insert(0, module_dir)
    from hicexplorer_gui.main_window import MainWindow
    from hicexplorer_gui.settings import Settings
    qsettings = QtCore.QSettings(str(tmp_path / "settings.ini"), QtCore.QSettings.IniFormat)
    win = MainWindow(Settings(qsettings))
    qtbot.addWidget(win)
    win.set_tools_dir(CPP_BIN)
    win.create_project(str(tmp_path / "project"))
    tab = win.active_project_tab()
    if tab.browser is None:
        pytest.skip("the matrix browser needs hicx_matrix (HICX_PYTHON_MODULE_DIR)")
    return tab


def open_in(window, record):
    view = window.open_analysis_view(record)
    assert view is not None, window.window.statusBar().currentMessage()
    return view


def same_bytes(a, b):
    assert os.path.isfile(a) and os.path.isfile(b), (a, b)
    assert filecmp.cmp(a, b, shallow=False), "{} differs from the CLI figure {}".format(a, b)


def test_qc_report_view_and_export(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import QCReportView
    record = run_tool(tmp_path, "hicQC", ["--logfiles", "{data}/QC/QC.log", "--outputFolder", "qc"])
    view = QCReportView(record, loader)
    qtbot.addWidget(view)
    assert set(view.tables) == {"QC", "discarded", "distance", "read_orientation", "unmapable"}
    assert view.tabs.count() == 5
    export = tmp_path / "export"
    written = view.export_figure(str(export))
    assert written == [str(export)]
    cli_files = sorted(os.listdir(os.path.join(record.cwd, "qc")))
    assert sorted(os.listdir(str(export))) == cli_files
    uuid = re.compile(rb"T_[0-9a-f]{5}")
    for name in cli_files:
        a, b = os.path.join(record.cwd, "qc", name), str(export / name)
        if name.endswith(".html"):
            assert uuid.sub(b"T_xxxxx", open(a, "rb").read()) == uuid.sub(b"T_xxxxx", open(b, "rb").read())
        else:
            same_bytes(b, a)


def test_qc_report_of_a_matrix_build_exports_like_hicqc(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import QCReportView
    record = run_tool(tmp_path, "hicBuildMatrix", [
        "-s", "{data}/R1_1000.bam", "{data}/R2_1000.bam", "--outFileName", "matrix.h5", "-bs", "100000",
        "--QCfolder", "qc", "--restrictionSequence", "AAGCTT", "--danglingSequence", "AGCT",
        "-rs", "{data}/hicFindRestSite/hindIII.bed"])
    view = QCReportView(record, loader)
    qtbot.addWidget(view)
    assert "QC" in view.tables
    written = view.export_figure(str(tmp_path / "export"))
    for png in sorted(glob.glob(os.path.join(record.cwd, "qc", "*.png"))):
        same_bytes(os.path.join(written[0], os.path.basename(png)), png)


def test_distance_decay_view_overlays_and_exports(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import DistanceDecayView
    record = run_tool(tmp_path, "hicPlotDistVsCounts", [
        "--matrices", "{data}/small_test_matrix_50kb_res.h5", "{data}/small_test_matrix_50kb_res.cool",
        "--labels", "h5", "cool", "--plotFile", "decay.png", "--outFileData", "decay.txt"])
    view = DistanceDecayView(record, loader)
    qtbot.addWidget(view)
    assert set(view.curves) == {("h5", "all"), ("cool", "all")}
    other = run_tool(tmp_path, "hicPlotDistVsCounts", [
        "--matrices", "{data}/Li_et_al_2015.h5", "--plotFile", "li.png", "--outFileData", "li.txt"], name="other")
    view.add_table(os.path.join(other.cwd, "li.txt"))
    assert len(view.curves) == 3
    without_table = run_tool(tmp_path, "hicPlotDistVsCounts", [
        "--matrices", "{data}/small_test_matrix_50kb_res.h5", "{data}/small_test_matrix_50kb_res.cool",
        "--labels", "h5", "cool", "--plotFile", "decay.png"], name="plotdata")
    from_plot_data = DistanceDecayView(without_table, loader)
    qtbot.addWidget(from_plot_data)
    assert set(from_plot_data.curves) == set(view.curves) - {("Li_et_al_2015.h5", "all")}
    same_bytes(view.export_figure(str(tmp_path / "exported.png"))[0], os.path.join(record.cwd, "decay.png"))


def test_hic_viewpoint_view_moves_the_browser_and_exports(qtbot, tmp_path, window):
    record = run_tool(tmp_path, "hicPlotViewpoint", [
        "--matrix", "{data}/Li_et_al_2015.h5", "--region", "X:3000000-3500000", "-rp", "X:3200000",
        "--outFileName", "viewpoint.png", "--interactionOutFileName", "interactions"])
    view = open_in(window, record)
    assert list(view.curves) == ["Li_et_al_2015.h5"]
    region = view.position_region(3200000)
    assert region == "X:2950000-3450000"
    with qtbot.waitSignal(window.browser.fetch_finished, timeout=60000):
        view.navigate.emit(view.matrix, region)
    assert window.tabs.currentWidget() is window.browser
    assert window.browser.sources[0].path == os.path.join(DATA, "Li_et_al_2015.h5")
    assert window.browser.chrom == "X"
    assert window.browser.view_x == pytest.approx((2950000, 3450000))
    same_bytes(view.export_figure(str(tmp_path / "exported.png"))[0], os.path.join(record.cwd, "viewpoint.png"))


def test_chic_viewpoint_view_and_export(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import ViewpointView
    record = run_tool(tmp_path, "chicPlotViewpoint", [
        "-if", "{data}/cHi-C/chicViewpoint/two_matrices.hdf5", "--range", "200000", "200000",
        "-o", "plots.tar.gz", "-t", "1", "--combinationMode", "oneGene", "--combinationName", "Eya1"])
    view = ViewpointView(record, loader)
    qtbot.addWidget(view)
    assert view.selector.count() >= 2 and "Eya1" in view.selector.itemText(0)
    assert view.findChild(type(view.status), "unavailable_chic_position") is not None
    exported = view.export_figure(str(tmp_path / "exported.tar.gz"))[0]
    proc = subprocess.run([REFERENCE_PYTHON, EQUIV, "compare", "--format", "tar_images", "--class", "E0",
                           os.path.join(record.cwd, "plots.tar.gz"), exported],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    assert proc.returncode == 0, proc.stdout


def test_aggregate_contacts_view_pairs_move_the_browser_and_export(qtbot, tmp_path, window):
    record = run_tool(tmp_path, "hicAggregateContacts", [
        "--matrix", "{data}/Li_et_al_2015.h5", "--numberOfBins", "30", "--BED",
        "{data}/hicAggregateContacts/test_regions.bed", "--mode", "intra-chr", "--range", "50000:900000",
        "--operationType", "sum", "--outFileName", "aggregate.png", "--outFilePrefixMatrix", "m",
        "--outFileContactPairs", "p"])
    view = open_in(window, record)
    assert [os.path.basename(p) for p in view.matrix_files] == ["m_genome.tab"]
    assert view.values.shape == (31, 31) or view.values.shape[0] == view.values.shape[1]
    assert view.pairs
    chrom, region = view.pair_region(0)
    with qtbot.waitSignal(window.browser.fetch_finished, timeout=60000):
        view.pair_table.selectRow(0)
    assert window.browser.chrom == chrom
    start, end = (int(v) for v in region.split(":")[1].split("-"))
    assert window.browser.view_x == pytest.approx((start, end))
    same_bytes(view.export_figure(str(tmp_path / "exported.png"))[0], os.path.join(record.cwd, "aggregate.png"))


def test_saddle_view_and_export(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import SaddleView
    record = run_tool(tmp_path, "hicCompartmentalization", [
        "-m", "{data}/hicPCA/obsexp_norm.h5", "--pca", "{data}/hicCompartmentalization/pca1.bedgraph",
        "-o", "ratio.png", "--outliers", "0.0", "--quantile", "30", "--outputMatrix", "matrix.npz"])
    view = SaddleView(record, loader)
    qtbot.addWidget(view)
    assert view.saddles.shape == (1, 30, 30)
    assert view.ratios.shape == (1, 29)
    same_bytes(view.export_figure(str(tmp_path / "exported.png"))[0], os.path.join(record.cwd, "ratio.png"))


def test_correlation_view_shows_hicrep_unavailable_and_exports_both_figures(qtbot, tmp_path, loader):
    from hicexplorer_gui.analysis_views import CorrelationView, UNAVAILABLE
    # The inputs and options of the harness case hicCorrelate.test_correlate.
    matrix = "{data}/hicCorrectMatrix/small_test_matrix_ICEcorrected_chrUextra_chr3LHet.h5"
    record = run_tool(tmp_path, "hicCorrelate", [
        "--matrices", matrix, matrix, "--labels", "first", "second", "--method", "spearman", "--log1p",
        "--colorMap", "jet", "--outFileNameHeatmap", "heatmap.png", "--outFileNameScatter", "scatter.png"])
    view = CorrelationView(record, loader)
    qtbot.addWidget(view)
    assert view.labels == ["first", "second"]
    assert view.results.shape == (2, 2) and view.results[0, 0] == pytest.approx(1.0)
    label = view.findChild(type(view.status), "unavailable_hicrep")
    assert label is not None and "PLAN 9.10" in label.text() and not label.isEnabled()
    written = view.export_figure(str(tmp_path / "exported.png"))
    assert len(written) == 2
    same_bytes(written[0], os.path.join(record.cwd, "heatmap.png"))
    same_bytes(written[1], os.path.join(record.cwd, "scatter.png"))


def test_differential_view_volcano_selects_a_tad_in_the_browser(qtbot, tmp_path, window):
    record = run_tool(tmp_path, "hicDifferentialTAD", [
        "-tm", "{data}/hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.cool",
        "-cm", "{data}/hicDifferentialTAD/GSM2644947_Auxin2days-R1.100000_chr1.cool",
        "-td", "{data}/hicDifferentialTAD/untreated_R1_domains.bed", "-m", "intra-TAD", "-mr", "one",
        "-t", "1", "-o", "diff", "--sharedMask", "--correctForMultipleTesting", "fdr"])
    view = open_in(window, record)
    accepted = [l for l in open(os.path.join(record.cwd, "diff_accepted.diff_tad")) if not l.startswith("#")]
    rejected = [l for l in open(os.path.join(record.cwd, "diff_rejected.diff_tad")) if not l.startswith("#")]
    assert len(view.rows) == len(accepted) + len(rejected)
    assert view.adjusted and "adjusted p-value intra-TAD" in view.columns
    points = view.volcano_points()
    assert points and all(y >= 0 for _, _, y in points)
    for key in ("differential_loops", "differential_compartments"):
        assert "PLAN 9.7" in view.findChild(type(view.status), "unavailable_" + key).text()
    row = points[0][0]
    with qtbot.waitSignal(window.browser.fetch_finished, timeout=60000):
        view.select_tad(row)
    chrom, span = view.region(row).split(":")
    start, end = (int(v) for v in span.split("-"))
    assert window.browser.sources[0].path.endswith("GSM2644945_Untreated-R1.100000_chr1.cool")
    assert window.browser.chrom == chrom
    assert window.browser.view_x == pytest.approx((start, end))
    exported = view.export_figure(str(tmp_path / "volcano.png"))[0]
    assert os.path.getsize(exported) > 0


def test_views_open_from_a_run_history_entry(qtbot, tmp_path, window):
    record = run_tool(tmp_path, "hicPlotDistVsCounts", [
        "--matrices", "{data}/small_test_matrix_50kb_res.h5", "--plotFile", "decay.png", "--outFileData", "d.txt"])
    history = tmp_path / "history"
    history.mkdir()
    summary = {"name": "t", "exit_code": 0, "workdir": record.cwd, "steps": [
        {"step": "decay", "tool": record.tool, "exit_code": 0, "argv": record.argv},
        {"step": "info", "tool": "hicInfo", "exit_code": 0, "argv": ["hicInfo", "-m", "x"]}]}
    (history / "summary.json").write_text(json.dumps(summary))
    views = window.open_analysis_views(str(history))
    assert [type(v).__name__ for v in views] == ["DistanceDecayView"]
    assert window.tabs.currentWidget() is views[0]
