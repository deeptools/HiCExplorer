"""The GUI shell (PLAN 10.4, redesigned 2026-09-18 into a tab per project):
project tabs, each project's own tool browser, runs from forms and the
workflow editor, data-driven tool filtering, and Settings moved into the menu.

Runs go through the workflow engine with the real C++ tools, and their
outputs are E0 against invoking the same command directly.
"""

import os
import re
import subprocess

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6 import QtCore  # noqa: E402

from gui_support import DATA, REPO, equiv_compare, needs_reference, needs_tools  # noqa: E402
from hicexplorer_gui.catalog import PLANNED_TOOLS  # noqa: E402
from hicexplorer_gui.main_window import MainWindow  # noqa: E402
from hicexplorer_gui.project import Project  # noqa: E402
from hicexplorer_gui.settings import Settings  # noqa: E402

pytestmark = needs_tools


def make_window(qtbot, tmp_path, name="settings.ini"):
    qsettings = QtCore.QSettings(str(tmp_path / name), QtCore.QSettings.IniFormat)
    win = MainWindow(Settings(qsettings))
    qtbot.addWidget(win)
    return win


@pytest.fixture
def window(qtbot, tmp_path):
    win = make_window(qtbot, tmp_path)
    win.create_project(str(tmp_path / "project"))
    return win


@pytest.fixture
def tab(window):
    """The ProjectTab of the single project the ``window`` fixture opened."""
    return window.active_project_tab()


def wait_finished(qtbot, controller, timeout=300000):
    with qtbot.waitSignal(controller.finished, timeout=timeout) as blocker:
        pass
    return blocker.args[0]


def ported_tool_names():
    # cpp/tools/<tool>.cpp plus the executables CMake builds from another
    # tool's source under a second name (hicQC from hicPrepareQCreport.cpp).
    tools_dir = os.path.join(REPO, "cpp", "tools")
    ported = {os.path.splitext(n)[0] for n in os.listdir(tools_dir)
              if n.endswith(".cpp") and n.startswith(("hic", "chic"))}
    with open(os.path.join(tools_dir, "CMakeLists.txt")) as handle:
        ported.update(re.findall(r"add_executable\(((?:hic|chic)\w+)", handle.read()))
    return ported


def test_tool_browser_lists_available_and_unavailable(tab):
    tree = tab.tool_tree
    available, missing = tree.topLevelItem(0), tree.topLevelItem(1)
    ported = ported_tool_names()
    assert available.text(0) == "Available ({})".format(len(ported))
    assert missing.childCount() == len(PLANNED_TOOLS) - len(ported)
    names = [missing.child(i).text(0) for i in range(missing.childCount())]
    assert "hicTADClassifier" in names
    item = missing.child(names.index("hicTADClassifier"))
    assert "not ported to C++ yet, PLAN tier 8" in item.toolTip(0)
    assert tab.open_tool_form("hicTADClassifier") is None


def test_invalid_form_does_not_run(qtbot, tab):
    page = tab.open_tool_form("hicInfo")
    assert page.run() is False
    assert "required" in page.form.fields["matrices"].error.text()
    assert not tab.controller.running


def _direct(argv, cwd):
    return subprocess.run(argv, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


@needs_reference
def test_run_from_form_is_e0_against_the_direct_command(qtbot, window, tab, tmp_path):
    project = tab.project
    matrix = os.path.join(DATA, "small_test_matrix.h5")

    page = tab.open_tool_form("hicMergeMatrixBins")
    out = os.path.join(project.path, "merged.h5")
    page.form.set_values({"matrix": matrix, "outFileName": out, "numBins": 5})
    assert page.run()
    assert wait_finished(qtbot, tab.controller) == 0
    direct_dir = tmp_path / "direct"
    direct_dir.mkdir()
    direct_out = str(direct_dir / "merged.h5")
    page.form.set_values({"matrix": matrix, "outFileName": direct_out, "numBins": 5})
    proc = _direct(page.form.argv(window.loader.executable("hicMergeMatrixBins")), str(direct_dir))
    assert proc.returncode == 0, proc.stderr
    code, text = equiv_compare("h5", out, direct_out)
    assert code == 0, text

    # An output outside the project, where users keep their results.
    elsewhere = tmp_path / "results elsewhere"
    elsewhere.mkdir()
    page = tab.open_tool_form("hicInfo")
    info = str(elsewhere / "info.txt")
    page.form.set_values({"matrices": [matrix], "outFileName": info})
    assert page.run()
    assert wait_finished(qtbot, tab.controller) == 0
    assert os.path.isfile(info)
    direct_info = str(direct_dir / "info.txt")
    page.form.set_values({"matrices": [matrix], "outFileName": direct_info})
    assert _direct(page.form.argv(window.loader.executable("hicInfo")), str(direct_dir)).returncode == 0
    code, text = equiv_compare("text", info, direct_info)
    assert code == 0, text

    history = project.history()
    assert [h["name"] for h in history[:2]] == ["hicInfo", "hicMergeMatrixBins"]
    step = history[0]["steps"][0]
    assert step["status"] == "succeeded" and step["exit_code"] == 0
    assert step["peak_rss_kb"] > 0 and step["cpu_seconds"] is not None
    # provenance: absolute output paths and the command line in the history
    assert step["outputs"] == {"outFileName": info}
    assert info in step["argv"] and os.path.isabs(step["argv"][0])
    assert history[1]["steps"][0]["outputs"] == {"outFileName": out}
    tab.run_view.refresh_history()
    assert tab.run_view.history.count() == 2
    tab.run_view.history.setCurrentRow(0)
    assert tab.run_view.table.item(0, 2).text() == "succeeded"


def test_output_in_a_missing_directory_is_refused(qtbot, tab, tmp_path):
    page = tab.open_tool_form("hicInfo")
    missing = str(tmp_path / "no such directory" / "info.txt")
    page.form.set_values({"matrices": [os.path.join(DATA, "small_test_matrix.h5")], "outFileName": missing})
    # validation fails at once, so wait for the signal around run()
    with qtbot.waitSignal(tab.controller.finished, timeout=60000) as blocker:
        assert page.run()
    assert blocker.args[0] == 1
    assert "does not exist" in tab.run_view.log.toPlainText()
    assert not os.path.exists(os.path.dirname(missing))


def test_cancel_stops_a_run(qtbot, tab):
    page = tab.open_tool_form("hicCorrectMatrix", "correct")
    matrix = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    page.form.set_values({"matrix": matrix, "outFileName": os.path.join(tab.project.path, "c.cool"),
                          "correctionMethod": "ICE"}, "correct")
    assert page.run()
    qtbot.waitUntil(lambda: tab.controller.step_records()
                    and tab.controller.step_records()[0].get("status") == "running", timeout=60000)
    tab.run_view.cancel()
    assert wait_finished(qtbot, tab.controller) == 130
    assert tab.project.history()[0]["exit_code"] == 130


@needs_reference
def test_workflow_editor_connects_saves_validates_and_runs(qtbot, tab, tmp_path):
    editor = tab.editor
    editor.name_edit.setText("merge-info")
    first = editor.add_step("hicMergeMatrixBins", "merge")
    first.form.set_values({"matrix": os.path.join(DATA, "small_test_matrix.h5"), "outFileName": "merged.h5",
                           "numBins": 5})
    second = editor.add_step("hicInfo", "info")
    second.form.refresh_links()
    link = second.form.fields["matrices"].link
    labels = [link.itemText(i) for i in range(link.count())]
    assert "merge: outFileName" in labels
    link.activated.emit(labels.index("merge: outFileName"))
    second.form.fields["outFileName"].edit.setText("info.txt")
    assert second.form.values()["matrices"] == ["${steps.merge.outputs.outFileName}"]

    path = editor.save()
    assert os.path.isfile(path)
    assert editor.validate() == [] or all(not m.startswith("error") for m in editor.validate())
    with qtbot.waitSignal(tab.controller.finished, timeout=300000) as blocker:
        editor.run()
    assert blocker.args[0] == 0
    workdir = tab.project.workdir_for("merge-info")
    assert os.path.isfile(os.path.join(workdir, "merged.h5"))
    assert os.path.isfile(os.path.join(workdir, "info.txt"))

    # reload the saved workflow into a fresh editor state
    editor.load(path)
    assert [e.step_id for e in editor.step_editors] == ["merge", "info"]
    assert editor.step_editors[1].form.values()["matrices"] == ["${steps.merge.outputs.outFileName}"]
    assert Project(tab.project.path).workflows() == [path]


# -- the redesign: a project per tab -----------------------------------------

def test_two_projects_open_at_once_keep_independent_state(qtbot, tmp_path):
    win = make_window(qtbot, tmp_path)
    win.create_project(str(tmp_path / "alpha"))
    win.create_project(str(tmp_path / "beta"))
    assert win.tabs.count() == 2
    assert [win.tabs.tabText(i) for i in range(2)] == ["alpha", "beta"]

    alpha, beta = win.tabs.widget(0), win.tabs.widget(1)
    assert alpha.project.path != beta.project.path
    assert alpha is not beta
    assert alpha.controller is not beta.controller
    assert alpha.tabs is not beta.tabs

    # opening a tool form in one project does not touch the other
    alpha.open_tool_form("hicInfo")
    assert any(alpha.tabs.tabText(i) == "hicInfo" for i in range(alpha.tabs.count()))
    assert not any(beta.tabs.tabText(i) == "hicInfo" for i in range(beta.tabs.count()))

    # loading data into one project does not filter the other's tool list
    alpha.load_data_file(os.path.join(DATA, "small_test_matrix.cool"))
    assert alpha.current_data is not None
    assert beta.current_data is None
    assert alpha.filtered_entries != beta.entries or alpha.current_data.format != beta.current_data


def test_closing_one_project_tab_leaves_the_other_untouched(qtbot, tmp_path):
    win = make_window(qtbot, tmp_path)
    win.create_project(str(tmp_path / "alpha"))
    win.create_project(str(tmp_path / "beta"))
    beta = win.tabs.widget(1)
    win._close_project_tab(0)
    assert win.tabs.count() == 1
    assert win.tabs.widget(0) is beta
    assert beta.project.name == "beta"


def test_closing_a_project_with_a_run_in_progress_asks_first(qtbot, tmp_path, monkeypatch):
    win = make_window(qtbot, tmp_path)
    win.create_project(str(tmp_path / "project"))
    tab = win.active_project_tab()
    project = tab.project  # kept: tab itself is deleteLater()'d by the confirmed close below
    page = tab.open_tool_form("hicCorrectMatrix", "correct")
    matrix = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    page.form.set_values({"matrix": matrix, "outFileName": os.path.join(tab.project.path, "c.cool"),
                          "correctionMethod": "ICE"}, "correct")
    assert page.run()
    qtbot.waitUntil(lambda: tab.controller.running, timeout=60000)

    from PySide6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *a, **k: QtWidgets.QMessageBox.No)
    win._close_project_tab(0)
    assert win.tabs.count() == 1  # declined: still open
    assert tab.controller.running

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", lambda *a, **k: QtWidgets.QMessageBox.Yes)
    win._close_project_tab(0)
    assert win.tabs.count() == 0  # confirmed: closed, and the run is told to cancel
    # The tab (and its controller, a QObject parented to it) is deleteLater()'d
    # by the close above, so what is left to check is on disk, not the Qt
    # object: the history entry the cancelled run's worker thread still writes.
    qtbot.waitUntil(lambda: bool(project.history()) and project.history()[0].get("exit_code") is not None,
                    timeout=60000)
    assert project.history()[0]["exit_code"] == 130


# -- data-driven tool filtering ------------------------------------------------

def test_loading_a_cool_file_filters_to_matrix_compatible_tools(tab):
    from hicexplorer_gui.catalog import entries_for_format
    matrix = os.path.join(DATA, "small_test_matrix.cool")
    entry = tab.load_data_file(matrix)
    assert entry.format == "cool"
    assert tab.current_data is entry
    expected = {e.name for e in entries_for_format(tab.entries, "cool")}
    assert expected  # the real spec data of the built tools gives a non-empty set
    assert {e.name for e in tab.filtered_entries} == expected
    assert "hicInfo" in expected  # hicInfo --help-json declares matrices: input, [h5, cool, mcool]
    assert "hicBuildMatrix" not in expected  # takes BAM/SAM, not cool
    shown = {tab.tool_tree.topLevelItem(0).child(i).text(0)
            for i in range(tab.tool_tree.topLevelItem(0).childCount())}
    assert shown == expected


def test_loading_a_bam_file_filters_to_bam_compatible_tools(tab):
    from hicexplorer_gui.catalog import entries_for_format
    bam = os.path.join(DATA, "small_test_R1_unsorted.bam")
    entry = tab.load_data_file(bam)
    assert entry.format == "bam"
    expected = {e.name for e in entries_for_format(tab.entries, "bam")}
    assert expected
    assert "hicBuildMatrix" in expected  # samFiles: input, [bam, sam]
    assert "hicInfo" not in expected  # hicInfo takes matrices, not alignments
    assert {e.name for e in tab.filtered_entries} == expected


def test_loading_a_fastq_file_filters_to_hicalignreads(tab, tmp_path):
    # PLAN tier 13: hicAlignReads wraps minibwa and declares fastq/fastq.gz
    # as an input format in its own --help-json, so it is found the same
    # way test_loading_a_cool_file_filters_to_matrix_compatible_tools finds
    # hicInfo for cool -- no FASTQ-specific code path left to test here.
    from hicexplorer_gui.catalog import entries_for_format
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")
    entry = tab.load_data_file(str(fastq))
    assert entry.format == "fastq"
    expected = {e.name for e in entries_for_format(tab.entries, "fastq")}
    assert expected  # the real spec data of the built tools gives a non-empty set
    assert "hicAlignReads" in expected  # hicAlignReads --help-json declares inFile: input, [fastq, fastq.gz]
    assert {e.name for e in tab.filtered_entries} == expected
    assert "tool(s) accept fastq input" in tab.data_status.text()
    shown = {tab.tool_tree.topLevelItem(0).child(i).text(0)
            for i in range(tab.tool_tree.topLevelItem(0).childCount())}
    assert shown == expected


def test_loading_a_fastq_file_says_no_tool_reads_it_directly_when_none_is_built(tab, tmp_path, monkeypatch):
    # The FASTQ_MESSAGE fallback still exists for a tool directory that
    # predates PLAN tier 13 (no hicAlignReads executable in it): simulate
    # that by filtering hicAlignReads out of the catalog entries the tab
    # already has, rather than needing a second, hicAlignReads-less build.
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")
    monkeypatch.setattr(tab, "entries", [e for e in tab.entries if e.name != "hicAlignReads"])
    entry = tab.load_data_file(str(fastq))
    assert entry.format == "fastq"
    assert tab.filtered_entries == []
    assert "no available tool reads FASTQ directly" in tab.data_status.text()
    assert "align it to BAM first" in tab.data_status.text()
    assert tab.tool_tree.topLevelItem(0).childCount() == 0


def test_loading_a_matrix_file_auto_opens_it_in_the_matrix_browser(tab):
    matrix = os.path.join(DATA, "small_test_matrix.cool")
    entry = tab.load_data_file(matrix)
    assert entry.format == "cool"
    assert tab.browser is not None
    assert tab.browser.sources[0] is not None
    assert os.path.abspath(tab.browser.sources[0].path) == os.path.abspath(matrix)
    # auto-open also switches to the Matrix browser tab, the same as clicking
    # the manual "Open in Matrix browser" button would
    assert tab.tabs.currentWidget() is tab.browser


def test_loading_a_non_matrix_file_does_not_touch_the_matrix_browser(tab):
    bam = os.path.join(DATA, "small_test_R1_unsorted.bam")
    entry = tab.load_data_file(bam)
    assert entry.format == "bam"
    assert tab.browser is not None
    assert tab.browser.sources[0] is None
    assert tab.tabs.currentWidget() is not tab.browser


def test_manual_open_in_browser_button_still_opens_the_selected_data_file(tab):
    # the "Open in Matrix browser" button (open_in_browser_button) stays for
    # files opened another way, e.g. reselecting an earlier load or a tool's
    # own output; auto-open on load must not replace it.
    a = os.path.join(DATA, "small_test_matrix.cool")
    b = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    tab.load_data_file(a)
    tab.load_data_file(b)
    assert os.path.abspath(tab.browser.sources[0].path) == os.path.abspath(b)
    tab.set_current_data(tab.data_files[0])
    assert tab.open_in_browser_button.isEnabled()
    tab.open_in_browser_button.click()
    assert os.path.abspath(tab.browser.sources[0].path) == os.path.abspath(a)


def test_without_data_the_full_tool_list_is_shown(tab):
    assert tab.current_data is None
    assert tab.filtered_entries == tab.entries
    assert "showing every available tool" in tab.data_status.text()


def test_loading_a_second_file_updates_the_filter(tab):
    cool = os.path.join(DATA, "small_test_matrix.cool")
    bam = os.path.join(DATA, "small_test_R1_unsorted.bam")
    tab.load_data_file(cool)
    cool_names = {e.name for e in tab.filtered_entries}
    tab.load_data_file(bam)
    bam_names = {e.name for e in tab.filtered_entries}
    assert cool_names != bam_names
    assert tab.current_data.path == bam
    assert len(tab.data_files) == 2
    # both loaded files stay visible, the current one is marked
    texts = [tab.data_list.item(i).text() for i in range(tab.data_list.count())]
    assert any("small_test_matrix.cool" in t and t.startswith("  ") for t in texts)
    assert any("small_test_R1_unsorted.bam" in t and t.startswith(">") for t in texts)


# -- settings moved into the menu ---------------------------------------------

def test_settings_is_not_a_project_tab(window, tab):
    titles_top = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert "Settings" not in titles_top
    titles_nested = [tab.tabs.tabText(i) for i in range(tab.tabs.count())]
    assert "Settings" not in titles_nested


def test_settings_opens_from_the_file_menu_and_applies(window, tmp_path):
    dialog = window.open_settings_dialog()
    assert dialog.isVisible()
    dialog.page.edit.setText(str(tmp_path / "missing"))
    assert dialog.page.apply() is False
    assert "is not a directory" in dialog.page.status.text()
    dialog.close()
