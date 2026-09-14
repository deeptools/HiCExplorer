"""The GUI shell (PLAN 10.4): projects, tool browser, runs from forms and the workflow editor.

Runs go through the workflow engine with the real C++ tools, and their
outputs are E0 against invoking the same command directly.
"""

import os
import subprocess

import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

from PySide6 import QtCore  # noqa: E402

from gui_support import DATA, equiv_compare, needs_reference, needs_tools  # noqa: E402
from hicexplorer_gui.catalog import PLANNED_TOOLS  # noqa: E402
from hicexplorer_gui.main_window import MainWindow  # noqa: E402
from hicexplorer_gui.project import Project  # noqa: E402
from hicexplorer_gui.settings import Settings  # noqa: E402

pytestmark = needs_tools


@pytest.fixture
def window(qtbot, tmp_path):
    qsettings = QtCore.QSettings(str(tmp_path / "settings.ini"), QtCore.QSettings.IniFormat)
    win = MainWindow(Settings(qsettings))
    qtbot.addWidget(win)
    win.create_project(str(tmp_path / "project"))
    return win


def wait_finished(qtbot, controller, timeout=300000):
    with qtbot.waitSignal(controller.finished, timeout=timeout) as blocker:
        pass
    return blocker.args[0]


def test_tool_browser_lists_available_and_unavailable(window):
    tree = window.tool_tree
    available, missing = tree.topLevelItem(0), tree.topLevelItem(1)
    assert available.text(0) == "Available (30)"
    assert missing.childCount() == len(PLANNED_TOOLS) - 30
    names = [missing.child(i).text(0) for i in range(missing.childCount())]
    assert "hicPlotMatrix" in names
    item = missing.child(names.index("hicPlotMatrix"))
    assert "not ported to C++ yet, PLAN tier 7" in item.toolTip(0)
    assert window.open_tool_form("hicPlotMatrix") is None


def test_settings_reject_a_missing_directory_inline(window, tmp_path):
    window.settings_page.edit.setText(str(tmp_path / "missing"))
    assert window.settings_page.apply() is False
    assert "is not a directory" in window.settings_page.status.text()


def test_invalid_form_does_not_run(qtbot, window):
    page = window.open_tool_form("hicInfo")
    assert page.run() is False
    assert "required" in page.form.fields["matrices"].error.text()
    assert not window.controller.running


def _direct(argv, cwd):
    return subprocess.run(argv, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


@needs_reference
def test_run_from_form_is_e0_against_the_direct_command(qtbot, window, tmp_path):
    project = window.project
    matrix = os.path.join(DATA, "small_test_matrix.h5")

    page = window.open_tool_form("hicMergeMatrixBins")
    out = os.path.join(project.path, "merged.h5")
    page.form.set_values({"matrix": matrix, "outFileName": out, "numBins": 5})
    assert page.run()
    assert wait_finished(qtbot, window.controller) == 0
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
    page = window.open_tool_form("hicInfo")
    info = str(elsewhere / "info.txt")
    page.form.set_values({"matrices": [matrix], "outFileName": info})
    assert page.run()
    assert wait_finished(qtbot, window.controller) == 0
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
    window.run_view.refresh_history()
    assert window.run_view.history.count() == 2
    window.run_view.history.setCurrentRow(0)
    assert window.run_view.table.item(0, 2).text() == "succeeded"


def test_output_in_a_missing_directory_is_refused(qtbot, window, tmp_path):
    page = window.open_tool_form("hicInfo")
    missing = str(tmp_path / "no such directory" / "info.txt")
    page.form.set_values({"matrices": [os.path.join(DATA, "small_test_matrix.h5")], "outFileName": missing})
    # validation fails at once, so wait for the signal around run()
    with qtbot.waitSignal(window.controller.finished, timeout=60000) as blocker:
        assert page.run()
    assert blocker.args[0] == 1
    assert "does not exist" in window.run_view.log.toPlainText()
    assert not os.path.exists(os.path.dirname(missing))


def test_cancel_stops_a_run(qtbot, window):
    page = window.open_tool_form("hicCorrectMatrix", "correct")
    matrix = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    page.form.set_values({"matrix": matrix, "outFileName": os.path.join(window.project.path, "c.cool"),
                          "correctionMethod": "ICE"}, "correct")
    assert page.run()
    qtbot.waitUntil(lambda: window.controller.step_records()
                    and window.controller.step_records()[0].get("status") == "running", timeout=60000)
    window.run_view.cancel()
    assert wait_finished(qtbot, window.controller) == 130
    assert window.project.history()[0]["exit_code"] == 130


@needs_reference
def test_workflow_editor_connects_saves_validates_and_runs(qtbot, window, tmp_path):
    editor = window.editor
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
    with qtbot.waitSignal(window.controller.finished, timeout=300000) as blocker:
        editor.run()
    assert blocker.args[0] == 0
    workdir = window.project.workdir_for("merge-info")
    assert os.path.isfile(os.path.join(workdir, "merged.h5"))
    assert os.path.isfile(os.path.join(workdir, "info.txt"))

    # reload the saved workflow into a fresh editor state
    editor.load(path)
    assert [e.step_id for e in editor.step_editors] == ["merge", "info"]
    assert editor.step_editors[1].form.values()["matrices"] == ["${steps.merge.outputs.outFileName}"]
    assert Project(window.project.path).workflows() == [path]
