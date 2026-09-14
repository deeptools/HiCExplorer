"""Running tools and workflows through the workflow engine, and the run view.

A RunController runs one engine Runner on a worker thread and reports
through Qt signals; a finished run is copied into the project's history. The
RunView shows the live log, one row per step (status, exit code, peak RSS,
CPU time), a cancel button and the project's history.
"""

import datetime
import json
import os
import shutil
import threading

import yaml
from PySide6 import QtCore, QtWidgets

from .forms import is_list_arg
from .workflow.model import WorkflowError, load_workflow
from .workflow.runner import Runner

STATUS_COLUMNS = ["Step", "Tool", "Status", "Exit", "Peak RSS (MB)", "CPU (s)", "Wall (s)"]


def single_tool_workflow(spec, subcommand, values, project_dir, threads=1):
    """A one-step workflow document for a form's values.

    Output files inside the project become declared outputs, so the engine
    tracks and hashes them; the project directory is the work directory.
    """
    args = {}
    outputs = {}
    by_dest = spec.arguments_by_dest(subcommand)
    thread_dest = spec.thread_dest(subcommand)
    for dest, value in values.items():
        arg = by_dest[dest]
        if dest == thread_dest:
            try:
                threads = max(1, int(value))
            except (TypeError, ValueError):
                args[dest] = value
            continue
        info = arg.get("file") or {}
        if info.get("role") == "output" and value is not None:
            items = value if isinstance(value, list) else [value]
            refs = []
            for index, path in enumerate(items):
                absolute = os.path.abspath(os.path.join(project_dir, path))
                inside = os.path.commonpath([absolute, project_dir]) == project_dir and absolute != project_dir
                if inside and isinstance(path, str) and not path.startswith("${"):
                    name = dest if len(items) == 1 else "{}_{}".format(dest, index + 1)
                    outputs[name] = os.path.relpath(absolute, project_dir)
                    refs.append("${outputs.%s}" % name)
                else:
                    refs.append(path)
            args[dest] = refs if is_list_arg(arg) else refs[0]
        else:
            args[dest] = value
    step = {"id": "run", "tool": spec.tool, "threads": threads, "args": args, "outputs": outputs}
    if subcommand is not None:
        step["subcommand"] = subcommand
    return {"version": 1, "name": spec.tool + (" " + subcommand if subcommand else ""), "threads": threads,
            "steps": [step]}


class RunController(QtCore.QObject):
    """Runs one workflow at a time on a worker thread."""

    log = QtCore.Signal(str)
    started = QtCore.Signal(str)
    finished = QtCore.Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.runner = None
        self.workflow = None
        self.history_dir = None
        self.thread = None
        self.exit_code = None
        self.name = ""

    @property
    def running(self):
        return self.thread is not None and self.thread.is_alive()

    def run_tool(self, project, loader, spec, subcommand, values, threads=1):
        document = single_tool_workflow(spec, subcommand, values, project.path, threads)
        history = project.new_history_entry(document["name"])
        path = os.path.join(history, "workflow.yaml")
        with open(path, "w") as handle:
            yaml.safe_dump(document, handle, sort_keys=False)
        return self._start(project, loader, path, project.path, document["name"], history, force=True)

    def run_workflow(self, project, loader, workflow_path, name, force=False):
        workdir = project.workdir_for(name)
        history = project.new_history_entry(name)
        shutil.copyfile(workflow_path, os.path.join(history, "workflow.yaml"))
        return self._start(project, loader, workflow_path, workdir, name, history, force=force)

    def _start(self, project, loader, path, workdir, name, history, force):
        if self.running:
            raise RuntimeError("a run is still in progress")
        self.name = name
        self.history_dir = history
        self.exit_code = None
        try:
            self.workflow = load_workflow(path, workdir)
        except WorkflowError as exc:
            self.workflow = None
            self._write_summary(1, str(exc), workdir, path)
            self.log.emit("error: {}".format(exc))
            self.finished.emit(1)
            return False
        self.runner = Runner(self.workflow, loader, force=force, log=self.log.emit)
        self._started_at = _now()
        self.thread = threading.Thread(target=self._run, name="hicexplorer-run", daemon=True)
        self.started.emit(name)
        self.thread.start()
        return True

    def _run(self):
        try:
            code = self.runner.run()
        except Exception as exc:  # noqa: BLE001 - reported in the log, never swallowed
            self.log.emit("engine error: {!r}".format(exc))
            code = 1
        self.exit_code = code
        self._write_summary(code, None, self.workflow.workdir, self.workflow.path)
        self.finished.emit(code)

    def cancel(self):
        if self.runner is not None:
            self.runner.cancel()

    def wait(self, timeout=None):
        if self.thread is not None:
            self.thread.join(timeout)

    def step_records(self):
        """run.json of every step of the current workflow that has one."""
        if self.workflow is None:
            return []
        out = []
        for step in self.workflow.steps:
            record = _read_json(os.path.join(self.workflow.state_dir, "runs", step.id, "run.json"))
            if self.runner is not None and step.id in self.runner.status:
                record = dict(record or {}, step=step.id, tool=step.tool)
                record["status"] = self.runner.status[step.id]
            elif record is None or record.get("start", "") < getattr(self, "_started_at", ""):
                record = {"step": step.id, "tool": step.tool, "status": "running" if self.running else "not run"}
            out.append(record)
        return out

    def _write_summary(self, code, error, workdir, path):
        steps = []
        if self.workflow is not None:
            for step in self.workflow.steps:
                run_dir = os.path.join(self.workflow.state_dir, "runs", step.id)
                target = os.path.join(self.history_dir, "steps", step.id)
                if os.path.isdir(run_dir):
                    shutil.copytree(run_dir, target, dirs_exist_ok=True)
                record = _read_json(os.path.join(target, "run.json")) or {}
                status = self.runner.status.get(step.id, "not run") if self.runner else "not run"
                steps.append({"step": step.id, "tool": step.tool, "status": status,
                              "exit_code": record.get("exit_code"), "peak_rss_kb": record.get("peak_rss_kb"),
                              "cpu_seconds": record.get("cpu_seconds"), "wall_seconds": record.get("wall_seconds")})
        summary = {"name": self.name, "exit_code": code, "error": error, "workflow": path, "workdir": workdir,
                   "started": getattr(self, "_started_at", _now()), "finished": _now(), "steps": steps}
        with open(os.path.join(self.history_dir, "summary.json"), "w") as handle:
            json.dump(summary, handle, indent=2)


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _fmt(value, scale=1.0, digits=2):
    if value is None:
        return ""
    return "{:.{}f}".format(value * scale, digits)


class RunView(QtWidgets.QWidget):
    """Live log, step table, cancel, and the project's run history."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.controller = None
        self.project = None
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.addWidget(splitter)

        history_box = QtWidgets.QWidget()
        history_layout = QtWidgets.QVBoxLayout(history_box)
        history_layout.setContentsMargins(0, 0, 0, 0)
        history_layout.addWidget(QtWidgets.QLabel("History"))
        self.history = QtWidgets.QListWidget()
        self.history.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.history.setWordWrap(True)
        self.history.currentItemChanged.connect(self._history_selected)
        history_layout.addWidget(self.history)
        splitter.addWidget(history_box)

        detail = QtWidgets.QWidget()
        detail_layout = QtWidgets.QVBoxLayout(detail)
        detail_layout.setContentsMargins(0, 0, 0, 0)
        row = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel("No run yet.")
        self.status.setWordWrap(True)
        row.addWidget(self.status, 1)
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        row.addWidget(self.cancel_button)
        detail_layout.addLayout(row)
        self.table = QtWidgets.QTableWidget(0, len(STATUS_COLUMNS))
        self.table.setHorizontalHeaderLabels(STATUS_COLUMNS)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setWordWrap(True)
        detail_layout.addWidget(self.table, 1)
        detail_layout.addWidget(QtWidgets.QLabel("Log"))
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setLineWrapMode(QtWidgets.QPlainTextEdit.WidgetWidth)
        detail_layout.addWidget(self.log, 2)
        splitter.addWidget(detail)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.refresh_table)

    def set_project(self, project):
        self.project = project
        self.refresh_history()

    def attach(self, controller):
        self.controller = controller
        controller.started.connect(self._started)
        controller.log.connect(self._append_log)
        controller.finished.connect(self._finished)

    def _started(self, name):
        self.log.clear()
        self.status.setText("Running {}".format(name))
        self.cancel_button.setEnabled(True)
        self.timer.start()

    def _append_log(self, line):
        self.log.appendPlainText(line)

    def _finished(self, code):
        self.timer.stop()
        self.cancel_button.setEnabled(False)
        self.refresh_table()
        state = {0: "finished", 130: "cancelled"}.get(code, "failed")
        self.status.setText("{} {} (exit status {})".format(self.controller.name, state, code))
        self.refresh_history()

    def cancel(self):
        if self.controller is not None:
            self.controller.cancel()
            self.status.setText("Cancelling {}".format(self.controller.name))

    def refresh_table(self):
        if self.controller is None:
            return
        self._fill_table(self.controller.step_records())

    def _fill_table(self, records):
        self.table.setRowCount(len(records))
        for row, record in enumerate(records):
            cells = [record.get("step", ""), record.get("tool", ""), record.get("status", ""),
                     "" if record.get("exit_code") is None else str(record.get("exit_code")),
                     _fmt(record.get("peak_rss_kb"), 1 / 1024.0, 1), _fmt(record.get("cpu_seconds")),
                     _fmt(record.get("wall_seconds"))]
            for column, text in enumerate(cells):
                self.table.setItem(row, column, QtWidgets.QTableWidgetItem(text))

    def refresh_history(self):
        self.history.blockSignals(True)
        self.history.clear()
        if self.project is not None:
            for entry in self.project.history():
                state = {0: "finished", 130: "cancelled"}.get(entry.get("exit_code"), "failed")
                item = QtWidgets.QListWidgetItem("{}\n{} ({})".format(entry.get("name", ""),
                                                                     entry.get("started", ""), state))
                item.setData(QtCore.Qt.UserRole, entry)
                self.history.addItem(item)
        self.history.blockSignals(False)

    def _history_selected(self, item, _previous=None):
        if item is None or (self.controller is not None and self.controller.running):
            return
        entry = item.data(QtCore.Qt.UserRole)
        self._fill_table(entry.get("steps") or [])
        lines = []
        for step in entry.get("steps") or []:
            for stream in ("stdout.txt", "stderr.txt"):
                path = os.path.join(entry["dir"], "steps", step["step"], stream)
                if os.path.isfile(path):
                    with open(path, errors="replace") as handle:
                        text = handle.read()
                    if text.strip():
                        lines.append("== {} {}\n{}".format(step["step"], stream, text.rstrip()))
        if entry.get("error"):
            lines.append("error: " + entry["error"])
        self.log.setPlainText("\n".join(lines))
        self.status.setText("{} (exit status {}), started {}".format(
            entry.get("name", ""), entry.get("exit_code"), entry.get("started", "")))
