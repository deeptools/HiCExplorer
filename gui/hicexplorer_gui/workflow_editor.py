"""A minimal workflow editor: add steps, connect outputs to inputs, save, validate, run.

Each step is a tool form. An output file field names a path relative to the
workflow's work directory and becomes a declared output of the step; an input
file field can be connected to any output of an earlier step, which writes a
``${steps.<id>.outputs.<name>}`` reference into the workflow file.
"""

import os
import re

import yaml
from PySide6 import QtCore, QtWidgets

from .forms import ToolForm, is_list_arg
from .project import safe_name
from .workflow.model import WorkflowError, load_workflow
from .workflow.validate import validate_workflow

_ID_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")


class StepEditor(QtWidgets.QWidget):
    def __init__(self, editor, step_id, spec, parent=None):
        super().__init__(parent)
        self.editor = editor
        self.spec = spec
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Step id"))
        self.id_edit = QtWidgets.QLineEdit(step_id)
        self.id_edit.textChanged.connect(editor.steps_changed)
        row.addWidget(self.id_edit, 1)
        row.addWidget(QtWidgets.QLabel("Threads"))
        self.threads = QtWidgets.QSpinBox()
        self.threads.setRange(1, 1024)
        row.addWidget(self.threads)
        layout.addLayout(row)
        self.form = ToolForm(spec, link_provider=lambda: editor.upstream_outputs(self))
        self.form.changed.connect(editor.outputs_changed)
        layout.addWidget(self.form, 1)

    @property
    def step_id(self):
        return self.id_edit.text().strip()

    def outputs(self):
        """{output name: relative path} from the output file fields that are set."""
        out = {}
        values = self.form.values()
        for arg in self.form.arguments():
            info = arg.get("file") or {}
            value = values.get(arg["dest"])
            if info.get("role") != "output" or value is None:
                continue
            items = value if isinstance(value, list) else [value]
            for index, path in enumerate(items):
                if isinstance(path, str) and not path.startswith("${"):
                    name = arg["dest"] if len(items) == 1 else "{}_{}".format(arg["dest"], index + 1)
                    out[name] = path
        return out

    def to_data(self):
        values = self.form.values()
        args = {}
        outputs = self.outputs()
        for arg in self.form.arguments():
            dest = arg["dest"]
            if dest not in values:
                continue
            value = values[dest]
            info = arg.get("file") or {}
            if info.get("role") == "output":
                items = value if isinstance(value, list) else [value]
                refs = []
                for index, path in enumerate(items):
                    name = dest if len(items) == 1 else "{}_{}".format(dest, index + 1)
                    refs.append("${outputs.%s}" % name if name in outputs else path)
                value = refs if is_list_arg(arg) else refs[0]
            args[dest] = value
        data = {"id": self.step_id, "tool": self.spec.tool}
        if self.form.subcommand() is not None:
            data["subcommand"] = self.form.subcommand()
        if self.threads.value() > 1:
            data["threads"] = self.threads.value()
        data["args"] = args
        if outputs:
            data["outputs"] = outputs
        return data


class WorkflowEditor(QtWidgets.QWidget):
    run_requested = QtCore.Signal(str, str)  # workflow path, name

    def __init__(self, parent=None):
        super().__init__(parent)
        self.project = None
        self.entries = []
        self.loader = None
        self.step_editors = []
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Workflow name"))
        self.name_edit = QtWidgets.QLineEdit("workflow")
        top.addWidget(self.name_edit, 1)
        top.addWidget(QtWidgets.QLabel("Thread budget"))
        self.budget = QtWidgets.QSpinBox()
        self.budget.setRange(1, 1024)
        self.budget.setValue(4)
        top.addWidget(self.budget)
        outer.addLayout(top)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(QtWidgets.QLabel("Steps"))
        self.step_list = QtWidgets.QListWidget()
        self.step_list.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.step_list.setWordWrap(True)
        self.step_list.setTextElideMode(QtCore.Qt.ElideNone)
        self.step_list.currentRowChanged.connect(self._select_step)
        left_layout.addWidget(self.step_list, 1)
        self.tool_combo = QtWidgets.QComboBox()
        self.tool_combo.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.tool_combo.setMinimumContentsLength(10)
        left_layout.addWidget(self.tool_combo)
        buttons = QtWidgets.QHBoxLayout()
        self.add_button = QtWidgets.QPushButton("Add step")
        self.add_button.clicked.connect(lambda: self.add_step(self.tool_combo.currentText()))
        buttons.addWidget(self.add_button)
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.remove_button.clicked.connect(self.remove_current)
        buttons.addWidget(self.remove_button)
        left_layout.addLayout(buttons)
        splitter.addWidget(left)
        self.stack = QtWidgets.QStackedWidget()
        self.placeholder = QtWidgets.QLabel("Add a step to edit its arguments.")
        self.placeholder.setWordWrap(True)
        self.stack.addWidget(self.placeholder)
        splitter.addWidget(self.stack)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)
        outer.addWidget(splitter, 1)

        actions = QtWidgets.QHBoxLayout()
        self.save_button = QtWidgets.QPushButton("Save")
        self.save_button.clicked.connect(self.save)
        actions.addWidget(self.save_button)
        self.validate_button = QtWidgets.QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate)
        actions.addWidget(self.validate_button)
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run)
        actions.addWidget(self.run_button)
        actions.addStretch(1)
        outer.addLayout(actions)
        self.messages = QtWidgets.QListWidget()
        self.messages.setWordWrap(True)
        self.messages.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.messages.setMaximumHeight(120)
        outer.addWidget(self.messages)

    # -- context --------------------------------------------------------
    def set_context(self, project, entries, loader):
        self.project = project
        self.entries = entries
        self.loader = loader
        self.tool_combo.clear()
        for entry in entries:
            if entry.available:
                self.tool_combo.addItem(entry.name)

    def _spec(self, tool):
        for entry in self.entries:
            if entry.name == tool and entry.available:
                return entry.spec
        return None

    # -- steps ----------------------------------------------------------
    def add_step(self, tool, step_id=None):
        spec = self._spec(tool)
        if spec is None:
            self.show_messages(["error: tool {} is not available".format(tool)])
            return None
        existing = {e.step_id for e in self.step_editors}
        if step_id is None:
            base = re.sub(r"[^A-Za-z0-9_]", "_", tool)
            step_id, n = base, 2
            while step_id in existing:
                step_id, n = "{}_{}".format(base, n), n + 1
        editor = StepEditor(self, step_id, spec)
        self.step_editors.append(editor)
        self.stack.addWidget(editor)
        self.steps_changed()
        self.step_list.setCurrentRow(len(self.step_editors) - 1)
        return editor

    def remove_current(self):
        row = self.step_list.currentRow()
        if row < 0:
            return
        editor = self.step_editors.pop(row)
        self.stack.removeWidget(editor)
        editor.deleteLater()
        self.steps_changed()

    def _select_step(self, row):
        if 0 <= row < len(self.step_editors):
            self.step_editors[row].form.refresh_links()
            self.stack.setCurrentWidget(self.step_editors[row])
        else:
            self.stack.setCurrentWidget(self.placeholder)

    def steps_changed(self, *_):
        row = self.step_list.currentRow()
        self.step_list.blockSignals(True)
        self.step_list.clear()
        for editor in self.step_editors:
            self.step_list.addItem("{} ({})".format(editor.step_id, editor.spec.tool))
        self.step_list.blockSignals(False)
        if self.step_editors:
            self.step_list.setCurrentRow(min(max(row, 0), len(self.step_editors) - 1))
        self.outputs_changed()

    def outputs_changed(self, *_):
        for editor in self.step_editors:
            editor.form.refresh_links()

    def upstream_outputs(self, editor):
        out = []
        for other in self.step_editors:
            if other is editor:
                break
            for name in other.outputs():
                out.append(("{}: {}".format(other.step_id, name),
                            "${steps.%s.outputs.%s}" % (other.step_id, name)))
        return out

    # -- document -------------------------------------------------------
    def to_data(self):
        return {"version": 1, "name": self.name_edit.text().strip() or "workflow",
                "threads": self.budget.value(), "steps": [e.to_data() for e in self.step_editors]}

    def workflow_path(self):
        return os.path.join(self.project.workflows_dir, safe_name(self.name_edit.text().strip() or "workflow") + ".yaml")

    def save(self):
        if self.project is None:
            self.show_messages(["error: open or create a project first"])
            return None
        problems = ["error: step id {!r} is not valid".format(e.step_id)
                    for e in self.step_editors if not _ID_RE.match(e.step_id)]
        if problems:
            self.show_messages(problems)
            return None
        os.makedirs(self.project.workflows_dir, exist_ok=True)
        path = self.workflow_path()
        with open(path, "w") as handle:
            yaml.safe_dump(self.to_data(), handle, sort_keys=False)
        self.show_messages(["saved {}".format(path)])
        return path

    def load(self, path):
        with open(path) as handle:
            data = yaml.safe_load(handle) or {}
        for editor in list(self.step_editors):
            self.stack.removeWidget(editor)
            editor.deleteLater()
        self.step_editors = []
        self.name_edit.setText(str(data.get("name") or os.path.splitext(os.path.basename(path))[0]))
        if isinstance(data.get("threads"), int):
            self.budget.setValue(data["threads"])
        for raw in data.get("steps") or []:
            editor = self.add_step(raw.get("tool"), raw.get("id"))
            if editor is None:
                continue
            editor.threads.setValue(int(raw.get("threads") or 1))
            outputs = raw.get("outputs") or {}
            values = {}
            for dest, value in (raw.get("args") or {}).items():
                values[dest] = _expand_outputs(value, outputs)
            editor.form.set_values(values, raw.get("subcommand"))
        self.steps_changed()

    def validate(self):
        path = self.save()
        if path is None:
            return None
        try:
            workflow = load_workflow(path, self.project.workdir_for(self.name_edit.text().strip() or "workflow"))
        except WorkflowError as exc:
            self.show_messages(["error: {}".format(exc)])
            return ["error: {}".format(exc)]
        messages = [str(m) for m in validate_workflow(workflow, self.loader)]
        self.show_messages(messages or ["workflow {} is valid ({} steps)".format(workflow.name, len(workflow.steps))])
        return messages

    def run(self):
        messages = self.validate()
        if messages is None or any(m.startswith("error") for m in messages):
            return
        self.run_requested.emit(self.workflow_path(), self.name_edit.text().strip() or "workflow")

    def show_messages(self, lines):
        self.messages.clear()
        for line in lines:
            item = QtWidgets.QListWidgetItem(line)
            if line.startswith("error"):
                item.setForeground(QtCore.Qt.red)
            self.messages.addItem(item)


def _expand_outputs(value, outputs):
    def one(v):
        if isinstance(v, str):
            match = re.fullmatch(r"\$\{outputs\.([^}]*)\}", v)
            if match and match.group(1) in outputs:
                return outputs[match.group(1)]
        return v
    if isinstance(value, list):
        return [one(v) for v in value]
    return one(value)
