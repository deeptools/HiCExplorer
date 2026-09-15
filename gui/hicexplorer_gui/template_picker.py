"""The "new workflow from template" dialog (cpp/PLAN.md 10.7)."""

import yaml
from PySide6 import QtCore, QtWidgets

from .templates import TemplateError, load_templates


class TemplatePicker(QtWidgets.QDialog):
    """Lists the templates, shows each one's parameters with their defaults and
    descriptions and its unavailable steps with the reason, and instantiates
    the chosen one. ``workflow`` holds the result after ``accept``."""

    def __init__(self, entries, parent=None, templates=None):
        super().__init__(parent)
        self.setWindowTitle("New workflow from template")
        self.entries = entries
        self.templates = templates if templates is not None else load_templates()
        self.workflow = None
        self.fields = {}

        layout = QtWidgets.QHBoxLayout(self)
        self.list = QtWidgets.QListWidget()
        self.list.addItems([t.title for t in self.templates])
        self.list.setMaximumWidth(220)
        self.list.currentRowChanged.connect(self._show)
        layout.addWidget(self.list)

        right = QtWidgets.QVBoxLayout()
        self.description = QtWidgets.QLabel()
        self.description.setWordWrap(True)
        right.addWidget(self.description)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.form_host = QtWidgets.QWidget()
        self.form = QtWidgets.QFormLayout(self.form_host)
        self.form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.form.setRowWrapPolicy(QtWidgets.QFormLayout.WrapLongRows)
        scroll.setWidget(self.form_host)
        right.addWidget(scroll, 1)
        self.unavailable = QtWidgets.QLabel()
        self.unavailable.setWordWrap(True)
        self.unavailable.setObjectName("unavailable_steps")
        right.addWidget(self.unavailable)
        name_row = QtWidgets.QHBoxLayout()
        name_row.addWidget(QtWidgets.QLabel("Workflow name:"))
        self.name = QtWidgets.QLineEdit()
        name_row.addWidget(self.name, 1)
        right.addLayout(name_row)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        right.addWidget(self.status)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel)
        self.create_button = buttons.addButton("Create", QtWidgets.QDialogButtonBox.AcceptRole)
        buttons.accepted.connect(self.create)
        buttons.rejected.connect(self.reject)
        right.addWidget(buttons)
        layout.addLayout(right, 1)
        if self.templates:
            self.list.setCurrentRow(0)

    @property
    def template(self):
        row = self.list.currentRow()
        return self.templates[row] if 0 <= row < len(self.templates) else None

    def _show(self, _row):
        template = self.template
        while self.form.rowCount():
            self.form.removeRow(0)
        self.fields = {}
        if template is None:
            return
        self.description.setText("<b>{}</b><br>{}".format(template.title, template.description))
        self.name.setText(template.name)
        for parameter in template.parameters:
            field = QtWidgets.QWidget()
            row = QtWidgets.QHBoxLayout(field)
            row.setContentsMargins(0, 0, 0, 0)
            edit = QtWidgets.QLineEdit()
            if "default" in parameter:
                edit.setText(self._text(parameter["default"]))
            else:
                edit.setPlaceholderText("required")
            edit.setToolTip(parameter.get("description", ""))
            row.addWidget(edit, 1)
            if parameter["kind"] in ("file", "files", "directory"):
                browse = QtWidgets.QPushButton("...")
                browse.setFixedWidth(32)
                browse.clicked.connect(lambda _=False, p=parameter, e=edit: self._browse(p, e))
                row.addWidget(browse)
            label = QtWidgets.QLabel("{}:".format(parameter["name"]))
            label.setToolTip(parameter.get("description", ""))
            self.form.addRow(label, field)
            hint = QtWidgets.QLabel(parameter.get("description", ""))
            hint.setWordWrap(True)
            hint.setStyleSheet("color: gray;")
            self.form.addRow("", hint)
            self.fields[parameter["name"]] = (parameter, edit)
        missing = template.unavailable_steps(self.entries)
        if missing:
            self.unavailable.setText("<b>Unavailable steps</b> (not part of the workflow):<br>" + "<br>".join(
                "{} ({}): {}".format(tool, purpose, reason) for tool, purpose, reason in missing))
        else:
            self.unavailable.setText("All steps are available.")
        self.status.clear()

    @staticmethod
    def _text(value):
        if isinstance(value, (list, dict)):
            return yaml.safe_dump(value, default_flow_style=True).strip()
        return str(value)

    def _browse(self, parameter, edit):
        if parameter["kind"] == "directory":
            path = QtWidgets.QFileDialog.getExistingDirectory(self, parameter["name"])
            paths = [path] if path else []
        elif parameter["kind"] == "files":
            paths, _ = QtWidgets.QFileDialog.getOpenFileNames(self, parameter["name"])
        else:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, parameter["name"])
            paths = [path] if path else []
        if paths:
            edit.setText(self._text(paths if parameter["kind"] == "files" else paths[0]))

    def values(self):
        out = {}
        for name, (parameter, edit) in self.fields.items():
            text = edit.text().strip()
            if not text:
                continue
            if parameter["kind"] in ("file", "directory"):
                out[name] = text
            else:
                try:
                    value = yaml.safe_load(text)
                except yaml.YAMLError:
                    value = text
                out[name] = [value] if parameter["kind"] == "files" and not isinstance(value, list) else value
        return out

    def create(self):
        template = self.template
        if template is None:
            return
        try:
            self.workflow = template.instantiate(self.values(), name=self.name.text().strip() or None)
        except TemplateError as exc:
            self.status.setText(str(exc))
            return
        self.accept()
