"""Forms generated from a tool specification (``--help-json``).

Every settable argument gets a field: a check box for store_true and
store_false, a drop-down for single-valued choices, a text box otherwise
(several values separated by spaces, quoted where needed), and a file field
with a picker filtered by the argument's role and formats. Required
arguments, types, choices, nargs and mutually exclusive groups are checked
against the specification; errors appear next to the field and in a summary
line, never in a popup.
"""

import shlex

from PySide6 import QtCore, QtWidgets

from .workflow.cmdline import argv_for_values, shell_join
from .workflow.spec import NON_SETTABLE_ACTIONS
from .workflow.validate import _check_value

FORMAT_GLOBS = {
    "cool": ["*.cool"], "mcool": ["*.mcool"], "h5": ["*.h5"], "hic": ["*.hic"],
    "bed": ["*.bed", "*.bed.gz"], "bedgraph": ["*.bedgraph", "*.bg"], "bigwig": ["*.bw", "*.bigwig"],
    "bam": ["*.bam"], "sam": ["*.sam"], "txt": ["*.txt"], "tsv": ["*.tsv"], "npz": ["*.npz"],
    "png": ["*.png"], "pdf": ["*.pdf"], "svg": ["*.svg"], "fasta": ["*.fa", "*.fasta", "*.fa.gz"],
    "gff": ["*.gff", "*.gtf"], "hicpro": ["*.matrix"], "homer": ["*.homer", "*.gz"],
    "ginteractions": ["*.tsv", "*.ginteractions"], "broadpeak": ["*.broadPeak"], "narrowpeak": ["*.narrowPeak"],
}

ERROR_STYLE = "color: #b00020;"
HELP_STYLE = "color: palette(mid);"


def name_filter(formats):
    if not formats:
        return "All files (*)"
    globs = []
    for fmt in formats:
        for glob in FORMAT_GLOBS.get(fmt.lower(), ["*." + fmt]):
            if glob not in globs:
                globs.append(glob)
    return "{} ({});;All files (*)".format(", ".join(formats), " ".join(globs))


def is_list_arg(arg):
    nargs = arg.get("nargs")
    return nargs in ("+", "*") or (isinstance(nargs, int) and not isinstance(nargs, bool))


def convert(text, type_name):
    """Python's int() and float() for int and float arguments; text otherwise."""
    if type_name == "int":
        return int(text)
    if type_name == "float":
        return float(text)
    return text


def display(value):
    if value is None:
        return "None"
    if isinstance(value, list):
        return " ".join(display(v) for v in value)
    return str(value)


class Field(QtWidgets.QWidget):
    """One argument: editor row, help text and an inline error line."""

    changed = QtCore.Signal()

    def __init__(self, arg, parent=None):
        super().__init__(parent)
        self.arg = arg
        self.dest = arg["dest"]
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 2, 0, 6)
        outer.setSpacing(2)
        flags = ", ".join(arg.get("flags") or []) or self.dest
        label_text = "<b>{}</b>".format(flags)
        if arg.get("required"):
            label_text += " (required)"
        if arg.get("cpp_only"):
            label_text += " <i>C++ only</i>"
        self.label = QtWidgets.QLabel(label_text)
        self.label.setWordWrap(True)
        self.label.setTextFormat(QtCore.Qt.RichText)
        outer.addWidget(self.label)
        self.row = QtWidgets.QHBoxLayout()
        self.row.setSpacing(4)
        outer.addLayout(self.row)
        help_text = (arg.get("help") or "").strip()
        if arg.get("note"):
            help_text = (help_text + " " if help_text else "") + "Note: " + arg["note"]
        self.help = QtWidgets.QLabel(help_text)
        self.help.setWordWrap(True)
        self.help.setStyleSheet(HELP_STYLE)
        self.help.setVisible(bool(help_text))
        outer.addWidget(self.help)
        self.error = QtWidgets.QLabel()
        self.error.setWordWrap(True)
        self.error.setStyleSheet(ERROR_STYLE)
        self.error.hide()
        outer.addWidget(self.error)

    def set_error(self, message):
        self.error.setText(message or "")
        self.error.setVisible(bool(message))

    # read() -> (present, value, error); set_value(value); clear()
    def read(self):
        raise NotImplementedError

    def set_value(self, value):
        raise NotImplementedError

    def clear(self):
        raise NotImplementedError


class FlagField(Field):
    def __init__(self, arg, parent=None):
        super().__init__(arg, parent)
        default = arg.get("default")
        text = "set (stores False)" if arg.get("action") == "store_false" else "set"
        self.box = QtWidgets.QCheckBox(text)
        self.box.toggled.connect(self.changed)
        self.row.addWidget(self.box)
        self.row.addStretch(1)
        self._default = default

    def read(self):
        if not self.box.isChecked():
            return False, None, None
        return True, (self.arg.get("action") != "store_false"), None

    def set_value(self, value):
        self.box.setChecked(value is not None and value != self._default)

    def clear(self):
        self.box.setChecked(False)


class ChoiceField(Field):
    def __init__(self, arg, parent=None):
        super().__init__(arg, parent)
        self.combo = QtWidgets.QComboBox()
        self.combo.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.combo.setMinimumContentsLength(8)
        self.combo.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        default = arg.get("default")
        self.combo.addItem("(not set: default {})".format(display(default)), userData=_UNSET)
        for choice in arg.get("choices") or []:
            if choice is None:
                continue
            self.combo.addItem(display(choice), userData=choice)
        self.combo.currentIndexChanged.connect(self.changed)
        self.row.addWidget(self.combo, 1)
        self._invalid = None

    def read(self):
        if self._invalid is not None:
            return True, self._invalid, None
        data = self.combo.currentData()
        if data is _UNSET:
            return False, None, None
        return True, data, None

    def set_value(self, value):
        self._invalid = None
        if value is None:
            self.combo.setCurrentIndex(0)
            return
        for index in range(1, self.combo.count()):
            if self.combo.itemData(index) == value:
                self.combo.setCurrentIndex(index)
                return
        # A value the drop-down cannot show (e.g. loaded from a file): keep it
        # so validation reports it.
        self._invalid = value
        self.changed.emit()

    def clear(self):
        self._invalid = None
        self.combo.setCurrentIndex(0)


class _Unset:
    def __repr__(self):
        return "<unset>"


_UNSET = _Unset()


class TextField(Field):
    def __init__(self, arg, parent=None):
        super().__init__(arg, parent)
        self.edit = QtWidgets.QLineEdit()
        self.edit.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.edit.setMinimumWidth(80)
        nargs = arg.get("nargs")
        hint = []
        if nargs == "+":
            hint.append("one or more values, separated by spaces")
        elif nargs == "*":
            hint.append("any number of values, separated by spaces")
        elif isinstance(nargs, int) and not isinstance(nargs, bool):
            hint.append("{} value{}, separated by spaces".format(nargs, "" if nargs == 1 else "s"))
        if arg.get("type") in ("int", "float"):
            hint.append(arg["type"])
        if arg.get("choices"):
            hint.append("one of " + ", ".join(display(c) for c in arg["choices"] if c is not None))
        default = arg.get("default")
        if default is not None:
            hint.append("default " + display(default))
        self.edit.setPlaceholderText("; ".join(hint))
        self.edit.textChanged.connect(self.changed)
        self.row.addWidget(self.edit, 1)

    def tokens(self):
        text = self.edit.text().strip()
        if not text:
            return []
        if is_list_arg(self.arg):
            return shlex.split(text)
        return [text]

    def read(self):
        text = self.edit.text().strip()
        if not text:
            return False, None, None
        try:
            tokens = self.tokens()
        except ValueError as exc:
            return True, None, "cannot split the values: {}".format(exc)
        type_name = self.arg.get("type")
        values = []
        for token in tokens:
            if token.startswith("${"):
                values.append(token)
                continue
            try:
                values.append(convert(token, type_name))
            except ValueError:
                return True, None, "invalid {} value: {!r}".format(type_name, token)
        if is_list_arg(self.arg):
            return True, values, None
        return True, values[0], None

    def set_value(self, value):
        if value is None:
            self.edit.clear()
        elif isinstance(value, list):
            self.edit.setText(" ".join(shlex.quote(display(v)) for v in value))
        else:
            self.edit.setText(display(value))

    def clear(self):
        self.edit.clear()


class FileField(TextField):
    """A text field with a picker filtered by role and formats, and, in the
    workflow editor, a drop-down of upstream step outputs to connect."""

    def __init__(self, arg, parent=None, link_provider=None, start_dir=None):
        super().__init__(arg, parent)
        info = arg.get("file") or {}
        self.role = info.get("role") or "input"
        self.kind = info.get("kind") or "file"
        self.formats = info.get("formats") or []
        self.start_dir = start_dir
        self.browse = QtWidgets.QPushButton("Browse")
        self.browse.clicked.connect(self.pick)
        self.row.addWidget(self.browse)
        self.link_provider = link_provider
        self.link = None
        if link_provider is not None and self.role == "input":
            self.link = QtWidgets.QComboBox()
            self.link.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
            self.link.setMinimumContentsLength(10)
            self.link.setToolTip("Connect the output of an earlier step")
            self.link.activated.connect(self._link_chosen)
            self.refresh_links()
            outer = self.layout()
            outer.insertWidget(2, self.link)
        if self.formats:
            self.label.setText(self.label.text() + " [{}: {}]".format(self.role, ", ".join(self.formats)))
        else:
            self.label.setText(self.label.text() + " [{}]".format(self.role))

    def file_filter(self):
        return name_filter(self.formats)

    def refresh_links(self):
        if self.link is None:
            return
        current = self.edit.text()
        self.link.blockSignals(True)
        self.link.clear()
        self.link.addItem("Connect an upstream output...", userData=None)
        for label, reference in self.link_provider():
            self.link.addItem(label, userData=reference)
        self.link.blockSignals(False)
        self.edit.setText(current)

    def _link_chosen(self, index):
        reference = self.link.itemData(index)
        if reference is None:
            return
        if is_list_arg(self.arg):
            tokens = self.tokens() if self.edit.text().strip() else []
            tokens.append(reference)
            self.set_value(tokens)
        else:
            self.edit.setText(reference)

    def pick(self):
        start = self.start_dir() if callable(self.start_dir) else (self.start_dir or "")
        dialog_filter = self.file_filter()
        multi = is_list_arg(self.arg)
        if self.kind == "directory":
            path = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose a directory", start)
            paths = [path] if path else []
        elif self.role == "output" or self.kind == "prefix":
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Choose an output file", start, dialog_filter)
            paths = [path] if path else []
        elif multi:
            paths, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "Choose input files", start, dialog_filter)
        else:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose an input file", start, dialog_filter)
            paths = [path] if path else []
        if not paths:
            return
        if multi:
            tokens = self.tokens() if self.edit.text().strip() else []
            self.set_value(tokens + paths)
        else:
            self.edit.setText(paths[0])


def make_field(arg, link_provider=None, start_dir=None):
    action = arg.get("action") or "store"
    if action in ("store_true", "store_false", "store_const"):
        return FlagField(arg)
    if arg.get("file"):
        return FileField(arg, link_provider=link_provider, start_dir=start_dir)
    choices = [c for c in (arg.get("choices") or []) if c is not None]
    if choices and not is_list_arg(arg):
        return ChoiceField(arg)
    return TextField(arg)


class ToolForm(QtWidgets.QWidget):
    """The form of one tool, or of one subcommand of a tool."""

    changed = QtCore.Signal()

    def __init__(self, spec, parent=None, link_provider=None, start_dir=None):
        super().__init__(parent)
        self.spec = spec
        self.link_provider = link_provider
        self.start_dir = start_dir
        self.fields = {}
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        title = QtWidgets.QLabel("<b>{}</b> {}".format(spec.tool, spec.version))
        layout.addWidget(title)
        if spec.description:
            description = QtWidgets.QLabel(spec.description.strip())
            description.setWordWrap(True)
            layout.addWidget(description)
        self.subcommand_combo = None
        if spec.has_subcommands:
            row = QtWidgets.QHBoxLayout()
            row.addWidget(QtWidgets.QLabel("Subcommand"))
            self.subcommand_combo = QtWidgets.QComboBox()
            if not spec.subcommand_required:
                self.subcommand_combo.addItem("(none)", userData=None)
            for name in spec.commands:
                self.subcommand_combo.addItem(name, userData=name)
            self.subcommand_combo.currentIndexChanged.connect(self._rebuild)
            row.addWidget(self.subcommand_combo, 1)
            layout.addLayout(row)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        self.summary.setStyleSheet(ERROR_STYLE)
        self.summary.hide()
        layout.addWidget(self.summary)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        layout.addWidget(self.scroll, 1)
        self._rebuild()

    # -- structure ------------------------------------------------------
    def subcommand(self):
        if self.subcommand_combo is None:
            return None
        return self.subcommand_combo.currentData()

    def set_subcommand(self, name):
        if self.subcommand_combo is None:
            return
        index = self.subcommand_combo.findData(name)
        if index >= 0:
            self.subcommand_combo.setCurrentIndex(index)

    def arguments(self):
        return [a for a in self.spec.arguments(self.subcommand())
                if (a.get("action") or "store") not in NON_SETTABLE_ACTIONS]

    def _groups(self):
        groups = [(g.get("title") or "", g.get("arguments") or []) for g in self.spec.data.get("groups") or []]
        sub = self.subcommand()
        if sub is not None:
            groups += [("{}: {}".format(sub, g.get("title") or ""), g.get("arguments") or [])
                       for g in self.spec.commands[sub].get("groups") or []]
        return groups

    def _rebuild(self, *_):
        old = self.values() if self.fields else {}
        container = QtWidgets.QWidget()
        column = QtWidgets.QVBoxLayout(container)
        column.setContentsMargins(4, 4, 4, 4)
        self.fields = {}
        for title, args in self._groups():
            settable = [a for a in args if (a.get("action") or "store") not in NON_SETTABLE_ACTIONS]
            if not settable:
                continue
            box = QtWidgets.QGroupBox(title)
            box_layout = QtWidgets.QVBoxLayout(box)
            for arg in settable:
                field = make_field(arg, self.link_provider, self.start_dir)
                field.changed.connect(self.changed)
                self.fields[arg["dest"]] = field
                box_layout.addWidget(field)
            column.addWidget(box)
        exclusive = self.spec.mutually_exclusive(self.subcommand())
        if exclusive:
            text = "Mutually exclusive: " + "; ".join(
                " / ".join(g.get("dests") or []) + (" (one required)" if g.get("required") else "")
                for g in exclusive)
            note = QtWidgets.QLabel(text)
            note.setWordWrap(True)
            column.addWidget(note)
        column.addStretch(1)
        self.scroll.setWidget(container)
        for dest, value in old.items():
            if dest in self.fields:
                self.fields[dest].set_value(value)
        self.changed.emit()

    def refresh_links(self):
        for field in self.fields.values():
            if isinstance(field, FileField):
                field.refresh_links()

    # -- values ---------------------------------------------------------
    def read(self):
        """(values of the set fields, errors of unreadable fields)."""
        values, errors = {}, {}
        for dest, field in self.fields.items():
            present, value, error = field.read()
            if error:
                errors[dest] = error
            elif present:
                values[dest] = value
        return values, errors

    def values(self):
        return self.read()[0]

    def set_values(self, values, subcommand=None):
        if subcommand is not None:
            self.set_subcommand(subcommand)
        for field in self.fields.values():
            field.clear()
        for dest, value in values.items():
            if dest in self.fields:
                self.fields[dest].set_value(value)

    def check_values(self, values):
        """Errors {dest: message} of a set of values against the specification."""
        errors = {}
        by_dest = self.spec.arguments_by_dest(self.subcommand())
        for dest, value in values.items():
            arg = by_dest.get(dest)
            if arg is None:
                errors[dest] = "unknown argument"
                continue
            messages = []
            if isinstance(value, str) and value.startswith("${"):
                continue
            if isinstance(value, list) and any(isinstance(v, str) and v.startswith("${") for v in value):
                continue
            _check_value(messages.append, "", dest, arg, arg.get("action") or "store", value)
            if messages:
                errors[dest] = messages[0].split(": ", 2)[-1]
        for dest, arg in by_dest.items():
            if (arg.get("action") or "store") in NON_SETTABLE_ACTIONS:
                continue
            if arg.get("required") and values.get(dest) is None and dest not in errors:
                errors[dest] = "required"
        for group in self.spec.mutually_exclusive(self.subcommand()):
            dests = group.get("dests") or []
            present = [d for d in dests if values.get(d) is not None]
            if len(present) > 1:
                for d in present:
                    others = ", ".join(x for x in present if x != d)
                    errors.setdefault(d, "not allowed together with {}".format(others))
            elif group.get("required") and not present:
                for d in dests:
                    errors.setdefault(d, "one of {} is required".format(", ".join(dests)))
        if self.spec.has_subcommands and self.spec.subcommand_required and self.subcommand() is None:
            errors["<subcommand>"] = "a subcommand is required"
        return errors

    def validate(self):
        """Checks the form, shows the errors inline and returns {dest: message}."""
        values, errors = self.read()
        for dest, message in self.check_values(values).items():
            errors.setdefault(dest, message)
        for dest, field in self.fields.items():
            field.set_error(errors.get(dest))
        if errors:
            self.summary.setText("{} problem{}: {}".format(
                len(errors), "" if len(errors) == 1 else "s",
                "; ".join("{}: {}".format(d, m) for d, m in sorted(errors.items()))))
            self.summary.show()
        else:
            self.summary.hide()
        return errors

    def argv(self, executable=None):
        return argv_for_values(self.spec, self.subcommand(), self.values(), executable)

    def command_line(self, executable=None):
        return shell_join(self.argv(executable))
