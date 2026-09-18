"""One open project: its own data, tool browser, runs, workflow editor,
matrix browser and dynamically opened tool-form and analysis-view tabs.

A ProjectTab is added as one top-level tab of the main window's project
QTabWidget (main_window.MainWindow). Everything a project needs is scoped to
its ProjectTab instance: nothing here is shared between two open projects
except the window-global C++ tool directory and the resulting tool catalog
(entries, loader), which the main window pushes into every open ProjectTab
with set_catalog() whenever Settings changes it.
"""

import os

from PySide6 import QtCore, QtGui, QtWidgets

from .catalog import entries_for_format
from .dataformats import FASTQ_MESSAGE, LOADABLE_FORMATS, detect_format, is_fastq_format, is_matrix_format
from .forms import ERROR_STYLE, FORMAT_GLOBS, ToolForm
from .runs import RunController, RunView
from .workflow_editor import WorkflowEditor


class DataFile:
    def __init__(self, path, fmt):
        self.path = path
        self.format = fmt

    @property
    def name(self):
        return os.path.basename(self.path)

    def __repr__(self):
        return "DataFile({!r}, format={!r})".format(self.path, self.format)


def _open_data_filter():
    globs = []
    for fmt in LOADABLE_FORMATS:
        for glob in FORMAT_GLOBS.get(fmt, ["*." + fmt]):
            if glob not in globs:
                globs.append(glob)
    return "Hi-C and alignment data ({});;All files (*)".format(" ".join(globs))


class ToolPage(QtWidgets.QWidget):
    """A tool form with validate, run and the command line it produces."""

    def __init__(self, project_tab, spec, parent=None):
        super().__init__(parent)
        self.owner = project_tab
        self.spec = spec
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        self.form = ToolForm(spec, start_dir=lambda: project_tab.project.path if project_tab.project else "")
        layout.addWidget(self.form, 1)
        self.command = QtWidgets.QLabel()
        self.command.setWordWrap(True)
        self.command.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        layout.addWidget(self.command)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        row = QtWidgets.QHBoxLayout()
        self.validate_button = QtWidgets.QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate)
        row.addWidget(self.validate_button)
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.clicked.connect(self.run)
        row.addWidget(self.run_button)
        row.addStretch(1)
        layout.addLayout(row)
        self.form.changed.connect(self.update_command)
        self.update_command()

    def update_command(self):
        self.command.setText("Command: " + self.form.command_line())

    def validate(self):
        errors = self.form.validate()
        self.status.setStyleSheet(ERROR_STYLE if errors else "")
        self.status.setText("" if not errors else "Fix the marked fields before running.")
        return errors

    def run(self):
        if self.validate():
            return False
        if self.owner.project is None:
            self.status.setStyleSheet(ERROR_STYLE)
            self.status.setText("Open or create a project first; runs are recorded in its history.")
            return False
        if self.owner.controller.running:
            self.status.setStyleSheet(ERROR_STYLE)
            self.status.setText("Another run is in progress in this project.")
            return False
        self.status.setStyleSheet("")
        self.status.setText("Running; see the Runs tab.")
        self.owner.show_runs()
        return self.owner.controller.run_tool(self.owner.project, self.owner.loader, self.spec,
                                              self.form.subcommand(), self.form.values())


class ProjectTab(QtWidgets.QWidget):
    """Everything one open project owns: data, tool catalog filter, runs,
    workflow editor, matrix browser and dynamic tool-form/analysis tabs."""

    def __init__(self, window, project, entries, loader, parent=None):
        super().__init__(parent)
        self.window = window
        self.project = project
        self.entries = entries
        self.filtered_entries = entries
        self.loader = loader
        self.data_files = []
        self.current_data = None
        self.controller = RunController(self)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(splitter)

        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        left_layout.setContentsMargins(2, 2, 2, 2)
        self._build_data_panel(left_layout)
        self._build_tools_panel(left_layout)
        self._build_project_panel(left_layout)
        splitter.addWidget(left)
        splitter.setStretchFactor(0, 0)
        left.setMinimumWidth(260)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabCloseRequested.connect(self._close_tab)
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(1, 1)

        self.run_view = RunView()
        self.run_view.attach(self.controller)
        self.run_view.set_project(project)
        self.editor = WorkflowEditor()
        self.editor.run_requested.connect(self.run_workflow)
        self.editor.set_context(project, entries, loader)
        self.browser = None
        self._fixed = []
        for widget, title in ((self.run_view, "Runs"), (self.editor, "Workflow editor")):
            self.tabs.addTab(widget, title)
            self._fixed.append(widget)
        self._add_browser()
        self.controller.finished.connect(lambda _code: self.refresh_project())

        self._apply_filter()
        self.refresh_project()

    # -- layout -----------------------------------------------------------
    def _add_browser(self):
        try:
            from .browser import MatrixBrowser
        except ImportError as exc:  # hicx_matrix or pyqtgraph missing
            page = QtWidgets.QLabel("The matrix browser needs pyqtgraph and the hicx_matrix module: {}".format(exc))
            page.setWordWrap(True)
            self.browser = None
        else:
            page = self.browser = MatrixBrowser()
        self.tabs.addTab(page, "Matrix browser")
        self._fixed.append(page)

    def _build_data_panel(self, layout):
        layout.addWidget(_bold_label("Data"))
        self.load_button = QtWidgets.QPushButton("Load data file...")
        self.load_button.clicked.connect(self.load_data_dialog)
        layout.addWidget(self.load_button)
        self.data_list = QtWidgets.QListWidget()
        self.data_list.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.data_list.setMaximumHeight(90)
        self.data_list.itemClicked.connect(self._data_item_clicked)
        layout.addWidget(self.data_list)
        row = QtWidgets.QHBoxLayout()
        self.remove_data_button = QtWidgets.QPushButton("Remove")
        self.remove_data_button.clicked.connect(self.remove_current_data)
        row.addWidget(self.remove_data_button)
        self.open_in_browser_button = QtWidgets.QPushButton("Open in Matrix browser")
        self.open_in_browser_button.clicked.connect(self._open_current_in_browser)
        row.addWidget(self.open_in_browser_button)
        layout.addLayout(row)
        self.data_status = QtWidgets.QLabel("No data loaded; showing every available tool.")
        self.data_status.setWordWrap(True)
        layout.addWidget(self.data_status)

    def _build_tools_panel(self, layout):
        layout.addWidget(_bold_label("Tools"))
        self.tool_filter = QtWidgets.QLineEdit()
        self.tool_filter.setPlaceholderText("Filter tools")
        self.tool_filter.textChanged.connect(self._fill_tools)
        layout.addWidget(self.tool_filter)
        self.tool_tree = QtWidgets.QTreeWidget()
        self.tool_tree.setHeaderHidden(True)
        self.tool_tree.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tool_tree.currentItemChanged.connect(self._tool_selected)
        self.tool_tree.itemDoubleClicked.connect(self._tool_opened)
        layout.addWidget(self.tool_tree, 3)
        self.tool_detail = QtWidgets.QLabel()
        self.tool_detail.setWordWrap(True)
        self.tool_detail.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        layout.addWidget(self.tool_detail, 1)

    def _build_project_panel(self, layout):
        layout.addWidget(_bold_label("Project"))
        self.project_tree = QtWidgets.QTreeWidget()
        self.project_tree.setHeaderHidden(True)
        self.project_tree.itemDoubleClicked.connect(self._project_item_opened)
        layout.addWidget(self.project_tree, 2)

    def _close_tab(self, index):
        if self.tabs.widget(index) in self._fixed:
            return
        widget = self.tabs.widget(index)
        self.tabs.removeTab(index)
        widget.deleteLater()

    def show_runs(self):
        self.tabs.setCurrentWidget(self.run_view)

    def show_editor(self):
        self.tabs.setCurrentWidget(self.editor)

    def show_browser(self):
        self.tabs.setCurrentIndex(self.tabs.indexOf(self._fixed[-1]))

    # -- catalog and data-driven filtering ---------------------------------
    def set_catalog(self, entries, loader):
        """Called by the main window when the window-global tool directory
        changes; the tool catalog is shared, the filter stays project-local."""
        self.entries = entries
        self.loader = loader
        self.editor.set_context(self.project, entries, loader)
        self._apply_filter()

    def load_data_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load a data file", self.project.path,
                                                         _open_data_filter())
        if path:
            self.load_data_file(path)

    def load_data_file(self, path):
        fmt = detect_format(path)
        entry = DataFile(path, fmt)
        self.data_files.append(entry)
        self.set_current_data(entry)
        if is_matrix_format(fmt):
            self._open_current_in_browser()
        return entry

    def set_current_data(self, entry):
        self.current_data = entry
        self._refresh_data_list()
        self._apply_filter()

    def remove_current_data(self):
        row = self.data_list.currentRow()
        if not (0 <= row < len(self.data_files)):
            return
        removed = self.data_files.pop(row)
        if self.current_data is removed:
            self.current_data = self.data_files[-1] if self.data_files else None
        self._refresh_data_list()
        self._apply_filter()

    def _data_item_clicked(self, item):
        index = self.data_list.row(item)
        if 0 <= index < len(self.data_files):
            self.set_current_data(self.data_files[index])

    def _refresh_data_list(self):
        self.data_list.blockSignals(True)
        self.data_list.clear()
        for entry in self.data_files:
            fmt_text = entry.format or "unrecognised"
            mark = "> " if entry is self.current_data else "  "
            item = QtWidgets.QListWidgetItem("{}{} ({})".format(mark, entry.name, fmt_text))
            item.setToolTip(entry.path)
            self.data_list.addItem(item)
            if entry is self.current_data:
                self.data_list.setCurrentItem(item)
        self.data_list.blockSignals(False)
        self.open_in_browser_button.setEnabled(
            self.current_data is not None and is_matrix_format(self.current_data.format) and self.browser is not None)

    def _open_current_in_browser(self):
        if self.current_data is None or self.browser is None:
            return
        if self.browser.open_matrix(self.current_data.path) is not None:
            self.show_browser()

    def _apply_filter(self):
        if self.current_data is None:
            self.filtered_entries = self.entries
            self.data_status.setText("No data loaded; showing every available tool.")
        elif self.current_data.format is None:
            self.filtered_entries = self.entries
            self.data_status.setText(
                "{} is not a recognised format; showing every available tool.".format(self.current_data.name))
        elif is_fastq_format(self.current_data.format):
            self.filtered_entries = []
            self.data_status.setText(
                "{}: {}.".format(self.current_data.name, FASTQ_MESSAGE))
        else:
            self.filtered_entries = entries_for_format(self.entries, self.current_data.format)
            if self.filtered_entries:
                self.data_status.setText("{}: {} tool(s) accept {} input.".format(
                    self.current_data.name, len(self.filtered_entries), self.current_data.format))
            else:
                self.data_status.setText("{}: no available tool accepts {} input.".format(
                    self.current_data.name, self.current_data.format))
        self._fill_tools()

    # -- tools --------------------------------------------------------------
    def _fill_tools(self, *_):
        text = self.tool_filter.text().strip().lower()
        self.tool_tree.clear()
        available = QtWidgets.QTreeWidgetItem(["Available"])
        missing = QtWidgets.QTreeWidgetItem(["Not available"])
        # With data loaded, the candidate list is already filtered to formats
        # the loaded file matches (entries_for_format only returns available
        # tools, so "Not available" stays empty); without data, fall back to
        # the full catalog, available and unavailable alike.
        candidates = self.filtered_entries if self.current_data is not None else self.entries
        for entry in candidates:
            if text and text not in entry.name.lower():
                continue
            item = QtWidgets.QTreeWidgetItem([entry.name])
            item.setData(0, QtCore.Qt.UserRole, entry.name)
            if entry.available:
                item.setToolTip(0, entry.spec.description or entry.name)
                available.addChild(item)
            else:
                item.setToolTip(0, entry.reason)
                item.setForeground(0, QtGui.QBrush(QtGui.QColor("gray")))
                missing.addChild(item)
        available.setText(0, "Available ({})".format(available.childCount()))
        missing.setText(0, "Not available ({})".format(missing.childCount()))
        self.tool_tree.addTopLevelItems([available, missing])
        available.setExpanded(True)
        missing.setExpanded(missing.childCount() > 0)
        metrics = self.tool_tree.fontMetrics()
        longest = max([metrics.horizontalAdvance(e.name) for e in self.entries] or [0])
        scroll = self.tool_tree.verticalScrollBar().sizeHint().width()
        self.tool_tree.setMinimumWidth(longest + 2 * self.tool_tree.indentation() + scroll + 16)

    def entry(self, name):
        for entry in self.entries:
            if entry.name == name:
                return entry
        return None

    def _tool_selected(self, item, _previous=None):
        name = item.data(0, QtCore.Qt.UserRole) if item is not None else None
        entry = self.entry(name) if name else None
        if entry is None:
            self.tool_detail.setText("")
        elif entry.available:
            self.tool_detail.setText("<b>{}</b><br>{}".format(entry.name, entry.spec.description))
        else:
            self.tool_detail.setText("<b>{}</b><br>Not available: {}".format(entry.name, entry.reason))

    def _tool_opened(self, item, _column=0):
        name = item.data(0, QtCore.Qt.UserRole)
        if name:
            self.open_tool_form(name)

    def open_tool_form(self, name, subcommand=None):
        entry = self.entry(name)
        if entry is None or not entry.available:
            self.window.statusBar().showMessage("{}: {}".format(name, entry.reason if entry else "unknown tool"))
            return None
        page = ToolPage(self, entry.spec)
        if subcommand is not None:
            page.form.set_subcommand(subcommand)
        index = self.tabs.addTab(page, name)
        self.tabs.setCurrentIndex(index)
        return page

    # -- project (workflows and history) -------------------------------------
    def refresh_project(self):
        self.project_tree.clear()
        if self.project is None:
            return
        root = QtWidgets.QTreeWidgetItem([self.project.name])
        workflows = QtWidgets.QTreeWidgetItem(["Workflows"])
        for path in self.project.workflows():
            item = QtWidgets.QTreeWidgetItem([os.path.basename(path)])
            item.setData(0, QtCore.Qt.UserRole, ("workflow", path))
            workflows.addChild(item)
        history = QtWidgets.QTreeWidgetItem(["History"])
        for entry in self.project.history():
            item = QtWidgets.QTreeWidgetItem(["{} ({})".format(entry.get("name"), entry.get("exit_code"))])
            item.setToolTip(0, entry.get("started", ""))
            item.setData(0, QtCore.Qt.UserRole, ("history", entry["dir"]))
            history.addChild(item)
        root.addChildren([workflows, history])
        self.project_tree.addTopLevelItem(root)
        root.setExpanded(True)
        workflows.setExpanded(True)
        self.run_view.refresh_history()

    def _project_item_opened(self, item, _column=0):
        data = item.data(0, QtCore.Qt.UserRole)
        if not data:
            return
        kind, path = data
        if kind == "workflow":
            self.editor.load(path)
            self.show_editor()
        elif not self.open_analysis_views(path):
            self.show_runs()

    # -- analysis views and templates (PLAN 10.6, 10.7) ----------------------
    def open_analysis_views(self, history_dir):
        """Opens a view for every step of a run whose tool has one."""
        from .analysis_views import records_from_history
        try:
            records = records_from_history(history_dir)
        except (OSError, ValueError) as exc:
            self.window.statusBar().showMessage("Cannot read the run: {}".format(exc))
            return []
        return [view for view in (self.open_analysis_view(record) for record in records) if view is not None]

    def open_analysis_view(self, record):
        from .analysis_views import open_view
        if self.loader is None:
            self.window.statusBar().showMessage("Set the C++ tool directory first (File > Settings).")
            return None
        try:
            view = open_view(record, self.loader)
        except Exception as exc:  # noqa: BLE001  a view that cannot read its data says why
            self.window.statusBar().showMessage("Cannot open the {} view: {}".format(record.tool, exc))
            return None
        view.navigate.connect(self.navigate_browser)
        self.tabs.addTab(view, "{}: {}".format(view.title, record.name))
        self.tabs.setCurrentWidget(view)
        return view

    def navigate_browser(self, matrix, region):
        """Shows region in this project's matrix browser, opening matrix
        first when it is not the matrix already open."""
        if self.browser is None:
            return False
        current = self.browser.sources[0]
        if matrix and (current is None or os.path.abspath(current.path) != os.path.abspath(matrix)):
            if self.browser.open_matrix(matrix) is None:
                return False
        ok = self.browser.goto(region)
        self.show_browser()
        return ok

    def new_from_template(self):
        from .template_picker import TemplatePicker
        dialog = TemplatePicker(self.entries, self)
        if dialog.exec() != QtWidgets.QDialog.Accepted or dialog.workflow is None:
            return None
        return self.create_workflow(dialog.workflow)

    def create_workflow(self, workflow):
        """Saves a workflow document in the project and opens it in the editor."""
        import yaml
        from .project import safe_name
        path = os.path.join(self.project.workflows_dir, safe_name(workflow["name"]) + ".yaml")
        with open(path, "w") as handle:
            yaml.safe_dump(workflow, handle, sort_keys=False)
        self.refresh_project()
        self.editor.load(path)
        self.show_editor()
        return path

    def run_workflow(self, path, name):
        if self.project is None or self.controller.running:
            return False
        self.show_runs()
        return self.controller.run_workflow(self.project, self.loader, path, name)


def _bold_label(text):
    label = QtWidgets.QLabel("<b>{}</b>".format(text))
    return label
