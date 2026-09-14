"""The main window: projects, the tool browser, tool forms, runs, the workflow editor and settings."""

import os

from PySide6 import QtCore, QtGui, QtWidgets

from .catalog import tool_entries
from .forms import ERROR_STYLE, ToolForm
from .project import Project, ProjectError
from .runs import RunController, RunView
from .workflow_editor import WorkflowEditor


class ToolPage(QtWidgets.QWidget):
    """A tool form with validate, run and the command line it produces."""

    def __init__(self, window, spec, parent=None):
        super().__init__(parent)
        self.window = window
        self.spec = spec
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        self.form = ToolForm(spec, start_dir=lambda: window.project.path if window.project else "")
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
        if self.window.project is None:
            self.status.setStyleSheet(ERROR_STYLE)
            self.status.setText("Open or create a project first; runs are recorded in its history.")
            return False
        if self.window.controller.running:
            self.status.setStyleSheet(ERROR_STYLE)
            self.status.setText("Another run is in progress.")
            return False
        self.status.setStyleSheet("")
        self.status.setText("Running; see the Runs tab.")
        self.window.show_runs()
        return self.window.controller.run_tool(self.window.project, self.window.loader, self.spec,
                                               self.form.subcommand(), self.form.values())


class SettingsPage(QtWidgets.QWidget):
    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(QtWidgets.QLabel("<b>C++ tool directory</b>"))
        hint = QtWidgets.QLabel("The directory holding the HiCExplorer C++ tool executables "
                                "(for example the tools directory of a build).")
        hint.setWordWrap(True)
        layout.addWidget(hint)
        row = QtWidgets.QHBoxLayout()
        self.edit = QtWidgets.QLineEdit(window.settings.tools_dir)
        row.addWidget(self.edit, 1)
        browse = QtWidgets.QPushButton("Browse")
        browse.clicked.connect(self.browse)
        row.addWidget(browse)
        apply_button = QtWidgets.QPushButton("Apply")
        apply_button.clicked.connect(self.apply)
        row.addWidget(apply_button)
        layout.addLayout(row)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        layout.addStretch(1)

    def browse(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "C++ tool directory", self.edit.text())
        if path:
            self.edit.setText(path)

    def apply(self):
        path = self.edit.text().strip()
        if path and not os.path.isdir(path):
            self.status.setStyleSheet(ERROR_STYLE)
            self.status.setText("{} is not a directory.".format(path))
            return False
        self.window.set_tools_dir(path)
        available = sum(1 for e in self.window.entries if e.available)
        self.status.setStyleSheet("" if available else ERROR_STYLE)
        self.status.setText("{} of {} tools available.".format(available, len(self.window.entries)))
        return True


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, settings, project_path=None, parent=None):
        super().__init__(parent)
        self.settings = settings
        self.project = None
        self.entries = []
        self.loader = None
        self.controller = RunController(self)
        self.setWindowTitle("HiCExplorer")

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabCloseRequested.connect(self._close_tab)
        self.setCentralWidget(self.tabs)

        self.run_view = RunView()
        self.run_view.attach(self.controller)
        self.editor = WorkflowEditor()
        self.editor.run_requested.connect(self.run_workflow)
        self.settings_page = SettingsPage(self)
        self.browser = None
        self._fixed = []
        for widget, title in ((self.run_view, "Runs"), (self.editor, "Workflow editor"),
                              (self.settings_page, "Settings")):
            self.tabs.addTab(widget, title)
            self._fixed.append(widget)
        self._add_browser()
        self.controller.finished.connect(lambda _code: self.refresh_project())

        self._build_docks()
        self._build_menus()
        self.statusBar().showMessage("Ready")
        self.set_tools_dir(settings.tools_dir)
        if project_path:
            self.open_project(project_path)

    # -- layout ---------------------------------------------------------
    def _add_browser(self):
        try:
            from .browser import MatrixBrowser
        except ImportError as exc:  # hicx_matrix or pyqtgraph missing
            page = QtWidgets.QLabel("The matrix browser needs pyqtgraph and the hicx_matrix module: {}".format(exc))
            page.setWordWrap(True)
            self.browser = None
        else:
            page = self.browser = MatrixBrowser()
        self.tabs.insertTab(2, page, "Matrix browser")
        self._fixed.append(page)

    def _build_docks(self):
        project_dock = QtWidgets.QDockWidget("Project", self)
        project_dock.setObjectName("project")
        self.project_tree = QtWidgets.QTreeWidget()
        self.project_tree.setHeaderHidden(True)
        self.project_tree.itemDoubleClicked.connect(self._project_item_opened)
        project_dock.setWidget(self.project_tree)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, project_dock)

        tools_dock = QtWidgets.QDockWidget("Tools", self)
        tools_dock.setObjectName("tools")
        panel = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(2, 2, 2, 2)
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
        tools_dock.setWidget(panel)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, tools_dock)
        self.splitDockWidget(project_dock, tools_dock, QtCore.Qt.Vertical)
        self.resizeDocks([project_dock, tools_dock], [1, 3], QtCore.Qt.Vertical)
        self.resizeDocks([project_dock], [260], QtCore.Qt.Horizontal)

    def _build_menus(self):
        menu = self.menuBar().addMenu("&File")
        new_action = QtGui.QAction("New project...", self)
        new_action.triggered.connect(self._new_project_dialog)
        menu.addAction(new_action)
        open_action = QtGui.QAction("Open project...", self)
        open_action.triggered.connect(self._open_project_dialog)
        menu.addAction(open_action)
        menu.addSeparator()
        quit_action = QtGui.QAction("Quit", self)
        quit_action.triggered.connect(self.close)
        menu.addAction(quit_action)
        view = self.menuBar().addMenu("&View")
        for title, method in (("Runs", self.show_runs), ("Workflow editor", self.show_editor),
                              ("Matrix browser", self.show_browser), ("Settings", self.show_settings)):
            action = QtGui.QAction(title, self)
            action.triggered.connect(method)
            view.addAction(action)

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

    def show_settings(self):
        self.tabs.setCurrentWidget(self.settings_page)

    def show_browser(self):
        self.tabs.setCurrentIndex(self.tabs.indexOf(self._fixed[-1]))

    # -- tools ----------------------------------------------------------
    def set_tools_dir(self, path):
        self.settings.tools_dir = path
        self.entries, self.loader = tool_entries(path)
        self._fill_tools()
        self.editor.set_context(self.project, self.entries, self.loader)
        available = sum(1 for e in self.entries if e.available)
        self.statusBar().showMessage("{} of {} tools available".format(available, len(self.entries)))

    def _fill_tools(self, *_):
        text = self.tool_filter.text().strip().lower()
        self.tool_tree.clear()
        available = QtWidgets.QTreeWidgetItem(["Available"])
        missing = QtWidgets.QTreeWidgetItem(["Not available"])
        for entry in self.entries:
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
        missing.setExpanded(True)
        # Wide enough for the longest tool name, so no name is elided.
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
            self.statusBar().showMessage("{}: {}".format(name, entry.reason if entry else "unknown tool"))
            return None
        page = ToolPage(self, entry.spec)
        if subcommand is not None:
            page.form.set_subcommand(subcommand)
        index = self.tabs.addTab(page, name)
        self.tabs.setCurrentIndex(index)
        return page

    # -- projects -------------------------------------------------------
    def _new_project_dialog(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "New project directory")
        if path:
            self.create_project(path)

    def _open_project_dialog(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Open project directory")
        if path:
            self.open_project(path)

    def create_project(self, path):
        return self._set_project(Project.create(path))

    def open_project(self, path):
        try:
            project = Project(path)
        except ProjectError as exc:
            self.statusBar().showMessage(str(exc))
            return None
        return self._set_project(project)

    def _set_project(self, project):
        self.project = project
        self.settings.add_recent_project(project.path)
        self.setWindowTitle("HiCExplorer - {}".format(project.name))
        self.run_view.set_project(project)
        self.editor.set_context(project, self.entries, self.loader)
        self.refresh_project()
        return project

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
        else:
            self.show_runs()

    def run_workflow(self, path, name):
        if self.project is None or self.controller.running:
            return False
        self.show_runs()
        return self.controller.run_workflow(self.project, self.loader, path, name)
