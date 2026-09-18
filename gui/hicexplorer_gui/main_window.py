"""The main window: a tab per open project, and the window-global C++ tool
directory (Settings, moved into the File menu, PLAN 10 GUI redesign
2026-09-18).

Each open project is its own ProjectTab (project_tab.py): its own data,
filtered tool browser, runs, workflow editor, matrix browser and dynamically
opened tool-form and analysis-view tabs. The tool directory is the one thing
every open project shares, since there is a single build of the C++ tools;
the resulting catalog (entries, loader) is pushed into every open ProjectTab
whenever it changes.
"""

import os

from PySide6 import QtCore, QtGui, QtWidgets

from .catalog import tool_entries
from .project import Project, ProjectError
from .project_tab import ProjectTab


class SettingsPage(QtWidgets.QWidget):
    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(QtWidgets.QLabel("<b>C++ tool directory</b>"))
        hint = QtWidgets.QLabel("The directory holding the HiCExplorer C++ tool executables "
                                "(for example the tools directory of a build). Shared by every "
                                "open project, since there is one build of the tools.")
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

        layout.addWidget(QtWidgets.QLabel("<b>Shared minibwa index cache</b>"))
        index_hint = QtWidgets.QLabel(
            "Where hicBuildIndex saves a custom reference index by default (PLAN tier 13), "
            "so several projects can reuse it instead of rebuilding it. A tool form's own file "
            "picker can still point anywhere, including inside a single project; this is only "
            "the starting suggestion. Leave empty to default to the current project's directory.")
        index_hint.setWordWrap(True)
        layout.addWidget(index_hint)
        index_row = QtWidgets.QHBoxLayout()
        self.index_edit = QtWidgets.QLineEdit(window.settings.index_cache_dir)
        index_row.addWidget(self.index_edit, 1)
        index_browse = QtWidgets.QPushButton("Browse")
        index_browse.clicked.connect(self.browse_index_cache)
        index_row.addWidget(index_browse)
        index_apply = QtWidgets.QPushButton("Apply")
        index_apply.clicked.connect(self.apply_index_cache)
        index_row.addWidget(index_apply)
        layout.addLayout(index_row)
        self.index_status = QtWidgets.QLabel()
        self.index_status.setWordWrap(True)
        layout.addWidget(self.index_status)
        layout.addStretch(1)

    def browse(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "C++ tool directory", self.edit.text())
        if path:
            self.edit.setText(path)

    def apply(self):
        path = self.edit.text().strip()
        if path and not os.path.isdir(path):
            self.status.setStyleSheet("color: #d0314b;")
            self.status.setText("{} is not a directory.".format(path))
            return False
        self.window.set_tools_dir(path)
        available = sum(1 for e in self.window.entries if e.available)
        self.status.setStyleSheet("" if available else "color: #d0314b;")
        self.status.setText("{} of {} tools available.".format(available, len(self.window.entries)))
        return True

    def browse_index_cache(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Shared minibwa index cache",
                                                           self.index_edit.text())
        if path:
            self.index_edit.setText(path)

    def apply_index_cache(self):
        path = self.index_edit.text().strip()
        if path and not os.path.isdir(path):
            self.index_status.setStyleSheet("color: #d0314b;")
            self.index_status.setText("{} is not a directory.".format(path))
            return False
        self.window.settings.index_cache_dir = path
        self.index_status.setStyleSheet("")
        self.index_status.setText("Saved." if path else "Cleared; defaults to the project directory.")
        return True


class SettingsDialog(QtWidgets.QDialog):
    """The C++ tool directory, as a non-modal dialog reached from File >
    Settings, rather than a permanently visible tab."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.page = SettingsPage(window)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.page)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)
        self.resize(520, 220)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, settings, project_path=None, parent=None):
        super().__init__(parent)
        self.settings = settings
        self.entries = []
        self.loader = None
        self.setWindowTitle("HiCExplorer")

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabCloseRequested.connect(self._close_project_tab)
        self.setCentralWidget(self.tabs)

        self._settings_dialog = None
        self._build_menus()
        self.statusBar().showMessage("Ready")
        self.set_tools_dir(settings.tools_dir)
        if project_path:
            self.open_project(project_path)

    # -- menus --------------------------------------------------------------
    def _build_menus(self):
        menu = self.menuBar().addMenu("&File")
        new_action = QtGui.QAction("New project...", self)
        new_action.triggered.connect(self._new_project_dialog)
        menu.addAction(new_action)
        open_action = QtGui.QAction("Open project...", self)
        open_action.triggered.connect(self._open_project_dialog)
        menu.addAction(open_action)
        self.recent_menu = menu.addMenu("Open recent")
        self._fill_recent_menu()
        menu.addSeparator()
        template_action = QtGui.QAction("New workflow from template...", self)
        template_action.triggered.connect(self.new_from_template)
        menu.addAction(template_action)
        menu.addSeparator()
        settings_action = QtGui.QAction("Settings...", self)
        settings_action.triggered.connect(self.open_settings_dialog)
        menu.addAction(settings_action)
        menu.addSeparator()
        quit_action = QtGui.QAction("Quit", self)
        quit_action.triggered.connect(self.close)
        menu.addAction(quit_action)

        analysis = self.menuBar().addMenu("&Analysis")
        views_action = QtGui.QAction("Open the analysis views of a run...", self)
        views_action.triggered.connect(self._open_views_dialog)
        analysis.addAction(views_action)

        view = self.menuBar().addMenu("&View")
        for title, method in (("Runs", self.show_runs), ("Workflow editor", self.show_editor),
                              ("Matrix browser", self.show_browser)):
            action = QtGui.QAction(title, self)
            action.triggered.connect(method)
            view.addAction(action)
        view.addSeparator()
        settings_view_action = QtGui.QAction("Settings...", self)
        settings_view_action.triggered.connect(self.open_settings_dialog)
        view.addAction(settings_view_action)

    def _fill_recent_menu(self):
        self.recent_menu.clear()
        recent = self.settings.recent_projects
        if not recent:
            empty = QtGui.QAction("(no recent projects)", self)
            empty.setEnabled(False)
            self.recent_menu.addAction(empty)
            return
        for path in recent:
            action = QtGui.QAction(path, self)
            action.triggered.connect(lambda checked=False, p=path: self.open_project(p))
            self.recent_menu.addAction(action)

    def open_settings_dialog(self):
        if self._settings_dialog is None:
            self._settings_dialog = SettingsDialog(self, self)
        else:
            self._settings_dialog.page.edit.setText(self.settings.tools_dir)
            self._settings_dialog.page.index_edit.setText(self.settings.index_cache_dir)
        self._settings_dialog.show()
        self._settings_dialog.raise_()
        self._settings_dialog.activateWindow()
        return self._settings_dialog

    # -- the window-global tool catalog --------------------------------------
    def set_tools_dir(self, path):
        self.settings.tools_dir = path
        self.entries, self.loader = tool_entries(path)
        for index in range(self.tabs.count()):
            self.tabs.widget(index).set_catalog(self.entries, self.loader)
        available = sum(1 for e in self.entries if e.available)
        self.statusBar().showMessage("{} of {} tools available".format(available, len(self.entries)))

    # -- project tabs ---------------------------------------------------------
    def _new_project_dialog(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "New project directory")
        if path:
            self.create_project(path)

    def _open_project_dialog(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Open project directory")
        if path:
            self.open_project(path)

    def create_project(self, path):
        return self._add_project_tab(Project.create(path))

    def open_project(self, path):
        existing = self._find_open_tab(path)
        if existing is not None:
            self.tabs.setCurrentWidget(existing)
            return existing.project
        try:
            project = Project(path)
        except ProjectError as exc:
            self.statusBar().showMessage(str(exc))
            return None
        return self._add_project_tab(project)

    def _find_open_tab(self, path):
        target = os.path.abspath(path)
        for index in range(self.tabs.count()):
            tab = self.tabs.widget(index)
            if os.path.abspath(tab.project.path) == target:
                return tab
        return None

    def _add_project_tab(self, project):
        tab = ProjectTab(self, project, self.entries, self.loader)
        index = self.tabs.addTab(tab, project.name)
        self.tabs.setCurrentIndex(index)
        self.settings.add_recent_project(project.path)
        self._fill_recent_menu()
        return tab.project

    def _close_project_tab(self, index):
        tab = self.tabs.widget(index)
        if tab.controller.running:
            answer = QtWidgets.QMessageBox.question(
                self, "Close project",
                "A run is in progress in project \"{}\". Close it anyway?".format(tab.project.name),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if answer != QtWidgets.QMessageBox.Yes:
                return
            tab.controller.cancel()
        self.tabs.removeTab(index)
        tab.deleteLater()

    def active_project_tab(self):
        return self.tabs.currentWidget()

    # -- view menu, acting on the active project -----------------------------
    def show_runs(self):
        tab = self.active_project_tab()
        if tab is not None:
            tab.show_runs()

    def show_editor(self):
        tab = self.active_project_tab()
        if tab is not None:
            tab.show_editor()

    def show_browser(self):
        tab = self.active_project_tab()
        if tab is not None:
            tab.show_browser()

    def new_from_template(self):
        tab = self.active_project_tab()
        if tab is None:
            self.statusBar().showMessage("Open or create a project first.")
            return None
        return tab.new_from_template()

    def _open_views_dialog(self):
        tab = self.active_project_tab()
        if tab is None:
            self.statusBar().showMessage("Open or create a project first.")
            return None
        start = tab.project.history_dir if tab.project is not None else ""
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Run history entry", start)
        if path:
            tab.open_analysis_views(path)
