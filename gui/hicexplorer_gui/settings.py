"""User settings, kept by QSettings (never in the repository)."""

import os

from PySide6 import QtCore


class Settings:
    """The C++ tool directory, the shared minibwa index cache directory, and
    recent projects.

    ``HICX_CPP_BIN`` fills in the tool directory when none is stored, so a
    session can be pointed at a build without touching the stored settings.
    ``HICX_INDEX_CACHE`` does the same for the index cache directory
    (PLAN.md tier 13: a minibwa index built with hicBuildIndex can be saved
    per project, using the project directory the same way every other tool
    output can, or in this one shared directory so several projects reuse
    the same custom reference index instead of rebuilding it).
    """

    def __init__(self, qsettings=None):
        self._s = qsettings if qsettings is not None else QtCore.QSettings("HiCExplorer", "hicexplorer-gui")

    @property
    def tools_dir(self):
        value = self._s.value("tools_dir", "", type=str)
        return value or os.environ.get("HICX_CPP_BIN", "")

    @tools_dir.setter
    def tools_dir(self, value):
        self._s.setValue("tools_dir", value or "")
        self._s.sync()

    @property
    def index_cache_dir(self):
        value = self._s.value("index_cache_dir", "", type=str)
        return value or os.environ.get("HICX_INDEX_CACHE", "")

    @index_cache_dir.setter
    def index_cache_dir(self, value):
        self._s.setValue("index_cache_dir", value or "")
        self._s.sync()

    @property
    def recent_projects(self):
        value = self._s.value("recent_projects", [], type=list)
        return [p for p in value if isinstance(p, str)]

    def add_recent_project(self, path):
        items = [path] + [p for p in self.recent_projects if p != path]
        self._s.setValue("recent_projects", items[:10])
        self._s.sync()
