"""Analysis views (cpp/PLAN.md 10.6).

Each view shows the data outputs of one tool run, taken from a step record of
the run history: the tool, its logged command line and its working directory.
Where a result has a genomic position, the view emits ``navigate(matrix,
region)``, and the main window moves the matrix browser there.

Figure export reruns the logged command line with only the figure options
pointing at the export target, in a temporary directory that receives every
other output, so the exported figure is the figure the command-line tool
writes for those parameters, drawn by ``hicexplorer_plot``. A view whose tool
writes no figure (hicDifferentialTAD) exports its own plot and says so.

Analyses the port does not provide are shown, with the reason, rather than
left out (``UNAVAILABLE``).
"""

import glob
import json
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

UNAVAILABLE = {
    "hicrep": "HiCRep (stratum-adjusted correlation coefficient): not implemented, PLAN 9.10",
    "differential_loops": "Differential loops: the count-based differential engine is not "
                          "implemented yet, PLAN 9.7",
    "differential_compartments": "Differential compartments: the count-based differential "
                                 "engine is not implemented yet, PLAN 9.7",
    "chic_position": "The chicPlotViewpoint data holds distances to the reference point, not "
                     "genomic positions, so this view does not move the matrix browser",
}

# The dests of each tool's figure outputs.
FIGURE_DESTS = {
    "hicPlotDistVsCounts": ["plotFile"],
    "hicPlotViewpoint": ["outFileName"],
    "chicPlotViewpoint": ["outFileName"],
    "hicAggregateContacts": ["outFileName", "diagnosticHeatmapFile"],
    "hicCompartmentalization": ["outputFileName"],
    "hicCorrelate": ["outFileNameHeatmap", "outFileNameScatter"],
    "hicQC": ["outputFolder"],
    "hicPrepareQCreport": ["outputFolder"],
}
QC_BUILD_TOOLS = ("hicBuildMatrix", "hicBuildMatrixMicroC", "hicQuickQC")


class StepRecord:
    """One tool run: the tool, its argv (argv[0] the executable) and working directory."""

    def __init__(self, tool, argv, cwd, name=None):
        self.tool, self.argv, self.cwd = tool, list(argv), cwd
        self.name = name or tool


def records_from_history(history_dir):
    """The step records of a run history entry whose tools have a view."""
    with open(os.path.join(history_dir, "summary.json")) as handle:
        summary = json.load(handle)
    out = []
    for step in summary.get("steps") or []:
        if step.get("tool") in VIEW_CLASSES and step.get("exit_code") == 0 and step.get("argv"):
            out.append(StepRecord(step["tool"], step["argv"], summary.get("workdir") or history_dir,
                                  name="{} ({})".format(step.get("step"), step["tool"])))
    return out


class Argv:
    """Option values of a logged command line, found by any of the spec's flags."""

    def __init__(self, argv, spec):
        self.argv = list(argv)
        self.spec = spec
        sub = None
        if spec.has_subcommands and len(self.argv) > 1 and self.argv[1] in spec.commands:
            sub = self.argv[1]
        self.subcommand = sub
        self.by_dest = spec.arguments_by_dest(sub)
        self.flags = {flag for arg in self.by_dest.values() for flag in arg.get("flags") or []}

    def _find(self, dest):
        flags = set(self.by_dest.get(dest, {}).get("flags") or [])
        for index, token in enumerate(self.argv):
            if token in flags:
                end = index + 1
                while end < len(self.argv) and self.argv[end] not in self.flags:
                    end += 1
                return index, end
            if "=" in token and token.split("=", 1)[0] in flags:
                return index, index + 1
        return None

    def values(self, dest):
        found = self._find(dest)
        if found is None:
            return []
        index, end = found
        token = self.argv[index]
        if "=" in token and token.split("=", 1)[0] in self.flags and end == index + 1:
            return [token.split("=", 1)[1]]
        return self.argv[index + 1:end]

    def value(self, dest, default=None):
        values = self.values(dest)
        return values[0] if values else default

    def given(self, dest):
        return self._find(dest) is not None

    def replaced(self, replacements, removed=()):
        """A copy with dest -> [values] replaced (appended with the option's
        long flag when the command line does not have it) and the removed dests
        dropped."""
        from .workflow.spec import long_flag
        argv = list(self.argv)
        spans, appended = [], []
        for dest in list(replacements) + list(removed):
            found = self._find(dest)
            if found is not None:
                spans.append((found[0], found[1], dest))
            elif dest in replacements and dest in self.by_dest:
                appended += [long_flag(self.by_dest[dest])] + list(replacements[dest])
        for start, end, dest in sorted(spans, reverse=True):
            flag = self.argv[start].split("=", 1)[0]
            new = [flag] + list(replacements[dest]) if dest in replacements else []
            argv[start:end] = new
        return argv + appended


class _DecadeAxis(pg.AxisItem):
    """In log mode, labels the decades only: the minor ticks (2..9 per
    decade) keep their lines but no text, which otherwise runs together."""

    def tickStrings(self, values, scale, spacing):
        strings = super().tickStrings(values, scale, spacing)
        if not self.logMode:
            return strings
        return [s if abs(v - round(v)) < 1e-6 else "" for v, s in zip(values, strings)]


def _plot_widget(title=None, log_x=False, log_y=False):
    axes = {name: _DecadeAxis(name) for name, log in (("bottom", log_x), ("left", log_y)) if log}
    widget = pg.PlotWidget(title=title, axisItems=axes)
    widget.setBackground("w")
    widget.setMenuEnabled(False)
    widget.showGrid(x=True, y=True, alpha=0.2)
    widget.setLogMode(x=log_x, y=log_y)
    return widget


def _unavailable_label(key):
    label = QtWidgets.QLabel("Unavailable. " + UNAVAILABLE[key])
    label.setWordWrap(True)
    label.setEnabled(False)
    label.setObjectName("unavailable_" + key)
    return label


def _table(headers, rows):
    table = QtWidgets.QTableWidget(len(rows), len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    for r, row in enumerate(rows):
        for c, value in enumerate(row):
            table.setItem(r, c, QtWidgets.QTableWidgetItem(str(value)))
    table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Interactive)
    table.resizeColumnsToContents()
    for c in range(len(headers)):
        table.setColumnWidth(c, min(table.columnWidth(c), 260))
    return table


PALETTE = [(31, 119, 180), (214, 39, 40), (44, 160, 44), (148, 103, 189), (255, 127, 14),
           (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)]


class AnalysisView(QtWidgets.QWidget):
    """Base of the views: header with the run and the export button."""

    navigate = QtCore.Signal(str, str)   # matrix path ("" keeps the open matrix), region
    title = "Analysis"

    def __init__(self, record, loader, parent=None):
        super().__init__(parent)
        self.record = record
        self.loader = loader
        self.spec = loader.load(record.tool)
        self.args = Argv(record.argv, self.spec)
        self.messages = []
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        header = QtWidgets.QHBoxLayout()
        self.caption = QtWidgets.QLabel("<b>{}</b>: {}".format(self.title, record.name))
        self.caption.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.caption.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Preferred)
        header.addWidget(self.caption, 1)
        self.export_button = QtWidgets.QPushButton("Export figure...")
        self.export_button.clicked.connect(self._export_dialog)
        header.addWidget(self.export_button)
        outer.addLayout(header)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        outer.addWidget(self.status)
        self.body = QtWidgets.QVBoxLayout()
        outer.addLayout(self.body, 1)
        self.load()

    # -- paths and runs -------------------------------------------------
    def path(self, value):
        value = os.path.expanduser(value)
        return value if os.path.isabs(value) else os.path.normpath(os.path.join(self.record.cwd, value))

    def set_status(self, text):
        self.status.setText(text)

    def executable(self, tool=None):
        return self.loader.executable(tool or self.record.tool)

    def rerun(self, replacements, removed=(), workdir=None):
        """Runs the logged command with replaced options in workdir (a new
        temporary directory when None): every other output the tool writes is
        redirected there, inputs keep their absolute paths."""
        workdir = workdir or tempfile.mkdtemp(prefix="hicx-view-")
        argv = self.args.replaced(replacements, removed)
        redirected = Argv(argv, self.spec)
        extra = {}
        for dest, arg in redirected.by_dest.items():
            info = arg.get("file") or {}
            if dest in replacements or not redirected.given(dest):
                continue
            values = redirected.values(dest)
            if info.get("role") == "output":
                extra[dest] = [os.path.join(workdir, os.path.basename(v.rstrip("/")) or "out")
                               for v in values]
            elif info.get("role") == "input":
                extra[dest] = [self.path(v) for v in values]
        argv = Argv(argv, self.spec).replaced(extra)
        argv[0] = self.executable()
        proc = subprocess.run(argv, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
        return proc, workdir, argv

    def plot_data(self):
        """The figure data (C++-only --plotData) of the logged command."""
        tmp = tempfile.mkdtemp(prefix="hicx-view-")
        target = os.path.join(tmp, "plot_data.json")
        try:
            figure = {dest: [os.path.join(tmp, "figure_{}".format(dest))]
                      for dest in FIGURE_DESTS.get(self.record.tool, []) if self.args.given(dest)}
            figure["plotData"] = [target]
            proc, _, _ = self.rerun(figure, workdir=tmp)
            if proc.returncode != 0:
                raise RuntimeError("{} --plotData exited {}: {}".format(
                    self.record.tool, proc.returncode, proc.stderr.strip()[-400:]))
            with open(target) as handle:
                return json.load(handle)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    # -- export ---------------------------------------------------------
    def figure_targets(self, target):
        """dest -> export path: the first figure at target, further ones beside it."""
        dests = [d for d in FIGURE_DESTS.get(self.record.tool, []) if self.args.given(d)]
        stem, ext = os.path.splitext(target)
        out = {}
        for index, dest in enumerate(dests):
            out[dest] = target if index == 0 else "{}_{}{}".format(stem, dest, ext)
        return out

    def export_figure(self, target):
        """Writes the figure(s) the CLI tool writes for the logged parameters.
        Returns the written paths; raises RuntimeError with the tool's error."""
        targets = self.figure_targets(os.path.abspath(target))
        if not targets:
            raise RuntimeError("{} writes no figure for this run".format(self.record.tool))
        tmp = tempfile.mkdtemp(prefix="hicx-export-")
        try:
            proc, _, _ = self.rerun({d: [p] for d, p in targets.items()}, removed=("plotData",),
                                    workdir=tmp)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if proc.returncode != 0:
            raise RuntimeError("{} exited {}: {}".format(self.record.tool, proc.returncode,
                                                         proc.stderr.strip()[-600:]))
        return list(targets.values())

    def _export_dialog(self):
        folder = FIGURE_DESTS.get(self.record.tool) == ["outputFolder"] or self.record.tool in QC_BUILD_TOOLS
        if folder:
            target = QtWidgets.QFileDialog.getExistingDirectory(self, "Export the QC report into")
        else:
            target, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export figure", "",
                                                              "Images (*.png *.pdf *.svg)")
        if not target:
            return
        try:
            written = self.export_figure(target)
        except RuntimeError as exc:
            self.set_status("Export failed: {}".format(exc))
            return
        self.set_status("Exported: " + ", ".join(written))

    def load(self):
        raise NotImplementedError


# ---------------------------------------------------------------------------
# QC report


def read_tsv(path):
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    return (lines[0], lines[1:]) if lines else ([], [])


class QCReportView(AnalysisView):
    title = "QC report"

    def folder(self):
        if self.record.tool in QC_BUILD_TOOLS:
            return self.path(self.args.value("QCfolder", "."))
        return self.path(self.args.value("outputFolder", "."))

    def load(self):
        folder = self.folder()
        self.tabs = QtWidgets.QTabWidget()
        self.body.addWidget(self.tabs)
        self.tables = {}
        for path in sorted(glob.glob(os.path.join(folder, "*_table.txt"))):
            headers, rows = read_tsv(path)
            name = os.path.basename(path)[:-len("_table.txt")]
            page = QtWidgets.QSplitter(QtCore.Qt.Vertical)
            table = _table(headers, rows)
            page.addWidget(table)
            page.addWidget(self._chart(headers, rows))
            self.tabs.addTab(page, name.replace("_", " "))
            self.tables[name] = (headers, rows)
        if not self.tables:
            self.set_status("No QC tables in {}".format(folder))

    @staticmethod
    def _chart(headers, rows):
        plot = _plot_widget()
        columns = [i for i, h in enumerate(headers) if i > 0 and not h.rstrip().endswith("%")
                   and all(_is_number(row[i]) for row in rows if i < len(row))]
        # Horizontal bars: the column names are long, and on the y-axis each
        # gets its own line instead of overlapping under the bars.
        thickness = 0.8 / max(1, len(rows))
        for r, row in enumerate(rows):
            y = np.arange(len(columns)) + r * thickness
            values = [float(row[i]) for i in columns]
            plot.addItem(pg.BarGraphItem(x0=0, y=y, width=values, height=thickness * 0.9,
                                         brush=pg.mkBrush(*PALETTE[r % len(PALETTE)]),
                                         name=row[0] if row else None))
        plot.getAxis("left").setTicks([[(i + 0.4 - thickness / 2, headers[c]) for i, c in enumerate(columns)]])
        plot.invertY(True)
        return plot

    def export_figure(self, target):
        if self.record.tool not in QC_BUILD_TOOLS:
            return super().export_figure(target)
        # The matrix builders draw their QC report as hicQC does from QC.log.
        argv = [self.executable("hicQC"), "--logfiles", os.path.join(self.folder(), "QC.log"),
                "--outputFolder", os.path.abspath(target)]
        proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if proc.returncode != 0:
            raise RuntimeError("hicQC exited {}: {}".format(proc.returncode, proc.stderr.strip()[-600:]))
        return [os.path.abspath(target)]


def _is_number(text):
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


# ---------------------------------------------------------------------------
# distance decay


def read_distance_table(path):
    """(label, chromosome) -> (distances, contacts) of a --outFileData table.
    The tool appends one table per curve, each with its header line."""
    headers, rows = read_tsv(path)
    index = {name: i for i, name in enumerate(headers)}
    curves = {}
    for row in rows:
        if len(row) <= index["Contacts"] or not _is_number(row[index["Distance"]]):
            continue
        key = (row[index["Matrix"]], row[index["Chromosome"]])
        distances, contacts = curves.setdefault(key, ([], []))
        distances.append(float(row[index["Distance"]]))
        contacts.append(float(row[index["Contacts"]]))
    return curves


class DistanceDecayView(AnalysisView):
    title = "Distance decay"

    def load(self):
        self.plot = _plot_widget(log_x=True, log_y=True)
        self.plot.setLabel("bottom", "genomic distance (bp)")
        self.plot.setLabel("left", "mean contacts")
        self.legend = self.plot.addLegend(offset=(-10, 10))
        row = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("Overlay another distance table...")
        add.clicked.connect(self._add_dialog)
        row.addWidget(add)
        row.addStretch(1)
        self.body.addLayout(row)
        self.body.addWidget(self.plot, 1)
        self.curves = {}
        data = self.args.value("outFileData")
        if data and os.path.isfile(self.path(data)):
            self.add_table(self.path(data))
        else:
            payload = self.plot_data()
            for matrix in payload.get("matrices") or []:
                for series in matrix.get("series") or []:
                    self._add_curve((matrix.get("label"), series.get("chrom")), series["x"], series["y"])

    def add_table(self, path):
        for key, (x, y) in read_distance_table(path).items():
            self._add_curve(key, x, y)

    def _add_curve(self, key, x, y):
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        keep = (x > 0) & (y > 0)
        pen = pg.mkPen(PALETTE[len(self.curves) % len(PALETTE)], width=2)
        name = "{} {}".format(*key) if key[1] not in (None, "all") else str(key[0])
        self.curves[key] = self.plot.plot(x[keep], y[keep], pen=pen, name=name)

    def _add_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Distance table (--outFileData)")
        if path:
            self.add_table(path)


# ---------------------------------------------------------------------------
# viewpoints


class ViewpointView(AnalysisView):
    title = "Viewpoint"

    def load(self):
        self.plot = _plot_widget()
        self.legend = self.plot.addLegend(offset=(-10, 10))
        self.curves = {}
        self.chrom = None
        self.matrix = None
        if self.record.tool == "hicPlotViewpoint":
            self._load_hic()
        else:
            self._load_chic()
        self.body.addWidget(self.plot, 1)

    def _load_hic(self):
        self.matrix = self.path(self.args.value("matrix", ""))
        reference = self.args.values("referencePoint")
        self.chrom = reference[0].split(":")[0] if reference else None
        self.plot.setLabel("bottom", "position on {} (bp); click to show it in the matrix browser".format(self.chrom))
        self.plot.setLabel("left", "contacts")
        prefix = self.args.value("interactionOutFileName")
        files = sorted(glob.glob(self.path(prefix) + "_*.bedgraph")) if prefix else []
        if files:
            for path in files:
                x, y = [], []
                with open(path) as handle:
                    for line in handle:
                        fields = line.split()
                        if len(fields) >= 7:
                            x.append((float(fields[4]) + float(fields[5])) / 2.0)
                            y.append(float(fields[6]))
                label = os.path.basename(path)[len(os.path.basename(self.path(prefix))) + 1:-len(".bedgraph")]
                self._curve(label, x, y)
        else:
            payload = self.plot_data()
            start, end = payload["region_start"], payload["region_end"]
            for label, values in zip(payload.get("legend") or [], payload.get("data") or []):
                x = np.linspace(start, end, len(values), endpoint=False)
                self._curve(label, x, values)
        self.plot.scene().sigMouseClicked.connect(self._clicked)

    def _load_chic(self):
        payload = self.plot_data()
        upstream, downstream = payload["range"]
        resolution = payload.get("resolution") or 1
        self.plot.setLabel("bottom", "distance to the reference point (bp)")
        self.plot.setLabel("left", "normalised contacts")
        self.chic_items = [item for chunk in payload.get("chunks") or [] for group in chunk
                           for item in group.get("items") or [] if not item.get("skip")]
        self.selector = QtWidgets.QComboBox()
        self.selector.addItems([item["label"] for item in self.chic_items])
        self.selector.currentIndexChanged.connect(self._show_chic)
        self.body.addWidget(self.selector)
        self.body.addWidget(_unavailable_label("chic_position"))
        self._chic_x = lambda n: np.linspace(-upstream, downstream, n)
        self._show_chic(0)

    def _show_chic(self, index):
        self.plot.clear()
        self.legend.clear()
        self.curves = {}
        if not self.chic_items:
            return
        item = self.chic_items[index]
        data = np.asarray(item["data"], dtype=float)
        self._curve(item["label"], self._chic_x(len(data)), data)
        if item.get("background"):
            background = np.asarray(item["background"], dtype=float)
            self._curve("background model", self._chic_x(len(background)), background)

    def _curve(self, label, x, y):
        pen = pg.mkPen(PALETTE[len(self.curves) % len(PALETTE)], width=2)
        self.curves[label] = self.plot.plot(np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                            pen=pen, name=label)

    def position_region(self, x, half_width=250000):
        start = max(0, int(x - half_width))
        return "{}:{}-{}".format(self.chrom, start, int(x + half_width))

    def _clicked(self, event):
        if self.chrom is None:
            return
        point = self.plot.plotItem.vb.mapSceneToView(event.scenePos())
        self.navigate.emit(self.matrix or "", self.position_region(point.x()))


# ---------------------------------------------------------------------------
# aggregate contacts and saddle plots


def _colormap(name):
    if name == "blue-white-red":
        return pg.ColorMap([0.0, 0.5, 1.0], [(30, 60, 200), (255, 255, 255), (200, 20, 20)])
    return pg.colormap.get(name)


def _heatmap(values, title, colormap="viridis"):
    plot = _plot_widget(title)
    plot.invertY(True)
    plot.setAspectLocked(True)
    image = pg.ImageItem(np.asarray(values, dtype=float), axisOrder="row-major")
    image.setColorMap(_colormap(colormap))
    plot.addItem(image)
    return plot, image


class AggregateContactsView(AnalysisView):
    title = "Aggregate contacts"

    def load(self):
        self.matrix = self.path(self.args.value("matrix", ""))
        prefix = self.args.value("outFilePrefixMatrix")
        self.matrix_files = sorted(glob.glob(self.path(prefix) + "_*.tab")) if prefix else []
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        self.selector = QtWidgets.QComboBox()
        self.selector.addItems([os.path.basename(p) for p in self.matrix_files])
        self.selector.currentIndexChanged.connect(self._show_matrix)
        left_layout.addWidget(self.selector)
        self.heatmap, self.image = _heatmap(np.zeros((1, 1)), "aggregate", "magma")
        left_layout.addWidget(self.heatmap, 1)
        splitter.addWidget(left)
        pairs_prefix = self.args.value("outFileContactPairs")
        self.pairs = []
        for path in sorted(glob.glob(self.path(pairs_prefix) + "_*.tab")) if pairs_prefix else []:
            with open(path) as handle:
                for line in handle:
                    fields = line.split("\t")
                    if len(fields) >= 7:
                        self.pairs.append([os.path.basename(path)] + [f.strip() for f in fields[:7]])
        self.pair_table = _table(["file", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "value"],
                                 self.pairs)
        self.pair_table.itemSelectionChanged.connect(self._pair_selected)
        right = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right)
        right_layout.addWidget(QtWidgets.QLabel("Contact pairs: select one to show it in the matrix browser"))
        right_layout.addWidget(self.pair_table, 1)
        splitter.addWidget(right)
        self.body.addWidget(splitter, 1)
        if not self.matrix_files:
            self.set_status("The run wrote no --outFilePrefixMatrix tables.")
        self._show_matrix(0)

    def _show_matrix(self, index):
        if 0 <= index < len(self.matrix_files):
            self.values = np.loadtxt(self.matrix_files[index], delimiter="\t", ndmin=2)
            self.image.setImage(self.values, autoLevels=True)

    def pair_region(self, row):
        _, c1, s1, e1, c2, s2, e2, _ = self.pairs[row]
        if c1 != c2:
            return c1, "{}:{}-{}".format(c1, s1, e1)
        lo, hi = min(int(s1), int(s2)), max(int(e1), int(e2))
        pad = max(50000, (hi - lo) // 2)
        return c1, "{}:{}-{}".format(c1, max(0, lo - pad), hi + pad)

    def _pair_selected(self):
        rows = {index.row() for index in self.pair_table.selectedIndexes()}
        if rows:
            self.navigate.emit(self.matrix, self.pair_region(min(rows))[1])


class SaddleView(AnalysisView):
    title = "Compartment saddle"

    def load(self):
        matrix = self.args.value("outputMatrix")
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.saddles = np.zeros((0, 1, 1))
        if matrix:
            path = self.path(matrix)
            path = path if path.endswith(".npz") else path + ".npz"
            if os.path.isfile(path):
                with np.load(path) as archive:
                    self.saddles = np.asarray(archive["arr_0"], dtype=float)
        left = QtWidgets.QWidget()
        left_layout = QtWidgets.QVBoxLayout(left)
        self.selector = QtWidgets.QComboBox()
        labels = [os.path.basename(p) for p in self.args.values("obsexp_matrices")]
        self.selector.addItems(labels[:len(self.saddles)])
        self.selector.currentIndexChanged.connect(self._show)
        left_layout.addWidget(self.selector)
        self.heatmap, self.image = _heatmap(np.zeros((1, 1)), "log2 normalised interactions per quantile pair",
                                            "blue-white-red")
        left_layout.addWidget(self.heatmap, 1)
        splitter.addWidget(left)
        self.ratio_plot = _plot_widget("polarization ratio")
        self.ratio_plot.setLabel("bottom", "quantile")
        dat = self.path(self.args.value("outputFileName", "")) + "_dat"
        self.ratios = np.loadtxt(dat, ndmin=2) if os.path.isfile(dat) else np.zeros((0, 0))
        for index, row in enumerate(self.ratios):
            finite = np.isfinite(row)
            self.ratio_plot.plot(np.arange(len(row))[finite], row[finite],
                                 pen=pg.mkPen(PALETTE[index % len(PALETTE)], width=2), symbol="o",
                                 symbolSize=5)
        splitter.addWidget(self.ratio_plot)
        self.body.addWidget(splitter, 1)
        if not len(self.saddles):
            self.set_status("The run wrote no --outputMatrix archive.")
        self._show(0)

    def _show(self, index):
        if 0 <= index < len(self.saddles):
            values = self.saddles[index]
            with np.errstate(divide="ignore", invalid="ignore"):
                shown = np.log2(values)
            shown[~np.isfinite(shown)] = 0.0
            self.image.setImage(shown, autoLevels=True)


# ---------------------------------------------------------------------------
# correlation


class CorrelationView(AnalysisView):
    title = "Correlation"

    def load(self):
        payload = self.plot_data()
        self.labels = [str(label).strip("'") for label in payload.get("labels") or []]
        self.results = np.asarray(payload.get("results") or [[]], dtype=float)
        if self.results.ndim == 2 and self.results.shape[0] == self.results.shape[1] > 1:
            # The tool fills the upper triangle only (its figure masks the
            # lower one); a zero there is no correlation, so mirror it.
            lower = np.tril_indices(self.results.shape[0], -1)
            if not self.results[lower].any():
                self.results[lower] = self.results.T[lower]
        info = QtWidgets.QLabel("{} correlation{}".format(payload.get("method", ""),
                                                         " of log1p values" if payload.get("log1p") else ""))
        self.body.addWidget(info)
        self.body.addWidget(_unavailable_label("hicrep"))
        plot, self.image = _heatmap(self.results, "correlation", "viridis")
        n = len(self.labels)
        ticks = [[(i + 0.5, self.labels[i]) for i in range(n)]]
        plot.getAxis("bottom").setTicks(ticks)
        plot.getAxis("left").setTicks(ticks)
        for i in range(n):
            for j in range(n):
                text = pg.TextItem("{:.3f}".format(self.results[i, j]), color="w", anchor=(0.5, 0.5))
                text.setPos(j + 0.5, i + 0.5)
                plot.addItem(text)
        self.body.addWidget(plot, 1)


# ---------------------------------------------------------------------------
# differential TADs


def read_diff_tad(path):
    """(columns, rows) of a hicDifferentialTAD accepted or rejected file."""
    columns, rows = [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith("# Chromosome"):
                columns = line[2:].rstrip("\n").split("\t")
            elif line.startswith("#") or not line.strip():
                continue
            else:
                rows.append(line.rstrip("\n").split("\t"))
    return columns, rows


class DifferentialView(AnalysisView):
    title = "Differential TADs"
    TESTS = ("intra-TAD", "left-inter-TAD", "right-inter-TAD")

    def load(self):
        self.matrix = self.path(self.args.value("targetMatrix", ""))
        prefix = self.path(self.args.value("outFileNamePrefix", ""))
        self.columns, self.rows = [], []
        for status in ("rejected", "accepted"):
            path = "{}_{}.diff_tad".format(prefix, status)
            if os.path.isfile(path):
                columns, rows = read_diff_tad(path)
                self.columns = ["status"] + columns
                self.rows.extend([[status] + row for row in rows])
        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(QtWidgets.QLabel("Test:"))
        self.test = QtWidgets.QComboBox()
        self.test.addItems(self.TESTS)
        self.test.currentIndexChanged.connect(self._draw)
        controls.addWidget(self.test)
        self.adjusted = any(c.startswith("adjusted p-value") for c in self.columns)
        controls.addWidget(QtWidgets.QLabel("y: -log10 {}p-value".format("adjusted " if self.adjusted else "")))
        controls.addStretch(1)
        self.body.addLayout(controls)
        for key in ("differential_loops", "differential_compartments"):
            self.body.addWidget(_unavailable_label(key))
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.volcano = _plot_widget("volcano")
        self.volcano.setLabel("bottom", "effect: rank-sum statistic W")
        self.volcano.setLabel("left", "-log10 p")
        self.scatter = pg.ScatterPlotItem(size=7, pen=None)
        self.scatter.sigClicked.connect(self._point_clicked)
        self.volcano.addItem(self.scatter)
        splitter.addWidget(self.volcano)
        self.table = _table(self.columns, self.rows)
        self.table.itemSelectionChanged.connect(self._row_selected)
        splitter.addWidget(self.table)
        splitter.setSizes([500, 700])
        self.body.addWidget(splitter, 1)
        if not self.rows:
            self.set_status("No accepted or rejected file with prefix {}".format(prefix))
        self._draw()

    def volcano_points(self, test=None):
        test = test or self.test.currentText()
        index = {name: i for i, name in enumerate(self.columns)}
        w = index.get("W " + test)
        p = index.get(("adjusted p-value " if self.adjusted else "p-value ") + test)
        points = []
        for row_index, row in enumerate(self.rows):
            try:
                effect, pvalue = float(row[w]), float(row[p])
            except (TypeError, ValueError, IndexError):
                continue
            if not (math.isfinite(effect) and math.isfinite(pvalue)) or pvalue <= 0:
                continue
            points.append((row_index, effect, -math.log10(pvalue)))
        return points

    def _draw(self, *_):
        spots = [{"pos": (effect, y), "data": row,
                  "brush": pg.mkBrush(214, 39, 40, 200) if self.rows[row][0] == "rejected"
                  else pg.mkBrush(90, 90, 90, 150)} for row, effect, y in self.volcano_points()]
        self.scatter.setData(spots)

    def region(self, row):
        index = {name: i for i, name in enumerate(self.columns)}
        values = self.rows[row]
        chrom = values[index["Chromosome"]]
        start, end = int(values[index["start"]]), int(values[index["end"]])
        pad = (end - start) // 2
        return "{}:{}-{}".format(chrom, max(0, start - pad), end + pad)

    def select_tad(self, row):
        self.table.selectRow(row)
        self.navigate.emit(self.matrix, self.region(row))

    def _row_selected(self):
        rows = {index.row() for index in self.table.selectedIndexes()}
        if rows:
            self.navigate.emit(self.matrix, self.region(min(rows)))

    def _point_clicked(self, _item, points, *_):
        if points:
            self.select_tad(points[0].data())

    def figure_targets(self, target):
        return {}

    def export_figure(self, target):
        """hicDifferentialTAD writes no figure; the volcano plot is exported as shown."""
        from pyqtgraph.exporters import ImageExporter
        exporter = ImageExporter(self.volcano.plotItem)
        exporter.export(os.path.abspath(target))
        return [os.path.abspath(target)]


VIEW_CLASSES = {
    "hicQC": QCReportView, "hicPrepareQCreport": QCReportView,
    "hicBuildMatrix": QCReportView, "hicBuildMatrixMicroC": QCReportView, "hicQuickQC": QCReportView,
    "hicPlotDistVsCounts": DistanceDecayView,
    "hicPlotViewpoint": ViewpointView, "chicPlotViewpoint": ViewpointView,
    "hicAggregateContacts": AggregateContactsView,
    "hicCompartmentalization": SaddleView,
    "hicCorrelate": CorrelationView,
    "hicDifferentialTAD": DifferentialView,
}


def open_view(record, loader, parent=None):
    return VIEW_CLASSES[record.tool](record, loader, parent)
