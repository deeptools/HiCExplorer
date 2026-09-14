"""Matrix browser (PLAN 10.5): cool, mcool, h5 and .hic through hicx_matrix.

Only the visible region is read and held. Navigation (the region box, zoom
and pan) schedules a fetch of the bins the view shows; the fetch runs on a
worker thread, a newer request replaces a waiting one, and the result of a
request that has been superseded is dropped, so the interface never blocks
and never shows a stale region. For mcool and .hic files the resolution is
chosen from the span on screen: the finest one that keeps the view within
MAX_BINS bins. A single-resolution file (cool, h5) whose view would exceed
MAX_BINS is not fetched; the browser asks to zoom in instead.
"""

import math
import os
import threading

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

import hicx_matrix

MAX_BINS = 1200
MODES = ("single", "side by side", "difference")
COLORMAPS = ("white-red", "viridis", "inferno", "magma", "cividis")


def colormap(name):
    if name == "white-red":
        return pg.ColorMap([0.0, 0.5, 1.0], [(255, 255, 255), (250, 130, 90), (140, 0, 0)])
    if name == "blue-white-red":
        return pg.ColorMap([0.0, 0.5, 1.0], [(30, 60, 200), (255, 255, 255), (200, 20, 20)])
    return pg.colormap.get(name)


def parse_region_text(text, lengths):
    """(chrom, start, end) of "chrom" or "chrom:start-end" (commas allowed)."""
    text = text.strip().replace(",", "")
    if ":" in text:
        chrom, _, span = text.rpartition(":")
        if "-" not in span:
            raise ValueError("expected chrom:start-end, found {!r}".format(text))
        start_text, end_text = span.split("-", 1)
        try:
            start, end = int(float(start_text)), int(float(end_text))
        except ValueError:
            raise ValueError("start and end of {!r} are not numbers".format(text))
    else:
        chrom, start, end = text, 0, None
    if chrom not in lengths:
        raise ValueError("unknown chromosome {!r}".format(chrom))
    end = lengths[chrom] if end is None else end
    start, end = max(0, start), min(lengths[chrom], end)
    if end <= start:
        raise ValueError("the region {!r} is empty".format(text))
    return chrom, start, end


def choose_resolution(resolutions, span, max_bins=MAX_BINS):
    """The finest resolution that shows ``span`` bp in at most max_bins bins."""
    if not resolutions:
        return None
    for resolution in sorted(resolutions):
        if span / resolution <= max_bins:
            return resolution
    return max(resolutions)


class MatrixSource:
    def __init__(self, path):
        self.path = path
        self.file = hicx_matrix.open(path)
        self.format = self.file.format
        self.chromosomes = self.file.chromosomes()
        self.lengths = dict(self.chromosomes)
        self.resolutions = list(self.file.resolutions())
        self.normalizations = list(self.file.normalizations())

    @property
    def name(self):
        return os.path.basename(self.path)


class Request:
    def __init__(self, generation, fetches, chrom, rows, cols):
        self.generation = generation
        self.fetches = fetches    # [(source, region1, region2, resolution, normalization)]
        self.chrom = chrom
        self.rows = rows          # (start, end) bp of the fetched rows, bin aligned
        self.cols = cols


class Fetcher(QtCore.QObject):
    """One fetch at a time off the UI thread; newer requests supersede older ones."""

    done = QtCore.Signal(object, object, str)   # request, [arrays], error
    _finished = QtCore.Signal(object, object, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.generation = 0
        self.running = None
        self.pending = None
        self.stale = 0
        self._finished.connect(self._on_finished)

    def submit(self, fetches, chrom, rows, cols):
        self.generation += 1
        request = Request(self.generation, fetches, chrom, rows, cols)
        if self.running is None:
            self._start(request)
        else:
            if self.pending is not None:
                self.stale += 1
            self.pending = request
        return request

    @property
    def busy(self):
        return self.running is not None or self.pending is not None

    def _start(self, request):
        self.running = request
        threading.Thread(target=self._work, args=(request,), name="hicx-fetch", daemon=True).start()

    def _work(self, request):
        arrays, error = [], ""
        try:
            for source, region1, region2, resolution, normalization in request.fetches:
                arrays.append(source.file.fetch(region1, region2, resolution=resolution,
                                                normalization=normalization))
        except Exception as exc:  # noqa: BLE001 - shown inline in the browser
            arrays, error = [], str(exc)
        self._finished.emit(request, arrays, error)

    def _on_finished(self, request, arrays, error):
        self.running = None
        if request.generation == self.generation:
            self.done.emit(request, arrays, error)
        else:
            self.stale += 1
        if self.pending is not None:
            pending, self.pending = self.pending, None
            self._start(pending)


def _is_int(text):
    try:
        int(text)
    except ValueError:
        return False
    return True


def track_kind(path):
    """"bigwig" by extension; otherwise from the first data line: two
    coordinate pairs (chrom start end chrom start end ..., as hicDetectLoops
    and BEDPE write them) are loops, chrom start end value is a bedGraph, and
    anything else with chrom start end is a BED of domains."""
    if path.lower().endswith((".bw", ".bigwig")):
        return "bigwig"
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 6 and all(_is_int(f) for f in (fields[1], fields[2], fields[4], fields[5])):
                return "loops"
            if len(fields) == 4 and _is_int(fields[1]) and _is_int(fields[2]):
                try:
                    float(fields[3])
                    return "bedgraph"
                except ValueError:
                    return "tads"
            return "tads"
    return "tads"


class Track:
    def __init__(self, path, kind):
        self.path = path
        self.kind = kind          # "tads", "loops", "bedgraph", "bigwig"
        self.data = {}
        self.handle = None
        if kind == "bigwig":
            import pyBigWig
            self.handle = pyBigWig.open(path)
        else:
            self._read()

    @property
    def name(self):
        return os.path.basename(self.path)

    def _read(self):
        opener = open
        if self.path.endswith(".gz"):
            import gzip
            opener = gzip.open
        with opener(self.path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                try:
                    if self.kind == "loops":
                        if fields[0] != fields[3]:
                            continue
                        item = (int(fields[1]), int(fields[2]), int(fields[4]), int(fields[5]))
                    elif self.kind == "bedgraph":
                        item = (int(fields[1]), int(fields[2]), float(fields[3]))
                    else:
                        item = (int(fields[1]), int(fields[2]))
                except (IndexError, ValueError):
                    continue
                self.data.setdefault(fields[0], []).append(item)

    def signal(self, chrom, start, end, bins):
        """(x, y) of a 1D track over [start, end)."""
        if self.kind == "bigwig":
            if chrom not in self.handle.chroms():
                return np.array([]), np.array([])
            end = min(end, self.handle.chroms()[chrom])
            if end <= start:
                return np.array([]), np.array([])
            bins = max(1, min(bins, end - start))
            values = np.array(self.handle.stats(chrom, int(start), int(end), nBins=bins), dtype=float)
            edges = np.linspace(start, end, bins + 1)
            return edges, values
        items = [i for i in self.data.get(chrom, []) if i[1] > start and i[0] < end]
        if not items:
            return np.array([]), np.array([])
        items.sort()
        x = np.array([v for i in items for v in (i[0], i[1])], dtype=float)
        y = np.array([i[2] for i in items for _ in (0, 1)], dtype=float)
        return x, y


class GenomeAxis(pg.AxisItem):
    """Tick labels in kb or Mb instead of scientific notation."""

    def tickStrings(self, values, scale, spacing):
        out = []
        for value in values:
            position = value * scale
            if spacing * scale >= 1e6 or abs(position) >= 1e7:
                out.append("{:g} Mb".format(round(position / 1e6, 3)))
            else:
                out.append("{:g} kb".format(round(position / 1e3, 1)))
        return out


def genome_axes():
    return {"left": GenomeAxis("left"), "bottom": GenomeAxis("bottom")}


class Panel:
    """One matrix view: a plot with the image and the overlays."""

    def __init__(self, layout, title):
        self.plot = layout.addPlot(title=title, axisItems=genome_axes())
        self.plot.setAspectLocked(True)
        self.plot.invertY(True)
        self.plot.showGrid(False, False)
        self.image = pg.ImageItem(axisOrder="row-major")
        self.plot.addItem(self.image)
        self.tads = pg.PlotDataItem(pen=pg.mkPen((0, 90, 200), width=1.5), connect="finite")
        self.loops = pg.ScatterPlotItem(symbol="s", size=10, pen=pg.mkPen((0, 150, 0), width=1.5),
                                        brush=None, pxMode=True)
        self.plot.addItem(self.tads)
        self.plot.addItem(self.loops)
        self.message = pg.TextItem("", color=(180, 0, 30), anchor=(0, 0))
        self.plot.addItem(self.message)


class MatrixBrowser(QtWidgets.QWidget):
    fetch_finished = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sources = [None, None]
        self.tracks = []
        self.chrom = None
        self.shown = []           # the arrays of the last applied fetch, as fetched
        self.shown_request = None
        self.error = ""
        self.fetcher = Fetcher(self)
        self.fetcher.done.connect(self._fetched)
        self._navigating = False

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        controls = QtWidgets.QGridLayout()
        controls.setHorizontalSpacing(6)
        self.open_a = QtWidgets.QPushButton("Open matrix A")
        self.open_a.clicked.connect(lambda: self._pick_matrix(0))
        self.open_b = QtWidgets.QPushButton("Open matrix B")
        self.open_b.clicked.connect(lambda: self._pick_matrix(1))
        self.add_track_button = QtWidgets.QPushButton("Add track")
        self.add_track_button.clicked.connect(self._pick_track)
        self.region = QtWidgets.QLineEdit()
        self.region.setPlaceholderText("chr:start-end")
        self.region.returnPressed.connect(lambda: self.goto(self.region.text()))
        self.go = QtWidgets.QPushButton("Go")
        self.go.clicked.connect(lambda: self.goto(self.region.text()))
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems(MODES)
        self.mode.currentTextChanged.connect(self._mode_changed)
        self.normalization = QtWidgets.QComboBox()
        self.normalization.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.normalization.setMinimumContentsLength(8)
        self.normalization.currentTextChanged.connect(lambda _t: self.refresh())
        self.log_scale = QtWidgets.QCheckBox("log")
        self.log_scale.setChecked(True)
        self.log_scale.toggled.connect(lambda _c: self._redraw())
        self.cmap = QtWidgets.QComboBox()
        self.cmap.addItems(COLORMAPS)
        self.cmap.currentTextChanged.connect(lambda _t: self._redraw())
        row1 = [self.open_a, self.open_b, self.add_track_button, QtWidgets.QLabel("Mode"), self.mode]
        row2 = [QtWidgets.QLabel("Region"), self.region, self.go, QtWidgets.QLabel("Normalization"),
                self.normalization, self.log_scale, QtWidgets.QLabel("Colors"), self.cmap]
        for column, widget in enumerate(row1):
            controls.addWidget(widget, 0, column)
        controls.addWidget(self.region, 1, 1, 1, 3)
        controls.addWidget(row2[0], 1, 0)
        controls.addWidget(self.go, 1, 4)
        second = QtWidgets.QHBoxLayout()
        for widget in row2[3:]:
            second.addWidget(widget)
        second.addStretch(1)
        outer.addLayout(controls)
        outer.addLayout(second)
        self.status = QtWidgets.QLabel("Open a cool, mcool, h5 or .hic file.")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)

        self.graphics = pg.GraphicsLayoutWidget()
        self.graphics.setBackground("w")
        outer.addWidget(self.graphics, 1)
        self.panels = []
        self.track_plot = None
        self.track_curves = []
        self._build_panels()
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(150)
        self.timer.timeout.connect(self._view_changed)

    # -- files ----------------------------------------------------------
    def _pick_matrix(self, slot):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open matrix", "", "Hi-C matrices (*.cool *.mcool *.h5 *.hic);;All files (*)")
        if path:
            self.open_matrix(path, slot)

    def _pick_track(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Add track", "", "Tracks (*.bed *.bedgraph *.bg *.bw *.bigwig *.bed.gz);;All files (*)")
        if path:
            self.add_track(path)

    def open_matrix(self, path, slot=0):
        try:
            source = MatrixSource(path)
        except Exception as exc:  # noqa: BLE001
            self._set_status("Cannot open {}: {}".format(path, exc), error=True)
            return None
        self.sources[slot] = source
        if slot == 0:
            self.normalization.blockSignals(True)
            self.normalization.clear()
            self.normalization.addItems(source.normalizations)
            self.normalization.blockSignals(False)
            if self.chrom is None or self.chrom not in source.lengths:
                chrom, length = source.chromosomes[0]
                self.region.setText("{}:0-{}".format(chrom, length))
        self._build_panels()
        if self.region.text():
            self.goto(self.region.text())
        return source

    def add_track(self, path, kind=None):
        if kind is None:
            try:
                kind = track_kind(path)
            except OSError as exc:
                self._set_status("Cannot read track {}: {}".format(path, exc), error=True)
                return None
        try:
            track = Track(path, kind)
        except Exception as exc:  # noqa: BLE001
            self._set_status("Cannot read track {}: {}".format(path, exc), error=True)
            return None
        self.tracks.append(track)
        self._build_panels()
        self._redraw()
        return track

    # -- layout ---------------------------------------------------------
    def _mode_changed(self, _text):
        self._build_panels()
        self.refresh()

    def _build_panels(self):
        self.graphics.clear()
        mode = self.mode.currentText()
        titles = ["A"]
        if mode == "side by side":
            titles = ["A", "B"]
        elif mode == "difference":
            titles = ["A - B"]
        self.panels = []
        for column, title in enumerate(titles):
            source = self.sources[0 if title != "B" else 1]
            label = title + (": " + source.name if source is not None and title != "A - B" else "")
            panel = Panel(self.graphics, label)
            panel.plot.sigRangeChanged.connect(self._range_changed)
            if self.panels:
                panel.plot.setXLink(self.panels[0].plot)
                panel.plot.setYLink(self.panels[0].plot)
            self.panels.append(panel)
        one_d = [t for t in self.tracks if t.kind in ("bedgraph", "bigwig")]
        self.track_plot = None
        self.track_curves = []
        if one_d:
            self.graphics.nextRow()
            self.track_plot = self.graphics.addPlot(colspan=len(self.panels),
                                                    axisItems={"bottom": GenomeAxis("bottom")})
            self.track_plot.setMaximumHeight(140)
            self.track_plot.setXLink(self.panels[0].plot)
            self.track_plot.addLegend(offset=(5, 5))
            palette = [(200, 60, 0), (0, 100, 180), (60, 140, 60), (120, 60, 160)]
            for index, track in enumerate(one_d):
                curve = self.track_plot.plot(stepMode="center", pen=pg.mkPen(palette[index % 4], width=1.2),
                                             name=track.name)
                self.track_curves.append((track, curve))
        if self.shown_request is not None:
            self._redraw()

    # -- navigation -----------------------------------------------------
    def goto(self, text):
        source = self.sources[0]
        if source is None:
            self._set_status("Open matrix A first.", error=True)
            return False
        try:
            chrom, start, end = parse_region_text(text, source.lengths)
        except ValueError as exc:
            self._set_status(str(exc), error=True)
            return False
        self.chrom = chrom
        self._goto_target = (chrom, start, end)
        self._apply_goto()
        # Panels that were just built are laid out later; apply the region
        # again then, so the view shows what was typed at any window size.
        QtCore.QTimer.singleShot(0, self._apply_goto)
        self._schedule(chrom, (start, end), (start, end))
        return True

    def _apply_goto(self):
        target = getattr(self, "_goto_target", None)
        if target is None or not self.panels or target[0] != self.chrom:
            return
        _, start, end = target
        self._navigating = True
        self.panels[0].plot.setRange(xRange=(start, end), yRange=(start, end), padding=0)
        self._navigating = False

    def zoom(self, factor):
        """Zooms the view about its centre (factor < 1 zooms in)."""
        if self.chrom is None:
            return
        self._goto_target = None
        (x0, x1), (y0, y1) = self.panels[0].plot.viewRange()
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        hw, hh = (x1 - x0) * factor / 2, (y1 - y0) * factor / 2
        self.panels[0].plot.setRange(xRange=(cx - hw, cx + hw), yRange=(cy - hh, cy + hh), padding=0)
        self._view_changed()

    def pan(self, fraction):
        if self.chrom is None:
            return
        self._goto_target = None
        (x0, x1), (y0, y1) = self.panels[0].plot.viewRange()
        dx, dy = (x1 - x0) * fraction, (y1 - y0) * fraction
        self.panels[0].plot.setRange(xRange=(x0 + dx, x1 + dx), yRange=(y0 + dy, y1 + dy), padding=0)
        self._view_changed()

    def _range_changed(self, *_):
        if not self._navigating and self.chrom is not None:
            if QtWidgets.QApplication.mouseButtons() != QtCore.Qt.NoButton:
                self._goto_target = None  # the user drags: follow the view
            self.timer.start()

    def _view_changed(self):
        if self.chrom is None or self.sources[0] is None:
            return
        length = self.sources[0].lengths[self.chrom]
        (x0, x1), (y0, y1) = self.panels[0].plot.viewRange()
        # The aspect ratio is locked, so a panel that is wider (or taller)
        # than square shows more of one axis. Fetch the centred square whose
        # side is the shorter axis: the region the user navigated to.
        side = min(x1 - x0, y1 - y0)
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        x0, x1, y0, y1 = cx - side / 2, cx + side / 2, cy - side / 2, cy + side / 2
        cols = (int(max(0, math.floor(x0))), int(min(length, math.ceil(x1))))
        rows = (int(max(0, math.floor(y0))), int(min(length, math.ceil(y1))))
        if cols[1] <= cols[0] or rows[1] <= rows[0]:
            return
        self.region.setText("{}:{:,}-{:,}".format(self.chrom, cols[0], cols[1]))
        self._schedule(self.chrom, rows, cols)

    def refresh(self):
        if self.shown_request is not None or self.chrom is not None:
            self._view_changed()

    def resolution_for(self, source, span, others=()):
        resolutions = source.resolutions
        for other in others:
            resolutions = [r for r in resolutions if r in other.resolutions]
        return choose_resolution(resolutions, span)

    def _schedule(self, chrom, rows, cols):
        mode = self.mode.currentText()
        a, b = self.sources
        normalization = self.normalization.currentText() or "none"
        span = max(rows[1] - rows[0], cols[1] - cols[0])
        if mode != "single" and b is None:
            self._set_status("Open matrix B for the {} view.".format(mode), error=True)
            return None
        involved = [a] if mode == "single" else [a, b]
        for source in involved:
            if chrom not in source.lengths:
                self._set_status("{} has no chromosome {}.".format(source.name, chrom), error=True)
                return None
            if normalization not in source.normalizations:
                self._set_status("{} has no normalization {!r} (available: {}).".format(
                    source.name, normalization, ", ".join(source.normalizations)), error=True)
                return None
        fetches = []
        resolution_used = None
        for index, source in enumerate(involved):
            others = [s for s in involved if s is not source] if mode == "difference" else []
            resolution = self.resolution_for(source, span, others)
            if mode == "difference" and not resolution and (a.resolutions or b.resolutions):
                self._set_status("A and B share no resolution.", error=True)
                return None
            if resolution is None:
                self._set_status("{} has variable bin sizes, which the browser does not show.".format(
                    source.name), error=True)
                return None
            if span / resolution > MAX_BINS:
                self._set_status("Zoom in: this view needs {:,} bins at {} bp, the browser reads at most {:,}."
                                 .format(int(span / resolution), resolution, MAX_BINS), error=True)
                return None
            r0, r1 = (rows[0] // resolution) * resolution, -(-rows[1] // resolution) * resolution
            c0, c1 = (cols[0] // resolution) * resolution, -(-cols[1] // resolution) * resolution
            length = source.lengths[chrom]
            region1 = "{}:{}-{}".format(chrom, r0, min(r1, length))
            region2 = "{}:{}-{}".format(chrom, c0, min(c1, length))
            fetches.append((source, region1, region2, resolution, normalization))
            resolution_used = resolution
        resolution = fetches[0][3]
        aligned_rows = ((rows[0] // resolution) * resolution, min(-(-rows[1] // resolution) * resolution,
                                                                  a.lengths[chrom]))
        aligned_cols = ((cols[0] // resolution) * resolution, min(-(-cols[1] // resolution) * resolution,
                                                                  a.lengths[chrom]))
        if (self.shown_request is not None and not self.error and self.shown_request.chrom == chrom
                and [f[1:] for f in self.shown_request.fetches] == [f[1:] for f in fetches]
                and [f[0] for f in self.shown_request.fetches] == [f[0] for f in fetches]):
            return self.shown_request
        self._set_status("Loading {} at {} bp".format(fetches[0][1], resolution_used))
        return self.fetcher.submit(fetches, chrom, aligned_rows, aligned_cols)

    # -- results --------------------------------------------------------
    def _fetched(self, request, arrays, error):
        self.error = error
        if error:
            self._set_status(error, error=True)
            self.fetch_finished.emit()
            return
        self.shown = arrays
        self.shown_request = request
        source, region1, region2, resolution, normalization = request.fetches[0]
        shape = arrays[0].shape
        self._set_status("{} x {} at {:,} bp ({} bins x {} bins, {})".format(
            region1, region2, resolution, shape[0], shape[1], normalization))
        self._redraw()
        self.fetch_finished.emit()

    def _redraw(self):
        request = self.shown_request
        if request is None or not self.shown:
            return
        mode = self.mode.currentText()
        images = list(self.shown)
        cmap_name = self.cmap.currentText()
        if mode == "difference" and len(images) == 2 and images[0].shape == images[1].shape:
            left, right = images
            if self.log_scale.isChecked():
                left, right = np.log1p(np.clip(left, 0, None)), np.log1p(np.clip(right, 0, None))
            images = [left - right]
            cmap_name = "blue-white-red"
            transform_log = False
        else:
            transform_log = self.log_scale.isChecked()
        for index, panel in enumerate(self.panels):
            if index >= len(images):
                break
            array = images[index]
            display = np.log1p(np.clip(array, 0, None)) if transform_log else array
            finite = display[np.isfinite(display)]
            display = np.where(np.isfinite(display), display, 0.0)
            if finite.size:
                if mode == "difference":
                    bound = float(np.percentile(np.abs(finite), 99)) or 1.0
                    levels = (-bound, bound)
                else:
                    high = float(np.percentile(finite, 99.5))
                    levels = (float(finite.min()), high if high > finite.min() else float(finite.min()) + 1)
            else:
                levels = (0, 1)
            panel.image.setImage(display, levels=levels, autoLevels=False)
            panel.image.setColorMap(colormap(cmap_name))
            _, region1, region2, resolution, _ = request.fetches[min(index, len(request.fetches) - 1)]
            r0, c0 = request.rows[0], request.cols[0]
            panel.image.setRect(QtCore.QRectF(c0, r0, array.shape[1] * resolution, array.shape[0] * resolution))
            panel.message.setText("")
            self._draw_overlays(panel, request)
        self._draw_tracks(request)

    def _draw_overlays(self, panel, request):
        chrom = request.chrom
        x, y, lx, ly = [], [], [], []
        for track in self.tracks:
            if track.kind == "tads":
                for start, end in track.data.get(chrom, []):
                    if end < request.cols[0] or start > request.cols[1]:
                        continue
                    x += [start, end, end, start, start, np.nan]
                    y += [start, start, end, end, start, np.nan]
            elif track.kind == "loops":
                for s1, e1, s2, e2 in track.data.get(chrom, []):
                    for cx, cy in (((s2 + e2) / 2, (s1 + e1) / 2), ((s1 + e1) / 2, (s2 + e2) / 2)):
                        if request.cols[0] <= cx <= request.cols[1] and request.rows[0] <= cy <= request.rows[1]:
                            lx.append(cx)
                            ly.append(cy)
        panel.tads.setData(np.array(x, dtype=float), np.array(y, dtype=float))
        panel.loops.setData(lx, ly)

    def _draw_tracks(self, request):
        if self.track_plot is None:
            return
        start, end = request.cols
        width = max(50, int(self.graphics.width()))
        for track, curve in self.track_curves:
            x, values = track.signal(request.chrom, start, end, width)
            if len(values):
                if track.kind == "bedgraph":
                    xs = x.reshape(-1, 2)
                    edges = np.concatenate([xs[:, 0], xs[-1:, 1]])
                    curve.setData(edges, values.reshape(-1, 2)[:, 0])
                else:
                    curve.setData(x, np.nan_to_num(values))
            else:
                curve.setData([], [])

    def _set_status(self, text, error=False):
        self.status.setText(text)
        self.status.setStyleSheet("color: #b00020;" if error else "")
        if error:
            self.error = text
            self.fetch_finished.emit()
