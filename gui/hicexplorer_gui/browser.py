"""Matrix browser (PLAN 10.5): cool, mcool, h5 and .hic through hicx_matrix.

Only the visible region is read and held. Navigation (the region box, zoom
and pan) schedules a fetch of the bins the view shows; the fetch runs on a
worker thread, a newer request replaces a waiting one, and the result of a
request that has been superseded is dropped, so the interface never blocks
and never shows a stale region. For mcool and .hic files the resolution is
chosen from the span on screen: the finest one that keeps the view within
MAX_BINS bins. A single-resolution file (cool, h5) whose view would exceed
MAX_BINS is not fetched; the browser asks to zoom in instead.

The axis ranges of every matrix view and of the track under it are the
region, in every mode: the views are sized square instead of locking the
aspect ratio (which would widen one axis), their ranges are synchronised
explicitly, navigation stays inside the chromosome, and overlays are clipped
to the region.
"""

import math
import os
import threading

import numpy as np
import pyqtgraph as pg
import shiboken6
from PySide6 import QtCore, QtGui, QtWidgets

import hicx_matrix

MAX_BINS = 1200
MODES = ("single", "side by side", "difference")
COLORMAPS = ("white-red", "viridis", "inferno", "magma", "cividis")
AXIS_WIDTH = 72
TRACK_HEIGHT = 130


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
    def __init__(self, generation, fetches, chrom, rows, cols, view_rows, view_cols):
        self.generation = generation
        self.fetches = fetches    # [(source, region1, region2, resolution, normalization)]
        self.chrom = chrom
        self.rows = rows          # (start, end) bp of the fetched rows, bin aligned
        self.cols = cols
        self.view_rows = view_rows  # (start, end) bp of the region on screen
        self.view_cols = view_cols


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

    def submit(self, fetches, chrom, rows, cols, view_rows, view_cols):
        self.generation += 1
        request = Request(self.generation, fetches, chrom, rows, cols, view_rows, view_cols)
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
        try:
            self._finished.emit(request, arrays, error)
        except RuntimeError:
            # The browser was closed while this fetch ran; nobody needs it.
            pass

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
                        score = None
                        if len(fields) > 6:
                            try:
                                score = float(fields[6])
                            except ValueError:
                                score = None
                        item = (int(fields[1]), int(fields[2]), int(fields[4]), int(fields[5]), score)
                    elif self.kind == "bedgraph":
                        item = (int(fields[1]), int(fields[2]), float(fields[3]))
                    else:
                        item = (int(fields[1]), int(fields[2]))
                except (IndexError, ValueError):
                    continue
                self.data.setdefault(fields[0], []).append(item)

    def signal(self, chrom, start, end, bins):
        """(edges, values) of a 1D track over [start, end): len(edges) is
        len(values) + 1 and every edge lies inside the interval."""
        if self.kind == "bigwig":
            if chrom not in self.handle.chroms():
                return np.array([]), np.array([])
            end = min(end, self.handle.chroms()[chrom])
            if end <= start:
                return np.array([]), np.array([])
            bins = max(1, min(bins, int(end - start)))
            values = np.array(self.handle.stats(chrom, int(start), int(end), nBins=bins), dtype=float)
            return np.linspace(start, end, bins + 1), values
        items = sorted(i for i in self.data.get(chrom, []) if i[1] > start and i[0] < end)
        if not items:
            return np.array([]), np.array([])
        edges = np.array([i[0] for i in items] + [items[-1][1]], dtype=float)
        return np.clip(edges, start, end), np.array([i[2] for i in items], dtype=float)


def clipped_tad_lines(domains, rows, cols):
    """(x, y) with NaN breaks of the domain squares' edges inside the region."""
    r0, r1 = rows
    c0, c1 = cols
    x, y = [], []

    def horizontal(at, a, b):
        if r0 <= at <= r1:
            a, b = max(a, c0), min(b, c1)
            if a < b:
                x.extend([a, b, np.nan])
                y.extend([at, at, np.nan])

    def vertical(at, a, b):
        if c0 <= at <= c1:
            a, b = max(a, r0), min(b, r1)
            if a < b:
                x.extend([at, at, np.nan])
                y.extend([a, b, np.nan])

    for start, end in domains:
        horizontal(start, start, end)
        horizontal(end, start, end)
        vertical(start, start, end)
        vertical(end, start, end)
    return np.array(x, dtype=float), np.array(y, dtype=float)


ARC_SAMPLES = 24


def loop_arc_path(loops, cols):
    """(x, y) with NaN breaks tracing one parabolic arc per loop whose two
    anchor midpoints both lie inside ``cols`` (a single linear axis, matching
    the matrix's x-axis). Apex height is scaled by the loop's score when one
    was parsed, falling back to its span (e2 - s1); heights are normalised to
    the tallest arc drawn so the track always uses its own available height."""
    kept, heights = [], []
    c0, c1 = cols
    for s1, e1, s2, e2, score in loops:
        x1, x2 = (s1 + e1) / 2.0, (s2 + e2) / 2.0
        if c0 <= x1 <= c1 and c0 <= x2 <= c1:
            kept.append((x1, x2))
            heights.append(abs(score) if score is not None else float(e2 - s1))
    if not kept:
        return np.array([]), np.array([])
    peak = max(heights) or 1.0
    x, y = [], []
    t = np.linspace(-1.0, 1.0, ARC_SAMPLES)
    parabola = 1.0 - t ** 2
    for (x1, x2), height in zip(kept, heights):
        left, right = min(x1, x2), max(x1, x2)
        mid, half = (left + right) / 2.0, max(right - left, 1.0) / 2.0
        x.extend((mid + t * half).tolist())
        x.append(np.nan)
        y.extend((parabola * (height / peak)).tolist())
        y.append(np.nan)
    return np.array(x, dtype=float), np.array(y, dtype=float)


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


def _plot_item(title=None, left=True):
    axes = {"bottom": GenomeAxis("bottom")}
    if left:
        axes["left"] = GenomeAxis("left")
    plot = pg.PlotItem(title=title, axisItems=axes)
    plot.hideButtons()
    plot.setMenuEnabled(False)
    plot.getAxis("left").setWidth(AXIS_WIDTH)
    plot.vb.setDefaultPadding(0.0)
    return plot


class Panel:
    """One matrix view: a plot with the image and the overlays."""

    def __init__(self, title):
        self.title = title
        self.plot = _plot_item(title)
        self.plot.invertY(True)
        self.image = pg.ImageItem(axisOrder="row-major")
        self.plot.addItem(self.image)
        self.tads = pg.PlotDataItem(pen=pg.mkPen((0, 90, 200), width=1.5), connect="finite")
        self.plot.addItem(self.tads)


class MatrixBrowser(QtWidgets.QWidget):
    fetch_finished = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sources = [None, None]
        self.tracks = []
        self.chrom = None
        self.view_x = None        # (start, end) bp of the columns on screen
        self.view_y = None        # (start, end) bp of the rows on screen
        self.shown = []           # the arrays of the last applied fetch, as fetched
        self.shown_request = None
        self.error = ""
        self.fetcher = Fetcher(self)
        self.fetcher.done.connect(self._fetched)
        self._navigating = False
        self._fit_attempts = 0

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
        for column, widget in enumerate([self.open_a, self.open_b, self.add_track_button,
                                         QtWidgets.QLabel("Mode"), self.mode]):
            controls.addWidget(widget, 0, column)
        controls.addWidget(QtWidgets.QLabel("Region"), 1, 0)
        controls.addWidget(self.region, 1, 1, 1, 3)
        controls.addWidget(self.go, 1, 4)
        second = QtWidgets.QHBoxLayout()
        for widget in (QtWidgets.QLabel("Normalization"), self.normalization, self.log_scale,
                       QtWidgets.QLabel("Colors"), self.cmap):
            second.addWidget(widget)
        second.addStretch(1)
        outer.addLayout(controls)
        outer.addLayout(second)
        self.status = QtWidgets.QLabel("Open a cool, mcool, h5 or .hic file.")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)

        self.graphics = pg.GraphicsLayoutWidget()
        self.graphics.setBackground("w")
        self.graphics.installEventFilter(self)
        outer.addWidget(self.graphics, 1)
        self.panels = []
        self.track_plot = None
        self.track_curves = []
        self.arc_plot = None
        self.arc_curves = []
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(150)
        self.timer.timeout.connect(self._view_changed)
        self._build_panels()

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
                self.chrom, self.view_x, self.view_y = None, None, None
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
            panel = Panel(label)
            self.graphics.addItem(panel.plot, row=0, col=column)
            panel.plot.sigRangeChanged.connect(lambda *_args, p=panel: self._range_changed(p))
            self.panels.append(panel)
        one_d = [t for t in self.tracks if t.kind in ("bedgraph", "bigwig")]
        loop_tracks = [t for t in self.tracks if t.kind == "loops"]
        self.track_plot = None
        self.track_curves = []
        self.arc_plot = None
        self.arc_curves = []
        row = 1
        if one_d:
            self.track_plot = _plot_item(left=False)
            self.graphics.addItem(self.track_plot, row=row, col=0)
            self.track_plot.setMouseEnabled(x=False, y=False)
            self.track_plot.addLegend(offset=(5, 5))
            palette = [(200, 60, 0), (0, 100, 180), (60, 140, 60), (120, 60, 160)]
            for index, track in enumerate(one_d):
                curve = self.track_plot.plot(stepMode="center", pen=pg.mkPen(palette[index % 4], width=1.2),
                                             name=track.name)
                self.track_curves.append((track, curve))
            row += 1
        if loop_tracks:
            # Arc track: replaces the old point overlay of loop anchors drawn
            # directly on the matrix (PLAN 10.5); one parabolic arc per loop,
            # laid out as its own panel row so it shares the matrix's x-axis
            # like the bedgraph/bigwig tracks above.
            self.arc_plot = _plot_item(left=False)
            self.graphics.addItem(self.arc_plot, row=row, col=0)
            self.arc_plot.setMouseEnabled(x=False, y=False)
            self.arc_plot.addLegend(offset=(5, 5))
            palette = [(160, 20, 130), (0, 130, 130), (180, 110, 0)]
            for index, track in enumerate(loop_tracks):
                curve = self.arc_plot.plot(pen=pg.mkPen(palette[index % len(palette)], width=1.5),
                                           connect="finite", name=track.name)
                self.arc_curves.append((track, curve))
        self._apply_limits()
        if self.view_x is not None:
            self._set_view(self.view_x, self.view_y)
        if self.shown_request is not None:
            self._redraw()
        self._fit_attempts = 0
        QtCore.QTimer.singleShot(0, self._fit)

    def eventFilter(self, watched, event):
        if watched is self.graphics and event.type() == QtCore.QEvent.Resize:
            self._fit_attempts = 0
            QtCore.QTimer.singleShot(0, self._fit)
        return False

    def _fit(self):
        """Sizes the matrix views square and the track as wide as view A."""
        # Deferred through QTimer.singleShot(0, ...), so it can fire after the
        # browser (e.g. a project tab being closed) was already deleted.
        if not shiboken6.isValid(self) or not self.panels:
            return
        width = self.graphics.width() - 16
        height = self.graphics.height() - 16
        if width < 60 or height < 60:
            return
        plot = self.panels[0].plot
        # Space the axes and the title take around the view. A view wider than
        # its plot (a title not yet elided) would make the measure too small,
        # so it never goes below the fixed left axis width.
        dw, dh = AXIS_WIDTH + 12, 64
        if shiboken6.isValid(plot.vb) and plot.vb.width() > 1 and plot.vb.height() > 1:
            dw = max(AXIS_WIDTH + 2, plot.size().width() - plot.vb.width())
            dh = max(40, plot.size().height() - plot.vb.height())
        n = len(self.panels)
        extra_rows = (self.track_plot is not None) + (self.arc_plot is not None)
        track = (TRACK_HEIGHT + 10) * extra_rows
        side = int(max(40, min(width / n - dw - 6 * n, height - track - dh)))
        for panel in self.panels:
            self._elide_title(panel, side)
            size = QtCore.QSizeF(side + dw, side + dh)
            panel.plot.setMinimumSize(size)
            panel.plot.setMaximumSize(size)
        for extra in (self.track_plot, self.arc_plot):
            if extra is not None:
                size = QtCore.QSizeF(side + dw, TRACK_HEIGHT)
                extra.setMinimumSize(size)
                extra.setMaximumSize(size)
        self._fit_attempts += 1
        if self._fit_attempts < 4:
            QtCore.QTimer.singleShot(0, self._check_square)

    @staticmethod
    def _elide_title(panel, width):
        """A plot title is the minimum width of its layout column, so a long
        file name would widen the view past the square. The title is elided
        until its rendered width fits; the tooltip keeps the full name."""
        label = panel.plot.titleLabel
        label.item.setToolTip(panel.title)
        metrics = QtGui.QFontMetrics(label.item.font())
        budget = int(width * 0.95)
        while True:
            text = metrics.elidedText(panel.title, QtCore.Qt.ElideMiddle, max(budget, 1))
            if text != label.text:
                panel.plot.setTitle(text)
            if label.itemRect().width() <= width or budget <= 20:
                return
            budget = int(budget * 0.85)

    def _check_square(self):
        if not shiboken6.isValid(self) or not self.panels:
            return
        vb = self.panels[0].plot.vb
        if shiboken6.isValid(vb) and abs(vb.width() - vb.height()) > 1:
            self._fit()

    def _apply_limits(self):
        if self.chrom is None or self.sources[0] is None or self.chrom not in self.sources[0].lengths:
            return
        length = self.sources[0].lengths[self.chrom]
        self._navigating = True
        try:
            for panel in self.panels:
                panel.plot.vb.setLimits(xMin=0, xMax=length, yMin=0, yMax=length,
                                        maxXRange=length, maxYRange=length)
            if self.track_plot is not None:
                self.track_plot.vb.setLimits(xMin=0, xMax=length, maxXRange=length)
            if self.arc_plot is not None:
                self.arc_plot.vb.setLimits(xMin=0, xMax=length, maxXRange=length)
        finally:
            self._navigating = False

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
        self._apply_limits()
        self._set_view((start, end), (start, end))
        self._view_changed()
        return True

    def _clamp(self, lo, hi):
        length = float(self.sources[0].lengths[self.chrom])
        span = min(max(hi - lo, 1.0), length)
        lo = min(max(0.0, lo), length - span)
        return (lo, lo + span)

    def _set_view(self, x_range, y_range):
        """Shows exactly this region in every matrix view and on the track."""
        x_range = self._clamp(*x_range)
        y_range = self._clamp(*y_range)
        self.view_x, self.view_y = x_range, y_range
        self._navigating = True
        try:
            for panel in self.panels:
                panel.plot.vb.setRange(xRange=x_range, yRange=y_range, padding=0)
            if self.track_plot is not None:
                self.track_plot.vb.setXRange(*x_range, padding=0)
            if self.arc_plot is not None:
                self.arc_plot.vb.setXRange(*x_range, padding=0)
        finally:
            self._navigating = False
        self.region.setText("{}:{:,}-{:,}".format(self.chrom, int(round(x_range[0])), int(round(x_range[1]))))

    def zoom(self, factor):
        """Zooms the view about its centre (factor < 1 zooms in)."""
        if self.view_x is None:
            return
        cx, cy = sum(self.view_x) / 2, sum(self.view_y) / 2
        hx = (self.view_x[1] - self.view_x[0]) * factor / 2
        hy = (self.view_y[1] - self.view_y[0]) * factor / 2
        self._set_view((cx - hx, cx + hx), (cy - hy, cy + hy))
        self._view_changed()

    def pan(self, fraction):
        if self.view_x is None:
            return
        dx = (self.view_x[1] - self.view_x[0]) * fraction
        dy = (self.view_y[1] - self.view_y[0]) * fraction
        self._set_view((self.view_x[0] + dx, self.view_x[1] + dx), (self.view_y[0] + dy, self.view_y[1] + dy))
        self._view_changed()

    def _range_changed(self, panel):
        """Mouse zoom or pan in one view: every view and the track follow."""
        if self._navigating or self.chrom is None or self.sources[0] is None:
            return
        (x0, x1), (y0, y1) = panel.plot.vb.viewRange()
        self._set_view((x0, x1), (y0, y1))
        self.timer.start()

    def _view_changed(self):
        if self.chrom is None or self.sources[0] is None or self.view_x is None:
            return
        length = self.sources[0].lengths[self.chrom]
        cols = (int(max(0, math.floor(self.view_x[0]))), int(min(length, math.ceil(self.view_x[1]))))
        rows = (int(max(0, math.floor(self.view_y[0]))), int(min(length, math.ceil(self.view_y[1]))))
        if cols[1] <= cols[0] or rows[1] <= rows[0]:
            return
        self._schedule(self.chrom, rows, cols)

    def refresh(self):
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
        for source in involved:
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
            length = source.lengths[chrom]
            r0, r1 = (rows[0] // resolution) * resolution, -(-rows[1] // resolution) * resolution
            c0, c1 = (cols[0] // resolution) * resolution, -(-cols[1] // resolution) * resolution
            region1 = "{}:{}-{}".format(chrom, r0, min(r1, length))
            region2 = "{}:{}-{}".format(chrom, c0, min(c1, length))
            fetches.append((source, region1, region2, resolution, normalization))
        resolution = fetches[0][3]
        length = a.lengths[chrom]
        aligned_rows = ((rows[0] // resolution) * resolution, min(-(-rows[1] // resolution) * resolution, length))
        aligned_cols = ((cols[0] // resolution) * resolution, min(-(-cols[1] // resolution) * resolution, length))
        current = self.shown_request
        if (current is not None and not self.error and current.chrom == chrom
                and [f[1:] for f in current.fetches] == [f[1:] for f in fetches]
                and [f[0] for f in current.fetches] == [f[0] for f in fetches]):
            current.view_rows, current.view_cols = rows, cols
            self._redraw()
            return current
        self._set_status("Loading {} at {:,} bp".format(fetches[0][1], resolution))
        return self.fetcher.submit(fetches, chrom, aligned_rows, aligned_cols, rows, cols)

    # -- results --------------------------------------------------------
    def _fetched(self, request, arrays, error):
        self.error = error
        if error:
            self._set_status(error, error=True)
            self.fetch_finished.emit()
            return
        self.shown = arrays
        self.shown_request = request
        _, region1, region2, resolution, normalization = request.fetches[0]
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
        transform_log = self.log_scale.isChecked()
        if mode == "difference" and len(images) == 2 and images[0].shape == images[1].shape:
            left, right = images
            if transform_log:
                left, right = np.log1p(np.clip(left, 0, None)), np.log1p(np.clip(right, 0, None))
            images = [left - right]
            cmap_name = "blue-white-red"
            transform_log = False
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
                    low = float(finite.min())
                    levels = (low, high if high > low else low + 1)
            else:
                levels = (0, 1)
            panel.image.setImage(display, levels=levels, autoLevels=False)
            panel.image.setColorMap(colormap(cmap_name))
            resolution = request.fetches[min(index, len(request.fetches) - 1)][3]
            panel.image.setRect(QtCore.QRectF(request.cols[0], request.rows[0],
                                              array.shape[1] * resolution, array.shape[0] * resolution))
            self._draw_overlays(panel, request)
        self._draw_tracks(request)
        self._draw_arcs(request)

    def _draw_overlays(self, panel, request):
        rows, cols = request.view_rows, request.view_cols
        domains = []
        for track in self.tracks:
            if track.kind == "tads":
                domains.extend(track.data.get(request.chrom, []))
        x, y = clipped_tad_lines(domains, rows, cols)
        panel.tads.setData(x, y)

    def _draw_arcs(self, request):
        if self.arc_plot is None:
            return
        cols = request.view_cols
        for track, curve in self.arc_curves:
            x, y = loop_arc_path(track.data.get(request.chrom, []), cols)
            curve.setData(x, y)
        if self.view_x is not None:
            self._navigating = True
            try:
                self.arc_plot.vb.setXRange(*self.view_x, padding=0)
            finally:
                self._navigating = False

    def _draw_tracks(self, request):
        if self.track_plot is None:
            return
        start, end = request.view_cols
        width = max(50, int(self.track_plot.vb.width()) or 50)
        for track, curve in self.track_curves:
            edges, values = track.signal(request.chrom, start, end, width)
            if len(values):
                curve.setData(edges, np.nan_to_num(values))
            else:
                curve.setData([], [])
        if self.view_x is not None:
            self._navigating = True
            try:
                self.track_plot.vb.setXRange(*self.view_x, padding=0)
            finally:
                self._navigating = False

    def _set_status(self, text, error=False):
        self.status.setText(text)
        self.status.setStyleSheet("color: #d0314b;" if error else "")
        if error:
            self.error = text
            self.fetch_finished.emit()
