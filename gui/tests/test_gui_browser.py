"""Matrix browser (PLAN 10.5), offscreen.

- The arrays shown equal hicx_matrix fetches of the regions the browser
  reports, and are E2 against the reference readers (cooler, hicmatrix,
  hicstraw; the oracle of cpp/python/tests).
- Zooming switches resolution on mcool and .hic.
- Side-by-side and difference views, stale fetches, tracks.
- A scripted browsing session stays within a peak RSS bound.
"""

import json
import os
import resource
import subprocess
import sys
import time

import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")
pytest.importorskip("pyqtgraph")

from gui_support import DATA, GUI_DIR, REPO, import_hicx_matrix  # noqa: E402

hicx_matrix = import_hicx_matrix()

from hicexplorer_gui import browser as browser_module  # noqa: E402
from hicexplorer_gui.browser import MatrixBrowser, Track  # noqa: E402

sys.path.insert(0, os.path.join(REPO, "cpp", "python", "tests"))
import support  # noqa: E402

LARGE_HIC = os.environ.get("HICX_LARGE_HIC")
needs_large = pytest.mark.skipif(not LARGE_HIC or not os.path.isfile(LARGE_HIC),
                                 reason="HICX_LARGE_HIC is not set or missing")
# Peak RSS bound of the scripted session (browse_session.py), in kB.
SESSION_RSS_BOUND_KB = 400 * 1024


@pytest.fixture
def browser(qtbot):
    widget = MatrixBrowser()
    qtbot.addWidget(widget)
    widget.resize(800, 800)
    widget.show()
    return widget


def settle(qtbot, widget, timeout=120000):
    def idle():
        return not widget.fetcher.busy and not widget.timer.isActive()
    qtbot.waitUntil(idle, timeout=timeout)
    qtbot.wait(50)
    qtbot.waitUntil(idle, timeout=timeout)


def parse(region):
    chrom, span = region.replace(",", "").rsplit(":", 1)
    start, end = span.split("-")
    return chrom, int(start), int(end)


def check_against_module(widget):
    """Every shown array equals a fresh hicx_matrix fetch of its request."""
    assert widget.shown_request is not None and not widget.error, widget.error
    for array, (source, region1, region2, resolution, normalization) in zip(
            widget.shown, widget.shown_request.fetches):
        fresh = hicx_matrix.open(source.path).fetch(region1, region2, resolution=resolution,
                                                    normalization=normalization)
        assert np.array_equal(array, fresh, equal_nan=True)
    return widget.shown_request.fetches


def show(qtbot, widget, path, region, normalization="none"):
    widget.open_matrix(path)
    settle(qtbot, widget)
    widget.normalization.setCurrentText(normalization)
    settle(qtbot, widget)
    assert widget.goto(region)
    settle(qtbot, widget)
    return check_against_module(widget)


def test_cool_region_matches_module_and_cooler(qtbot, browser, tmp_path):
    path = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    (source, region1, region2, resolution, _), = show(qtbot, browser, path, "1:10,000,000-12,000,000")
    assert resolution == 10000
    oracle = support.run_oracle("cooler", [{"id": "q", "what": "matrix", "uri": path, "balance": False,
                                            "region1": region1, "region2": region2}], tmp_path)
    support.assert_e2(browser.shown[0], oracle["q"]["matrix"], "cooler")


def test_mcool_switches_resolution_and_matches_cooler(qtbot, browser, tmp_path):
    path = os.path.join(DATA, "hicBuildMatrix", "multi_small_test_matrix.mcool")
    queries, arrays, used = [], [], []
    for region in ("chr2L", "chr2L:1,000,000-2,000,000"):
        (source, region1, region2, resolution, _), = show(qtbot, browser, path, region)
        used.append(resolution)
        queries.append({"id": "q{}".format(len(queries)), "what": "matrix", "balance": False,
                        "uri": "{}::/resolutions/{}".format(path, resolution),
                        "region1": region1, "region2": region2})
        arrays.append(browser.shown[0])
    assert used[0] > used[1], used
    assert used == [20000, 5000], used
    oracle = support.run_oracle("cooler", queries, tmp_path)
    for query, array in zip(queries, arrays):
        support.assert_e2(array, oracle[query["id"]]["matrix"], "cooler")


def test_h5_region_matches_module_and_hicmatrix(qtbot, browser, tmp_path):
    path = os.path.join(DATA, "hicDifferentialTAD", "GSM2644945_Untreated-R1.100000_chr1.h5")
    (source, region1, region2, resolution, _), = show(qtbot, browser, path, "chr1:50,000,000-80,000,000")
    oracle = support.run_oracle("hicmatrix", [{"id": "q", "what": "matrix", "path": path,
                                               "region1": list(parse(region1)), "region2": list(parse(region2))}],
                                tmp_path)
    support.assert_e2(browser.shown[0], oracle["q"]["matrix"], "hicmatrix")


def hicstraw_query(query_id, path, fetch, norm):
    source, region1, region2, resolution, _ = fetch
    return {"id": query_id, "what": "matrix", "path": path, "resolution": resolution, "norm": norm,
            "matrix_type": "observed", "region1": list(parse(region1)), "region2": list(parse(region2))}


@pytest.mark.parametrize("normalization,norm", [("none", "NONE"), ("KR", "KR")])
def test_hic_region_matches_module_and_hicstraw(qtbot, browser, tmp_path, normalization, norm):
    path = os.path.join(DATA, "hicHyperoptDetectLoopsHiCCUPS", "SRR1791297_30.hic")
    fetch, = show(qtbot, browser, path, "NC_001134.8:100,000-700,000", normalization)
    oracle = support.run_oracle("hicstraw", [hicstraw_query("q", path, fetch, norm)], tmp_path)
    support.assert_e2(browser.shown[0], oracle["q"]["matrix"], "hicstraw")


@needs_large
def test_large_hic_region_and_zoom_resolution_switching(qtbot, browser, tmp_path):
    fetch, = show(qtbot, browser, LARGE_HIC, "1:50,000,000-52,000,000", "KR")
    oracle = support.run_oracle("hicstraw", [hicstraw_query("q", LARGE_HIC, fetch, "KR")], tmp_path)
    support.assert_e2(browser.shown[0], oracle["q"]["matrix"], "hicstraw")
    assert browser.goto("1")
    settle(qtbot, browser)
    used = [browser.shown_request.fetches[0][3]]
    for _ in range(2):
        browser.zoom(0.1)
        settle(qtbot, browser)
        check_against_module(browser)
        used.append(browser.shown_request.fetches[0][3])
    assert used[0] > used[1] > used[2], used
    for array in browser.shown:
        assert max(array.shape) <= browser_module.MAX_BINS


def test_side_by_side_and_difference(qtbot, browser):
    # The chr1_chr2 coolers store their 'format' attribute as a fixed length
    # string, which cooler.fileops.is_cooler rejects and cooler.Cooler opens.
    a = os.path.join(DATA, "hicDifferentialTAD", "GSM2644945_Untreated-R1.100000_chr1_chr2.cool")
    b = os.path.join(DATA, "hicDifferentialTAD", "GSM2644947_Auxin2days-R1.100000_chr1_chr2.cool")
    browser.open_matrix(a, 0)
    browser.open_matrix(b, 1)
    browser.mode.setCurrentText("side by side")
    assert browser.goto("chr1:20,000,000-60,000,000")
    settle(qtbot, browser)
    fetches = check_against_module(browser)
    assert len(browser.panels) == 2 and [f[0].path for f in fetches] == [a, b]
    browser.mode.setCurrentText("difference")
    settle(qtbot, browser)
    fetches = check_against_module(browser)
    assert len(browser.panels) == 1 and len(fetches) == 2
    image = browser.panels[0].image.image
    assert image.shape == browser.shown[0].shape
    expected = np.log1p(np.clip(browser.shown[0], 0, None)) - np.log1p(np.clip(browser.shown[1], 0, None))
    assert np.allclose(image, np.where(np.isfinite(expected), expected, 0.0))


def test_navigation_never_blocks_and_stale_fetches_are_dropped(qtbot, browser):
    path = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    browser.open_matrix(path)
    settle(qtbot, browser)
    stale_before = browser.fetcher.stale
    regions = ["1:{},000,000-{},000,000".format(s, s + 8) for s in range(10, 110, 10)]
    durations = []
    for region in regions:
        start = time.perf_counter()
        browser.goto(region)
        durations.append(time.perf_counter() - start)
    assert max(durations) < 0.2, durations
    settle(qtbot, browser)
    assert browser.fetcher.stale - stale_before >= len(regions) - 2
    chrom, start, end = parse(regions[-1])
    request = browser.shown_request
    assert request.rows[0] <= start and request.rows[1] >= end - 10000
    check_against_module(browser)


def test_single_resolution_file_refuses_a_too_large_view(qtbot, browser):
    browser.open_matrix(os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool"))
    settle(qtbot, browser)
    browser.goto("1")
    settle(qtbot, browser)
    assert "Zoom in" in browser.status.text()


def test_tracks_overlay_and_signal(qtbot, browser, tmp_path):
    import pyBigWig
    tads = tmp_path / "tads.bed"
    tads.write_text("1\t112000000\t114000000\n1\t114000000\t117000000\n2\t0\t100\n")
    graph = tmp_path / "signal.bedgraph"
    graph.write_text("".join("1\t{}\t{}\t{}\n".format(s, s + 100000, (s // 100000) % 7)
                             for s in range(111000000, 119000000, 100000)))
    wig = str(tmp_path / "signal.bw")
    handle = pyBigWig.open(wig, "w")
    handle.addHeader([("1", 249250621)])
    handle.addEntries(["1"] * 3, [112000000, 114000000, 116000000], ends=[113000000, 115000000, 117000000],
                      values=[1.0, 2.0, 3.0])
    handle.close()
    browser.open_matrix(os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool"))
    for path in (str(tads), os.path.join(DATA, "hicDetectLoops", "loops.bedgraph"), str(graph), wig):
        assert browser.add_track(path) is not None
    assert [t.kind for t in browser.tracks] == ["tads", "loops", "bedgraph", "bigwig"]
    browser.goto("1:111,000,000-119,000,000")
    settle(qtbot, browser)
    panel = browser.panels[0]
    assert np.count_nonzero(np.isfinite(panel.tads.xData)) >= 10
    assert len(panel.loops.data) >= 1
    curves = dict((t.kind, c) for t, c in browser.track_curves)
    assert len(curves["bedgraph"].yData) > 0 and np.nanmax(curves["bedgraph"].yData) == 6
    assert np.nanmax(curves["bigwig"].yData) == 3.0
    assert Track(wig, "bigwig").signal("1", 112000000, 113000000, 4)[1].tolist() == [1.0, 1.0, 1.0, 1.0]


def assert_view_is_region(widget, start, end):
    """Axis ranges of every view and of the track equal the region, views are
    square, and no overlay or track item lies outside the region."""
    for panel in widget.panels:
        (x0, x1), (y0, y1) = panel.plot.vb.viewRange()
        assert (x0, x1) == pytest.approx((start, end), abs=1.0)
        assert (y0, y1) == pytest.approx((start, end), abs=1.0)
        vb = panel.plot.vb
        assert abs(vb.width() - vb.height()) <= 2, (vb.width(), vb.height())
        # the view lies inside its own plot, so neighbouring views never overlap
        inner = vb.mapToScene(vb.rect()).boundingRect()
        outer = panel.plot.geometry()
        assert outer.left() - 1 <= inner.left() and inner.right() <= outer.right() + 1, (inner, outer)
        xs, ys = panel.tads.xData, panel.tads.yData
        if xs is not None and len(xs):
            finite = np.isfinite(xs)
            assert start <= xs[finite].min() and xs[finite].max() <= end
            assert start <= ys[finite].min() and ys[finite].max() <= end
        for spot in panel.loops.points():
            assert start <= spot.pos().x() <= end and start <= spot.pos().y() <= end
    if widget.track_plot is not None:
        assert tuple(widget.track_plot.vb.viewRange()[0]) == pytest.approx((start, end), abs=1.0)
        track_vb, matrix_vb = widget.track_plot.vb, widget.panels[0].plot.vb
        left_track = track_vb.mapToScene(track_vb.rect().topLeft()).x()
        left_matrix = matrix_vb.mapToScene(matrix_vb.rect().topLeft()).x()
        assert abs(left_track - left_matrix) <= 2 and abs(track_vb.width() - matrix_vb.width()) <= 2
        for _track, curve in widget.track_curves:
            if curve.xData is not None and len(curve.xData):
                assert start <= curve.xData.min() and curve.xData.max() <= end


@pytest.mark.parametrize("mode", ["single", "side by side", "difference"])
def test_axes_and_track_follow_the_region_in_every_mode(qtbot, tmp_path, mode):
    widget = MatrixBrowser()
    qtbot.addWidget(widget)
    widget.resize(1280, 720)
    widget.show()
    tads = tmp_path / "tads.bed"
    tads.write_text("1\t108000000\t111000000\n1\t111000000\t114000000\n1\t118000000\t121000000\n")
    graph = tmp_path / "signal.bedgraph"
    graph.write_text("".join("1\t{}\t{}\t{}\n".format(s, s + 250000, (s // 250000) % 5)
                             for s in range(100000000, 130000000, 250000)))
    # a file name wider than a small view: titles must not widen the views
    cool = str(tmp_path / "a_matrix_file_name_much_wider_than_a_small_view_gm12878_chr1.cool")
    os.symlink(os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool"), cool)
    widget.open_matrix(cool, 0)
    widget.open_matrix(cool, 1)
    for path in (str(tads), os.path.join(DATA, "hicDetectLoops", "loops.bedgraph"), str(graph)):
        assert widget.add_track(path) is not None
    widget.mode.setCurrentText(mode)
    assert widget.goto("1:110,000,000-120,000,000")
    settle(qtbot, widget)
    qtbot.wait(100)
    assert not widget.error, widget.error
    assert widget.shown_request.view_cols == (110000000, 120000000)
    assert_view_is_region(widget, 110000000, 120000000)
    widget.zoom(0.5)
    settle(qtbot, widget)
    assert_view_is_region(widget, 112500000, 117500000)
    widget.pan(0.4)
    settle(qtbot, widget)
    assert_view_is_region(widget, 114500000, 119500000)
    for size in ((1920, 1080), (960, 520)):
        widget.resize(*size)
        settle(qtbot, widget)
        qtbot.wait(100)
        assert_view_is_region(widget, 114500000, 119500000)


def test_scripted_session_peak_rss(tmp_path):
    # The session's fetch count and its memory bound are defined for the
    # scripted walk over both the cool and the large .hic file.
    if not (LARGE_HIC and os.path.isfile(LARGE_HIC)):
        pytest.skip("HICX_LARGE_HIC does not name an existing .hic file; the scripted "
                    "session needs it for its fetch count and memory bound")
    cool = os.path.join(DATA, "hicTADClassifier", "gm12878_chr1.cool")
    argv = [sys.executable, os.path.join(GUI_DIR, "tests", "browse_session.py"), cool]
    if LARGE_HIC and os.path.isfile(LARGE_HIC):
        argv.append(LARGE_HIC)
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen",
               PYTHONPATH=os.pathsep.join(p for p in (GUI_DIR, os.environ.get("HICX_PYTHON_MODULE_DIR", ""),
                                                      os.environ.get("PYTHONPATH", "")) if p))
    proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                          env=env, timeout=1200)
    assert proc.returncode == 0, proc.stderr[-3000:]
    result = json.loads(proc.stdout.strip().splitlines()[-1])
    child_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    print("scripted session:", result, "children peak kB:", child_peak)
    assert result["fetches"] >= 20
    assert result["peak_rss_kb"] <= SESSION_RSS_BOUND_KB, result
