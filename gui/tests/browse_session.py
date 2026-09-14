"""A scripted matrix browser session for the memory measurement of PLAN 10.5.

    QT_QPA_PLATFORM=offscreen python browse_session.py COOL [LARGE_HIC]

Browses gm12878_chr1.cool (regions, zooms, pans at its single resolution)
and, when given, a large .hic file (whole chromosome down to a few Mb, with
resolution switching). Prints one JSON object: the number of fetches, the
largest array shape and the peak RSS of this process in kB.
"""

import json
import resource
import sys
import time

from PySide6 import QtWidgets

from hicexplorer_gui.browser import MatrixBrowser


def wait(app, browser, timeout=600):
    deadline = time.time() + timeout
    while time.time() < deadline:
        app.processEvents()
        if not browser.fetcher.busy and not browser.timer.isActive():
            app.processEvents()
            if not browser.fetcher.busy and not browser.timer.isActive():
                return
        time.sleep(0.01)
    raise RuntimeError("browser did not settle")


def main():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    browser = MatrixBrowser()
    browser.resize(1000, 1000)
    browser.show()
    fetches, largest = 0, (0, 0)

    def record():
        nonlocal fetches, largest
        fetches += 1
        for array in browser.shown:
            if array.size > largest[0] * largest[1]:
                largest = array.shape

    browser.fetch_finished.connect(record)
    browser.open_matrix(sys.argv[1])
    wait(app, browser)
    for region in ("1:1,000,000-11,000,000", "1:100,000,000-108,000,000", "1:200,000,000-202,000,000"):
        browser.goto(region)
        wait(app, browser)
        for _ in range(3):
            browser.pan(0.3)
            wait(app, browser)
        browser.zoom(0.5)
        wait(app, browser)
        browser.zoom(1.5)
        wait(app, browser)
    if len(sys.argv) > 2:
        browser.open_matrix(sys.argv[2])
        wait(app, browser)
        browser.goto("1")
        wait(app, browser)
        for _ in range(6):
            browser.zoom(0.35)
            wait(app, browser)
            browser.pan(0.2)
            wait(app, browser)
        browser.normalization.setCurrentText("KR")
        wait(app, browser)
        for region in ("1:150,000,000-151,000,000", "2:10,000,000-30,000,000", "X:50,000,000-52,000,000"):
            browser.goto(region)
            wait(app, browser)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({"fetches": fetches, "largest_array": list(largest), "peak_rss_kb": peak,
                      "stale_dropped": browser.fetcher.stale}))
    browser.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
