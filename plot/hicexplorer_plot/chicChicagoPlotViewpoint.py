"""The figure of chicChicagoPlotViewpoint.

chicChicagoPlotViewpoint has no Python original. --style scatter reproduces
R Chicago::plotBaits() (checked directly against the real installed
package): x is distSign, y is the raw N column, score only picks a point's
colour (two significance tiers, chicChicagoPlotViewpoint.cpp's own
--plevel1/2, R's own defaults 5 and 3), a grey vertical line marks the bait
at x = 0, and an optional Brownian mean line with a dashed upper 95% band
overlays the background model when the C++ side was given one. --style arcs
(this project's own convention, no R precedent) draws one arc per kept
interaction, in either --baitID's own relative frame or --region's absolute
genomic coordinates.

Data (all styles):
    plotFile, style, title, xLabel, plevel1, plevel2, dpi, score
Data (style == "scatter"):
    x, y, hasBackground, bmeanX, bmeanY, bmeanUpperY (last three only when
    hasBackground)
Data (style == "arcs"):
    arcX1, arcX2 (the two ends of each arc), arcHeight (log1p(N)), anchors
    (bait positions to mark with a small tick)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402
import numpy as np  # noqa: E402


def _tier_masks(score, plevel1, plevel2):
    score = np.asarray(score, dtype=np.float64)
    background = score < plevel2
    lev2 = (score >= plevel2) & (score < plevel1)
    lev1 = score >= plevel1
    return background, lev2, lev1


def _draw_scatter(ax, data):
    x = np.array(data["x"], dtype=np.float64)
    y = np.array(data["y"], dtype=np.float64)
    score = data["score"]
    plevel1, plevel2 = data["plevel1"], data["plevel2"]

    if data.get("hasBackground"):
        bx = np.array(data["bmeanX"], dtype=np.float64)
        by = np.array(data["bmeanY"], dtype=np.float64)
        bu = np.array(data["bmeanUpperY"], dtype=np.float64)
        ax.plot(bx, by, color='darkgrey', linewidth=1, label='Brownian mean')
        ax.plot(bx, bu, color='darkgrey', linewidth=1, linestyle='--',
                label='Brownian 95% upper')

    background, lev2, lev1 = _tier_masks(score, plevel1, plevel2)
    if np.any(background):
        ax.scatter(x[background], y[background], s=10, c='black', label='background')
    if np.any(lev2):
        ax.scatter(x[lev2], y[lev2], s=20, c='blue', label='score >= {}'.format(plevel2))
    if np.any(lev1):
        ax.scatter(x[lev1], y[lev1], s=32, c='red', label='score >= {}'.format(plevel1))

    ax.axvline(0, color='grey', linewidth=1)
    ax.set_ylabel('N')


def _arc_points(x1, x2, height, n=60):
    """A quadratic Bezier from (x1, 0) to (x2, 0) peaking at `height` at the
    midpoint, the usual genome-browser arc shape."""
    t = np.linspace(0.0, 1.0, n)
    mx = (x1 + x2) / 2.0
    xs = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * mx + t ** 2 * x2
    ys = 2 * (1 - t) * t * height
    return xs, ys


def _draw_arcs(ax, data):
    x1 = np.array(data["arcX1"], dtype=np.float64)
    x2 = np.array(data["arcX2"], dtype=np.float64)
    heights = np.array(data["arcHeight"], dtype=np.float64)
    score = data["score"]
    plevel1, plevel2 = data["plevel1"], data["plevel2"]
    background, lev2, lev1 = _tier_masks(score, plevel1, plevel2)
    colour = np.where(lev1, 'red', np.where(lev2, 'blue', 'black'))
    zorder = np.where(lev1, 3, np.where(lev2, 2, 1))

    order = np.argsort(zorder)
    for i in order:
        xs, ys = _arc_points(x1[i], x2[i], heights[i])
        ax.plot(xs, ys, color=colour[i], linewidth=1.2, alpha=0.8, zorder=int(zorder[i]))

    for anchor in data.get("anchors", []):
        ax.plot([anchor, anchor], [0, -0.03 * max(heights.max(initial=1.0), 1.0)],
                color='grey', linewidth=1.5)

    for label, c in (('background', 'black'), ('score >= {}'.format(plevel2), 'blue'),
                     ('score >= {}'.format(plevel1), 'red')):
        ax.plot([], [], color=c, label=label)
    ax.axhline(0, color='lightgrey', linewidth=0.8, zorder=0)
    ax.set_ylabel('log1p(N)')
    ax.set_ylim(bottom=min(-0.05 * max(heights.max(initial=1.0), 1.0), ax.get_ylim()[0]))


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(9, 5))

    if data["style"] == "arcs":
        _draw_arcs(ax, data)
    else:
        _draw_scatter(ax, data)

    ax.set_xlabel(data["xLabel"])
    ax.set_title(data["title"])
    ax.legend(prop={'size': 'small'}, loc='upper right')
    plt.tight_layout()
    plt.savefig(data["plotFile"], dpi=data.get("dpi"))
    plt.close(fig)
