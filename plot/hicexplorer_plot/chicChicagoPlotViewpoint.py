"""The figure of chicChicagoPlotViewpoint.

chicChicagoPlotViewpoint has no Python original; this reproduces R
Chicago::plotBaits() (checked directly against the real installed package):
x is distSign, y is the raw N column, score only picks a point's colour
(two significance tiers, chicChicagoPlotViewpoint.cpp's own --plevel1/2,
R's own defaults 5 and 3), a grey vertical line marks the bait at x = 0, and
an optional Brownian mean line with a dashed upper 95% band overlays the
background model when the C++ side was given one.

Data:
    plotFile      the figure file
    baitLabel     the plot title
    plevel1       more significant score threshold (red)
    plevel2       less significant score threshold (blue)
    x             distSign per kept row
    y             N per kept row
    score         score per kept row, same order as x/y
    dpi           --dpi
    hasBackground whether bmeanX/bmeanY/bmeanUpperY are present
    bmeanX        distSign, sorted ascending
    bmeanY        Bmean at each bmeanX
    bmeanUpperY   Bmean + 1.96 * sqrt(Bmean + Bmean^2 / dispersion) at each bmeanX
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402
import numpy as np  # noqa: E402


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42

    x = np.array(data["x"], dtype=np.float64)
    y = np.array(data["y"], dtype=np.float64)
    score = np.array(data["score"], dtype=np.float64)
    plevel1 = data["plevel1"]
    plevel2 = data["plevel2"]

    fig, ax = plt.subplots(figsize=(8, 5))

    if data.get("hasBackground"):
        bx = np.array(data["bmeanX"], dtype=np.float64)
        by = np.array(data["bmeanY"], dtype=np.float64)
        bu = np.array(data["bmeanUpperY"], dtype=np.float64)
        ax.plot(bx, by, color='darkgrey', linewidth=1, label='Brownian mean')
        ax.plot(bx, bu, color='darkgrey', linewidth=1, linestyle='--',
                label='Brownian 95% upper')

    background = score < plevel2
    lev2 = (score >= plevel2) & (score < plevel1)
    lev1 = score >= plevel1

    if np.any(background):
        ax.scatter(x[background], y[background], s=10, c='black', label='background')
    if np.any(lev2):
        ax.scatter(x[lev2], y[lev2], s=20, c='blue', label='score >= {}'.format(plevel2))
    if np.any(lev1):
        ax.scatter(x[lev1], y[lev1], s=32, c='red', label='score >= {}'.format(plevel1))

    ax.axvline(0, color='grey', linewidth=1)
    ax.set_xlabel('distance from bait (bp)')
    ax.set_ylabel('N')
    ax.set_title(data["baitLabel"])
    ax.legend(prop={'size': 'small'}, loc='upper right')
    plt.tight_layout()
    plt.savefig(data["plotFile"], dpi=data.get("dpi"))
    plt.close(fig)
