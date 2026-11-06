"""The ratio plot of hicInterIntraTAD.

The calls are those at the end of hicexplorer/hicInterIntraTAD.py main(),
unchanged. The C++ binary writes the table and passes the two ratio columns
in table order.

Data:
    plotFile   --outFileNameRatioPlot
    fontsize   --fontsize
    dpi        --dpi
    x          inter_left_intra_ratio per TAD
    y          inter_right_intra_ratio per TAD
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


def draw(data):
    plt.scatter(data["x"], data["y"], s=20, alpha=0.7)
    plt.xlabel('Inter-left/intra TAD contact ratio', fontsize=data["fontsize"])
    plt.ylabel('Inter-right/intra TAD contact ratio', fontsize=data["fontsize"])
    plt.tight_layout()
    plt.savefig(data["plotFile"], dpi=data["dpi"])
    plt.close()
