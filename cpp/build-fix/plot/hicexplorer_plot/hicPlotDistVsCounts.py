"""The figure of hicPlotDistVsCounts.

The calls are those of hicexplorer/hicPlotDistVsCounts.py main() from
`mpl.rcParams['pdf.fonttype'] = 42` to `plt.close(fig)`, unchanged. The C++
binary computes compute_distance_mean, the scale factors and the scaled
curves, and writes, per positional matrix in command line order, the curves
the loop draws (a curve with at most one point is skipped there and not
written).

Data:
    plotFile     the figure file
    perchr       --perchr
    num_chroms   len(chroms), the chromosomes with more than one distance
    maxdepth     --maxdepth
    plotsize     --plotsize or null
    matrices     [{"label": labels[matrix_file],
                   "series": [{"chrom", "x": distances, "y": scaled means}]}]
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib as mpl  # noqa: E402
import numpy as np  # noqa: E402


def draw(data):
    mpl.rcParams['pdf.fonttype'] = 42
    matrices = data["matrices"]
    perchr = data["perchr"]
    maxdepth = data["maxdepth"]

    if len(matrices) > 1 and perchr:
        max_cols = 4
        num_rows = int(np.ceil(float(data["num_chroms"]) / max_cols))
        num_cols = min(data["num_chroms"], max_cols)
    else:
        num_cols = num_rows = 1

    if data["plotsize"] is None:
        width = 6
        height = 4
    else:
        width, height = data["plotsize"]
    fig = plt.figure(figsize=(width * num_cols, height * num_rows))

    axs = np.empty((num_rows, num_cols), dtype='object')
    for matrix in matrices:
        idx = 0
        for series in matrix["series"]:
            chrom = series["chrom"]
            x = tuple(series["x"])
            if perchr and len(matrices) == 1:
                col = 0
                row = 0
            else:
                col = idx % num_cols
                row = idx // num_cols
            if axs[row, col] is None:
                ax = plt.subplot2grid((num_rows, num_cols), (row, col))
                ax.set_xlabel('genomic distance')
                ax.set_ylabel('corrected Hi-C counts')
                try:
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                except ValueError:
                    continue
            else:
                ax = axs[row, col]
            y = np.array(series["y"], dtype=np.float64)
            if perchr and len(matrices) > 1:
                label = matrix["label"]
                ax.set_title(chrom)
            elif perchr:
                label = chrom
            else:
                label = matrix["label"]

            ax.plot(x, y, label=label)
            axs[row, col] = ax
            idx += 1

    for ax in axs.reshape(-1):
        if ax is None:
            continue
        ax.legend(prop={'size': 'small'})
        ax.set_xlim(0, maxdepth)
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(data["plotFile"], bbox_inches='tight', bbox_extra_artists=(lgd,))
    plt.close(fig)
