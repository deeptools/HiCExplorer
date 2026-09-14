"""The figure of hicCorrectMatrix diagnostic_plot.

The calls are those of hicexplorer/hicCorrectMatrix.py plot_total_contact_dist
and the rcParams line of main(), unchanged. The C++ binary loads the matrix,
masks the zero coverage bins, and computes for every panel the coverage per
bin without the diagonal, its MAD (median of the positive values and median
absolute deviation) and the modified z-score filter below 5.

Data:
    plotName   --plotName
    xMax       --xMax or null
    perchr     --perchr
    panels     [{"title": chromosome name or null, "row_sum": the kept coverage,
                 "median", "med_abs_deviation"}]
"""

import matplotlib
from matplotlib import use
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.ticker import MultipleLocator, FormatStrFormatter  # noqa: E402
import numpy as np  # noqa: E402


class MAD(object):
    """hicCorrectMatrix.MAD reduced to what the figure reads, value_to_mad,
    with the median and the median absolute deviation the C++ computed."""

    def __init__(self, median, med_abs_deviation):
        self.mad_b_value = 0.6745
        self.median = np.float64(median)
        self.med_abs_deviation = np.float64(med_abs_deviation)

    def value_to_mad(self, value):
        diff = value - self.median
        if self.med_abs_deviation == 0.0:
            return self.mad_b_value * diff

        return self.mad_b_value * diff / self.med_abs_deviation


class _Args(object):
    pass


def plot_total_contact_dist(panels, args):
    use('Agg')

    majorlocator = MultipleLocator(1)
    majorformatter = FormatStrFormatter('%d')
    minorlocator = MultipleLocator(0.2)

    def plot_histogram(row_sum_values, mad_values, ax1, title=None):

        if args.xMax:
            ax1.set_xlim(ax1.get_xlim()[0], args.xMax)
            row_sum_values = row_sum_values[row_sum_values < args.xMax]

        ax1.set_xlabel("total counts per bin")
        ax1.set_ylabel("frequency")
        ax1.patch.set_visible(False)
        dist, bin_s, __ = ax1.hist(row_sum_values, 100, color='green')

        # add second axis on top
        ax2 = ax1.twiny()
        ax2.set_xlabel("modified z-score")
        ax2.xaxis.set_major_locator(majorlocator)
        ax2.xaxis.set_major_formatter(majorformatter)
        ax2.xaxis.grid(True, which='minor')
        # for the minor ticks, use no labels; default NullFormatter
        ax2.xaxis.set_minor_locator(minorlocator)

        ax2.set_xlim(mad_values.value_to_mad(np.array(ax1.get_xlim())))

        # get first local mininum value
        local_min = [x for x, y in enumerate(dist) if 1 <= x < len(
            dist) - 1 and dist[x - 1] > y < dist[x + 1]]

        if len(local_min) > 0:
            threshold = bin_s[local_min[0]]
        else:
            threshold = None

        if threshold:
            mad_threshold = mad_values.value_to_mad(threshold)
            ymin, ymax = ax2.get_ylim()
            ax2.vlines(mad_threshold, ymin, ymax)

    if args.perchr:
        chroms = [panel["title"] for panel in panels]
        num_rows = int(np.ceil(float(len(chroms)) / 5))
        num_cols = min(len(chroms), 5)
        grids = gridspec.GridSpec(num_rows, num_cols)
        fig = plt.figure(figsize=(6 * num_cols, 5 * num_rows))
        ax = {}
        for plot_num, panel in enumerate(panels):
            chrname = panel["title"]
            row_sum = np.array(panel["row_sum"], dtype=np.float64)
            mad = MAD(panel["median"], panel["med_abs_deviation"])
            col = plot_num % num_cols
            row = plot_num // num_cols
            ax[chrname] = fig.add_subplot(grids[row, col])

            plot_histogram(row_sum, mad, ax[chrname], title=chrname)
            ax[chrname].set_title(chrname)
    else:
        fig = plt.figure()
        panel = panels[0]
        row_sum = np.array(panel["row_sum"], dtype=np.float64)
        mad = MAD(panel["median"], panel["med_abs_deviation"])
        ax = fig.add_subplot(111)
        plot_histogram(row_sum, mad, ax)

    plt.tight_layout()
    plt.savefig(args.plotName)
    plt.close()


def draw(data):
    matplotlib.rcParams['pdf.fonttype'] = 42
    args = _Args()
    args.plotName = data["plotName"]
    args.xMax = data["xMax"]
    args.perchr = data["perchr"]
    plot_total_contact_dist(data["panels"], args)
