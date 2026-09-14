"""The polarization plot of hicCompartmentalization.

The calls are those of hicexplorer/hicCompartmentalization.py
plot_polarization_ratio and the rcParams line of main(), unchanged. The C++
binary writes --outputMatrix and the _dat file and passes the polarization
ratios of every matrix.

Data:
    outputFileName   --outputFileName, the figure
    labels           per matrix, its file name without the last extension
    quantile         --quantile
    ratios           per matrix, within / between per quantile
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def plot_polarization_ratio(polarization_ratio, plotName, labels,
                            number_of_quantiles):
    for i, r in enumerate(polarization_ratio):
        plt.plot(r, marker="o", label=labels[i])
    plt.axhline(1, c='grey', ls='--', lw=1)
    plt.axvline(number_of_quantiles / 2, c='grey', ls='--', lw=1)
    plt.legend(loc='best')
    plt.xlabel('Quantiles')
    plt.ylabel('signal within comp. / signla between comp.')
    plt.title('compartment polarization ratio')
    plt.savefig(plotName)


def draw(data):
    matplotlib.rcParams['pdf.fonttype'] = 42
    ratios = [[np.float64(value) for value in row] for row in data["ratios"]]
    plot_polarization_ratio(ratios, data["outputFileName"], data["labels"],
                            data["quantile"])
