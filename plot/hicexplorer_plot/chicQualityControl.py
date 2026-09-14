"""The sparsity and histogram figures of chicQualityControl.

The calls are those at the end of hicexplorer/chicQualityControl.py main(),
unchanged. The C++ binary writes the tables and passes, per matrix, the
sparsity of every reference point that is faulty (-1.0) in no matrix, in
file order.

Data:
    outFileNameSparsity    the sparsity distribution figure
    outFileNameHistogram   the histogram figure
    dpi                    --dpi
    sparsity               the --sparsity threshold
    labels                 the basename of every matrix
    x                      per matrix, the sparsity values of the kept points
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def draw(data):
    labels = data["labels"]
    x = [np.array(values, dtype=np.float64) for values in data["x"]]
    y = [[i] * len(x[i]) for i in range(len(labels))]

    for i in range(len(labels)):
        plt.plot(x[i], y[i], 'o', mfc='none', markersize=0.3,
                 label=labels[i])
    plt.yticks([])
    plt.xlabel("Sparsity level")

    plt.axvline(x=data["sparsity"], c='r', label='sparsity threshold', linewidth=0.3)
    plt.xscale('log')
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    plt.legend(loc='center', bbox_to_anchor=(1.4, 0.5))
    plt.savefig(data["outFileNameSparsity"], dpi=data["dpi"])

    plt.close()
    for i in range(len(labels)):
        plt.hist(x[i], bins=100, alpha=0.5, label=labels[i])
    plt.xlabel("Sparsity level")
    plt.ylabel("Number of counts")

    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    plt.legend(loc='center', bbox_to_anchor=(1.4, 0.5))
    plt.savefig(data["outFileNameHistogram"], dpi=data["dpi"])
