"""Console scripts named after the plotting tools.

Each runs the C++ binary of the same name, which parses the command line,
computes the figure's data and then draws it through this package with the
interpreter these scripts run in. The binary is looked up in HICX_CPP_BIN,
then on PATH.
"""

import os
import shutil
import sys


def _run(tool):
    directory = os.environ.get("HICX_CPP_BIN")
    binary = os.path.join(directory, tool) if directory else shutil.which(tool)
    if not binary or not os.access(binary, os.X_OK):
        sys.stderr.write("{}: the C++ binary was not found (set HICX_CPP_BIN to the "
                         "directory of the HiCExplorer C++ tools)\n".format(tool))
        return 2
    env = dict(os.environ)
    env.setdefault("HICX_PLOT_PYTHON", sys.executable)
    package_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env["PYTHONPATH"] = os.pathsep.join(p for p in (package_root, env.get("PYTHONPATH")) if p)
    os.execve(binary, [binary] + sys.argv[1:], env)
    return 1


def hicPlotDistVsCounts():
    return _run("hicPlotDistVsCounts")


def hicPlotViewpoint():
    return _run("hicPlotViewpoint")
