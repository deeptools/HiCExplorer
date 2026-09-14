"""python -m hicexplorer_plot TOOL DATA.json [--remove-data]
python -m hicexplorer_plot --check TOOL

Draws the figure of TOOL from the data file its C++ binary wrote. With
--remove-data the data file is deleted once it has been read; the binaries
pass it for the temporary file they create, and a file written with the C++
option --plotData is kept.

Before drawing, and alone with --check, the renderer versions are checked
against the pins of hicexplorer_plot.environment; a mismatch or a missing
package exits with status 3 before any figure is written.
"""

import importlib
import json
import os
import sys

from hicexplorer_plot import environment


def _remove(paths):
    for path in paths:
        try:
            os.remove(path)
        except OSError:
            pass


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) == 2 and argv[0] == "--check":
        return 0 if environment.check(argv[1]) else environment.REFUSAL_EXIT
    remove = "--remove-data" in argv
    argv = [a for a in argv if a != "--remove-data"]
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    tool, data_path = argv
    try:
        with open(data_path) as handle:
            data = json.load(handle)
    finally:
        if remove:
            _remove([data_path])
    # Arrays too large for JSON travel in files the document names; with
    # --remove-data they are temporary too and go once the figure is drawn.
    temporary = data.get("temporary_files", []) if remove else []
    if not environment.check(tool):
        _remove(temporary)
        return environment.REFUSAL_EXIT
    try:
        module = importlib.import_module("hicexplorer_plot." + tool)
    except ImportError as error:
        _remove(temporary)
        sys.stderr.write("hicexplorer_plot: no drawing module for {}: {}\n".format(tool, error))
        return 2
    try:
        result = module.draw(data)
    finally:
        _remove(temporary)
    return 0 if result is None else int(result)


if __name__ == "__main__":
    sys.exit(main())
