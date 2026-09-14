"""python -m hicexplorer_plot TOOL DATA.json [--remove-data]

Draws the figure of TOOL from the data file its C++ binary wrote. With
--remove-data the data file is deleted once it has been read; the binaries
pass it for the temporary file they create, and a file written with the C++
option --plotData is kept.
"""

import importlib
import json
import os
import sys


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    remove = "--remove-data" in argv
    argv = [a for a in argv if a != "--remove-data"]
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    tool, data_path = argv
    try:
        module = importlib.import_module("hicexplorer_plot." + tool)
    except ImportError as error:
        sys.stderr.write("hicexplorer_plot: no drawing module for {}: {}\n".format(tool, error))
        return 2
    try:
        with open(data_path) as handle:
            data = json.load(handle)
    finally:
        if remove:
            try:
                os.remove(data_path)
            except OSError:
                pass
    result = module.draw(data)
    return 0 if result is None else int(result)


if __name__ == "__main__":
    sys.exit(main())
