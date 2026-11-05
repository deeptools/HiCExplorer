import logging
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.getLogger('hicmatrix').setLevel(logging.ERROR)
# logging.getLogger('hicmatrix').setLevel(logging.DEBUG)

logging.getLogger('numexpr').setLevel(logging.ERROR)


import warnings
import sys

if not sys.warnoptions:
    warnings.simplefilter("ignore")

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import multiprocessing as mp
mp.set_start_method('fork')


def _find_extension():
    """Make hicexplorer._hicx, the C++ implementation of the public API,
    importable from a source checkout.

    An installed package has the extension next to this file. A source
    checkout has it in the build tree, <build>/python/hicexplorer, which is
    added to the package path: from HICX_PYTHON_EXT_DIR, from the build tree
    of HICX_BIN_DIR (<build>/tools), or from cpp/build of this checkout."""
    import os

    def has_extension(directory):
        try:
            return any(name.startswith('_hicx.') for name in os.listdir(directory))
        except OSError:
            return False

    here = os.path.dirname(os.path.abspath(__file__))
    if has_extension(here):
        return
    candidates = [os.environ.get('HICX_PYTHON_EXT_DIR')]
    bin_dir = os.environ.get('HICX_BIN_DIR')
    if bin_dir:
        candidates.append(os.path.join(os.path.dirname(os.path.abspath(bin_dir)), 'python', 'hicexplorer'))
    candidates.append(os.path.join(os.path.dirname(here), 'cpp', 'build', 'python', 'hicexplorer'))
    for directory in candidates:
        if directory and has_extension(directory):
            __path__.append(os.path.abspath(directory))
            return


_find_extension()
