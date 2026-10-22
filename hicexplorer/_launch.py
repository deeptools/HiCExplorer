"""Console-script targets that start a C++ executable without importing the Python stack."""
from hicexplorer._cpp import entry_point


def __getattr__(name):
    # Only tool names resolve; any other lookup (pytest_plugins, __path__ ...) must fail normally.
    if not name.startswith(('hic', 'chic')):
        raise AttributeError(name)
    return entry_point(name)
