"""Console-script targets that start a C++ executable without importing the Python stack."""
from hicexplorer._cpp import entry_point


def __getattr__(name):
    if name.startswith('_'):
        raise AttributeError(name)
    return entry_point(name)
