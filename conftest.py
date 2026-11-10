# numpy 2 writes scalars as np.float64(0.75) where numpy 1 wrote 0.75; the
# doctests show the numpy 1 form, which is also what HiCExplorer 3 printed.
import numpy


def pytest_configure(config):
    if int(numpy.__version__.split('.')[0]) >= 2:
        numpy.set_printoptions(legacy='1.25')
