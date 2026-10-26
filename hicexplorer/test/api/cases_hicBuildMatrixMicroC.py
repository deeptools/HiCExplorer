"""Equivalence cases for hicexplorer.hicBuildMatrixMicroC (apieq.py).

The module's own public name is the class Tester; the functions it takes
from hicexplorer.lib.buildMatrixMethods by its star import are covered by
cases_lib_buildMatrixMethods.py, and one of them is called here through this
module to show that it is the same function.
"""
import os

from apieq import Case, File

M = 'hicexplorer.hicBuildMatrixMicroC.'
SUFFIX = os.path.join('hicexplorer', 'test', 'test_data') + '/'


def tester(Tester, environment=None):
    """root and bam_file_1, the package location cut off; with
    HIC_TEST_DATA_DIR set to `environment` when given."""
    saved = os.environ.get('HIC_TEST_DATA_DIR')
    if environment is None:
        os.environ.pop('HIC_TEST_DATA_DIR', None)
    else:
        os.environ['HIC_TEST_DATA_DIR'] = environment
    try:
        t = Tester()
    finally:
        if saved is None:
            os.environ.pop('HIC_TEST_DATA_DIR', None)
        else:
            os.environ['HIC_TEST_DATA_DIR'] = saved
    if environment:
        return t.root, t.bam_file_1, sorted(vars(t))
    return (t.root.endswith(SUFFIX), t.root[-len(SUFFIX):], t.bam_file_1[-len(SUFFIX) - 7:], sorted(vars(t)))


def star_imported(get_bins):
    return get_bins(1000, [('chr1', 2500), ('chr2', 999)])


CASES = [
    Case('Tester/default', M + 'Tester', driver='tester'),
    Case('Tester/environment', M + 'Tester', args=(File(),), driver='tester'),
    Case('Tester/environment_trailing_slash', M + 'Tester', args=('/data/hic/',), driver='tester'),
    Case('Tester/environment_empty', M + 'Tester', args=('',), driver='tester'),
    Case('get_bins/star_import', M + 'get_bins', driver='star_imported'),
]
