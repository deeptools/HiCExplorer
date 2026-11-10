import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import apieq  # noqa: E402
import reference  # noqa: E402


@pytest.fixture(scope='session')
def reference_results():
    """The 3.7.6 records of every case, computed once per session (cached)."""
    try:
        return reference.results_for(apieq.cases_files(HERE))
    except reference.Unavailable as reason:
        pytest.skip('3.7.6 reference unavailable: {}'.format(reason))
