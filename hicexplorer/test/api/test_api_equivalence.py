"""Each public function against its HiCExplorer 3.7.6 reference, one test per
case of the cases_*.py files in this directory (apieq.py explains the
comparison, cpp/python/API_BINDINGS.md how to add cases)."""
import os

import pytest

import apieq

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = [(path, case) for path in apieq.cases_files(HERE)
         for case in apieq.load_cases_module(path).CASES]


@pytest.mark.parametrize('path,case', CASES,
                         ids=['{}::{}'.format(os.path.basename(p)[6:-3], c.id) for p, c in CASES])
def test_matches_3_7_6(path, case, reference_results, tmp_path):
    expected = reference_results[(os.path.basename(path), case.id)]
    module = apieq.load_cases_module(path)
    actual = apieq.run_case(module, case, str(tmp_path / apieq.safe_name(case.id)))
    message = apieq.compare_records(expected, actual, case)
    assert message is None, message
